#!/usr/bin/env python3
"""Run the same md-flexible case with reference md-flexible and AutoPasDM.

The comparison does not require identical particle owner ranks. Inverted Pressure
uses measured local computation times, so exact adaptive boundaries can differ
between implementations and runs. This regression checks that

1. the final owned-particle state agrees by particle ID within tolerance, and
2. optionally, both runs actually move at least one subdomain boundary.

Example:
    python3 scripts/run_md_flexible_comparison.py \
        --reference /path/to/AutoPas/build/examples/md-flexible/md-flexible \
        --distributed ./md-flexible \
        --config input/autopasdm_load_balancing_compare.yaml \
        --ranks 8 --require-load-balancing

For SLURM, e.g. use:
    --launcher 'srun -n {ranks}'
"""

from __future__ import annotations

import argparse
import glob
import math
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


PARTICLE_FILE_RE = re.compile(r"_Particles_(\d+)_(\d+)\.vtu$")
RANK_FILE_RE = re.compile(r"_Ranks_(\d+)_(\d+)\.vtu$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True, type=Path, help="Original md-flexible executable")
    parser.add_argument("--distributed", required=True, type=Path, help="AutoPasDM md-flexible executable")
    parser.add_argument("--config", required=True, type=Path, help="YAML input shared by both runs")
    parser.add_argument("--ranks", type=int, default=8, help="MPI ranks (default: 8)")
    parser.add_argument("--reference-ranks", type=int, help="Override the reference MPI rank count")
    parser.add_argument("--distributed-ranks", type=int, help="Override the AutoPasDM MPI rank count")
    parser.add_argument(
        "--require-load-balancing",
        action="store_true",
        help="require the decomposition boundaries to move in both runs",
    )
    parser.add_argument(
        "--launcher",
        default="mpirun -np {ranks}",
        help="MPI launcher command. {ranks} is replaced with --ranks (default: 'mpirun -np {ranks}')",
    )
    parser.add_argument("--absolute-tolerance", type=float, default=1e-3)
    parser.add_argument("--relative-tolerance", type=float, default=1e-4)
    parser.add_argument("--boundary-move-tolerance", type=float, default=1e-6)
    parser.add_argument(
        "--keep-output",
        type=Path,
        help=(
            "Keep both run directories below this directory instead of deleting temporary output. "
            "The destination must be absent or empty; existing contents are never removed."
        ),
    )
    return parser.parse_args()


def absolute_existing_file(path: Path, description: str) -> Path:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise RuntimeError(f"{description} does not exist: {path}")
    return path


def launcher_command(template: str, ranks: int) -> list[str]:
    if ranks <= 0:
        raise RuntimeError("--ranks must be greater than zero")
    return shlex.split(template.replace("{ranks}", str(ranks)))


def run_case(label: str, executable: Path, config: Path, run_directory: Path, launcher: list[str]) -> None:
    run_directory.mkdir(parents=True, exist_ok=True)
    command = [*launcher, str(executable), "--yaml-filename", str(config)]
    log_file = run_directory / "run.log"

    print(f"Running {label}:")
    print("  " + shlex.join(command))
    print(f"  cwd: {run_directory}")

    with log_file.open("w", encoding="utf-8") as log:
        result = subprocess.run(command, cwd=run_directory, stdout=log, stderr=subprocess.STDOUT, check=False)

    if result.returncode != 0:
        tail = log_file.read_text(encoding="utf-8", errors="replace").splitlines()[-40:]
        raise RuntimeError(
            f"{label} failed with exit code {result.returncode}. Last log lines:\n" + "\n".join(tail)
        )


def data_array(piece: ET.Element, name: str) -> list[str]:
    for element in piece.iter("DataArray"):
        if element.attrib.get("Name") == name:
            return (element.text or "").split()
    raise RuntimeError(f'DataArray "{name}" not found')


def optional_data_array(piece: ET.Element, name: str) -> list[str] | None:
    for element in piece.iter("DataArray"):
        if element.attrib.get("Name") == name:
            return (element.text or "").split()
    return None


def point_coordinates(piece: ET.Element) -> list[float]:
    points = next(piece.iter("Points"), None)
    if points is None:
        raise RuntimeError("Points element not found")
    array = next(points.iter("DataArray"), None)
    if array is None:
        raise RuntimeError("Points DataArray not found")
    return [float(value) for value in (array.text or "").split()]


def matching_files(root: Path, pattern: re.Pattern[str]) -> list[tuple[int, int, Path]]:
    matches: list[tuple[int, int, Path]] = []
    for filename in glob.glob(str(root / "**" / "*.vtu"), recursive=True):
        match = pattern.search(filename)
        if match:
            rank = int(match.group(1))
            iteration = int(match.group(2))
            matches.append((iteration, rank, Path(filename)))
    if not matches:
        raise RuntimeError(f"No matching VTK files found below {root}")
    return matches


def first_and_final_iterations(root: Path, pattern: re.Pattern[str]) -> tuple[int, int]:
    iterations = [iteration for iteration, _, _ in matching_files(root, pattern)]
    return min(iterations), max(iterations)


def files_at_iteration(root: Path, pattern: re.Pattern[str], iteration: int) -> list[tuple[int, Path]]:
    return sorted(
        (rank, filename)
        for file_iteration, rank, filename in matching_files(root, pattern)
        if file_iteration == iteration
    )


def read_particle_state(root: Path) -> tuple[int, dict[int, dict[str, object]]]:
    _, final_iteration = first_and_final_iterations(root, PARTICLE_FILE_RE)
    result: dict[int, dict[str, object]] = {}

    for rank, filename in files_at_iteration(root, PARTICLE_FILE_RE, final_iteration):
        tree = ET.parse(filename)
        piece = next(tree.getroot().iter("Piece"))

        ids = [int(value) for value in data_array(piece, "ids")]
        type_ids = [int(value) for value in data_array(piece, "typeIds")]
        if len(type_ids) != len(ids):
            raise RuntimeError(f"Type-ID/ID length mismatch in {filename}")

        fields: dict[str, list[float]] = {"position": point_coordinates(piece)}
        required_vector_arrays = {
            "velocity": ("velocities", 3),
            "force": ("forces", 3),
        }
        optional_vector_arrays = {
            "quaternion": ("quaternions", 4),
            "angularVelocity": ("angularVelocities", 3),
            "torque": ("torques", 3),
        }
        for field, (array_name, _) in required_vector_arrays.items():
            fields[field] = [float(value) for value in data_array(piece, array_name)]
        for field, (array_name, _) in optional_vector_arrays.items():
            raw_values = optional_data_array(piece, array_name)
            if raw_values is not None:
                fields[field] = [float(value) for value in raw_values]

        for field, values in fields.items():
            components = 4 if field == "quaternion" else 3
            if len(values) != components * len(ids):
                raise RuntimeError(f"{field}/ID length mismatch in {filename}")

        for index, particle_id in enumerate(ids):
            if particle_id in result:
                raise RuntimeError(f"Particle ID {particle_id} occurs more than once in final owned output")
            particle_state: dict[str, object] = {"rank": rank, "type": type_ids[index]}
            for field, values in fields.items():
                components = 4 if field == "quaternion" else 3
                begin = components * index
                particle_state[field] = tuple(values[begin : begin + components])
            result[particle_id] = particle_state

    return final_iteration, result


def read_boxes(root: Path, iteration: int) -> dict[int, tuple[tuple[float, ...], tuple[float, ...]]]:
    boxes: dict[int, tuple[tuple[float, ...], tuple[float, ...]]] = {}
    for rank, filename in files_at_iteration(root, RANK_FILE_RE, iteration):
        tree = ET.parse(filename)
        piece = next(tree.getroot().iter("Piece"))
        coordinates = point_coordinates(piece)
        if len(coordinates) != 24:
            raise RuntimeError(f"Expected 8 three-dimensional rank-box points in {filename}")
        points = [coordinates[index : index + 3] for index in range(0, len(coordinates), 3)]
        box_min = tuple(min(point[d] for point in points) for d in range(3))
        box_max = tuple(max(point[d] for point in points) for d in range(3))
        boxes[rank] = (box_min, box_max)
    return boxes


def maximum_boundary_movement(root: Path) -> tuple[int, int, float]:
    initial_iteration, final_iteration = first_and_final_iterations(root, RANK_FILE_RE)
    initial_boxes = read_boxes(root, initial_iteration)
    final_boxes = read_boxes(root, final_iteration)

    if set(initial_boxes) != set(final_boxes):
        raise RuntimeError("Rank IDs differ between initial and final domain output")

    movement = 0.0
    for rank in initial_boxes:
        for initial, final in zip(initial_boxes[rank][0] + initial_boxes[rank][1], final_boxes[rank][0] + final_boxes[rank][1]):
            movement = max(movement, abs(initial - final))
    return initial_iteration, final_iteration, movement


def close(a: float, b: float, absolute_tolerance: float, relative_tolerance: float) -> bool:
    return math.isclose(a, b, abs_tol=absolute_tolerance, rel_tol=relative_tolerance)


def compare_states(
    reference: dict[int, dict[str, object]],
    distributed: dict[int, dict[str, object]],
    absolute_tolerance: float,
    relative_tolerance: float,
) -> tuple[bool, int, float]:
    reference_ids = set(reference)
    distributed_ids = set(distributed)
    if reference_ids != distributed_ids:
        print("FAIL: final owned particle IDs differ.")
        print("  only in reference  :", sorted(reference_ids - distributed_ids))
        print("  only in distributed:", sorted(distributed_ids - reference_ids))
        return False, 0, math.inf

    ownership_changes = 0
    maximum_difference = 0.0
    failures: list[str] = []

    for particle_id in sorted(reference_ids):
        expected = reference[particle_id]
        actual = distributed[particle_id]

        if expected["rank"] != actual["rank"]:
            ownership_changes += 1
        if expected["type"] != actual["type"]:
            failures.append(f"id={particle_id}: type {expected['type']} != {actual['type']}")

        for field in ("position", "velocity", "force", "quaternion", "angularVelocity", "torque"):
            if (field in expected) != (field in actual):
                failures.append(f"id={particle_id}: {field} is present in only one output")
                continue
            if field not in expected:
                continue
            expected_vector = expected[field]
            actual_vector = actual[field]
            if field == "quaternion":
                direct_difference = max(abs(a - b) for a, b in zip(expected_vector, actual_vector))
                negated_difference = max(abs(a + b) for a, b in zip(expected_vector, actual_vector))
                if negated_difference < direct_difference:
                    actual_vector = tuple(-component for component in actual_vector)
            for component, (a, b) in enumerate(zip(expected_vector, actual_vector)):
                maximum_difference = max(maximum_difference, abs(a - b))
                if not close(a, b, absolute_tolerance, relative_tolerance):
                    failures.append(
                        f"id={particle_id}: {field}[{component}] reference={a:.12g} distributed={b:.12g}"
                    )

    if failures:
        print("FAIL: final particle states differ beyond tolerance.")
        for failure in failures[:20]:
            print("  " + failure)
        if len(failures) > 20:
            print(f"  ... and {len(failures) - 20} more")
        return False, ownership_changes, maximum_difference

    return True, ownership_changes, maximum_difference


def main() -> int:
    args = parse_args()

    try:
        reference = absolute_existing_file(args.reference, "Reference executable")
        distributed = absolute_existing_file(args.distributed, "Distributed executable")
        config = absolute_existing_file(args.config, "Configuration")
        reference_ranks = args.reference_ranks if args.reference_ranks is not None else args.ranks
        distributed_ranks = args.distributed_ranks if args.distributed_ranks is not None else args.ranks
        reference_launcher = launcher_command(args.launcher, reference_ranks)
        distributed_launcher = launcher_command(args.launcher, distributed_ranks)

        if args.keep_output:
            work_root = args.keep_output.expanduser().resolve()
            if work_root.exists():
                if not work_root.is_dir():
                    raise RuntimeError(f"Output destination exists and is not a directory: {work_root}")
                if next(work_root.iterdir(), None) is not None:
                    raise RuntimeError(f"Output destination is not empty; refusing to overwrite: {work_root}")
            else:
                work_root.mkdir(parents=True)
            cleanup = None
        else:
            cleanup = tempfile.TemporaryDirectory(prefix="autopasdm_md_flexible_comparison_")
            work_root = Path(cleanup.name)

        reference_root = work_root / "reference"
        distributed_root = work_root / "distributed"

        run_case("reference md-flexible", reference, config, reference_root, reference_launcher)
        run_case("AutoPasDM md-flexible", distributed, config, distributed_root, distributed_launcher)

        boundaries_ok = True
        if args.require_load_balancing:
            ref_initial, ref_final, ref_movement = maximum_boundary_movement(reference_root)
            autopasdm_initial, autopasdm_final, autopasdm_movement = maximum_boundary_movement(distributed_root)

            print("\nAdaptive-domain check:")
            print(f"  reference   iterations {ref_initial}->{ref_final}, max boundary movement {ref_movement:.12g}")
            print(f"  distributed iterations {autopasdm_initial}->{autopasdm_final}, max boundary movement {autopasdm_movement:.12g}")

            if ref_movement <= args.boundary_move_tolerance:
                print("FAIL: reference md-flexible did not measurably change its decomposition.")
                boundaries_ok = False
            if autopasdm_movement <= args.boundary_move_tolerance:
                print("FAIL: AutoPasDM did not measurably change its decomposition.")
                boundaries_ok = False

        ref_iteration, ref_particles = read_particle_state(reference_root)
        autopasdm_iteration, autopasdm_particles = read_particle_state(distributed_root)

        print("\nFinal-state comparison:")
        print(f"  reference final iteration   : {ref_iteration}")
        print(f"  distributed final iteration : {autopasdm_iteration}")
        print(f"  reference owned particles   : {len(ref_particles)}")
        print(f"  distributed owned particles : {len(autopasdm_particles)}")

        if ref_iteration != autopasdm_iteration:
            print("FAIL: final iteration numbers differ; refusing to compare particle states.")
            states_ok = False
            ownership_changes = 0
            maximum_difference = math.inf
        else:
            states_ok, ownership_changes, maximum_difference = compare_states(
                ref_particles,
                autopasdm_particles,
                args.absolute_tolerance,
                args.relative_tolerance,
            )
        print(f"  particles owned by different ranks: {ownership_changes} (diagnostic only)")
        print(f"  maximum absolute component diff   : {maximum_difference:.12g}")
        print(f"  absolute tolerance                : {args.absolute_tolerance:.12g}")
        print(f"  relative tolerance                : {args.relative_tolerance:.12g}")

        if boundaries_ok and states_ok:
            if args.require_load_balancing:
                print("\nPASS: both decompositions adapt and the final particle state agrees within tolerance.")
            else:
                print("\nPASS: the final particle state agrees within tolerance.")
            result = 0
        else:
            print("\nFAIL: md-flexible comparison failed.")
            result = 1

        if args.keep_output:
            print(f"Output kept in: {work_root}")
        if cleanup is not None:
            cleanup.cleanup()
        return result

    except (OSError, RuntimeError, subprocess.SubprocessError, ET.ParseError) as error:
        print(f"ERROR: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())