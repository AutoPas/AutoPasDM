# AutoPasDM

AutoPasDM is a distributed particle-container layer around AutoPas.
It keeps distributed ownership, migration, halo exchange, adaptive load balancing,
and communication out of the simulator while reusing AutoPas for node-level particle
storage, interaction traversal, and auto-tuning.

The current implementation uses a Cartesian process grid in x, y, and z. It starts
from a regular decomposition and can adapt the subdomain boundaries at runtime with
inverted-pressure load balancing while keeping the process-grid topology fixed.
Boundary conditions are selected per dimension. `periodic`, `reflective`, and `none`
are represented in the distributed topology. Reflective boundaries are non-periodic
for migration and halo exchange; the bundled md-flexible example applies reflective
mirror-wall forces at application level for both single-site and multi-site particles.

The simulator-facing API deliberately exposes operations on the distributed particle
system instead of the node-level AutoPas container. Existing AutoPas particle types
and interaction functors can therefore be reused while AutoPasDM owns the
distributed execution mechanisms.

## Repository layout

```text
AutoPasDM/
├── include/autopasdm/   # public headers and template code
├── src/                           # compiled non-template implementation
├── examples/
│   ├── minimal_mpi.cpp            # low-level runtime / integration example
│   ├── lj/                        # minimal LJ simulator using AutoPasDM
│   └── md-flexible/               # md-flexible integration and regression inputs
├── docs/
└── tests/
```

The expected development workspace is:

```text
<workspace>/
├── AutoPas/
└── AutoPasDM/
```

## Build

From the AutoPasDM repository:

```bash
cmake -S . -B build \
  -DAUTOPAS_SOURCE_DIR=../AutoPas \
  -DMD_FLEXIBLE_USE_MPI=ON \
  -DMD_FLEXIBLE_FUNCTOR_AUTOVEC=ON

cmake --build build -j
```

The adapted md-flexible executable is then available at:

```text
build/examples/md-flexible/md-flexible
```

The minimal LJ example is available at:

```text
build/examples/lj/lj
```

and a small md-flexible-particle integration example at:

```text
build/examples/md-flexible/autopasdm-lj-example
```

### Quick MPI test

After building with MPI support, a small md-flexible example can be started directly
from the AutoPasDM repository root with four MPI ranks:

```bash
mpirun -np 4 ./build/examples/md-flexible/md-flexible \
  --yaml-filename ./examples/md-flexible/input/autopasdm_performance_baseline.yaml
```

This runs a short periodic Lennard-Jones simulation and is a convenient check that
the AutoPasDM md-flexible integration starts correctly with multiple ranks.

## Current feature set

- Cartesian domain decomposition in one, two, or three dimensions.
- Per-dimension `periodic`, `reflective`, and `none` boundary topology.
- Root-owned and already-distributed initial particle insertion.
- Staged x -> y -> z particle migration through direct face neighbors.
- Staged x -> y -> z halo exchange with propagation to edge and corner neighbors.
- Split-phase non-blocking MPI payload transport for timestep migration and halo exchange.
- Adaptive regular-grid decomposition with inverted-pressure load balancing, plus an
  optional ALL/TENSOR load-balancing path when compiled with ALL support.
- Adaptive ownership lookup using the complete process-plane split positions.
- A compact operation-based public API for particle updates, reductions, diagnostics,
  distributed synchronization, and AutoPas interaction functors.
- Integration with the bundled md-flexible example in single-site and multi-site modes,
  including `None`, `InvertedPressure`, and optional `ALL` load-balancer settings and
  `load-balancing-interval`.

The ALL/TENSOR path is optional and requires building with
`-DMD_FLEXIBLE_ENABLE_ALLLBL=ON` and an available ALL dependency. It is implemented,
and when enabled is covered by a two-rank AutoPasDM balancing test and a short
two-rank md-flexible application test.
Communication/computation overlap is not used yet, and MPI is the only communication
backend. See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for details and design constraints.

## Tests

AutoPasDM tests use GoogleTest and are enabled by default for a standalone
checkout. CTest launches MPI tests with their required rank counts. The suite covers,
among other things:

- Cartesian and adaptive domain geometry,
- initial particle redistribution,
- staged migration and halo exchange in 3D,
- periodic, reflective, and open (`none`) topology,
- split-phase non-blocking neighbor communication,
- the public AutoPasDM interaction and operation API,
- inverted-pressure boundary movement and minimum-width enforcement,
- repeated adaptive load-balancing steps,
- adaptive `targetRank()` ownership, and
- an eight-rank 3D end-to-end load-balancing, migration, halo, and force check.

All tests can be run with:

```bash
cmake --build build -j
ctest --test-dir build --output-on-failure
```

### md-flexible 3D force regression

A small end-to-end regression input is available at
`examples/md-flexible/tests/yamlTestFiles/autopasdm_force_compare_3d.yaml`. It exercises internal
face, edge, and corner interactions as well as periodic boundaries in all three
dimensions. The force output can be compared across AutoPasDM decompositions as shown
below. The same input is also used in CI to compare original md-flexible and AutoPasDM
with eight MPI ranks.

From `build/examples/md-flexible`:

```bash
rm -rf force_compare_3d/np1 force_compare_3d/np8
mkdir -p force_compare_3d/np1 force_compare_3d/np8

./md-flexible \
  --yaml-filename ../../../examples/md-flexible/tests/yamlTestFiles/autopasdm_force_compare_3d.yaml \
  --vtk-output-folder force_compare_3d/np1

mpirun -np 8 ./md-flexible \
  --yaml-filename ../../../examples/md-flexible/tests/yamlTestFiles/autopasdm_force_compare_3d.yaml \
  --vtk-output-folder force_compare_3d/np8

python3 ../../../examples/md-flexible/tests/scripts/compare_vtk_forces.py \
  force_compare_3d/np1 \
  force_compare_3d/np8
```

The decomposition comparison should report the same owned-particle IDs and forces
within the configured tolerance. For an original-AutoPas-vs-AutoPasDM comparison, run
from the repository root after building both MPI-enabled md-flexible executables:

```bash
python3 examples/md-flexible/tests/scripts/run_md_flexible_comparison.py \
  --reference /path/to/AutoPas/build/examples/md-flexible/md-flexible \
  --distributed /path/to/AutoPasDM/build/examples/md-flexible/md-flexible \
  --config examples/md-flexible/tests/yamlTestFiles/autopasdm_force_compare_3d.yaml \
  --ranks 8
```

The comparison checks IDs, types, positions, velocities, and forces in the final
output, without requiring the same owner rank.

### md-flexible load-balancing regression

`examples/md-flexible/tests/yamlTestFiles/autopasdm_load_balancing_compare.yaml` and
`examples/md-flexible/tests/scripts/run_load_balancing_regression.py` compare the bundled AutoPasDM
integration against an MPI-enabled original md-flexible build. Both decompositions
must adapt and the final global particle state is compared independently of the exact
owner rank, because runtime-based work measurements can produce different adaptive
split positions.

For example:

```bash
python3 examples/md-flexible/tests/scripts/run_load_balancing_regression.py \
  --reference /path/to/AutoPas/build/examples/md-flexible/md-flexible \
  --distributed ./md-flexible \
  --config examples/md-flexible/tests/yamlTestFiles/autopasdm_load_balancing_compare.yaml \
  --ranks 8
```

The regression uses the general `run_md_flexible_comparison.py` runner and additionally
requires both load balancers to move at least one subdomain boundary.

### md-flexible thermostat comparison

CI compares the original one-rank thermostat run with a four-rank AutoPasDM run using
`examples/md-flexible/tests/yamlTestFiles/autopasdm_thermostat_compare.yaml`. It compares particle velocities and checks that
both results reach the expected global temperature. The runner accepts different rank
counts with `--reference-ranks` and `--distributed-ranks`.

### md-flexible dynamic simulation comparison

`examples/md-flexible/tests/yamlTestFiles/autopasdm_dynamic_compare.yaml` advances a small 3D system
for 100 steps with a non-zero time step. Three particles cross internal MPI boundaries
in x, y, and z, while another pair exercises Lennard-Jones force evaluation. CI runs
this input with eight ranks in both implementations and compares the final global
particle states by ID.

### md-flexible reflective-boundary regression

The single-site md-flexible example contains a reflective-wall regression input at
`examples/md-flexible/tests/yamlTestFiles/autopasdm_reflective_compare.yaml`. It places particles near
all six global faces and two opposite corners. With sigma = epsilon = 1 and a wall
distance of 0.5, the expected force component from each reflective wall is exactly
24. The same setup can be compared between one rank and an eight-rank `2 x 2 x 2`
decomposition.

From `build/examples/md-flexible`:

```bash
rm -rf reflective_compare/np1 reflective_compare/np8
mkdir -p reflective_compare/np1 reflective_compare/np8

./md-flexible \
  --yaml-filename ../../../examples/md-flexible/tests/yamlTestFiles/autopasdm_reflective_compare.yaml \
  --vtk-output-folder reflective_compare/np1

mpirun -np 8 ./md-flexible \
  --yaml-filename ../../../examples/md-flexible/tests/yamlTestFiles/autopasdm_reflective_compare.yaml \
  --vtk-output-folder reflective_compare/np8

python3 ../../../examples/md-flexible/tests/scripts/compare_vtk_forces.py \
  reflective_compare/np1 \
  reflective_compare/np8

python3 ../../../examples/md-flexible/tests/scripts/check_reflective_forces.py reflective_compare/np8
```

The first script checks decomposition independence. The second checks the forces
against the analytical values of the reflective regression setup.

The reflective Lennard-Jones wall force is an AutoPasDM application feature and is not
implemented by the original md-flexible reference. CI therefore checks this path within
AutoPasDM: it compares one-rank and eight-rank forces and checks the eight-rank result
against the analytical wall forces.

### md-flexible multi-site comparison

`examples/md-flexible/tests/yamlTestFiles/autopasdm_multisite_compare.yaml` contains two tri-site
molecules interacting across an internal rank corner. The short run exercises the
multi-site force calculation, torque generation, and rotational integration. CI builds
the original md-flexible reference in multi-site mode and compares the final molecular
state with an eight-rank AutoPasDM run, including quaternions, angular velocities, and
torques.

The test build can be disabled with `-DAUTOPASDM_BUILD_TESTS=OFF`.

See [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) for the current architecture,
design constraints, and remaining limitations.
