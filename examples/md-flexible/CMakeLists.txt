# md-flexible adapted as a bundled AutoPasDM example.
# AutoPas and AutoPasDM are configured by the repository root.

list(APPEND CMAKE_MODULE_PATH "${AUTOPAS_SOURCE_DIR}/cmake/modules")

include(autopas_yaml-cpp)
include(autopas_all)

option(MD_FLEXIBLE_USE_MPI "Enable MPI support for the AutoPasDM md-flexible example." ON)
option(MD_FLEXIBLE_FUNCTOR_AUTOVEC "Compile AutoVec Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_AVX "If instruction set is available, compile AVX Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_SVE "If instruction set is available, compile SVE Functor." OFF)
option(MD_FLEXIBLE_FUNCTOR_HWY "Compile Highway Functor." ON)
option(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC "Compile auto vectorized Axilrod-Teller-Muto Functor." OFF)
option(MD_FLEXIBLE_CALC_GLOBALS "Compile md-flexible functors with global calculations." OFF)

file(
    GLOB_RECURSE MDFlex_SRC
    "src/*.cpp"
    "src/*.h"
)

add_executable(md-flexible ${MDFlex_SRC})
target_include_directories(md-flexible PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}")

# md-flexible exposes the optional ALL dependency through its configuration, while the
# corresponding implementation is part of the application-independent AutoPasDM library.
if(MD_FLEXIBLE_ENABLE_ALLLBL)
    target_compile_definitions(
        md-flexible
        PRIVATE
            AUTOPASDM_ENABLE_ALL_LOAD_BALANCER
            MD_FLEXIBLE_ENABLE_ALLLBL
    )
endif()

target_link_libraries(
    md-flexible
    PUBLIC
        autopas
        autopasdm
        molecularDynamicsLibrary
        yaml-cpp
        ${ALL_LIB}
        $<$<BOOL:${MD_FLEXIBLE_USE_MPI}>:MPI::MPI_CXX>
)

target_compile_definitions(
    md-flexible
    PUBLIC
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AVX}>:MD_FLEXIBLE_FUNCTOR_AVX>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_SVE}>:MD_FLEXIBLE_FUNCTOR_SVE>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_HWY}>:MD_FLEXIBLE_FUNCTOR_HWY>
        $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC>
        $<$<BOOL:${MD_FLEXIBLE_CALC_GLOBALS}>:MD_FLEXIBLE_CALC_GLOBALS>
)

# Keep upstream's build-time notice for functors that are unavailable in multisite mode.
string(TOLOWER "${MD_FLEXIBLE_MODE}" MD_FLEXIBLE_MODE_lower)
if(MD_FLEXIBLE_MODE_lower MATCHES "multisite")
    if(MD_FLEXIBLE_FUNCTOR_AVX)
        message(WARNING "AVX Lennard-Jones functor has not been implemented for Multi-Site Molecules")
    endif()
    if(MD_FLEXIBLE_FUNCTOR_SVE)
        message(WARNING "SVE Lennard-Jones functor has not been implemented for Multi-Site Molecules")
    endif()
    if(MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC)
        message(WARNING "Axilrod-Teller-Muto functor has not been implemented for Multi-Site Molecules")
    endif()
    if(MD_FLEXIBLE_FUNCTOR_HWY)
        message(WARNING "Highway functor has not been implemented for Multi-Site Molecules")
    endif()
endif()

# Copy input and helper files next to the executable as in the original
# md-flexible build.
file(
    GLOB_RECURSE MDFLEX_AUX_FILES
    RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}"
    "*.sh"
    "*.gp"
    "*.py"
    "*.yaml"
    "*.rule"
    "*.frule"
)

foreach(aux_file ${MDFLEX_AUX_FILES})
    configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/${aux_file}"
        "${CMAKE_CURRENT_BINARY_DIR}/${aux_file}"
        COPYONLY
    )
endforeach()

# Small integration example: md-flexible's MoleculeLJ + md-flexible's ordinary
# AutoPas LJFunctor running through AutoPasDM. No functor adapter.
if(MD_FLEXIBLE_MODE STREQUAL "SINGLESITE")
    add_executable(autopasdm-lj-example examples/AutoPasDMLJExample.cpp)
    target_include_directories(autopasdm-lj-example PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(
        autopasdm-lj-example
        PRIVATE
            autopasdm
            molecularDynamicsLibrary
    )
endif()

if(AUTOPASDM_BUILD_TESTS)
    add_subdirectory(tests)

    set(MDFLEX_TEST_FUNCTOR "Lennard-Jones")

    # The AVX functor needs both an AVX-enabled md-flexible build and an AutoPas
    # vectorization configuration that can run the corresponding traversals.
    set(HAS_AVX false)
    if(MD_FLEXIBLE_FUNCTOR_AVX AND MD_FLEXIBLE_MODE STREQUAL "SINGLESITE")
        if(NOT AUTOPAS_USE_AUTOVEC)
            message(STATUS "No AutoPas vectorization used; not adding md-flexible AVX tests.")
        elseif(AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SSE")
            message(STATUS "Only SSE specified; not adding md-flexible AVX tests.")
        elseif(AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "SVE")
            message(STATUS "SVE vectorization level; not adding md-flexible AVX tests.")
        elseif(AUTOPAS_VECTOR_INSTRUCTIONS MATCHES "NATIVE")
            message(STATUS "Native vectorization level; detecting AVX capability.")
            try_compile(
                HAS_AVX
                "${CMAKE_BINARY_DIR}"
                "${AUTOPAS_SOURCE_DIR}/cmake/tests/has_avx_test.cpp"
                COMPILE_DEFINITIONS "-march=native"
                OUTPUT_VARIABLE HAS_AVX_ERROR
            )
        else()
            message(STATUS "Adding md-flexible AVX tests for the selected vectorization level.")
            set(HAS_AVX true)
        endif()
    endif()

    if(HAS_AVX)
        set(MDFLEX_TEST_FUNCTOR "Lennard-Jones-AVX")
    endif()

    # Small single-rank runs exercise the same command-line paths as upstream md-flexible.
    add_test(
        NAME md-flexible.test-static
        COMMAND md-flexible
            --container linked
            --cutoff 1.
            --distribution-mean 5.0
            --distribution-stddeviation 2.0
            --data-layout soa
            --functor "Lennard-Jones"
            --iterations 10
            --particle-generator gauss
            --particles-total 10
            --traversal c08,sliced
            --verlet-rebuild-frequency 5
            --verlet-skin-radius 0
            --boundary-type none
            --deltaT 0.
            --no-end-config
            --no-progress-bar
        CONFIGURATIONS checkExamples
    )
    set_tests_properties(md-flexible.test-static PROPERTIES LABELS "md-flexible;integration")
    if(MD_FLEXIBLE_USE_MPI)
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-static)
    endif()

    add_test(
        NAME md-flexible.test-sim
        COMMAND md-flexible
            --container all
            --traversal all
            --cutoff 1.5
            --functor "${MDFLEX_TEST_FUNCTOR}"
            --vectorization-pattern 1xVectorLength
            --tuning-phases 1
            --particle-generator grid
            --particles-per-dimension 10
            --particle-spacing 1.1225
            --verlet-rebuild-frequency 4
            --verlet-skin-radius 0.2
            --boundary-type periodic
            --deltaT 0.005
            --no-end-config
            --no-progress-bar
        CONFIGURATIONS checkExamples
    )
    set_tests_properties(md-flexible.test-sim PROPERTIES LABELS "md-flexible;integration")
    if(MD_FLEXIBLE_USE_MPI)
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-sim)
    endif()

    # Keep the upstream lightweight performance smoke test when its Highway
    # single-site functor is available.
    if(MD_FLEXIBLE_MODE STREQUAL "SINGLESITE" AND MD_FLEXIBLE_FUNCTOR_HWY)
        add_test(
            NAME md-flexMeasurePerf
            COMMAND
                bash "${CMAKE_CURRENT_SOURCE_DIR}/scripts/measurePerf.sh"
                "$<TARGET_FILE:md-flexible>"
                -silent
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexMeasurePerf
            PROPERTIES
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                LABELS "md-flexible;performance"
                ENVIRONMENT "${LSAN_OPTIONS_STR}"
        )
        if(MD_FLEXIBLE_USE_MPI)
            autopasdm_disable_lsan_for_mpi_test(md-flexMeasurePerf)
        endif()
    endif()

    # Exercise the complete YAML -> Simulation -> AutoPasDM ALL/TENSOR path on two ranks.
    if(MD_FLEXIBLE_ENABLE_ALLLBL AND MD_FLEXIBLE_USE_MPI)
        add_test(
            NAME md-flexible.test-all-load-balancer
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/tests/yamlTestFiles/autopasdm_all_load_balancer.yaml"
                ${MPIEXEC_POSTFLAGS}
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexible.test-all-load-balancer
            PROPERTIES
                PROCESSORS 2
                LABELS "md-flexible;integration;load-balancing"
                PASS_REGULAR_EXPRESSION "Owned: 126"
        )
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-all-load-balancer)
    endif()

    # Exercise AutoPas' MPI divide-and-conquer tuner through the AutoPasDM
    # md-flexible integration on two ranks.
    if(MD_FLEXIBLE_USE_MPI)
        add_test(
            NAME md-flexible.test-mpi-tuning
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/tests/yamlTestFiles/autopasdm_mpi_tuning.yaml"
                ${MPIEXEC_POSTFLAGS}
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexible.test-mpi-tuning
            PROPERTIES
                PROCESSORS 2
                LABELS "md-flexible;integration;tuning"
                PASS_REGULAR_EXPRESSION "Owned: 64"
        )
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-mpi-tuning)
    endif()

    # Exercise the full CLI checkpoint path: one run writes a two-rank VTK checkpoint,
    # and a second run restores it while using a deliberately different input particle count.
    if(MD_FLEXIBLE_USE_MPI)
        set(MDFLEX_CHECKPOINT_TEST_OUTPUT_FOLDER "autopasdm-cli-checkpoint-output")
        set(MDFLEX_CHECKPOINT_TEST_FILE
            "${CMAKE_CURRENT_BINARY_DIR}/${MDFLEX_CHECKPOINT_TEST_OUTPUT_FOLDER}/autopasdm-cli-checkpoint/autopasdm-cli-checkpoint_Particles_1.pvtu"
        )

        add_test(
            NAME md-flexible.test-checkpoint-cli-write
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/tests/yamlTestFiles/autopasdm_cli_checkpoint.yaml"
                --vtk-filename autopasdm-cli-checkpoint
                --vtk-output-folder "${MDFLEX_CHECKPOINT_TEST_OUTPUT_FOLDER}"
                --vtk-write-frequency 1
                ${MPIEXEC_POSTFLAGS}
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexible.test-checkpoint-cli-write
            PROPERTIES
                PROCESSORS 2
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                FIXTURES_SETUP MDFlexCLICheckpoint
                LABELS "md-flexible;integration;checkpoint;cli"
                PASS_REGULAR_EXPRESSION "Owned: 8"
        )
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-checkpoint-cli-write)

        add_test(
            NAME md-flexible.test-checkpoint-cli-restart
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/tests/yamlTestFiles/autopasdm_cli_checkpoint_restart.yaml"
                --checkpoint "${MDFLEX_CHECKPOINT_TEST_FILE}"
                --vtk-filename autopasdm-cli-checkpoint-restart
                --vtk-output-folder "${MDFLEX_CHECKPOINT_TEST_OUTPUT_FOLDER}"
                --vtk-write-frequency 1
                ${MPIEXEC_POSTFLAGS}
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexible.test-checkpoint-cli-restart
            PROPERTIES
                PROCESSORS 2
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                FIXTURES_REQUIRED MDFlexCLICheckpoint
                LABELS "md-flexible;integration;checkpoint;cli"
                PASS_REGULAR_EXPRESSION "Owned: 8"
        )
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-checkpoint-cli-restart)
    endif()

    if(MD_FLEXIBLE_USE_MPI)
        add_test(
            NAME md-flexible.test-checkpoint-cli-restart-different-rank-count
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 1
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/tests/yamlTestFiles/autopasdm_cli_checkpoint_restart.yaml"
                --checkpoint "${MDFLEX_CHECKPOINT_TEST_FILE}"
                --vtk-filename autopasdm-cli-checkpoint-restart-one-rank
                --vtk-output-folder "${MDFLEX_CHECKPOINT_TEST_OUTPUT_FOLDER}"
                --vtk-write-frequency 1
                ${MPIEXEC_POSTFLAGS}
            CONFIGURATIONS checkExamples
        )
        set_tests_properties(
            md-flexible.test-checkpoint-cli-restart-different-rank-count
            PROPERTIES
                PROCESSORS 1
                WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
                FIXTURES_REQUIRED MDFlexCLICheckpoint
                LABELS "md-flexible;integration;checkpoint;cli"
                PASS_REGULAR_EXPRESSION "Owned: 8"
        )
        autopasdm_disable_lsan_for_mpi_test(md-flexible.test-checkpoint-cli-restart-different-rank-count)
    endif()

    if(HAS_AVX)
        add_test(
            NAME md-flexible-avx.test-unaligned
            COMMAND md-flexible
                --no-end-config
                --functor "Lennard-Jones-AVX"
                --vectorization-pattern 1xVectorLength
                --deltaT 0
                --particle-generator uniform
                --log-level debug
                --traversal lc_c04_combined_SoA
                --particles-total 71
                --no-progress-bar
                CONFIGURATIONS checkExamples
        )
        set_tests_properties(md-flexible-avx.test-unaligned PROPERTIES LABELS "md-flexible;integration")
        if(MD_FLEXIBLE_USE_MPI)
            autopasdm_disable_lsan_for_mpi_test(md-flexible-avx.test-unaligned)
        endif()
    endif()
endif()

# In a larger AutoPas build this target collects example tests;
# standalone AutoPasDM builds do not necessarily define it.
if(TARGET checkExamples)
    add_dependencies(checkExamples md-flexible)
endif()
