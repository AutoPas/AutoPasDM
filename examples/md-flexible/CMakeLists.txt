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
    if(MD_FLEXIBLE_FUNCTOR_AVX AND MD_FLEXIBLE_MODE STREQUAL "SINGLESITE")
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

    # Exercise the complete YAML -> Simulation -> AutoPasDM ALL/TENSOR path on two ranks.
    if(MD_FLEXIBLE_ENABLE_ALLLBL AND MD_FLEXIBLE_USE_MPI)
        add_test(
            NAME md-flexible.test-all-load-balancer
            COMMAND
                ${MPIEXEC_EXECUTABLE}
                ${MPIEXEC_NUMPROC_FLAG} 2
                ${MPIEXEC_PREFLAGS}
                $<TARGET_FILE:md-flexible>
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/input/autopasdm_all_load_balancer.yaml"
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
                --yaml-filename "${CMAKE_CURRENT_SOURCE_DIR}/input/autopasdm_mpi_tuning.yaml"
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
    endif()

    if(MD_FLEXIBLE_FUNCTOR_AVX AND MD_FLEXIBLE_MODE STREQUAL "SINGLESITE")
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
    endif()
endif()
