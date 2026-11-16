if (NOT AUTOPASDM_BUILD_TESTS)
    message(STATUS "AUTOPASDM_BUILD_TESTS is disabled. Not building md-flexible tests!")
    return()
endif ()
message(STATUS "mdflexTests - adding tests.")

file(
    GLOB_RECURSE
    MDFlexTest_SRC
    "../src/*.cpp"
    "../src/*.h"
    "./*.cpp"
    "./*.h"
)

# remove main to avoid duplicate two main functions
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*main.cpp")
# remove files which are not needed for tests and expensive to compile
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*templateInstantiations/computeInteractions.*cpp")
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*Simulation.((cpp)|h)")
# Particle serialization is now implemented and tested by AutoPasDM; this legacy md-flexible test has no implementation.
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*ParticleSerializationToolsTest.((cpp)|h)")
# This regression is built as a separate MPI test below, not as part of the single-process gtest executable.
list(FILTER MDFlexTest_SRC EXCLUDE REGEX ".*AllLoadBalancerTest.cpp$")

add_executable(mdFlexTests ${MDFlexTest_SRC})

target_compile_definitions(
    mdFlexTests PRIVATE
    YAMLDIRECTORY=\"${CMAKE_CURRENT_SOURCE_DIR}/yamlTestFiles/\"
    $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_AUTOVEC>
    $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_AVX}>:MD_FLEXIBLE_FUNCTOR_AVX>
    $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_SVE}>:MD_FLEXIBLE_FUNCTOR_SVE>
    $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_HWY}>:MD_FLEXIBLE_FUNCTOR_HWY>
    $<$<BOOL:${MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC}>:MD_FLEXIBLE_FUNCTOR_ATM_AUTOVEC>
    $<$<BOOL:${MD_FLEXIBLE_CALC_GLOBALS}>:MD_FLEXIBLE_CALC_GLOBALS>
)

include(autopas_all)

target_include_directories(
    mdFlexTests
    PUBLIC ${AUTOPAS_SOURCE_DIR}/tests ${AUTOPAS_SOURCE_DIR}/tests/testAutopas ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(
  mdFlexTests
  PUBLIC
    autopasdm
    autopas
    molecularDynamicsLibrary
    GTest::gmock
    yaml-cpp
    ${ALL_LIB} 
    MPI::MPI_CXX
)

if(MD_FLEXIBLE_ENABLE_ALLLBL)
    target_compile_definitions(mdFlexTests PRIVATE AUTOPASDM_ENABLE_ALL_LOAD_BALANCER MD_FLEXIBLE_ENABLE_ALLLBL)
endif ()

if(MD_FLEXIBLE_ENABLE_ALLLBL)
    add_executable(
        autopasdm_all_load_balancer_tests
        AllLoadBalancerTest.cpp
        "${CMAKE_CURRENT_SOURCE_DIR}/../../../tests/MPITestMain.cpp"
    )

    target_include_directories(
        autopasdm_all_load_balancer_tests
        PRIVATE
            "${PROJECT_SOURCE_DIR}/tests"
    )

    target_link_libraries(
        autopasdm_all_load_balancer_tests
        PRIVATE
            autopasdm
            GTest::gtest
            ${ALL_LIB}
            MPI::MPI_CXX
    )

    target_compile_definitions(autopasdm_all_load_balancer_tests PRIVATE AUTOPASDM_ENABLE_ALL_LOAD_BALANCER)
    target_compile_features(autopasdm_all_load_balancer_tests PRIVATE cxx_std_20)

    add_test(
        NAME AutoPasDM/AllLoadBalancerTest.RebalancesMigratesAndComputesAcrossUpdatedBoundary
        COMMAND
            ${MPIEXEC_EXECUTABLE}
            ${MPIEXEC_NUMPROC_FLAG} 2
            ${MPIEXEC_PREFLAGS}
            $<TARGET_FILE:autopasdm_all_load_balancer_tests>
            "--gtest_filter=AllLoadBalancerTest.RebalancesMigratesAndComputesAcrossUpdatedBoundary"
            ${MPIEXEC_POSTFLAGS}
    )
    set_tests_properties(
        AutoPasDM/AllLoadBalancerTest.RebalancesMigratesAndComputesAcrossUpdatedBoundary
        PROPERTIES
            PROCESSORS 2
            LABELS "md-flexible;integration;load-balancing"
    )
endif ()

# this cmake module was only introduced in 3.10
include(GoogleTest)
# more robust, queries the compiled executable
gtest_discover_tests(
        mdFlexTests TEST_PREFIX "AutoPasDM/mdFlex/"
        # Use PRE_TEST discovery mode to avoid running MPI executable during build
        # This prevents SIGPIPE errors in CI when sanitizers are enabled
        DISCOVERY_MODE PRE_TEST
        PROPERTIES
            ENVIRONMENT "${LSAN_OPTIONS_STR}"
)
