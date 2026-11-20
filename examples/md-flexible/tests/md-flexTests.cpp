/**
 * @file md-flexTests.cpp
 * @author F. Gratl
 * @date 06.11.19
 */

#include <gtest/gtest.h>

#include <cstring>

#include "autopas/utils/WrapMPI.h"

int main(int argc, char **argv) {
  // GoogleTest discovery only needs the registered tests and does not need MPI. Avoid initializing
  // MPI here because some MPI runtimes leave allocations that LeakSanitizer reports at process exit.
  bool isTestDiscovery = false;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--gtest_list_tests") == 0) {
      isTestDiscovery = true;
      break;
    }
  }

  testing::InitGoogleTest(&argc, argv);
  // set the gtest death test style to threadsafe
  testing::FLAGS_gtest_death_test_style = "threadsafe";

  if (not isTestDiscovery) {
    autopas::AutoPas_MPI_Init(&argc, &argv);
  }
  int output = RUN_ALL_TESTS();
  if (not isTestDiscovery) {
    autopas::AutoPas_MPI_Finalize();
  }
  return output;
}
