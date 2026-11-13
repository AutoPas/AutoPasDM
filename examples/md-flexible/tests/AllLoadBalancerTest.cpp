#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "TestForceFunctor.h"
#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Particle.h"
#include "autopasdm/Runtime.h"

namespace {

using Particle = autopasdm::Particle;

constexpr std::array<double, 3> globalMin{0., 0., 0.};
constexpr std::array<double, 3> globalMax{4., 4., 4.};
constexpr double cutoff = 0.6;
constexpr std::array<double, 10> initialXPositions{0.5, 1.0, 1.5, 1.9, 1.99, 2.01, 2.1, 2.5, 3.0, 3.5};

Particle makeParticle(unsigned long id, double x, double y = 2., double z = 2.) {
  Particle particle;
  particle.setID(id);
  particle.setR({x, y, z});
  particle.setF({0., 0., 0.});
  return particle;
}

auto makeConfigurator() {
  return [](auto &autoPas) {
    autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
    autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    // Each rank computes forces for its owned particles, including against halo particles.
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
    autoPas.setVerletSkin(0.0);
  };
}

std::vector<Particle> makeInitialParticles() {
  std::vector<Particle> particles;
  particles.reserve(initialXPositions.size());
  for (std::size_t id = 0; id < initialXPositions.size(); ++id) {
    particles.push_back(makeParticle(static_cast<unsigned long>(id), initialXPositions[id]));
  }
  return particles;
}

std::vector<std::array<double, 3>> calculateExpectedForces(const std::vector<Particle> &particles) {
  std::vector<std::array<double, 3>> expectedForces(particles.size(), {0., 0., 0.});
  for (std::size_t i = 0; i < particles.size(); ++i) {
    for (std::size_t j = i + 1; j < particles.size(); ++j) {
      std::array<double, 3> displacement{};
      double distanceSquared = 0.;
      for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        displacement[dimension] = particles[j].getR()[dimension] - particles[i].getR()[dimension];
        distanceSquared += displacement[dimension] * displacement[dimension];
      }
      if (distanceSquared >= cutoff * cutoff) {
        continue;
      }
      for (std::size_t dimension = 0; dimension < 3; ++dimension) {
        expectedForces[i][dimension] += displacement[dimension];
        expectedForces[j][dimension] -= displacement[dimension];
      }
    }
  }
  return expectedForces;
}

int ownerOfId(autopasdm::AutoPasDM<Particle> &particles, const autopasdm::Runtime &runtime, unsigned long id) {
  int locallyOwned = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) { locallyOwned |= (particle.getID() == id); });
  // Add one so rank 0 is distinguishable from the no-owner value 0.
  return static_cast<int>(particles.globalSum(locallyOwned * (runtime.rank() + 1))) - 1;
}

}  // namespace

TEST(AllLoadBalancerTest, RebalancesMigratesAndComputesAcrossUpdatedBoundary) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  particles.addParticlesFromRoot(makeInitialParticles());

  const double initialSharedBoundary = runtime.rank() == 0 ? particles.localBoxMax()[0] : particles.localBoxMin()[0];
  ASSERT_DOUBLE_EQ(initialSharedBoundary, 2.0);

  // Rank 0 reports substantially more work, so ALL/TENSOR should move the shared
  // boundary away from the initial equal-sized split.
  const double localWork = runtime.rank() == 0 ? 100. : 1.;
  particles.prepareInteractionsWithAllLoadBalancer(localWork);

  const double localSharedBoundary = runtime.rank() == 0 ? particles.localBoxMax()[0] : particles.localBoxMin()[0];
  const double sharedBoundary = particles.globalSum(localSharedBoundary) / 2.;
  EXPECT_NEAR(localSharedBoundary, sharedBoundary, 1e-12);
  EXPECT_GT(std::abs(sharedBoundary - initialSharedBoundary), 1e-6);

  // The particles nearest the old split must remain globally unique, and at least
  // one of them must change owners as the split moves.
  constexpr std::size_t initialParticleCount = initialXPositions.size();
  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), initialParticleCount);
  for (std::size_t id = 0; id < initialParticleCount; ++id) {
    int locallyOwned = 0;
    particles.forEachOwnedParticle(
        [&](const Particle &particle) { locallyOwned |= (particle.getID() == static_cast<unsigned long>(id)); });
    EXPECT_EQ(particles.globalSum(locallyOwned), 1) << "particle id " << id;
  }
  const int ownerNearLeft = ownerOfId(particles, runtime, 4);
  const int ownerNearRight = ownerOfId(particles, runtime, 5);
  EXPECT_TRUE(ownerNearLeft != 0 || ownerNearRight != 1);

  // Insert a pair on opposite sides of the new split. Their force check below
  // exercises the updated ownership, halo exchange, and interaction path.
  constexpr double pairOffset = 0.2;
  auto leftParticle = makeParticle(10, sharedBoundary - pairOffset, 3., 3.);
  auto rightParticle = makeParticle(11, sharedBoundary + pairOffset, 3., 3.);
  particles.addParticlesFromRoot(std::vector<Particle>{leftParticle, rightParticle});

  EXPECT_NE(ownerOfId(particles, runtime, 10), ownerOfId(particles, runtime, 11));
  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), initialParticleCount + 2);

  auto allParticles = makeInitialParticles();
  allParticles.push_back(leftParticle);
  allParticles.push_back(rightParticle);
  const auto expectedForces = calculateExpectedForces(allParticles);

  autopasdm::testing::TestForceFunctor<Particle> forceFunctor(cutoff);
  particles.computeInteractions(&forceFunctor);

  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = static_cast<std::size_t>(particle.getID());
    ASSERT_LT(id, expectedForces.size());
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getF()[dimension], expectedForces[id][dimension], 1e-12)
          << "force component " << dimension << " for particle " << id;
    }
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, particles.getLocalNumberOfOwnedParticles());

  particles.finalize();
}
