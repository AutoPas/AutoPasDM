#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Runtime.h"
#include "src/distributed/MoleculeLJParticleSerializer.h"

namespace {

using Particle = mdLib::MultisiteMoleculeLJ;

constexpr std::array<double, 3> globalMin{0., 0., 0.};
constexpr std::array<double, 3> globalMax{4., 1., 1.};
constexpr double cutoff = 0.2;

auto makeConfigurator() {
  return [](auto &autoPas) {
    autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
    autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::enabled});
    autoPas.setVerletSkin(0.0);
  };
}

TEST(MultisiteParticleMigrationTest, PreservesStateWhenParticleMigratesAcrossRankBoundary) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());

  Particle original;
  original.setID(42);
  original.setTypeId(7);
  original.setR({1.75, 0.5, 0.5});
  original.setV({1., 2., 3.});
  original.setF({4., 5., 6.});
  original.setOldF({7., 8., 9.});
  original.setQuaternion({0.1, 0.2, 0.3, 0.4});
  original.setAngularVel({10., 11., 12.});
  original.setTorque({13., 14., 15.});

  particles.addParticlesFromRoot(std::vector<Particle>{original});

  particles.applyToOwnedParticles([](auto &particle) {
    if (particle.getID() == 42) {
      particle.setR({2.05, 0.5, 0.5});
    }
  });

  // All ranks participate in migration and halo exchange.
  particles.prepareInteractions();

  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), 1);
  EXPECT_EQ(particles.getLocalNumberOfOwnedParticles(), runtime.rank() == 1 ? 1 : 0);

  if (runtime.rank() == 0) {
    // The migrated particle is within the halo width of rank 0's right boundary.
    EXPECT_EQ(particles.getLocalNumberOfHaloParticles(), 1);
  }

  if (runtime.rank() == 1) {
    std::size_t found = 0;
    particles.forEachOwnedParticle([&](const Particle &particle) {
      ++found;
      EXPECT_EQ(particle.getID(), 42);
      EXPECT_EQ(particle.getTypeId(), 7);
      EXPECT_EQ(particle.getR(), (std::array<double, 3>{2.05, 0.5, 0.5}));
      EXPECT_EQ(particle.getV(), original.getV());
      EXPECT_EQ(particle.getF(), original.getF());
      EXPECT_EQ(particle.getOldF(), original.getOldF());
      EXPECT_EQ(particle.getQuaternion(), original.getQuaternion());
      EXPECT_EQ(particle.getAngularVel(), original.getAngularVel());
      EXPECT_EQ(particle.getTorque(), original.getTorque());
    });
    EXPECT_EQ(found, 1);
  }

  particles.finalize();
}

}  // namespace