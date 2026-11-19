#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <string>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/InteractionTypeOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Runtime.h"
#include "src/ParallelVtkWriter.h"
#include "src/configuration/MDFlexConfig.h"
#include "src/distributed/MoleculeLJParticleSerializer.h"

namespace {

using Particle = ParticleType;
using Container = autopasdm::AutoPasDM<Particle>;

constexpr double cutoff = 0.5;
constexpr std::array<double, 3> globalMin{0., 0., 0.};
constexpr std::array<double, 3> globalMax{4., 4., 4.};
constexpr char outputFolder[] = "autopasdm-checkpoint-regression";
constexpr char sessionName[] = "distributed-checkpoint";
constexpr std::size_t checkpointIteration = 17;

auto makeConfigurator() {
  return [](auto &autoPas) {
    autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
    autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
    autoPas.setAllowedInteractionTypeOptions({autopas::InteractionTypeOption::pairwise});
    autoPas.setVerletSkin(0.2);
  };
}

std::vector<Particle> makeParticles() {
  Particle first;
  first.setID(10);
  first.setTypeId(0);
  first.setR({0.75, 1.25, 2.25});
  first.setV({0.125, -0.25, 0.5});
  first.setF({1.5, 2.5, -3.5});
#if MD_FLEXIBLE_MODE == MULTISITE
  first.setQuaternion({0.5, 0.5, 0.5, 0.5});
  first.setAngularVel({0.25, -0.5, 0.75});
  first.setTorque({-1.25, 2.25, -3.25});
#endif

  Particle second;
  second.setID(11);
  second.setTypeId(0);
  second.setR({3.25, 2.75, 1.75});
  second.setV({-0.375, 0.625, -0.875});
  second.setF({-4.5, 5.5, 6.5});
#if MD_FLEXIBLE_MODE == MULTISITE
  second.setQuaternion({0.5, -0.5, 0.5, -0.5});
  second.setAngularVel({-0.125, 0.375, -0.625});
  second.setTorque({4.25, -5.25, 6.25});
#endif

  return {first, second};
}

std::string checkpointFilename() {
  return std::string(outputFolder) + "/" + sessionName + "/" + sessionName + "_Particles_00017.pvtu";
}

void expectRestoredParticles(Container &container) {
  EXPECT_EQ(container.getGlobalNumberOfOwnedParticles(), 2);
  std::size_t checkedParticles = 0;
  container.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_TRUE(id == 10 or id == 11);
    const auto expected = makeParticles()[id == 10 ? 0 : 1];

    EXPECT_EQ(particle.getTypeId(), expected.getTypeId());
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getR()[dimension], expected.getR()[dimension], 1e-5);
      EXPECT_NEAR(particle.getV()[dimension], expected.getV()[dimension], 1e-5);
      EXPECT_NEAR(particle.getF()[dimension], expected.getF()[dimension], 1e-5);
    }
#if MD_FLEXIBLE_MODE == MULTISITE
    for (std::size_t component = 0; component < 4; ++component) {
      EXPECT_NEAR(particle.getQuaternion()[component], expected.getQuaternion()[component], 1e-5);
    }
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getAngularVel()[dimension], expected.getAngularVel()[dimension], 1e-5);
      EXPECT_NEAR(particle.getTorque()[dimension], expected.getTorque()[dimension], 1e-5);
    }
#endif
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, container.numberOfRanks() == 1 ? 2u : 1u);
}

}  // namespace

TEST(DistributedCheckpointTest, WritesTwoRankVtkCheckpoint) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  Container source(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  source.addParticlesFromRoot(makeParticles());

  ParallelVtkWriter writer(sessionName, outputFolder, 5, runtime);
  writer.recordTimestep(checkpointIteration, source);

  // Synchronize after each rank has closed its piece file before restart tests read it.
  EXPECT_EQ(source.globalSum(1), 2);
  source.finalize();
}

TEST(DistributedCheckpointTest, RestoresCheckpointWithSameRankCount) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  MDFlexConfig configuration;
  configuration.checkpointfile.value = checkpointFilename();
  configuration.loadParticlesFromCheckpoint(runtime.rank(), runtime.size());
  ASSERT_EQ(configuration.particles.size(), 1);

  Container restored(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  restored.addDistributedParticles(configuration.particles);
  expectRestoredParticles(restored);
  restored.finalize();
}

TEST(DistributedCheckpointTest, RestoresAndRedistributesCheckpointWithDifferentRankCount) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 1);

  MDFlexConfig configuration;
  configuration.checkpointfile.value = checkpointFilename();
  configuration.loadParticlesFromCheckpoint(runtime.rank(), runtime.size());
  ASSERT_EQ(configuration.particles.size(), 2);

  Container restored(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  restored.addDistributedParticles(configuration.particles);
  expectRestoredParticles(restored);
  restored.finalize();
}
