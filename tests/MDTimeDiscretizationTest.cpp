#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Runtime.h"
#include "src/TimeDiscretization.h"
#include "src/distributed/MoleculeLJParticleSerializer.h"

namespace {

using Particle = ParticleType;
using Properties = ParticlePropertiesLibraryType;

constexpr double cutoff = 0.5;
constexpr std::array<double, 3> globalMin{0., 0., 0.};
constexpr std::array<double, 3> globalMax{4., 4., 4.};

void configureAutoPas(auto &autoPas) {
  autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
  autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
  autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
  autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
  autoPas.setVerletSkin(0.2);
  autoPas.setVerletRebuildFrequency(1);
}

Properties makeTranslationProperties() {
  Properties properties(cutoff);
  properties.addSiteType(0, 1.0);
  properties.addLJParametersToSite(0, 1.0, 1.0);
  properties.addSiteType(1, 2.5);
  properties.addLJParametersToSite(1, 1.0, 1.0);
#if MD_FLEXIBLE_MODE == MULTISITE
  properties.addMolType(0, {0}, {{0., 0., 0.}}, {1., 1., 1.});
  properties.addMolType(1, {1}, {{0., 0., 0.}}, {1., 1., 1.});
#endif
  properties.calculateMixingCoefficients();
  return properties;
}

Particle makeParticle(unsigned long id, std::size_t typeId, const std::array<double, 3> &position,
                      const std::array<double, 3> &force) {
  Particle particle;
  particle.setID(id);
  particle.setTypeId(typeId);
  particle.setR(position);
  particle.setV({0., 0., 0.});
  particle.setF(force);
#if MD_FLEXIBLE_MODE == MULTISITE
  particle.setQuaternion({1., 0., 0., 0.});
  particle.setTorque({0., 0., 0.});
  particle.setAngularVel({0., 0., 0.});
#endif
  return particle;
}

auto makeConfigurator() {
  return [](auto &autoPas) { configureAutoPas(autoPas); };
}

}  // namespace

TEST(MDTimeDiscretizationTest, AppliesMassScaledGravityAndGlobalForceToDistributedParticles) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  auto properties = makeTranslationProperties();
  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());

  constexpr std::array<double, 3> globalForce{1., 0., 0.};
  constexpr std::array<double, 3> gravity{0., 0., -9.81};
  constexpr double deltaT = 0.1;
  const std::vector<Particle> inputParticles{
      makeParticle(0, 0, {0.75, 1., 1.}, {1., 0., -9.81}),
      makeParticle(1, 1, {3.25, 1., 1.}, {1., 0., -24.525}),
  };
  particles.addParticlesFromRoot(inputParticles);

  TimeDiscretization::calculatePositionsAndResetForces(particles, properties, deltaT, globalForce, gravity, false);

  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), 2);
  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_LT(id, 2ul);
    const auto mass = id == 0 ? 1.0 : 2.5;
    const auto initialX = id == 0 ? 0.75 : 3.25;

    EXPECT_NEAR(particle.getR()[0], initialX + deltaT * deltaT / (2.0 * mass), 1e-12);
    EXPECT_NEAR(particle.getR()[1], 1.0, 1e-12);
    EXPECT_NEAR(particle.getR()[2], 1.0 + 0.5 * gravity[2] * deltaT * deltaT, 1e-12);
    EXPECT_NEAR(particle.getF()[0], globalForce[0], 1e-12);
    EXPECT_NEAR(particle.getF()[1], globalForce[1], 1e-12);
    EXPECT_NEAR(particle.getF()[2], gravity[2] * mass, 1e-12);
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, 1);
  particles.finalize();
}

TEST(MDTimeDiscretizationTest, CalculatesGravityTorqueOnDistributedMultisiteMolecules) {
#if MD_FLEXIBLE_MODE == MULTISITE
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  Properties properties(cutoff);
  properties.addSiteType(0, 1.0);
  properties.addLJParametersToSite(0, 1.0, 1.0);
  properties.addSiteType(1, 3.0);
  properties.addLJParametersToSite(1, 1.0, 1.0);
  properties.addMolType(0, {0, 1}, {{-1., 0., 0.}, {1., 0., 0.}}, {1., 1., 1.});
  properties.calculateMixingCoefficients();

  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  std::vector<Particle> inputParticles;
  for (unsigned long id = 0; id < 2; ++id) {
    inputParticles.push_back(makeParticle(id, 0, {id == 0 ? 0.75 : 3.25, 1., 1.}, {0., 0., 0.}));
  }
  particles.addParticlesFromRoot(inputParticles);

  constexpr std::array<double, 3> globalForce{0., 0., 0.};
  constexpr std::array<double, 3> gravity{0., 0., -10.};
  TimeDiscretization::calculateQuaternionsAndResetTorques(particles, properties, 0.01, globalForce, gravity);

  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), 2);
  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    EXPECT_NEAR(particle.getTorque()[0], 0., 1e-12);
    EXPECT_NEAR(particle.getTorque()[1], 20., 1e-12);
    EXPECT_NEAR(particle.getTorque()[2], 0., 1e-12);
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, 1);
  particles.finalize();
#else
  GTEST_SKIP() << "Rotational gravity torque is only available for multisite particles.";
#endif
}
TEST(MDTimeDiscretizationTest, CalculatesVelocitiesForDistributedParticles) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  auto properties = makeTranslationProperties();
  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());

  constexpr double deltaT = 0.2;
  const std::vector<Particle> inputParticles{
      makeParticle(0, 0, {0.75, 1., 1.}, {2., 4., 6.}),
      makeParticle(1, 1, {3.25, 1., 1.}, {-1., 3., 5.}),
  };
  auto particlesWithOldForces = inputParticles;
  particlesWithOldForces[0].setOldF({4., 2., 0.});
  particlesWithOldForces[1].setOldF({3., -1., 1.});
  particlesWithOldForces[0].setV({0.5, -1., 2.});
  particlesWithOldForces[1].setV({-0.5, 0.25, 1.});
  particles.addParticlesFromRoot(particlesWithOldForces);

  TimeDiscretization::calculateVelocities(particles, properties, deltaT);

  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), 2);
  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_LT(id, 2ul);
    const auto mass = id == 0 ? 1.0 : 2.5;
    const auto initialVelocity = id == 0 ? std::array<double, 3>{0.5, -1., 2.} : std::array<double, 3>{-0.5, 0.25, 1.};
    const auto force = id == 0 ? std::array<double, 3>{2., 4., 6.} : std::array<double, 3>{-1., 3., 5.};
    const auto oldForce = id == 0 ? std::array<double, 3>{4., 2., 0.} : std::array<double, 3>{3., -1., 1.};

    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      const auto expectedVelocity =
          initialVelocity[dimension] + (force[dimension] + oldForce[dimension]) * deltaT / (2.0 * mass);
      EXPECT_NEAR(particle.getV()[dimension], expectedVelocity, 1e-12)
          << "velocity component " << dimension << " for particle " << id;
    }
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, 1);
  particles.finalize();
}

TEST(MDTimeDiscretizationTest, CalculatesAngularVelocitiesForDistributedMultisiteMolecules) {
#if MD_FLEXIBLE_MODE == MULTISITE
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  Properties properties(cutoff);
  properties.addSiteType(0, 1.0);
  properties.addLJParametersToSite(0, 1.0, 1.0);
  properties.addSiteType(1, 1.0);
  properties.addLJParametersToSite(1, 1.0, 1.0);
  properties.addMolType(0, {0, 1}, {{-0.1, 0., 0.}, {0.1, 0., 0.}}, {1., 2., 4.});
  properties.calculateMixingCoefficients();

  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  auto first = makeParticle(0, 0, {0.75, 1., 1.}, {0., 0., 0.});
  first.setAngularVel({0.5, -1., 0.25});
  first.setTorque({2., 4., 8.});
  auto second = makeParticle(1, 0, {3.25, 1., 1.}, {0., 0., 0.});
  second.setAngularVel({-0.5, 0.25, 1.});
  second.setTorque({-2., 6., 4.});
  particles.addParticlesFromRoot(std::vector<Particle>{first, second});

  constexpr double deltaT = 0.2;
  TimeDiscretization::calculateAngularVelocities(particles, properties, deltaT);

  EXPECT_EQ(particles.getGlobalNumberOfOwnedParticles(), 2);
  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_LT(id, 2ul);
    const auto initialAngularVelocity =
        id == 0 ? std::array<double, 3>{0.5, -1., 0.25} : std::array<double, 3>{-0.5, 0.25, 1.};
    const auto torque = id == 0 ? std::array<double, 3>{2., 4., 8.} : std::array<double, 3>{-2., 6., 4.};
    constexpr std::array<double, 3> momentOfInertia{1., 2., 4.};

    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      const auto expectedAngularVelocity =
          initialAngularVelocity[dimension] + torque[dimension] / momentOfInertia[dimension] * deltaT / 2.0;
      EXPECT_NEAR(particle.getAngularVel()[dimension], expectedAngularVelocity, 1e-12)
          << "angular velocity component " << dimension << " for particle " << id;
    }
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, 1);
  particles.finalize();
#else
  GTEST_SKIP() << "Angular velocity integration is only available for multisite particles.";
#endif
}
