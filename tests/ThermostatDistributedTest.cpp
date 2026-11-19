#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstddef>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Runtime.h"
#include "src/Thermostat.h"
#include "src/TypeDefinitions.h"
#include "src/distributed/MoleculeLJParticleSerializer.h"

namespace {

using Particle = ParticleType;
using Properties = ParticlePropertiesLibraryType;

constexpr double cutoff = 0.5;
constexpr std::array<double, 3> globalMin{0., 0., 0.};
constexpr std::array<double, 3> globalMax{4., 4., 4.};

auto makeConfigurator() {
  return [](auto &autoPas) {
    autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
    autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
    autoPas.setVerletSkin(0.2);
  };
}

Properties makeProperties() {
  Properties properties(cutoff);
#if MD_FLEXIBLE_MODE == MULTISITE
  properties.addSiteType(0, 1.0);
#else
  properties.addSiteType(0, 2.0);
#endif
  properties.addLJParametersToSite(0, 1.0, 1.0);
#if MD_FLEXIBLE_MODE == MULTISITE
  properties.addSiteType(1, 1.0);
  properties.addLJParametersToSite(1, 1.0, 1.0);
  properties.addMolType(0, {0, 1}, {{-0.1, 0., 0.}, {0.1, 0., 0.}}, {1., 2., 4.});
#endif
  properties.calculateMixingCoefficients();
  return properties;
}

Particle makeParticle(unsigned long id, const std::array<double, 3> &velocity,
                      const std::array<double, 3> &angularVelocity) {
  Particle particle;
  particle.setID(id);
  particle.setTypeId(0);
  particle.setR({id == 0 ? 0.75 : 3.25, 1., 1.});
  particle.setV(velocity);
#if MD_FLEXIBLE_MODE == MULTISITE
  particle.setQuaternion({1., 0., 0., 0.});
  particle.setAngularVel(angularVelocity);
  particle.setTorque({0., 0., 0.});
#else
  (void)angularVelocity;
#endif
  return particle;
}

}  // namespace

TEST(ThermostatDistributedTest, CalculatesAndAppliesGlobalTemperatureAcrossRanks) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  auto properties = makeProperties();
  autopasdm::AutoPasDM<Particle> particles(runtime, globalMin, globalMax, cutoff, makeConfigurator());

  const std::array<double, 3> firstVelocity{1., 0., 0.};
  const std::array<double, 3> secondVelocity{0., 2., 0.};
#if MD_FLEXIBLE_MODE == MULTISITE
  const std::array<double, 3> firstAngularVelocity{1., 0., 0.};
  const std::array<double, 3> secondAngularVelocity{0., 1., 0.};
#else
  const std::array<double, 3> firstAngularVelocity{};
  const std::array<double, 3> secondAngularVelocity{};
#endif
  particles.addParticlesFromRoot(std::vector<Particle>{makeParticle(0, firstVelocity, firstAngularVelocity),
                                                       makeParticle(1, secondVelocity, secondAngularVelocity)});

#if MD_FLEXIBLE_MODE == MULTISITE
  constexpr double expectedInitialTemperature = 13. / 12.;
#else
  constexpr double expectedInitialTemperature = 5. / 3.;
#endif
  const auto initialTemperature = Thermostat::calcTemperature(particles, properties);
  EXPECT_NEAR(initialTemperature, expectedInitialTemperature, 1e-12);

  const auto initialTemperatureByType = Thermostat::calcTemperatureComponent(particles, properties);
  ASSERT_EQ(initialTemperatureByType.size(), 1u);
  EXPECT_NEAR(initialTemperatureByType.at(0), expectedInitialTemperature, 1e-12);

  constexpr double targetTemperature = 3.0;
  constexpr double deltaTemperature = 0.5;
  const auto expectedTemperatureAfterApply = expectedInitialTemperature + deltaTemperature;
  const auto expectedVelocityScale = std::sqrt(expectedTemperatureAfterApply / expectedInitialTemperature);
  Thermostat::apply(particles, properties, targetTemperature, deltaTemperature);

  EXPECT_NEAR(Thermostat::calcTemperature(particles, properties), expectedTemperatureAfterApply, 1e-12);
  const auto temperatureByTypeAfterApply = Thermostat::calcTemperatureComponent(particles, properties);
  ASSERT_EQ(temperatureByTypeAfterApply.size(), 1u);
  EXPECT_NEAR(temperatureByTypeAfterApply.at(0), expectedTemperatureAfterApply, 1e-12);

  std::size_t checkedParticles = 0;
  particles.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_LT(id, 2ul);
    const auto &initialVelocityForParticle = id == 0 ? firstVelocity : secondVelocity;
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getV()[dimension], initialVelocityForParticle[dimension] * expectedVelocityScale, 1e-12)
          << "velocity component " << dimension << " for particle " << id;
    }
#if MD_FLEXIBLE_MODE == MULTISITE
    const auto &initialAngularVelocityForParticle = id == 0 ? firstAngularVelocity : secondAngularVelocity;
    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getAngularVel()[dimension],
                  initialAngularVelocityForParticle[dimension] * expectedVelocityScale, 1e-12)
          << "angular velocity component " << dimension << " for particle " << id;
    }
#endif
    ++checkedParticles;
  });
  EXPECT_EQ(checkedParticles, 1);
  particles.finalize();
}
