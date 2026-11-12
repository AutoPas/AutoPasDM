#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "autopas/options/ContainerOption.h"
#include "autopas/options/DataLayoutOption.h"
#include "autopas/options/Newton3Option.h"
#include "autopas/options/TraversalOption.h"
#include "autopasdm/AutoPasDM.h"
#include "autopasdm/Runtime.h"
#include "molecularDynamicsLibrary/LJMultisiteFunctor.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "src/distributed/MoleculeLJParticleSerializer.h"

namespace {

using Particle = mdLib::MultisiteMoleculeLJ;
using ParticlePropertiesLibraryType = ::ParticlePropertiesLibrary<double, std::size_t>;
using LJFunctor = mdLib::LJMultisiteFunctor<Particle, true, true, autopas::FunctorN3Modes::Both, false, false>;

constexpr double cutoff = 1.0;
constexpr std::array<double, 3> globalMin{0.0, 0.0, 0.0};
constexpr std::array<double, 3> globalMax{4.0, 4.0, 4.0};

auto makeConfigurator() {
  return [](auto &autoPas) {
    autoPas.setAllowedContainers({autopas::ContainerOption::directSum});
    autoPas.setAllowedTraversals({autopas::TraversalOption::ds_sequential});
    autoPas.setAllowedDataLayouts({autopas::DataLayoutOption::aos});
    autoPas.setAllowedNewton3Options({autopas::Newton3Option::disabled});
    autoPas.setVerletSkin(0.0);
  };
}

ParticlePropertiesLibraryType makeParticlePropertiesLibrary() {
  ParticlePropertiesLibraryType properties(cutoff);

  properties.addSiteType(0, 1.0);
  properties.addLJParametersToSite(0, 1.0, 0.4);

  // Use an asymmetric molecular geometry so that torques are tested as well.
  properties.addMolType(0, {0, 0, 0}, {{-0.1, -0.05, 0.0}, {0.1, -0.05, 0.0}, {0.0, 0.1, 0.0}}, {1.0, 1.0, 1.0});

  properties.calculateMixingCoefficients();
  return properties;
}

std::vector<Particle> makeParticles() {
  Particle first;
  first.setID(0);
  first.setTypeId(0);
  first.setR({1.65, 2.0, 2.0});
  first.setQuaternion({0.0, 0.0, 0.0, 1.0});
  first.setF({0.0, 0.0, 0.0});
  first.setTorque({0.0, 0.0, 0.0});
  first.setOldF({0.0, 0.0, 0.0});

  Particle second;
  second.setID(1);
  second.setTypeId(0);
  second.setR({2.35, 2.0, 2.0});
  // Rotate around the z-axis; quaternion components are ordered as x, y, z, w.
  second.setQuaternion({0.0, 0.0, std::sin(M_PI / 8.0), std::cos(M_PI / 8.0)});
  second.setF({0.0, 0.0, 0.0});
  second.setTorque({0.0, 0.0, 0.0});
  second.setOldF({0.0, 0.0, 0.0});

  return {first, second};
}

}  // namespace

TEST(MultisiteInteractionTest, ComputesSameForcesAndTorquesAcrossRankBoundary) {
  int argc = 0;
  char **argv = nullptr;
  autopasdm::Runtime runtime(argc, argv);
  ASSERT_EQ(runtime.size(), 2);

  auto properties = makeParticlePropertiesLibrary();
  const auto inputParticles = makeParticles();

  // Compute a local reference with both molecules in the same container.
  autopas::AutoPas<Particle> reference;
  reference.setBoxMin(globalMin);
  reference.setBoxMax(globalMax);
  reference.setCutoff(cutoff);
  makeConfigurator()(reference);
  reference.init();

  for (const auto &particle : inputParticles) {
    reference.addParticle(particle);
  }

  LJFunctor referenceFunctor(cutoff, properties);
  reference.computeInteractions(&referenceFunctor);

  std::array<std::array<double, 3>, 2> expectedForces{};
  std::array<std::array<double, 3>, 2> expectedTorques{};

  for (auto iter = reference.begin(); iter.isValid(); ++iter) {
    const auto id = iter->getID();
    ASSERT_LT(id, 2u);
    expectedForces[id] = iter->getF();
    expectedTorques[id] = iter->getTorque();
  }

  // Distribute the molecules across the x-rank boundary at x = 2.
  autopasdm::AutoPasDM<Particle> distributed(runtime, globalMin, globalMax, cutoff, makeConfigurator());
  distributed.addParticlesFromRoot(inputParticles);
  distributed.computeInteractions(new LJFunctor(cutoff, properties));

  EXPECT_EQ(distributed.getGlobalNumberOfOwnedParticles(), 2);

  // Each rank's local container should contain the other molecule as a halo particle.
  EXPECT_GE(distributed.getLocalNumberOfHaloParticles(), 1);

  std::size_t checkedParticles = 0;
  distributed.forEachOwnedParticle([&](const Particle &particle) {
    const auto id = particle.getID();
    ASSERT_LT(id, 2u);

    for (std::size_t dimension = 0; dimension < 3; ++dimension) {
      EXPECT_NEAR(particle.getF()[dimension], expectedForces[id][dimension], 1e-10)
          << "force component " << dimension << " for particle " << id;
      EXPECT_NEAR(particle.getTorque()[dimension], expectedTorques[id][dimension], 1e-10)
          << "torque component " << dimension << " for particle " << id;
    }
    ++checkedParticles;
  });

  EXPECT_EQ(checkedParticles, 1);
  distributed.finalize();
}