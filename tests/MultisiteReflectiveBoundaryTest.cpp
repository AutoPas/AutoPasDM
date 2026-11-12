#include <gtest/gtest.h>

#include <array>
#include <cstddef>
#include <vector>

#include "autopasdm/BoundaryType.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"
#include "molecularDynamicsLibrary/ParticlePropertiesLibrary.h"
#include "src/ReflectiveBoundary.h"

namespace {

using Particle = mdLib::MultisiteMoleculeLJ;
using ParticlePropertiesLibraryType = ::ParticlePropertiesLibrary<double, std::size_t>;

class MultisiteTestContainer {
 public:
  std::array<double, 3> globalMin{0.0, 0.0, 0.0};
  std::array<double, 3> globalMax{4.0, 4.0, 4.0};
  std::array<double, 3> localMin{0.0, 0.0, 0.0};
  std::array<double, 3> localMax{4.0, 4.0, 4.0};
  std::array<autopasdm::BoundaryType, 3> boundaries{autopasdm::BoundaryType::reflective, autopasdm::BoundaryType::none,
                                                    autopasdm::BoundaryType::none};
  std::vector<Particle> particles;

  [[nodiscard]] const auto &boundaryTypes() const { return boundaries; }
  [[nodiscard]] const auto &globalBoxMin() const { return globalMin; }
  [[nodiscard]] const auto &globalBoxMax() const { return globalMax; }
  [[nodiscard]] const auto &localBoxMin() const { return localMin; }
  [[nodiscard]] const auto &localBoxMax() const { return localMax; }

  template <class Kernel>
  void applyToOwnedParticlesInRegion(const std::array<double, 3> &regionMin, const std::array<double, 3> &regionMax,
                                     Kernel &&kernel) {
    for (auto &particle : particles) {
      const auto &position = particle.getR();
      const bool insideRegion = position[0] >= regionMin[0] and position[0] < regionMax[0] and
                                position[1] >= regionMin[1] and position[1] < regionMax[1] and
                                position[2] >= regionMin[2] and position[2] < regionMax[2];
      if (insideRegion) {
        kernel(particle);
      }
    }
  }
};

ParticlePropertiesLibraryType makeParticlePropertiesLibrary(
    std::vector<std::array<double, 3>> sitePositions = {{0.0, -0.05, 0.0}, {0.2, -0.05, 0.0}, {0.3, 0.1, 0.0}}) {
  ParticlePropertiesLibraryType properties(1.0);

  properties.addSiteType(0, 1.0);
  properties.addLJParametersToSite(0, 1.0, 1.0);

  properties.addMolType(0, std::vector<std::size_t>(sitePositions.size(), 0), sitePositions,
                        std::array<double, 3>{1.0, 1.0, 1.0});

  properties.calculateMixingCoefficients();
  return properties;
}

TEST(MultisiteReflectiveBoundaryTest, MatchesOriginalMirroredMoleculeInteraction) {
  MultisiteTestContainer container;

  Particle molecule;
  molecule.setID(0);
  molecule.setTypeId(0);
  molecule.setR({0.5, 2.0, 2.0});
  molecule.setQuaternion({1.0, 0.0, 0.0, 0.0});
  molecule.setF({0.0, 0.0, 0.0});
  molecule.setTorque({0.0, 0.0, 0.0});
  container.particles.push_back(molecule);

  auto properties = makeParticlePropertiesLibrary();

  ReflectiveBoundary::apply(container, properties, 1.0);

  const auto &force = container.particles[0].getF();
  const auto &torque = container.particles[0].getTorque();

  // The original implementation sums the pair interactions between all real and mirrored sites.
  EXPECT_NEAR(force[0], 10.511438613856308, 1e-12);
  EXPECT_NEAR(force[1], 0.0, 1e-12);
  EXPECT_NEAR(force[2], 0.0, 1e-12);

  // The site forces also produce the same torque as the original implementation.
  EXPECT_NEAR(torque[0], 0.0, 1e-12);
  EXPECT_NEAR(torque[1], 0.0, 1e-12);
  EXPECT_NEAR(torque[2], 1.0523572813630935, 1e-12);
}

TEST(MultisiteReflectiveBoundaryTest, DiscardsNetAttractiveMirrorInteraction) {
  MultisiteTestContainer container;

  Particle molecule;
  molecule.setID(0);
  molecule.setTypeId(0);
  molecule.setR({0.55, 2.0, 2.0});
  molecule.setQuaternion({1.0, 0.0, 0.0, 0.0});
  molecule.setF({3.0, 4.0, 5.0});
  molecule.setTorque({6.0, 7.0, 8.0});
  container.particles.push_back(molecule);

  auto properties = makeParticlePropertiesLibrary({{0.1, -0.1, 0.0}, {0.1, 0.1, 0.0}});

  ReflectiveBoundary::apply(container, properties, 1.0);

  const auto &force = container.particles[0].getF();
  const auto &torque = container.particles[0].getTorque();
  EXPECT_EQ(force, (std::array<double, 3>{3.0, 4.0, 5.0}));
  EXPECT_EQ(torque, (std::array<double, 3>{6.0, 7.0, 8.0}));
}

TEST(MultisiteReflectiveBoundaryTest, MatchesOriginalMirroredMoleculeInteractionAtUpperWall) {
  MultisiteTestContainer container;

  Particle molecule;
  molecule.setID(0);
  molecule.setTypeId(0);
  molecule.setR({3.5, 2.0, 2.0});
  molecule.setQuaternion({1.0, 0.0, 0.0, 0.0});
  molecule.setF({0.0, 0.0, 0.0});
  molecule.setTorque({0.0, 0.0, 0.0});
  container.particles.push_back(molecule);

  auto properties = makeParticlePropertiesLibrary({{0.0, -0.05, 0.0}, {-0.2, -0.05, 0.0}, {-0.3, 0.1, 0.0}});

  ReflectiveBoundary::apply(container, properties, 1.0);

  const auto &force = container.particles[0].getF();
  const auto &torque = container.particles[0].getTorque();

  EXPECT_NEAR(force[0], -10.511438613856308, 1e-12);
  EXPECT_NEAR(force[1], 0.0, 1e-12);
  EXPECT_NEAR(force[2], 0.0, 1e-12);
  EXPECT_NEAR(torque[0], 0.0, 1e-12);
  EXPECT_NEAR(torque[1], 0.0, 1e-12);
  EXPECT_NEAR(torque[2], -1.0523572813630935, 1e-12);
}

TEST(MultisiteReflectiveBoundaryTest, AddsIndependentWallContributionsAtReflectiveCorner) {
  MultisiteTestContainer container;
  container.boundaries[1] = autopasdm::BoundaryType::reflective;

  Particle molecule;
  molecule.setID(0);
  molecule.setTypeId(0);
  molecule.setR({0.5, 0.5, 2.0});
  molecule.setQuaternion({1.0, 0.0, 0.0, 0.0});
  molecule.setF({0.0, 0.0, 0.0});
  molecule.setTorque({0.0, 0.0, 0.0});
  container.particles.push_back(molecule);

  auto properties = makeParticlePropertiesLibrary();

  ReflectiveBoundary::apply(container, properties, 1.0);

  const auto &force = container.particles[0].getF();
  const auto &torque = container.particles[0].getTorque();

  EXPECT_NEAR(force[0], 10.511438613856308, 1e-12);
  EXPECT_NEAR(force[1], 481.2245621061274, 1e-12);
  EXPECT_NEAR(force[2], 0.0, 1e-12);
  EXPECT_NEAR(torque[0], 0.0, 1e-12);
  EXPECT_NEAR(torque[1], 0.0, 1e-12);
  EXPECT_NEAR(torque[2], 50.98734597954814, 1e-12);
}

}  // namespace
