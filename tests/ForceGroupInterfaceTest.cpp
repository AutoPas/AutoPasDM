#include <gtest/gtest.h>

#include <array>
#include <vector>

#include "autopasdm/timeintegration/ForceGroup.h"
#include "autopasdm/timeintegration/RespaIntegrator.h"
#include "autopasdm/timeintegration/TimeIntegrator.h"

namespace {

using autopasdm::timeintegration::ForceGroup;

struct TestParticle {
  double mass{1.0};
  std::array<double, 3> position{0.0, 0.0, 0.0};
  std::array<double, 3> velocity{0.0, 0.0, 0.0};
  std::array<double, 3> fastForce{0.0, 0.0, 0.0};
  std::array<double, 3> slowForce{0.0, 0.0, 0.0};
};

struct TestParticleSystem {
  std::vector<TestParticle> particles;

  template <class Function>
  void applyToOwnedParticles(Function &&function) {
    for (auto &particle : particles) {
      function(particle);
    }
  }
};

struct TestParticleOperations {
  void resetForce(TestParticle &particle, ForceGroup group) const {
    if (group == ForceGroup::fast) {
      particle.fastForce = {0.0, 0.0, 0.0};
    } else {
      particle.slowForce = {0.0, 0.0, 0.0};
    }
  }

  std::array<double, 3> getForce(const TestParticle &particle, ForceGroup group) const {
    return group == ForceGroup::fast ? particle.fastForce : particle.slowForce;
  }

  void setForce(TestParticle &particle, ForceGroup group, const std::array<double, 3> &force) const {
    if (group == ForceGroup::fast) {
      particle.fastForce = force;
    } else {
      particle.slowForce = force;
    }
  }

  double getMass(const TestParticle &particle) const { return particle.mass; }

  auto getPosition(const TestParticle &particle) const { return particle.position; }

  void setPosition(TestParticle &particle, const std::array<double, 3> &position) const {
    particle.position = position;
  }

  auto getVelocity(const TestParticle &particle) const { return particle.velocity; }

  void setVelocity(TestParticle &particle, const std::array<double, 3> &velocity) const {
    particle.velocity = velocity;
  }
};

struct TestForceEvaluator {
  std::array<double, 3> fastForce{1.0, 0.0, 0.0};
  std::array<double, 3> slowForce{0.0, 2.0, 0.0};
  int fastEvaluations{0};
  int slowEvaluations{0};

  template <class ParticleSystem, class ParticleOperations>
  void evaluate(ParticleSystem &particles, ForceGroup group, ParticleOperations &operations) {
    const auto &force = group == ForceGroup::fast ? fastForce : slowForce;

    if (group == ForceGroup::fast) {
      ++fastEvaluations;
    } else {
      ++slowEvaluations;
    }

    particles.applyToOwnedParticles([&](auto &particle) { operations.setForce(particle, group, force); });
  }
};

static_assert(autopasdm::timeintegration::ForceGroupOperationsFor<TestParticleOperations, TestParticle>);

static_assert(
    autopasdm::timeintegration::ForceGroupEvaluatorFor<TestForceEvaluator, TestParticleSystem, TestParticleOperations>);

TEST(ForceGroupInterfaceTest, KeepsFastAndSlowForcesSeparate) {
  TestParticleSystem particles{{TestParticle{}}};
  TestParticleOperations operations;
  TestForceEvaluator forceEvaluator;

  forceEvaluator.evaluate(particles, ForceGroup::fast, operations);
  forceEvaluator.evaluate(particles, ForceGroup::slow, operations);

  EXPECT_EQ(operations.getForce(particles.particles[0], ForceGroup::fast), (std::array<double, 3>{1.0, 0.0, 0.0}));
  EXPECT_EQ(operations.getForce(particles.particles[0], ForceGroup::slow), (std::array<double, 3>{0.0, 2.0, 0.0}));

  operations.resetForce(particles.particles[0], ForceGroup::fast);

  EXPECT_EQ(operations.getForce(particles.particles[0], ForceGroup::fast), (std::array<double, 3>{0.0, 0.0, 0.0}));
  EXPECT_EQ(operations.getForce(particles.particles[0], ForceGroup::slow), (std::array<double, 3>{0.0, 2.0, 0.0}));
}

TEST(RespaIntegratorTest, AdvancesOneOuterStepWithTwoInnerSteps) {
  TestParticle particle;
  particle.mass = 2.0;
  particle.position = {1.0, 0.0, 0.0};

  TestParticleSystem particles{{particle}};
  TestParticleOperations operations;
  TestForceEvaluator forceEvaluator;
  forceEvaluator.fastForce = {2.0, 0.0, 0.0};
  forceEvaluator.slowForce = {4.0, 0.0, 0.0};

  autopasdm::timeintegration::RespaIntegrator integrator{2};

  EXPECT_FALSE(integrator.initialize(particles, operations, forceEvaluator));
  EXPECT_FALSE(integrator.step(particles, 0.4, operations, forceEvaluator));

  const auto &result = particles.particles[0];

  // Two inner drifts with the velocity-Verlet kick sequence give x = 1.24.
  EXPECT_NEAR(result.position[0], 1.24, 1e-12);

  // The fast and slow force contributions both affect the final velocity.
  EXPECT_NEAR(result.velocity[0], 1.2, 1e-12);

  // Initialization evaluates both groups once. The outer step evaluates the
  // fast group once per inner step and the slow group once at its end.
  EXPECT_EQ(forceEvaluator.fastEvaluations, 3);
  EXPECT_EQ(forceEvaluator.slowEvaluations, 2);
}

}  // namespace