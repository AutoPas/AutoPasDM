#pragma once

#include <cstddef>
#include <functional>
#include <stdexcept>
#include <type_traits>

#include "autopas/utils/ArrayMath.h"
#include "autopasdm/timeintegration/ForceGroup.h"

namespace autopasdm::timeintegration {

class RespaIntegrator {
 public:
  explicit RespaIntegrator(std::size_t innerSteps) : _innerSteps(innerSteps) {
    if (_innerSteps == 0) {
      throw std::invalid_argument("RespaIntegrator requires at least one inner step.");
    }
  }

  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  bool initialize(ParticleSystem &particles, ParticleOperations &operations, ForceEvaluator &forceEvaluator) const {
    // Compute both force groups for the initial particle configuration.
    const bool fastTuning = evaluateForceGroup(particles, operations, forceEvaluator, ForceGroup::fast);
    const bool slowTuning = evaluateForceGroup(particles, operations, forceEvaluator, ForceGroup::slow);

    return fastTuning || slowTuning;
  }

  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  bool step(ParticleSystem &particles, double dt, ParticleOperations &operations,
            ForceEvaluator &forceEvaluator) const {
    bool tuning = false;

    // Apply the first half of the slow-force velocity update.
    updateVelocities(particles, operations, ForceGroup::slow, 0.5 * dt);

    const double innerDt = dt / static_cast<double>(_innerSteps);

    for (std::size_t innerStep = 0; innerStep < _innerSteps; ++innerStep) {
      // Apply the first half of the fast-force velocity update.
      updateVelocities(particles, operations, ForceGroup::fast, 0.5 * innerDt);

      // Drift positions using the updated velocities.
      updatePositions(particles, operations, innerDt);

      // Recompute fast forces at the updated positions.
      const bool fastTuning = evaluateForceGroup(particles, operations, forceEvaluator, ForceGroup::fast);
      tuning = fastTuning || tuning;

      // Complete the fast-force velocity update.
      updateVelocities(particles, operations, ForceGroup::fast, 0.5 * innerDt);
    }

    // Recompute slow forces once at the end of the outer step.
    const bool slowTuning = evaluateForceGroup(particles, operations, forceEvaluator, ForceGroup::slow);
    tuning = slowTuning || tuning;

    // Complete the slow-force velocity update.
    updateVelocities(particles, operations, ForceGroup::slow, 0.5 * dt);

    return tuning;
  }

 private:
  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  static bool evaluateForceGroup(ParticleSystem &particles, ParticleOperations &operations,
                                 ForceEvaluator &forceEvaluator, ForceGroup group) {
    // Clear only the force buffer for the group being evaluated.
    particles.applyToOwnedParticles([&](auto &particle) { operations.resetForce(particle, group); });

    using EvaluationResult = decltype(forceEvaluator.evaluate(particles, group, operations));

    if constexpr (std::is_void_v<EvaluationResult>) {
      forceEvaluator.evaluate(particles, group, operations);
      return false;
    } else {
      // Non-void evaluation results must be convertible to bool.
      return static_cast<bool>(forceEvaluator.evaluate(particles, group, operations));
    }
  }

  template <class ParticleSystem, class ParticleOperations>
  static void updateVelocities(ParticleSystem &particles, ParticleOperations &operations, ForceGroup group,
                               double timeFactor) {
    particles.applyToOwnedParticles([&](auto &particle) {
      const double mass = operations.getMass(particle);
      const auto force = operations.getForce(particle, group);
      auto velocity = operations.getVelocity(particle);

      velocity =
          autopas::utils::ArrayMath::add(velocity, autopas::utils::ArrayMath::mulScalar(force, timeFactor / mass));

      operations.setVelocity(particle, velocity);
    });
  }

  template <class ParticleSystem, class ParticleOperations>
  static void updatePositions(ParticleSystem &particles, ParticleOperations &operations, double dt) {
    particles.applyToOwnedParticles([&](auto &particle) {
      const auto velocity = operations.getVelocity(particle);
      auto position = operations.getPosition(particle);

      position = autopas::utils::ArrayMath::add(position, autopas::utils::ArrayMath::mulScalar(velocity, dt));

      operations.setPosition(particle, position);
    });
  }

  std::size_t _innerSteps;
};

}  // namespace autopasdm::timeintegration