#pragma once

#include <functional>
#include <type_traits>

#include "autopas/utils/ArrayMath.h"
#include "autopasdm/timeintegration/TimeIntegrator.h"

namespace autopasdm::timeintegration {

struct StormerVerlet {
  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  decltype(auto) initialize(ParticleSystem &particles, ParticleOperations &operations,
                            ForceEvaluator &forceEvaluator) const {
    // Clear the initial forces before evaluating them.
    particles.applyToOwnedParticles([&](auto &particle) { operations.resetForce(particle); });

    return std::invoke(forceEvaluator, particles);
  }

  template <class ParticleSystem, class ParticleOperations>
  void beginStep(ParticleSystem &particles, double dt, ParticleOperations &operations) const {
    // Update positions using the current forces, then reset the force buffer
    // for the upcoming force evaluation.
    particles.applyToOwnedParticles([&](auto &particle) {
      const double mass = operations.getMass(particle);
      const auto force = operations.getForce(particle);
      const auto velocity = operations.getVelocity(particle);

      operations.setOldForce(particle, force);
      operations.resetForce(particle);

      auto displacement =
          autopas::utils::ArrayMath::add(autopas::utils::ArrayMath::mulScalar(velocity, dt),
                                         autopas::utils::ArrayMath::mulScalar(force, dt * dt / (2.0 * mass)));

      auto position = operations.getPosition(particle);
      position = autopas::utils::ArrayMath::add(position, displacement);
      operations.setPosition(particle, position);
    });
  }

  template <class ParticleSystem, class ParticleOperations>
  void finishStep(ParticleSystem &particles, double dt, ParticleOperations &operations) const {
    // Complete the velocity update using the old and newly computed forces.
    particles.applyToOwnedParticles([&](auto &particle) {
      const double mass = operations.getMass(particle);
      const auto oldForce = operations.getOldForce(particle);
      const auto force = operations.getForce(particle);
      auto velocity = operations.getVelocity(particle);

      const auto forceSum = autopas::utils::ArrayMath::add(oldForce, force);
      const auto velocityChange = autopas::utils::ArrayMath::mulScalar(forceSum, dt / (2.0 * mass));

      velocity = autopas::utils::ArrayMath::add(velocity, velocityChange);
      operations.setVelocity(particle, velocity);
    });
  }

  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  auto step(ParticleSystem &particles, double dt, ParticleOperations &operations,
            ForceEvaluator &forceEvaluator) const {
    beginStep(particles, dt, operations);

    using ForceResult = std::invoke_result_t<ForceEvaluator &, ParticleSystem &>;
    if constexpr (std::is_void_v<ForceResult>) {
      std::invoke(forceEvaluator, particles);
      finishStep(particles, dt, operations);
    } else {
      auto result = std::invoke(forceEvaluator, particles);
      finishStep(particles, dt, operations);
      return result;
    }
  }
};

}  // namespace autopasdm::timeintegration