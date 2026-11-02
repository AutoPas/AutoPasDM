#pragma once

#include <concepts>
#include <functional>

#include "autopasdm/timeintegration/ForceGroup.h"

namespace autopasdm::timeintegration {

template <class ForceEvaluator, class ParticleSystem>
concept ForceEvaluatorFor = requires(ForceEvaluator &forceEvaluator, ParticleSystem &particles) {
  std::invoke(forceEvaluator, particles);
};

template <class Integrator, class ParticleSystem, class ParticleOperations, class ForceEvaluator>
concept TimeIntegrator = ForceEvaluatorFor<ForceEvaluator, ParticleSystem> and
    requires(Integrator &integrator, ParticleSystem &particles, ParticleOperations &operations,
             ForceEvaluator &forceEvaluator, double dt) {
  integrator.initialize(particles, operations, forceEvaluator);
  integrator.step(particles, dt, operations, forceEvaluator);
};

template <class ForceEvaluator, class ParticleSystem, class ParticleOperations>
concept ForceGroupEvaluatorFor = requires(ForceEvaluator &forceEvaluator, ParticleSystem &particles,
                                          ParticleOperations &operations) {
  forceEvaluator.evaluate(particles, ForceGroup::fast, operations);
  forceEvaluator.evaluate(particles, ForceGroup::slow, operations);
};

template <class ParticleOperations, class Particle>
concept ForceGroupOperationsFor = requires(ParticleOperations &operations, Particle &particle) {
  operations.resetForce(particle, ForceGroup::fast);
  operations.resetForce(particle, ForceGroup::slow);

  operations.getForce(particle, ForceGroup::fast);
  operations.getForce(particle, ForceGroup::slow);

  operations.setForce(particle, ForceGroup::fast, operations.getForce(particle, ForceGroup::fast));
  operations.setForce(particle, ForceGroup::slow, operations.getForce(particle, ForceGroup::slow));
};

template <class ParticleOperations, class Particle>
concept RespaParticleOperationsFor = ForceGroupOperationsFor<ParticleOperations, Particle> &&
    requires(ParticleOperations &operations, Particle &particle) {
  { operations.getMass(particle) } -> std::convertible_to<double>;

  operations.getPosition(particle);
  operations.setPosition(particle, operations.getPosition(particle));

  operations.getVelocity(particle);
  operations.setVelocity(particle, operations.getVelocity(particle));
};

}  // namespace autopasdm::timeintegration