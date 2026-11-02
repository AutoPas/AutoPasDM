
# Adding a Custom Time Integrator

AutoPasDM provides an interface for time integration that is independent of the
particle model. Users can implement their own integration scheme and the
particle operations required by that scheme in their own simulator project.

An integrator defines the order of the integration stages and the equations
used to update the simulation state. Particle operations provide access to
particle data, while a force evaluator computes forces for the current particle
configuration.

## Integration interface

A time integrator provides two methods:

- `initialize(particles, operations, forceEvaluator)` prepares the initial
  state before the first time step.
- `step(particles, dt, operations, forceEvaluator)` advances the simulation
  by one time step.

A typical `initialize()` implementation resets the forces and invokes the
force evaluator. A typical `step()` implementation updates the particle state,
prepares for a new force calculation, and invokes the force evaluator again.
The exact sequence depends on the integration scheme.

The integrator may return the result of the force evaluator. For example,
`AutoPasForceEvaluator` returns a status indicating whether AutoPas is tuning.
The application can use this information to decide how to proceed.

## Implementing an integrator

Implement your integrator in your simulator project. For example, you could
create:

```text
my-simulator/
└── src/
    └── integration/
        └── EulerIntegrator.h
```

The following example shows a simple explicit Euler integrator for particles
whose acceleration is computed as force divided by mass. The type belongs to
the application's namespace; it does not need to inherit from an AutoPasDM
base class.

```cpp
#pragma once

#include <functional>
#include <type_traits>

#include "autopas/utils/ArrayMath.h"

namespace mySimulation {

struct EulerIntegrator {
  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  decltype(auto) initialize(ParticleSystem &particles, ParticleOperations &operations,
                            ForceEvaluator &forceEvaluator) const {
    // Clear the initial forces before evaluating them.
    particles.applyToOwnedParticles(
        [&](auto &particle) { operations.resetForce(particle); });

    return std::invoke(forceEvaluator, particles);
  }

  template <class ParticleSystem, class ParticleOperations, class ForceEvaluator>
  auto step(ParticleSystem &particles, double dt, ParticleOperations &operations,
            ForceEvaluator &forceEvaluator) const {
    // Explicit Euler uses the current velocity and force to update the state.
    particles.applyToOwnedParticles([&](auto &particle) {
      const double mass = operations.getMass(particle);
      const auto velocity = operations.getVelocity(particle);
      const auto force = operations.getForce(particle);

      auto position = operations.getPosition(particle);
      position = autopas::utils::ArrayMath::add(
          position, autopas::utils::ArrayMath::mulScalar(velocity, dt));
      operations.setPosition(particle, position);

      auto updatedVelocity = autopas::utils::ArrayMath::add(
          velocity, autopas::utils::ArrayMath::mulScalar(force, dt / mass));
      operations.setVelocity(particle, updatedVelocity);
    });

    // Recompute forces for the updated particle configuration.
    particles.applyToOwnedParticles(
        [&](auto &particle) { operations.resetForce(particle); });

    using ForceResult = std::invoke_result_t<ForceEvaluator &, ParticleSystem &>;
    if constexpr (std::is_void_v<ForceResult>) {
      std::invoke(forceEvaluator, particles);
    } else {
      return std::invoke(forceEvaluator, particles);
    }
  }
};

}  // namespace mySimulation
```

For this example, the update equations are:

$$
\mathbf{r}_{n+1} = \mathbf{r}_n + \Delta t\,\mathbf{v}_n
$$

$$
\mathbf{v}_{n+1} = \mathbf{v}_n +
\Delta t\,\frac{\mathbf{F}_n}{m}
$$

The force evaluator is called after the update so that the stored forces
correspond to the new particle positions and are ready for the next time step.

This example illustrates the interface. Explicit Euler is a first-order
method; choose an integration scheme appropriate for your physical model and
accuracy requirements.

## Providing particle operations

The integrator calls methods on the `ParticleOperations` object to access
model-specific data. The Euler example requires these operations:

```cpp
resetForce(particle);
getMass(particle);
getPosition(particle);
setPosition(particle, position);
getVelocity(particle);
setVelocity(particle, velocity);
getForce(particle);
```

These operations should be implemented in the simulator project as well. For
example, an MD operations object can obtain a particle's mass from a particle
properties library and read or update position, velocity, and force through
the particle's methods.

`resetForce()` should set the force to the correct starting value before each
force evaluation. This may be zero, or it may include an external force that
must remain present while interaction forces are added.

The required operations depend on the integrator. Another integrator may use
different data or require additional operations. The `TimeIntegrator` concept
checks the integrator's `initialize()` and `step()` interface. Operations used
inside those methods are checked when the corresponding templates are
instantiated.

## Using the integrator

Include the integrator from your simulator project and construct it alongside
your particle operations and force evaluator. Call `initialize()` once before
the first time step:

```cpp
#include "autopasdm/integration/TimeIntegrator.h"
#include "integration/EulerIntegrator.h"

// The particle operations and force evaluator are defined by the application.
mySimulation::EulerIntegrator integrator;

static_assert(autopasdm::integration::TimeIntegrator<
              decltype(integrator), decltype(particles),
              decltype(particleOperations), decltype(forceEvaluator)>);

const auto initializationStatus =
    integrator.initialize(particles, particleOperations, forceEvaluator);

// Handle the status if required by the application.
(void)initializationStatus;

for (int step = 0; step < numberOfSteps; ++step) {
  const auto status =
      integrator.step(particles, deltaT, particleOperations, forceEvaluator);

  // Handle the status if required by the application.
  (void)status;
}
```

If the force evaluator returns no value, `initialize()` and `step()` also
return no value. If it returns a status, the integrator can forward that status
to the caller.

The `TimeIntegrator` concept checks the basic integrator interface at compile
time. It does not require inheritance from an AutoPasDM type.

## Particle and force handling

The integrator updates owned particles through
`particles.applyToOwnedParticles()`. Halo particles are managed by the
particle system and the force evaluation workflow; an integrator should not
update halos as if they were locally owned particles.

The force evaluator is responsible for computing forces and performing any
interaction preparation required by the particle system. The integrator
controls when the evaluator is called.