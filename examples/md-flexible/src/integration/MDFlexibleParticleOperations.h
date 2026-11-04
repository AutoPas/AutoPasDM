

#pragma once

#include <array>
#include <atomic>
#include <cmath>
#include <iostream>
#include <stdexcept>

#include "autopas/utils/ArrayMath.h"
#include "autopas/utils/WrapOpenMP.h"
#include "src/TypeDefinitions.h"

struct MDFlexibleParticleOperations {
  const ParticlePropertiesLibraryType &particlePropertiesLibrary;
  std::array<double, 3> globalForce;
  double verletSkinRadius;
  unsigned int verletRebuildFrequency;
  bool fastParticlesThrow;

  // Set if any particle exceeds the Verlet displacement limit during the parallel update.
  mutable std::atomic_bool _particleTooFast{false};

  template <class Particle>
  void resetForce(Particle &particle) const {
    // Set the force to the configured external force before pair forces are added.
    particle.setF(globalForce);
  }

  template <class Particle>
  double getMass(const Particle &particle) const {
    return particlePropertiesLibrary.getMolMass(particle.getTypeId());
  }

  template <class Particle>
  auto getPosition(const Particle &particle) const {
    return particle.getR();
  }

  template <class Particle, class Vector>
  void setPosition(Particle &particle, const Vector &position) const {
    const auto displacement = autopas::utils::ArrayMath::sub(position, particle.getR());

#ifndef AUTOPAS_ENABLE_DYNAMIC_CONTAINERS
    const double maxAllowedDistanceMoved = verletSkinRadius / verletRebuildFrequency / 2.;
    const double maxAllowedDistanceMovedSquared = maxAllowedDistanceMoved * maxAllowedDistanceMoved;

    if (not particle.addRDistanceCheck(displacement, maxAllowedDistanceMovedSquared) and
        maxAllowedDistanceMovedSquared > 0.) {
      const auto distanceMoved = std::sqrt(autopas::utils::ArrayMath::dot(displacement, displacement));
      AUTOPAS_OPENMP(critical)
      std::cerr << "A particle moved farther than verletSkinPerTimestep/2: " << distanceMoved << " > "
                << maxAllowedDistanceMoved << "\n"
                << particle << "\nNew Position: [" << position[0] << ", " << position[1] << ", " << position[2] << "]"
                << std::endl;
      _particleTooFast.store(true, std::memory_order_relaxed);
    }
#else
    particle.addR(displacement);
#endif
  }

  void checkForFastParticles() const {
    const bool particleTooFast = _particleTooFast.exchange(false, std::memory_order_relaxed);
    if (particleTooFast and fastParticlesThrow) {
      throw std::runtime_error("At least one particle was too fast!");
    }
  }

  template <class Particle>
  auto getVelocity(const Particle &particle) const {
    return particle.getV();
  }

  template <class Particle, class Vector>
  void setVelocity(Particle &particle, const Vector &velocity) const {
    particle.setV(velocity);
  }

  template <class Particle>
  auto getForce(const Particle &particle) const {
    return particle.getF();
  }

  template <class Particle>
  auto getOldForce(const Particle &particle) const {
    return particle.getOldF();
  }

  template <class Particle, class Vector>
  void setOldForce(Particle &particle, const Vector &force) const {
    particle.setOldF(force);
  }
};
