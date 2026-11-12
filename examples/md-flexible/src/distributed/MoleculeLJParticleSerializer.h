#pragma once

#include <array>
#include <cstddef>

#include "autopas/AutoPasDecl.h"
#include "autopasdm/ParticleSerialization.h"
#include "molecularDynamicsLibrary/MoleculeLJ.h"
#include "molecularDynamicsLibrary/MultisiteMoleculeLJ.h"

namespace autopasdm {

/**
 * MPI serialization for md-flexible's single-site particle type.
 *
 * The application defines which particle state has to be transferred. This
 * keeps AutoPasDM independent of md-flexible and of MoleculeLJ.
 */
template <>
struct ParticleSerializer<mdLib::MoleculeLJ> {
  /** Serialized state required to reconstruct a MoleculeLJ on another rank. */
  struct Message {
    std::array<double, 3> r{};
    std::array<double, 3> v{};
    std::array<double, 3> f{};
    std::array<double, 3> oldF{};
    unsigned long id{0};
    std::size_t typeId{0};
  };

  /**
   * Serialize one MoleculeLJ particle for distributed communication.
   *
   * @param particle Particle to serialize.
   */
  static Message pack(const mdLib::MoleculeLJ &particle) {
    return Message{particle.getR(),    particle.getV(),  particle.getF(),
                   particle.getOldF(), particle.getID(), particle.getTypeId()};
  }

  /**
   * Reconstruct one MoleculeLJ particle from its serialized state.
   *
   * @param message Serialized particle state.
   */
  static mdLib::MoleculeLJ unpack(const Message &message) {
    mdLib::MoleculeLJ particle;
    particle.setR(message.r);
    particle.setV(message.v);
    particle.setF(message.f);
    particle.setOldF(message.oldF);
    particle.setID(message.id);
    particle.setTypeId(message.typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    return particle;
  }
};

/**
 * MPI serialization for md-flexible's multi-site particle type.
 */
template <>
struct ParticleSerializer<mdLib::MultisiteMoleculeLJ> {
  struct Message {
    std::array<double, 3> r{};
    std::array<double, 3> v{};
    std::array<double, 3> f{};
    std::array<double, 3> oldF{};
    std::array<double, 4> quaternion{};
    std::array<double, 3> angularVelocity{};
    std::array<double, 3> torque{};
    unsigned long id{0};
    std::size_t typeId{0};
  };

  static Message pack(const mdLib::MultisiteMoleculeLJ &particle) {
    return {particle.getR(),      particle.getV(),          particle.getF(),
            particle.getOldF(),   particle.getQuaternion(), particle.getAngularVel(),
            particle.getTorque(), particle.getID(),         particle.getTypeId()};
  }

  static mdLib::MultisiteMoleculeLJ unpack(const Message &message) {
    mdLib::MultisiteMoleculeLJ particle;
    particle.setR(message.r);
    particle.setV(message.v);
    particle.setF(message.f);
    particle.setOldF(message.oldF);
    particle.setQuaternion(message.quaternion);
    particle.setAngularVel(message.angularVelocity);
    particle.setTorque(message.torque);
    particle.setID(message.id);
    particle.setTypeId(message.typeId);
    particle.setOwnershipState(autopas::OwnershipState::owned);
    return particle;
  }
};

}  // namespace autopasdm
