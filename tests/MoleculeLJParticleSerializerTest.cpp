#include <gtest/gtest.h>

#include <type_traits>

#include "src/distributed/MoleculeLJParticleSerializer.h"

TEST(MoleculeLJParticleSerializerTest, MultiSiteRoundTripPreservesState) {
  using Particle = mdLib::MultisiteMoleculeLJ;
  using Serializer = autopasdm::ParticleSerializer<Particle>;

  static_assert(std::is_trivially_copyable_v<typename Serializer::Message>);

  Particle original;
  original.setR({1., 2., 3.});
  original.setV({4., 5., 6.});
  original.setF({7., 8., 9.});
  original.setOldF({10., 11., 12.});
  original.setQuaternion({0.1, 0.2, 0.3, 0.4});
  original.setAngularVel({13., 14., 15.});
  original.setTorque({16., 17., 18.});
  original.setID(42);
  original.setTypeId(7);

  const auto restored = Serializer::unpack(Serializer::pack(original));

  EXPECT_EQ(restored.getR(), original.getR());
  EXPECT_EQ(restored.getV(), original.getV());
  EXPECT_EQ(restored.getF(), original.getF());
  EXPECT_EQ(restored.getOldF(), original.getOldF());
  EXPECT_EQ(restored.getQuaternion(), original.getQuaternion());
  EXPECT_EQ(restored.getAngularVel(), original.getAngularVel());
  EXPECT_EQ(restored.getTorque(), original.getTorque());
  EXPECT_EQ(restored.getID(), original.getID());
  EXPECT_EQ(restored.getTypeId(), original.getTypeId());
}