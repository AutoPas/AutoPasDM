#pragma once

namespace autopasdm::timeintegration {

/**
 * Identifies the timescale of a force contribution for a multiple-time-step
 * integrator such as r-RESPA.
 */
enum class ForceGroup {
  fast,
  slow,
};

}  // namespace autopasdm::timeintegration