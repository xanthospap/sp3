#ifndef __SP3C_IGS_FLAGS__
#define __SP3C_IGS_FLAGS__

#include "sp3_flag_details.hpp"
#include <cstdint>

namespace dso {

/// @class Sp3Flag A flag to hold all events recorded in the Sp3 for a field
struct Sp3Flag {
  /// Initialize unmarked
  sp3::uitype bits_{0};

  /// @brief Mark flag with an Sp3Event (aka, set the Sp3Event)
  /// @param[in] e The Sp3Event to turn on (aka set)
  void set(Sp3Event e) noexcept { bits_ |= (1 << static_cast<sp3::uitype>(e)); }

  /// @brief Enable bitwise multiple Sp3Event's set.
  /// This function actually enables the following code:
  /// Sp3Flag flag;
  /// flag.set(Sp3Event::bad_abscent_position|Sp3Event::bad_abscent_clock|
  ///           Sp3Event::clock_event|Sp3Event::has_clk_rate_stdev);
  /// @param[in] wf An Sp3FlagWrapper; note that actually we do not need for
  ///            this parameter to be an explicit Sp3FlagWrapper instance; e.g.
  ///            flag.set(Sp3Event::bad_abscent_position|
  ///                     Sp3Event::bad_abscent_clock|
  ///                     Sp3Event::clock_event);
  ///            will call this function because the input parameters will
  ///            result in a Sp3FlagWrapper instance.
  void set(sp3::Sp3FlagWrapper wf) noexcept { bits_ = wf.bits_; }

  /// @brief Un-Mark flag with an Sp3Event (aka, unset the Sp3Event)
  void clear(Sp3Event e) noexcept {
    bits_ &= (~(1 << static_cast<sp3::uitype>(e)));
  }

  /// @brief Clear all Sp3Event's and reset flag to empty/clean
  void reset() noexcept { bits_ = 0; }

  /// @brief Trigger, aka check if an Sp3Event is set
  bool is_set(Sp3Event e) const noexcept {
    return ((bits_ >> static_cast<sp3::uitype>(e)) & 1);
  }

  /// @brief Check if Sp3Flag is clean (no Sp3Event is set)
  bool is_clean() const noexcept { return !bits_; }

  /// @brief Set to reasonable default values
  void set_defaults() noexcept;
}; // Sp3Flag

} // dso

#endif
