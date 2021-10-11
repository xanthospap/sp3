#ifndef __SP3C_SATELLITE_FILE__
#define __SP3C_SATELLITE_FILE__

#include "ggdatetime/dtcalendar.hpp"

namespace dso {

/// @brief Number of characters used to describe a given sattelite vehicle
constexpr int SAT_ID_CHARS = 3;

/// @brief Number of characters used to describe a given sattelite vehicle as a
/// null terminating string
constexpr int SAT_ID_MAX_CHARS = SAT_ID_CHARS + 1;

/// @class Satellite ID as denoted in an Sp3 file
struct SatelliteId {
  /// @brief The id of the space vehicle (3chars plus the null terminating char)
  char id[SAT_ID_MAX_CHARS] = {'\0'};
  
  /// @brief Constructor from a c-string; this will only copy the first 
  /// SAT_ID_CHARS (from the input string) to the instance's id
  explicit SatelliteId(const char *str = nullptr) noexcept {
    if (str) set_id(str);
  }
  
  /// @brief Set id from a c-string; this will only copy the first 
  /// SAT_ID_CHARS (from the input string) to the instance's id
  void set_id(const char* str) noexcept {
    std::memcpy(id, str, SAT_ID_CHARS);
  }

  /// @brief Compare satellite id's
  bool operator==(const SatelliteId &s) const noexcept {
    return !std::strncmp(id, s.id, SAT_ID_CHARS);
  }

  /// @brief Compare satellite id's
  bool operator!=(const SatelliteId &s) const noexcept {
    return !this->operator==(s);
  }

  /// @brief Cast the id to an std::string
  std::string to_string() const noexcept { return std::string(id, SAT_ID_CHARS); }
};

}//dso
#endif