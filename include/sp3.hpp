/** @file
 * Define class to handle interaction with Sp3-C ephemerides files.
 */

#ifndef __SP3C_IGS_FILE__
#define __SP3C_IGS_FILE__

#include "datetime/calendar.hpp"
#include "satellite.hpp"
#include "sp3flag.hpp"
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <vector>
#ifdef DEBUG
#include "datetime/datetime_write.hpp"
#endif

namespace dso {

/** @class Sp3DataBlock
 * Instances of this class, hold Sp3 data records for one block (aka one
 * epoch) and one satellite.
 */
struct Sp3DataBlock {
  dso::datetime<dso::nanoseconds> t{dso::datetime<dso::nanoseconds>::min()};
  double state[8];      /** [ X, Y, Z, clk, Vx, Vy, Vz, Vc ] */
  double state_sdev[8]; /** following state__ */
  Sp3Flag flag;         /** flag for state */
};                      /* Sp3DataBlock */

class Sp3c {
public:
  /** Let's not write this more than once. */
  typedef std::ifstream::pos_type pos_type;

  /** @brief Constructor from filename */
  explicit Sp3c(const char *fn);

  /** @brief Copy not allowed ! */
  Sp3c(const Sp3c &) = delete;

  /** @brief Assignment not allowed ! */
  Sp3c &operator=(const Sp3c &) = delete;

  /** @brief Move Constructor. */
  Sp3c(Sp3c &&a) noexcept(
      std::is_nothrow_move_constructible<std::ifstream>::value) = default;

  /** @brief Move assignment operator. */
  Sp3c &operator=(Sp3c &&a) noexcept(
      std::is_nothrow_move_assignable<std::ifstream>::value) = default;

  /** get the Sp3 interval */
  auto interval() const noexcept { return interval__; }

  /** get the number of epochs included in the file */
  auto num_epochs() const noexcept { return num_epochs__; }

  /** get the initial epoch (datetime) in the Sp3 file */
  auto start_epoch() const noexcept { return start_epoch__; }

  /** Rewind to the start of data blocks (i.e. just after the header) */
  void rewind() noexcept { __istream.seekg(__end_of_head, std::ios::beg); }

  /** @brief Read the next data block and parse holding for a given SV
   * @param[in] satid The SV to collect records for
   * @param[out] block An Sp3DataBlock instance; if we encounter records
   *            for SV satid, block will be filled with the parsed values.
   *            Use the block's flag member to check which values where
   *            actually parsed (if any at all).
   *            Units are :
   *            * State Vector: [km] and [dm/sec]
   *            * Std. Deviations: [mm] and 10e-4[mm/sec]
   * @return -1: EOF encountered
   *          0: All ok
   *         >0: ERROR
   */
  int get_next_data_block(sp3::SatelliteId satid, Sp3DataBlock &block) noexcept;

  /** Assuming we are in a positio in the file where the next line to be read
   * is an epoch header line; resolve the date, but do not progress the
   * stream position.
   */
  int peak_next_data_block(dso::datetime<dso::nanoseconds> &t) noexcept;

  /** @brief Check if a given SV in included in the Sp3 (i.e. is included in
   *        the instance's sat_vec__ member).
   *
   * @note  It is assumed that the header of the file is already parsed; this
   *        should always be the case, since the file's header is parsed when
   *        it is constructed.
   *
   * @param[in] satid The satelliteid we are searching for, as recorded in Sp3
   *            (aka a 3-char id, as'G01', 'R27', etc)
   * @return True if satellite is included in the instance's sat_vec__; false
   *         otherwise.
   */
  bool has_sv(sp3::SatelliteId satid) const noexcept {
    using sp3::SatelliteId;
    return std::find_if(sat_vec__.cbegin(), sat_vec__.cend(),
                        [satid](const SatelliteId &s) { return s == satid; }) !=
           sat_vec__.cend();
  }

  /** @brief Number of satellites in sp3 file
   *
   * @note  It is assumed that the header of the file is already parsed; this
   *        should always be the case, since the file's header is parsed when
   *        it is constructed.
   */
  int num_sats() const noexcept { return sat_vec__.size(); }

  /** @brief Return the vector of satellites included in the sp3 file */
  std::vector<sp3::SatelliteId> sattellite_vector() const noexcept {
    return sat_vec__;
  }

  /** @brief Return the vector of satellites included in the sp3 file */
  std::vector<sp3::SatelliteId> &sattellite_vector() noexcept {
    return sat_vec__;
  }

#ifdef DEBUG
  void print_members() const noexcept;
#endif

private:
  /** @brief Read sp3c header; assign info */
  int read_header() noexcept;

  /** @brief Resolve an Epoch Header Record line */
  int resolve_epoch_line(dso::datetime<dso::nanoseconds> &t) noexcept;

  /** @brief Get and resolve the next Position and Clock Record */
  int get_next_position(sp3::SatelliteId &sat, double &xkm, double &ykm,
                        double &zkm, double &clk, double &xstdv, double &ystdv,
                        double &zstdv, double &cstdv, Sp3Flag &flag,
                        const sp3::SatelliteId *wsat = nullptr) noexcept;

  /** @brief Get and resolve the next Velocity and ClockRate-of-Change Record */
  int get_next_velocity(sp3::SatelliteId &sat, double &xkm, double &ykm,
                        double &zkm, double &clk, double &xstdv, double &ystdv,
                        double &zstdv, double &cstdv, Sp3Flag &flag,
                        const sp3::SatelliteId *wsat = nullptr) noexcept;

  /** The name of the file */
  std::string __filename;
  /** The infput (file) stream */
  std::ifstream __istream;
  /** the version 'c' or 'd' */
  char version__;
  /** Start epoch */
  dso::datetime<dso::nanoseconds> start_epoch__;
  /** Number of epochs in file */
  int num_epochs__,
      /** Number od SVs in file */
      num_sats__;
  /** Coordinate system (last char always '\0') */
  char crd_sys__[6] = {'\0'},
       /** Orbit type (last char always '\0') */
      orb_type__[4] = {'\0'},
       /** Agency (last char always '\0') */
      agency__[5] = {'\0'},
       /** Time system (last char always '\0') */
      time_sys__[4] = {'\0'};
  /** Epoch interval */
  dso::nanoseconds interval__;
  /** Mark the 'END OF HEADER' field */
  pos_type __end_of_head;
  /** Vector of satellite id's */
  std::vector<sp3::SatelliteId> sat_vec__;
  /** floating point base for position std. dev (mm or 10**-4 mm/sec) */
  double fpb_pos__,
      /** floating point base for clock std. dev (psec or 10**-4 psec/sec) */
      fpb_clk__;
}; /* class Sp3c */

/** Utility class, to iterate through the data blocks of an Sp3 file */
class Sp3Iterator {
  Sp3c *sp3_;
  sp3::SatelliteId id_;
  Sp3DataBlock block_;

public:
  Sp3Iterator(Sp3c &sp3) : sp3_(&sp3) {
    sp3_->rewind();
    if (sp3_->get_next_data_block(id_, block_)) {
      throw std::runtime_error(
          "ERROR Failed to create Sp3Iterator instance!\n");
    }
  };

  const Sp3DataBlock &data_block() const noexcept { return block_; }

  void begin() {
    sp3_->rewind();
    if (sp3_->get_next_data_block(id_, block_)) {
      throw std::runtime_error(
          "ERROR Failed to create Sp3Iterator instance!\n");
    }
    return;
  }

  int advance() noexcept { return sp3_->get_next_data_block(id_, block_); }

  dso::datetime<dso::nanoseconds> current_time() const noexcept {
    return block_.t;
  }

  int peak_next_epoch(dso::datetime<dso::nanoseconds> &t) const noexcept {
    return sp3_->peak_next_data_block(t);
  }

  int goto_epoch(const dso::datetime<dso::nanoseconds> &t) noexcept {
    int error = 0, advance_er = 0;
    dso::datetime<dso::nanoseconds> ct = block_.t;

    if (block_.t < t) {
      // peak next epoch from next header
      while (!advance_er && !(error = peak_next_epoch(ct))) {
        // if next epoch <  requested, read it in
        if (ct < t) {
          advance_er = advance();
        } else {
          break;
        }
      }
      if (error < 0) { // EOF encountered
        return -1;
      }
      if (error + advance_er)
        return error + advance_er;
    } else {
      this->begin();
      if (block_.t > t)
        return 10;
      return goto_epoch(t);
    }

    return error + advance_er;
  }

}; /* Sp3Iterator */

} /* namespace dso */

#endif
