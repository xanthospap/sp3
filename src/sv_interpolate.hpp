#ifndef __SV_SP3_INTERPOLATION_HPP__
#define __SV_SP3_INTERPOLATION_HPP__

#include "sp3.hpp"
#include <stdexcept>

namespace dso {
namespace sp3 {
int neville_interpolation(double x, double &y, double &dy, const double *xx,
                          const double *yy, int array_size, int mm,
                          int from_index = 0, double *cws = nullptr,
                          double *dws = nullptr) noexcept;
} // namespace sp3

constexpr dso::datetime_interval<dso::microseconds> MAX_MICROSEC_{dso::modified_julian_day(0), dso::microseconds{60 * 60 * 2 * 1'000'000L}}; // 2 hours
constexpr int MIN_INTERPOLATION_PTS = 4;
class SvInterpolator {
public:
  sp3::SatelliteId svid;
  int num_dpts{0};
  Sp3c *sp3{nullptr};
  Sp3DataBlock *data{nullptr};
  int last_index{0};

  // fill in the data instance using an sp3 instance
  int feed_from_sp3() noexcept {
    if (!sp3)
      return 1;

    if (!sp3->num_epochs()) {
      fprintf(stderr,
              "[ERROR] Sp3 instance has no epochs stored! Did you forget to "
              "read it's header? (traceback: %s)\n",
              __func__);
      return 1;
    }

    if (!sp3->has_sv(svid)) {
      fprintf(stderr,
              "[ERROR] Sp3 instance has no data records for the requested SV "
              "(traceback: %s)\n",
              __func__);
      return 2;
    }

    // allocate enough space; to be safe, use the number of epochs in the sp3
    // file, even though some records may be missing
    if (data)
      delete[] data;
    num_dpts = 0;
    data = new Sp3DataBlock[sp3->num_epochs()];

    // read the sp3 file through and grap data for the sv
    Sp3DataBlock block;
    sp3->rewind();
    int error, idx = 0;
    while (!(error = sp3->get_next_data_block(svid, block))) {
      // do not include data point if position and clock are missing
      if (!(block.flag.is_set(Sp3Event::bad_abscent_position) &&
            block.flag.is_set(Sp3Event::bad_abscent_clock)))
        data[idx++] = block;
    }

    // check for error while parsing
    if (error > 0) {
      fprintf(stderr,
              "[ERROR] Failed parsing sp3 file for the requested SV data "
              "(traceback: %s)\n",
              __func__);
      delete[] data;
    }

    num_dpts = (error > 0) ? 0 : idx;

    return error > 0;
  }

public:
  SvInterpolator(sp3::SatelliteId sid) noexcept : svid(sid){};
  SvInterpolator(sp3::SatelliteId sid, Sp3c &sp3obj) : svid(sid), sp3(&sp3obj) {
    if (int error = feed_from_sp3(); error) {
      throw std::runtime_error("[ERROR] Failed creating SvInterpolator "
                               "instance from Sp3 Error Code: " +
                               std::to_string(error));
    }
  };

  int index_hunt(dso::datetime<dso::microseconds> t) noexcept {
    int index = -1;
    // find the interval within the instance's block
    // |------|------|---X---|X------|-------|-------|-------|  X  > t-axis
    //               ^   1    2                                ^
    //               | last_index (A)                        | last_index (B)
    if (data[last_index].t <= t) {
      if (last_index == num_dpts - 1) { // case (B)
        return (index = last_index);
      } else if (data[last_index + 1].t > t) { // case (A1)
        return (index = last_index);
      } else { // try next block (if any)
        if (last_index + 2 < num_dpts - 1) {
          if (data[last_index + 1].t <= t && data[last_index + 2].t > t) {
            ++last_index;
            return (index = last_index);
          }
        } else if (last_index + 2 == num_dpts) {
          if (data[last_index + 1].t <= t)
            ++last_index;
          return (index = last_index);
        }
      }
    }

    auto it = std::lower_bound(
        data, data + num_dpts, t,
        [](const Sp3DataBlock &block, const dso::datetime<dso::microseconds>& tt) {
          return block.t < tt;
        });
    return static_cast<int>(data - it);
  }

  int interpolate_at(dso::datetime<dso::microseconds> t, double *result) noexcept {
    int index = index_hunt(t);
    // int lpts = index;
    // int rpts = num_dpts - index - 1;

    int start = index;
    while (start > 0 && t - data[start].t < MAX_MICROSEC_)
      --start;

    int stop = index;
    while (stop < num_dpts - 1 && data[stop].t - t < MAX_MICROSEC_)
      ++stop;

    if (stop - start + 1 < MIN_INTERPOLATION_PTS) {
      fprintf(stderr,
              "[ERROR] Cannot interpolate due to too few data points "
              "(traceback: %s)\n",
              __func__);
      return 1;
    }

    int size = stop - start + 1;
    double *td = new double[size];
    double *xd = new double[size];
    double *yd = new double[size];
    double *zd = new double[size];

    auto start_t = sp3->start_epoch();
    for (int i = 0; i < size; i++) {
      td[i] = data[start + i].t.delta_date(start_t).as_mjd();
      xd[i] = data[start + i].state[0];
      yd[i] = data[start + i].state[1];
      zd[i] = data[start + i].state[2];
    }
    double tx = t.delta_date(start_t).as_mjd();

    double y, dy;
    if (sp3::neville_interpolation(tx, y, dy, td, xd, size, size)) {
      fprintf(stderr, "[ERROR] Neville algorithm failed!\n");
      return 5;
    }

    *result = y;

    return 0;
  }
}; // SvInterpolator

} // namespace dso

#endif