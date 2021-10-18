#ifndef __SV_SP3_INTERPOLATION_HPP__
#define __SV_SP3_INTERPOLATION_HPP__

#include "sp3.hpp"
#include <stdexcept>
#ifdef DEBUG
#include <chrono>
#endif

namespace dso {
namespace sp3 {
int neville_interpolation(double x, double &y, double &dy, const double *xx,
                          const double *yy, int array_size, int mm,
                          int from_index = 0, double *cws = nullptr,
                          double *dws = nullptr) noexcept;

int neville_interpolation3(double t, double *estimates, double *destimates,
                           const double *__restrict__ tt,
                           const double *__restrict__ xx,
                           const double *__restrict__ yy,
                           const double *__restrict__ zz, int array_size,
                           int mm, int from_index, double *workspace) noexcept;
} // namespace sp3

constexpr dso::datetime_interval<dso::microseconds> MAX_MICROSEC_{
    dso::modified_julian_day(0),
    dso::microseconds{60 * 60 * 2 * 1'000'000L}}; // 1 hours

constexpr int MIN_INTERPOLATION_PTS = 4;

class SvInterpolator {
private:
  sp3::SatelliteId svid; ///< SV to interpolate
  int num_dpts{0}; ///< num of data points (blocks) available for SV; read from Sp3
  Sp3c *sp3{nullptr};///< Sp3 instance providing data values
  Sp3DataBlock *data{nullptr}; ///< data points/blocks to be collected from the Sp3
  int last_index{0};           ///< last index of data used in the interpolation
  double *txyz{nullptr}, ///< time, x, y and z data arrays used in interpolation
  *workspace{nullptr}; ///< workspace arena (allocate once) used in interpolation

  /// @brief Compute workspace arena size
  /// This function will compute the maximum number of data points to be used
  /// in the interpolation, based on the options available (aka MAX_MICROSEC_)
  /// and the data interval of the Sp3.
  /// We are considering points up to MAX_MICROSEC_ on the right and points
  /// up to MAX_MICROSEC_ on the left.
  /// @return Maximum number of points around a central point, with time tags
  ///         less than MAX_MICROSEC_ apart
  int compute_workspace_size() noexcept;

  /// @brief fill in the data array using an sp3 instance (aka collect SV blocks
  ///        from Sp3)
  int feed_from_sp3() noexcept;

  /// @brief Return the index of the data block in the data array, so that
  ///        bloc[i].t <= t < block[i+1].t
  ///        Try to make an educated guess...
  int index_hunt(const dso::datetime<dso::microseconds> &t) noexcept {
    int start_index = (data[last_index].t <= t) ? last_index : 0;
    auto it = std::lower_bound(
        data+start_index, data + num_dpts, t,
        [](const Sp3DataBlock &block, const dso::datetime<dso::microseconds>& tt) {
        return block.t < tt;
        });
    return static_cast<int>(it-data);
  }

public:
  SvInterpolator(sp3::SatelliteId sid) noexcept : svid(sid){};
 
  /// @brief Constructor from a SatelliteId and an Sp3c instance; this function
  ///        will:
  ///        1. call the feed_from_sp3() function, to allocate needed memory
  ///           and read in the satellite blocks off from the sp3 file
  ///        2. allocate enough workspace (a memory arena) for the t, x, y, and
  ///           z arrays (for later calls to interpolate at) and also the
  ///           c, d arrays that are needed for neville interpolation 
  SvInterpolator(sp3::SatelliteId sid, Sp3c &sp3obj);

  /// @brief Destructor (free memmory)
  ~SvInterpolator() noexcept {
    if (data && num_dpts) delete[] data;
    if (txyz) delete[] txyz;
    if (workspace) delete[] workspace;
    num_dpts = 0;
  }

  int num_data_points() const noexcept { return num_dpts; }

  int interpolate_at(dso::datetime<dso::microseconds> t, double *result,
                     double *error) noexcept;
}; // SvInterpolator

} // dso

#endif
