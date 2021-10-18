#include "sv_interpolate.hpp"

int dso::SvInterpolator::compute_workspace_size() noexcept {
  dso::milliseconds max_ms{MAX_MICROSEC_.to_sec_type<dso::milliseconds>()};
  dso::milliseconds data_every_ms{
      dso::cast_to<dso::microseconds, dso::milliseconds>(sp3->interval())};
  int lr_points =
      max_ms.as_underlying_type() / data_every_ms.as_underlying_type() + 1;
  return lr_points * 2 + 1;
}

int dso::SvInterpolator::feed_from_sp3() noexcept {
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
    data = nullptr;
  }

  // set num_dpts to 0 in case of error; else num_dpts is the actual number of
  // data blocks read in from the Sp3 instance, excluding the ones with bad
  // flags
  num_dpts = (error > 0) ? 0 : idx;

  return error > 0;
}

dso::SvInterpolator::SvInterpolator(sp3::SatelliteId sid, Sp3c &sp3obj)
    : svid(sid), sp3(&sp3obj) {
  if (int error = feed_from_sp3(); error) {
    throw std::runtime_error("[ERROR] Failed creating SvInterpolator "
                             "instance from Sp3 Error Code: " +
                             std::to_string(error));
  }
  int workspace_size = compute_workspace_size();
  txyz = new double[workspace_size * 4];
  workspace = new double[workspace_size * 6];
}

int dso::SvInterpolator::interpolate_at(dso::datetime<dso::microseconds> t,
                                        double *result,
                                        double *error) noexcept {

  int index = index_hunt(t);
#ifdef DEBUG
  if (index < 0 || index > num_dpts - 1) {
    fprintf(stderr, "[DEBUG] Invalid index! hunt returned index=%d\n", index);
  }
#endif

  // start point on the left ....
  int start = index;
  while (start > 0 && t - data[start].t < MAX_MICROSEC_)
    --start;

  // end point at the right (inclusive)
  int stop = index;
  while (stop < num_dpts - 1 && data[stop].t - t < MAX_MICROSEC_)
    ++stop;

  // number of data points to be used in interpolation
  int size = stop - start + 1;

  if (size < MIN_INTERPOLATION_PTS) {
    fprintf(stderr,
            "[ERROR] Cannot interpolate due to too few data points "
            "(traceback: %s)\n",
            __func__);
    return 1;
  }

  // seperate the workspace arena to arrays of x, y, z and time
  int wsz = compute_workspace_size();
  double *td = txyz + 0 * wsz;
  double *xd = txyz + 1 * wsz;
  double *yd = txyz + 2 * wsz;
  double *zd = txyz + 3 * wsz;

  // fill in arays for each component
  auto start_t = sp3->start_epoch();
  for (int i = 0; i < size; i++) {
    td[i] = data[start + i].t.delta_date(start_t).as_mjd();
    xd[i] = data[start + i].state[0];
    yd[i] = data[start + i].state[1];
    zd[i] = data[start + i].state[2];
  }

  // point to interpolate at (as fractional days)
  double tx = t.delta_date(start_t).as_mjd();

  // perform the interpolation for all components
  if (sp3::neville_interpolation3(tx, result, error, td, xd, yd, zd, size, size,
                                  0, workspace)) {
    fprintf(stderr, "[ERROR] Neville algorithm failed (traceback: %s)\n",
            __func__);
    return 5;
  }

  return 0;
}
