#include "sv_interpolate.hpp"
#include "datetime/calendar.hpp"

int dso::SvInterpolator::compute_workspace_size() noexcept {
  dso::nanoseconds lr_intrvl =
      dso::cast_to<dso::milliseconds, dso::nanoseconds>(max_millisec);
  int one_side_pts =
      lr_intrvl.as_underlying_type() / sp3->interval().as_underlying_type();
  ++one_side_pts;
  return one_side_pts * 2 + 1;
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

dso::SvInterpolator::SvInterpolator(sp3::SatelliteId sid, Sp3c &sp3obj,
                                    dso::milliseconds max_allowed_millisec)
    : svid(sid), sp3(&sp3obj), max_millisec(max_allowed_millisec) {
  if (int error = feed_from_sp3(); error) {
    throw std::runtime_error("[ERROR] Failed creating SvInterpolator "
                             "instance from Sp3 Error Code: " +
                             std::to_string(error));
  }
  int workspace_size = compute_workspace_size();
  txyz = new double[workspace_size * 4];
  workspace = new double[workspace_size * 6];
}

/*
int dso::SvInterpolator::interpolate_at(dso::datetime<dso::nanoseconds> t,
                                        double *result,
                                        double *error) noexcept {

  int index = index_hunt(t);
  last_index = index;
#ifdef DEBUG
  if (index < 0 || index > num_dpts - 1) {
    fprintf(stderr, "[DEBUG] Invalid index! hunt returned index=%d\n", index);
  }
#endif

  // the max allowed interval as datetime_interval to symplify computations
  dso::datetime_interval<dso::nanoseconds> max_t{
      dso::modified_julian_day{0},
      dso::cast_to<dso::milliseconds, dso::nanoseconds>(max_millisec)};

  // start point on the left ....
  int start = index;
  while (start > 0 && t - data[start].t < max_t)
    --start;
  if (index - start < min_dpts_on_each_side) {
    fprintf(
        stderr,
        "[ERROR] Cannot interpolate due to too few data points (on the left) "
        "(traceback: %s)\n",
        __func__);
#ifdef DEBUG
    fprintf(stderr,
            "[ERROR] Index=%d (start=%d), index mjd=%.9f, data start=%.9f\n",
            index, start, data[index].t.as_mjd(), sp3->start_epoch().as_mjd());
    dso::seconds masec =
        dso::cast_to<dso::milliseconds, dso::seconds>(max_millisec);
    dso::seconds sisec =
        dso::cast_to<dso::nanoseconds, dso::seconds>(this->sp3->interval());
    fprintf(stderr,
            "[ERROR] Interpolator allowed to go back %ld seconds, Sp3 file has "
            "positions every %ld seconds\n",
            masec.as_underlying_type(), sisec.as_underlying_type());
#endif
    return 1;
  }

  // end point at the right (inclusive)
  int stop = index;
  while (stop < num_dpts - 1 && data[stop].t - t < max_t)
    ++stop;
  if (stop - index < min_dpts_on_each_side) {
    fprintf(
        stderr,
        "[ERROR] Cannot interpolate due to too few data points (on the right) "
        "(traceback: %s)\n",
        __func__);
    return 1;
  }

  // number of data points to be used in interpolation
  int size = stop - start + 1;

  // seperate the workspace arena to arrays of x, y, z and time
  int wsz = compute_workspace_size();
  double *__restrict__ td = txyz + 0 * wsz;
  double *__restrict__ xd = txyz + 1 * wsz;
  double *__restrict__ yd = txyz + 2 * wsz;
  double *__restrict__ zd = txyz + 3 * wsz;

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
*/

int dso::SvInterpolator::interpolate_at(dso::datetime<dso::nanoseconds> t,
                                        double *pos, double *erpos,
                                        double *vel, double *ervel) noexcept {

  int index = index_hunt(t);
  last_index = index;
#ifdef DEBUG
  if (index < 0 || index > num_dpts - 1) {
    fprintf(stderr, "[DEBUG] Invalid index! hunt returned index=%d\n", index);
  }
#endif

  // the max allowed interval as datetime_interval to simplify computations
  dso::datetime_interval<dso::nanoseconds> max_t{
      0, dso::cast_to<dso::milliseconds, dso::nanoseconds>(max_millisec)};

  // start point on the left ....
  int start = index;
  while (start > 0 && t - data[start].t < max_t)
    --start;
  if (index - start < min_dpts_on_each_side) {
    fprintf(
        stderr,
        "[ERROR] Cannot interpolate due to too few data points (on the left) "
        "(traceback: %s)\n",
        __func__);
//#ifdef DEBUG
//    fprintf(stderr,
//            "[ERROR] Index=%d (start=%d), index mjd=%.9f, data start=%.9f\n",
//            index, start, data[index].t.as_mjd(), sp3->start_epoch().as_mjd());
//    dso::seconds masec =
//        dso::cast_to<dso::milliseconds, dso::seconds>(max_millisec);
//    dso::seconds sisec =
//        dso::cast_to<dso::nanoseconds, dso::seconds>(this->sp3->interval());
//    fprintf(stderr,
//            "[ERROR] Interpolator allowed to go back %ld seconds, Sp3 file has "
//            "positions every %ld seconds\n",
//            masec.as_underlying_type(), sisec.as_underlying_type());
//#endif
    return 1;
  }

  // end point at the right (inclusive)
  int stop = index;
  while (stop < num_dpts - 1 && data[stop].t - t < max_t)
    ++stop;
  if (stop - index < min_dpts_on_each_side) {
    fprintf(
        stderr,
        "[ERROR] Cannot interpolate due to too few data points (on the right) "
        "(traceback: %s)\n",
        __func__);
    return 1;
  }

  // number of data points to be used in interpolation
  int size = stop - start + 1;

  // seperate the workspace arena to arrays of x, y, z and time
  int wsz = compute_workspace_size();
  double *__restrict__ td = txyz + 0 * wsz;
  double *__restrict__ xd = txyz + 1 * wsz;
  double *__restrict__ yd = txyz + 2 * wsz;
  double *__restrict__ zd = txyz + 3 * wsz;

  // fill in arays for each component
  auto start_t = sp3->start_epoch();
  for (int i = 0; i < size; i++) {
    // td[i] = data[start + i].t.delta_date(start_t).as_mjd();
    td[i] =
        data[start + i].t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(
            start_t).seconds();
    xd[i] = data[start + i].state[0];
    yd[i] = data[start + i].state[1];
    zd[i] = data[start + i].state[2];
  }

  // point to interpolate at (as fractional days)
  double tx = // t.delta_date(start_t).as_mjd();
    t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(start_t).seconds();

  // perform the interpolation for all components
  if (sp3::neville_interpolation3(tx, pos, erpos, td, xd, yd, zd, size, size,
                                  0, workspace)) {
    fprintf(stderr, "[ERROR] Neville algorithm failed (traceback: %s)\n",
            __func__);
    return 5;
  }
  
  // interpolate velocity if needed
  if (vel && ervel) {
    // fill in arays for each component
    for (int i = 0; i < size; i++) {
      td[i] =
          data[start + i]
              .t.diff<dso::DateTimeDifferenceType::FractionalSeconds>(start_t).seconds();
      xd[i] = data[start + i].state[4];
      yd[i] = data[start + i].state[5];
      zd[i] = data[start + i].state[6];
    }

    // perform the interpolation for all components
    if (sp3::neville_interpolation3(tx, vel, ervel, td, xd, yd, zd, size, size,
                                    0, workspace)) {
      fprintf(stderr, "[ERROR] Neville algorithm failed (traceback: %s)\n",
              __func__);
      return 6;
    }
  }

  return 0;
}
