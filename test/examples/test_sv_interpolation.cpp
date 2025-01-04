#include "sv_interpolate.hpp"
#include <chrono>
#include <cstdio>

using namespace dso;
using dso::sp3::SatelliteId;

int main(int argc, char *argv[]) {
  if (argc < 2 || argc > 3) {
    fprintf(stderr, "Usage: %s <SP3c FILE> [SV]\n", argv[0]);
    return 1;
  }

  Sp3c sp3(argv[1]);
#ifdef DEBUG
  sp3.print_members();
#endif

  SatelliteId sv("L27");
  if (argc==3) sv.set_id(argv[2]);

  if (sp3.num_sats() == 1) {
    printf("\nSp3 file only includes one satellite; perorming interpolation "
           "for %s\n",
           sp3.sattellite_vector()[0].id);
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else if (!sp3.has_sv(sv)) {
    printf("\nSatellite %s not included in sp3 file\n", sv.id);
    return 0;
  }

  SvInterpolator sv_intrp(sv, sp3);
  printf("Fed interpolator with %d data points\n", sv_intrp.num_data_points());

  /*for (int i = 0; i < sv_intrp.num_dpts; i++)
    printf("Sp3 %.7f %.7f %.7f %.7f\n", sv_intrp.data[i].t.as_mjd(),
           sv_intrp.data[i].state[0], sv_intrp.data[i].state[1],
           sv_intrp.data[i].state[2]);*/

  auto start_t = sp3.start_epoch();
  auto stop_t =
      start_t + dso::datetime_interval<nanoseconds>(1, nanoseconds(0));
  auto every_t = dso::datetime_interval<nanoseconds>(
      0, nanoseconds(180 * dso::nanoseconds::sec_factor<long>()));

  // interpolate for the whole time span ...
  /*dso::seconds total_sec{sp3.num_epochs() *
                         sp3.interval().as_underlying_type() /
                         dso::microseconds::sec_factor<long>()};
  int total_days = total_sec.remove_days();
  stop_t = start_t + dso::datetime_interval<nanoseconds>(
                         modified_julian_day(total_days),
                         nanoseconds(total_sec.as_underlying_type() *
                                      dso::nanoseconds::sec_factor<long>()));*/

  auto start_timer = std::chrono::high_resolution_clock::now();
  double xyz[3] = {0}, dxdydz[3] = {0};
  printf("performing interpolation for interval %.5f to %.5f every %ld sec.\n",
         start_t.imjd().as_underlying_type() + start_t.fractional_days().days(),
         stop_t.imjd().as_underlying_type() + stop_t.fractional_days().days(),
         every_t.sec().as_underlying_type() /
             dso::nanoseconds::sec_factor<long>());
  while (start_t < stop_t) {
    sv_intrp.interpolate_at(start_t, xyz, dxdydz);
    // printf("Intrp %.7f %.7f\n", start_t.as_mjd(), *xyz);
    start_t += every_t;
  }
  auto stop_timer = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      stop_timer - start_timer);
  printf("Interpolation took about %ld milliseconds\n", duration.count());

  return 0;
}
