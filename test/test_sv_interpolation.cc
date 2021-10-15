#include "sv_interpolate.hpp"
#include <cstdio>

using namespace dso;
using dso::sp3::SatelliteId;

int main(int argc, char *argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: %s <SP3c FILE>\n", argv[0]);
    return 1;
  }

  Sp3c sp3(argv[1]);
  #ifdef DEBUG
  sp3.print_members();
  #endif

  SatelliteId sv("L27");

  if (sp3.num_sats() == 1) {
    printf("Sp3 file only includes one satellite; perorming interpolation for %s\n", sp3.sattellite_vector()[0].id);
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else if (!sp3.has_sv(sv)) {
    printf("Satellite %s not included in sp3 file\n", sv.id);
    return 0;
  }

  SvInterpolator sv_intrp(sv, sp3);
    printf("Fed interpolator with %d data points\n", sv_intrp.num_dpts);

    for (int i = 0; i < sv_intrp.num_dpts; i++)
      printf("Sp3 %.7f %.7f %.7f %.7f\n", sv_intrp.data[i].t.as_mjd(),
             sv_intrp.data[i].state[0]/1e3, sv_intrp.data[i].state[1]/1e3,
             sv_intrp.data[i].state[2]/1e3);

    auto start_t = sp3.start_epoch();
    auto stop_t = start_t + dso::datetime_interval<microseconds>(
                                modified_julian_day(5), microseconds(0));
    auto every_t = dso::datetime_interval<microseconds>(
        modified_julian_day(0), microseconds(1'000'000 * 30));

    double xyz[3] = {0};
    while (start_t < stop_t) {
      sv_intrp.interpolate_at(start_t, xyz);
      printf("Intrp %.7f %.7f\n", start_t.as_mjd(), *xyz/1e3);
      start_t += every_t;
    }

    return 0;
}