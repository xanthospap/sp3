#include "sp3.hpp"
#include <bits/c++config.h>
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
  Sp3DataBlock block;


  if (sp3.num_sats() == 1) {
    printf("Sp3 file only includes one satellite; extracting records for %s\n", sp3.sattellite_vector()[0].id);
    sv.set_id(sp3.sattellite_vector()[0].id);
  } else if (!sp3.has_sv(sv)) {
    printf("Satellite %s not included in sp3 file\n", sv.id);
    return 0;
  }

  // let's try reading the records; note that -1 denotes EOF
  int j;
  std::size_t rec_count = 0;
  do {
    // printf("reading record #%6lu", rec_count);
    j = sp3.get_next_data_block(sv, block);
    if (j > 0) {
      printf("Something went wrong ....status = %3d\n", j);
      return 1;
    } else if (j == -1) {
      printf("EOF encountered; Sp3 file read through!\n");
    }
    bool position_ok = !block.flag.is_set(Sp3Event::bad_abscent_position);
    // bool velocity_ok = !block.flag.is_set(Sp3Event::bad_abscent_velocity);
    if (position_ok) printf("%15.6f %15.7f %15.7f %15.7f\n", block.t.as_mjd(), block.state[0], block.state[1], block.state[2]);
    ++rec_count;
  } while (!j);

  printf("Num of records read: %6lu\n", rec_count);

  return 0;
}
