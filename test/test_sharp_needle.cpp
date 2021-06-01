#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("Sharp needle test", "[THIN_NEEDLE]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/sharp_needle1.json";
  const std::string frestart1 =
      HOME +
      "/Workspace/dnn/test/benchmark_results/dump_sharp_needle1_10000.txt";
  const std::string frestart2 =
      HOME +
      "/Workspace/dnn/test/benchmark_results/dump_sharp_needle1_20000.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  // Check the parameters essential to the simulation
  REQUIRE(dnn.get_omega() == 0.05);
  REQUIRE(dnn.get_dt() == 1e-2);
  REQUIRE(dnn.get_d() == Approx(628.3185));
  REQUIRE(dnn.get_grid_Fo() == Approx(0.06283185));
  REQUIRE(dnn.get_h() == 3);

  // Check the Needle
  REQUIRE(dnn.get_needles().size() == 1);
  REQUIRE(dnn.get_needles()[0].x0 == 0);
  REQUIRE(dnn.get_needles()[0].xf == 3);
  REQUIRE(dnn.get_needles()[0].y0 == 750);
  REQUIRE(dnn.get_needles()[0].yf == 750);

  // Check the BCs
  REQUIRE(dnn.get_BC(0).type == FLUX);
  REQUIRE(dnn.get_BC(1).type == FLUX);
  REQUIRE(dnn.get_BC(2).type == FLUX);
  REQUIRE(dnn.get_BC(3).type == FLUX);
  REQUIRE(dnn.get_BC(0).value == Approx(0));
  REQUIRE(dnn.get_BC(1).value == Approx(0));
  REQUIRE(dnn.get_BC(2).value == Approx(0));
  REQUIRE(dnn.get_BC(3).value == Approx(0));

  dnn.set_boundary_values();
  // Check set boundary values
  REQUIRE(dnn.get_grid().at(0, 0) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, 0) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(0, dnn.get_y_size() - 1) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, dnn.get_y_size() - 1) ==
          Approx(0.05));

  dnn.set_needle_values();

  // iteration 0 to 1e4
  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (int i = 0; i < 1e4; ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
    dnn.set_needle_values();
    dnn.grow_needles();
    dnn.set_needle_values();
    dnn.fix_concentration();
  }

  REQUIRE(dnn.get_needles()[0].rad == Approx(0.7261591607));
  REQUIRE(dnn.get_needles()[0].vel == Approx(1.896428001));
  REQUIRE(dnn.get_needles()[0].r == Approx(0.33434505));

  REQUIRE_NOTHROW(odnn.read_restart_file(frestart1));
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }

  // iteration 1e4 to 2e4
  for (int i = 0; i < 1e4; ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
    dnn.set_needle_values();
    dnn.grow_needles();
    dnn.set_needle_values();
    dnn.fix_concentration();
  }

  REQUIRE_NOTHROW(odnn.read_restart_file(frestart2));
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}
