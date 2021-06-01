#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>

TEST_CASE("2D parabolic needle test", "[2D_PARABOLIC]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/thick_needle2.json";
  const std::string frestart1 =
      HOME +
      "/Workspace/dnn/test/benchmark_results/dump_parabolic_needle2_1000.txt";
  const std::string frestart2 =
      HOME +
      "/Workspace/dnn/test/benchmark_results/dump_parabolic_needle2_2000.txt";
  const std::string frestart3 =
      HOME +
      "/Workspace/dnn/test/benchmark_results/dump_parabolic_needle2_10000.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  // Check the parameters essential to the simulation
  REQUIRE(dnn.get_omega() == 0.15);
  /* REQUIRE(dnn.get_dt() == 1e-3); */
  REQUIRE(dnn.get_d() == Approx(56.84926438));
  /* REQUIRE(dnn.get_grid_Fo() == Approx(0.05684926438)); */
  REQUIRE(dnn.get_A() == 5);
  REQUIRE(dnn.get_B() == 2);
  REQUIRE(dnn.get_C() == 1);

  // Check the Needle
  REQUIRE(dnn.get_needles().size() == 1);
  REQUIRE(dnn.get_needles()[0].x0 == 0);
  REQUIRE(dnn.get_needles()[0].xf == 6);
  REQUIRE(dnn.get_needles()[0].y0 == 700);
  REQUIRE(dnn.get_needles()[0].yf == 700);

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
  REQUIRE(dnn.get_grid().at(0, 0) == Approx(0.15));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, 0) == Approx(0.15));
  REQUIRE(dnn.get_grid().at(0, dnn.get_y_size() - 1) == Approx(0.15));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, dnn.get_y_size() - 1) ==
          Approx(0.15));

  dnn.set_needle_values();

  // iteration 0 to 1e3
  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (int i = 0; i < 1e3; ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
    dnn.set_needle_values();
    dnn.grow_needles();
    dnn.set_needle_values();
    /* dnn.fix_concentration(); */
  }

  /* REQUIRE(dnn.get_needles()[0].rad == Approx(0.4782198828)); */
  /* REQUIRE(dnn.get_needles()[0].vel == Approx(4.372650263)); */
  /* REQUIRE(dnn.get_needles()[0].r == Approx(0.8110392)); */

  REQUIRE_NOTHROW(odnn.read_restart_file(frestart1));
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }

  // iteration 1e3 to 2e3
  for (int i = 0; i < 1e3; ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
    dnn.set_needle_values();
    dnn.grow_needles();
    dnn.set_needle_values();
    /* dnn.fix_concentration(); */
  }

  REQUIRE_NOTHROW(odnn.read_restart_file(frestart2));
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }

  // iteration 2e3 to 10e3
  for (int i = 0; i < 8e3; ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
    dnn.set_needle_values();
    dnn.grow_needles();
    dnn.set_needle_values();
    /* dnn.fix_concentration(); */
  }

  REQUIRE_NOTHROW(odnn.read_restart_file(frestart3));
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}
