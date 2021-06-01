#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("fif contour calculation SS H3", "[FIF_SS_H3]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/ss_fif_h3.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_fif_h3.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));
  REQUIRE(dnn.get_h() == 3);

  dnn.set_boundary_values();
  dnn.set_needle_values();

  for (int i = 0; i < 237088; ++i) {
    dnn.run_steady_state();
  }

  dnn.set_needle_values();

  /* REQUIRE(dnn.get_citer() == 237088); */

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }

  REQUIRE(dnn.calculate_flux_intensity_factor_sq_j_integral(
              dnn.get_needles()[0].xf, dnn.get_needles()[0].yf) ==
          Approx(0.9853254975));
}

TEST_CASE("fif contour calculation SS H10", "[FIF_SS_H10]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/ss_fif_h10.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_fif_h10.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));
  REQUIRE(dnn.get_h() == 3);

  dnn.set_boundary_values();
  dnn.set_needle_values();

  for (int i = 0; i < 237088; ++i) {
    dnn.run_steady_state();
  }

  dnn.set_needle_values();

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }

  REQUIRE(dnn.calculate_flux_intensity_factor_sq_j_integral(
              dnn.get_needles()[0].xf, dnn.get_needles()[0].yf) ==
          Approx(0.9990276254));
}

/* TEST_CASE("fif contour calculation transient H3", "[FIF_T_H3]") { */
/*   std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
 */
/*   const std::string fn = HOME +
 * "/Workspace/dnn/test/test_inputs/transient_fif_h3.json"; */
/*   const std::string frestart = HOME +
 * "/Workspace/dnn/test/benchmark_results/dump_transient_fif_h3.txt"; */
/*   DNN dnn, odnn; */
/*   REQUIRE_NOTHROW(dnn.init(fn)); */
/*   REQUIRE_NOTHROW(odnn.init(fn)); */
/*   REQUIRE_NOTHROW(odnn.read_restart_file(frestart)); */
/*   REQUIRE(dnn.get_h() == 3); */
/*   dnn.run(); */
/*   for(std::size_t idx = 0; idx < dnn.get_x_size()*dnn.get_y_size(); ++idx) {
 */
/*     REQUIRE(dnn.get_grid().grid[idx].u ==
 * Approx(odnn.get_grid().grid[idx].u)); */
/*   } */
/*   REQUIRE(dnn.calculate_flux_intensity_factor_sq_j_integral(dnn.get_needles()[0].xf,
 * dnn.get_needles()[0].yf) == Approx(0.9853254975)); */
/* } */

/* TEST_CASE("fif contour calculation transient H10", "[FIF_T_H10]") { */
/*   std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
 */
/*   const std::string fn = HOME +
 * "/Workspace/dnn/test/test_inputs/transient_fif_h10.json"; */
/*   const std::string frestart = HOME +
 * "/Workspace/dnn/test/benchmark_results/dump_transient_fif_h10.txt"; */
/*   DNN dnn, odnn; */
/*   REQUIRE_NOTHROW(dnn.init(fn)); */
/*   REQUIRE_NOTHROW(odnn.init(fn)); */
/*   REQUIRE_NOTHROW(odnn.read_restart_file(frestart)); */
/*   REQUIRE(dnn.get_h() == 10); */
/*   dnn.run(); */
/*   for(std::size_t idx = 0; idx < dnn.get_x_size()*dnn.get_y_size(); ++idx) {
 */
/*     REQUIRE(dnn.get_grid().grid[idx].u ==
 * Approx(odnn.get_grid().grid[idx].u)); */
/*   } */
/*   REQUIRE(dnn.calculate_flux_intensity_factor_sq_j_integral(dnn.get_needles()[0].xf,
 * dnn.get_needles()[0].yf) == Approx(0.9990276253)); */
/* } */
