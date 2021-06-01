#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <random>
#include <stdexcept>

#define REQUIRE CHECK

// Random number generator instance
std::random_device dev;
std::mt19937 rng(dev());

TEST_CASE("Thick needle hand calculation", "[THICK_NEEDLE]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/thick_needle1.json";

  DNN dnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());

  std::uniform_int_distribution<std::mt19937::result_type> distx(
      1, dnn.get_x_size() - 2); // distribution in range [1, x_size()-2]
  std::uniform_int_distribution<std::mt19937::result_type> disty(
      1, dnn.get_y_size() - 2); // distribution in range [1, y_size()-2]

  /* SECTION("Check simulation parameters") { */
  REQUIRE(dnn.get_omega() == 0.05);
  REQUIRE(dnn.get_dt() == Approx(1e-4));
  REQUIRE(dnn.get_d() == Approx(588.683856));
  REQUIRE(dnn.get_grid_Fo() == Approx(0.0588683856));
  REQUIRE(dnn.get_wintegral() == 2);
  REQUIRE(dnn.get_A() == 3);
  REQUIRE(dnn.get_B() == 2);
  REQUIRE(dnn.get_C() == 1);
  /* } */

  /* SECTION("Check the needle") { */
  REQUIRE(dnn.get_needles().size() == 1);
  REQUIRE(dnn.get_needles()[0].x0 == 0);
  REQUIRE(dnn.get_needles()[0].xf == 4);
  REQUIRE(dnn.get_needles()[0].y0 == 25);
  REQUIRE(dnn.get_needles()[0].yf == 25);
  /* } */

  /* SECTION("Check the BCs") { */
  REQUIRE(dnn.get_BC(0).type == FLUX);
  REQUIRE(dnn.get_BC(1).type == FLUX);
  REQUIRE(dnn.get_BC(2).type == FLUX);
  REQUIRE(dnn.get_BC(3).type == FLUX);
  REQUIRE(dnn.get_BC(0).value == Approx(0));
  REQUIRE(dnn.get_BC(1).value == Approx(0));
  REQUIRE(dnn.get_BC(2).value == Approx(0));
  REQUIRE(dnn.get_BC(3).value == Approx(0));
  /* } */

  /* SECTION("Check inner domain") { */
  REQUIRE(dnn.get_grid().at(distx(rng), disty(rng)) == Approx(0.05));
  /* } */

  dnn.set_boundary_values();
  /* SECTION("Check set boundary values") { */
  REQUIRE(dnn.get_grid().at(0, 0) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, 0) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(0, dnn.get_y_size() - 1) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, dnn.get_y_size() - 1) ==
          Approx(0.05));
  REQUIRE(dnn.get_grid().at(distx(rng), 0) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(0, disty(rng)) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(distx(rng), dnn.get_y_size() - 1) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(dnn.get_x_size() - 1, disty(rng)) == Approx(0.05));
  /* } */

  dnn.set_needle_values();
  /* SECTION("Check needle values") { */
  REQUIRE(dnn.get_grid().at(0, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(5, 25) != Approx(0.0));
  REQUIRE(dnn.get_grid().at(6, 25) != Approx(0.0));
  /* } */

  // First iteration
  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  dnn.run_transient(temp);
  REQUIRE(temp[dnn.get_grid().index(5, 25)] == Approx(0.046858));
  dnn.get_grid().swap(temp);
  dnn.set_needle_values();
  /* SECTION("Check values after first transient run") { */
  REQUIRE(dnn.get_grid().at(0, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(5, 25) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(0, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(1, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(4, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(0, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(1, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(4, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(5, 26) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 24) == Approx(0.05));
  /* } */

  /* SECTION("Grow Needle") { */
  REQUIRE(dnn.calculate_line_integral(1, 6, 24, 26, 25, 25) ==
          Approx(0.0251328));
  REQUIRE(dnn.calculate_surface_integral(1, 6, 24, 26) == Approx(0.032854));
  REQUIRE_NOTHROW(dnn.grow_needles());
  /* REQUIRE(dnn.grow_needles() == Approx(1.315953e-5)); */
  REQUIRE(dnn.get_needles()[0].vel == Approx(4.761606));
  REQUIRE(dnn.get_needles()[0].rad == Approx(0.4582725));
  REQUIRE(dnn.get_needles()[0].r == Approx(0.00047616));
  /* } */

  dnn.set_needle_values();
  /* SECTION("Check thick needle growth") { */
  REQUIRE(dnn.get_grid().at(0, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(5, 25) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(0, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(4, 26) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(0, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(4, 24) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(5, 26) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 24) == Approx(0.05));
  /* } */

  dnn.fix_concentration();
  /* SECTION("Check needle tip concentration after fix") { */
  REQUIRE(dnn.get_grid().at(5, 25) == Approx(0.03749255));
  /* } */

  /*
   * SECOND ITERATION
   */

  dnn.run_transient(temp);
  dnn.get_grid().swap(temp);
  dnn.set_needle_values();
  // Check values after first transient run
  REQUIRE(dnn.get_grid().at(0, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(5, 25) == Approx(0.0374945));
  REQUIRE(dnn.get_grid().at(6, 25) == Approx(0.049214));
  REQUIRE(dnn.get_grid().at(7, 25) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 26) == Approx(0.041167));
  REQUIRE(dnn.get_grid().at(4, 26) == Approx(0.044308991));
  REQUIRE(dnn.get_grid().at(5, 26) == Approx(0.04901674));
  REQUIRE(dnn.get_grid().at(6, 26) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 26) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 24) == Approx(0.041167));
  REQUIRE(dnn.get_grid().at(4, 24) == Approx(0.044308991));
  REQUIRE(dnn.get_grid().at(5, 24) == Approx(0.04901674));
  REQUIRE(dnn.get_grid().at(6, 24) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 24) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 27) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(1, 27) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 27) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 27) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(4, 27) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(5, 27) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 27) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 27) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 23) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(1, 23) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 23) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 23) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(4, 23) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(5, 23) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 23) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 23) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(1, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(2, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(3, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(4, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 28) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(1, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(2, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(3, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(4, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 22) == Approx(0.05));

  // Grow Needle
  REQUIRE(dnn.calculate_line_integral(1, 6, 23, 27, 24, 26) ==
          Approx(0.014141804742));
  REQUIRE(dnn.calculate_surface_integral(1, 6, 23, 27) == Approx(0.095975));
  REQUIRE_NOTHROW(dnn.grow_needles());
  /* REQUIRE(dnn.grow_needles() == Approx(4.6061e-6)); */
  REQUIRE(dnn.get_needles()[0].vel == Approx(2.3649032387));
  REQUIRE(dnn.get_needles()[0].rad == Approx(0.6502693914));
  /* REQUIRE(dnn.get_needles()[0].r == Approx(0.00047616)); */

  dnn.set_needle_values();
  dnn.fix_concentration();
  // Check new needle values
  REQUIRE(dnn.get_grid().at(0, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 25) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(5, 25) == Approx(0.03680099217));
  REQUIRE(dnn.get_grid().at(6, 25) == Approx(0.049214));
  REQUIRE(dnn.get_grid().at(7, 25) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 26) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 26) == Approx(0.044308991));
  REQUIRE(dnn.get_grid().at(5, 26) == Approx(0.04901674));
  REQUIRE(dnn.get_grid().at(6, 26) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 26) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(2, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(3, 24) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(4, 24) == Approx(0.044308991));
  REQUIRE(dnn.get_grid().at(5, 24) == Approx(0.04901674));
  REQUIRE(dnn.get_grid().at(6, 24) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 24) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 27) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 27) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 27) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 27) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(4, 27) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(5, 27) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 27) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 27) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 23) == Approx(0.0));
  REQUIRE(dnn.get_grid().at(1, 23) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(2, 23) == Approx(0.046858));
  REQUIRE(dnn.get_grid().at(3, 23) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(4, 23) == Approx(0.0498026079));
  REQUIRE(dnn.get_grid().at(5, 23) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 23) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 23) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(1, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(2, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(3, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(4, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 28) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 28) == Approx(0.05));

  REQUIRE(dnn.get_grid().at(0, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(1, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(2, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(3, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(4, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(5, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(6, 22) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(7, 22) == Approx(0.05));
}
