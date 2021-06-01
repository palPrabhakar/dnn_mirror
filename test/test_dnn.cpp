#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("DNN_READ_INPUT", "[READ_INPUT]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/test.json";
  const std::string fn2 = HOME + "/Workspace/dnn/test/test_inputs/test2.json";
  const std::string fn3 = HOME + "/Workspace/dnn/test/test_inputs/test3.json";
  const std::string fn4 = HOME + "/Workspace/dnn/test/test_inputs/test4.json";
  DNN dnn;

  SECTION("Read Input file test failure") {
    REQUIRE_THROWS_AS(dnn.init("foo.json"), exit_exception);
    REQUIRE_THROWS_AS(dnn.init(fn2), exit_exception);
    REQUIRE_THROWS_AS(dnn.init(fn3), exit_exception);
    REQUIRE_THROWS_AS(dnn.init(fn4), exit_exception);
    /* REQUIRE_THROWS_WITH(dnn.read_input_file(fn2),
     * Catch::Matchers::Contains("Key not found")); */
  }

  dnn.init(fn);

  SECTION("Check file name.") { REQUIRE(dnn.get_input_file_name() == fn); }

  SECTION("Read Input file success") {
    /* REQUIRE_NOTHROW(dnn.read_input_file(fn)); */
    REQUIRE(dnn.get_x_size() == 100);
    REQUIRE(dnn.get_y_size() == 100);
    REQUIRE(dnn.get_dx() == Approx(1));
    REQUIRE(dnn.get_dy() == Approx(1));
    /* REQUIRE(dnn.get_coeff() == Approx(1.0)); */
    REQUIRE(dnn.get_niter() == 100);
    REQUIRE(dnn.get_h() == 3);
    REQUIRE(dnn.get_BC(BC_LEFT_WALL).type == "FLUX");
    REQUIRE(dnn.get_BC(BC_LEFT_WALL).value == Approx(0));
    REQUIRE(dnn.get_BC(BC_BOTTOM_WALL).type == "CONST");
    REQUIRE(dnn.get_BC(BC_BOTTOM_WALL).value == Approx(0));
    REQUIRE(dnn.get_BC(BC_RIGHT_WALL).type == "FLUX");
    REQUIRE(dnn.get_BC(BC_RIGHT_WALL).value == Approx(0.3236));
    REQUIRE(dnn.get_BC(BC_TOP_WALL).type == "CONST");
    REQUIRE(dnn.get_BC(BC_TOP_WALL).value == Approx(1));
    REQUIRE(dnn.get_needles().size() == 1);
    REQUIRE(dnn.get_needles()[0].x0 == Approx(0));
    REQUIRE(dnn.get_needles()[0].y0 == Approx(0));
    REQUIRE(dnn.get_needles()[0].xf == Approx(5));
    REQUIRE(dnn.get_needles()[0].yf == Approx(0));
    REQUIRE(dnn.get_epsilon() == Approx(1e-5));
    REQUIRE(dnn.get_dcl() == Approx(3e-9));
    REQUIRE(dnn.get_sigma() == Approx(0.083));
    REQUIRE(dnn.get_c_inf() == Approx(7));
    REQUIRE(dnn.get_c0() == Approx(7));
    REQUIRE(dnn.get_k() == Approx(0.13));
    REQUIRE(dnn.get_m() == Approx(-6.5));
    REQUIRE(dnn.get_gamma() == Approx(1.96e-7));
  }

  SECTION("Dump Output") {
    REQUIRE_NOTHROW(dnn.dump_output());
    /* REQUIRE_NOTHROW(dnn2.dump_output()); */
  }
}

TEST_CASE("DNN_LAPLACIAN_DIRICHLET") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/data.json";
  DNN dnn;
  dnn.init(fn);
  dnn.set_boundary_values();
  dnn.dump_output();
  for (int i = 0; i < 100; ++i) {
    dnn.run_laplacian();
  }
  dnn.dump_output();

  CHECK(dnn.get_sor() == 1.5);

  REQUIRE(dnn.get_grid().at(0, 0) == Approx(37.5));
  REQUIRE(dnn.get_grid().at(4, 0) == Approx(25));
  REQUIRE(dnn.get_grid().at(0, 4) == Approx(87.5));
  REQUIRE(dnn.get_grid().at(4, 4) == Approx(75));
  CHECK(dnn.get_grid().at(1, 1) == Approx(43).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 1) == Approx(33.30).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 1) == Approx(33.89).epsilon(0.01));
  CHECK(dnn.get_grid().at(1, 2) == Approx(63.21).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 2) == Approx(56.11).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 2) == Approx(52.34).epsilon(0.01));
  CHECK(dnn.get_grid().at(1, 3) == Approx(78.59).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 3) == Approx(76.06).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 3) == Approx(69.71).epsilon(0.01));
}

TEST_CASE("DNN_LAPLACIAN_NEUMANN") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/data2.json";
  DNN dnn;
  dnn.init(fn);
  dnn.dump_output();
  for (int i = 0; i < 100; ++i) {
    dnn.run_laplacian();
  }
  dnn.dump_output();

  CHECK(dnn.get_sor() == 1);

  CHECK(dnn.get_grid().at(1, 0) == Approx(70.06).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 0) == Approx(63.25).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 0) == Approx(50.85).epsilon(0.01));
  CHECK(dnn.get_grid().at(1, 1) == Approx(71.61).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 1) == Approx(66.03).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 1) == Approx(57.55).epsilon(0.01));
  CHECK(dnn.get_grid().at(1, 2) == Approx(76.0).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 2) == Approx(71.8).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 2) == Approx(63.28).epsilon(0.01));
  CHECK(dnn.get_grid().at(1, 3) == Approx(83.4).epsilon(0.01));
  CHECK(dnn.get_grid().at(2, 3) == Approx(82.6).epsilon(0.01));
  CHECK(dnn.get_grid().at(3, 3) == Approx(74.3).epsilon(0.01));
}

TEST_CASE("DNN_READ_RESTART") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/restart.json";
  const std::string restart_file =
      HOME + "/Workspace/dnn/dump/dump_ss_1m_iterm/dump_990000.txt";
  DNN dnn;
  dnn.init(fn);
  Needle needle;
  needle.x0 = 0;
  needle.xf = 160;
  needle.y0 = 30;
  needle.yf = 30;
  REQUIRE_NOTHROW(dnn.read_restart_file(restart_file));
  REQUIRE(dnn.calculate_flux_intensity_factor_sq_j_integral(160, 30) ==
          Approx(1.0).epsilon(0.1));
  REQUIRE(dnn.calculate_flux_intensity_factor(needle, 3, 2, 1) == Approx(1.0));
}

TEST_CASE("DNN_TRANSIENT_DIRICHLET") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/data3.json";
  DNN dnn;
  dnn.init(fn);
  dnn.set_boundary_values();
  for (int i = 0; i < 10000; ++i) {
    dnn.run_transient();
  }
  dnn.dump_output();

  /* CHECK(dnn.get_sor() == 1.5); */

  /* CHECK(dnn.get_grid().at(1, 1) == Approx(43).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(2, 1) == Approx(33.30).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(3, 1) == Approx(33.89).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(1, 2) == Approx(63.21).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(2, 2) == Approx(56.11).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(3, 2) == Approx(52.34).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(1, 3) == Approx(78.59).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(2, 3) == Approx(76.06).epsilon(0.01)); */
  /* CHECK(dnn.get_grid().at(3, 3) == Approx(69.71).epsilon(0.01)); */
}

TEST_CASE("DNN_ADD_NEEDLE") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/data3.json";
  DNN dnn;
  dnn.init(fn);
  dnn.set_boundary_values();
  dnn.set_needle_values();
  REQUIRE(dnn.get_grid().at(0, 0) == 0);
  REQUIRE(dnn.get_grid().at(1, 0) == 0);
  REQUIRE(dnn.get_grid().at(2, 0) == 0);
  REQUIRE(dnn.get_grid().at(3, 0) == 0);
  REQUIRE(dnn.get_grid().at(4, 0) == 200);
  dnn.get_needles()[0].xf = 10;
  dnn.set_needle_values();
  REQUIRE(dnn.get_grid().at(4, 0) == 0);
  REQUIRE(dnn.get_grid().at(5, 0) == 0);
  REQUIRE(dnn.get_grid().at(6, 0) == 0);
  REQUIRE(dnn.get_grid().at(7, 0) == 0);
  REQUIRE(dnn.get_grid().at(8, 0) == 0);
  REQUIRE(dnn.get_grid().at(9, 0) == 0);
  REQUIRE(dnn.get_grid().at(10, 0) == 0);
  REQUIRE(dnn.get_grid().at(0, 1) == 0);
  REQUIRE(dnn.get_grid().at(0, 2) == 0);
  REQUIRE(dnn.get_grid().at(0, 3) == 0);
  dnn.get_needles()[1].yf = 10;
  dnn.set_needle_values();
  REQUIRE(dnn.get_grid().at(0, 4) == 0);
  REQUIRE(dnn.get_grid().at(0, 5) == 0);
  REQUIRE(dnn.get_grid().at(0, 6) == 0);
  REQUIRE(dnn.get_grid().at(0, 7) == 0);
  REQUIRE(dnn.get_grid().at(0, 8) == 0);
  REQUIRE(dnn.get_grid().at(0, 9) == 0);
  REQUIRE(dnn.get_grid().at(0, 10) == 0);
  dnn.dump_output();
}

TEST_CASE("DNN_GROW_NEEDLE") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/restart.json";
  const std::string restart_file =
      HOME + "/Workspace/dnn/dump_ss_1m_iterm/dump_990000.txt";
  DNN dnn;
  dnn.init(fn);
  REQUIRE_NOTHROW(dnn.read_restart_file(restart_file));
  dnn.grow_needles();
  REQUIRE(dnn.get_needles()[0].x0 == 0);
  REQUIRE(dnn.get_needles()[0].xf == 161);
  REQUIRE(dnn.get_needles()[0].y0 == 30);
  REQUIRE(dnn.get_needles()[0].yf == 30);
  REQUIRE(dnn.get_needles()[0].r == Approx(0.422).epsilon(0.01));
}
