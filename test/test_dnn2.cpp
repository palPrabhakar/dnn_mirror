#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("DNN_INIT_THIN_NEEDLE", "[SUCCESS]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/success_thin_needle.json";
  DNN dnn;
  REQUIRE_NOTHROW(dnn.init(fn));

  SECTION("check read input file values thin needle") {
    REQUIRE(dnn.get_x_size() == 100);
    REQUIRE(dnn.get_y_size() == 100);
    REQUIRE(dnn.get_dx() == Approx(1));
    REQUIRE(dnn.get_dy() == Approx(1));
    REQUIRE(dnn.get_niter() == 100);
    REQUIRE(dnn.get_nprint() == 0);

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
    REQUIRE(dnn.get_needles()[0].y0 == Approx(50));
    REQUIRE(dnn.get_needles()[0].xf == Approx(5));
    REQUIRE(dnn.get_needles()[0].yf == Approx(50));

    REQUIRE(dnn.get_dcl() == Approx(3e-9));
    REQUIRE(dnn.get_sigma() == Approx(0.083));
    REQUIRE(dnn.get_c_inf() == Approx(7));
    REQUIRE(dnn.get_c0() == Approx(7));
    REQUIRE(dnn.get_k() == Approx(0.13));
    REQUIRE(dnn.get_m() == Approx(-6.5));
    REQUIRE(dnn.get_gamma() == Approx(1.96e-7));
    REQUIRE(dnn.get_omega() == Approx(0.05));
    REQUIRE(dnn.get_d() == Approx(628.318));

    // Test for thin needle
    REQUIRE(dnn.get_wintegral() == JINT);

    // Test for transient case
    REQUIRE(dnn.get_dt() == Approx(1e-3));

    // Specifi to thin needle
    REQUIRE(dnn.get_h() == 3);
    REQUIRE(dnn.get_grid_Fo() == Approx(.628318));
  }
}

TEST_CASE("DNN_INIT_THICK_NEEDLE", "[SUCCESS]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/success_thick_needle.json";
  DNN dnn;
  REQUIRE_NOTHROW(dnn.init(fn));

  SECTION("check read input file values thick needle") {
    REQUIRE(dnn.get_x_size() == 100);
    REQUIRE(dnn.get_y_size() == 100);
    REQUIRE(dnn.get_dx() == Approx(1));
    REQUIRE(dnn.get_dy() == Approx(1));
    REQUIRE(dnn.get_niter() == 100);
    REQUIRE(dnn.get_nprint() == 10);

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
    REQUIRE(dnn.get_needles()[0].xf == Approx(7));
    REQUIRE(dnn.get_needles()[0].yf == Approx(0));

    REQUIRE(dnn.get_dcl() == Approx(3e-9));
    REQUIRE(dnn.get_sigma() == Approx(0.083));
    REQUIRE(dnn.get_c_inf() == Approx(7));
    REQUIRE(dnn.get_c0() == Approx(7));
    REQUIRE(dnn.get_k() == Approx(0.13));
    REQUIRE(dnn.get_m() == Approx(-6.5));
    REQUIRE(dnn.get_gamma() == Approx(1.96e-7));
    REQUIRE(dnn.get_omega() == Approx(0.05));
    REQUIRE(dnn.get_d() == Approx(588.68));

    // Test for thick needle
    REQUIRE(dnn.get_wintegral() == NORMAL);
    REQUIRE(dnn.get_A() == 5);
    REQUIRE(dnn.get_B() == 2);
    REQUIRE(dnn.get_C() == 1);

    // Test for transient case
    REQUIRE(dnn.get_dt() == Approx(1e-3));
    REQUIRE(dnn.get_grid_Fo() == Approx(0.58868));
  }
}

TEST_CASE("DNN_INIT_STEADY_STATE", "[SUCCESS]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/success_steady_state.json";
  DNN dnn;
  REQUIRE_NOTHROW(dnn.init(fn));

  SECTION("check read input file values steady state") {
    REQUIRE(dnn.get_x_size() == 250);
    REQUIRE(dnn.get_y_size() == 250);
    REQUIRE(dnn.get_dx() == Approx(10));
    REQUIRE(dnn.get_dy() == Approx(10));
    REQUIRE(dnn.get_niter() == 500);
    REQUIRE(dnn.get_nprint() == 25);

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
    REQUIRE(dnn.get_needles()[0].xf == Approx(50));
    REQUIRE(dnn.get_needles()[0].yf == Approx(0));

    REQUIRE(dnn.get_dcl() == Approx(3e-9));
    REQUIRE(dnn.get_sigma() == Approx(0.083));
    REQUIRE(dnn.get_c_inf() == Approx(7));
    REQUIRE(dnn.get_c0() == Approx(7));
    REQUIRE(dnn.get_k() == Approx(0.13));
    REQUIRE(dnn.get_m() == Approx(-6.5));
    REQUIRE(dnn.get_gamma() == Approx(1.96e-7));
    REQUIRE(dnn.get_omega() == Approx(1));
    REQUIRE(dnn.get_d() == Approx(3e-9));

    // Steady state specific test
    REQUIRE(dnn.get_epsilon() == Approx(1e-9));
    REQUIRE(dnn.get_sor() == Approx(1.5));
  }
}

TEST_CASE("DNN_INIT", "[FAILURE]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";

  SECTION("Invalid input file name") {
    const std::string fn =
        HOME + "/Workspace/dnn/test/test_inputs/thisFileDoesNotExist.json";
    DNN dnn;
    REQUIRE_THROWS_AS(dnn.init(fn), exit_exception);
  }

  SECTION("Invalid file format") {
    const std::string fn =
        HOME + "/Workspace/dnn/test/test_inputs/invalid_input.json";
    DNN dnn;
    REQUIRE_THROWS_AS(dnn.init(fn), exit_exception);
  }

  SECTION("Key missing in input file") {
    const std::string fn1 = HOME + "/Workspace/dnn/test/test_inputs/test2.json";
    DNN dnn;
    REQUIRE_THROWS_AS(dnn.init(fn1), exit_exception);
  }
}

TEST_CASE("DNN_READ_RESTART", "[THIN_NEEDLE]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/restart.json";
  const std::string restart_file =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_restart_test.txt";
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
}

TEST_CASE("DNN_ADD_THIN_NEEDLE", "[SUCCESS]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/success_thin_needle.json";
  DNN dnn;
  dnn.init(fn);
  dnn.set_boundary_values();
  dnn.set_needle_values();
  REQUIRE(dnn.get_grid().at(0, 50) == 0);
  REQUIRE(dnn.get_grid().at(1, 50) == 0);
  REQUIRE(dnn.get_grid().at(2, 50) == 0);
  REQUIRE(dnn.get_grid().at(3, 50) == 0);
  REQUIRE(dnn.get_grid().at(4, 50) == 0);
  REQUIRE(dnn.get_grid().at(5, 50) == 0);
  REQUIRE(dnn.get_grid().at(6, 50) == Approx(0.05));

  // Increase needle length
  dnn.get_needles()[0].xf = 10;
  dnn.set_needle_values();

  REQUIRE(dnn.get_grid().at(0, 50) == 0);
  REQUIRE(dnn.get_grid().at(1, 50) == 0);
  REQUIRE(dnn.get_grid().at(2, 50) == 0);
  REQUIRE(dnn.get_grid().at(3, 50) == 0);
  REQUIRE(dnn.get_grid().at(4, 50) == 0);
  REQUIRE(dnn.get_grid().at(5, 50) == 0);
  REQUIRE(dnn.get_grid().at(6, 50) == 0);
  REQUIRE(dnn.get_grid().at(7, 50) == 0);
  REQUIRE(dnn.get_grid().at(8, 50) == 0);
  REQUIRE(dnn.get_grid().at(9, 50) == 0);
  REQUIRE(dnn.get_grid().at(10, 50) == 0);

  REQUIRE(dnn.get_grid().at(11, 50) == Approx(0.05));
  REQUIRE(dnn.get_grid().at(75, 75) == Approx(0.05));
}
