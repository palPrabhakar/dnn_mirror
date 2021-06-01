#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("CONST_TEMPERATURE_BC", "[SS1]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test1.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test1.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));
  for (int i = 0; i < 280; ++i) {
    dnn.run_steady_state();
  }
  /* REQUIRE(dnn.get_citer() == 280); */
  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX", "[SS2]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test2.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test2.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 491; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX2", "[SS3]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test3.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test3.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 317; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX3", "[SS4]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test4.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test4.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 317; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX4", "[SS5]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test5.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test5.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 355; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX5", "[SS6]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test6.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test6.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 355; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_FLUX", "[SS7]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test7.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test7.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 275; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_FLUX2", "[SS8]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test8.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test8.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 275; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_FLUX3", "[SS9]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn = HOME + "/Workspace/dnn/test/test_inputs/ss_test9.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test9.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 777; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_FLUX4", "[SS10]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/ss_test10.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_ss_test10.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  for (int i = 0; i < 495; ++i) {
    dnn.run_steady_state();
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    REQUIRE(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}
