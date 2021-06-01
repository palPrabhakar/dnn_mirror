#include "../src/dnn.hpp"
#include "catch.hpp"
#include <cstddef>
#include <iostream>
#include <stdexcept>

TEST_CASE("CONST_TEMPERATURE_BC", "[T1]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test1.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test1.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX1", "[T2]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test2.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test2.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX2", "[T3]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test3.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test3.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX3", "[T4]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test4.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test4.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX4", "[T5]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test5.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test5.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX5", "[T6]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test6.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test6.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX6", "[T7]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test7.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test7.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX7", "[T8]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test8.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test8.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX8", "[T9]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test9.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test9.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_TEMPERATURE_NO_FLUX9", "[T10]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test10.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test10.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}

TEST_CASE("CONST_FLUX_ADD", "[T11]") {
  std::string const HOME = std::getenv("HOME") ? std::getenv("HOME") : ".";
  const std::string fn =
      HOME + "/Workspace/dnn/test/test_inputs/transient_test11.json";
  const std::string frestart =
      HOME + "/Workspace/dnn/test/benchmark_results/dump_transient_test11.txt";
  DNN dnn(fn), odnn(fn);
  REQUIRE_NOTHROW(dnn.read_input_file());
  REQUIRE_NOTHROW(odnn.read_input_file());
  REQUIRE_NOTHROW(dnn.init());
  REQUIRE_NOTHROW(odnn.init());
  REQUIRE_NOTHROW(odnn.read_restart_file(frestart));

  std::vector<double> temp(dnn.get_x_size() * dnn.get_y_size());
  for (size_t i = 0; i < dnn.get_niter(); ++i) {
    dnn.run_transient(temp);
    dnn.get_grid().swap(temp);
  }

  for (std::size_t idx = 0; idx < dnn.get_x_size() * dnn.get_y_size(); ++idx) {
    CHECK(dnn.get_grid().at(idx) == Approx(odnn.get_grid().at(idx)));
  }
}
