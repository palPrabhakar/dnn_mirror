#include "../src/grid.hpp"
#include "catch.hpp"
#include <cstddef>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <vector>

TEST_CASE("CREATE_GRID", "[GRID]") {
  Grid grid;
  Grid grid1(4, 5, 0);

  SECTION("Check x size.") {
    REQUIRE(grid.size() == 0);
    REQUIRE(grid1.size() == 4 * 5);
  }

  SECTION("Set and Get values.") {
    REQUIRE_NOTHROW(grid1.set_value(5, 5, 5));
    REQUIRE(grid1.at(5, 5) == 5);
  }
}

TEST_CASE("SHIFT_GRID", "[GRID]") {
  Grid grid1(3, 3, 0);
  /* std::vector<double> d = {1, 2, 3, 4, 5, 6, 7, 8, 9}; */
  std::vector<Point> d1 = {Point(1, false), Point(2, false), Point(3, false),
                           Point(4, false), Point(5, false), Point(6, false),
                           Point(7, false), Point(8, false), Point(9, false)};

  grid1.set_values(d1);

  SECTION("Check Assignment") {
    REQUIRE(grid1.at(0, 0) == 1);
    REQUIRE(grid1.at(1, 0) == 2);
    REQUIRE(grid1.at(2, 0) == 3);
    REQUIRE(grid1.at(0, 1) == 4);
    REQUIRE(grid1.at(1, 1) == 5);
    REQUIRE(grid1.at(2, 1) == 6);
    REQUIRE(grid1.at(0, 2) == 7);
    REQUIRE(grid1.at(1, 2) == 8);
    REQUIRE(grid1.at(2, 2) == 9);
  }

  grid1.shift_values();

  SECTION("Check Shift") {
    CHECK(grid1.at(0, 0) == 2);
    CHECK(grid1.at(1, 0) == 3);
    CHECK(grid1.at(2, 0) == 3);
    CHECK(grid1.at(0, 1) == 5);
    CHECK(grid1.at(1, 1) == 6);
    CHECK(grid1.at(2, 1) == 6);
    CHECK(grid1.at(0, 2) == 8);
    CHECK(grid1.at(1, 2) == 9);
    CHECK(grid1.at(2, 2) == 9);
  }
}

TEST_CASE("SHIFT_GRID_2", "[GRID]") {
  Grid grid1(5, 4, 0);
  /* std::vector<double> d = {1, 2, 3, 4, 5, 6, 7, 8, 9}; */
  /* std::vector<double> d1 = {1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 1,
   * 2, 3, 4, 5}; */
  std::vector<Point> d1 = {
      Point(1, false), Point(2, false), Point(3, false), Point(4, false),
      Point(5, false), Point(1, false), Point(2, false), Point(3, false),
      Point(4, false), Point(5, false), Point(1, false), Point(2, false),
      Point(3, false), Point(4, false), Point(5, false), Point(1, false),
      Point(2, false), Point(3, false), Point(4, false), Point(5, false)};

  grid1.set_values(d1);

  SECTION("Check Assignments") {
    CHECK(grid1.at(0, 0) == 1);
    CHECK(grid1.at(1, 0) == 2);
    CHECK(grid1.at(2, 0) == 3);
    CHECK(grid1.at(3, 0) == 4);
    CHECK(grid1.at(4, 0) == 5);
    CHECK(grid1.at(0, 1) == 1);
    CHECK(grid1.at(1, 1) == 2);
    CHECK(grid1.at(2, 1) == 3);
    CHECK(grid1.at(3, 1) == 4);
    CHECK(grid1.at(4, 1) == 5);
    CHECK(grid1.at(0, 2) == 1);
    CHECK(grid1.at(1, 2) == 2);
    CHECK(grid1.at(2, 2) == 3);
    CHECK(grid1.at(3, 2) == 4);
    CHECK(grid1.at(4, 2) == 5);
    CHECK(grid1.at(0, 3) == 1);
    CHECK(grid1.at(1, 3) == 2);
    CHECK(grid1.at(2, 3) == 3);
    CHECK(grid1.at(3, 3) == 4);
    CHECK(grid1.at(4, 3) == 5);
  }

  grid1.shift_values(d1);

  SECTION("Check Shift") {
    CHECK(grid1.at(0, 0) == 2);
    CHECK(grid1.at(1, 0) == 3);
    CHECK(grid1.at(2, 0) == 4);
    CHECK(grid1.at(3, 0) == 5);
    CHECK(grid1.at(4, 0) == 5);
    CHECK(grid1.at(0, 1) == 2);
    CHECK(grid1.at(1, 1) == 3);
    CHECK(grid1.at(2, 1) == 4);
    CHECK(grid1.at(3, 1) == 5);
    CHECK(grid1.at(4, 1) == 5);
    CHECK(grid1.at(0, 2) == 2);
    CHECK(grid1.at(1, 2) == 3);
    CHECK(grid1.at(2, 2) == 4);
    CHECK(grid1.at(3, 2) == 5);
    CHECK(grid1.at(4, 2) == 5);
    CHECK(grid1.at(0, 3) == 2);
    CHECK(grid1.at(1, 3) == 3);
    CHECK(grid1.at(2, 3) == 4);
    CHECK(grid1.at(3, 3) == 5);
    CHECK(grid1.at(4, 3) == 5);
  }
}
