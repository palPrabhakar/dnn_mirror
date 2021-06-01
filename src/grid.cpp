#include "grid.hpp"
#include <algorithm>
#include <vector>

Grid::Grid() {}

Grid::Grid(std::size_t x, std::size_t y, double c0) : h_sz(x), v_sz(y) {
  ufield = std::vector<double>(x * y, c0);
  needles = std::vector<int>(x * y, 0);
}

void Grid::set_uvalues(std::vector<double> &vec) { ufield = vec; }

/* void Grid::set_needle_values(std::vector<bool>& vec) { */
/*   needles = vec; */
/* } */

void Grid::shift_values() {
  for (std::size_t i = 0; i < h_sz * v_sz; i += h_sz) {
    std::copy(ufield.begin() + i + 1, ufield.begin() + i + h_sz,
              ufield.begin() + i);
    std::copy(needles.begin() + i + 1, needles.begin() + i + h_sz,
              needles.begin() + i);
  }
}
