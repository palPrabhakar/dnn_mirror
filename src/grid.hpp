#ifndef GRID_H
#define GRID_H

#include <cstdlib>
#include <vector>

class Grid {
public:
  Grid();

  Grid(std::size_t, std::size_t, double c0);

  // Get
  inline const double at(std::size_t x, std::size_t y) const {
    return ufield[index(x, y)];
  }

  inline const double at(std::size_t idx) const { return ufield[idx]; }

  inline const bool is_needle(std::size_t x, std::size_t y) const {
    return needles[index(x, y)] == 1 ? true : false;
  }

  // Set
  inline void set_value(std::size_t x, std::size_t y, double val) {
    ufield[index(x, y)] = val;
  }

  inline void set_needle(std::size_t x, std::size_t y, bool val) {
    needles[index(x, y)] = val == true ? 1 : 0;
  }

  const std::size_t size() const { return h_sz * v_sz; }

  inline const std::size_t index(std::size_t x, std::size_t y) const {
    return x + h_sz * y;
  }

  void swap(std::vector<double> &vec) { ufield.swap(vec); }

  double *data() { return ufield.data(); }

  int *needle_data() { return needles.data(); }

  void set_uvalues(std::vector<double> &vec);

  void shift_values();

private:
  std::size_t h_sz, v_sz; // h_sz == xsize
  std::vector<double> ufield;
  std::vector<int> needles;
};

#endif /* ifndef GRID_H */
