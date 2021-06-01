#include "dnn.hpp"
#include "helper.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <exception>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <typeinfo>
#include <vector>
/* #include "data_types.hpp" */

/* DNN::DNN() {} */

void DNN::init() {
  /* read_input_file(); */
  grid = Grid(xsize, ysize, omega);

  set_d_value();

  if (dt < 0) {
    dt = dx * dy / (4 * d);
  }
  grid_Fo = (d * dt) / (dx * dx);

#if GPU

  init_GPU();
  bytes = xsize * ysize * sizeof(double);
  allocate_device_memory(reinterpret_cast<void **>(&d_grid),
                         static_cast<int>(bytes));
  allocate_device_memory(reinterpret_cast<void **>(&d_tgrid),
                         static_cast<int>(bytes));
  allocate_device_memory(reinterpret_cast<void **>(&d_gneedle),
                         xsize * ysize * sizeof(int));
  allocate_device_memory(reinterpret_cast<void **>(&d_needles),
                         needles.size() * sizeof(Needle));

#endif

  std::cout << std::setprecision(10) << "\nomega: " << omega
            << "\ngrid_Fo: " << grid_Fo << "\nd: " << d << "\n\n";

  // Just a work around
  if (nprint) {
    arr_time.reserve(n_iter / nprint + 1);
    arr_vel.reserve(n_iter / nprint + 1);
    arr_rad.reserve(n_iter / nprint + 1);
    arr_len.reserve(n_iter / nprint + 1);
  }
}

void DNN::init_ss() {
  grid = Grid(xsize, ysize, omega);

  set_d_value();

  if (!rfile_name.empty()) {
    read_restart_file(rfile_name);
  }

  needles[0].rad = _rad;
  needles[0].vel = _vel;
}

void DNN::set_d_value() {
  if (non_dim) {
    if (wintegral == JINT) {
      // Sharp needle
      d = M_PI / (2 * omega * omega);
    } else {
      // Thick needle
      double peclet = calculate_peclet_number();
      if (peclet == -1) {
        throw exit_exception("Unable to calcuate the peclet number\n");
      }
      d = 1 / (2 * peclet);
    }
  } else {
    d = dcl;
  }
}

template <typename T>
T DNN::get_json_value(const nlohmann::json jobj, std::string key) {
  // TODO
  // Type issue, implicity converts the float to int
  // Unable to check for it right now
  if (!jobj.contains(key)) {
    if (key.compare("sor") == 0) {
      return SOR;
    } else if (key.compare("transient") == 0) {
      return IS_TRANSIENT;
    } else if (key.compare("epsilon") == 0) {
      return static_cast<T>(EPSILON);
    } else if (key.compare("integral") == 0) {
      return static_cast<T>(JINT);
    } else if (key.compare("dt") == 0) {
      return static_cast<T>(-1.0);
    }

    else if (key.compare("non_dim") == 0) {
      return true;
    } else {
      char buffer[100];
      sprintf(buffer, "Key not found: %s", key.c_str());
      throw std::invalid_argument(buffer);
    }
  }

  T rval = jobj[key].get<T>();

  return rval;
}

std::vector<Needle> DNN::add_needle(const nlohmann::json jobj,
                                    std::string key) {
  if (!jobj.contains(key)) {
    char buffer[100];
    sprintf(buffer, "Key not found: %s", key.c_str());
    throw std::invalid_argument(buffer);
  }

  std::vector<Needle> needles;
  auto jarr = jobj["needle"];

  for (auto &val : jarr) {
    Needle needle;
    needle.x0 = val["x0"].get<std::size_t>();
    needle.y0 = val["y0"].get<std::size_t>();
    needle.xf = val["xf"].get<std::size_t>();
    needle.yf = val["yf"].get<std::size_t>();
    needle.vel = 0;
    needle.rad = 1;
    needle.r = 0;

    if ((needle.x0 != needle.xf) && (needle.y0 != needle.yf)) {
      throw std::invalid_argument("Invalid Needle co-ordinates, needle can "
                                  "only be aligned along the x or y axis.\n");
    }

    needles.push_back(needle);
  }

  return needles;
}

BC DNN::get_bc(const nlohmann::json jobj, std::string key) {
  if (!jobj.contains(key)) {
    char buffer[100];
    sprintf(buffer, "Key not found: %s", key.c_str());
    throw std::invalid_argument(buffer);
  }

  BC rval;

  std::string type = jobj[key].at("type").get<std::string>();
  double val = jobj[key].at("value").get<double>();

  std::transform(type.begin(), type.end(), type.begin(), ::toupper);

  if (type.compare(FLUX) != 0 && type.compare(CONSTANT) != 0) {
    throw std::invalid_argument(
        "Boundary condition type does not match FLUX or CONSTANT: " + key);
  }

  rval.set(type, val);

  return rval;
}

void DNN::read_input_file() {
  using json = nlohmann::json;
  json input;

  std::ifstream ifile(ifile_name);

  if (ifile.is_open()) {
    try {
      ifile >> input;

      // Values required for every case
      xsize = get_json_value<std::size_t>(input, "nx");
      ysize = get_json_value<std::size_t>(input, "ny");
      dx = get_json_value<double>(input, "dx");
      dy = get_json_value<double>(input, "dy");
      n_iter = get_json_value<unsigned long long>(input, "iterations");
      nprint = get_json_value<unsigned int>(
          input, "nprint"); // Value of 0 implies no printing

      // Material Paramters
      dcl = get_json_value<double>(input, "D");
      sigma = get_json_value<double>(input, "sigma");
      c_inf = get_json_value<double>(input, "c_inf");
      c0 = get_json_value<double>(input, "c0");
      k = get_json_value<double>(input, "k");
      m = get_json_value<double>(input, "m");
      gamma = get_json_value<double>(input, "gamma");
      omega = get_json_value<double>(input, "omega");

      // Boundary Conditions
      bc[LW] = get_bc(input, "BC_LEFT_WALL");
      bc[BW] = get_bc(input, "BC_BOTTOM_WALL");
      bc[RW] = get_bc(input, "BC_RIGHT_WALL");
      bc[TW] = get_bc(input, "BC_TOP_WALL");

      // Needles
      needles = add_needle(input, "needle");

      is_transient = get_json_value<bool>(input, "transient");
      non_dim = get_json_value<bool>(input, "non_dim");

      if (is_transient) {
        dt = get_json_value<double>(input, "dt");
      } else {
        sor = get_json_value<double>(input, "sor");
        _epsilon = get_json_value<double>(input, "epsilon");
        _rad = get_json_value<double>(input, "radius");
        _vel = get_json_value<double>(input, "velocity");
      }

      wintegral = get_json_value<int>(input, "integral");
      if (wintegral == JINT) {
        h = get_json_value<int>(input, "h");
      } else {
        A = get_json_value<int>(input, "A");
        B = get_json_value<int>(input, "B");
        C = get_json_value<int>(input, "C");
      }

    } catch (const std::invalid_argument &e) {
      std::string err_msg("Invalid parameter, ");
      err_msg.append(e.what());
      throw exit_exception(err_msg);
    } catch (const std::exception &e) {
      std::string err_msg("Exception while reading the input file: ");
      err_msg.append(ifile_name).append(e.what());
      throw exit_exception(err_msg);
    }
  } else {
    std::string err_msg("Unable to open file: ");
    err_msg.append(ifile_name);
    throw exit_exception(err_msg);
  }

  ifile.close();
}

void DNN::dump_output() {
  std::ofstream fp;
  char buffer[32];

  sprintf(buffer, "dump_%llu.txt", c_iter);
  fp.open(buffer);

  if (fp.is_open()) {
    for (std::size_t i = 0; i < ysize; ++i) {
      for (std::size_t j = 0; j < xsize; ++j) {
        /* printf("%f ", grid.at(j, i)); */
        fp << std::setprecision(10) << grid.at(j, i) << " ";
      }
      fp << "\n";
    }
    fp.close();
  }
}

void DNN::dump_vel_rad_data() {
  std::ofstream fp;
  std::size_t len = arr_time.size();

  fp.open("dump_vel_rad.dat");

  if (fp.is_open()) {
    for (std::size_t i = 0; i < len; ++i) {
      fp << std::setprecision(10) << arr_time[i] << " " << arr_vel[i] << " "
         << arr_rad[i] << " " << arr_len[i] << "\n";
    }
  }

  fp.close();
}

void DNN::run_dnn_model() {
  set_boundary_values();

#if GPU

  copy_to_device_memory(static_cast<void *>(d_grid),
                        static_cast<void *>(grid.data()),
                        static_cast<int>(bytes));
  copy_to_device_memory(static_cast<void *>(d_gneedle),
                        static_cast<void *>(grid.needle_data()),
                        xsize * ysize * sizeof(int));
  copy_to_device_memory(static_cast<void *>(d_needles),
                        static_cast<void *>(needles.data()),
                        needles.size() * sizeof(Needle));

  set_parabolic_needle(d_grid, d_gneedle, d_needles,
                       static_cast<int>(needles.size()));

#else

  double err = 99.99;
  std::vector<double> temp(grid.size());

  set_needle_values();

#endif

  for (c_iter = 0; c_iter < n_iter; ++c_iter) {

    if (nprint != 0 && c_iter % nprint == 0) {

#if GPU

      copy_to_host_memory(static_cast<void *>(grid.data()),
                          static_cast<void *>(d_grid), static_cast<int>(bytes));
      copy_to_host_memory(static_cast<void *>(grid.needle_data()),
                          static_cast<void *>(d_gneedle),
                          grid.size() * sizeof(int));
      copy_to_host_memory(static_cast<void *>(needles.data()),
                          static_cast<void *>(d_needles),
                          needles.size() * sizeof(Needle));

#endif

      dump_output();

      double len = (needles[0].xf - needles[0].x0) * dx + needles[0].r;

      arr_time.push_back(c_iter * dt);
      arr_vel.push_back(needles[0].vel);
      arr_rad.push_back(needles[0].rad);
      arr_len.push_back(len);

      std::cout << "citer: " << c_iter << ", vel: " << needles[0].vel
                << ", rad: " << needles[0].rad << ", needle.len: " << len
                << "\n";
    }

#if GPU
    // TODO:
    run_transient_CUDA(d_grid, d_tgrid, grid.size());
    swap_device_pointers(d_grid, d_tgrid);
    set_parabolic_needle(d_grid, d_gneedle, d_needles,
                         static_cast<int>(needles.size()));
    run_grow_needles(d_grid, d_gneedle, d_needles,
                     static_cast<int>(needles.size()));
    set_parabolic_needle(d_grid, d_gneedle, d_needles,
                         static_cast<int>(needles.size()));
    run_shift_domain(d_grid, d_gneedle, d_needles, static_cast<int>(ysize),
                     xsize / 5);
#else

    err = run_transient(temp);
    grid.swap(temp);
    set_needle_values();
    grow_needles();
    set_needle_values();
    /* fix_concentration(); */
    /* shift_domain(); */

#endif
  }

#if GPU

  copy_to_host_memory(static_cast<void *>(grid.data()),
                      static_cast<void *>(d_grid), static_cast<int>(bytes));
  copy_to_host_memory(static_cast<void *>(grid.needle_data()),
                      static_cast<void *>(d_gneedle),
                      grid.size() * sizeof(int));
  copy_to_host_memory(static_cast<void *>(needles.data()),
                      static_cast<void *>(d_needles),
                      needles.size() * sizeof(Needle));

  cuda_release_memory(static_cast<void *>(d_grid));
  cuda_release_memory(static_cast<void *>(d_tgrid));
  cuda_release_memory(static_cast<void *>(d_gneedle));
  cuda_release_memory(static_cast<void *>(d_needles));

#endif

  dump_output();
  dump_vel_rad_data();
}

void DNN::run_dnn_ss() {
  double err = 99;
  double fif;
  double fif_factor;
  double vel = _vel;
  double rad = _rad;

  set_needle_values();

  for (c_iter = 0; c_iter < n_iter && err > EPSILON; ++c_iter) {

    if (c_iter % nprint == 0) {
      std::cout << "Vel: " << vel << " Rad: " << rad << " err: " << err << "\n";
    }

    err = run_steady_state();

    fif = calculate_flux_intensity_factor(needles[0], A, B, C);

    fif_factor = 2 * d * d * fif * fif;
    vel = pow(fif_factor, 0.6666666666);
    rad = pow(fif_factor, -0.3333333333);

    destroy_needle();

    needles[0].vel = vel;
    needles[0].rad = rad;

    set_needle_values();
  }

  std::cout << "Simulation Done\n";
  std::cout << "Vel: " << vel << " Rad: " << rad << "\n";

  dump_output();
}

void DNN::grow_needles() {
  double fif_sq, vel, rad, fif_factor, fif;

  for (auto &needle : needles) {
    if (wintegral == JINT) {
      fif_sq =
          calculate_flux_intensity_factor_sq_j_integral(needle.xf, needle.yf);
    } else {
      fif = calculate_flux_intensity_factor(needle, A, B, C);
      fif_sq = fif * fif;
      destroy_needle();
    }

#if GROW

    fif_factor = 2 * d * d * fif_sq;
    vel = pow(fif_factor, 0.6666666666);
    rad = pow(fif_factor, -0.3333333333);

    needle.r += (vel * dt) / dx;
    needle.vel = vel;
    needle.rad = rad;

    if (needle.r > 1.0) {
      if (needle.yf == needle.y0) {
        needle.xf = needle.xf + 1;
        needle.r = needle.r - 1;
      } else {
        needle.yf = needle.yf + 1;
        needle.r = needle.r - 1;
      }
    }

#endif
  }
}

void DNN::fix_concentration() {
  std::size_t x, y;
  double u_corr;

  for (auto &needle : needles) {
    x = needle.xf;
    y = needle.yf;

    if (y == needle.y0) {
      double r = needle.r;
      double f1 = 1 - r;
      double f2 = 2 - r;
      double f3 = r * r - 2 * r + 2;
      double f4 = f3 * f1 / 2;
      double f5 = f3 * f1 + f2;
      /* u_corr = f1*grid.at(x+2, y)/(f2*f2) + (f1/(2*f2))*(grid.at(x+1, y+1) +
       * grid.at(x+1, y-1)); */
      u_corr = f1 * grid.at(x + 2, y) +
               f4 * (grid.at(x + 1, y + 1) + grid.at(x + 1, y - 1));
      grid.set_value(x + 1, y, u_corr / f5);
    } else {
      u_corr = ((1 - needle.r) / (4 - 3 * needle.r)) *
               (grid.at(x + 1, y + 1) + grid.at(x - 1, y + 1) +
                grid.at(x, y + 2) + grid.at(x, y));
      grid.set_value(x, y + 1, u_corr);
    }
  }
}

void DNN::fix_concentration_tourret() {
  std::size_t x, y;
  double u_corr;

  for (auto &needle : needles) {
    x = needle.xf;
    y = needle.yf;

    if (y == needle.y0) {
      u_corr = ((1 - needle.r) / (4 - 3 * needle.r)) *
               (grid.at(x + 2, y) + grid.at(x, y) + grid.at(x + 1, y + 1) +
                grid.at(x + 1, y - 1));
      grid.set_value(x + 1, y, u_corr);
    } else {
      u_corr = ((1 - needle.r) / (4 - 3 * needle.r)) *
               (grid.at(x + 1, y + 1) + grid.at(x - 1, y + 1) +
                grid.at(x, y + 2) + grid.at(x, y));
      grid.set_value(x, y + 1, u_corr);
    }
  }
}

double DNN::run_steady_state() {
  double max_error = -1.0f, error, nval;
  for (std::size_t i = 0; i < ysize; ++i) {   // y-coordinate
    for (std::size_t j = 0; j < xsize; ++j) { // x-coordinate

      if (grid.is_needle(j, i)) {
        continue;
      }

      if (j == 0 && i == 0) {
        // left bottom corner
        if (bc[LW].type == FLUX && bc[BW].type == FLUX) {
          /*
           * T(0,0) = 2*T(1, 0) - 2*dx*flux_x + 2*T(0, 1) - 2*dy*flux_y/4
           */
          nval = (2 * grid.at(j + 1, i) + 2 * dx * bc[LW].value / d +
                  2 * grid.at(j, i + 1) + 2 * dy * bc[BW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else if (bc[LW].type == FLUX && bc[BW].type == CONSTANT) {
          nval = bc[BW].value;
        } else if (bc[LW].type == CONSTANT && bc[BW].type == FLUX) {
          nval = bc[LW].value;
        } else {
          /*
           * T(0,0) = 2*const_val_x + 2*const_val_y/4
           */
          nval = (2 * bc[LW].value + 2 * bc[BW].value) / 4;
        }
      } else if (j == xsize - 1 && i == 0) {
        // right bottom corner
        if (bc[RW].type == FLUX && bc[BW].type == FLUX) {
          nval = (2 * grid.at(j - 1, i) - 2 * dx * bc[RW].value / d +
                  2 * grid.at(j, i + 1) + 2 * dy * bc[BW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else if (bc[RW].type == FLUX && bc[BW].type == CONSTANT) {
          nval = bc[BW].value;
        } else if (bc[RW].type == CONSTANT && bc[BW].type == FLUX) {
          nval = bc[RW].value;
        } else {
          nval = (2 * bc[RW].value + 2 * bc[BW].value) / 4;
        }
      } else if (j == 0 && i == ysize - 1) {
        // left top corner
        if (bc[LW].type == FLUX && bc[TW].type == FLUX) {
          nval = (2 * grid.at(j + 1, i) + 2 * dx * bc[LW].value / d +
                  2 * grid.at(j, i - 1) - 2 * dy * bc[TW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else if (bc[LW].type == FLUX && bc[TW].type == CONSTANT) {
          nval = bc[TW].value;
        } else if (bc[LW].type == CONSTANT && bc[TW].type == FLUX) {
          nval = bc[LW].value;
        } else {
          nval = (2 * bc[LW].value + 2 * bc[TW].value) / 4;
        }
      } else if (j == xsize - 1 && i == ysize - 1) {
        // right top corner
        if (bc[RW].type == FLUX && bc[TW].type == FLUX) {
          nval = (2 * grid.at(j - 1, i) - 2 * dx * bc[RW].value / d +
                  2 * grid.at(j, i - 1) - 2 * dy * bc[TW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else if (bc[RW].type == FLUX && bc[TW].type == CONSTANT) {
          nval = bc[TW].value;
        } else if (bc[RW].type == CONSTANT && bc[TW].type == FLUX) {
          nval = bc[RW].value;
        } else {
          nval = (2 * bc[RW].value + 2 * bc[TW].value) / 4;
        }
      } else if (j > 0 && j < xsize - 1 && i == 0) {
        // bottom wall
        if (bc[BW].type == FLUX) {
          nval = (grid.at(j + 1, i) + grid.at(j - 1, i) +
                  2 * grid.at(j, i + 1) + 2 * dy * bc[BW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else {
          nval = bc[BW].value;
        }
      } else if (j == 0 && i > 0 && i < ysize - 1) {
        // left wall
        if (bc[LW].type == FLUX) {
          nval = (2 * grid.at(j + 1, i) + 2 * dx * bc[LW].value / d +
                  grid.at(j, i + 1) + grid.at(j, i - 1)) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else {
          nval = bc[LW].value;
        }
      } else if (j > 0 && j < xsize - 1 && i == ysize - 1) {
        // top wall
        if (bc[TW].type == FLUX) {
          nval = (grid.at(j + 1, i) + grid.at(j - 1, i) +
                  2 * grid.at(j, i - 1) - 2 * dy * bc[TW].value / d) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else {
          nval = bc[TW].value;
        }
      } else if (j == xsize - 1 && i > 0 && i < ysize - 1) {
        // right wall
        if (bc[RW].type == FLUX) {
          nval = (2 * grid.at(j - 1, i) - 2 * dx * bc[RW].value / d +
                  grid.at(j, i + 1) + grid.at(j, i - 1)) /
                 4;
          nval = sor * nval + (1 - sor) * grid.at(j, i);
        } else {
          nval = bc[RW].value;
        }
      } else {
        // inner domain
        nval = (grid.at(j + 1, i) + grid.at(j - 1, i) + grid.at(j, i + 1) +
                grid.at(j, i - 1)) /
               4;
        nval = sor * nval + (1 - sor) * grid.at(j, i);
      }
      error = std::abs((nval - grid.at(j, i)) / nval);
      grid.set_value(j, i, nval);
      max_error = error > max_error ? error : max_error;
    }
  }
  return max_error;
}

double DNN::run_transient(std::vector<double> &out) {
  double max_error = -1.0f, error, nval;
  std::size_t xp, xn, yp, yn;

  for (std::size_t i = 1; i < ysize - 1; ++i) {   // y-coordinate
    for (std::size_t j = 1; j < xsize - 1; ++j) { // x-coordinate
      // inner domain
      xp = j - 1;
      xn = j + 1;
      yp = i - 1;
      yn = i + 1;
      nval = grid.at(j, i) +
             grid_Fo * (grid.at(xn, i) + grid.at(xp, i) + grid.at(j, yn) +
                        grid.at(j, yp) - 4 * grid.at(j, i));
      error = (std::abs(nval - grid.at(j, i)) / nval);
      out[grid.index(j, i)] = nval;
      max_error = error > max_error ? error : max_error;
    }
  }

  // left bottom corner
  xn = 1;
  yn = 1;
  if (bc[LW].type == FLUX && bc[BW].type == FLUX) {
    /*
     * T(0,0) = 2*T(1, 0) - 2*dx*flux_x + 2*T(0, 1) - 2*dy*flux_y/4
     */
    nval = grid.at(0, 0) +
           grid_Fo * (2 * grid.at(xn, 0) + (2 * dx * bc[LW].value) / d +
                      2 * grid.at(0, yn) + (2 * dy * bc[BW].value) / d -
                      4 * grid.at(0, 0));
  } else if (bc[LW].type == FLUX && bc[BW].type == CONSTANT) {
    nval = bc[BW].value;
  } else if (bc[LW].type == CONSTANT && bc[BW].type == FLUX) {
    nval = bc[LW].value;
  } else {
    nval = (bc[LW].value + bc[BW].value) / 2;
  }
  out[grid.index(0, 0)] = nval;

  // right bottom corner
  xp = xsize - 2;
  yn = 1;
  if (bc[RW].type == FLUX && bc[BW].type == FLUX) {
    nval = grid.at(xsize - 1, 0) +
           grid_Fo * (2 * grid.at(xp, 0) - (2 * dx * bc[RW].value) / d +
                      2 * grid.at(xsize - 1, yn) + (2 * dy * bc[BW].value) / d -
                      4 * grid.at(xsize - 1, 0));
  } else if (bc[RW].type == FLUX && bc[BW].type == CONSTANT) {
    nval = bc[BW].value;
  } else if (bc[RW].type == CONSTANT && bc[BW].type == FLUX) {
    nval = bc[RW].value;
  } else {
    nval = (bc[RW].value + bc[BW].value) / 2;
  }
  out[grid.index(xsize - 1, 0)] = nval;

  // left top corner
  xn = 1;
  yp = ysize - 2;
  if (bc[LW].type == FLUX && bc[TW].type == FLUX) {
    nval = grid.at(0, ysize - 1) +
           grid_Fo * (2 * grid.at(xn, ysize - 1) + (2 * dx * bc[LW].value) / d +
                      2 * grid.at(0, yp) - (2 * dy * bc[TW].value) / d -
                      4 * grid.at(0, ysize - 1));
  } else if (bc[LW].type == FLUX && bc[TW].type == CONSTANT) {
    nval = bc[TW].value;
  } else if (bc[LW].type == CONSTANT && bc[TW].type == FLUX) {
    nval = bc[LW].value;
  } else {
    nval = (bc[LW].value + bc[TW].value) / 2;
  }
  out[grid.index(0, ysize - 1)] = nval;

  // right top corner
  xp = xsize - 2;
  yp = ysize - 2;
  if (bc[RW].type == FLUX && bc[TW].type == FLUX) {
    nval = grid.at(xsize - 1, ysize - 1) +
           grid_Fo * (2 * grid.at(xp, ysize - 1) - (2 * dx * bc[RW].value) / d +
                      2 * grid.at(xsize - 1, yp) - (2 * dy * bc[TW].value) / d -
                      4 * grid.at(xsize - 1, ysize - 1));
  } else if (bc[RW].type == FLUX && bc[TW].type == CONSTANT) {
    nval = bc[TW].value;
  } else if (bc[RW].type == CONSTANT && bc[TW].type == FLUX) {
    nval = bc[RW].value;
  } else {
    nval = (bc[RW].value + bc[TW].value) / 2;
  }
  out[grid.index(xsize - 1, ysize - 1)] = nval;

  // bottom wall
  yn = 1;
  for (std::size_t j = 1; j < xsize - 1; ++j) {
    xp = j - 1;
    xn = j + 1;
    if (bc[BW].type == FLUX) {
      nval = grid.at(j, 0) +
             grid_Fo * (grid.at(xn, 0) + grid.at(xp, 0) + 2 * grid.at(j, yn) +
                        (2 * dy * bc[BW].value) / d - 4 * grid.at(j, 0));
    } else {
      nval = bc[BW].value;
    }
    out[grid.index(j, 0)] = nval;
  }

  // left wall
  xn = 1;
  for (std::size_t i = 1; i < ysize - 1; ++i) {
    yp = i - 1;
    yn = i + 1;
    if (bc[LW].type == FLUX) {
      nval = grid.at(0, i) +
             grid_Fo * (2 * grid.at(xn, i) + (2 * dx * bc[LW].value) / d +
                        grid.at(0, yn) + grid.at(0, yp) - 4 * grid.at(0, i));
    } else {
      nval = bc[LW].value;
    }
    out[grid.index(0, i)] = nval;
  }

  // top wall
  yp = ysize - 2;
  for (std::size_t j = 1; j < xsize - 1; ++j) {
    xp = j - 1;
    xn = j + 1;
    if (bc[TW].type == FLUX) {
      nval = grid.at(j, ysize - 1) +
             grid_Fo * (grid.at(xn, ysize - 1) + grid.at(xp, ysize - 1) +
                        2 * grid.at(j, yp) - (2 * dy * bc[TW].value) / d -
                        4 * grid.at(j, ysize - 1));
    } else {
      nval = bc[TW].value;
    }
    out[grid.index(j, ysize - 1)] = nval;
  }

  // right wall
  xp = xsize - 2;
  for (std::size_t i = 1; i < ysize - 1; ++i) {
    yp = i - 1;
    yn = i + 1;
    if (bc[RW].type == FLUX) {
      nval = grid.at(xsize - 1, i) +
             grid_Fo * (2 * grid.at(xp, i) - (2 * dx * bc[RW].value) / d +
                        grid.at(xsize - 1, yn) + grid.at(xsize - 1, yp) -
                        4 * grid.at(xsize - 1, i));
    } else {
      nval = bc[RW].value;
    }
    out[grid.index(xsize - 1, i)] = nval;
  }

  return max_error;
}

void DNN::set_boundary_values() {
  /*
   * The purpose of this method to change the initial values at the boundary
   */
  // left bottom corner
  if (bc[LW].type == CONSTANT && bc[BW].type == CONSTANT) {
    grid.set_value(0, 0, (bc[LW].value + bc[BW].value) / 2);
  } else if (bc[LW].type == CONSTANT || bc[BW].type == CONSTANT) {
    grid.set_value(0, 0, bc[LW].type == CONSTANT ? bc[LW].value : bc[BW].value);
  }

  // left top corner
  if (bc[LW].type == CONSTANT && bc[TW].type == CONSTANT) {
    grid.set_value(0, ysize - 1, (bc[LW].value + bc[TW].value) / 2);
  } else if (bc[LW].type == CONSTANT || bc[TW].type == CONSTANT) {
    grid.set_value(0, ysize - 1,
                   bc[LW].type == CONSTANT ? bc[LW].value : bc[TW].value);
  }

  // right top corner
  if (bc[TW].type == CONSTANT && bc[RW].type == CONSTANT) {
    grid.set_value(xsize - 1, ysize - 1, (bc[TW].value + bc[RW].value) / 2);
  } else if (bc[TW].type == CONSTANT || bc[RW].type == CONSTANT) {
    grid.set_value(xsize - 1, ysize - 1,
                   bc[TW].type == CONSTANT ? bc[TW].value : bc[RW].value);
  }

  // right bottom corner
  if (bc[RW].type == CONSTANT && bc[BW].type == CONSTANT) {
    grid.set_value(xsize - 1, 0, (bc[RW].value + bc[BW].value) / 2);
  } else if (bc[RW].type == CONSTANT || bc[BW].type == CONSTANT) {
    grid.set_value(xsize - 1, 0,
                   bc[RW].type == CONSTANT ? bc[RW].value : bc[BW].value);
  }

  if (bc[LW].type == CONSTANT) {
    for (std::size_t i = 1; i < ysize - 1; ++i) {
      grid.set_value(0, i, bc[LW].value);
    }
  }

  if (bc[TW].type == CONSTANT) {
    for (std::size_t i = 1; i < xsize - 1; ++i) {
      grid.set_value(i, ysize - 1, bc[TW].value);
    }
  }

  if (bc[RW].type == CONSTANT) {
    for (std::size_t i = 1; i < ysize - 1; ++i) {
      grid.set_value(xsize - 1, i, bc[RW].value);
    }
  }

  if (bc[BW].type == CONSTANT) {
    for (std::size_t i = 1; i < xsize - 1; ++i) {
      grid.set_value(i, 0, bc[BW].value);
    }
  }
}

void DNN::set_needle_values() {
  for (auto &needle : needles) {
    if (wintegral == JINT) {
      if (needle.y0 == needle.yf) {
        for (std::size_t i = needle.x0; i <= needle.xf; ++i) {
          grid.set_value(i, needle.y0, 0);
          grid.set_needle(i, needle.y0, true);
        }
      } else {
        for (std::size_t i = needle.y0; i <= needle.yf; ++i) {
          grid.set_value(needle.x0, i, 0);
          grid.set_needle(needle.x0, i, true);
        }
      }
    } else {
      // Assuming needle only aligned along the x-axis
      for (std::size_t i = needle.x0; i <= needle.xf; ++i) {
        std::size_t wd = static_cast<std::size_t>(
            sqrt(2 * needle.rad * (needle.xf - i) * dx) / dy);

        grid.set_value(i, needle.y0, 0);
        grid.set_needle(i, needle.y0, true);

        /* wd = wd < 3 ? wd : 3; */
        for (std::size_t j = 1; j <= wd; ++j) {
          grid.set_value(i, needle.y0 + j, 0);
          grid.set_value(i, needle.y0 - j, 0);
          grid.set_needle(i, needle.y0 + j, true);
          grid.set_needle(i, needle.y0 - j, true);
        }
      }
    }
  }
}

void DNN::destroy_needle() {
  for (auto &needle : needles) {
    for (std::size_t i = needle.x0; i <= needle.xf; ++i) {
      std::size_t wd = static_cast<std::size_t>(
          sqrt(2 * needle.rad * (needle.xf - i) * dx) / dy);

      grid.set_needle(i, needle.y0, false);

      /* wd = wd < 3 ? wd : 3; */
      for (std::size_t j = 1; j <= wd; ++j) {
        grid.set_needle(i, needle.y0 + j, false);
        grid.set_needle(i, needle.y0 - j, false);
      }
    }
  }
}

double DNN::calculate_flux_intensity_factor(const Needle &needle,
                                            const std::size_t A,
                                            const std::size_t B,
                                            const std::size_t C) {
  /*
   * A, B, C -> Tourret 2016
   * xt, yt -> tip coordinates
   * VALID ONLY FOR HORIZONTAL TIP
   */
  std::size_t yf, yb;
  const std::size_t xb = needle.xf - A;
  const std::size_t xf = needle.xf + B;
  static double factor = 1 / (4 * sqrt((A + needle.r + 1 / 2) * dx));
  double fval = 0.0;
  std::size_t y_nt = needle.y0;
  std::size_t y_nb = needle.y0;

  for (std::size_t i = 0; i < ysize; ++i) {
    if (!grid.is_needle(xb, needle.y0 + i)) {
      y_nt = needle.y0 + i - 1;
      y_nb = needle.y0 - i + 1;
      break;
    }
  }

  yf = y_nt + C;
  yb = y_nb - C;

  double line = calculate_line_integral(xb, xf, yb, yf, y_nb, y_nt);
  double surface = calculate_surface_integral(xb, xf, yb, yf);
  surface = (surface * needle.vel) / d;
  fval = factor * (line + surface);

  return fval;
}

double DNN::calculate_line_integral(const std::size_t xb, const std::size_t xf,
                                    const std::size_t yb, const std::size_t yf,
                                    const std::size_t y_nb,
                                    const std::size_t y_nt) {
  double flux_acc = 0.0;
  double bflux = 0.0;
  double rflux = 0.0;
  double tflux = 0.0;
  double lflux = 0.0;
  double dyu, dxu;

  // Flux intensity in bottom contour line
  for (std::size_t i = xb - 1; i <= xf; ++i) {
    std::size_t xp{i}, xn{i + 1}, yp{yb - 1}, yn{yb};
    dyu = (grid.at(xp, yn) + grid.at(xn, yn) - grid.at(xp, yp) -
           grid.at(xn, yp)) /
          (2 * dy);
    if (i == xb - 1 || i == xf) { // half contour
      bflux += -dyu * dx / 2;
    } else {
      bflux += -dyu * dx;
    }
  }

  // flux intensity in front contour line
  for (std::size_t i = yb - 1; i <= yf; ++i) {
    std::size_t xp{xf}, xn{xf + 1}, yp{i}, yn{i + 1};
    dxu = (grid.at(xn, yp) + grid.at(xn, yn) - grid.at(xp, yp) -
           grid.at(xp, yn)) /
          (2 * dx);
    if (i == yb - 1 || i == yf) { // half contour
      rflux += dxu * dy / 2;
    } else {
      rflux += dxu * dy;
    }
  }

  // Flux intensity in top contour line
  for (std::size_t i = xb - 1; i <= xf; ++i) {
    std::size_t xp{i}, xn{i + 1}, yp{yf}, yn{yf + 1};
    dyu = ((grid.at(xp, yn) + grid.at(xn, yn)) -
           (grid.at(xp, yp) + grid.at(xn, yp))) /
          (2 * dy);
    if (i == xb - 1 || i == xf) { // half contour
      tflux += dyu * dx / 2;
    } else {
      tflux += dyu * dx;
    }
  }

  // Flux intensity in back contour line top half
  for (std::size_t i = y_nt; i <= yf; ++i) {
    std::size_t xp{xb - 1}, xn{xb}, yp{i}, yn{i + 1};
    dxu = (grid.at(xn, yp) + grid.at(xn, yn) - grid.at(xp, yp) -
           grid.at(xp, yn)) /
          (2 * dx);
    if (i == yf) { // half contour
      lflux += -dxu * dy / 2;
    } else {
      lflux += -dxu * dy;
    }
  }

  // Flux intensity in back contour line bottom half
  for (std::size_t i = yb - 1; i < y_nb; ++i) {
    std::size_t xp{xb - 1}, xn{xb}, yp{i}, yn{i + 1};
    dxu = (grid.at(xn, yp) + grid.at(xn, yn) - grid.at(xp, yp) -
           grid.at(xp, yn)) /
          (2 * dx);
    if (i == yb - 1) { // half contour
      lflux += -dxu * dy / 2;
    } else {
      lflux += -dxu * dy;
    }
  }

  flux_acc = bflux + rflux + tflux + lflux;

  return flux_acc;
}

double DNN::calculate_surface_integral(const std::size_t hbc,
                                       const std::size_t hfc,
                                       const std::size_t vbc,
                                       const std::size_t vfc) {
  double integral = 0.0;
  double dxu;

  // internal
  for (std::size_t y = vbc; y <= vfc; ++y) {
    for (std::size_t x = hbc; x <= hfc; ++x) {
      // The flux inside the needle will be zero
      if (!grid.is_needle(x, y)) {
        dxu = (grid.at(x + 1, y) - grid.at(x - 1, y)) / (2 * dx);
        integral += dxu * dx * dy;
      }
    }
  }

  return integral;
}

double DNN::calculate_flux_intensity_factor_sq_j_integral(std::size_t xpos,
                                                          std::size_t ypos) {
  /*
   * Valid for sharp needle
   * (xpos, ypos) -> needle tip coordiantes
   * dx = dy = dx (modify it to reflect that)
   */
  std::size_t hfc = xpos + h; // horizontal front contour pos
  std::size_t hbc = xpos - h; // horizontal back contour pos
  std::size_t vfc = ypos + h; // vertical front pos
  std::size_t vbc = ypos - h; // vertical bac pos
  double flux_sq_acc = 0.0;   // c++11 uniform initilization
  double dxU = 0.0;
  double dyU = 0.0;
  double nx, ny;

  // flux intensity in front contour line
  for (std::size_t i = vbc; i < vfc; ++i) {
    std::size_t xp{hfc - 1}, xn{hfc}, yp{i}, yn{i + 1};
    nx = 1.0;
    ny = 0.0;
    dxU = (grid.at(xn, yp) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xp, yn));
    dyU = (grid.at(xp, yn) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xn, yp));
    if (i == vbc || i == vfc - 1) { // half contour
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (8 * dx);
    } else {
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (4 * dx);
    }
  }

  // Flux intensity in back contour line
  for (std::size_t i = vbc; i < vfc; ++i) {
    std::size_t xp{hbc}, xn{hbc + 1}, yp{i}, yn{i + 1};
    nx = -1.0;
    ny = 0.0;
    dxU = (grid.at(xn, yp) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xp, yn));
    dyU = (grid.at(xp, yn) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xn, yp));
    if (i == vbc || i == vfc - 1) { // half contour
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (8 * dx);
    } else {
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (4 * dx);
    }
  }

  // Flux intensity in top contour line
  for (std::size_t i = hbc; i < hfc; ++i) {
    std::size_t xp{i}, xn{i + 1}, yp{vfc - 1}, yn{vfc};
    nx = 0.0;
    ny = 1.0;
    dxU = (grid.at(xn, yp) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xp, yn));
    dyU = (grid.at(xp, yn) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xn, yp));
    if (i == hbc || i == hfc - 1) { // half contour
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (8 * dx);
    } else {
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (4 * dx);
    }
  }

  // Flux intensity in bottom contour line
  for (std::size_t i = hbc; i < hfc; ++i) {
    std::size_t xp{i}, xn{i + 1}, yp{vbc}, yn{vbc + 1};
    nx = 0.0;
    ny = -1.0;
    dxU = (grid.at(xn, yp) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xp, yn));
    dyU = (grid.at(xp, yn) + grid.at(xn, yn)) -
          (grid.at(xp, yp) + grid.at(xn, yp));
    if (i == hbc || i == hfc - 1) { // half contour
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (8 * dx);
    } else {
      flux_sq_acc +=
          ((dxU * dxU - dyU * dyU) * nx + 2 * dxU * dyU * ny) / (4 * dx);
    }
  }

  flux_sq_acc = flux_sq_acc / (2 * M_PI);

  return flux_sq_acc;
}

void DNN::read_restart_file(std::string restart_file) {
  std::ifstream ifile(restart_file);

  if (ifile.is_open()) {
    double val;

    for (std::size_t i = 0; i < ysize; ++i) {
      for (std::size_t j = 0; j < xsize; ++j) {
        ifile >> val;
        grid.set_value(j, i, val);
      }
    }

  } else {
    throw std::runtime_error(
        "Unable to open the restart file: " + restart_file + "\n");
  }

  ifile.close();
}

double DNN::calculate_peclet_number() {
  int status;
  int iter = 0, max_iter = 100;
  const gsl_root_fsolver_type *T;
  gsl_root_fsolver *s;
  double root = 0;
  double x_lo = 0.0, x_hi = 1.0;
  gsl_function F;

  peclet_params params;
  params.omega = omega;

  F.function = &peclet_number;
  F.params = &params;

  T = gsl_root_fsolver_brent;
  s = gsl_root_fsolver_alloc(T);
  gsl_root_fsolver_set(s, &F, x_lo, x_hi);

  do {
    iter++;
    status = gsl_root_fsolver_iterate(s);
    root = gsl_root_fsolver_root(s);
    x_lo = gsl_root_fsolver_x_lower(s);
    x_hi = gsl_root_fsolver_x_upper(s);
    status = gsl_root_test_interval(x_lo, x_hi, 0, 0.001);

    /* if (status == GSL_SUCCESS) */
    /*     printf("Converged:\n"); */

    /* printf("%5d [%.7f, %.7f] %.7f %.7f\n", iter, x_lo, x_hi, root, x_hi -
     * x_lo); */
  } while (status == GSL_CONTINUE && iter < max_iter);

  gsl_root_fsolver_free(s);

  return status == GSL_SUCCESS ? root : -1.0;
}

inline void DNN::shift_domain() {
  // This is just a temporary fix to run the simulations quickly
  if (needles[0].xf > xsize / 5) {
    grid.shift_values();
    needles[0].xf = needles[0].xf - 1;
  }
}

void DNN::init_GPU() const {
  int bc_type[4];
  double bc_values[4];

  for (std::size_t i = 0; i < 4; ++i) {
    bc_type[i] =
        bc[i].type.compare(static_cast<std::string>(CONSTANT)) == 0 ? 0 : 1;
    bc_values[i] = bc[i].value;
  }

#if GPU
  init_CUDA(grid_Fo, d, dx, dy, dt, xsize, ysize, A, B, C, bc_type, bc_values);
#endif
}

/* void DNN::calcualte_integral_thin_needle() { */
/*   double flux_acc{0.0}; */
/*   // For thin needle */
/*   for(std::size_t i = hbc; i < xt; ++i) { */
/*     std::size_t xp{i}, xn{i+1}, y{yt}, yp{yt-1}, yn{yt+1}; */
/*     double dyUn, dyUp; */
/*     dyUn = (grid.at(xp, yn) + grid.at(xn, yn) - grid.at(xp, y) - grid.at(xn,
 * y))/(2*dy); */
/*     /1* if (i == xt-1) { *1/ */
/*     dyUp = -1*(grid.at(xp, y) + grid.at(xn, y) - grid.at(xp, yp) -
 * grid.at(xn, yp))/(2*dy); */
/*     /1*   flux_acc += dyUn*(dx+needle.r) + dyUp*(dx+needle.r); *1/ */
/*     /1* } else { *1/ */
/*       flux_acc += (dyUn + dyUp)*dx; */
/*     /1* } *1/ */
/*   } */

/*   /1* flux_acc += (grid.at(xt+1, yt) - grid.at(xt, yt)/((1-needle.r)*dx))*dy;
 * *1/ */
/*   /1* flux_acc += ((grid.at(xt+1, yt) + grid.at(xt+1, yt-1) - grid.at(xt, yt)
 * - grid.at(xt, yt-1))/dx)*dy/2; *1/ */
/*   /1* flux_acc += ((grid.at(xt+1, yt) + grid.at(xt+1, yt+1) - grid.at(xt, yt)
 * - grid.at(xt, yt+1))/dx)*dy/2; *1/ */
/*   /1* flux_acc += (grid.at(xt+1, yt) - grid.at(xt, yt))/((1-needle.r)*dx);
 * *1/ */
/*   flux_acc += ((1-needle.r)*(1-needle.r)*grid.at(xt+2, yt) - grid.at(xt, yt)
 * + needle.r*grid.at(xt+1, yt))/((1-needle.r)*(2-needle.r)); */
/*   fval = factor*flux_acc; */
/* } */
