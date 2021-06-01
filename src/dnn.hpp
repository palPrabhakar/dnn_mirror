#ifndef DNN_H
#define DNN_H

#include "data_types.hpp"
#include "grid.hpp"
#include "json.hpp"
#include <cstddef>
#include <cstdlib>
#include <exception>
#include <fstream>
#include <string>
#include <vector>

/* #include <cuda_runtime.h> */

/*
 * Both constructors are almost similar
 * Is there a need for no value constructor??
 * Since the vector will be default initialized on creation
 * Something to think about later!!
 */

#ifndef GPU
#define GPU 1
#endif

#ifndef GROW
#define GROW 1
#endif

#ifndef NON_DIM
#define NON_DIM 1
#endif

#define FLUX "FLUX"
#define CONSTANT "CONST"
#define LW 0
#define BW 1
#define RW 2
#define TW 3
#define JINT 1
#define NORMAL 2 // Thick needle contour integral

// Default simulation values
const bool IS_TRANSIENT = true;
const double SOR = 1.0;
const double EPSILON = 1e-5;
const double RADIUS = 1.0;

class DNN {
public:
  DNN(std::string file_name) : ifile_name(file_name), bc(std::vector<BC>(4)) {}
  DNN(std::string file_name, std::string _rfile_name)
      : ifile_name(file_name), rfile_name(_rfile_name), bc(std::vector<BC>(4)) {
  }

  void init();

  void init_ss();

  void init_GPU() const;

  // Read input file and set simulation parameters
  void read_input_file();

  const std::string get_input_file_name() const { return ifile_name; }

  const std::size_t get_x_size() const { return xsize; }

  const std::size_t get_y_size() const { return ysize; }

  const double get_dx() const { return dx; }

  const double get_dy() const { return dy; }

  const double get_dt() const { return dt; }

  const unsigned long long get_niter() const { return n_iter; }

  const unsigned long long get_citer() const { return c_iter; }

  const unsigned int get_nprint() const { return nprint; }

  const double get_epsilon() const { return _epsilon; }

  const double get_sor() const { return sor; }

  const double get_grid_Fo() const { return grid_Fo; }

  const double get_dcl() const { return dcl; }

  const double get_sigma() const { return sigma; }

  const double get_c_inf() const { return c_inf; }

  const double get_c0() const { return c0; }

  const double get_k() const { return k; }

  const double get_m() const { return m; }

  const double get_gamma() const { return gamma; }

  const double get_omega() const { return omega; }

  const double get_d() const { return d; }

  const int get_wintegral() const { return wintegral; }

  const unsigned int get_h() const { return h; }

  const unsigned int get_A() const { return A; }

  const unsigned int get_B() const { return B; }

  const unsigned int get_C() const { return C; }

  const bool get_is_transient() const { return is_transient; }

  std::vector<Needle> &get_needles() { return needles; }

  const BC get_BC(int bwall) const { return bc[bwall]; }

  Grid &get_grid() { return grid; }

  inline void shift_domain();

  void set_d_value();

  double calculate_peclet_number();

  // Dump output to a file
  void dump_output();

  void dump_vel_rad_data();

  // Read restart file
  void read_restart_file(std::string);

  // Run Simulations
  void run_dnn_model();

  void run_dnn_ss();

  double run_steady_state();

  double run_transient(std::vector<double> &);

  void set_boundary_values();

  void set_needle_values();

  void destroy_needle();

  void grow_needles();

  void fix_concentration();

  void fix_concentration_tourret();

  // FIF calculation for sharp needle
  double calculate_flux_intensity_factor_sq_j_integral(std::size_t,
                                                       std::size_t);

  double calculate_flux_intensity_factor(const Needle &, const std::size_t,
                                         const std::size_t, const std::size_t);

  double calculate_line_integral(const std::size_t, const std::size_t,
                                 const std::size_t, const std::size_t,
                                 const std::size_t, const std::size_t);

  double calculate_surface_integral(const std::size_t, const std::size_t,
                                    const std::size_t, const std::size_t);

private:
  std::string ifile_name;   // input file name
  std::string rfile_name;   // restart file
  std::size_t xsize, ysize; // no. of grid points in x and y dir
  unsigned int nprint;      // print_every n iterations
  double dx, dy;            // grid spacing

  unsigned long long n_iter; // total no. of max iterations
  unsigned long long c_iter; // current iteration

  bool is_transient = false; // transient problem
  bool non_dim = true;
  double dt;       // time step
  double sor;      // successive over relaxation factor
  double _epsilon; // residual threshold

  int wintegral;  // FIF integration type 1 --> J Integral, 2 --> Thick needle
                  // integral
  unsigned int h; // half contour length thin needle
  unsigned int A, B, C; // thick needle contour parameters

  double dcl;   // dcl - diffusion coeff liq
  double sigma; // tip selection parameter
  double c_inf; // far field concentration
  double c0;    // eqmb (planar interface) liquidus concentration at referecnce
                // temprature T0
  double k;     // partition coeff.
  double m;     // liquidus slope
  double gamma; // gibbs-thomson coefficient
  double omega; // supersaturation
  double grid_Fo; // grid fourier number
  double d;       // reduced dcl
  double d_o;     // capillary length
  double rho_s;   // steady state radius

  // for ss run
  double _rad;
  double _vel;

  // Grid
  Grid grid;

  std::vector<BC> bc;
  std::vector<Needle> needles;

  // Device Resources
  size_t bytes;
  double *d_grid;
  double *d_tgrid;
  int *d_gneedle;
  Needle *d_needles;

  // Dump store
  std::vector<double> arr_time; // time array
  std::vector<double> arr_vel;  // vel array
  std::vector<double> arr_rad;  // rad array
  std::vector<double> arr_len;  // len array

  template <typename T> T get_json_value(const nlohmann::json, std::string);

  BC get_bc(const nlohmann::json, std::string);

  std::vector<Needle> add_needle(const nlohmann::json, std::string);
};

#endif
