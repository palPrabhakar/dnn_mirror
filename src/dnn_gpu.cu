/*
  This file contains the gpu version of the DNN functions
 */
#include "data_types.hpp"
#include <assert.h>
#include <helper_cuda.h>
#include <math.h>
#include <stdio.h>
#include <thrust/copy.h>

#define LW 0 // left wall
#define BW 1 // bottom wall
#define RW 2 // right wall
#define TW 3 // top wall

#define CONST 0
#define FLUX 1

#define BLOCK_SIZE 128

__constant__ double gridFo;
__constant__ double d;
__constant__ double dx;
__constant__ double dy;
__constant__ double dt;
__constant__ unsigned int xsize;
__constant__ unsigned int ysize;
__constant__ int bc_type[4];
__constant__ double bc_values[4];
__constant__ int A;
__constant__ int B;
__constant__ int C;

__device__ __forceinline__ int index(int x, int y) { return x + y * xsize; }

__device__ void destroy_needle(int *is_needle, Needle needle) {
  for (int i = needle.x0; i <= needle.xf; ++i) {
    /* int idx = i + needle.y0*xsize; */
    /* int idx = index(i, needle.y0); */

    is_needle[index(i, needle.y0)] = 0;

    int wd = (int)sqrt(2 * needle.rad * (needle.xf - i) * dx) / dy;
    /* wd = wd < 3 ? wd : 3; */
    for (int j = 1; j <= wd; ++j) {
      /* is_needle[idx + j*xsize] = 0; */
      /* is_needle[idx - j*xsize] = 0; */
      is_needle[index(i, needle.y0 + j)] = 0;
      is_needle[index(i, needle.y0 - j)] = 0;
    }
  }
}

__device__ double calculate_line_integral(double *grid, int xb, int xf, int yb,
                                          int yf, int y_nb, int y_nt) {
  double dxu;
  double dyu;
  double flux_acc = 0.0;
  double xfactor = 1.0 / (2.0 * dx);
  double yfactor = 1.0 / (2.0 * dy);

  // Bottom contour line
  // Edge case bottom contour
  dyu = yfactor * (grid[index(xb - 1, yb)] + grid[index(xb, yb)] -
                   grid[index(xb - 1, yb - 1)] - grid[index(xb, yb - 1)]);
  flux_acc -= 0.5 * dyu * dx;
  for (int i = xb; i < xf; ++i) {
    dyu = yfactor * (grid[index(i, yb)] + grid[index(i + 1, yb)] -
                     grid[index(i, yb - 1)] - grid[index(i + 1, yb - 1)]);
    flux_acc -= dyu * dx;
  }
  dyu = yfactor * (grid[index(xf, yb)] + grid[index(xf + 1, yb)] -
                   grid[index(xf, yb - 1)] - grid[index(xf + 1, yb - 1)]);
  flux_acc -= 0.5 * dyu * dx;

  // Front contour line
  dxu = xfactor * (grid[index(xf + 1, yb - 1)] + grid[index(xf + 1, yb)] -
                   grid[index(xf, yb - 1)] - grid[index(xf, yb)]);
  flux_acc += 0.5 * dxu * dy;
  for (int i = yb; i < yf; ++i) {
    dxu = xfactor * (grid[index(xf + 1, i)] + grid[index(xf + 1, i + 1)] -
                     grid[index(xf, i)] - grid[index(xf, i + 1)]);
    flux_acc += dxu * dy;
  }
  dxu = xfactor * (grid[index(xf + 1, yf)] + grid[index(xf + 1, yf + 1)] -
                   grid[index(xf, yf)] - grid[index(xf, yf + 1)]);
  flux_acc += 0.5 * dxu * dy;

  // Top contour line
  dyu = yfactor * (grid[index(xb - 1, yf + 1)] + grid[index(xb, yf + 1)] -
                   grid[index(xb - 1, yf)] - grid[index(xb, yf)]);
  flux_acc += 0.5 * dyu * dx;
  for (int i = xb; i < xf; ++i) {
    dyu = yfactor * (grid[index(i, yf + 1)] + grid[index(i + 1, yf + 1)] -
                     grid[index(i, yf)] - grid[index(i + 1, yf)]);
    flux_acc += dyu * dx;
  }
  dyu = yfactor * (grid[index(xf, yf + 1)] + grid[index(xf + 1, yf + 1)] -
                   grid[index(xf, yf)] - grid[index(xf + 1, yf)]);
  flux_acc += 0.5 * dyu * dx;

  // Back contour line
  // Split in two, due to do needle
  for (int i = y_nt; i < yf; ++i) {
    dxu = xfactor * (grid[index(xb, i)] + grid[index(xb, i + 1)] -
                     grid[index(xb - 1, i)] - grid[index(xb - 1, i + 1)]);
    flux_acc -= dxu * dy;
  }
  dxu = xfactor * (grid[index(xb, yf)] + grid[index(xb, yf + 1)] -
                   grid[index(xb - 1, yf)] - grid[index(xb - 1, yf + 1)]);
  flux_acc -= 0.5 * dxu * dy;

  dxu = xfactor * (grid[index(xb, yb - 1)] + grid[index(xb, yb)] -
                   grid[index(xb - 1, yb - 1)] - grid[index(xb - 1, yb)]);
  flux_acc -= 0.5 * dxu * dy;
  for (int i = yb; i < y_nb; ++i) {
    dxu = xfactor * (grid[index(xb, i)] + grid[index(xb, i + 1)] -
                     grid[index(xb - 1, i)] - grid[index(xb - 1, i + 1)]);
    flux_acc -= dxu * dy;
  }

  return flux_acc;
}

__device__ double calculate_surface_integral(double *grid, int *is_needle,
                                             int xb, int xf, int yb, int yf) {
  double xfactor = 1.0 / (2.0 * dx);
  double integral = 0.0;
  double dxu;

  for (int y = yb; y <= yf; ++y) {
    for (int x = xb; x <= xf; ++x) {
      if (is_needle[index(x, y)] == 0) {
        dxu = xfactor * (grid[index(x + 1, y)] - grid[index(x - 1, y)]);
        integral += dxu * dx * dy;
      }
    }
  }

  return integral;
}

__device__ double calculate_flux_intensity_factor(double *grid, int *is_needle,
                                                  Needle needle) {
  double factor = 1 / (4 * sqrt((A + needle.r + 1 / 2) * dx));

  int xb = needle.xf - A;
  int xf = needle.xf + B;

  int y_nt, y_nb;
  for (int i = 0; i < ysize; ++i) {
    if (is_needle[index(xb, needle.y0 + i)] == 0) {
      y_nt = needle.y0 + i - 1;
      y_nb = needle.y0 - i + 1;
      break;
    }
  }

  int yf = y_nt + C;
  int yb = y_nb - C;

  double line = calculate_line_integral(grid, xb, xf, yb, yf, y_nb, y_nt);
  double surface = calculate_surface_integral(grid, is_needle, xb, xf, yb, yf);
  /* printf("\n\nIntegral values line: %f, surface: %f\n\n", line, surface); */
  surface = (needle.vel * surface) / d;

  return (line + surface) * factor;
}

__device__ void grow_needle(double *grid, int *is_needle, Needle *needle) {
  double fif = calculate_flux_intensity_factor(grid, is_needle, *needle);

  /* printf("\nFIF calcuated: %f\n\n", fif); */

  double fif_factor = fif * fif * 2 * d * d;
  double vel = pow(fif_factor, 0.6666666666);
  double rad = pow(fif_factor, -0.3333333333);

  if (rad < needle->rad) {
    destroy_needle(is_needle, *needle);
  }

  /* printf("\nNeedle rad: %f, vel: %f, r: %f\n\n", needle->rad, needle->vel,
   * needle->r); */
  /* printf("\nNew needle rad: %f, vel: %f, r: %f\n\n", vel, rad, 0.0); */

  needle->r += (vel * dt) / dx;
  needle->rad = rad;
  needle->vel = vel;

  if (needle->r > 1.0) {
    needle->xf += 1;
    needle->r -= 1.0;
  }
}

__global__ void print_device_constants() {
  int idx = threadIdx.x;

  if (idx == 0) {
    printf("device_gridFo: %f\n", gridFo);
    printf("device_d: %f\n", d);
    printf("device_dx: %f\n", dx);
    printf("device_dy: %f\n", dy);
    printf("device_xsize: %d\n", xsize);
    printf("device_ysize: %d\n", ysize);

    for (int i = 0; i < 4; ++i) {
      printf("BC_TYPE[%d]: %d, BC_VAL[%d], %f\n", i, bc_type[i], i,
             bc_values[i]);
    }
  }
}

__global__ void print_needles(Needle *needles, int size) {
  int idx = threadIdx.x;

  if (idx == 0) {

    for (int i = 0; i < size; ++i) {
      printf("Device: (x0: %d, y0: %d), (xf: %d, yf: %d), rad: %f, vel: %f, r: "
             "%f\n",
             needles[i].x0, needles[i].y0, needles[i].xf, needles[i].yf,
             needles[i].rad, needles[i].vel, needles[i].r);
    }
  }
}

__global__ void grow_needles(double *grid, int *is_needle, Needle *needles,
                             int size) {
  // Assuming needle is alinged along the x-axis
  int idx = threadIdx.x + blockDim.x * blockIdx.x;

  if (idx < size) {
    grow_needle(grid, is_needle, &needles[idx]);
  }
}

__device__ void set_needle(double *grid, int *is_needle, Needle needle) {
  for (int i = needle.x0; i <= needle.xf; ++i) {
    /* int idx = i + needle.y0*xsize; */
    int idx = index(i, needle.y0);

    grid[idx] = 0;
    is_needle[idx] = 1;

    int wd = (int)sqrt(2 * needle.rad * (needle.xf - i) * dx) / dy;
    /* wd = wd < 3 ? wd : 3; */
    for (int j = 1; j <= wd; ++j) {
      /* grid[idx + j*xsize] = 0; */
      /* grid[idx - j*xsize] = 0; */
      /* is_needle[idx + j*xsize] = 1; */
      /* is_needle[idx - j*xsize] = 1; */
      int idx_t = index(i, needle.y0 + j);
      int idx_b = index(i, needle.y0 - j);
      grid[idx_t] = 0;
      grid[idx_b] = 0;
      is_needle[idx_t] = 1;
      is_needle[idx_b] = 1;
    }
  }
}

__device__ void shift_row(double *grid, int *is_needle, int start_index) {
  thrust::copy(thrust::device, grid + start_index + 1,
               grid + start_index + xsize, grid + start_index);
  thrust::copy(thrust::device, is_needle + start_index + 1,
               is_needle + start_index + xsize, is_needle + start_index);
}

__global__ void set_needles(double *grid, int *is_needle, Needle *needles,
                            int size) {
  // Assuming needly only aligned along the x-axis
  int idx = threadIdx.x + blockDim.x * blockIdx.x;

  if (idx < size) {
    set_needle(grid, is_needle, needles[idx]);
  }
}

__global__ void shift_domain(double *grid, int *is_needle, Needle *needles,
                             double pos) {
  if (needles[0].xf > pos) {
    int idx = threadIdx.x + blockDim.x * blockIdx.x;

    if (idx < ysize) {
      shift_row(grid, is_needle, idx * xsize);
      if (idx == 0) {
        needles[0].xf = needles[0].xf - 1;
      }
    }
  }
}

__global__ void explicit_euler(double *grid, double *ogrid) {
  // I do not see any reason why the solver shoulde be aware of the needle
  int idx = threadIdx.x + blockDim.x * blockIdx.x;
  if (idx < xsize * ysize) {
    int idxl = idx - 1;
    int idxr = idx + 1;
    int idxt = idx + xsize;
    int idxb = idx - xsize;
    double nval;
    if (idx == 0) {
      // Left Bottom corner
      int bc = bc_type[LW] + bc_type[BW];
      if (bc == 2 * FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxr] + 2 * dx * bc_values[LW] / d +
                                  2 * grid[idxt] + 2 * dy * bc_values[BW] / d -
                                  4 * grid[idx]);
      } else if (bc == FLUX) {
        if (bc_type[LW] == 0) {
          nval = bc_values[LW];
        } else {
          nval = bc_values[BW];
        }
      } else {
        nval = (bc_values[LW] + bc_values[BW]) / 2;
      }
    } else if (idx == xsize - 1) {
      // Right bottom corner
      int bc = bc_type[RW] + bc_type[BW];
      if (bc == 2 * FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxl] - 2 * dx * bc_values[RW] / d +
                                  2 * grid[idxt] + 2 * dy * bc_values[BW] / d -
                                  4 * grid[idx]);
      } else if (bc == FLUX) {
        if (bc_type[RW] == 0) {
          nval = bc_values[RW];
        } else {
          nval = bc_values[BW];
        }
      } else {
        nval = (bc_values[RW] + bc_values[BW]) / 2;
      }
    } else if (idx == (ysize - 1) * xsize) {
      // Left top corner
      int bc = bc_type[LW] + bc_type[TW];
      if (bc == 2 * FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxr] + 2 * dx * bc_values[LW] / d +
                                  2 * grid[idxb] - 2 * dy * bc_values[TW] / d -
                                  4 * grid[idx]);
      } else if (bc == FLUX) {
        if (bc_type[LW] == 0) {
          nval = bc_values[LW];
        } else {
          nval = bc_values[TW];
        }
      } else {
        nval = (bc_values[LW] + bc_values[TW]) / 2;
      }
    } else if (idx == xsize * ysize - 1) {
      // Right top corner
      int bc = bc_type[RW] + bc_type[TW];
      if (bc == 2 * FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxl] - 2 * dx * bc_values[RW] / d +
                                  2 * grid[idxb] - 2 * dy * bc_values[TW] / d -
                                  4 * grid[idx]);
      } else if (bc == FLUX) {
        if (bc_type[LW] == 0) {
          nval = bc_values[RW];
        } else {
          nval = bc_values[TW];
        }
      } else {
        nval = (bc_values[RW] + bc_values[TW]) / 2;
      }
    } else if (idx < xsize) {
      // Bottom wall
      if (bc_type[BW] == FLUX) {
        nval =
            grid[idx] + gridFo * (grid[idxr] + grid[idxl] + 2 * grid[idxt] +
                                  2 * dy * bc_values[BW] / d - 4 * grid[idx]);
      } else {
        nval = bc_values[BW];
      }
    } else if ((idx + 1) % xsize == 0) {
      // Right wall
      if (bc_type[RW] == FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxl] - 2 * dx * bc_values[RW] / d +
                                  grid[idxt] + grid[idxb] - 4 * grid[idx]);
      } else {
        nval = bc_values[RW];
      }
    } else if (idx % xsize == 0) {
      // Left Wall
      if (bc_type[LW] == FLUX) {
        nval =
            grid[idx] + gridFo * (2 * grid[idxr] + 2 * dx * bc_values[LW] / d +
                                  grid[idxt] + grid[idxb] - 4 * grid[idx]);
      } else {
        nval = bc_values[LW];
      }
    } else if (idx > (ysize - 1) * xsize) {
      // Top wall
      if (bc_type[TW] == FLUX) {
        nval =
            grid[idx] + gridFo * (grid[idxr] + grid[idxl] + 2 * grid[idxb] -
                                  2 * dy * bc_values[TW] / d - 4 * grid[idx]);
      } else {
        nval = bc_values[TW];
      }
    } else {
      // Inner domain
      nval = grid[idx] + gridFo * (grid[idxl] + grid[idxr] + grid[idxt] +
                                   grid[idxb] - 4 * grid[idx]);
    }
    ogrid[idx] = nval;
  }
}

void run_transient_CUDA(double *d_grid, double *d_ogrid, int size) {
  int nthreads = BLOCK_SIZE;
  int nblocks = ceil((double)size / (double)nthreads);
  explicit_euler<<<nblocks, nthreads>>>(d_grid, d_ogrid);
}

void set_parabolic_needle(double *d_grid, int *d_gneedle, Needle *d_needles,
                          int size) {
  int block_size = 32;
  int nblocks =
      ceil(static_cast<double>(size) / static_cast<double>(block_size));
  set_needles<<<nblocks, block_size>>>(d_grid, d_gneedle, d_needles, size);
  /* print_needles<<<1, 32>>>(needles, size); */
  /* cudaDeviceSynchronize(); */
}

void run_grow_needles(double *d_grid, int *d_gneedle, Needle *d_needles,
                      int size) {
  int block_size = 32;
  int nblocks =
      ceil(static_cast<double>(size) / static_cast<double>(block_size));
  grow_needles<<<nblocks, block_size>>>(d_grid, d_gneedle, d_needles, size);
}

void run_shift_domain(double *d_grid, int *d_gneedle, Needle *d_needles,
                      int ysize, double pos) {
  int block_size = 32;
  int nblocks =
      ceil(static_cast<double>(ysize) / static_cast<double>(block_size));
  shift_domain<<<nblocks, block_size>>>(d_grid, d_gneedle, d_needles, pos);
}

void init_CUDA(double _gridFo, double _d, double _dx, double _dy, double _dt,
               unsigned int _xsize, unsigned int _ysize, int _A, int _B, int _C,
               int _bc_type[4], double _bc_values[4]) {
  printf("INIT_CUDA Function\n");
  checkCudaErrors(cudaMemcpyToSymbol(gridFo, &_gridFo, sizeof(double)));
  checkCudaErrors(cudaMemcpyToSymbol(d, &_d, sizeof(double)));
  checkCudaErrors(cudaMemcpyToSymbol(dx, &_dx, sizeof(double)));
  checkCudaErrors(cudaMemcpyToSymbol(dy, &_dy, sizeof(double)));
  checkCudaErrors(cudaMemcpyToSymbol(dt, &_dt, sizeof(double)));
  checkCudaErrors(cudaMemcpyToSymbol(xsize, &_xsize, sizeof(unsigned int)));
  checkCudaErrors(cudaMemcpyToSymbol(ysize, &_ysize, sizeof(unsigned int)));
  checkCudaErrors(cudaMemcpyToSymbol(A, &_A, sizeof(int)));
  checkCudaErrors(cudaMemcpyToSymbol(B, &_B, sizeof(int)));
  checkCudaErrors(cudaMemcpyToSymbol(C, &_C, sizeof(int)));
  checkCudaErrors(cudaMemcpyToSymbol(bc_type, _bc_type, 4 * sizeof(int)));
  checkCudaErrors(
      cudaMemcpyToSymbol(bc_values, _bc_values, 4 * sizeof(double)));
  /* print_device_constants<<<1, 32>>>(); */
  /* cudaDeviceSynchronize(); */
}

void copy_to_device_memory(void *d_ptr, void *h_ptr, int bytes) {
  /* printf("\n\nCopy to device memory called.\n\n"); */
  checkCudaErrors(cudaMemcpy(d_ptr, h_ptr, bytes, cudaMemcpyHostToDevice));
}

void copy_to_host_memory(void *h_ptr, void *d_ptr, int bytes) {
  /* printf("\n\nCopy to host memory called.\n\n"); */
  checkCudaErrors(cudaMemcpy(h_ptr, d_ptr, bytes, cudaMemcpyDeviceToHost));
}

void allocate_device_memory(void **d_ptr_ptr, int bytes) {
  /* printf("\n\nAllocate to device memory called.\n\n"); */
  checkCudaErrors(cudaMalloc(d_ptr_ptr, bytes));
}

void cuda_synchronize_device() { checkCudaErrors(cudaDeviceSynchronize()); }

void cuda_release_memory(void *d_ptr) { checkCudaErrors(cudaFree(d_ptr)); }
