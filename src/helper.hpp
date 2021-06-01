#ifndef HELPER_H
#define HELPER_H

#include "data_types.hpp"

struct peclet_params {
  double omega;
};

double peclet_number(double, void *);

int is_device_pointer(const void *);

void swap_device_pointers(double *&, double *&);

// GPU Related methods
void init_CUDA(double, double, double, double, double, unsigned int,
               unsigned int, int, int, int, int[4], double[4]);
void run_transient_CUDA(double *, double *, int);
void copy_to_device_memory(void *, void *, int);
void copy_to_host_memory(void *, void *, int);
void allocate_device_memory(void **, int);
void set_parabolic_needle(double *, int *, Needle *, int);
void run_grow_needles(double *, int *, Needle *, int);
void run_shift_domain(double *, int *, Needle *, int, double);
void cuda_synchronize_device();
void cuda_release_memory(void *ptr);

#endif
