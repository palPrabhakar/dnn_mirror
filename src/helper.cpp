#include <cmath>
#include <iostream>
/* #include <cuda_runtime.h> */
#include "dnn.hpp"
#include "helper.hpp"

#if GPU
#include "helper_cuda.h"
#endif

/* int is_device_pointer(const void *ptr) */
/* { */
/*   int is_device_ptr = 0; */
/*   cudaPointerAttributes attributes; */

/*   checkCudaErrors(cudaPointerGetAttributes(&attributes, ptr)); */

/*   if(attributes.devicePointer != NULL) */
/*   { */
/*     is_device_ptr = 1; */
/*     std::cout<<"\nYes *ptr is device pointer\n"; */
/*   } */
/*   else { */
/*     std::cout<<"\n LOL!! *ptr is not device pointer\n"; */
/*   } */

/*   return is_device_ptr; */
/* } */

double peclet_number(double x, void *params) {
  double omega = ((peclet_params *)params)->omega;
  return sqrt(M_PI * x) * exp(x) * erfc(sqrt(x)) - omega;
}

void swap_device_pointers(double *&p1, double *&p2) {
  double *temp;
  temp = p1;
  p1 = p2;
  p2 = temp;
}
