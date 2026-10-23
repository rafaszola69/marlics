#ifndef GPU_H
#define GPU_H
// Entry point to the CUDA code (src/gpu/), so that the rest of the program needs no CUDA headers.
// The CPU-only build links gpu_stub.cpp instead, which reports that the GPU is not available.
#include "integrator.h"

//Integrator running on the GPU, or NULL if the integrator has no GPU version.
Integrator * make_gpu_integrator(GEOMETRY * geometry, const struct Simulation_Parameters * sim_param);

#endif
