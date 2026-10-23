// GPU versions of the integrators: the same templates as the CPU ones, run by CudaExec.
#include "exec_cuda.h"
#include "fire_integrator_impl.h"
#include "integrator_dp5_impl.h"
#include "integrator_rk2_impl.h"
#include "integrator_euler_impl.h"
#include "integrator_factory.h"
#include "gpu.h"

Integrator * make_gpu_integrator(GEOMETRY * geometry, const struct Simulation_Parameters * sim_param)
{
  return make_integrator<CudaExec>(geometry, sim_param);
}
