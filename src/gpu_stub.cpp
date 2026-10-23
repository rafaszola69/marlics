// CPU-only build: no GPU support (see "make gpu").
#include "driver.h"
#include "gpu.h"
#include <iostream>
#include <cstdlib>

Integrator * make_gpu_integrator(GEOMETRY * geometry, const struct Simulation_Parameters * sim_param)
{
  std::cout << "This Marlics was compiled without GPU support (device gpu).\n"
            << "Build it with \"make gpu\" or use \"device cpu\".\nAborting the program.\n";
  exit(2);
}
