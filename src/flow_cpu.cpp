// FFTs of the CPU backflow (pocketfft). Kept out of the headers: old nvcc versions cannot compile it.
#define POCKETFFT_NO_MULTITHREADING
#include "third_party/pocketfft_hdronly.h"
#include "flow_cpu.h"

void FlowCpu::fft(bool forward)
{
  const pocketfft::shape_t shape={(size_t)L.Nz,(size_t)L.Ny,(size_t)L.Nx};
  const ptrdiff_t c=sizeof(std::complex<double>);
  const pocketfft::stride_t stride={(ptrdiff_t)L.Ny*L.Nx*c,(ptrdiff_t)L.Nx*c,c};
  pocketfft::shape_t axes={0,1,2};
  double scale=1./n;
  if(walls)
    {
      //shape axis of x, y, z: 2, 1, 0
      axes={(size_t)(2-walls->W.pb),(size_t)(2-walls->W.pa)};
      scale=1./n_plane;
    }
  for(int a=0; a<3; a++)
    pocketfft::c2c(shape,stride,stride,axes,forward,fk[a].data(),fk[a].data(),forward ? 1. : scale);
}
