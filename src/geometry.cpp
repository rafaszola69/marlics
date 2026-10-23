#include "geometry.h"
#include "driver.h"
#include "energy.h"
#include "energy_ldg.h"
#include "boundary.h"
#include "boundary_strong.h"
#include "boundary_rp.h"
#include "boundary_fournier_galatola.h"
#include "boundary_homeotropic.h"
#include <stdio.h> 
#include <vector>


double Pi=3.14159265359;
//Base class to latter geometry definition
GEOMETRY::GEOMETRY(const struct Simulation_Parameters * lc) :
  bulk_energy(lc),
  Nx(lc->Nx),
  Ny(lc->Ny),
  Nz(lc->Nz),
  dx_1(1/lc->dx),
  dy_1(1/lc->dy),
  dz_1(1/lc->dz),
  dV(lc->dx*lc->dy*lc->dz),
  dA(lc->dx*lc->dy)
{

  std::cout << "Nx=" << Nx << " \n";
  std::cout << "Ny=" << Ny << " \n";
  std::cout << "Nz=" << Nz << " \n\n";


};

int GEOMETRY::surface_site(int i, int j, int k, double v[3], double & area) const
{
  std::cout << "The " << geometry_name << " geometry does not define its surface sites.\nAborting the program.\n";
  exit(1);
}

//Default: no derivatives at surface sites (used only by geometries that override fill_ki).
void GEOMETRY::surface_stencil(int i, int j, int k, int nm[3], int np[3]) const
{
  const int s=(k*Ny+j)*Nx+i;
  for(int a=0; a<3; a++) nm[a]=np[a]=s;
}

void GEOMETRY::build_lattice( void )
{
  surface.clear();
  for(int s=0; s<Nx*Ny*Nz; s++)
    {
      if(point_type[s]<2) continue;
      const int i=s%Nx, j=(s/Nx)%Ny, k=s/(Nx*Ny);
      SurfaceSite S;
      S.site=s;
      S.bc=surface_site(i,j,k,S.v,S.area);
      surface_stencil(i,j,k,S.nm,S.np);
      for(int a=0; a<3; a++)
        S.w[a]= (S.nm[a]!=s && S.np[a]!=s) ? 0.5 : 1.0;
      surface.push_back(S);
    }

  boundary_params.clear();
  for(BOUNDARY * bc : bc_conditions) boundary_params.push_back(bc->params());
}

LatticeView GEOMETRY::view( void ) const
{
  LatticeView L;
  L.Nx=Nx; L.Ny=Ny; L.Nz=Nz;
  for(int a=0; a<3; a++) L.periodic[a]=is_periodic(a);
  L.d_1[0]=dx_1; L.d_1[1]=dy_1; L.d_1[2]=dz_1;
  L.dV=dV;
  L.point_type=point_type;
  L.n_surface=surface.size();
  L.surface=surface.data();
  L.bulk=bulk_energy.bulk();
  L.n_bc=boundary_params.size();
  L.bc=boundary_params.data();
  return L;
}

//Called inside a parallel region. Dynamic schedule: bulk sites, the expensive ones, are unevenly
//spread (a droplet fills only the middle of the box). Bulk sites first (surface and outside sites are zeroed), then
//the surface sites; the implicit barrier between the two loops keeps them from overlapping.
void GEOMETRY::fill_ki(double * k_i, const double * Qij) const
{
  const LatticeView L=view();

  #pragma omp for schedule(dynamic,Nx)
  for(int s=0; s<Nx*Ny*Nz; s++)
    {
      if(point_type[s]==1) bulk_site_force(L,Qij,s,&k_i[5*s]);
      else for(int ll=0; ll<=4; ll++) k_i[5*s+ll]=0.;
    }

  #pragma omp for schedule(static)
  for(int n=0; n<L.n_surface; n++)
    surface_site_force(L,Qij,L.surface[n],&k_i[5*L.surface[n].site]);
}

//Every liquid crystal site carries bulk energy over its volume dV; surface sites also carry
//the anchoring energy over the area they represent. Bulk sites therefore see exactly the
//energy whose gradient fill_ki uses, including the links to surface sites.
void GEOMETRY::Energy_calc(double * k_i, const double * Qij) const
{
  const LatticeView L=view();

  #pragma omp for schedule(static)
  for(int s=0; s<Nx*Ny*Nz; s++)
    {
      k_i[2*s+0]= point_type[s]==0 ? 0. : dV*bulk_site_energy(L,Qij,s);
      k_i[2*s+1]=0.;
    }

  #pragma omp for schedule(static)
  for(int n=0; n<L.n_surface; n++)
    k_i[2*L.surface[n].site+1]=surface_site_energy(L,Qij,L.surface[n]);
}

void GEOMETRY::boundary_init(struct Simulation_Parameters * sim_param)
{
    std::cout <<"\nInitiating boundary conditions:\n\n";
    int ii=0;
    std::string anc_type;
    for (ii=0; ii< number_of_boundaries; ii++)
    {
        try
        {
            anc_type = sim_param->anchoring_type.at(ii);
            std::cout << "Boundary " << ii << ":" << anc_type <<".\n";
        }
        catch(std::out_of_range dummy_var)
        {
            std::cout<< "You must define boundaries "<< boundary_needed_to_be_defined << " in the " << geometry_name << " geometry.\nPlease review your input file.\nAborting the program.\n\n";
            exit(0);    
        }

        if (strcasecmp(anc_type.c_str(),"planar_strong") == 0)
        {
            // Salvar os valores originais (se existirem)
            double original_Wo1 = 0.0;
            double original_Wo2 = 0.0;
            bool has_Wo1 = false;
            bool has_Wo2 = false;

            try {
                original_Wo1 = sim_param->Wo1.at(ii);
                has_Wo1 = true;
            } catch(std::out_of_range dummy_var) {}

            try {
                original_Wo2 = sim_param->Wo2.at(ii);
                has_Wo2 = true;
            } catch(std::out_of_range dummy_var) {}

            // Definir valores altos para Wo1 e Wo2 para este ancoramento específico
            sim_param->Wo1[ii] = 10000.0;
            sim_param->Wo2[ii] = 50000.0;

            // Criar o objeto Boundary_Fg com os valores modificados
            bc_conditions[ii] = new Boundary_Fg(sim_param, ii);

            // Restaurar os valores originais (se necessário)
            if (has_Wo1) {
                sim_param->Wo1[ii] = original_Wo1;
            } else {
                sim_param->Wo1.erase(ii);
            }

            if (has_Wo2) {
                sim_param->Wo2[ii] = original_Wo2;
            } else {
                sim_param->Wo2.erase(ii);
            }

            std::cout << "Using Planar Strong Boundary condition (implemented as FG with high Wo1=10.0 and Wo2=50.0)" << std::endl;
        }
        else if (strcasecmp(anc_type.c_str(),"strong") == 0 || strcasecmp(anc_type.c_str(),"fixed") == 0)
        {
            bc_conditions[ii] = new Strong_Boundary(sim_param, ii);
        }
        else if (strcasecmp(anc_type.c_str(),"rp") == 0 || strcasecmp(anc_type.c_str(),"Rapini-Papoular") == 0)
        {
            bc_conditions[ii] = new Boundary_Rp(sim_param, ii);
        }
        else if (strcasecmp(anc_type.c_str(),"fg") == 0 || strcasecmp(anc_type.c_str(),"fournier-galatola") == 0)
        {
            bc_conditions[ii] = new Boundary_Fg(sim_param, ii);
        }
        else if (strcasecmp(anc_type.c_str(),"homeotropic") == 0)
        {
            bc_conditions[ii] = new Boundary_Homeotropic(sim_param, ii);
        }
        else
        {
            std::cout<<"There is no anchoring type named " << anc_type <<".\nPlease check your input file 'anchoring type' fields.\n\n Aborting the program.'\n";
            exit(0);
        } 
    }  
    build_lattice();
};

