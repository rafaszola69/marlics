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

//Bulk energy density at site s, averaged over the 8 octants of one-sided differences.
//Its exact gradient is the compact stencil used by fill_ki (3-point second derivatives and
//4-corner mixed derivatives), so the time integrators decrease it, unlike a central-difference energy.
double GEOMETRY::bulk_energy_octant(const double * Qij, int s, int xm, int xp, int ym, int yp, int zm, int zp, const double v[3]) const
{
  const int nx[2]={xm,xp}, ny[2]={ym,yp}, nz[2]={zm,zp};
  const double sg[2]={-1.,1.};
  double QN[5], dQ[15], f=0;

  for(int ll=0; ll<=4;ll++) QN[ll]=Qij[5*s+ll];
  for(int a=0; a<2; a++)
    for(int b=0; b<2; b++)
      for(int c=0; c<2; c++)
        {
          for(int ll=0; ll<=4;ll++)
            {
              dQ[ll]   =sg[a]*(Qij[5*nx[a]+ll]-QN[ll])*dx_1;
              dQ[5+ll] =sg[b]*(Qij[5*ny[b]+ll]-QN[ll])*dy_1;
              dQ[10+ll]=sg[c]*(Qij[5*nz[c]+ll]-QN[ll])*dz_1;
            }
          f+=bulk_energy.energy_calculation(QN,dQ,v);
        }
  return 0.125*f;
}

int GEOMETRY::neighbor(int i, int j, int k, int axis, int dir) const
{
  int c[3]={i,j,k};
  const int n[3]={Nx,Ny,Nz};
  const int s=(k*Ny+j)*Nx+i;

  c[axis]+=dir;
  if(c[axis]<0 || c[axis]>=n[axis])
    {
      if(!is_periodic(axis)) return s;
      c[axis]=(c[axis]+n[axis])%n[axis];
    }
  const int t=(c[2]*Ny+c[1])*Nx+c[0];
  return point_type[t]==0 ? s : t;
}

int GEOMETRY::surface_site(int i, int j, int k, double v[3], double & area) const
{
  std::cout << "The " << geometry_name << " geometry does not define its surface sites.\nAborting the program.\n";
  exit(1);
}

//Every liquid crystal site carries bulk energy over its volume dV; surface sites also carry
//the anchoring energy over the area they represent. Bulk sites therefore see exactly the
//energy whose gradient fill_ki uses, including the links to surface sites.
void GEOMETRY::Energy_calc(double * k_i, const double * Qij) const
{
  #pragma omp for schedule(static)
  for(int s=0; s<Nx*Ny*Nz; s++)
    {
      const int i=s%Nx, j=(s/Nx)%Ny, k=s/(Nx*Ny);
      double v[3]={0.,0.,0.};
      const double dQ[15]={0.};

      k_i[2*s+0]=0.;
      k_i[2*s+1]=0.;
      if(point_type[s]==0) continue;

      k_i[2*s+0]=dV*bulk_energy_octant(Qij,s,
                                       neighbor(i,j,k,0,-1),neighbor(i,j,k,0,1),
                                       neighbor(i,j,k,1,-1),neighbor(i,j,k,1,1),
                                       neighbor(i,j,k,2,-1),neighbor(i,j,k,2,1),v);
      if(point_type[s]>=2)
        {
          double area;
          const int bc=surface_site(i,j,k,v,area);
          k_i[2*s+1]=area*bc_conditions[bc]->energy_calculation(&Qij[5*s],dQ,v);
        }
    }
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
};

