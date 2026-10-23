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
  bulk_force(lc),
  Nx(lc->Nx),
  Ny(lc->Ny),
  Nz(lc->Nz),
  dx_1(1/lc->dx),
  dy_1(1/lc->dy),
  dz_1(1/lc->dz)
{

  std::cout << "Nx=" << Nx << " \n";
  std::cout << "Ny=" << Ny << " \n";
  std::cout << "Nz=" << Nz << " \n\n";


};


  
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

