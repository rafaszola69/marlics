#ifndef WEDGE_SLAB2_
#define WEDGE_SLAB2_


class Geometry_Wedge2 : public GEOMETRY
{
  
  
 public:


  virtual void fill_ki(double *, const double * ) const ;
    
  Geometry_Wedge2(const struct Simulation_Parameters *);
  ~Geometry_Wedge2(void) {};
  virtual void Energy_calc(double * k_i, const double * Qij) const ;
  virtual int backwardPTP(double *, const double * ) const  ;


  
   private:
	int Nz2;
  double wTheta;
  virtual int * fill_point_type( void ) const ;  
  

  //void check_bulk_limits(int i, int j, int k) const ;
  //void check_surface_limits(int i, int j, int k) const ;
  
};

#endif
