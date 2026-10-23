#ifndef WEDGE_SLAB_
#define WEDGE_SLAB_


class Geometry_Wedge : public GEOMETRY
{
  
  
 public:


  virtual void fill_ki(double *, const double * ) const ;
    
  Geometry_Wedge(const struct Simulation_Parameters *);
  ~Geometry_Wedge(void) {};
  virtual int surface_site(int, int, int, double [3], double &) const ;
  virtual bool is_periodic(int axis) const { return axis<2; }
  virtual int backwardPTP(double *, const double * ) const  ;


  
   private:
	int Nz2;
  virtual int * fill_point_type( void ) const ;  
  

  //void check_bulk_limits(int i, int j, int k) const ;
  //void check_surface_limits(int i, int j, int k) const ;
  
};

#endif
