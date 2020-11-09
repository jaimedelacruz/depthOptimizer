#ifndef READMODELHPP
#define READMODELHPP

#include "Arrays.hpp"
#include "model.hpp"
//#include "io_hdf5.hpp"

namespace ml{

  // **************************************************************************** //

  template<typename T, size_t N>
  void readAndFill(std::vector<hsize_t> const& dims, std::string const &vname, int const idx, Model<T> &model)
  {
    T const nElements = ch5n::getTotalSize(dims);
    T* __restrict__ buffer = new T [nElements]();
    
    
    
  }
  
  
  
  // **************************************************************************** //

  template <typename T>
  ml::Model<T> readModel(std::string const &filename, int &nt, int const tstep = 0){

    //static std::vector<std::string> const VNAMES = {"z", "y", "x", "electron_density", "hydrogen_populations", "temperature", "velocity_z", "velocity_turbulent","density", "B_x", "B_y", "B_z"};
    
    static std::vector<std::string> const VNAMES = {"z",  "temperature", "velocity_z", "velocity_turbulent", "B_x", "B_y", "B_z", "electron_density", "density", "hydrogen_populations"};
    constexpr static const int vidx[10] = {0,1,2,3,4,5,6,7,8,9};
    
    ch5n::ch5 IO(filename, ch5n::r);
    int ny=0, nx=0, ndep=0, nhydrogen = 0;

    std::vector<std::string> var_names = IO.getVarNames();
    std::vector<hsize_t> dims = IO.getVariableDimensions("hydrogen_populations");

    nt = 1;
    if(dims.size() == 5){
      nt = dims[0];
      nhydrogen = dims[1];
      ny = dims[2];
      nx = dims[3];
      ndep = dims[4];
    }else if(dims.size() == 4){
      nhydrogen = dims[0];
      ny = dims[1];
      nx = dims[2];
      ndep = dims[3];
    }else{
      std::cerr<<"ml::readModel: ERROR, number of dimensions for hydrogen must be 4 or 5"<<std::endl;
      exit(1);
    }

    ml::Model<T> model(ny, nx, ndep, nhydrogen);

    int const nvar = VNAMES.size();
    for(int ii=0; ii<nvar; ++ii){
      for(auto &ot: var_names){
	
	if(VNAMES[ii] == ot){
	  std::cerr<<"ml::readModel: reading ["<<ot<<"]\n";

	  std::vector<hsize_t> iDim = IO.getVariableDimensions(ot);
	  size_t nElements = 0;

	  // --- Read variable --- //
	  IO.read_raw_external<T>(ot, nElements, &model.d(vidx[ii],0,0,0), tstep);
	  
	}

		     
      }

    }
    

    
    return model;
  }
}

#endif

