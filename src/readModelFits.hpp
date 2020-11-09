#ifndef RMODELFITS
#define RMODELFITS

#include <array>
#include <string>
#include <iostream>
#include <fitsio.h>
#include <fstream>
#include <sys/stat.h>
#include <typeinfo>


#include "Arrays.hpp"
#include "model.hpp"

namespace fits{

  // *************************************************** //

  template<typename T> size_t getTotalSize(std::vector<T> const& dim)
  {
    size_t res = 1;
    int const nEl = dim.size();
    for(int ii=0; ii<nEl;++ii) res *= dim[ii];
    return res;
  }
  
  // *************************************************** //

  inline bool file_exists(std::string const &name) {
    std::ifstream f(name.c_str(),std::ifstream::in);
    if (f.good()) {
      f.close();
      return true;
    } else {
      f.close();
      return false;
    }
  }

  // *************************************************** //

  template<typename T>
  inline std::string dim2string(std::vector<T> const& dim)
  {
    std::string res = "(";
    res.append(std::to_string(dim[0]));
    for(int ii=1; ii<int(dim.size()); ++ii) res.append(std::string(",")+std::to_string(dim[ii]));

    res+= ")";
    return res;
  }

  // *************************************************** //

  template<typename T>
  struct getFitsType{
    static int run(){return 0;}
  };

  template<> struct getFitsType<double>{
    static int run(){return TDOUBLE; }
  };

  template<> struct getFitsType<float>{
    static int run(){return TFLOAT; }
  };
  
  template<> struct getFitsType<int>{
    static int run(){return TINT; }
  };

  template<> struct getFitsType<long>{
    static int run(){return TLONGLONG; }
  };
  
  template<> struct getFitsType<short>{
    static int run(){return TSHORT; }
  };
  
  template<> struct getFitsType<char>{
    static int run(){return TBYTE; }
  };
  // *************************************************** //

  template<typename T>
  void read_fits(std::string const& filename, size_t const nbuffer, T* buffer)
  {
    
    if(!file_exists(filename)){
      std::cerr<<"fits::read_fits: ERROR, file not found -> "<<filename<<std::endl;
      exit(1);
    }
    
    
    // --- Open file and get dimensions --- //
        
    fitsfile *afptr;
    long npixel = 1;
    int naxis = 0, status = 0, type = 0;
    long naxes[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1}, firstpix[9] = {1, 1, 1, 1, 1, 1, 1, 1, 1};
    
    fits_open_file(&afptr, filename.c_str(), READONLY, &status); 
    if(status){
      std::cerr<<"fits::read_fits: ERROR, file not open -> "<<filename<<std::endl;
      return;
    }

    fits_get_img_param(afptr, 9, &type, &naxis, naxes, &status);
   
    npixel = naxes[0] * naxes[1] * naxes[2] * naxes[3] * naxes[4] * naxes[5] * naxes[6] * naxes[7] * naxes[8];
    
    if(npixel > nbuffer){
      std::cerr<< "fits::read_fits: ERROR, buffer size != data size\n";
      exit(1);
    }
    
    int datatype = getFitsType<T>::run();
    
    //std::cout<<"fits::read_fits: info, reading file -> "<<filename<<std::endl;
    fits_read_pix(afptr, datatype, firstpix, npixel, NULL, buffer, NULL, &status);
    fits_close_file(afptr,  &status);
    
  } // Function
    // *************************************************** //

  std::vector<long> fits_getDims(std::string const& name)
  {
    
    std::vector<long> dim;
    
    
    if(!file_exists(name)){
      std::cerr<<"fits::read_fits: ERROR, file not found -> "<<name<<std::endl;
      exit(1);
    }
    
    
    /* --- open file --- */
    
    fitsfile *afptr;
    int status = 0, naxis = 0;
    
    fits_open_file(&afptr, name.c_str(), READONLY, &status);
    if(status){
      fprintf(stderr,"fits_getDims: Error, could not open file -> %s\n", name.c_str());
      return dim;
    }
    
    
    /* --- Get info --- */
    
    fits_get_img_dim(afptr, &naxis, &status);  /* read dimensions */
    dim.resize(naxis);
    fits_get_img_size(afptr, naxis, &dim[0], &status);
    
    fits_close_file(afptr,  &status);

    std::vector<long> out(naxis);
    for(int ii=0; ii<naxis;++ii) out[ii] = dim[naxis-1-ii];
    
    return out;
    
  }
  
  // *************************************************** //

  template<typename T>
  ml::Model<T> readModel(std::string const& filename){

    std::vector<long> dims = fits_getDims(filename);
    size_t const nElements = getTotalSize(dims);
    int const nDims = dims.size();

    if(nDims != 4){
      std::cerr<<"readModel: ERROR, input file must have 4 dimensions\n";
      exit(1);
    }

    int const nHydrogen = dims[2] - ml::NVAR;
    
    ml::Model<T> model(dims[0], dims[1], dims[3], nHydrogen);
    
    std::string sdims = dim2string(dims);
    std::cout<<"fits::read_model: reading file -> "<<filename<<" "<<sdims<<std::endl;
    
    read_fits(filename, nElements, &model.d(0,0,0,0));
    
    return model;
  }
  
  // *************************************************** //

  template<typename T, typename U>
  void write_fits(std::string const& filename, std::vector<long> const& dims,  T* data, bool const verbose)
  {
    bool over = false;
    if(file_exists(filename)){
      remove(filename.c_str());
      over = true;
    }

    if(verbose){
      if(over) fprintf(stdout,"write_fits: saving (re-writing) -> %s\n",
		       filename.c_str());
      else fprintf(stdout,"write_fits: saving (creating) -> %s\n",
		   filename.c_str());  
    }
    

    int status= 0, dtype = 0, naxes = dims.size();
    long int fpixel = 1, npix = 1;
    for(int ii=0; ii<naxes; ++ii) npix *= dims[ii];

    fitsfile *fptr = NULL;
    fits_create_file(&fptr, filename.c_str(), &status);

    int const naxes2 = naxes;
    long int idim[naxes2];

    for(int ii=0;ii<naxes2;ii++) idim[naxes-ii-1] = (long int)dims[ii];

  
  if     (typeid(U) == typeid(unsigned char)){
    dtype = TBYTE;
    fits_create_img(fptr, BYTE_IMG, naxes, &idim[0], &status);
  }else if(typeid(U) == typeid(short int)){
    dtype = TSHORT;
    fits_create_img(fptr, SHORT_IMG, naxes, &idim[0], &status);
  }else if(typeid(U) == typeid(int)){
    dtype = TINT;
    fits_create_img(fptr, LONG_IMG, naxes, &idim[0], &status);
  }else if(typeid(U) == typeid(long int)){
    dtype = TLONGLONG;
    fits_create_img(fptr, LONGLONG_IMG, naxes, &idim[0], &status);
  }else if(typeid(U) == typeid(float)){
    dtype = TFLOAT;
    fits_create_img(fptr, FLOAT_IMG, naxes, &idim[0], &status);
  }else if(typeid(U) == typeid(double)){
    dtype = TDOUBLE;
    fits_create_img(fptr, DOUBLE_IMG, naxes, &idim[0], &status);
  }else{
    fprintf(stderr,"fits_save_one: Error, unsupported type\n");
    return;
  }
    
  //dtype = getFitsType<T>::run();

    fits_write_img(fptr, dtype, fpixel, npix, data, &status); // Write the array of integers to the image
    fits_close_file(fptr, &status);
  }



  // *************************************************** //

  template<typename T> void writeModel(ml::Model<T> &m, std::string const& nameout)
  {   
    std::vector<long> dim = {m.ny, m.nx, ml::NVAR+m.nHydrogen, m.nDep};
    write_fits<T,float>(nameout, dim, &m.d(0,0,0,0), true);
  }

  // *************************************************** //

  

  
}



#endif
