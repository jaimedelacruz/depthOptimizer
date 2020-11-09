// **************************************************************************** //

/*
  Depth optimization routines
  Gradient definitions taken from MULTI (Carlsson 1986)
  but with the adition of a search in velocity too.

  Coded by J. de la Cruz Rodriguez (ISP-SU 2020)
  
 */

// **************************************************************************** //

#include <vector>
#include <cstdio>
#include <iostream>
#include <string>
#include <omp.h>

#include "interpolation.hpp"
#include "model.hpp"
#include "readModelFits.hpp"

using namespace std;

// **************************************************************************** //

int main(int const narg, char *argv[])
{
  using fp = float;
  
  // --- read input parameters --- //
  
  if(narg < 8){
    cerr << "USAGE: ./optimizer.x filein.h5 fileout.h5 nthreads convert_units wsize tempmax taumax" << endl;
    exit(0);
  }

  string const filein   = string(argv[1]);
  string const fileout  = string(argv[2]);
  int    const nthreads = std::max<int>(1, stoi(argv[3]));
  int    const units    = std::max<int>(0, stoi(argv[4]));
  int    const wsize    = std::max<int>(1, stoi(argv[5]));
  fp     const tempmax  = std::max<int>(1, stod(argv[6]));
  fp     const taumax   = std::max<int>(1, stod(argv[7]));

  cout << "main: in="<<filein<<", out="<<fileout<<", nthreads="<<nthreads<<endl;
  cout << "main: convert_to_SI="<<units<<", smooth_window="<<wsize<<", Tg_max="<<tempmax<<", ltau_max="<<taumax<<endl<<endl;


  
  // --- Read input model --- //

  int nt = 0;
  ml::Model<fp> model =  fits::readModel<fp>(filein);



  

  // --- do I need to fill hydrogen? --- //

  bool const fill_hydrogen = (((model.d(0,0,ml::NVAR,0) + model.d(0,0,ml::NVAR,1)) < 1.e-32) && (model.nHydrogen <= 6)) ? true : false;


  
  
  // --- Prepare variables for parallel loop --- //
  
  vector<line> dum;
  eos::witt *eos = NULL;  
 
  long ipix = 0 , tid = 0, xx=0, yy=0;
  
  int const ny = model.ny;
  int const nx = model.nx;
  long const npix = nx*ny;
  int oper = -1, per = 0;


  
  // --- parallel loop --- //
#pragma omp parallel default(shared) firstprivate(ipix, tid, xx, yy, eos) num_threads(nthreads)
  {
    tid = omp_get_thread_num();
    eos =  new eos::witt(dum);

#pragma omp for schedule(dynamic, 2)
    for(ipix=0; ipix<npix; ipix++){


      // --- get pixel coordinates --- //
      yy = ipix / npix;
      xx = ipix - yy*npix;


      // --- extract pixel and optimize --- //
      model.template operator()<double>(yy,xx).Optimize(fill_hydrogen, units, *eos, tempmax, taumax, wsize);


      // --- progress count --- //
      if(tid == 0){
	per = (ipix * 100) / std::max<int>(1,npix-1);
	if(per != oper){
	  oper = per;
         fprintf(stderr, "\rprocessing -> %3d%s", per,"%");
	}
      }
    }

    delete eos;
  }// parallel block
  fprintf(stderr, "\rprocessing -> %3d%s\n", 100,"%");


  // --- write to disk --- //

  fits::writeModel(model, fileout);
  
  
}

// **************************************************************************** //
