#ifndef MODELHPP
#define MODELHPP

/* ---
   Model tools for depth optimization

   Coded by J. de la Cruz Rodriguez (ISP-SU 2020)
   --- */

#include <cstring>
#include <typeinfo>

#include "Arrays.hpp"
#include "partition.h"
#include "witt.h"

namespace ml{
  
  // ------------------------------------------------------------------------ //

  // --- z, temp, vlos, vturb, ne, dens, bx, by, bz, nH --- //
  constexpr static const int NVAR = 10;
  
  // ------------------------------------------------------------------------ //

  // ******************************************************************************* //
  
  template<typename Z, typename T>
  class model1D{
  public:
    int const nDep;
    int const nHydrogen;
    T *z, *tau, *temp, *vz, *vturb, *Bx, *By, *Bz, *Ne, *rho, *nH; 
    Z *data_in;
    
    model1D(int const ndep, int const nhyd, Z* data): nDep(ndep), nHydrogen(nhyd), data_in(data)
    {

      // --- We might want to have the input data as float but perform
      // --- all the operations in double.
      // --- Make a copy on the fly and copy back in the destructor

      if(typeid(Z) == typeid(T)){
	z = (T*)data_in; 
      }else{
	int const ntot = nDep*(ml::NVAR+nHydrogen);
	z = new T [ntot]();
	for(int ii=0; ii<ntot; ++ii) z[ii] = T(data_in[ii]);
      }
      
      tau   = z + ndep*1;
      temp  = z + ndep*2;
      vz    = z + ndep*3;
      vturb = z + ndep*4;
      Bx    = z + ndep*5;
      By    = z + ndep*6;
      Bz    = z + ndep*7;
      Ne    = z + ndep*8;
      rho   = z + ndep*9;
      nH    = z + ndep*10;
    }

    // ------------------------------------------------------------------------ //

    ~model1D(){
      if(typeid(T) != typeid(Z)){
	
	int const ntot = nDep*(ml::NVAR+nHydrogen);
	for(int ii=0; ii<ntot; ++ii) data_in[ii] = Z(z[ii]);
	delete [] z;
	
      }
    }
    
    // ------------------------------------------------------------------------ //

    void to_CGS()const{
      
      constexpr static const T M_TO_CM        = 100;
      constexpr static const T TESLA_TO_GAUSS = 10000;
      constexpr static const T M3_TO_CM3      = M_TO_CM * M_TO_CM * M_TO_CM;
      constexpr static const T KG_TO_GR       = 1000;
      constexpr static const T RHO_TR         = KG_TO_GR / M3_TO_CM3;  

      
      for(int ii=0; ii<nDep; ++ii) z[ii] *= M_TO_CM;	
      for(int ii=0; ii<nDep; ++ii) vz[ii]*= M_TO_CM;
      for(int ii=0; ii<nDep; ++ii) vturb[ii]*= M_TO_CM;
      for(int ii=0; ii<nDep; ++ii) Bx[ii]*= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) By[ii]*= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) Bz[ii]*= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) Ne[ii] /= M3_TO_CM3;
      for(int ii=0; ii<nDep; ++ii) rho[ii]*= RHO_TR;

      for(int jj=0; jj<nHydrogen; ++jj){
	int const off = jj*nDep;
	for(int ii=0; ii<nDep; ++ii) nH[off+ii]/=  M3_TO_CM3;
      }
    }
    
    
    // ------------------------------------------------------------------------ //

    void to_SI()const{
      
      constexpr static const T M_TO_CM        = 100;
      constexpr static const T TESLA_TO_GAUSS = 10000;
      constexpr static const T M3_TO_CM3      = M_TO_CM * M_TO_CM * M_TO_CM;
      constexpr static const T KG_TO_GR       = 1000;
      constexpr static const T RHO_TR         = KG_TO_GR / M3_TO_CM3;  

      
      for(int ii=0; ii<nDep; ++ii) z[ii] /= M_TO_CM;	
      for(int ii=0; ii<nDep; ++ii) vz[ii]/= M_TO_CM;
      for(int ii=0; ii<nDep; ++ii) vturb[ii]/= M_TO_CM;
      for(int ii=0; ii<nDep; ++ii) Bx[ii]/= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) By[ii]/= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) Bz[ii]/= TESLA_TO_GAUSS;
      for(int ii=0; ii<nDep; ++ii) Ne[ii] *= M3_TO_CM3;
      for(int ii=0; ii<nDep; ++ii) rho[ii]/= RHO_TR;

      for(int jj=0; jj<nHydrogen; ++jj){
	int const off = jj*nDep;
	for(int ii=0; ii<nDep; ++ii) nH[off+ii]*=  M3_TO_CM3;
      }
    }
    
    // ------------------------------------------------------------------------ //

    template<typename U>
    static void _fill_nH_LTE_6(int const nLev, U const& temp,  U const& Pg,
			       U const& Pe, T (&H)[6] , eos::witt &eos)
    {
      constexpr static const U BK = 1.3806488E-16;
      
      
      if(nLev > 6){
	std::cerr<<"model1D::_fill_nH_LTE_6: your H populations must have at most 6 levels in order for this routine to fill the pupulations, exiting! nLev = "<<nLev <<std::endl;
	  exit(1);
      }
      
      // --- Energy [ergs]: {0.000, 82258.211, 97491.219, 102822.766, 105290.508, 109677.617} * H * C
      constexpr static const U ELEV[6] = {0.0, 1.63401468e-11, 1.93661011e-11, 2.04251840e-11,
					  2.09153875e-11, 2.17868629e-11};
      
      constexpr static const U G[6]    = {2,8,18,32,50,1};
      
      U const BKT = BK * temp;
      U pp[6] = {};//, ein[6] = {}, pf[6] = {};

      
      // --- solve EOS for H --- //

      U iPe = Pe, iPg = Pg,  itemp = temp;
      eos.gasc<U>(itemp, iPe, iPg, pp); 


      // --- copy total H I and H II --- //

      U const nH1 = (pp[0] * pp[4]) / BKT;
      U const nH2 = (pp[1] * pp[4]) / BKT;

      
      
      // --- Solve Boltzmann eq. for bound levels --- //
      
      double nH1tot = 1.0;
      int const nLev1 = nLev-1;
      H[0]     = 1.0;
      H[nLev1] = nH2;

      for(int ii = 1; ii<nLev1; ++ii){
	H[ii] = (G[ii] / G[0]) * exp(-ELEV[ii]/BKT); // nH_i / nH_0
	nH1tot += H[ii]; // nHI / nHI_0
      }

      
      // --- Calculate (nHI_i / nHI_0) x nHI_0 to get the actual level populations --- //
      
      nH1tot = nH1 / nH1tot;
      for(int ii=0; ii<nLev1; ++ii) H[ii] *=nH1tot;
    }

    // ------------------------------------------------------------------------ //

    static void _get_tau_scale(int const ndep, int const nLev, eos::witt &eos, const T* const __restrict__ z,
			       const T* const __restrict__ temp,
			       T* __restrict__ rho, T* __restrict__ Ne, T* __restrict__ ltau,
			       T* __restrict__ nH, bool const solve_nh)
    {
      constexpr static const double BK = 1.3806488E-16;
      double Wav = 5000.;

      T* __restrict__ alpha_500 = new T [ndep]();

      // --- Loop through all depth and compute the relevant quantities --- //

      for(int kk=0; kk<ndep; ++kk){
	double Pe  = Ne[kk]*BK*temp[kk];
	double Pg = eos.pg_from_rho<double>(double(temp[kk]), double(rho[kk]), Pe);

	// --- Solve nH ? --- //
	
	if(solve_nh){
	  double H[6] = {};
	  _fill_nH_LTE_6<double>(nLev, double(temp[kk]), Pg, Pe, H, eos);
	  for(int ii=0; ii<nLev; ++ii) nH[ii*ndep + kk] = H[ii];
	}


	// --- Fill background opacity at 500 nmm --- //

	double opac = 0;
	eos.contOpacity<double>(double(temp[kk]), Pg, Pe, 1, &Wav, &opac);
	alpha_500[kk] = opac;
      }


      
      // --- Compute Tau-scale --- //
      
      ipol::bezier_integral(ndep, z, alpha_500, ltau);
      delete [] alpha_500;
    }

    
    // ------------------------------------------------------------------------ //

    
    void Optimize(bool const fill_Hydrogen, int const convert_units, eos::witt &eos,
		  T const temp_max, T const ltau_max, int const wsmooth, int const nDep2)const
    {
      static T const log11           = 1 / log10(1.1);
      constexpr static const T vscal = 1.0E-5 / 0.5;

      
      int const ndep = nDep;
      
      // --- convert to CGS? --- //
      
      if(convert_units)
	to_CGS();



      
      // --- calculate tau scale and fill nH if neccesary --- //

      _get_tau_scale(ndep, nHydrogen, eos, z, temp, rho, Ne, tau, nH, bool(fill_Hydrogen));

      
      // --- Find limits --- //

      int k0 = 0, k1 = nDep-1;
      for(int kk=1; kk<nDep; ++kk){
	if((temp[kk] > temp_max) && (k0 == kk-1)) k0 = kk;
	if( tau[kk] <= ltau_max)                  k1 = kk;
      }
      
      k1 = std::min<int>(k1+1, nDep-1);
      if((k1-k0) < nDep/3) k1 = nDep-1;
      
      int const kk0 = k0+1;
      int const kk1 = k1; 

      
      // --- compute gradients --- //
      
      T* __restrict__ aind = new T [nDep]();

      for(int kk=kk0; kk<= kk1; ++kk){
	T const tdiv = std::abs(log10(temp[kk]) - log10(temp[kk-1])) * log11;
	T const rdiv = std::abs(log10(rho[kk])  - log10( rho[kk-1])) * log11;
	T const vdiv = std::abs(vz[kk]  -  vz[kk-1]) * vscal;
	T const ldiv = std::abs(tau[kk] - tau[kk-1]) * 10;

	// --- take the largest --- //
	
	aind[kk] = aind[kk-1] + std::max<T>(std::max<T>(std::max<T>(tdiv, rdiv), vdiv), ldiv);
      }


      
      // --- smooth gradients --- //
      
      int const nAind = k1 - k0 + 1;
      ipol::smooth_and_scale_gradients<T>(nAind, &aind[k0], wsmooth, nDep2-1);

      T* __restrict__ index = ipol::arange<T>(nDep2);
      T* __restrict__ buffer = new T [nDep2]();


      
      // --- interpolate all variables --- //

      int const nTot  = NVAR + nHydrogen;
      
      for(int ii = 0; ii<nTot; ++ii){
	T* __restrict__ iVar = z + ii*nDep; 
	ipol::linear<T>(nAind, &aind[k0], &iVar[k0], nDep2, index, buffer);
	std::memcpy(iVar, buffer, nDep2*sizeof(T));
      }
      
      
      // --- clean up --- //
      
      delete [] buffer;
      delete [] index;
      delete [] aind;


      if(convert_units)
	to_SI();

    }
    
    // ------------------------------------------------------------------------ //

  };
  
  // ******************************************************************************* //
  
  template<typename T> class Model{
  public:
    int ny, nx, nDep, nHydrogen;
    mem::Array<T,4> d;
   
    ~Model(){};
    
    Model(): ny(0), nx(0), nDep(0), nHydrogen(0), d(){};
    
    // ------------------------------------------------------------------------ //

    Model(int const ny_in, int const nx_in, int const ndep_in, int const nHyd):
      ny(ny_in), nx(nx_in), nDep(ndep_in), nHydrogen(nHyd),d(ny_in, nx_in, NVAR+nHyd, ndep_in){d.Zero();}
    
    // ------------------------------------------------------------------------ //

    Model(Model<T> const &in):
      ny(in.ny), nx(in.nx), nDep(in.nDep), nHydrogen(in.nHydrogen), d(in.d){};
    
    // ------------------------------------------------------------------------ //
    template<typename U>
    model1D<T,U> operator()(long const yy, long const xx)
    {
      return model1D<T,U>(nDep, nHydrogen, &d(yy,xx,0,0));
    }

    // ------------------------------------------------------------------------ //
    
  };
  
  // ******************************************************************************* //

  


  
}
#endif
