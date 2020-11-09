#ifndef INTER_HPP
#define INTER_HPP

#include <algorithm>
#include <cmath>
#include <cstring>

//
// Coded by J. de la Cruz Rodriguez (ISP-SU, 2020)
//


namespace ipol{

  // ********************************************************************************* //

  template<class T> inline T signFortran(T const &val)
  {return ((val < static_cast<T>(0)) ? static_cast<T>(-1) : static_cast<T>(1));}
  

  // ********************************************************************************* //

  template<typename T> inline
  T getDerSteffen(T const &xu, T const &xc, T const &xd, T const &yu, T const &yc, T const &yd)
  {
    //
    // High order harmonic derivatives
    // Ref: Steffen (1990), A&A..239..443S
    //
    
    T const odx = xc - xu;
    T const dx  = xd - xc;
        
    T const S0 = (yd - yc) / dx;
    T const Su = (yc - yu) / odx;
    
    T const P0 = std::abs((Su*dx + S0*odx) / (odx+dx)) / 2;
    
    return (signFortran(S0) + signFortran(Su)) * std::min<T>(std::abs(Su),std::min<T>(std::abs(S0), P0));

  }
  
  // ********************************************************************************* //

  template<typename T>
  void getDerivatives(int const n, const T* const __restrict__ x, const T* const __restrict__ y, T* __restrict__ yp)
  {
    yp[0] = 0;
    yp[n-1] = 0;

    int const n1 = n-1;
    
    for(int ii=1; ii<n1; ++ii)
      yp[ii] = getDerSteffen(x[ii-1], x[ii], x[ii+1],y[ii-1], y[ii], y[ii+1]);
  }
  
  // ********************************************************************************* //

  template<typename T>
  void hermitian(int const n, const T* const __restrict__ x, const T* const __restrict__ y,
	       int const nn, const T* const __restrict__ xx,  T*  __restrict__ yy)
  {

    // --- Get derivatives --- //
    
    T* __restrict__ yp = new T [n]();
    getDerivatives<T>(n, x, y, yp);

    
    // --- limits and order --- //
    
    int k0=0, k1=0, dk=0;
    if((x[1]-x[0]) > 0){
      dk = 1, k0=-1, k1=n;
    }else{
      dk = -1, k0=n, k1=-1;
    }

    int j0=0, j1=0, dj=0;
    if((xx[1]-xx[0]) > 0){
      j0 = -1, j1=nn; dj=1;
    }else{
      j1 = -1, j0=nn; dj=-1;
    }
			 
    
    int const kk0 = k0;
    int const kk1 = k1;
    int const dkk = dk;
    
    int const jj1 = j1;
    int const jj0 = j0;
    int const djj = dj;
    
    int off = jj0+djj;

    for(int kk=kk0+2*dkk; kk != kk1; kk+=dkk){
      T const dx = x[kk] - x[kk-dkk];

      int const off1 = off;
      
      for(int jj=off1; jj != jj1; jj += djj){
	
	// --- clip value to the existing data range, it will repeat edge values automatically --- //
	T const ixx = std::max<T>(std::min<T>(xx[jj], x[kk1-dkk]), x[kk0+dkk]);
	  
	if((ixx >= x[kk-dkk]) && (ixx <= x[kk])){
	  T const u = (ixx-x[kk-dkk])/dx;
	  T const uu = u*u;
	  T const uuu = u*uu;
	  
	  yy[jj] = y[kk-dkk] * (1 - 3*uu + 2*uuu) + (3*uu - 2*uuu) * y[kk]
	    + (uuu - 2*uu + u) * dx * yp[kk-dkk] + (uuu - uu) * dx * yp[kk];
	    
	  off += djj;
	} // if
      } // jj
    } // kk

    
    delete [] yp;
  }
  
  // ********************************************************************************* //

  template<typename T>
  void linear(int const n, const T* const __restrict__ x, const T* const __restrict__ y,
	      int const nn, const T* const __restrict__ xx,  T*  __restrict__ yy)
  {

    // --- limits and order --- //
    
    int k0=0, k1=0, dk=0;
    if((x[1]-x[0]) > 0){
      dk = 1, k0=-1, k1=n;
    }else{
      dk = -1, k0=n, k1=-1;
    }

    int j0=0, j1=0, dj=0;
    if((xx[1]-xx[0]) > 0){
      j0 = -1, j1=nn; dj=1;
    }else{
      j1 = -1, j0=nn; dj=-1;
    }
			 
    
    int const kk0 = k0;
    int const kk1 = k1;
    int const dkk = dk;
    
    int const jj1 = j1;
    int const jj0 = j0;
    int const djj = dj;
    
    int off = jj0+djj;

    for(int kk=kk0+2*dkk; kk != kk1; kk+=dkk){
      T const dx = x[kk] - x[kk-dkk];
      
      T const a = (y[kk] - y[kk-dkk]) / dx;
      T const b = y[kk-dkk] - a*x[kk-dkk];

      int const off1 = off;
      
      for(int jj=off1; jj != jj1; jj += djj){
	
	// --- clip value to the existing data range, it will repeat edge values automatically --- //
	T const ixx = std::max<T>(std::min<T>(xx[jj], x[kk1-dkk]), x[kk0+dkk]);

	
	if((ixx >= x[kk-dkk]) && (ixx <= x[kk])){	  
	  yy[jj] = a*ixx + b;
	    
	  off += djj;
	} // if
      } // jj
    } // kk
  }

  // ********************************************************************************* //

  template<typename T> void bezier_integral(int const N, const T* const __restrict__ x,
				       const T* const __restrict__ y, T* __restrict__ res)
  {

    // --- Compute derivative of y --- //
    
    T* __restrict__ yp = new T [N]();
    getDerivatives<T>(N, x, y, yp);
    
    
    // --- integrate --- //

    res[0] = 0;
    for(int ii=1; ii<N; ++ii){
      double const dx = x[ii]-x[ii-1];

      double const Cu = y[ii-1] + (dx*yp[ii-1])/3.0;
      double const C0 = y[ii]   - (dx*yp[ii]  )/3.0;

      // --- Bezier integral ---//
      
      res[ii] =  (std::abs(dx) * (y[ii-1] + y[ii] + C0 + Cu) / 4);
    }


    // --- extrapolate tau_500 at the top boundary --- //
    
    T const y1 = res[1];
    T const y2 = y1 + res[2];
    res[0] = exp(2*log(y1) - log(y2));

    double iTau = 0;
    
    for(int ii=0; ii<N; ++ii){
      iTau += res[ii];
      res[ii] =  log10(iTau);
    }
    
    
    delete [] yp;
  }


  // ********************************************************************************* //

  template<typename T> void smooth_and_scale_gradients(int const N, T* __restrict__ d, int const wsize, T const scaling)
  {
    int const N1 = N-1;

    
    // --- standard smoothing with top-hat PSF --- //
    
    if(wsize > 1){

      int const w2 = wsize / 2;
      int const Nw2 = N+w2;
      
      T* __restrict__ d_orig = new T [N+2*w2];
      std::memcpy(&d_orig[w2], d, sizeof(T)*N);
      for(int ii=0; ii<w2; ++ii){
	d_orig[ii]      = d[0];
	d_orig[ii+Nw2]  = d[N1]; 
      }
      
      
      for(int kk = w2; kk < Nw2; ++kk){
	int const j0 = kk-w2;
	int const j1 = kk+w2;
	
	T sum = 0;
	for(int jj=j0; jj<=j1; ++jj)
	  sum += d_orig[jj];
	
	d[kk-w2] = sum / (j1-j0+1); 
      }
      
      delete [] d_orig;   
    }

    
    // --- make sure that the scaling of the array is correct --- //
    
    T const off = d[0];
    T const range = scaling / (d[N1] - d[0]);
    
    for(int ii=0; ii<N; ++ii){
      d[ii] = (d[ii]-off) * range;
    }
    
    
  }
  
  // ********************************************************************************* //

  template<typename T> inline T* arange(int const N)
  {
    T* __restrict__ res = new T [N]();
    for(int ii=0; ii<N; ++ii) res[ii] = T(ii);
    return res;
  }
  
  // ********************************************************************************* //


}


#endif

