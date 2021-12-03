#ifndef TIME_EVOLUTION_H
#define TIME_EVOLUTION_H


#include <iostream>
#include <iterator>
#include <algorithm>
#include "boost/multi_array.hpp"
#include <cassert>
#include <complex>
#include <mpi.h>
#include <cmath>

using namespace std;


class time_evolution
{
private:
  
public:
  
  
  time_evolution()
  {
  }
  virtual ~time_evolution() { }
  
  
   int64_t site(int64_t n)       //implement periodic bounndary condition
  {
    int64_t nn;
    if(n>-1 && n<N1) {nn=n;}
    else if (n>N1-1){nn=n-N1;}
    else if (n<0) {nn=N1+n;}
    return (nn);
  } 
  
  void operator() (complex<double> *&dxdtF11,complex<double> *&dxdtF12,complex<double> *&dxdtF13,complex<double> *&dxdtF21,complex<double> *&dxdtF22,complex<double> *&dxdtF23,complex<double> *&dxdtF31,complex<double> *&dxdtF32,complex<double> *&dxdtF33, complex<double> *&dxdtrho11,complex<double> *&dxdtrho12,complex<double> *&dxdtrho13,complex<double> *&dxdtrho21,complex<double> *&dxdtrho22,complex<double> *&dxdtrho23,complex<double> *&dxdtrho31,complex<double> *&dxdtrho32,complex<double> *&dxdtrho33, FandIni &F, GetMSigma &getmsigma, MemoryInt &memint, const int64_t t, const int64_t tt, int Direction)
  {
  if(tt>1)
  {
   if(Direction==0)//reverse direction
   {
   #pragma omp parallel for 
   for(int64_t z1=0; z1<N1; z1++)
   {
     for(int64_t z2=0; z2<N1; z2++)
     {
       for(int64_t qperp=0; qperp<Nt; qperp++)
       {
     dxdtrho11[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot11_access(z1, site(z2+1), qperp, tt)+F.rhot11_access(z1, site(z2-1), qperp, tt)-2.*F.rhot11_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rhot13_access(z1, site(z2+1), qperp, tt)-F.rhot13_access(z1, z2, qperp, tt));
      dxdtrho21[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot21_access(z1, site(z2+1), qperp, tt)+F.rhot21_access(z1, site(z2-1), qperp, tt)-2.*F.rhot21_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rhot23_access(z1, site(z2+1), qperp, tt)-F.rhot23_access(z1, z2, qperp, tt)); 
      dxdtrho31[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot31_access(z1, site(z2+1), qperp, tt)+F.rhot31_access(z1, site(z2-1), qperp, tt)-2.*F.rhot31_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rhot33_access(z1, site(z2+1), qperp, tt)-F.rhot33_access(z1, z2, qperp, tt)); 
      dxdtrho12[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot12_access(z1, site(z2+1), qperp, tt)+F.rhot12_access(z1, site(z2-1), qperp, tt)-2.*F.rhot12_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot12_access(z1, z2, qperp, tt);
      dxdtrho22[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot22_access(z1, site(z2+1), qperp, tt)+F.rhot22_access(z1, site(z2-1), qperp, tt)-2.*F.rhot22_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot22_access(z1, z2, qperp, tt);
      dxdtrho32[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rhot32_access(z1, site(z2+1), qperp, tt)+F.rhot32_access(z1, site(z2-1), qperp, tt)-2.*F.rhot32_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot32_access(z1, z2, qperp, tt);
 
     complex<double> bufferrho13;
     complex<double> bufferrho23;
     complex<double> bufferrho33;
   
     bufferrho13=0.;
     bufferrho23=0.;
     bufferrho33=0.;
     #if (ORDER==0)
     bufferrho13=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.rhot13_access(z1, z2, qperp, tt);
     bufferrho23=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.rhot23_access(z1, z2, qperp, tt);
     bufferrho33=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.rhot33_access(z1, z2, qperp, tt);
     #elif (ORDER==1)
     int i=0;
     for(int64_t k=1;k<3;k++)
     {
       for(int64_t n=0; n<k;n++)
       {
         for(int64_t m=0; m<k; m++)
         {
     bufferrho13+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot13_access(z1, site(z2+n-m), qperp, tt);
     bufferrho23+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot23_access(z1, site(z2+n-m), qperp, tt);
     bufferrho33+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot33_access(z1, site(z2+n-m), qperp, tt);
     i++;
         }
       }
     }
     #elif (ORDER==2)
     int i=0;
     for(int64_t k=1;k<4;k++)
     {
       for(int64_t n=0; n<k;n++)
       {
         for(int64_t m=0; m<k; m++)
         {
     bufferrho13+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot13_access(z1, site(z2+n-m), qperp, tt);
     bufferrho23+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot23_access(z1, site(z2+n-m), qperp, tt);
     bufferrho33+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.rhot33_access(z1, site(z2+n-m), qperp, tt);
     i++;
         }
       }
     }
     #endif     
     dxdtrho13[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot13_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rhot11_access(z1, z2, qperp, tt)-F.rhot11_access(z1, site(z2-1), qperp, tt))+bufferrho13+memint.SigmarhoP31_access(z1, z2, qperp);
     dxdtrho23[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot23_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rhot21_access(z1, z2, qperp, tt)-F.rhot21_access(z1, site(z2-1), qperp, tt))+bufferrho23+memint.SigmarhoP32_access(z1, z2, qperp);
     dxdtrho33[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rhot33_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rhot31_access(z1, z2, qperp, tt)-F.rhot31_access(z1, site(z2-1), qperp, tt))+bufferrho33+memint.SigmarhoP33_access(z1, z2, qperp);
  //if((z1==10)&&(z2==19)&&(qperp==1)){cout<<"M:"<<buffer31<<endl; cout<<"sigmarhoF31:"<<memint.SigmarhoF31_access(z1, z2, qperp)<<endl; cout<<"sigmaFP31:"<<memint.SigmaFP31_access(z1, z2, qperp)<<endl;}
     }
     
    }
   }
   }
   else//positive direction
   {
   #pragma omp parallel for 
   for(int64_t z1=0; z1<N1; z1++)
   {
     for(int64_t z2=0; z2<N1; z2++)
     {
       for(int64_t qperp=0; qperp<Nt; qperp++)
       {
      
  
     dxdtF11[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft11_access(z1, site(z2+1), qperp, tt)+F.Ft11_access(z1, site(z2-1), qperp, tt)-2.*F.Ft11_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.Ft13_access(z1, site(z2+1), qperp, tt)-F.Ft13_access(z1, z2, qperp, tt));
      dxdtF21[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft21_access(z1, site(z2+1), qperp, tt)+F.Ft21_access(z1, site(z2-1), qperp, tt)-2.*F.Ft21_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.Ft23_access(z1, site(z2+1), qperp, tt)-F.Ft23_access(z1, z2, qperp, tt)); 
      dxdtF31[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft31_access(z1, site(z2+1), qperp, tt)+F.Ft31_access(z1, site(z2-1), qperp, tt)-2.*F.Ft31_access(z1, z2, qperp, tt))-1./a_t*(1.-exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.Ft33_access(z1, site(z2+1), qperp, tt)-F.Ft33_access(z1, z2, qperp, tt)); 
      dxdtF12[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft12_access(z1, site(z2+1), qperp, tt)+F.Ft12_access(z1, site(z2-1), qperp, tt)-2.*F.Ft12_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft12_access(z1, z2, qperp, tt);
      dxdtF22[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft22_access(z1, site(z2+1), qperp, tt)+F.Ft22_access(z1, site(z2-1), qperp, tt)-2.*F.Ft22_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft22_access(z1, z2, qperp, tt);
      dxdtF32[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.Ft32_access(z1, site(z2+1), qperp, tt)+F.Ft32_access(z1, site(z2-1), qperp, tt)-2.*F.Ft32_access(z1, z2, qperp, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft32_access(z1, z2, qperp, tt);

     complex<double> buffer13;
     complex<double> buffer23;
     complex<double> buffer33;
     
     buffer13=0.;
     buffer23=0.;
     buffer33=0.;
  
     #if (ORDER==0)
     buffer13=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.Ft13_access(z1, z2, qperp, tt);
     buffer23=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.Ft23_access(z1, z2, qperp, tt);
     buffer33=1./(a_t*a_t)*getmsigma.M_ptr[z2]*F.Ft33_access(z1, z2, qperp, tt);
  
     #elif (ORDER==1)
     int i=0;
     for(int64_t k=1;k<3;k++)
     {
       for(int64_t n=0; n<k;n++)
       {
         for(int64_t m=0; m<k; m++)
         {
     buffer13+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft13_access(z1, site(z2+n-m), qperp, tt);
     buffer23+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft23_access(z1, site(z2+n-m), qperp, tt);
     buffer33+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft33_access(z1, site(z2+n-m), qperp, tt);   
     i++;
         }
       }
     }
     #elif (ORDER==2)
     int i=0;
     for(int64_t k=1;k<4;k++)
     {
       for(int64_t n=0; n<k;n++)
       {
         for(int64_t m=0; m<k; m++)
         {
     buffer13+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft13_access(z1, site(z2+n-m), qperp, tt);
     buffer23+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft23_access(z1, site(z2+n-m), qperp, tt);
     buffer33+=1./(a_t*a_t)*getmsigma.M_ptr[z2+N1*i]*F.Ft33_access(z1, site(z2+n-m), qperp, tt);   
     i++;
         }
       }
     }
     #endif     
  

     dxdtF13[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft13_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.Ft11_access(z1, z2, qperp, tt)-F.Ft11_access(z1, site(z2-1), qperp, tt))+buffer13+memint.SigmarhoF31_access(z1, z2, qperp)+memint.SigmaFP31_access(z1, z2, qperp);
     dxdtF23[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft23_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.Ft21_access(z1, z2, qperp, tt)-F.Ft21_access(z1, site(z2-1), qperp, tt))+buffer23+memint.SigmarhoF32_access(z1, z2, qperp)+memint.SigmaFP32_access(z1, z2, qperp);
     dxdtF33[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.Ft33_access(z1, z2, qperp, tt)-1./a_t*(exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.Ft31_access(z1, z2, qperp, tt)-F.Ft31_access(z1, site(z2-1), qperp, tt))+buffer33+memint.SigmarhoF33_access(z1, z2, qperp)+memint.SigmaFP33_access(z1, z2, qperp);
  //if((z1==10)&&(z2==19)&&(qperp==1)){cout<<"M:"<<buffer31<<endl; cout<<"sigmarhoF31:"<<memint.SigmarhoF31_access(z1, z2, qperp)<<endl; cout<<"sigmaFP31:"<<memint.SigmaFP31_access(z1, z2, qperp)<<endl;}
     }
     
    }
   }
   }
  
  }
  else
  {
   #pragma omp parallel for
   for(int64_t z1=0; z1<N1; z1++)
   {
     for(int64_t z2=0; z2<N1; z2++)
     {
       for(int64_t qperp=0; qperp<Nt; qperp++)
       {
      
     dxdtF11[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F11_access(site(z1+1), z2, qperp, t, tt)+F.F11_access(site(z1-1), z2, qperp, t, tt)-2.*F.F11_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.F31_access(site(z1+1), z2, qperp, t, tt)-F.F31_access(z1, z2, qperp, t, tt));
      dxdtF12[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F12_access(site(z1+1), z2, qperp, t, tt)+F.F12_access(site(z1-1), z2, qperp, t, tt)-2.*F.F12_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.F32_access(site(z1+1), z2, qperp, t, tt)-F.F32_access(z1, z2, qperp, t, tt)); 
      dxdtF13[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F13_access(site(z1+1), z2, qperp, t, tt)+F.F13_access(site(z1-1), z2, qperp, t, tt)-2.*F.F13_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.F33_access(site(z1+1), z2, qperp, t, tt)-F.F33_access(z1, z2, qperp, t, tt)); 
      dxdtF21[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F21_access(site(z1+1), z2, qperp, t, tt)+F.F21_access(site(z1-1), z2, qperp, t, tt)-2.*F.F21_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F21_access(z1, z2, qperp, t, tt);
      dxdtF22[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F22_access(site(z1+1), z2, qperp, t, tt)+F.F22_access(site(z1-1), z2, qperp, t, tt)-2.*F.F22_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F22_access(z1, z2, qperp, t, tt);
      dxdtF23[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.F23_access(site(z1+1), z2, qperp, t, tt)+F.F23_access(site(z1-1), z2, qperp, t, tt)-2.*F.F23_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F23_access(z1, z2, qperp, t, tt);
 
     dxdtrho11[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho11_access(site(z1+1), z2, qperp, t, tt)+F.rho11_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho11_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rho31_access(site(z1+1), z2, qperp, t, tt)-F.rho31_access(z1, z2, qperp, t, tt));
      dxdtrho12[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho12_access(site(z1+1), z2, qperp, t, tt)+F.rho12_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho12_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rho32_access(site(z1+1), z2, qperp, t, tt)-F.rho32_access(z1, z2, qperp, t, tt)); 
      dxdtrho13[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho13_access(site(z1+1), z2, qperp, t, tt)+F.rho13_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho13_access(z1, z2, qperp, t, tt))-1./a_t*(1.-exp(-M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*1./a*(F.rho33_access(site(z1+1), z2, qperp, t, tt)-F.rho33_access(z1, z2, qperp, t, tt)); 
      dxdtrho21[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho21_access(site(z1+1), z2, qperp, t, tt)+F.rho21_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho21_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho21_access(z1, z2, qperp, t, tt);
      dxdtrho22[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho22_access(site(z1+1), z2, qperp, t, tt)+F.rho22_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho22_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho22_access(z1, z2, qperp, t, tt);
      dxdtrho23[z1+ z2*N1+qperp*N1*N1]=1./(a*a)*(F.rho23_access(site(z1+1), z2, qperp, t, tt)+F.rho23_access(site(z1-1), z2, qperp, t, tt)-2.*F.rho23_access(z1, z2, qperp, t, tt))-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho23_access(z1, z2, qperp, t, tt);
 
     complex<double> buffer31;
     complex<double> buffer32;
     complex<double> buffer33;
     complex<double> bufferrho31;
     complex<double> bufferrho32;
     complex<double> bufferrho33;
     buffer31=0.;
     buffer32=0.;
     buffer33=0.;
     bufferrho31=0.;
     bufferrho32=0.;
     bufferrho33=0.;
     #if (ORDER==0)
     buffer31=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F31_access(z1, z2, qperp, t, tt);
     buffer32=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F32_access(z1, z2, qperp, t, tt);
     buffer33=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F33_access(z1, z2, qperp, t, tt);
     bufferrho31=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho31_access(z1, z2, qperp, t, tt);
     bufferrho32=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho32_access(z1, z2, qperp, t, tt);
     bufferrho33=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho33_access(z1, z2, qperp, t, tt);
     #elif (ORDER==1)
     buffer31=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F31_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.F31_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.F31_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.F31_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.F31_access(z1, z2, qperp, t, tt);
     buffer32=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F32_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.F32_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.F32_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.F32_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.F32_access(z1, z2, qperp, t, tt);
     buffer33=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.F33_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.F33_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.F33_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.F33_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.F33_access(z1, z2, qperp, t, tt);
 
     bufferrho31=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho31_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.rho31_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.rho31_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.rho31_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.rho31_access(z1, z2, qperp, t, tt);
     bufferrho32=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho32_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.rho32_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.rho32_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.rho32_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.rho32_access(z1, z2, qperp, t, tt);
     bufferrho33=1./(a_t*a_t)*getmsigma.M_ptr[z1]*F.rho33_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+N1]*F.rho33_access(z1, z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+2*N1]*F.rho33_access(site(z1-1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+3*N1]*F.rho33_access(site(z1+1), z2, qperp, t, tt)+1./(a_t*a_t)*getmsigma.M_ptr[z1+4*N1]*F.rho33_access(z1, z2, qperp, t, tt);
 
     #elif (ORDER==2)
     int i=0;
     for(int64_t k=1;k<4;k++)
     {
       for(int64_t n=0; n<k;n++)
       {
         for(int64_t m=0; m<k; m++)
         {
     buffer31+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.F31_access(site(z1+n-m), z2, qperp, t, tt);
     buffer32+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.F32_access(site(z1+n-m), z2, qperp, t, tt);
     buffer33+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.F33_access(site(z1+n-m), z2, qperp, t, tt);      
     bufferrho31+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.rho31_access(site(z1+n-m), z2, qperp, t, tt);
     bufferrho32+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.rho32_access(site(z1+n-m), z2, qperp, t, tt);
     bufferrho33+=1./(a_t*a_t)*getmsigma.M_ptr[z1+N1*i]*F.rho33_access(site(z1+n-m), z2, qperp, t, tt);
     i++;
         }
       }
     }
     #endif
        
     dxdtF31[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F31_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.F11_access(z1, z2, qperp, t, tt)-F.F11_access(site(z1-1), z2, qperp, t, tt))+buffer31+memint.SigmarhoF31_access(z1, z2, qperp)+memint.SigmaFP31_access(z1, z2, qperp);
     dxdtF32[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F32_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.F12_access(z1, z2, qperp, t, tt)-F.F12_access(site(z1-1), z2, qperp, t, tt))+buffer32+memint.SigmarhoF32_access(z1, z2, qperp)+memint.SigmaFP32_access(z1, z2, qperp);
     dxdtF33[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.F33_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.F13_access(z1, z2, qperp, t, tt)-F.F13_access(site(z1-1), z2, qperp, t, tt))+buffer33+memint.SigmarhoF33_access(z1, z2, qperp)+memint.SigmaFP33_access(z1, z2, qperp);
     
     dxdtrho31[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho31_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rho11_access(z1, z2, qperp, t, tt)-F.rho11_access(site(z1-1), z2, qperp, t, tt))+bufferrho31+memint.SigmarhoP31_access(z1, z2, qperp);
     dxdtrho32[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho32_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rho12_access(z1, z2, qperp, t, tt)-F.rho12_access(site(z1-1), z2, qperp, t, tt))+bufferrho32+memint.SigmarhoP32_access(z1, z2, qperp);
     dxdtrho33[z1+ z2*N1+qperp*N1*N1]=-4./(a_t*a_t)*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))*F.rho33_access(z1, z2, qperp, t, tt)-1./a_t*(exp(M_I*2.*M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-1.)*1./a*(F.rho13_access(z1, z2, qperp, t, tt)-F.rho13_access(site(z1-1), z2, qperp, t, tt))+bufferrho33+memint.SigmarhoP33_access(z1, z2, qperp);
  
  //if((z1==10)&&(z2==19)&&(qperp==1)){cout<<"M:"<<buffer31<<endl; cout<<"sigmarhoF31:"<<memint.SigmarhoF31_access(z1, z2, qperp)<<endl; cout<<"sigmaFP31:"<<memint.SigmaFP31_access(z1, z2, qperp)<<endl;}
     }
     
    }
   }
  }
   
  
  }  
  
 
  
  




};

#endif		// TIME_EVOLUTION_H
