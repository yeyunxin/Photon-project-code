#ifndef EULERMETHOD_H
#define EULERMETHOD_H


#include <complex.h>

using namespace std;

template < typename Func >
class MyEulerSolver
{
    public:
    
       int64_t siteqperp(int64_t n)       //implement periodic bounndary condition
       {
         int64_t nn;
         if(n>-1 && n<Nt) {nn=n;}
         else if (n>Nt-1){nn=n-Nt;}
         else if (n<0) {nn=Nt+n;}
         return (nn);
       } 
       void make_step(FandIni &F, GetMSigma &getmsigma, MemoryInt &memint, const int64_t t, const int64_t tt, int Direction)
       {
       rhs(dxdtF11, dxdtF12, dxdtF13, dxdtF21, dxdtF22, dxdtF23, dxdtF31, dxdtF32, dxdtF33,dxdtrho11, dxdtrho12, dxdtrho13, dxdtrho21, dxdtrho22, dxdtrho23, dxdtrho31, dxdtrho32, dxdtrho33, F, getmsigma, memint, t, tt, Direction);
       #pragma omp parallel for 
       for(int64_t z1=0; z1<N1; z1++)
       {
         for(int64_t z2=0; z2<N1; z2++)
         {
           for(int64_t qperp=0; qperp<Nt; qperp++)
           {
           
           (F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF11[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           (F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF12[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF13[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF21[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF22[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF23[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF31[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF32[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtF33[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           
           (F.rho11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho11[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           (F.rho12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho12[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho13[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho21[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho22[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho23[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho31[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho32[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rho33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*2])=2.*(F.rho33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+t*Nt*N1*N1*2])-(F.rho33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t-1)*Nt*N1*N1*2])+dxdtrho33[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
             
           
           }
         }
       }
      
      /* for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
          
           if((t+1)>tt)
           {
          (F.F11[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F21[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F31[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F12[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F22[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F32[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F13[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F23[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F33[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           
           }
         
           
           }
         }
       }*/
    
       }
       
       void make_stept_rev(FandIni &F, GetMSigma &getmsigma, MemoryInt &memint, const int64_t t, const int64_t tt, int Direction)
       {
       rhs(dxdtF11, dxdtF12, dxdtF13, dxdtF21, dxdtF22, dxdtF23, dxdtF31, dxdtF32, dxdtF33,dxdtrho11, dxdtrho12, dxdtrho13, dxdtrho21, dxdtrho22, dxdtrho23, dxdtrho31, dxdtrho32, dxdtrho33, F, getmsigma, memint, t, tt, Direction);
       #pragma omp parallel for
       for(int64_t z1=0; z1<N1; z1++)
       {
         for(int64_t z2=0; z2<N1; z2++)
         {
           for(int64_t qperp=0; qperp<Nt; qperp++)
           {
           
           (F.rhot11[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot11[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho11[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           (F.rhot12[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot12[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho12[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot13[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot13[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho13[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot21[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot21[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho21[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot22[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot22[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho22[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot23[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot23[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho23[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot31[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot31[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho31[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot32[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot32[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho32[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.rhot33[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])=2.*(F.rhot33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.rhot33[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])+dxdtrho33[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           }
         }
       }
      
      /* for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
          
           if((t+1)>tt)
           {
          (F.F11[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F21[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F31[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F12[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F22[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F32[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F13[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F23[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F33[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           
           }
         
           
           }
         }
       }*/
    
       }

       void make_stept(FandIni &F, GetMSigma &getmsigma, MemoryInt &memint, const int64_t t, const int64_t tt, int Direction)
       {
       rhs(dxdtF11, dxdtF12, dxdtF13, dxdtF21, dxdtF22, dxdtF23, dxdtF31, dxdtF32, dxdtF33, dxdtrho11, dxdtrho12, dxdtrho13, dxdtrho21, dxdtrho22, dxdtrho23, dxdtrho31, dxdtrho32, dxdtrho33, F, getmsigma, memint, t, tt, Direction);
       #pragma omp parallel for
       for(int64_t z1=0; z1<N1; z1++)
       {
         for(int64_t z2=0; z2<N1; z2++)
         {
           for(int64_t qperp=0; qperp<Nt; qperp++)
           {
           
           (F.Ft11[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft11[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF11[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           (F.Ft12[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft12[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF12[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft13[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft13[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF13[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft21[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft21[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF21[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft22[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft22[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF22[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft23[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft23[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF23[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft31[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft31[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF31[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft32[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft32[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF32[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
           (F.Ft33[z1+ z2*N1+qperp*N1*N1+(tt+1)*Nt*N1*N1])=2.*(F.Ft33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1])-(F.Ft33[z1+ z2*N1+qperp*N1*N1+(tt-1)*Nt*N1*N1])+dxdtF33[z1+ z2*N1+qperp*N1*N1]*dt * double(n_out_selfenergies)*dt * double(n_out_selfenergies);
        
           }
         }
       }
      
      /* for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
          
           if((t+1)>tt)
           {
          (F.F11[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F11[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F21[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F12[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F31[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F13[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F12[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F21[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F22[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F22[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F32[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F23[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F13[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F31[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F23[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F32[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           (F.F33[z2+ z1*N1+siteqperp(Nt-qperp)*N1*N1+(t+1)*Nt*N1*N1+tt*Nt*N1*N1*n_times_selfenergies])=(F.F33[z1+ z2*N1+qperp*N1*N1+tt*Nt*N1*N1+(t+1)*Nt*N1*N1*n_times_selfenergies]);
           
           }
         
           
           }
         }
       }*/
    
       }
	 
	MyEulerSolver()
        {

        dxdtF11=new complex<double>[N1*N1*Nt];
        dxdtF12=new complex<double>[N1*N1*Nt];
        dxdtF13=new complex<double>[N1*N1*Nt];
        dxdtF21=new complex<double>[N1*N1*Nt];
        dxdtF22=new complex<double>[N1*N1*Nt];
        dxdtF23=new complex<double>[N1*N1*Nt];
        dxdtF31=new complex<double>[N1*N1*Nt];
        dxdtF32=new complex<double>[N1*N1*Nt];
        dxdtF33=new complex<double>[N1*N1*Nt];
        dxdtrho11=new complex<double>[N1*N1*Nt];
        dxdtrho12=new complex<double>[N1*N1*Nt];
        dxdtrho13=new complex<double>[N1*N1*Nt];
        dxdtrho21=new complex<double>[N1*N1*Nt];
        dxdtrho22=new complex<double>[N1*N1*Nt];
        dxdtrho23=new complex<double>[N1*N1*Nt];
        dxdtrho31=new complex<double>[N1*N1*Nt];
        dxdtrho32=new complex<double>[N1*N1*Nt];
        dxdtrho33=new complex<double>[N1*N1*Nt];

        }
	
        ~MyEulerSolver() { }
      
      private:
	
      complex<double> *dxdtF11;
      complex<double> *dxdtF12;
      complex<double> *dxdtF13;
      complex<double> *dxdtF21;
      complex<double> *dxdtF22;
      complex<double> *dxdtF23;
      complex<double> *dxdtF31;
      complex<double> *dxdtF32;
      complex<double> *dxdtF33;
      complex<double> *dxdtrho11;
      complex<double> *dxdtrho12;
      complex<double> *dxdtrho13;
      complex<double> *dxdtrho21;
      complex<double> *dxdtrho22;
      complex<double> *dxdtrho23;
      complex<double> *dxdtrho31;
      complex<double> *dxdtrho32;
      complex<double> *dxdtrho33;
	Func rhs;
	

};

#endif // EULERMETHOD_H
