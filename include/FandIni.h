#ifndef FANDINI_H
#define FANDINI_H

#include <complex>
#include <iostream>
#include <cmath>
using namespace std;

class FandIni
{
public:


complex<double>* F11;
complex<double>* F12;
complex<double>* F13;
complex<double>* F21;
complex<double>* F22;
complex<double>* F23;
complex<double>* F31;
complex<double>* F32;
complex<double>* F33;

complex<double>* Ft11;
complex<double>* Ft12;
complex<double>* Ft13;
complex<double>* Ft21;
complex<double>* Ft22;
complex<double>* Ft23;
complex<double>* Ft31;
complex<double>* Ft32;
complex<double>* Ft33;

complex<double>* rho11;
complex<double>* rho12;
complex<double>* rho13;
complex<double>* rho21;
complex<double>* rho22;
complex<double>* rho23;
complex<double>* rho31;
complex<double>* rho32;
complex<double>* rho33;

complex<double>* rhot11;
complex<double>* rhot12;
complex<double>* rhot13;
complex<double>* rhot21;
complex<double>* rhot22;
complex<double>* rhot23;
complex<double>* rhot31;
complex<double>* rhot32;
complex<double>* rhot33;

complex<double>* Fa11;
complex<double>* Fa12;
complex<double>* Fa13;
complex<double>* Fa21;
complex<double>* Fa22;
complex<double>* Fa23;
complex<double>* Fa31;
complex<double>* Fa32;
complex<double>* Fa33;

complex<double>* Fb11;
complex<double>* Fb12;
complex<double>* Fb13;
complex<double>* Fb21;
complex<double>* Fb22;
complex<double>* Fb23;
complex<double>* Fb31;
complex<double>* Fb32;
complex<double>* Fb33;

complex<double>* Fc11;
complex<double>* Fc12;
complex<double>* Fc13;
complex<double>* Fc21;
complex<double>* Fc22;
complex<double>* Fc23;
complex<double>* Fc31;
complex<double>* Fc32;
complex<double>* Fc33;

complex<double>* Fd11;
complex<double>* Fd12;
complex<double>* Fd13;
complex<double>* Fd21;
complex<double>* Fd22;
complex<double>* Fd23;
complex<double>* Fd31;
complex<double>* Fd32;
complex<double>* Fd33;

complex<double>* Fe11;
complex<double>* Fe12;
complex<double>* Fe13;
complex<double>* Fe21;
complex<double>* Fe22;
complex<double>* Fe23;
complex<double>* Fe31;
complex<double>* Fe32;
complex<double>* Fe33;

complex<double>* Ff11;
complex<double>* Ff12;
complex<double>* Ff13;
complex<double>* Ff21;
complex<double>* Ff22;
complex<double>* Ff23;
complex<double>* Ff31;
complex<double>* Ff32;
complex<double>* Ff33;

FandIni()
{


F11=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F12=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F13=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F21=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F22=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F23=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F31=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F32=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
F33=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];

Ft11=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft12=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft13=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft21=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft22=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft23=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft31=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft32=new complex<double>[N1*N1*Nt*n_times_selfenergies];
Ft33=new complex<double>[N1*N1*Nt*n_times_selfenergies];

rho11=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho12=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho13=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho21=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho22=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho23=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho31=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho32=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];
rho33=new complex<double>[N1*N1*Nt*n_times_selfenergies*2];

rhot11=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot12=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot13=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot21=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot22=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot23=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot31=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot32=new complex<double>[N1*N1*Nt*n_times_selfenergies];
rhot33=new complex<double>[N1*N1*Nt*n_times_selfenergies];

Fa11=new complex<double>[N1*N1*Nt];
Fa12=new complex<double>[N1*N1*Nt];
Fa13=new complex<double>[N1*N1*Nt];
Fa21=new complex<double>[N1*N1*Nt];
Fa22=new complex<double>[N1*N1*Nt];
Fa23=new complex<double>[N1*N1*Nt];
Fa31=new complex<double>[N1*N1*Nt];
Fa32=new complex<double>[N1*N1*Nt];
Fa33=new complex<double>[N1*N1*Nt];

Fb11=new complex<double>[N1*N1*Nt];
Fb12=new complex<double>[N1*N1*Nt];
Fb13=new complex<double>[N1*N1*Nt];
Fb21=new complex<double>[N1*N1*Nt];
Fb22=new complex<double>[N1*N1*Nt];
Fb23=new complex<double>[N1*N1*Nt];
Fb31=new complex<double>[N1*N1*Nt];
Fb32=new complex<double>[N1*N1*Nt];
Fb33=new complex<double>[N1*N1*Nt];

Fc11=new complex<double>[N1*N1*Nt];
Fc12=new complex<double>[N1*N1*Nt];
Fc13=new complex<double>[N1*N1*Nt];
Fc21=new complex<double>[N1*N1*Nt];
Fc22=new complex<double>[N1*N1*Nt];
Fc23=new complex<double>[N1*N1*Nt];
Fc31=new complex<double>[N1*N1*Nt];
Fc32=new complex<double>[N1*N1*Nt];
Fc33=new complex<double>[N1*N1*Nt];

Fd11=new complex<double>[N1*N1*Nt];
Fd12=new complex<double>[N1*N1*Nt];
Fd13=new complex<double>[N1*N1*Nt];
Fd21=new complex<double>[N1*N1*Nt];
Fd22=new complex<double>[N1*N1*Nt];
Fd23=new complex<double>[N1*N1*Nt];
Fd31=new complex<double>[N1*N1*Nt];
Fd32=new complex<double>[N1*N1*Nt];
Fd33=new complex<double>[N1*N1*Nt];

Fe11=new complex<double>[N1*N1*Nt];
Fe12=new complex<double>[N1*N1*Nt];
Fe13=new complex<double>[N1*N1*Nt];
Fe21=new complex<double>[N1*N1*Nt];
Fe22=new complex<double>[N1*N1*Nt];
Fe23=new complex<double>[N1*N1*Nt];
Fe31=new complex<double>[N1*N1*Nt];
Fe32=new complex<double>[N1*N1*Nt];
Fe33=new complex<double>[N1*N1*Nt];

Ff11=new complex<double>[N1*N1*Nt];
Ff12=new complex<double>[N1*N1*Nt];
Ff13=new complex<double>[N1*N1*Nt];
Ff21=new complex<double>[N1*N1*Nt];
Ff22=new complex<double>[N1*N1*Nt];
Ff23=new complex<double>[N1*N1*Nt];
Ff31=new complex<double>[N1*N1*Nt];
Ff32=new complex<double>[N1*N1*Nt];
Ff33=new complex<double>[N1*N1*Nt];
}

~FandIni()
{
  delete [] F11;
  delete [] F12;
  delete [] F13;
  delete [] F21;
  delete [] F22;
  delete [] F23;
  delete [] F31;
  delete [] F32;
  delete [] F33;
  
  delete [] Ft11;
  delete [] Ft12;
  delete [] Ft13;
  delete [] Ft21;
  delete [] Ft22;
  delete [] Ft23;
  delete [] Ft31;
  delete [] Ft32;
  delete [] Ft33;

  delete [] rho11;
  delete [] rho12;
  delete [] rho13;
  delete [] rho21;
  delete [] rho22;
  delete [] rho23;
  delete [] rho31;
  delete [] rho32;
  delete [] rho33;
  
  delete [] rhot11;
  delete [] rhot12;
  delete [] rhot13;
  delete [] rhot21;
  delete [] rhot22;
  delete [] rhot23;
  delete [] rhot31;
  delete [] rhot32;
  delete [] rhot33;
  
  delete [] Fa11;
  delete [] Fa12;
  delete [] Fa13;
  delete [] Fa21;
  delete [] Fa22;
  delete [] Fa23;
  delete [] Fa31;
  delete [] Fa32;
  delete [] Fa33;
  
  delete [] Fb11;
  delete [] Fb12;
  delete [] Fb13;
  delete [] Fb21;
  delete [] Fb22;
  delete [] Fb23;
  delete [] Fb31;
  delete [] Fb32;
  delete [] Fb33;
  
  delete [] Fc11;
  delete [] Fc12;
  delete [] Fc13;
  delete [] Fc21;
  delete [] Fc22;
  delete [] Fc23;
  delete [] Fc31;
  delete [] Fc32;
  delete [] Fc33;
  
  delete [] Fd11;
  delete [] Fd12;
  delete [] Fd13;
  delete [] Fd21;
  delete [] Fd22;
  delete [] Fd23;
  delete [] Fd31;
  delete [] Fd32;
  delete [] Fd33;
  
  delete [] Fe11;
  delete [] Fe12;
  delete [] Fe13;
  delete [] Fe21;
  delete [] Fe22;
  delete [] Fe23;
  delete [] Fe31;
  delete [] Fe32;
  delete [] Fe33;
  
  delete [] Ff11;
  delete [] Ff12;
  delete [] Ff13;
  delete [] Ff21;
  delete [] Ff22;
  delete [] Ff23;
  delete [] Ff31;
  delete [] Ff32;
  delete [] Ff33;
}
complex<double> abstildep(int64_t qperp, int64_t q3)
{
 complex<double> buffer;
 buffer=sqrt(((2./a_t)*(cos(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-M_I*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*((2./a_t)*(cos(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))+M_I*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))+((2./a)*(cos(M_PI*(double(q3)-double(N1)/2.)/double(N1))-M_I*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1)))*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1)))*((2./a)*(cos(M_PI*(double(q3)-double(N1)/2.)/double(N1))+M_I*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1)))*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1))));
 return buffer;
}
complex<double> tildept(int64_t qperp)
{
 complex<double> buffer;
 buffer=(2./a_t)*(cos(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt))-M_I*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt)))*sin(M_PI*(double(qperp)-double(Nt)/2.)/double(Nt));
 return buffer;
}
complex<double> tildep3(int64_t q3)
{
 complex<double> buffer;
 buffer=(2./a)*(cos(M_PI*(double(q3)-double(N1)/2.)/double(N1))-M_I*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1)))*sin(M_PI*(double(q3)-double(N1)/2.)/double(N1));
 return buffer;
}

double omega(int64_t qperp, int64_t q3)
{
 double buffer;
 buffer=2./(dt*double(n_out_selfenergies))*asin(sqrt(dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)/(a*a)*sin(a*(2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1)*a))/2.)*sin(a*(2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1)*a))/2.)+dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)/(a_t*a_t)*sin(a_t*(2.*M_PI*(double(qperp)-double(Nt)/2.)/(double(Nt)*a_t))/2.)*sin(a_t*(2.*M_PI*(double(qperp)-double(Nt)/2.)/(double(Nt)*a_t))/2.)));
 return buffer;
}
double nq_ini_box(int64_t qperp, int64_t q3)
{
 double buffer;
 if(omega(qperp, q3)<Q)
 {
 buffer=1./(e*e);
 }
 else
 {
 buffer=0.;
 }
 return buffer;
}
       int64_t siteqperp(int64_t n)       //implement periodic bounndary condition
       {
         int64_t nn;
         if(n>-1 && n<Nt) {nn=n;}
         else if (n>Nt-1){nn=n-Nt;}
         else if (n<0) {nn=Nt+n;}
         return (nn);
       } 

void initialize()
{
for(int64_t n=0; n<N1*N1*Nt*n_times_selfenergies*2; n++)
{
  
      F11[n]=0.;
      F12[n]=0.;
      F13[n]=0.;
      F21[n]=0.;
      F22[n]=0.;
      F23[n]=0.;
      F31[n]=0.;
      F32[n]=0.;
      F33[n]=0.;
      
      rho11[n]=0.;
      rho12[n]=0.;
      rho13[n]=0.;
      rho21[n]=0.;
      rho22[n]=0.;
      rho23[n]=0.;
      rho31[n]=0.;
      rho32[n]=0.;
      rho33[n]=0.;

}

for(int64_t qperp=0; qperp<Nt; qperp++)
{
  for(int64_t z1=0; z1<N1; z1++)
  {
    for(int64_t z2=0; z2<N1; z2++)
    {
      for(int64_t q3=0; q3<N1; q3++)
      {
      if((qperp!=(Nt/2))||(q3!=(N1/2)))
      {
      #if (INIT==0)
      //t=0 and t'=0
   F11[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
   rho11[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));

      //t=1 and t'=0
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
      
   
     //t=1 and t'=1 
     
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
 
      //t=0 and t'=1
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     #elif (INIT==1)
      //t=0 and t'=0
   F11[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
    
     F12[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F21[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F23[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F32[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
   rho11[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
    
     rho12[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho21[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho23[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho32[z1+ z2*N1+qperp*N1*N1]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));

      //t=1 and t'=0
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
      
   
     //t=1 and t'=1 
     
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*(exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))))*(1.+2.*nq_ini_box(qperp, q3));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]+=((double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))));
 
      //t=0 and t'=1
     F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(4.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))+(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies))+2.*nq_ini_box(qperp, q3)*cos(omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
    
     rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildept(qperp)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildep3(q3)*conj(tildept(siteqperp(Nt-qperp)))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*(M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((-tildep3(q3)*conj(tildept(qperp))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(-tildept(siteqperp(Nt-qperp))*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     
     rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=0.;
     
     rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1]+=(double(Nt*Nt)*a_t*a_t)/(2.*double(N1)*a)*(1./abstildep(qperp, q3))*((1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2))-(1.-tildep3(q3)*conj(tildep3(q3))/(abstildep(qperp, q3)*abstildep(qperp, q3)))*M_I*(exp(-M_I*omega(qperp,q3)*dt*double(n_out_selfenergies)))*exp(-M_I*2.*M_PI*(double(q3)-double(N1)/2.)/(double(N1))*double(z1-z2)));
     #endif
     
      }     
      }
    }
  }
 }
 
for(int64_t qperp=0; qperp<Nt; qperp++)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
        
           Fc11[z1+ z2*N1+qperp*N1*N1]=F11[z1+ z2*N1+qperp*N1*N1];
           Fc12[z1+ z2*N1+qperp*N1*N1]=F12[z1+ z2*N1+qperp*N1*N1];
           Fc13[z1+ z2*N1+qperp*N1*N1]=F13[z1+ z2*N1+qperp*N1*N1];
           Fc21[z1+ z2*N1+qperp*N1*N1]=F21[z1+ z2*N1+qperp*N1*N1];
           Fc22[z1+ z2*N1+qperp*N1*N1]=F22[z1+ z2*N1+qperp*N1*N1];
           Fc23[z1+ z2*N1+qperp*N1*N1]=F23[z1+ z2*N1+qperp*N1*N1];
           Fc31[z1+ z2*N1+qperp*N1*N1]=F31[z1+ z2*N1+qperp*N1*N1];
           Fc32[z1+ z2*N1+qperp*N1*N1]=F32[z1+ z2*N1+qperp*N1*N1];
           Fc33[z1+ z2*N1+qperp*N1*N1]=F33[z1+ z2*N1+qperp*N1*N1];
           
           Fe11[z1+ z2*N1+qperp*N1*N1]=F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe12[z1+ z2*N1+qperp*N1*N1]=F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe13[z1+ z2*N1+qperp*N1*N1]=F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe21[z1+ z2*N1+qperp*N1*N1]=F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe22[z1+ z2*N1+qperp*N1*N1]=F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe23[z1+ z2*N1+qperp*N1*N1]=F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe31[z1+ z2*N1+qperp*N1*N1]=F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe32[z1+ z2*N1+qperp*N1*N1]=F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           Fe33[z1+ z2*N1+qperp*N1*N1]=F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*2];
           
           Ff11[z1+ z2*N1+qperp*N1*N1]=F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff12[z1+ z2*N1+qperp*N1*N1]=F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff13[z1+ z2*N1+qperp*N1*N1]=F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff21[z1+ z2*N1+qperp*N1*N1]=F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff22[z1+ z2*N1+qperp*N1*N1]=F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff23[z1+ z2*N1+qperp*N1*N1]=F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff31[z1+ z2*N1+qperp*N1*N1]=F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff32[z1+ z2*N1+qperp*N1*N1]=F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
           Ff33[z1+ z2*N1+qperp*N1*N1]=F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2];
        
           }
         }
       }

}
void acess_ini(int64_t tt)
{
for(int64_t n=0; n<N1*N1*Nt*n_times_selfenergies; n++)
{
  
      Ft11[n]=0.;
      Ft12[n]=0.;
      Ft13[n]=0.;
      Ft21[n]=0.;
      Ft22[n]=0.;
      Ft23[n]=0.;
      Ft31[n]=0.;
      Ft32[n]=0.;
      Ft33[n]=0.;
      
        
      rhot11[n]=0.;
      rhot12[n]=0.;
      rhot13[n]=0.;
      rhot21[n]=0.;
      rhot22[n]=0.;
      rhot23[n]=0.;
      rhot31[n]=0.;
      rhot32[n]=0.;
      rhot33[n]=0.;

}
for(int64_t t=0; t<2; t++)
{
for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
        
           (Ft11[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F11[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft12[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F12[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft13[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F13[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft21[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F21[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft22[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F22[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft23[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F23[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft31[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F31[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft32[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F32[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (Ft33[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(F33[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           
           (rhot11[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho11[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot12[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho12[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot13[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho13[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot21[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho21[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot22[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho22[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot23[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho23[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot31[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho31[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot32[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho32[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);
           (rhot33[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1])=(rho33[z1+ z2*N1+qperp*N1*N1+t*Nt*N1*N1+(tt)*Nt*N1*N1*2]);           
        
           }
         }
       }
}
if(tt==2)
{

for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
           (rhot11[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot12[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot13[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot21[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot22[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot23[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot31[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot32[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);
           (rhot33[z1+ z2*N1+qperp*N1*N1+(tt)*Nt*N1*N1])=(rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1+Nt*N1*N1*2]);           
        
           }
         }
       }

}
else
{
for(int64_t t=0; t<2; t++)
{
for(int64_t qperp=0; qperp<Nt; qperp++)// providing value for F with tt>t using the the symmetry property of F: Fij(px,py,z,t,z',t')=Fji(-px,-py,z',t',z,t)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
           
           (rhot11[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho11[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot12[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho12[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot13[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho13[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot21[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho21[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot22[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho22[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot23[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho23[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot31[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho31[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot32[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho32[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);
           (rhot33[z1+ z2*N1+qperp*N1*N1+(tt-t)*Nt*N1*N1])=(rho33[z1+ z2*N1+qperp*N1*N1+(1-t)*Nt*N1*N1+Nt*N1*N1*2]);           
        
           }
         }
       }
}
}
}

void diagupdate(int64_t t)
{

for(int64_t qperp=0; qperp<Nt; qperp++)
       {
         for(int64_t z1=0; z1<N1; z1++)
         {
           for(int64_t z2=0; z2<N1; z2++)
           {  
           Fa11[z1+ z2*N1+qperp*N1*N1]=Fc11[z1+ z2*N1+qperp*N1*N1];
           Fa12[z1+ z2*N1+qperp*N1*N1]=Fc12[z1+ z2*N1+qperp*N1*N1];
           Fa13[z1+ z2*N1+qperp*N1*N1]=Fc13[z1+ z2*N1+qperp*N1*N1];
           Fa21[z1+ z2*N1+qperp*N1*N1]=Fc21[z1+ z2*N1+qperp*N1*N1];
           Fa22[z1+ z2*N1+qperp*N1*N1]=Fc22[z1+ z2*N1+qperp*N1*N1];
           Fa23[z1+ z2*N1+qperp*N1*N1]=Fc23[z1+ z2*N1+qperp*N1*N1];
           Fa31[z1+ z2*N1+qperp*N1*N1]=Fc31[z1+ z2*N1+qperp*N1*N1];
           Fa32[z1+ z2*N1+qperp*N1*N1]=Fc32[z1+ z2*N1+qperp*N1*N1];
           Fa33[z1+ z2*N1+qperp*N1*N1]=Fc33[z1+ z2*N1+qperp*N1*N1];
           
           Fb11[z1+ z2*N1+qperp*N1*N1]=Fe11[z1+ z2*N1+qperp*N1*N1];
           Fb12[z1+ z2*N1+qperp*N1*N1]=Fe12[z1+ z2*N1+qperp*N1*N1];
           Fb13[z1+ z2*N1+qperp*N1*N1]=Fe13[z1+ z2*N1+qperp*N1*N1];
           Fb21[z1+ z2*N1+qperp*N1*N1]=Fe21[z1+ z2*N1+qperp*N1*N1];
           Fb22[z1+ z2*N1+qperp*N1*N1]=Fe22[z1+ z2*N1+qperp*N1*N1];
           Fb23[z1+ z2*N1+qperp*N1*N1]=Fe23[z1+ z2*N1+qperp*N1*N1];
           Fb31[z1+ z2*N1+qperp*N1*N1]=Fe31[z1+ z2*N1+qperp*N1*N1];
           Fb32[z1+ z2*N1+qperp*N1*N1]=Fe32[z1+ z2*N1+qperp*N1*N1];
           Fb33[z1+ z2*N1+qperp*N1*N1]=Fe33[z1+ z2*N1+qperp*N1*N1];
           
           Fc11[z1+ z2*N1+qperp*N1*N1]=Ff11[z1+ z2*N1+qperp*N1*N1];
           Fc12[z1+ z2*N1+qperp*N1*N1]=Ff12[z1+ z2*N1+qperp*N1*N1];
           Fc13[z1+ z2*N1+qperp*N1*N1]=Ff13[z1+ z2*N1+qperp*N1*N1];
           Fc21[z1+ z2*N1+qperp*N1*N1]=Ff21[z1+ z2*N1+qperp*N1*N1];
           Fc22[z1+ z2*N1+qperp*N1*N1]=Ff22[z1+ z2*N1+qperp*N1*N1];
           Fc23[z1+ z2*N1+qperp*N1*N1]=Ff23[z1+ z2*N1+qperp*N1*N1];
           Fc31[z1+ z2*N1+qperp*N1*N1]=Ff31[z1+ z2*N1+qperp*N1*N1];
           Fc32[z1+ z2*N1+qperp*N1*N1]=Ff32[z1+ z2*N1+qperp*N1*N1];
           Fc33[z1+ z2*N1+qperp*N1*N1]=Ff33[z1+ z2*N1+qperp*N1*N1];
           
           Fd11[z1+ z2*N1+qperp*N1*N1]=Ft11[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd12[z1+ z2*N1+qperp*N1*N1]=Ft12[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd13[z1+ z2*N1+qperp*N1*N1]=Ft13[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd21[z1+ z2*N1+qperp*N1*N1]=Ft21[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd22[z1+ z2*N1+qperp*N1*N1]=Ft22[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd23[z1+ z2*N1+qperp*N1*N1]=Ft23[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd31[z1+ z2*N1+qperp*N1*N1]=Ft31[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd32[z1+ z2*N1+qperp*N1*N1]=Ft32[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           Fd33[z1+ z2*N1+qperp*N1*N1]=Ft33[z1+ z2*N1+qperp*N1*N1+(t-2)*Nt*N1*N1];
           
           Fe11[z1+ z2*N1+qperp*N1*N1]=Ft11[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe12[z1+ z2*N1+qperp*N1*N1]=Ft12[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe13[z1+ z2*N1+qperp*N1*N1]=Ft13[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe21[z1+ z2*N1+qperp*N1*N1]=Ft21[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe22[z1+ z2*N1+qperp*N1*N1]=Ft22[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe23[z1+ z2*N1+qperp*N1*N1]=Ft23[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe31[z1+ z2*N1+qperp*N1*N1]=Ft31[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe32[z1+ z2*N1+qperp*N1*N1]=Ft32[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           Fe33[z1+ z2*N1+qperp*N1*N1]=Ft33[z1+ z2*N1+qperp*N1*N1+(t-1)*Nt*N1*N1];
           
           Ff11[z1+ z2*N1+qperp*N1*N1]=Ft11[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff12[z1+ z2*N1+qperp*N1*N1]=Ft12[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff13[z1+ z2*N1+qperp*N1*N1]=Ft13[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff21[z1+ z2*N1+qperp*N1*N1]=Ft21[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff22[z1+ z2*N1+qperp*N1*N1]=Ft22[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff23[z1+ z2*N1+qperp*N1*N1]=Ft23[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff31[z1+ z2*N1+qperp*N1*N1]=Ft31[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff32[z1+ z2*N1+qperp*N1*N1]=Ft32[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];
           Ff33[z1+ z2*N1+qperp*N1*N1]=Ft33[z1+ z2*N1+qperp*N1*N1+(t)*Nt*N1*N1];          
           }
         }
       }
}

       
complex<double> F11_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F12_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F13_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F21_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F22_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F23_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F31_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F32_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> F33_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return F33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
}  

complex<double> Ft11_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft12_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft13_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft21_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft22_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft23_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft31_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft32_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> Ft33_access(int z1, int z2, int qperp, int64_t t1)
{
return Ft33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
}  

complex<double> rho11_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho12_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho13_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho21_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho22_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho23_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho31_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho32_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
} 

complex<double> rho33_access(int z1, int z2, int qperp, int64_t t1, int64_t t2)
{
return rho33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t2+Nt*N1*N1*2*t1];
}  

complex<double> rhot11_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot11[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot12_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot12[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot13_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot13[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot21_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot21[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot22_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot22[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot23_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot23[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot31_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot31[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot32_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot32[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
} 

complex<double> rhot33_access(int z1, int z2, int qperp, int64_t t1)
{
return rhot33[z1+ z2*N1+qperp*N1*N1+Nt*N1*N1*t1];
}  

complex<double> Fa11_access(int z1, int z2, int qperp)
{
return Fa11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa12_access(int z1, int z2, int qperp)
{
return Fa12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa13_access(int z1, int z2, int qperp)
{
return Fa13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa21_access(int z1, int z2, int qperp)
{
return Fa21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa22_access(int z1, int z2, int qperp)
{
return Fa22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa23_access(int z1, int z2, int qperp)
{
return Fa23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa31_access(int z1, int z2, int qperp)
{
return Fa31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa32_access(int z1, int z2, int qperp)
{
return Fa32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fa33_access(int z1, int z2, int qperp)
{
return Fa33[z1+ z2*N1+qperp*N1*N1];
}  

complex<double> Fb11_access(int z1, int z2, int qperp)
{
return Fb11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb12_access(int z1, int z2, int qperp)
{
return Fb12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb13_access(int z1, int z2, int qperp)
{
return Fb13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb21_access(int z1, int z2, int qperp)
{
return Fb21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb22_access(int z1, int z2, int qperp)
{
return Fb22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb23_access(int z1, int z2, int qperp)
{
return Fb23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb31_access(int z1, int z2, int qperp)
{
return Fb31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb32_access(int z1, int z2, int qperp)
{
return Fb32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fb33_access(int z1, int z2, int qperp)
{
return Fb33[z1+ z2*N1+qperp*N1*N1];
}  

complex<double> Fc11_access(int z1, int z2, int qperp)
{
return Fc11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc12_access(int z1, int z2, int qperp)
{
return Fc12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc13_access(int z1, int z2, int qperp)
{
return Fc13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc21_access(int z1, int z2, int qperp)
{
return Fc21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc22_access(int z1, int z2, int qperp)
{
return Fc22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc23_access(int z1, int z2, int qperp)
{
return Fc23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc31_access(int z1, int z2, int qperp)
{
return Fc31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc32_access(int z1, int z2, int qperp)
{
return Fc32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fc33_access(int z1, int z2, int qperp)
{
return Fc33[z1+ z2*N1+qperp*N1*N1];
}  

complex<double> Fd11_access(int z1, int z2, int qperp)
{
return Fd11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd12_access(int z1, int z2, int qperp)
{
return Fd12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd13_access(int z1, int z2, int qperp)
{
return Fd13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd21_access(int z1, int z2, int qperp)
{
return Fd21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd22_access(int z1, int z2, int qperp)
{
return Fd22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd23_access(int z1, int z2, int qperp)
{
return Fd23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd31_access(int z1, int z2, int qperp)
{
return Fd31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd32_access(int z1, int z2, int qperp)
{
return Fd32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fd33_access(int z1, int z2, int qperp)
{
return Fd33[z1+ z2*N1+qperp*N1*N1];
}  

complex<double> Fe11_access(int z1, int z2, int qperp)
{
return Fe11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe12_access(int z1, int z2, int qperp)
{
return Fe12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe13_access(int z1, int z2, int qperp)
{
return Fe13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe21_access(int z1, int z2, int qperp)
{
return Fe21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe22_access(int z1, int z2, int qperp)
{
return Fe22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe23_access(int z1, int z2, int qperp)
{
return Fe23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe31_access(int z1, int z2, int qperp)
{
return Fe31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe32_access(int z1, int z2, int qperp)
{
return Fe32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Fe33_access(int z1, int z2, int qperp)
{
return Fe33[z1+ z2*N1+qperp*N1*N1];
}  

complex<double> Ff11_access(int z1, int z2, int qperp)
{
return Ff11[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff12_access(int z1, int z2, int qperp)
{
return Ff12[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff13_access(int z1, int z2, int qperp)
{
return Ff13[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff21_access(int z1, int z2, int qperp)
{
return Ff21[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff22_access(int z1, int z2, int qperp)
{
return Ff22[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff23_access(int z1, int z2, int qperp)
{
return Ff23[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff31_access(int z1, int z2, int qperp)
{
return Ff31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff32_access(int z1, int z2, int qperp)
{
return Ff32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> Ff33_access(int z1, int z2, int qperp)
{
return Ff33[z1+ z2*N1+qperp*N1*N1];
}  

private:



};

#endif // FANDINI_H

