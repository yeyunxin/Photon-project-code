#ifndef OBSERVABLES_H
#define OBSERVABLES_H

#include <complex>
#include <iostream>
#include <cmath>

using namespace std;
class Observables
{
public:
complex<double>* dNoverVd3p;
complex<double>* F11qqa;
complex<double>* F12qqa;
complex<double>* F13qqa;
complex<double>* F21qqa;
complex<double>* F22qqa;
complex<double>* F23qqa;
complex<double>* F31qqa;
complex<double>* F32qqa;
complex<double>* F33qqa;

complex<double>* F11qqb;
complex<double>* F12qqb;
complex<double>* F13qqb;
complex<double>* F21qqb;
complex<double>* F22qqb;
complex<double>* F23qqb;
complex<double>* F31qqb;
complex<double>* F32qqb;
complex<double>* F33qqb;

complex<double>* F11qqc;
complex<double>* F12qqc;
complex<double>* F13qqc;
complex<double>* F21qqc;
complex<double>* F22qqc;
complex<double>* F23qqc;
complex<double>* F31qqc;
complex<double>* F32qqc;
complex<double>* F33qqc;

complex<double>* F11qqd;
complex<double>* F12qqd;
complex<double>* F13qqd;
complex<double>* F21qqd;
complex<double>* F22qqd;
complex<double>* F23qqd;
complex<double>* F31qqd;
complex<double>* F32qqd;
complex<double>* F33qqd;

complex<double>* F11qqe;
complex<double>* F12qqe;
complex<double>* F13qqe;
complex<double>* F21qqe;
complex<double>* F22qqe;
complex<double>* F23qqe;
complex<double>* F31qqe;
complex<double>* F32qqe;
complex<double>* F33qqe;

complex<double>* F11qqf;
complex<double>* F12qqf;
complex<double>* F13qqf;
complex<double>* F21qqf;
complex<double>* F22qqf;
complex<double>* F23qqf;
complex<double>* F31qqf;
complex<double>* F32qqf;
complex<double>* F33qqf;


Observables()
{

dNoverVd3p = new complex<double>[Nt*N1];
F11qqa = new complex<double>[Nt*N1];
F12qqa = new complex<double>[Nt*N1];
F13qqa = new complex<double>[Nt*N1];
F21qqa = new complex<double>[Nt*N1];
F22qqa = new complex<double>[Nt*N1];
F23qqa = new complex<double>[Nt*N1];
F31qqa = new complex<double>[Nt*N1];
F32qqa = new complex<double>[Nt*N1];
F33qqa = new complex<double>[Nt*N1];

F11qqb = new complex<double>[Nt*N1];
F12qqb = new complex<double>[Nt*N1];
F13qqb = new complex<double>[Nt*N1];
F21qqb = new complex<double>[Nt*N1];
F22qqb = new complex<double>[Nt*N1];
F23qqb = new complex<double>[Nt*N1];
F31qqb = new complex<double>[Nt*N1];
F32qqb = new complex<double>[Nt*N1];
F33qqb = new complex<double>[Nt*N1];

F11qqc = new complex<double>[Nt*N1];
F12qqc = new complex<double>[Nt*N1];
F13qqc = new complex<double>[Nt*N1];
F21qqc = new complex<double>[Nt*N1];
F22qqc = new complex<double>[Nt*N1];
F23qqc = new complex<double>[Nt*N1];
F31qqc = new complex<double>[Nt*N1];
F32qqc = new complex<double>[Nt*N1];
F33qqc = new complex<double>[Nt*N1];

F11qqd = new complex<double>[Nt*N1];
F12qqd = new complex<double>[Nt*N1];
F13qqd = new complex<double>[Nt*N1];
F21qqd = new complex<double>[Nt*N1];
F22qqd = new complex<double>[Nt*N1];
F23qqd = new complex<double>[Nt*N1];
F31qqd = new complex<double>[Nt*N1];
F32qqd = new complex<double>[Nt*N1];
F33qqd = new complex<double>[Nt*N1];

F11qqe = new complex<double>[Nt*N1];
F12qqe = new complex<double>[Nt*N1];
F13qqe = new complex<double>[Nt*N1];
F21qqe = new complex<double>[Nt*N1];
F22qqe = new complex<double>[Nt*N1];
F23qqe = new complex<double>[Nt*N1];
F31qqe = new complex<double>[Nt*N1];
F32qqe = new complex<double>[Nt*N1];
F33qqe = new complex<double>[Nt*N1];

F11qqf = new complex<double>[Nt*N1];
F12qqf = new complex<double>[Nt*N1];
F13qqf = new complex<double>[Nt*N1];
F21qqf = new complex<double>[Nt*N1];
F22qqf = new complex<double>[Nt*N1];
F23qqf = new complex<double>[Nt*N1];
F31qqf = new complex<double>[Nt*N1];
F32qqf = new complex<double>[Nt*N1];
F33qqf = new complex<double>[Nt*N1];


}

~Observables()
{

delete [] dNoverVd3p;
delete [] F11qqa;
delete [] F12qqa;
delete [] F13qqa;
delete [] F21qqa;
delete [] F22qqa;
delete [] F23qqa;
delete [] F31qqa;
delete [] F32qqa;
delete [] F33qqa;

delete [] F11qqb;
delete [] F12qqb;
delete [] F13qqb;
delete [] F21qqb;
delete [] F22qqb;
delete [] F23qqb;
delete [] F31qqb;
delete [] F32qqb;
delete [] F33qqb;

delete [] F11qqc;
delete [] F12qqc;
delete [] F13qqc;
delete [] F21qqc;
delete [] F22qqc;
delete [] F23qqc;
delete [] F31qqc;
delete [] F32qqc;
delete [] F33qqc;

delete [] F11qqd;
delete [] F12qqd;
delete [] F13qqd;
delete [] F21qqd;
delete [] F22qqd;
delete [] F23qqd;
delete [] F31qqd;
delete [] F32qqd;
delete [] F33qqd;

delete [] F11qqe;
delete [] F12qqe;
delete [] F13qqe;
delete [] F21qqe;
delete [] F22qqe;
delete [] F23qqe;
delete [] F31qqe;
delete [] F32qqe;
delete [] F33qqe;

delete [] F11qqf;
delete [] F12qqf;
delete [] F13qqf;
delete [] F21qqf;
delete [] F22qqf;
delete [] F23qqf;
delete [] F31qqf;
delete [] F32qqf;
delete [] F33qqf;

}
void ini_class_variables()
{
  for(int64_t i=0; i<N1*Nt; i++)
  {
    F11qqa[i]=0.;
    F12qqa[i]=0.;
    F13qqa[i]=0.;
    F21qqa[i]=0.;
    F22qqa[i]=0.;
    F23qqa[i]=0.;
    F31qqa[i]=0.;
    F32qqa[i]=0.;
    F33qqa[i]=0.;
    
    F11qqb[i]=0.;
    F12qqb[i]=0.;
    F13qqb[i]=0.;
    F21qqb[i]=0.;
    F22qqb[i]=0.;
    F23qqb[i]=0.;
    F31qqb[i]=0.;
    F32qqb[i]=0.;
    F33qqb[i]=0.;
    
    F11qqc[i]=0.;
    F12qqc[i]=0.;
    F13qqc[i]=0.;
    F21qqc[i]=0.;
    F22qqc[i]=0.;
    F23qqc[i]=0.;
    F31qqc[i]=0.;
    F32qqc[i]=0.;
    F33qqc[i]=0.;
   
    F11qqd[i]=0.;
    F12qqd[i]=0.;
    F13qqd[i]=0.;
    F21qqd[i]=0.;
    F22qqd[i]=0.;
    F23qqd[i]=0.;
    F31qqd[i]=0.;
    F32qqd[i]=0.;
    F33qqd[i]=0.;
    
    F11qqe[i]=0.;
    F12qqe[i]=0.;
    F13qqe[i]=0.;
    F21qqe[i]=0.;
    F22qqe[i]=0.;
    F23qqe[i]=0.;
    F31qqe[i]=0.;
    F32qqe[i]=0.;
    F33qqe[i]=0.;
    
    F11qqf[i]=0.;
    F12qqf[i]=0.;
    F13qqf[i]=0.;
    F21qqf[i]=0.;
    F22qqf[i]=0.;
    F23qqf[i]=0.;
    F31qqf[i]=0.;
    F32qqf[i]=0.;
    F33qqf[i]=0.;  
  }
}
       int64_t siteqperp(int64_t n)       //implement periodic bounndary condition
       {
         int64_t nn;
         if(n>-1 && n<Nt) {nn=n;}
         else if (n>Nt-1){nn=n-Nt;}
         else if (n<0) {nn=Nt+n;}
         return (nn);
       } 
       int64_t siteqz(int64_t n)       //implement periodic bounndary condition
       {
         int64_t nn;
         if(n>-1 && n<N1) {nn=n;}
         else if (n>N1-1){nn=n-N1;}
         else if (n<0) {nn=N1+n;}
         return (nn);
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
void compute_observables(FandIni &F, const int64_t t)
{


  ini_class_variables();
  #pragma omp parallel for 
  for(int64_t qz=0; qz<N1; qz++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      if((qperp!=(Nt/2))||(qz!=(N1/2)))
      {
      for(int64_t z1=0; z1<N1; z1++)
      {
        for(int64_t z2=0; z2<N1; z2++)
        {
   
          
       
          F11qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqa[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fa33_access(z1, z2, siteqperp(Nt-qperp));
          
          F11qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqb[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fb33_access(z1, z2, siteqperp(Nt-qperp));
          
          F11qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqc[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fc33_access(z1, z2, siteqperp(Nt-qperp));
          
          F11qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqd[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fd33_access(z1, z2, siteqperp(Nt-qperp));
          
          F11qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqe[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Fe33_access(z1, z2, siteqperp(Nt-qperp));
          
          F11qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff11_access(z1, z2, siteqperp(Nt-qperp));
          F12qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff12_access(z1, z2, siteqperp(Nt-qperp));
          F13qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff13_access(z1, z2, siteqperp(Nt-qperp));
          F21qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff21_access(z1, z2, siteqperp(Nt-qperp));
          F22qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff22_access(z1, z2, siteqperp(Nt-qperp));
          F23qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff23_access(z1, z2, siteqperp(Nt-qperp));
          F31qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff31_access(z1, z2, siteqperp(Nt-qperp));
          F32qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff32_access(z1, z2, siteqperp(Nt-qperp));
          F33qqf[qz+qperp*N1]+=a*a*(cos(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1))+M_I*sin(double(z1-z2)*2.*M_PI*(double(qz)-double(N1)/2.)/double(N1)))*F.Ff33_access(z1, z2, siteqperp(Nt-qperp));
        }
      }
      dNoverVd3p[qz+qperp*N1]=-1.+ 1./(2.*abstildep(qperp, qz)*a*double(N1)*a_t*a_t*double(Nt*Nt))*((abstildep(qperp, qz)*abstildep(qperp, qz)*F11qqc[qz+qperp*N1]+M_I*abstildep(qperp, qz)*(1./(2.*dt * double(n_out_selfenergies)))*(F11qqe[qz+qperp*N1]-F11qqb[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]-F11qqe[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]+F11qqb[qz+qperp*N1])+(1./(4.*dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)))*(F11qqf[qz+qperp*N1]+F11qqa[qz+qperp*N1]-F11qqd[qz+qperp*N1]-F11qqd[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]))*(1.-tildept(qperp)*conj(tildept(qperp))/(abstildep(qperp, qz)*abstildep(qperp, qz)))+(abstildep(qperp, qz)*abstildep(qperp, qz)*F13qqc[qz+qperp*N1]+M_I*abstildep(qperp, qz)*(1./(2.*dt * double(n_out_selfenergies)))*(F13qqe[qz+qperp*N1]-F31qqb[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]-F31qqe[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]+F13qqb[qz+qperp*N1])+(1./(4.*dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)))*(F13qqf[qz+qperp*N1]+F13qqa[qz+qperp*N1]-F13qqd[qz+qperp*N1]-F31qqd[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]))*(-tildept(qperp)*conj(tildep3(qz))/(abstildep(qperp, qz)*abstildep(qperp, qz)))+(abstildep(qperp, qz)*abstildep(qperp, qz)*F22qqc[qz+qperp*N1]+M_I*abstildep(qperp, qz)*(1./(2.*dt * double(n_out_selfenergies)))*(F22qqe[qz+qperp*N1]-F22qqb[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]-F22qqe[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]+F22qqb[qz+qperp*N1])+(1./(4.*dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)))*(F22qqf[qz+qperp*N1]+F22qqa[qz+qperp*N1]-F22qqd[qz+qperp*N1]-F22qqd[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]))*1.+(abstildep(qperp, qz)*abstildep(qperp, qz)*F31qqc[qz+qperp*N1]+M_I*abstildep(qperp, qz)*(1./(2.*dt * double(n_out_selfenergies)))*(F31qqe[qz+qperp*N1]-F13qqb[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]-F13qqe[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]+F31qqb[qz+qperp*N1])+(1./(4.*dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)))*(F31qqf[qz+qperp*N1]+F31qqa[qz+qperp*N1]-F31qqd[qz+qperp*N1]-F13qqd[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]))*(-tildep3(qz)*conj(tildept(qperp))/(abstildep(qperp, qz)*abstildep(qperp, qz)))+(abstildep(qperp, qz)*abstildep(qperp, qz)*F33qqc[qz+qperp*N1]+M_I*abstildep(qperp, qz)*(1./(2.*dt * double(n_out_selfenergies)))*(F33qqe[qz+qperp*N1]-F33qqb[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]-F33qqe[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]+F33qqb[qz+qperp*N1])+(1./(4.*dt*double(n_out_selfenergies)*dt*double(n_out_selfenergies)))*(F33qqf[qz+qperp*N1]+F33qqa[qz+qperp*N1]-F33qqd[qz+qperp*N1]-F33qqd[siteqz(N1-qz)+siteqperp(Nt-qperp)*N1]))*(1.-tildep3(qz)*conj(tildep3(qz))/(abstildep(qperp, qz)*abstildep(qperp, qz))));
    }
    }
  }
  

}
complex<double> F11qqa_access(int qz, int qperp)
{
return F11qqa[qz+qperp*N1];
} 
complex<double> F12qqa_access(int qz, int qperp)
{
return F12qqa[qz+qperp*N1];
} 
complex<double> F13qqa_access(int qz, int qperp)
{
return F13qqa[qz+qperp*N1];
} 
complex<double> F21qqa_access(int qz, int qperp)
{
return F21qqa[qz+qperp*N1];
} 
complex<double> F22qqa_access(int qz, int qperp)
{
return F22qqa[qz+qperp*N1];
} 
complex<double> F23qqa_access(int qz, int qperp)
{
return F23qqa[qz+qperp*N1];
} 
complex<double> F31qqa_access(int qz, int qperp)
{
return F31qqa[qz+qperp*N1];
} 
complex<double> F32qqa_access(int qz, int qperp)
{
return F32qqa[qz+qperp*N1];
} 
complex<double> F33qqa_access(int qz, int qperp)
{
return F33qqa[qz+qperp*N1];
} 

complex<double> F11qqb_access(int qz, int qperp)
{
return F11qqb[qz+qperp*N1];
} 
complex<double> F12qqb_access(int qz, int qperp)
{
return F12qqb[qz+qperp*N1];
} 
complex<double> F13qqb_access(int qz, int qperp)
{
return F13qqb[qz+qperp*N1];
} 
complex<double> F21qqb_access(int qz, int qperp)
{
return F21qqb[qz+qperp*N1];
} 
complex<double> F22qqb_access(int qz, int qperp)
{
return F22qqb[qz+qperp*N1];
} 
complex<double> F23qqb_access(int qz, int qperp)
{
return F23qqb[qz+qperp*N1];
} 
complex<double> F31qqb_access(int qz, int qperp)
{
return F31qqb[qz+qperp*N1];
} 
complex<double> F32qqb_access(int qz, int qperp)
{
return F32qqb[qz+qperp*N1];
} 
complex<double> F33qqb_access(int qz, int qperp)
{
return F33qqb[qz+qperp*N1];
} 

complex<double> F11qqc_access(int qz, int qperp)
{
return F11qqc[qz+qperp*N1];
} 
complex<double> F12qqc_access(int qz, int qperp)
{
return F12qqc[qz+qperp*N1];
} 
complex<double> F13qqc_access(int qz, int qperp)
{
return F13qqc[qz+qperp*N1];
} 
complex<double> F21qqc_access(int qz, int qperp)
{
return F21qqc[qz+qperp*N1];
} 
complex<double> F22qqc_access(int qz, int qperp)
{
return F22qqc[qz+qperp*N1];
} 
complex<double> F23qqc_access(int qz, int qperp)
{
return F23qqc[qz+qperp*N1];
} 
complex<double> F31qqc_access(int qz, int qperp)
{
return F31qqc[qz+qperp*N1];
} 
complex<double> F32qqc_access(int qz, int qperp)
{
return F32qqc[qz+qperp*N1];
} 
complex<double> F33qqc_access(int qz, int qperp)
{
return F33qqc[qz+qperp*N1];
} 

complex<double> F11qqd_access(int qz, int qperp)
{
return F11qqd[qz+qperp*N1];
} 
complex<double> F12qqd_access(int qz, int qperp)
{
return F12qqd[qz+qperp*N1];
} 
complex<double> F13qqd_access(int qz, int qperp)
{
return F13qqd[qz+qperp*N1];
} 
complex<double> F21qqd_access(int qz, int qperp)
{
return F21qqd[qz+qperp*N1];
} 
complex<double> F22qqd_access(int qz, int qperp)
{
return F22qqd[qz+qperp*N1];
} 
complex<double> F23qqd_access(int qz, int qperp)
{
return F23qqd[qz+qperp*N1];
} 
complex<double> F31qqd_access(int qz, int qperp)
{
return F31qqd[qz+qperp*N1];
} 
complex<double> F32qqd_access(int qz, int qperp)
{
return F32qqd[qz+qperp*N1];
} 
complex<double> F33qqd_access(int qz, int qperp)
{
return F33qqd[qz+qperp*N1];
} 

complex<double> F11qqe_access(int qz, int qperp)
{
return F11qqe[qz+qperp*N1];
} 
complex<double> F12qqe_access(int qz, int qperp)
{
return F12qqe[qz+qperp*N1];
} 
complex<double> F13qqe_access(int qz, int qperp)
{
return F13qqe[qz+qperp*N1];
} 
complex<double> F21qqe_access(int qz, int qperp)
{
return F21qqe[qz+qperp*N1];
} 
complex<double> F22qqe_access(int qz, int qperp)
{
return F22qqe[qz+qperp*N1];
} 
complex<double> F23qqe_access(int qz, int qperp)
{
return F23qqe[qz+qperp*N1];
} 
complex<double> F31qqe_access(int qz, int qperp)
{
return F31qqe[qz+qperp*N1];
} 
complex<double> F32qqe_access(int qz, int qperp)
{
return F32qqe[qz+qperp*N1];
} 
complex<double> F33qqe_access(int qz, int qperp)
{
return F33qqe[qz+qperp*N1];
} 

complex<double> F11qqf_access(int qz, int qperp)
{
return F11qqf[qz+qperp*N1];
} 
complex<double> F12qqf_access(int qz, int qperp)
{
return F12qqf[qz+qperp*N1];
} 
complex<double> F13qqf_access(int qz, int qperp)
{
return F13qqf[qz+qperp*N1];
} 
complex<double> F21qqf_access(int qz, int qperp)
{
return F21qqf[qz+qperp*N1];
} 
complex<double> F22qqf_access(int qz, int qperp)
{
return F22qqf[qz+qperp*N1];
} 
complex<double> F23qqf_access(int qz, int qperp)
{
return F23qqf[qz+qperp*N1];
} 
complex<double> F31qqf_access(int qz, int qperp)
{
return F31qqf[qz+qperp*N1];
} 
complex<double> F32qqf_access(int qz, int qperp)
{
return F32qqf[qz+qperp*N1];
} 
complex<double> F33qqf_access(int qz, int qperp)
{
return F33qqf[qz+qperp*N1];
} 

complex<double> dNoverVd3p_access(int qz, int qperp)
{
return dNoverVd3p[qz+qperp*N1];
} 
private:

};






#endif // OBSERVABLES_H

