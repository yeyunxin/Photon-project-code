#ifndef MEMORYINT_H
#define MEMORYINT_H

#include <complex>
#include <iostream>
#include <cmath>
#include <fftw3.h>
using namespace std;

class MemoryInt
{
public:
complex<double>* SigmarhoF31;
complex<double>* SigmarhoF32;
complex<double>* SigmarhoF33;
complex<double>* SigmaFP31;
complex<double>* SigmaFP32;
complex<double>* SigmaFP33;
complex<double>* SigmarhoP31;
complex<double>* SigmarhoP32;
complex<double>* SigmarhoP33;


MemoryInt()
{

SigmarhoF31=new complex<double>[N1*N1*Nt];
SigmarhoF32=new complex<double>[N1*N1*Nt];
SigmarhoF33=new complex<double>[N1*N1*Nt];
SigmaFP31=new complex<double>[N1*N1*Nt];
SigmaFP32=new complex<double>[N1*N1*Nt];
SigmaFP33=new complex<double>[N1*N1*Nt];
SigmarhoP31=new complex<double>[N1*N1*Nt];
SigmarhoP32=new complex<double>[N1*N1*Nt];
SigmarhoP33=new complex<double>[N1*N1*Nt];


}

~MemoryInt()
{
  delete [] SigmarhoF31;
  delete [] SigmarhoF32;
  delete [] SigmarhoF33;
  delete [] SigmaFP31;
  delete [] SigmaFP32;
  delete [] SigmaFP33;
  delete [] SigmarhoP31;
  delete [] SigmarhoP32;
  delete [] SigmarhoP33;
  
}

void process_mem_usage(double& vm_usage, double& resident_set)
{
   using std::ios_base;
   using std::ifstream;
   using std::string;

   vm_usage     = 0.0;
   resident_set = 0.0;

   // 'file' stat seems to give the most reliable results
   //
   ifstream stat_stream("/proc/self/stat",ios_base::in);

   // dummy vars for leading entries in stat that we don't care about
   //
   string pid, comm, state, ppid, pgrp, session, tty_nr;
   string tpgid, flags, minflt, cminflt, majflt, cmajflt;
   string utime, stime, cutime, cstime, priority, nice;
   string O, itrealvalue, starttime;

   // the two fields we want
   //
   unsigned long vsize;
   long rss;

   stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr
               >> tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt
               >> utime >> stime >> cutime >> cstime >> priority >> nice
               >> O >> itrealvalue >> starttime >> vsize >> rss; // don't care about the rest

   stat_stream.close();

   long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024; // in case x86-64 is configured to use 2MB pages
   vm_usage     = vsize / 1024.0;
   resident_set = rss * page_size_kb;
}

  int64_t site(int64_t n)       //implement periodic bounndary condition
  {
    int64_t nn;
    if(n>-1 && n<N1) {nn=n;}
    else if (n>N1-1){nn=n-N1;}
    else if (n<0) {nn=N1+n;}
    return (nn);
  } 

void integral(int64_t t,int64_t tt, FandIni &F, GetMSigma &getmsigma)
{


  using std::cout;
  using std::endl;

  double vm, rss;
  
  
for(int64_t n=0; n<N1*N1*Nt; n++)
{
SigmarhoF31[n]=0.;
SigmarhoF32[n]=0.;
SigmarhoF33[n]=0.;
SigmaFP31[n]=0.;
SigmaFP32[n]=0.;
SigmaFP33[n]=0.;
SigmarhoP31[n]=0.;
SigmarhoP32[n]=0.;
SigmarhoP33[n]=0.;
}
     
if(t<(time_spets_storage-1))
{
//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=0; ttt<(t+1); ttt++)
        {
        if((ttt>0)&&(ttt<t))
        {
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F33_access(zpp, z2, qperp, ttt, tt));
        
       
        }
        else
        {
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F33_access(zpp, z2, qperp, ttt, tt));        
        
        }
        
        
        }
      }
     
    }
  }
}
        /*std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        } */ 

}
else
{
//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=(t-(time_spets_storage-2)); ttt<(t+1); ttt++)
        {
        if((ttt>(t-(time_spets_storage-2)))&&(ttt<t))
        {
      
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F33_access(zpp, z2, qperp, ttt, tt));
       
        }
        else
        {
       
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.F33_access(zpp, z2, qperp, ttt, tt));        
        
       
        }
        
        
        
        
        }
      }
   
    }
  }
}
       /* std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        }  */

}


//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=0; ttt<(tt+1); ttt++)
        {
      
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z1, zpp, ttt))*(F.rho31_access(zpp, z2, qperp, ttt, tt));
        SigmaFP32[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z1, zpp, ttt))*(F.rho32_access(zpp, z2, qperp, ttt, tt));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z1, zpp, ttt))*(F.rho33_access(zpp, z2, qperp, ttt, tt));        
      
        }
      }
     
    }
  }
}

if((t-tt)<(time_spets_storage-1))
{
//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=tt; ttt<(t+1); ttt++)
        {
        if((ttt>tt)&&(ttt<t))
        {
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho33_access(zpp, z2, qperp, ttt, tt));
        
       
        }
        else
        {
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho33_access(zpp, z2, qperp, ttt, tt));        
        
        }
        
        
        }
      }
     
    }
  }
}
        /*std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        } */ 

}
else
{
//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=(t-(time_spets_storage-2)); ttt<(t+1); ttt++)
        {
        if((ttt>(t-(time_spets_storage-2)))&&(ttt<t))
        {
      
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho33_access(zpp, z2, qperp, ttt, tt));
       
        }
        else
        {
       
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho31_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho32_access(zpp, z2, qperp, ttt, tt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z1, zpp, ttt))*(F.rho33_access(zpp, z2, qperp, ttt, tt));        
        
       
        }
        
        
        
        
        }
      }
   
    }
  }
}
       /* std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        }  */

}

        /*std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        } */ 




        /*fftw_complex *in, *out;
        fftw_complex *in1, *out1;
        fftw_complex *INa, *OUTa;
        fftw_complex *IN1a, *OUT1a; 
        fftw_complex *INb, *OUTb;
        fftw_complex *IN1b, *OUT1b;        
        in=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        in1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        fftw_plan p;
        p = fftw_plan_dft_1d(N1, in, out,
                         FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan Pa;
        Pa = fftw_plan_dft_1d(N1, INa, OUTa,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan Pb;
        Pb = fftw_plan_dft_1d(N1, INb, OUTb,
                         FFTW_FORWARD, FFTW_ESTIMATE);                         
        fftw_plan p1;
        p1 = fftw_plan_dft_1d(N1, in1, out1,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan P1a;
        P1a = fftw_plan_dft_1d(N1, IN1a, OUT1a,
                         FFTW_BACKWARD, FFTW_ESTIMATE);       
        fftw_plan P1b;
        P1b = fftw_plan_dft_1d(N1, IN1b, OUT1b,
                         FFTW_BACKWARD, FFTW_ESTIMATE);                                                  
std::clock_t c_start1 = std::clock();
auto t_start1 = std::chrono::high_resolution_clock::now();        
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t ttt=0; ttt<(tt+1); ttt++)
  {

        for(int64_t n=0;n<N1;n++)
        {
        in[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in[n][1] 	= 0.;
        in1[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in1[n][1] 	= 0.;
        }
       
       
        fftw_execute(p);  
        
       
        fftw_execute(p1);  
        

        
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      
        for(int64_t z2=0; z2<N1; z2++)
        {
         
        for(int64_t q3=0;q3<N1;q3++)
        {
        if((qperp!=(Nt/2))||(q3!=(N1/2)))
        {
        INa[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INa[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        }
        else
        {
        INa[q3][0] 	= 0.;
        INa[q3][1] 	= 0.;
        INb[q3][0] 	= 0.;
        INb[q3][1] 	= 0.;
        IN1a[q3][0] 	= 0.;
        IN1a[q3][1] 	= 0.;
        IN1b[q3][0] 	= 0.;
        IN1b[q3][1] 	= 0.;        
        }
        }
       
      
        fftw_execute(Pa);  
        
     
        fftw_execute(Pb); 
        
        fftw_execute(P1a); 
        
        fftw_execute(P1b);              
       
        if((ttt>0)&&(ttt<tt))
        {
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));        
        
        }
        else
        {
        
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));
        }        
       
        }
      
    }

  }
}
        std::clock_t c_end1 = std::clock();
        auto t_end1 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used second memint: "
              << 1000.0 * (c_end1 - c_start1) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed second memint: "
              << std::chrono::duration<double, std::milli>(t_end1-t_start1).count()
              << " ms\n";
        }  

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    
    fftw_destroy_plan(p1);
    fftw_free(in1);
    fftw_free(out1);
    fftw_cleanup();
    
    
    fftw_destroy_plan(Pa);
    fftw_free(INa);
    fftw_free(OUTa);
    fftw_cleanup();
    
   
    fftw_destroy_plan(Pb);
    fftw_free(INb);
    fftw_free(OUTb);
    fftw_cleanup();    
    
    fftw_destroy_plan(P1a);
    fftw_free(IN1a);
    fftw_free(OUT1a);
    fftw_cleanup();
    
   
    fftw_destroy_plan(P1b);
    fftw_free(IN1b);
    fftw_free(OUT1b);
    fftw_cleanup();     */ 


}

void integral1(int64_t t,int64_t tt, FandIni &F, GetMSigma &getmsigma)
{

for(int64_t n=0; n<N1*N1*Nt; n++)
{
SigmarhoF31[n]=0.;
SigmarhoF32[n]=0.;
SigmarhoF33[n]=0.;
SigmaFP31[n]=0.;
SigmaFP32[n]=0.;
SigmaFP33[n]=0.;
SigmarhoP31[n]=0.;
SigmarhoP32[n]=0.;
SigmarhoP33[n]=0.;
}

if(tt<(time_spets_storage-1))
{
//std::clock_t c_start2 = std::clock();
//auto t_start2 = std::chrono::high_resolution_clock::now();       
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=0; ttt<(tt+1); ttt++)
        {
        if((ttt>0)&&(ttt<tt))
        {
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft13_access(z1, zpp, qperp, ttt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft23_access(z1, zpp, qperp, ttt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft33_access(z1, zpp, qperp, ttt));
        }
        else
        {
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft13_access(z1, zpp, qperp, ttt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft23_access(z1, zpp, qperp, ttt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft33_access(z1, zpp, qperp, ttt));        
        
        }
        
        
        }
      }
     
    }
  }
}
        /*std::clock_t c_end2 = std::clock();
        auto t_end2 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end2 - c_start2) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end2-t_start2).count()
              << " ms\n";
        }  */
}
else
{
//std::clock_t c_start2 = std::clock();
//auto t_start2 = std::chrono::high_resolution_clock::now();   
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=(tt-(time_spets_storage-2)); ttt<(tt+1); ttt++)
        {
        if((ttt>(tt-(time_spets_storage-2)))&&(ttt<tt))
        {
      
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft13_access(z1, zpp, qperp, ttt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft23_access(z1, zpp, qperp, ttt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft33_access(z1, zpp, qperp, ttt));
       
        }
        else
        {
       
        SigmarhoF31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft13_access(z1, zpp, qperp, ttt));
        SigmarhoF32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft23_access(z1, zpp, qperp, ttt));
        SigmarhoF33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.Ft33_access(z1, zpp, qperp, ttt));        
        
       
        }
        
        
        
        
        }
      }
   
    }
  }
}
       /* std::clock_t c_end2 = std::clock();
        auto t_end2 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end2 - c_start2) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end2-t_start2).count()
              << " ms\n";
        }  */
}

if(t<(time_spets_storage-1))
{
//std::clock_t c_start3 = std::clock();
//auto t_start3 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=0; ttt<(t+1); ttt++)
        {
        if((ttt>0)&&(ttt<t))
        {
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmaFP32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));
        
       
        }
        else
        {
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmaFP32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));        
        
        }
        
        
        }
      }
     
    }
  }
}
       /* std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        }  */

}
else
{
//std::clock_t c_start3 = std::clock();
//auto t_start3 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=(t-(time_spets_storage-2)); ttt<(t+1); ttt++)
        {
        if((ttt>(t-(time_spets_storage-2)))&&(ttt<t))
        {
      
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmaFP32[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=-(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));
         
        }
        else
        {
     
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmaFP32[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=-1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaF(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));        
        
        }
        
        
        
        
        }
      }
   
    }
  }
}
        /*std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        } */ 
}

   
    

        /*fftw_complex *in, *out;
        fftw_complex *in1, *out1;
        fftw_complex *INa, *OUTa;
        fftw_complex *IN1a, *OUT1a; 
        fftw_complex *INb, *OUTb;
        fftw_complex *IN1b, *OUT1b;        
        in=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        in1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        fftw_plan p;
        p = fftw_plan_dft_1d(N1, in, out,
                         FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan Pa;
        Pa = fftw_plan_dft_1d(N1, INa, OUTa,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan Pb;
        Pb = fftw_plan_dft_1d(N1, INb, OUTb,
                         FFTW_FORWARD, FFTW_ESTIMATE);                         
        fftw_plan p1;
        p1 = fftw_plan_dft_1d(N1, in1, out1,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan P1a;
        P1a = fftw_plan_dft_1d(N1, IN1a, OUT1a,
                         FFTW_BACKWARD, FFTW_ESTIMATE);       
        fftw_plan P1b;
        P1b = fftw_plan_dft_1d(N1, IN1b, OUT1b,
                         FFTW_BACKWARD, FFTW_ESTIMATE);                                                  

std::clock_t c_start3 = std::clock();
auto t_start3 = std::chrono::high_resolution_clock::now();            
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t ttt=0; ttt<(tt+1); ttt++)
  {

        for(int64_t n=0;n<N1;n++)
        {
        in[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in[n][1] 	= 0.;
        in1[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in1[n][1] 	= 0.;
        }
       
       
        fftw_execute(p);  
       
        fftw_execute(p1);  
        

        
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      
        for(int64_t z2=0; z2<N1; z2++)
        {
         
        for(int64_t q3=0;q3<N1;q3++)
        {
        if((qperp!=(Nt/2))||(q3!=(N1/2)))
        {
        INa[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INa[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        }
        else
        {
        INa[q3][0] 	= 0.;
        INa[q3][1] 	= 0.;
        INb[q3][0] 	= 0.;
        INb[q3][1] 	= 0.;
        IN1a[q3][0] 	= 0.;
        IN1a[q3][1] 	= 0.;
        IN1b[q3][0] 	= 0.;
        IN1b[q3][1] 	= 0.;        
        }
        }
      
        fftw_execute(Pa);  
        
      
        fftw_execute(Pb); 
        
        fftw_execute(P1a); 
        
        fftw_execute(P1b);              
       
        if((ttt>0)&&(ttt<tt))
        {
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));        
        
        }
        else
        {
        
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));
        }        
       
        }
      
    }

  }
}

        std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used second memint: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed second memint: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        }  
        
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    
    fftw_destroy_plan(p1);
    fftw_free(in1);
    fftw_free(out1);
    fftw_cleanup();
    
    
    fftw_destroy_plan(Pa);
    fftw_free(INa);
    fftw_free(OUTa);
    fftw_cleanup();
    
   
    fftw_destroy_plan(Pb);
    fftw_free(INb);
    fftw_free(OUTb);
    fftw_cleanup();    
    
    fftw_destroy_plan(P1a);
    fftw_free(IN1a);
    fftw_free(OUT1a);
    fftw_cleanup();
    
   
    fftw_destroy_plan(P1b);
    fftw_free(IN1b);
    fftw_free(OUT1b);
    fftw_cleanup();   */   

  
}
void integral1_rev(int64_t t,int64_t tt, FandIni &F, GetMSigma &getmsigma)
{

for(int64_t n=0; n<N1*N1*Nt; n++)
{
SigmarhoF31[n]=0.;
SigmarhoF32[n]=0.;
SigmarhoF33[n]=0.;
SigmaFP31[n]=0.;
SigmaFP32[n]=0.;
SigmaFP33[n]=0.;
SigmarhoP31[n]=0.;
SigmarhoP32[n]=0.;
SigmarhoP33[n]=0.;
}


//std::clock_t c_start0 = std::clock();
//auto t_start0 = std::chrono::high_resolution_clock::now(); 
#pragma omp parallel for 
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t z2=0; z2<N1; z2++)
  {
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      for(int64_t zpp=0; zpp<N1; zpp++)
      {
        for(int64_t ttt=tt; ttt<(t+1); ttt++)
        {
        if((ttt>tt)&&(ttt<t))
        {
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));
        
       
        }
        else
        {
        SigmarhoP31[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot13_access(z1, zpp, qperp, ttt));
        SigmarhoP32[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot23_access(z1, zpp, qperp, ttt));
        SigmarhoP33[z1+ z2*N1+qperp*N1*N1]+=1./2.*(a*dt * double(n_out_selfenergies)/(a_t*a_t))*(getmsigma.SigmaRho(z2, zpp, ttt))*(F.rhot33_access(z1, zpp, qperp, ttt));        
        
        }
        
        
        }
      }
     
    }
  }
}
        /*std::clock_t c_end0 = std::clock();
        auto t_end0 = std::chrono::high_resolution_clock::now();
        if((t==198)&&(tt==1))
        {
        std::cout << std::fixed << "CPU time used first memint: "
              << 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed first memint: "
              << std::chrono::duration<double, std::milli>(t_end0-t_start0).count()
              << " ms\n";
        } */ 

        /*fftw_complex *in, *out;
        fftw_complex *in1, *out1;
        fftw_complex *INa, *OUTa;
        fftw_complex *IN1a, *OUT1a; 
        fftw_complex *INb, *OUTb;
        fftw_complex *IN1b, *OUT1b;        
        in=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        in1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        out1=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTa=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1a=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        INb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUTb=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        IN1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        OUT1b=(fftw_complex*)fftw_malloc(N1*sizeof(fftw_complex));
        fftw_plan p;
        p = fftw_plan_dft_1d(N1, in, out,
                         FFTW_BACKWARD, FFTW_ESTIMATE);
        fftw_plan Pa;
        Pa = fftw_plan_dft_1d(N1, INa, OUTa,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan Pb;
        Pb = fftw_plan_dft_1d(N1, INb, OUTb,
                         FFTW_FORWARD, FFTW_ESTIMATE);                         
        fftw_plan p1;
        p1 = fftw_plan_dft_1d(N1, in1, out1,
                         FFTW_FORWARD, FFTW_ESTIMATE);
        fftw_plan P1a;
        P1a = fftw_plan_dft_1d(N1, IN1a, OUT1a,
                         FFTW_BACKWARD, FFTW_ESTIMATE);       
        fftw_plan P1b;
        P1b = fftw_plan_dft_1d(N1, IN1b, OUT1b,
                         FFTW_BACKWARD, FFTW_ESTIMATE);                                                  

std::clock_t c_start3 = std::clock();
auto t_start3 = std::chrono::high_resolution_clock::now();            
for(int64_t z1=0; z1<N1; z1++)
{
  for(int64_t ttt=0; ttt<(tt+1); ttt++)
  {

        for(int64_t n=0;n<N1;n++)
        {
        in[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in[n][1] 	= 0.;
        in1[n][0] 	= pow(-1.0,double(n))*(getmsigma.SigmaF(z1, n, ttt));
        in1[n][1] 	= 0.;
        }
       
       
        fftw_execute(p);  
       
        fftw_execute(p1);  
        

        
    for(int64_t qperp=0; qperp<Nt; qperp++)
    {
      
        for(int64_t z2=0; z2<N1; z2++)
        {
         
        for(int64_t q3=0;q3<N1;q3++)
        {
        if((qperp!=(Nt/2))||(q3!=(N1/2)))
        {
        INa[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INa[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(qperp))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][0] 	= real((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        INb[q3][1] 	= imag((out[q3][0]+M_I*out[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(-M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1a[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(-(F.tildep3(q3))*conj(F.tildept(F.siteqperp(Nt-qperp)))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][0] 	= real((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        IN1b[q3][1] 	= imag((out1[q3][0]+M_I*out1[q3][1])*((1./(F.abstildep(qperp,q3)))*(1.-(F.tildep3(q3))*conj(F.tildep3(q3))/(F.abstildep(qperp, q3)*F.abstildep(qperp, q3))))*M_I*exp(M_I*F.omega(qperp,q3)*double(ttt-tt)*dt * double(n_out_selfenergies)));
        }
        else
        {
        INa[q3][0] 	= 0.;
        INa[q3][1] 	= 0.;
        INb[q3][0] 	= 0.;
        INb[q3][1] 	= 0.;
        IN1a[q3][0] 	= 0.;
        IN1a[q3][1] 	= 0.;
        IN1b[q3][0] 	= 0.;
        IN1b[q3][1] 	= 0.;        
        }
        }
      
        fftw_execute(Pa);  
        
      
        fftw_execute(Pb); 
        
        fftw_execute(P1a); 
        
        fftw_execute(P1b);              
       
        if((ttt>0)&&(ttt<tt))
        {
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));        
        
        }
        else
        {
        
        SigmaFP31[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTa[z2][0]+M_I*OUTa[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1a[z2][0]+M_I*OUT1a[z2][1])*pow(-1.0,double(z2));
        SigmaFP33[z1+ z2*N1+qperp*N1*N1]+=1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUTb[z2][0]+M_I*OUTb[z2][1])*pow(-1.0,double(z2))-1./2.*(dt * double(n_out_selfenergies)*double(Nt*Nt)/(2.*double(N1)))*(OUT1b[z2][0]+M_I*OUT1b[z2][1])*pow(-1.0,double(z2));
        }        
       
        }
      
    }

  }
}

        std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t+1==tt)
        {
        std::cout << std::fixed << "CPU time used second memint: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed second memint: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        }  
        
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    fftw_cleanup();
    
    fftw_destroy_plan(p1);
    fftw_free(in1);
    fftw_free(out1);
    fftw_cleanup();
    
    
    fftw_destroy_plan(Pa);
    fftw_free(INa);
    fftw_free(OUTa);
    fftw_cleanup();
    
   
    fftw_destroy_plan(Pb);
    fftw_free(INb);
    fftw_free(OUTb);
    fftw_cleanup();    
    
    fftw_destroy_plan(P1a);
    fftw_free(IN1a);
    fftw_free(OUT1a);
    fftw_cleanup();
    
   
    fftw_destroy_plan(P1b);
    fftw_free(IN1b);
    fftw_free(OUT1b);
    fftw_cleanup();   */   

  
}
complex<double> SigmarhoF31_access(int z1, int z2, int qperp)
{
return SigmarhoF31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmarhoF32_access(int z1, int z2, int qperp)
{
return SigmarhoF32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmarhoF33_access(int z1, int z2, int qperp)
{
return SigmarhoF33[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmaFP31_access(int z1, int z2, int qperp)
{
return SigmaFP31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmaFP32_access(int z1, int z2, int qperp)
{
return SigmaFP32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmaFP33_access(int z1, int z2, int qperp)
{
return SigmaFP33[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmarhoP31_access(int z1, int z2, int qperp)
{
return SigmarhoP31[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmarhoP32_access(int z1, int z2, int qperp)
{
return SigmarhoP32[z1+ z2*N1+qperp*N1*N1];
} 

complex<double> SigmarhoP33_access(int z1, int z2, int qperp)
{
return SigmarhoP33[z1+ z2*N1+qperp*N1*N1];
} 
private:



};

#endif // MEMORYINT_H

