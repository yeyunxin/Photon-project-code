#include <boost/lambda/lambda.hpp>
#include <ios>
#include <iostream>
#include <omp.h>
#include <cstring>
#include <ctime>
#include <chrono>
#include <cmath>
#include <time.h>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <algorithm>
#include "boost/multi_array.hpp"
#include <boost/array.hpp> 
#include <cassert>
#include <complex>
#include <stdint.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_randist.h>
#include <sys/resource.h>




using namespace std;

#include "include/parameters.h"
#include "include/getMandSigma.h"
#include "include/FandIni.h"
#include "include/MemoryInt.h"
#include "include/time_evolution.h"
#include "include/EulerMethod.h"
#include "include/observables.h"
#include "include/output.h"


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

int main(int argc, char *argv[])
{
  /*using std::cout;
  using std::endl;

  double vm, rss;*/


  GetMSigma getmsigma;
  FandIni F;
  MemoryInt memint;
  MyEulerSolver<time_evolution> stepper;
  Observables obs;
  Output out;
  int Direction; //set 0 for reverse evolution, set 1 for postitive evolution.

  

  int run;
  if (argc != 2) {cout << "main() ERROR: argc != 2" << endl; return 1;}
  sscanf(argv[1], "%d", &run);
  
  char *basefolder = out.new_dir(run);
  
    getmsigma.open_file();
    F.initialize();       //get initial condition for F
 
    for(int64_t t=2; t<1000; t++)
    {
      //std::clock_t c_start = std::clock();
      //auto t_start = std::chrono::high_resolution_clock::now();
      getmsigma.readin(t-1);
      /*std::clock_t c_end = std::clock();
      auto t_end = std::chrono::high_resolution_clock::now();
      if(t==199)
      {std::cout << std::fixed << "CPU time used read in: "
              << 1000.0 * (c_end - c_start) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed read in: "
              << std::chrono::duration<double, std::milli>(t_end-t_start).count()
              << " ms\n";}*/
              
      //if(t==168){cout<<getmsigma.M_ptr[2+N1*8]<<endl; cout<<getmsigma.SigmaF(1, 3, 38)<<endl; cout<<getmsigma.SigmaRho(2, 0, 55)<<endl;}
    
   
      for(int64_t tt=0; tt<2; tt++)
      {
        //std::clock_t c_start1 = std::clock();
        //auto t_start1 = std::chrono::high_resolution_clock::now();        
        memint.integral(t-1,tt, F, getmsigma);
        /*std::clock_t c_end1 = std::clock();
        auto t_end1 = std::chrono::high_resolution_clock::now();  
        if((t==199)&&(tt==1))
        {std::cout << std::fixed << "CPU time used memint: "
              << 1000.0 * (c_end1 - c_start1) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed memint: "
              << std::chrono::duration<double, std::milli>(t_end1-t_start1).count()
              << " ms\n";}*/

        //std::clock_t c_start2 = std::clock();
        //auto t_start2 = std::chrono::high_resolution_clock::now();                            
        stepper.make_step(F ,getmsigma, memint, t-1, tt, Direction); 
        /*std::clock_t c_end2 = std::clock();
        auto t_end2 = std::chrono::high_resolution_clock::now();  
        if((t==199)&&(tt==1))
        {std::cout << std::fixed << "CPU time used make step: "
              << 1000.0 * (c_end2 - c_start2) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed make step: "
              << std::chrono::duration<double, std::milli>(t_end2-t_start2).count()
              << " ms\n";}  */      
         
      }
    }
    
    for(int64_t t=2; t<1000; t++)
    {
      F.acess_ini(t);
      if(t>3)
      {
      Direction=0;
      for(int64_t tt=t-2; tt>1; tt--)
      {
        //std::clock_t c_start3 = std::clock();
        //auto t_start3 = std::chrono::high_resolution_clock::now();
        getmsigma.readint_rev(t,tt+1);
        /*std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used read in: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed read in: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        }*/
        //if(tt==35&&t==23){cout<<getmsigma.M_ptr[3+N1*6]<<endl; cout<<getmsigma.SigmaF(0, 2, 33)<<endl; cout<<getmsigma.SigmaRho(1, 1, 27)<<endl;}

        //std::clock_t c_start4 = std::clock();
        //auto t_start4 = std::chrono::high_resolution_clock::now();      
        memint.integral1_rev(t,tt+1, F, getmsigma);
        /*std::clock_t c_end4 = std::clock();
        auto t_end4 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used memint: "
              << 1000.0 * (c_end4 - c_start4) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed memint: "
              << std::chrono::duration<double, std::milli>(t_end4-t_start4).count()
              << " ms\n";
        }*/
      
        //std::clock_t c_start5 = std::clock();
        //auto t_start5 = std::chrono::high_resolution_clock::now(); 
        stepper.make_stept_rev(F ,getmsigma, memint, t, tt+1, Direction); 
        /*std::clock_t c_end5 = std::clock();
        auto t_end5 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used make step: "
              << 1000.0 * (c_end5 - c_start5) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed make step: "
              << std::chrono::duration<double, std::milli>(t_end5-t_start5).count()
              << " ms\n";
        }*/
       
      }     
      }
      Direction=1; //positive evolution
      for(int64_t tt=2; tt<t+1; tt++)
      {
        //std::clock_t c_start3 = std::clock();
        //auto t_start3 = std::chrono::high_resolution_clock::now();
        getmsigma.readint(t,tt-1);
        /*std::clock_t c_end3 = std::clock();
        auto t_end3 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used read in: "
              << 1000.0 * (c_end3 - c_start3) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed read in: "
              << std::chrono::duration<double, std::milli>(t_end3-t_start3).count()
              << " ms\n";
        }*/
        //if(tt==35&&t==23){cout<<getmsigma.M_ptr[3+N1*6]<<endl; cout<<getmsigma.SigmaF(0, 2, 33)<<endl; cout<<getmsigma.SigmaRho(1, 1, 27)<<endl;}

        //std::clock_t c_start4 = std::clock();
        //auto t_start4 = std::chrono::high_resolution_clock::now();      
        memint.integral1(t,tt-1, F, getmsigma);
        /*std::clock_t c_end4 = std::clock();
        auto t_end4 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used memint: "
              << 1000.0 * (c_end4 - c_start4) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed memint: "
              << std::chrono::duration<double, std::milli>(t_end4-t_start4).count()
              << " ms\n";
        }*/
      
        //std::clock_t c_start5 = std::clock();
        //auto t_start5 = std::chrono::high_resolution_clock::now(); 
        stepper.make_stept(F ,getmsigma, memint, t, tt-1, Direction); 
        /*std::clock_t c_end5 = std::clock();
        auto t_end5 = std::chrono::high_resolution_clock::now();
        if(t==tt)
        {
        std::cout << std::fixed << "CPU time used make step: "
              << 1000.0 * (c_end5 - c_start5) / CLOCKS_PER_SEC << " ms\n"
              << "Wall clock time passed make step: "
              << std::chrono::duration<double, std::milli>(t_end5-t_start5).count()
              << " ms\n";
        }*/
       
      }
     F.diagupdate(t);
     obs.compute_observables(F,t); //compute photon occupation
     
   
     out.append_to_files(obs, basefolder, t); 
     

      
     //if(tt==299){out.append_to_filesF(F, basefolder); } 
    }
    getmsigma.close_file();
  return 0;
  
}
