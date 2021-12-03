#ifndef OUTPUT_H
#define OUTPUT_H

#include <complex>
#include <iostream>
#include "observables.h"
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;


// -------------------------------------------------------------------------------------
class Output
{
  public:
    char *new_dir(int run)
    {
      char *basefolder;
      basefolder = new char[1000];
      char subfolder[1000], print[1000], currentfolder[1000], codefolder[1000];
      char to_execute[1000];
      int buffer;
      std::string err = getcwd(currentfolder, sizeof(currentfolder));
      sprintf(basefolder, "%s/output_%d", currentfolder, run);	
      sprintf(to_execute, "rm -r %s", basefolder);
      buffer = system(to_execute);
      mkdir(basefolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      
      sprintf(print, "out.new_dir(): Created output folder %s...\n", basefolder);
      cout << print << endl;
      
      sprintf(subfolder, "%s/code_copy", basefolder);
      mkdir(subfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      sprintf(codefolder, "%s/include", currentfolder);
      sprintf(to_execute, "cp -R %s %s/", codefolder, subfolder);
      buffer = system(to_execute);
      sprintf(to_execute, "cp %s/main.cpp %s/", currentfolder, subfolder);
      buffer = system(to_execute);
      sprintf(to_execute, "cp %s/main.exe %s/", currentfolder, subfolder);
      buffer = system(to_execute);
      	
      sprintf(subfolder, "%s/photon_occ", basefolder);
      mkdir(subfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      
      sprintf(subfolder, "%s/Fij", basefolder);
      mkdir(subfolder, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  
      return basefolder;
    }
    char *get_photon_occ_filename(char *basefolder)
    {
      char *print;
      print = new char[1000];
     
      sprintf(print, "%s/photon_occ/test.dat", basefolder);
      
      return print;
    }
    
    char *get_Fij_filename(char *basefolder)
    {
      char *print;
      print = new char[1000];
     
      sprintf(print, "%s/Fij/test.dat", basefolder);
      
      return print;
    }
    
    void photon_occ_to_file(Observables &obs, int64_t t, ofstream &out_photon_occ)
    {
      char print[10000];
      for (int64_t qperp=0; qperp<Nt; qperp++)
      {
	for (int64_t qz=0; qz<N1; qz++)
	{
	  sprintf(print, "%e\t%.10e\t%.10e\t%.10e\t%.10e\n", double(t-1)*dt * double(n_out_selfenergies)*m, double(qperp), double(qz), obs.dNoverVd3p_access(qz, qperp).real(), obs.dNoverVd3p_access(qz, qperp).imag());
	  out_photon_occ << print;
	}
      }
    }

    void Fij_to_file(FandIni &F, ofstream &out_Fij)
    {
      char print[10000];
	for (int64_t t=0; t<300; t++)
	{
	  sprintf(print, "%.10e\n", F.Ft11[0+ 0*N1+0*N1*N1+t*Nt*N1*N1].real());
	 
	  out_Fij << print;
	}
    }
       
    
    void append_to_files(Observables &obs, char *basefolder, int64_t t)
    {
      
      ofstream out_photon_occ;
      out_photon_occ.open(get_photon_occ_filename(basefolder), std::ofstream::out | std::ofstream::app);
      photon_occ_to_file(obs, t, out_photon_occ);
      out_photon_occ.close();
      
    
    }
    
    void append_to_filesF(FandIni &F, char *basefolder)
    {
      
    
      ofstream out_Fij;
      out_Fij.open(get_Fij_filename(basefolder), std::ofstream::out | std::ofstream::app);
      Fij_to_file(F, out_Fij);
      out_Fij.close();
    }
    
    private:
    
    
    
    
};

#endif // OUTPUT_H   
    
