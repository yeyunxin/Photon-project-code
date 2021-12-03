#ifndef GETMSIGMA_H
#define GETMSIGMA_H

#include <complex>
#include <iostream>
#include <fstream>

using namespace std;

class GetMSigma
{
public:
  double* M_ptr;
  double* SigmaF_ptr;
  double* SigmaRho_ptr;

  std::ifstream myfile;
  std::ifstream mySigmaFfile;
  std::ifstream mySigmaRhofile; 
  
  GetMSigma()
  {
    
    #if (ORDER==0)
    M_ptr = new double[N1];
    #elif (ORDER==1)
    M_ptr = new double[5*N1];
    #elif (ORDER==2)
    M_ptr = new double[14*N1];
    #endif
    SigmaF_ptr = new double[N1*N1*n_times_selfenergies];
    SigmaRho_ptr = new double[N1*N1*n_times_selfenergies];
  }
  
  
  ~GetMSigma()
  {
    
    
    delete [] M_ptr;
    delete [] SigmaF_ptr;
    delete [] SigmaRho_ptr;
    
  }

  
  void ini_class_variables()
  {

    #if (ORDER==0)
    memset(M_ptr, 0, sizeof(double)*N1);
    #elif (ORDER==1)
    memset(M_ptr, 0, sizeof(double)*5*N1);
    #elif (ORDER==2)
    memset(M_ptr, 0, sizeof(double)*14*N1);
    #endif
    memset(SigmaF_ptr, 0, sizeof(double)*N1*N1*n_times_selfenergies);
    memset(SigmaRho_ptr, 0, sizeof(double)*N1*N1*n_times_selfenergies);
   
    

  }
  
  char *get_m_filename(char *currentfolder)
  {  
      char *print;
      print = new char[1000];
      sprintf(print, "%s/m/sample_0.dat", currentfolder);
      
      return print;
  }
  
  char *get_sigmaF_filename(char *currentfolder)
  {  
      char *print;
      print = new char[1000];
      sprintf(print, "%s/sigmaF/sample_0.dat", currentfolder);
      
      return print;
  }
  
  char *get_sigmaRho_filename(char *currentfolder)
  {  
      char *print;
      print = new char[1000];
      sprintf(print, "%s/sigmaRho/sample_0.dat", currentfolder);
      
      return print;
  }

  void open_file()
  {
  
  char currentfolder[1000];
  std::string err = getcwd(currentfolder, sizeof(currentfolder));
  
  
  myfile.open(get_m_filename(currentfolder),std::ifstream::binary);
  mySigmaFfile.open(get_sigmaF_filename(currentfolder),std::ifstream::binary);
  mySigmaRhofile.open(get_sigmaRho_filename(currentfolder),std::ifstream::binary);
  }
  void close_file()
  {
  
  myfile.close();
  mySigmaFfile.close();
  mySigmaRhofile.close();
  }
  
  
  void readin(int64_t t)
  {
  int64_t timestep1;
 

  
  timestep1=t;
  
 
  ini_class_variables();
  
  #if (ORDER==0)
  for(int64_t i=81+timestep1*N1*117;i<81+(timestep1+1)*N1*117;i+=117)
  {
  
  fseek ( myfile , i , SEEK_SET );
  fread (&buffer,sizeof(buffer),1,myfile);
  M_ptr[nz]=e*e*stod(buffer);
  nz++;
  }
  #elif (ORDER==1)
  myfile.seekg(sizeof(double)*timestep1*5*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*5*N1);
  for(int64_t nz=0; nz<5*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  } 
  #elif (ORDER==2)
  myfile.seekg(sizeof(double)*timestep1*14*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*14*N1);
  for(int64_t nz=0; nz<14*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  }
  #endif
  
 
  mySigmaRhofile.seekg(sizeof(double)*timestep1*N1*N1*tmax, std::ios_base::beg);
  mySigmaRhofile.read((char*) &SigmaRho_ptr[0], sizeof(double)*(timestep1+1)*N1*N1);
 
  mySigmaFfile.seekg(sizeof(double)*timestep1*N1*N1*tmax, std::ios_base::beg);
  mySigmaFfile.read((char*) &SigmaF_ptr[0], sizeof(double)*(timestep1+1)*N1*N1);
  
  for(int64_t nz=0; nz<(timestep1+1)*N1*N1; nz++)
  {
  SigmaRho_ptr[nz]=e*e*SigmaRho_ptr[nz];
  SigmaF_ptr[nz]=e*e*SigmaF_ptr[nz];
  
  }
  }
  
  void readint(int64_t t, int64_t tt)
  {
  int64_t timestep1;
  int64_t timestep2;
 
  timestep1=t;
  timestep2=tt;
 
 
  ini_class_variables();
  #if (ORDER==0)
  for(int64_t i=81+timestep1*N1*117;i<81+(timestep1+1)*N1*117;i+=117)
  {
  
  fseek ( myfile , i , SEEK_SET );
  fread (&buffer,sizeof(buffer),1,myfile);
  M_ptr[nz]=e*e*stod(buffer);
  nz++;
  }
  #elif (ORDER==1)
  myfile.seekg(sizeof(double)*timestep2*5*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*5*N1);
  for(int64_t nz=0; nz<5*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  } 
  #elif (ORDER==2)
  myfile.seekg(sizeof(double)*timestep2*14*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*14*N1);
  for(int64_t nz=0; nz<14*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  }
  #endif
  
 
  mySigmaRhofile.seekg(sizeof(double)*timestep2*N1*N1*tmax, std::ios_base::beg);
  mySigmaRhofile.read((char*) &SigmaRho_ptr[0], sizeof(double)*(timestep1+1)*N1*N1);
 
  mySigmaFfile.seekg(sizeof(double)*timestep2*N1*N1*tmax, std::ios_base::beg);
  mySigmaFfile.read((char*) &SigmaF_ptr[0], sizeof(double)*(timestep1+1)*N1*N1);
  
  for(int64_t nz=0; nz<(timestep1+1)*N1*N1; nz++)
  {
  SigmaRho_ptr[nz]=e*e*SigmaRho_ptr[nz];
  SigmaF_ptr[nz]=e*e*SigmaF_ptr[nz];
  
  }
 
  }
  
  void readint_rev(int64_t t, int64_t tt)
  {
  int64_t timestep1;
  int64_t timestep2;
 
  timestep1=t;
  timestep2=tt;
 
 
  ini_class_variables();
  #if (ORDER==0)
  for(int64_t i=81+timestep1*N1*117;i<81+(timestep1+1)*N1*117;i+=117)
  {
  
  fseek ( myfile , i , SEEK_SET );
  fread (&buffer,sizeof(buffer),1,myfile);
  M_ptr[nz]=e*e*stod(buffer);
  nz++;
  }
  #elif (ORDER==1)
  myfile.seekg(sizeof(double)*timestep2*5*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*5*N1);
  for(int64_t nz=0; nz<5*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  } 
  #elif (ORDER==2)
  myfile.seekg(sizeof(double)*timestep2*14*N1, std::ios_base::beg);
  myfile.read((char*) &M_ptr[0], sizeof(double)*14*N1);
  for(int64_t nz=0; nz<14*N1; nz++)
  {
  M_ptr[nz]=e*e*M_ptr[nz];

  }
  #endif
  
 
  mySigmaRhofile.seekg(sizeof(double)*timestep2*N1*N1*tmax, std::ios_base::beg);
  mySigmaRhofile.read((char*) &SigmaRho_ptr[0], sizeof(double)*(timestep1+1)*N1*N1);
  
  for(int64_t nz=0; nz<(timestep1+1)*N1*N1; nz++)
  {
  SigmaRho_ptr[nz]=e*e*SigmaRho_ptr[nz];
  }
 
  }
double SigmaF(int z1, int z2, int64_t t)
{
return SigmaF_ptr[z2+z1*N1+t*N1*N1];
}  

double SigmaRho(int z1, int z2, int64_t t)
{
return SigmaRho_ptr[z2+z1*N1+t*N1*N1];
}  
  
private:



};



#endif // GETMSIGMA_H

