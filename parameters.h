#ifndef PARAMS_H
#define PARAMS_H

// ----------------------------------------------------------------
const complex<double> M_I = sqrt(complex<double>(-1.));			// complex number i


#define ORDER 1								// set 0 for LO, 
									// set 1 for NLO, 
									// set 2 for NNLO
									
#define INIT 0								// set 0 for zero initial photon number, 
									// set 1 for box initial condition for photon, 
																		
									  
// ----------- Parameters with initial values ----------------------------------------

//  int			OMP_NUM_THREADS = 8;							// #threads to use.

  const double 		e = 0.1;							// electric charge
  
  const double		a = 1.;									// lattice spacing numerical value
  
  const double 		m_times_a = 0.1;							// fermion mass [1/a]
  
  const double 		m = m_times_a/a;							// fermion mass numerical value
  
  const double		a_t = 1.;									// lattice spacing numerical value in transversal directions

  
  const int64_t 	N1 = 64;								// longitudinal grid size
  
  const int64_t 	Nt = 16;								// grid size in the transversal directions
  
  const double 		r = 1.0;								// Wilson parameter
  
  const double 		dt = a/20.;								// time step size
  
  double 		tf = 50. * a;								// #time-steps
  
  double 		Q = 0.1/a;                 //momentum scale in box initial condition
  
  const int 		n_out_selfenergies = 1;

  const int 		n_times_selfenergies = 1000;
  
  const int 		tmax = int(floor(tf / (dt * double(n_out_selfenergies))));
  
  
  
 
  
  const int 		time_spets_storage = 300;
  

  
#endif //PARAMS_H
