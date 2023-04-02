#ifndef _COSMOLOGY_H_
#define _COSMOLOGY_H_

#include <math.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include <fstream> 
#include <gsl/gsl_integration.h>

using namespace std;

#define FRACT_FLOAT_ERR (float) (1e-3)

/* FLAGS */
#define window_function_flag (int) 0       // 0 for top hat, 1 for sharp k and 2 for gaussian
#define Tk_flag (int) 1                   // 0 for simple Tk and 1 for the complex Tk

/* COSMOLOGICAL PARAMETERS */

//#define OM0 (double) (1)
#define N_nu (int) 3.05
#define hlittle (double) (0.6774)//(0.6774)   //(0.6727)
#define OMnu (double) (0.01) // (0.0014)   		// Planck 15 value
#define OMb (double) (0.0486) //(0.0486) //(0.02237/hlittle/hlittle)
#define OMc (double) (0.2589) // (0.1188/hlittle/hlittle) //(0.1200/hlittle/hlittle)   // (OMm-OMb)
#define OMm (double) (0.3075-OMnu) //(0.307) //(0.3166) // (0.32)
//#define z_eq (double) (3400)
#define OMl (double) (1-OMm) // (0.691) //(0.6910)// (0.6847)
#define OMtot (double) (1)  //(OMc+OMl+OMb+OMnu)
#define SIGMA8 (double) (0.8159)// (0.82)
#define POWER_INDEX (double) (0.9667)
#define delta_c0 (double) (1.686)
//#define rho0 (double) (1.2735344*1e11)   // Msun/Mpc^3 
#define rho0 (double) (8.534251287809296*1e10)   // Msun/Mpc^3 
#define T_CMB (double) (2.7255)

//
#define Theta_T (double) (T_CMB/2.7)
#define f_nu (double) (OMnu/OMm)
#define f_cb (double) ((OMc+OMb)/OMm)
#define f_b (double) (OMb/OMm)
#define f_c (double) (OMc/OMm)
#define f_nub (double) ((OMnu+OMb)/OMm)
 
// Normalization constant 2014021896 

#define Gamma_bardeen (double) (OMm * hlittle * exp(-OMb-(OMb/OMm)))   //Eq 5 (from Liddle et al. 1996) Bardeen et al. (1986)
//#define Gamma_bardeen (double) (OMm * hlittle * exp(-OMb-(OMb*sqrt(OMb*hlittle)/OMm)))   //  (for  keyword : "use_sugiyama_baryons" in HMF  

#endif
