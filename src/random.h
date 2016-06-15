/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012,2013,2014 Elias Potapov. */
/* Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
   The GSL Team. */

/*****************************************************************************************/
/* This file is part of DINAMICA. */

/* DINAMICA is free software: you can redistribute it and/or modify */
/* it under the terms of the GNU General Public License as published by */
/* the Free Software Foundation, either version 3 of the License, or */
/* (at your option) any later version. */

/* DINAMICA is distributed in the hope that it will be useful, */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the */
/* GNU General Public License for more details. */

/* You should have received a copy of the GNU General Public License */
/* along with DINAMICA.  If not, see <http://www.gnu.org/licenses/>. */
/****************************************************************************************/
/****************************************************************************************/
/* Original author is Elias Potapov <elias.potapov@gmail.com>
   Lomonosov Moscow State University, Biophysics Dep..
   Lebedev Physical Inst., Dep. of Theoretical Physics.
   Moscow, Russia */
/****************************************************************************************/
#include <string.h>
#include <gsl/gsl_rng.h>
#define NUM_POI 1000
#define SEED_SIZE 10

double *u;/* helping variable, the output of a generator */
char *rnd_throw_fname; /* filename to place the results of random throwing */
unsigned long int rng_seed;/*SEED for the generator. We keep the global value since the seed
	      can be set manually and we need, in this case, to keep its value for
	      all sessions.*/
int short seed_flag;/*Flag for how to compute seed.*/
int dyn_check_flag;/* Flag whether to check the dynamics with T-system */
int only_init_flag;/* Flag whether to write out only initial conditions without actual data */
int num_poin;/* Number of points in the phase space to generate */
/* ************************************************************************** */
/*                           ENVIRONMENTAL STUFF                              */
/* GSL_RNG_TYPE=<type> */
/* GSL_RNG_SEED=<seed> */
const char *rng_seed_env_name;/*const string ="GSL_RNG_TYPE"*/
const char *rng_type_env_name;/*const string ="GSL_RNG_SEED"*/
char *rng_type_env_val;/* <type> name */
char *rng_seed_env_val;/* <seed> string */
/* ************************************************************************** */
/*GSL definitions*/
const gsl_rng_type * rng_type;/*See GSL reference*/
gsl_rng * rng;/*See GSL reference*/
/*Distribution*/
char * distribution_type;
int distribution_n_par;
double *distribution_par;
/**************/
/* double *MAX_RND; /\* Upper bound for random dispersion *\/ */
/* double *MIN_RND; /\* Lower bound for random dispersion *\/ */
double Rnd_Rad; /* Radius of DIM-dimensional sphere limiting the random throwing  */
double *Rnd_Rad_Bnd;/* 2xDIM long vector of min and max of random radius for throwing for
		       the variables */
int *Rnd_Rad_Bnd_Idx;/* array of logical indices showing which vars are to be limited
			with Rnd_Rad_Bnd constraints. */
double Rnd_Low; /* Lower bound for throwing */
