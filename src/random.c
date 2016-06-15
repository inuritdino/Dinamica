/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012,2013,2014 Elias Potapov. */
/* Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
   2008, 2009, 2010, 2011 The GSL Team. */
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
   Tampere University of Technology, Department of Signal Processing
   Moscow, Russia / Tampere, Finland */
/****************************************************************************************/

#include "init.h"
#include "random.h"
#include "errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#ifdef CLOCKS_PER_SEC /* if the time.h macro is defined */
double CLOCKS_IN_SEC = CLOCKS_PER_SEC / 10; /* We use deci-seconds,i.e. 0.1 sec */
#else
double CLOCKS_IN_SEC = 100000; /* Use 10^5, since usually there are 10^6 clocks in sec */
#endif


regStat *rnd_init_cond_run()
{
  int i,j,k,m;
  FILE *out,*in;
  double **frame;
  double *rand;
  regReport *report = (regReport *)malloc(num_poin*sizeof(regReport));
  regStat *stat;
  //    regReport *one;
  /* Deal with the graphics flag */
  m = graph_flag;/* Keep the graph_flag value */
  if(m)
    graph_flag = 0;
  /* Initiate the generator */
  if(seed_flag==2)
    rng_seed = generate_seed();
  /* Sample points */
  rand = rgen_uni(num_poin*DIM,rng_seed);
  /* Open the output file for writing the intial conditions */
  out = fopen(rnd_throw_fname,"w");
  /* We use DIM-dimensional sphere for random throwing */
  for(i=0;i<num_poin;i++){
    fprintf(stdout,"\t***** %d *****\n",i+1);
    /* Print vars and pars into the header of the output.
       Only for the 1st iteration.*/
    if(!i){
      fprintf(out,"#");
      for(j=0;j<DIM;j++)
	fprintf(out,"U(%d)=%s,",j+1,var_name[j]);
      fprintf(out,"\n#");
      for(j=0;j<PARDIM;j++)
	fprintf(out,"%s=%G,",par_name[j],mu.P[j]);
      fprintf(out,"\n\n");
    }
    /* Generate and print the initials to the output */
    fprintf(out,"#%d) ",i+1);
    for(j=0;j<DIM;j++){
      /* The sphere for throwing, accounting for the Rnd_Rad constraints */
      if(Rnd_Rad_Bnd_Idx[j])
	xin[j] = Rnd_Rad_Bnd[j] +
	  (Rnd_Rad_Bnd[j+DIM] - Rnd_Rad_Bnd[j])*rand[i*DIM+j];
      else
	xin[j] = Rnd_Low + 2*Rnd_Rad*rand[i*DIM+j];
      /* Take ceiling of the xin */
      yin[j] = (long int)ceil(xin[j]);
      /* Make the initials current state */
      x[j] = xin[j];
      y[j] = yin[j];
      /* Print the initials to the output */
      fprintf(out,"%G ",xin[j]);
    }
    fprintf(out,"\n");
    /* First, run transiently */
    run(method,mynum.trans_time,0,1,1);
    /* Then, seriously */
    run(method,mynum.total_time,1,0,1);
    /* Check the dynamics */
    if(dyn_check_flag){
      report = get_sys_dynamics(report,i);
      if(report == NULL){
	fprintf(stderr,"Error occurred. Exit.\n");
	return NULL;
      }
      fprintf(out,"# %s\n",report_translate(report[i]));
      if(only_init_flag)
	fprintf(out,"\n");
    }
    if(!only_init_flag){
      /* Load the frame of just generated time series: data_name - global */
      in = fopen(data_name,"r");
      frame = load_frame(in,frame,&k);
      fclose(in);
      /* Write the frame right after the init conditions in the output */
      write_frame(out,frame,k);
      free_frame(frame);
    }
  }
  if(dyn_check_flag){
    fprintf(stdout,"----------------\n");
    fprintf(stdout,"Dynamics report:\n");
    for(i=0;i<num_poin;i++)
      fprintf(stdout,"#%d) %s\n",i+1,report_translate(report[i]));
    fprintf(stdout,"Abs.tol = %G, rel.tol = %G\n",eps_abs_am,eps_rel_am);
    /* Get statistics of found regimes */
    stat = regime_stat(report,num_poin,NULL);
    //free_regime_stat(stat);
  }

  if(m)
    graph_flag = 1;
  /* Close the stream */
  fclose(out);
  /* Draw the results */
  if(graph_flag && !only_init_flag){
    /* Do not plot more than 10 frames */
    if(num_poin > 10)
      gplot_frames(rnd_throw_fname,10,0,graph.yInd[0]);
    else
      gplot_frames(rnd_throw_fname,num_poin,0,graph.yInd[0]);
  }
  return stat;
}

gsl_rng * init_rng(const long int seed_num, const gsl_rng_type *TYPE, gsl_rng *RNG)
{
  /* Put the SEED into the environment */
  memset(rng_seed_env_val,'\0',SEED_SIZE);
  sprintf(rng_seed_env_val,"%ld",seed_num);
  setenv(rng_seed_env_name,rng_seed_env_val,1);

  /* Put the RNG type to the envronment */
  setenv(rng_type_env_name,rng_type_env_val,1);
    
  gsl_rng_env_setup();

  TYPE = gsl_rng_default;
  RNG = gsl_rng_alloc(TYPE);
  if(RNG == NULL)
    fprintf(stderr,"Could not allocate RNG.\n");
    
  return RNG;
}

double * rgen_uni(const int num_points, const long int seed_num)
{
  /* This function generates the uniformly distributed random number in (0,1) */
  int i;
  double *rand;
  /* Init the generator: rng_type and rng are the globals */
  rng = init_rng(seed_num,rng_type,rng);
  
  /*Reallocating memory...*/
  rand = (double *)malloc(num_points*sizeof(double));
  if(rand == NULL){
    fprintf(stderr,"Error: cannot allocate memory (rgen_uni).\n");
    return NULL;
  }
  /* Finally, call to the GSL routing for getting uniform random numbers */
  for(i=0;i<num_points;i++)
    rand[i] = gsl_rng_uniform_pos(rng);

  gsl_rng_free(rng);
  
  return rand;
}

int random_dist(const int num_points, const long int seed_num)
{
  int i,j;
  printf("DONT USE THIS FUNCTION!");
  if(distribution_type==NULL){
    fprintf(stderr,"You should provide valid string for the generator.\n");
    return 10;
  }

  if(num_points>DIM*NUM_POI)
    u = (double *)realloc(u,num_points*sizeof(double));

  for(i=0;i<num_points;i++)
      u[i] = gsl_ran_ugaussian(rng);/*Default behaviour.*/
  
  return 0;
}

unsigned long int generate_seed()
{
  int static gen_period = 0;
  unsigned long int seed = (int)clock();

  if(seed > CLOCKS_IN_SEC)
    seed = (unsigned long int)(seed/CLOCKS_IN_SEC);
  
  return seed;
}

int random_init()
{
  rng_seed_env_name = "GSL_RNG_SEED";
  rng_type_env_name = "GSL_RNG_TYPE";
  rng_seed_env_val = (char *)malloc(SEED_SIZE*sizeof(char));
  rng_type_env_val = "ranlxd2";/*Default random generator.*/
  num_poin = 10;
  seed_flag = 2;/*Default: automatically computed seed.*/
  dyn_check_flag = 0;
  only_init_flag = 1;

  rnd_throw_fname = strdup("rand.throw");

  int j;

  u = (double *)malloc(DIM*NUM_POI*sizeof(double));
  /*Distribution*/
  distribution_type = "gauss";/*Default*/
  distribution_n_par = 2;
  distribution_par = (double *)realloc(distribution_par,distribution_n_par*sizeof(double));
  distribution_par[0] = 0;/*Mean*/
  distribution_par[1] = 1;/*Std. dev.*/
  /******/
  /* MAX_RND = (double *)malloc(DIM*sizeof(double)); */
  /* MIN_RND = (double *)malloc(DIM*sizeof(double)); */

  /* for(j=0;j<DIM;j++){ */
  /*   MAX_RND[j] = 100.0; */
  /*   MIN_RND[j] = 0.0; */
  /* } */
  /* Sphere radius for random throwing */
  Rnd_Rad = 10;
  Rnd_Low = 0;
  Rnd_Rad_Bnd = (double *)malloc(2*DIM*sizeof(double));
  Rnd_Rad_Bnd_Idx = (int *)malloc(DIM*sizeof(int));
  for(j=0;j<DIM;j++){
    Rnd_Rad_Bnd[j] = Rnd_Low;
    Rnd_Rad_Bnd[j+DIM] = Rnd_Rad;
    Rnd_Rad_Bnd_Idx[j] = 0;
  }
  
  return 0;
}

void random_free()
{/* Freeing everything related to Random menu when closing DINAMICA */
  free(rng_seed_env_val);
  free(rnd_throw_fname);
  free(u);
  free(distribution_par);
  free(Rnd_Rad_Bnd);
  free(Rnd_Rad_Bnd_Idx);
}
