/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2013 Elias Potapov. */
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

#include "init.h"
#include "singularity.h"
#include "random.h"
#include <math.h>
#include <string.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>


int func_multiroot(const gsl_vector *var, void *params, \
		   gsl_vector *fun)
{
  int i,j,k;

  /* Y - is another instance of X in multiroot system */
  double y[DIM];
  /* Setting all Ys to values of vector var, that must be assigned with */
/* 		   initial points for searching */
  for(i=0;i<DIM;i++)
    y[i] = gsl_vector_get (var,i);
  /* Retrieve RHS values from that Ys */
  func_odeiv (t, y, f, &mu);
  /* More generally, constant array of V is equal to F */
  double v[DIM];
  /* Assigning F to constant V */
  for(i=0;i<DIM;i++)
    v[i] = f[i];
  /* Assigning RHS to gsl_vector, following GSL Team instructions */
  for(i=0;i<DIM;i++)
    gsl_vector_set(fun,i,v[i]);
  
  return GSL_SUCCESS;
}

int multiroot_find(void)
{
  int i,j,k,l;
  double tmp[DIM];
  /* Algorithm for finding roots should be able to be changed whenever */
  const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
  gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc (T,DIM);

  /* Initializing various working entities */
  int status;/* Status of interation process */
  size_t iter = 0; /* Number of iterations */
  int sol_count = 0; /* Number of solutions */
  double x_sol[MAX_N_SOL]; /* For solution filter. Check only for 1
  element. */
  FILE *out_file; /* File for output writing */
  double x_rnd[DIM]; /* Storing random values, for not touching real X */

  /* Initializing random generator */
  const gsl_rng_type *RT = gsl_rng_default;
  gsl_rng *r = gsl_rng_alloc(RT);
  /* Generating SEED */
  for(k=0;k<DIM;k+=2)/* Heuristic rules to update the seed */
    /*all formulae here is from my mind:)*/
    rng_seed+=(int)0.01*x[k]*rng_seed-(int)0.01*x[k+1]*rng_seed+k;

  /* The following is useless, since rng_seed is defines as UNSIGNED INT */
  /* if(rng_seed<0) */
  /*   rng_seed=-rng_seed; */
  r = init_rng(rng_seed,RT,r);/* INIT RNG */

  /* Allocating solution array */
  gsl_vector *xroot = gsl_vector_alloc(DIM);

  /* GSL multiroot function */
  gsl_multiroot_function froot = {&func_multiroot, DIM, &mu};

  /* How many points to generate as initials */
  int pts;
  if(RND_FLAG) /*  If random then PTS points */
    pts=PTS;
  else /* If no random and one initial set(current Xs), then 1 point */
    pts=1;
  for(j=0;j<pts;j++){
    if(RND_FLAG){
      for(i=0;i<DIM;i++){
	x_rnd[i]=gsl_rng_uniform_pos(r);
	x_rnd[i]=x_rnd[i]*(MAX_SNG-MIN_SNG)+MIN_SNG;
      }
    }
    else {
      for(i=0;i<DIM;i++)
	x_rnd[i]=x[i]; /* Not touching real X ever */
    }
    for(i=0;i<DIM;i++)
      gsl_vector_set(xroot,i,x_rnd[i]);/* Current X is located, not X_init */

    gsl_multiroot_fsolver_set(s,&froot,xroot);

    do
      {
	iter++;
	status = gsl_multiroot_fsolver_iterate(s);

	/* print_state(iter,s); */
	if(status)
	  break;

	status = gsl_multiroot_test_residual(s->f,TEST);
      }
    while( status == GSL_CONTINUE && iter < MAX_ITER );

    /* SORTING THE SOLUTIONS. */
    if(strcmp(gsl_strerror(status),"success") == 0)
      {
	if(sol_count==0){
	  x_sol[sol_count]=gsl_vector_get(s->x,0);
	  printf("%d) ",j+1); 
	  //printf("status = %s\n", gsl_strerror(status));
	  print_state(iter,s);
	  printf("Solution N %d\n",sol_count+1);
	  if(FILE_FLAG!=0 && ss_name!=NULL){
	    out_file=fopen(ss_name,"w");
	    if(out_file == NULL){
	      printf("Cannot open file\n");
	      exit(8);
	    }
	    fprintf(out_file,"Solution N %d\n",sol_count+1);
	    for(k=0;k<DIM;k++){
	      /* Print to the file */
	      fprintf(out_file,"U(%d)=%lf\n",k+1,gsl_vector_get(s->x,k));
	    }
	  }
	  for(k=0;k<DIM;k++){
	    /* Print to the stdout */
	    fprintf(stdout,"U(%d)=%lf\n",k+1,gsl_vector_get(s->x,k));
	  }
	  for(k=0;k<DIM;k++)
	    tmp[k]=gsl_vector_get(s->x,k);
	  if((ss_stab(tmp))==0)
	    printf("Stable\n");
	  else
	    printf("Unstable\n");
	  sol_count++;
	}
	if(sol_count > 0){
	  for(l=0;l<sol_count;l++){
	    if(fabs((double)gsl_vector_get(s->x,0)-x_sol[l]) < SOL_ERROR ){
		
	      break;
	    }
	    if(l == sol_count-1){
	      x_sol[sol_count]=gsl_vector_get(s->x,0);
	      printf("%d) ",j+1); 
	      //printf("status = %s\n", gsl_strerror(status));
	      print_state(iter,s);
	      printf("Solution N %d\n",sol_count+1);
	      if(FILE_FLAG!=0 && ss_name!=NULL){
		if(out_file == NULL){
		  printf("Cannot open file\n");
		  exit(8);
		}
		fprintf(out_file,"Solution N %d\n",sol_count+1);
		for(k=0;k<DIM;k++){
		  /* Print to the file */
		  fprintf(out_file,"U(%d)=%lf\n",k+1,gsl_vector_get(s->x,k));
		}
	      }
	      for(k=0;k<DIM;k++){
		/* Print to stdout */
		fprintf(stdout,"U(%d)=%lf\n",k+1,gsl_vector_get(s->x,k));
	      }
	      if(sol_count+1 == MAX_N_SOL){
		printf("Got max number of solutions. Exit.\n");
		return 0;
	      }
	      for(k=0;k<DIM;k++)
		tmp[k]=gsl_vector_get(s->x,k);
	      if((ss_stab(tmp))==0)
		printf("Stable\n");
	      else
		printf("Unstable\n");
	      sol_count++;
	    }
	  }
	}
      }
	  
    iter=0;
  }
  
  if(FILE_FLAG != 0)
    fclose(out_file);
  gsl_multiroot_fsolver_free(s);
  gsl_vector_free(xroot);
  if(RND_FLAG)
    gsl_rng_free(r);

  return GSL_SUCCESS;
}

int
print_state(size_t iter, gsl_multiroot_fsolver *s)
{
  int k;
  double temp[DIM];
  double dfdx[DIM*DIM];
  double dfdt[DIM];
  printf("iter = %3u  F(x) = ",(unsigned)iter);
  for(k=0;k<DIM;k++){
    printf("%.3e ",gsl_vector_get(s->f,k));
  }
  printf("\n");
  /* Print Jacobian for the steady state value */
  /* for(k=0;k<DIM;k++) */
  /*   temp[k]=gsl_vector_get(s->x,k); */
  /* if(jac_flag) */
  /*   printf("Analytical J=\n"); */
  /* else */
  /*   printf("Numerical J=\n"); */
  /* jac(t,temp,dfdx,dfdt,&mu); */
  /* jac_flag=!jac_flag; */
  /* printf("====================\n"); */
  /* if(jac_flag) */
  /*   printf("Analytical J=\n"); */
  /* else */
  /*   printf("Numerical J=\n"); */
  /* jac(t,temp,dfdx,dfdt,&mu); */
  /* jac_flag=!jac_flag; */
  
  return 0;
}

int multiroot_init()
{
  TEST=1e-7; /* Error for test of whether iteration process is successful. */
  SOL_ERROR=1e-4; /* Error to distinguish two different roots */
  RND_FLAG=1; /* Flag for use of random generation of PTS points to start from */
  PTS=500; /* Number of points to generate */
  FILE_FLAG=0; /* Flag for writing file output */
  MAX_SNG=100; /* Max bound random generator points */
  MIN_SNG=0; /* Min bound random generator points */
  MAX_ITER=1000; /* max num of iterations of each initial point provided */
  MAX_N_SOL=20;  /* max num of solutions(roots) to find */

  return 0;
}
