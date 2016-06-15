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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland*/
/****************************************************************************************/

#include "init.h"
#include "random.h"
#include <math.h>
#include <stdio.h>

int solver_lang(double x[],double h,
		int(*func)(double,const double [],double [],void *),
		double Nstd[])
{
 /* NOTE: Nstd[] is an array of random numbers following standard gaussian distribution*/
  int i;
  double ksi[DIM];/* Brownian coefficients */
  double dksidx[DIM*DIM];
  double dksidt[DIM];
  double f[DIM];/* Deterministic RHS's of the equations */

  if(!strcmp(method,"eu")){
    func(t,x,f,&mu);/* Calculate deterministic RHS's */
    lang_amend(t,x,ksi,mu.P);/* calculate ksi */
    for(i=0;i<DIM;i++){
      x[i] = x[i] + f[i]*h + ksi[i]*sqrt(h)*Nstd[i];
    }
  }
  else if(!strcmp(method,"milst")){
    func(t,x,f,&mu);/* Calculate deterministic RHS's at the given point*/
    lang_amend(t,x,ksi,mu.P);/* Calculate noise terms ksi */
    jac_general(t,x,dksidx,dksidt,&mu,lang_amend);
    for(i=0;i<DIM;i++){
      x[i] = x[i] + f[i]*h + ksi[i]*sqrt(h)*Nstd[i] + \
	0.5*ksi[i]*dksidx[i*DIM+i]*(h*Nstd[i]*Nstd[i] - h);
    }
  }
  else{
    fprintf(stderr,"Error: I do not know the stochastic method.\n");
    return 100;
  }
  
  return 0;
}

int solver(double x[],double h,
	   int(*func)(double,const double [],double [],void *)) 
{
  int i;
  double k1[DIM],k2[DIM],k3[DIM],k4[DIM];
  double x_prev[DIM];

  if(strcmp(method,"run-kut4")==0){
    // 1-st STEP
    for(i=0;i<DIM;i++)
      x_prev[i] = x[i];

    func(t,x,f,&mu);
    for(i=0;i<DIM;i++)
      k1[i] = f[i];

    //2-nd STEP
    for(i=0;i<DIM;i++)
      x[i] = x_prev[i] + h*k1[i]/2;

    func(t,x,f,&mu);
    for(i=0;i<DIM;i++)
      k2[i] = f[i];

    //3-rd STEP
    for(i=0;i<DIM;i++)
      x[i] = x_prev[i] + h*k2[i]/2;

    func(t,x,f,&mu);
    for(i=0;i<DIM;i++)
      k3[i] = f[i];

    //4-th STEP
    for(i=0;i<DIM;i++)
      x[i] = x_prev[i] + h*k3[i];

    func(t,x,f,&mu);
    for(i=0;i<DIM;i++)
      k4[i] = f[i];

    //Finally, calculating the next point
    for(i=0;i<DIM;i++)
      x[i] = x_prev[i] + h*(k1[i] + 2*k2[i] + 2*k3[i] + k4[i])/6;

  }
  else if (strcmp(method,"eu")==0){
    func(t,x,f,&mu);/*Compute derivatives*/
    for(i=0;i<DIM;i++)/*Get x for next moment*/
      x[i] = x[i] + f[i]*h;
  }
  else {
    fprintf(stderr,"Error: the method must be either 'eu' or 'run-kut4'.\n");
    return 1;
  }
  return 0;

}

int solver_lin(const double x_0[],double x_lin[],
	       double f_lin[], double h,int(*func)(double,const double [],
	       double [], double [], double [], void *)) 
{
  /*REMEMBER: linear system always takes the current vector x_0, so call
    this function before solver() for nonlinear system which changes
    x_0. x_0 is an attractor values.*/
  int i,j;
  double k1[DIM],k2[DIM],k3[DIM],k4[DIM];
  double x_prev[DIM];
  double dfdx[DIM*DIM];
  double dfdt[DIM];

  if(strcmp(method,"run-kut4")==0){
    // 1-st STEP
    i=1;
    while(i <= DIM){
      x_prev[i-1]=*(x_lin+i-1);
      i++;
    }
    jac(t,x_0,dfdx,dfdt,&mu);
    func(t,x_lin,f_lin,dfdx,dfdt,&mu);i=1;
    while (i <= DIM){
      k1[i-1]=*(f_lin+i-1);
      i++;
    }

    //2-nd STEP
    i=1;
    while(i<=DIM){
      *(x_lin+i-1)=x_prev[i-1] + h*k1[i-1]/2;
      i++;
    }
    jac(t,x_0,dfdx,dfdt,&mu);
    func(t,x_lin,f_lin,dfdx,dfdt,&mu);i=1;
    while(i<=DIM){
      k2[i-1]=*(f_lin+i-1);
      i++;
    }

    //3-rd STEP
    i=1;
    while(i<=DIM){
      *(x_lin+i-1)=x_prev[i-1] + h*k2[i-1]/2;
      i++;
    }
    jac(t,x_0,dfdx,dfdt,&mu);
    func(t,x_lin,f_lin,dfdx,dfdt,&mu);i=1;
    while(i<=DIM){
      k3[i-1]=*(f_lin+i-1);
      i++;
    }

    //4-th STEP
    i=1;
    while(i<=DIM){
      *(x_lin+i-1)=x_prev[i-1] + h*k3[i-1];
      i++;
    }
    jac(t,x_0,dfdx,dfdt,&mu);
    func(t,x_lin,f_lin,dfdx,dfdt,&mu);i=1;
    while(i<=DIM){
      k4[i-1]=*(f_lin+i-1);
      i++;
    }
    /*Finally, calculating the next point*/
    i=1;
    while(i<=DIM){
      *(x_lin+i-1)=x_prev[i-1] + h*(k1[i-1] + 2*k2[i-1] + 2*k3[i-1] + k4[i-1])/6;
      i++;
    }
  }
  else if (strcmp(method,"eu")==0){
    i=1;
    while(i<=DIM){
      x_prev[i-1]=*(x_lin+i-1);/*Store x from previous moment*/
      i++;
    }
    jac(t,x_0,dfdx,dfdt,&mu);
    func(t,x_lin,f_lin,dfdx,dfdt,&mu);/*Compute derivatives*/
    i=1;
    while(i<=DIM){/*Get x for next moment*/
      *(x_lin+i-1) = x_prev[i-1] + *(f_lin+i-1)*h;
      i++;
    }
  }
  else {
    fprintf(stderr,"Error: I dont know the method\n");
    return 1;
  }
  return 0;

}

int gill_init(long int *X, double *A, double *A_sum)
{
  /*This function initialize the gillespie algorithm. X is array of
    species numbers. A is the array of propensities of
    reactions. A_sum is the sum of A.*/
  int i;
  *A_sum=0.0;
  double Y[DIM];
  for(i=0;i<DIM;i++)/*Doublify*/
    Y[i]=(double)X[i];
  /*Compute all propensities*/
  propensity(Y,mu.P,A);
  for(i=0;i<NREAC;i++)
    //printf("A[%d]=%lf\n",i,A[i]);
    *A_sum = *A_sum + A[i];
  
  return 0;
}


