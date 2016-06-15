/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012,2013 Elias Potapov. */
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
#include "errors.h"
#include "random.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_deriv.h>

int func_odeiv(double t, const double x[], double f[],
	       void *params)
{
  int i;
  MU nu = *(MU *)params;
	
  rhs(t,x,f,nu.P);

  return GSL_SUCCESS;
}

int jac_general(double t,const double x[],double *dfdx,double *dfdt, void *params,
		int(*func)(double,const double [],double [],void *))
{
  /* This function is a copy of the function jac(...) below, except for it accepts
     the additional argument int(*func) which extends the capability of the
     jac_general to any function. While jac(...) is needed for the GSL step function
     algorithms, jac_general(...) can be used to calculate Jacobian matrix for some
     methods (e.g. Milstein) for the Langevin SDE's, where the derivatives of the
     stochastic terms must be calculated. NOTE: the calculation of the Jacobian takes
     place only numercially.*/
  MU nu = *(MU *)params; /* Get the parameter values */
  int i,j;/* looping variables */
  double h = mynum.step;/*Step for computing derivatives, the same like*/
  double temp[DIM];/* Temporary array for copying var values */
  /* for the integration. Must be changed maybe. */
  for(i=0;i<DIM;i++){
    dfdt[i]=0.0;
    for(j=0;j<DIM;j++)
      dfdx[i*DIM+j]=0.0;
  }
  /* This function computes numerically the Jacobian, in the point
     specified in x[] array, if jac_flag=0. The method is five-point
     stencil approximation for computing derivatives. The step is
     taken from mynum.step, which is the step for the
     integration~(maybe will need to be fixed).*/
  for(j=0;j<DIM;j++){/* Loop over variables */
    /* Get function at  (X+2h) point */
    for(i=0;i<DIM;i++) temp[i] = x[i];/* Copy array */
    temp[j] = x[j]+2*h;/* Change the desired var */
    func(t,temp,f,&nu);/* Compute all RHS's */
    for(i=0;i<DIM;i++) dfdx[i*DIM+j] += -f[i];
    /* Get function at (X+h) point */
    for(i=0;i<DIM;i++) temp[i] = x[i];
    temp[j] = x[j]+h;
    func(t,temp,f,&nu);
    for(i=0;i<DIM;i++) dfdx[i*DIM+j] += 8*f[i];
    /* Get function at (X-h) point */
    for(i=0;i<DIM;i++) temp[i] = x[i];
    temp[j] = x[j]-h;
    func(t,temp,f,&nu);
    for(i=0;i<DIM;i++) dfdx[i*DIM+j] += -8*f[i];
    /* Get function at (X-2h) point */
    for(i=0;i<DIM;i++) temp[i] = x[i];
    temp[j] = x[j]-2*h;
    func(t,temp,f,&nu);
    for(i=0;i<DIM;i++) dfdx[i*DIM+j] += f[i];
    /* Dividing by 12h */
    for(i=0;i<DIM;i++) dfdx[i*DIM+j] /= 12*h;
  }

  return 0;
}

int jac (double t,const double x[], double *dfdx,
	 double *dfdt, void *params)
{
  int i,j;/* looping variables */
  double h = mynum.step;/*Step for computing derivatives, the same like*/
  double temp[DIM];/* Temporary array for copying var values */
  /* for the integration. Must be changed maybe. */
  MU nu = *(MU *)params;
  for(i=0;i<DIM;i++){
    dfdt[i]=0.0;
    for(j=0;j<DIM;j++)
      dfdx[i*DIM+j]=0.0;
  }
  if(jac_flag)
    jacobian(t,x,dfdx,dfdt,nu.P);
  else{
    /* This function computes numerically the Jacobian, in the point
       specified in x[] array, if jac_flag=0. The method is five-point
       stencil approximation for computing derivatives. The step is
       taken from mynum.step, which is the step for the
       integration~(maybe will need to be fixed).*/
    for(j=0;j<DIM;j++){/* Loop over variables */
      /* Get function at  (X+2h) point */
      for(i=0;i<DIM;i++) temp[i] = x[i];/* Copy array */
      temp[j] = x[j]+2*h;/* Change the desired var */
      func_odeiv(t,temp,f,&nu);/* Compute all RHS's */
      for(i=0;i<DIM;i++) dfdx[i*DIM+j] += -f[i];
      /* Get function at (X+h) point */
      for(i=0;i<DIM;i++) temp[i] = x[i];
      temp[j] = x[j]+h;
      func_odeiv(t,temp,f,&nu);
      for(i=0;i<DIM;i++) dfdx[i*DIM+j] += 8*f[i];
      /* Get function at (X-h) point */
      for(i=0;i<DIM;i++) temp[i] = x[i];
      temp[j] = x[j]-h;
      func_odeiv(t,temp,f,&nu);
      for(i=0;i<DIM;i++) dfdx[i*DIM+j] += -8*f[i];
      /* Get function at (X-2h) point */
      for(i=0;i<DIM;i++) temp[i] = x[i];
      temp[j] = x[j]-2*h;
      func_odeiv(t,temp,f,&nu);
      for(i=0;i<DIM;i++) dfdx[i*DIM+j] += f[i];
      /* Dividing by 12h */
      for(i=0;i<DIM;i++) dfdx[i*DIM+j] /= 12*h;
    }
  }
  /* Printing */
  if(jac_flag)
    printf("Analytical J=\n");
  else
    printf("Numerical J=\n");
  for(i=0;i<DIM;i++){
    for(j=0;j<DIM;j++){
      printf("%lf ",dfdx[i*DIM+j]);}
    printf("\n");
  }

  return GSL_SUCCESS;
}

int ode_lin(double t, const double x_lin[], double f_lin[], const double
	    dfdx[], double dfdt[], void *params)
{
  MU nu = *(MU *)params;
  int i,j;
  
  for(i=0;i<DIM;i++){
    f_lin[i]=0.0;
    for(j=0;j<DIM;j++)
      f_lin[i] += dfdx[i*DIM+j]*x_lin[j];
  }

  return 0;
}
