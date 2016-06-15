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
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
/* It is the core algorithms for computing the main Lyapunov exponent
   for a given time series (needed for exploring chaotic attractors)
   and other Lyapunov exponents computation realted algorithms.*/  

#include "init.h"
#include "lyapunov.h"
#include <math.h>
#include <string.h>


double * mle (double *MLE)
{
  int i,j;
  /*
    mle(...) function computes the Main Lyapunov Exponent (MLE) which
    is the indicator of chaos based on the algorithm from A.Wolf,
    J.Swift, H.Swinney and J.Vastano (1985), Physica 16D, "Determining
    Lyapunov exponents from a time series". Before calling one should
    get on the attractor, which is presumably chaotic.
  */
  /*****************************************************************/
  /*
    We do not need to reconstruct the attractor, since we have the
    whole the attractor's set of time series from the numerical
    integration (`xs' and `ts' storage arrays).
    One needs to compute the chaotic attractor trajectory before
    calling this function and to use Runge-Kutta (4,5) method with
    fixed time step for integration or any other fixed time step
    method. It is better to use the same method for both trajectories.
  */

  /*
    Set the initial distance DI between reference (fiducial) and test
    trajectories. It is defined in lyap_init(...)
    DF is the current(final) distance during the integration.
  */
  double DF=0.0;

  /*
    We want to set a new test initial point DI distant from the
    reference orbit initial point.  We want also the separation to be
    equal for all directions of the phase space. That is | x - x0 | =
    DI/sqrt(DIM) for the each direction, where x0 is reference
    init. point and x is test init. point. We assume, that absolute
    values | x - x0 | is resolved to the positive difference, i.e. x =
    x0 + DI/sqrt(DIM).
  */
  double z0[DIM];/*Previous point of the test orbit*/
  double z1[DIM];/*Next point of the test orbit*/
  for(i=0;i<DIM;i++)
    z0[i]=DI/sqrt(DIM) + xs[0][i];

  /*
    The threshold we allow two orbits to be separated at max is DTR.
    It is defined outside: lyap_init(...)
  */
  
  /*
    We determine time step for the integration, it must be the same
    used for getting the reference trajectory.
  */
  double DT=ts[1]-ts[0];

  /*
      Rescaling parameter: how big the distance became with time.
  */
  double A=0.0;
  
  /*
    We go through the whole trajectory
  */
  j=0;/*Number of point of ref. trajectory*/
  /*
    Next we check whether some there is anything that was integrated.
  */
  if (j>=(write_count-1)){
    fprintf(stderr,"No stored trajectory.\n");
    return MLE;
  }
  while(j<(write_count-1)){
    /*
      We integrate our test trajectory before the distance will overcome
      the threshold.
    */
    /*Copy initial vector*/
    for(i=0;i<DIM;i++)
      z1[i]=z0[i];
    while (DF <= DTR){
      j++;
      /*Integrator*/
      solver(z1,DT,func_odeiv);
      /*Compute the distance DF*/
      DF=0.0;
      for(i=0;i<DIM;i++)
	DF=DF+pow((z1[i]-xs[j][i]),2);
      DF=sqrt(DF);
      if(j==write_count-1)
	break;
    }
    if(j==write_count-1)
      break;
    /*Set the rescaling parameter*/
    A=DF/DI;

    /*
      Cumulative sum of the Lyapunov exponents.
    */
    *MLE = *MLE + log(A)/log(2);

    /*
      The rescaling: looking for the replacement point...
    */
    for(i=0;i<DIM;i++)
      z0[i]=DI/sqrt(DIM) + xs[j][i];
    DF=0.0;
   /*  for(i=0;i<DIM;i++) */
/*       z0[i]=xs[j][i]+(z1[i]-xs[j][i])/A; */
  }
  /*
    Final value for MLE
  */
  *MLE = *MLE / (ts[write_count-1]-ts[0]);
  printf("DI = %G; DTR = %G\n",DI,DTR);
  printf("|-----------------------------\n");
  printf("|\t MLE = %G\n",*MLE);
  printf("|-----------------------------\n");
  return MLE;
}

int mle_cont()
{
  int i,j;
  double tmp;
  FILE *output;
  double LT1=mynum.total_time;
  output=fopen(mle_file,"w");
  if( (mu.P[LPI] < LPARMIN) || (mu.P[LPI] > LPARMAX) )
    {
      printf("The parameter (%s) ",par_name[LPI]);
      printf("is out of the specified range");
      printf(": [%G,%G].\n",LPARMIN,LPARMAX);
      return 1;
    }
  /*
    Making head of the output file
  */
  fprintf(output,"# MLE computation summary:\n");
  fprintf(output,"# Integr. method: %s; ",mynum.method);
  fprintf(output,"dt = %G\n",mynum.step);
  fprintf(output,"# Total time: LT1 = %G\n",LT1);
  fprintf(output,"# Initial conditions: \n");
  fprintf(output,"# ");
  for(i=0;i<DIM;i++)
    fprintf(output,"%G ",x[i]);
  fprintf(output,"\n");
  fprintf(output,"# Par.: %s = %G\n",par_name[LPI],mu.P[LPI]);
  fprintf(output,"# Par. step: LDP = %G\n",LDP);
  if(LDP>0)
    fprintf(output,"# Range: [%G, %G]\n",mu.P[LPI],LPARMAX);
  if(LDP<0)
    fprintf(output,"# Range: [%G, %G]\n",LPARMIN,mu.P[LPI]);
  fprintf(output,"# Init. distance: DI = %G\n",DI);
  fprintf(output,"# Threshold distance: DTR = %G\n",DTR);
  fprintf(output,"\n\n");
  fprintf(output,"# %s\t",par_name[LPI]);
  fprintf(output,"MLE\n");
  while( (mu.P[LPI] >= (LPARMIN-1E-10)) && (mu.P[LPI] <= (LPARMAX+1E-10)) )
    {
      printf("%s = %G\n",par_name[LPI],mu.P[LPI]);
      /*
	Additional integraion to get the attractor. Any
	integrator. Time is LT1, no output to the file,
	transient. Using mynum.trans_time.
      */
      printf("Integrating transiently...\n");
      run(mynum.method,mynum.trans_time,0,1,1);
      /*
	Fixed step method integration on the attractor.
      */
      /*Changing method first*/
      if(strcmp(method,"run-kut4")!=0){
	printf("Changing the method to the Run-Kut 4...\n");
	strncpy(method,"run-kut4",(size_t)strlen("run-kut4"));
      }
      /*Changing writing step*/
      if(mynum.write_step != 1){
	printf("Changing the writing step to 1...\n");
	mynum.write_step=1;
      }
      /*Integrating: no writing data to file, not transient*/
      printf("Forming the fiducial trajectory...\n");
      run(mynum.method,LT1,0,0,1);
      /*
	Given the trajecory run mle(...) function
      */
      tmp=0.0;
      mle(&tmp);
      /*
	Writing output
      */
      fprintf(output,"%lf\t%lf\n",mu.P[LPI],tmp);
      /*
	Apply step for the parameter
      */
      mu.P[LPI] = mu.P[LPI] + LDP;
    }
  /*
    Get the parameter back to correspond the trajectory.
  */
  mu.P[LPI] = mu.P[LPI] - LDP;
  fclose(output);
  
  return 0;
}

double * lyap_spec(double *LCE)/*Lyapunov Characteristic Exponents,LCE is
			    parameter to the function*/
{
  /*The function computes the Lyapunov Characteristic Exponent
    spectrum based on the motion of nonlinear and linear systems and
    Gram-Schmidt reorthonormalization process at each step of the
    integrator. The algorithm described in A.Wolf, J.Swift, H.Swinney
    and J.Vastano (1985), Physica 16D, "Determining Lyapunov exponents
    from a time series".*/
  int i,j,k,npr;/*loop var*/
  double t1;/*Time*/
  /*Assume x[] on the attractor. x[] is globally defined.*/
  /*Create vectors for linear system, so called
   ORTHONORMAL FRAME. x[i] is i-th vector.*/
  double xlin[DIM][DIM];double f_lin[DIM];
  double znorm[DIM];/*Normalization factor*/
  double GSC[DIM-1];/*Gram-Schmidt coefficients*/
  /*Initial conditions for linear system*/
  for(i=0;i<DIM;i++)
    for(k=0;k<DIM;k++)
      xlin[i][k]=0.0;
  /*Nullify LCE and make frame orthonormal*/
  for(i=0;i<DIM;i++){
    xlin[i][i]=1.0;
    LCE[i]=0.0;
  }

  t=0.0;t1=mynum.total_time;
  double h=mynum.step;
  npr=1;
  if((strcmp(mynum.method,"run-kut4")!=0) &&	\
     (strcmp(mynum.method,"eu")!=0)){
    fprintf(stderr,"Change to fixed step integrator.\n");
    return LCE;
  }
  while(t <= t1){
    /*Calling integrator*/
    for(i=0;i<DIM;i++)
      solver_lin(x,xlin[i],f_lin,h,ode_lin);
    solver(x,h,func_odeiv);
    
    /*Normalize first vector*/
    znorm[0] = 0.0;
    for(i=0;i<DIM;i++)
      znorm[0] += pow(xlin[0][i],2);
    znorm[0] = sqrt(znorm[0]);
    for(i=0;i<DIM;i++)
      xlin[0][i] = xlin[0][i]/znorm[0];
    
    /*New orthonormal set of vectors*/
    for(j=1;j<DIM;j++){
      /*Generate j Gram-Schmidt reorthonormalization coefficients*/
      for(k=0;k<j;k++){/*over vectors*/
	GSC[k]=0.0;
	for(i=0;i<DIM;i++)/*over components in each vector*/
	  GSC[k] += xlin[j][i]*xlin[k][i];
      }

      /*Construct a new vector*/
      for(k=0;k<DIM;k++)/*over components*/
	for(i=0;i<j;i++)/*over vectors*/
	  xlin[j][k] = xlin[j][k] - GSC[i]*xlin[i][k];

      /*Calculate the vector's norm*/
      znorm[j]=0.0;
      for(k=0;k<DIM;k++)
	znorm[j] += pow(xlin[j][k],2);
      znorm[j] = sqrt(znorm[j]);
      /*Normalize the new vector*/
      for(k=0;k<DIM;k++)
	xlin[j][k] = xlin[j][k]/znorm[j];
    }
    /*Update running vector magnitudes*/
    for(k=0;k<DIM;k++)
      LCE[k] += log(znorm[k])/log(2);
    t += h;
    /*Printing*/
    if(t > (t1/10)*npr){/*10 printouts*/
      printf("t = %G: ",t);
      for(i=0;i<DIM;i++)
	printf("%G ",LCE[i]/t);
      printf("\n");
      npr++;
    }
  }
  printf("|-----------------------------\n");
  for(i=0;i<DIM;i++){
    LCE[i]=LCE[i]/t;
    printf("|\t LCE(%d) = %G\n",i+1,LCE[i]);}
  printf("|-----------------------------\n");
  
  return LCE;
}

int lyap_spec_cont()
{
  int i,j;
  double tmp[DIM];
  FILE *output;
  double LT1=mynum.total_time;
  output=fopen(lyap_spec_file,"w");
  if( (mu.P[LPI] < LPARMIN) || (mu.P[LPI] > LPARMAX) )
    {
      printf("The parameter (%s) ",par_name[LPI]);
      printf("is out of the specified range");
      printf(": [%G,%G].\n",LPARMIN,LPARMAX);
      return 1;
    }
  /*
    Making head of the output file
  */
  fprintf(output,"# LCE spectrum computation summary:\n");
  fprintf(output,"# Integr. method: %s; ",mynum.method);
  fprintf(output,"dt = %G\n",mynum.step);
  fprintf(output,"# Total time: LT1 = %G\n",LT1);
  fprintf(output,"# Initial conditions: \n");
  fprintf(output,"# ");
  for(i=0;i<DIM;i++)
    fprintf(output,"%G ",x[i]);
  fprintf(output,"\n");
  fprintf(output,"# Par.: %s = %G\n",par_name[LPI],mu.P[LPI]);
  fprintf(output,"# Par. step: LDP = %G\n",LDP);
  if(LDP>0)
    fprintf(output,"# Range: [%G, %G]\n",mu.P[LPI],LPARMAX);
  if(LDP<0)
    fprintf(output,"# Range: [%G, %G]\n",LPARMIN,mu.P[LPI]);
  fprintf(output,"\n\n");
  fprintf(output,"# %s\t",par_name[LPI]);
  for(i=0;i<DIM;i++)
    fprintf(output,"LCE(%d)\t",i+1);
  fprintf(output,"\n");
  while( (mu.P[LPI] >= (LPARMIN-1E-10)) && (mu.P[LPI] <= (LPARMAX+1E-10)) )
    {
      printf("%s = %G\n",par_name[LPI],mu.P[LPI]);
      /*Call lyap_spec algorithm. Initial point in x[] must be on the
	attractor. Be sure about it before the call.*/
      lyap_spec(tmp);
      /*
	Writing output
      */
      fprintf(output,"%lf ",mu.P[LPI]);
      for(i=0;i<DIM;i++)
	fprintf(output,"%lf ",tmp[i]);
      fprintf(output,"\n");
      /*
	Apply step for the parameter
      */
      mu.P[LPI] = mu.P[LPI] + LDP;
    }
  /*
    Get the parameter back to correspond the trajectory.
  */
  mu.P[LPI] = mu.P[LPI] - LDP;
  
  fclose(output);
  return 0;
}

int lyap_init()
{
  DI = 1e-7;
  DTR = 1e-5;
  LPARMAX = 30.0;
  LPARMIN = 15.0;
  LDP = 0.1;
  LPI = 0;
  mle_file = (char *)malloc(150*sizeof(char));
  lyap_spec_file = (char *)malloc(150*sizeof(char));
  mle_file = strcpy(mle_file,"dat.mle");
  lyap_spec_file = strcpy(lyap_spec_file,"dat.lces");
  
  return 0;
}

void lyap_free()
{
  free(mle_file);
  free(lyap_spec_file);
}
