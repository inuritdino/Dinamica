/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012,2013,2014 Elias Potapov. */
/* Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007,
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
#include "errors.h"
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_eigen.h>/* For eigenvalues */
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_permutation.h>
#define TRUE 1
#define FALSE 0
//size of the buffer defined in trajectory determination function traj()
#define MAX_N_CMPPAIRS 100
#define MIN_ISEC 6
#define MIN_PEAK 6

int traj_init()
{
  traj_verb = 0;
  eps_traj_abs = 0.001;
  eps_traj_rel = 0.05;
  eps_traj_rel_hg = 0.05;
  
  return 0;
}

double D(double eps_abs, double eps_rel, double value) {
  return(eps_abs + eps_rel*value);
}
/* int traj() */
/* { */
/*   int i,j,k,l,m; */
/*   int n=LDIM; */
/*   int p=(int)DIM/LDIM; */
/*   int min_n_isec[p];//Min amount of intersections, that all vars have. */
/*   int min_n_peak[p];//Min amount of peaks, that all vars have. */
	
/*   double E_tper[MAX_N_ISEC][MAX_N_CMPPAIRS][p];//Errors for intersections. */
/*   /\* You should uncomment next line if You want counting for peak height in regime *\/ */
/*   /\*+definition. *\/ */
/*   //double E_peak[MAX_N_ISEC][MAX_N_CMPPAIRS][MAX_LOCAL_DIM];//Errors for peaks. */
/*   double E_t_peak[MAX_N_PEAKS][MAX_N_CMPPAIRS][p];//Errors for time of peaks. Not used. */
/*   double E,E1;//E---errors for amplitude. E1---helping error. */
/*   int count[p]; */
/*   int count_a[p]; */
/*   int short cmppairs; //Number of compared variables pairs. */
/*   float period=0.0; */
/*   double regime_ratio[p]; */
/*   /\*We add numerical values for each dynamical state to the */
/*     value of `regime' variable each time particular regime is */
/*     determined. It would help to debug wrong determination of */
/*     dynamical state. So first we need to nullify `regime' and */
/*     `reg[i]'.*\/ */
/*   regime=0; */
/*   for(i=0;i<DIM;i++) */
/*     reg[i]=0; */
/*   /\******************\/ */
/*   /\**STEADY STATE**\/ */
/*   for(i=0;i<DIM;i++) */
/*     ss_factor[i]=0.0; */
/*   steady_state(ss_factor); */
/*   /\*Here comes criterion of the SS existence*\/ */
/*   j=0;k=0; */
/*   for(i=0;i<DIM;i++) */
/*     if(ss_factor[i]>0.5){/\*50% criterion*\/ */
/*       reg[i]=SS;/\*We have SS for the i-th var*\/ */
/*       j++; */
/*     } */
/*   E=0;E1=0;/\*errors nullifying*\/ */
/*   if(j==DIM){ */
/*     regime+=SS; */
/*     for(i=0;i<DIM;i++){ */
/*       for(l=i;l<DIM-p;l+=p){ */
/* 	E1=xs[n_steps-1][l]-xs[n_steps-1][l+p]; */
/* 	if(E1>0) E1=-E1; */
/* 	E+=E1; */
/*       } */
/*       E=E/LDIM;/\*In case LDIM=1 => E=0.*\/ */
/*       /\*Then we decide what is bigger: ihss variables or hss. In */
/* 	general, all vars should the same type(hss or ihss), but in */
/* 	some case of bad error comparing they can be different.*\/ */
/*       if(E<eps_abs_am) k--; */
/*       if(E>eps_abs_am) k++; */
/*     } */
/*   } */
/*   if(k>0) */
/*     regime+=IHSS; */
/*   if(k<0) */
/*     regime+=HSS; */
/*   if(j==DIM && LDIM==1)/\*j==DIM means all variables in SS, i.e. reg[i]s==SS*\/ */
/*     regime=SS;/\*There is no homogeneous or inhomogeneous states */
/* 		in single(not coupled to any other) system*\/ */
/*   print_regime(); */
/*   if(regime>0)/\*Regime is SS+...so SS is the dynamical state*\/ */
/*     return 0; */
/*   /\**END OF STEADY STATE**\/ */
/*   /\*========================================================*\/ */
/*   for(i=0; i<n; i++){ */
/*     count[i]=0; */
/*     count_a[i]=0; */
/*     min_n_isec[i]=n_isec[0]; */
/*     min_n_peak[i]=n_peaks[0]; */
/*     regime_ratio[i]=1.0;/\*This is necessary for correct */
/* 			  trajectory determination*\/ */
/*   } */
/*   /\*****Creating min of INTERSECTIONS.******\/ */
/*   for(k=0;k<n;k++){/\*...intersections.*\/ */
/*     for(i=k;i<DIM;i+=n){ */
/*       if(n_isec[i] < min_n_isec[k]) */
/* 	min_n_isec[k] = n_isec[i]; */
/*     } */
/*   } */
/*   /\**Forming counts (l), if l>0 then additional integration.**\/ */
/*   j=0; */
/*   while(per_method == 0){ */
/*     for(k=0;k<n;k++){ */
/*       l=0; */
/*       if(min_n_isec[k] < MIN_ISEC) */
/* 	l++; */
/*     } */
/*     m=0; */
/*     for(i=0;i<DIM;i++){ */
/*       if(per_ratio[i]>eps_per_ratio) */
/* 	m++; */
/*     } */
/*     if(l==0 && m==0) */
/*       break;/\* We have enough to form per_ratios(l==0) and */
/* 	       convergence to L.C(m==0). *\/ */

/*     /\* Additional integration, getting attractor. *\/ */
/*     if(j<traj_trans){ */
/*       if(l>0){ */
/* 	fprintf(stdout,"Too few intersections.\n"); */
/* 	fprintf(stdout,"Getting attractor...%d\n",j+1); */
/* 	run(mynum.method,mynum.trans_time,0,1,1); */
/* 	fprintf(stdout,"Done.\n"); */
/*       } */
/*       if(m>0){ */
/* 	fprintf(stdout,"Did not get L.C\n"); */
/* 	fprintf(stdout,"Getting attractor...%d\n",j+1); */
/* 	run(mynum.method,mynum.trans_time,0,1,1); */
/* 	fprintf(stdout,"Done.\n"); */
/*       } */
/*     } */
/*     j++; */
/*     if(j==traj_trans){ */
/*       not_get_flag=TRUE; */
/*       break; */
/*     } */
/*   } */
/*   /\*****Creating min of PEAKS.******\/ */
/*   for(k=0;k<n;k++){ */
/*     for(i=k;i<DIM;i+=n){ */
/*       if(n_peaks[i] < min_n_peak[k]) */
/* 	min_n_peak[k] = n_peaks[i]; */
/*     } */
/*   } */
/*   /\**Forming counts (l), if l>0 then additional integration.**\/ */
/*   j=0; */
/*   while(per_method == 1){ */
/*     l=0; */
/*     for(i=0;i<n;i++){ */
/*       if(min_n_peak[i] < MIN_PEAK) */
/* 	l++; */
/*     } */
/*     m=0; */
/*     for(i=0;i<DIM;i++){ */
/*       if(per_ratio[i]>eps_per_ratio) */
/* 	m++; */
/*     } */
/*     if(l==0 && m==0) */
/*       break;/\* We have enough to form per_ratios(l==0) and */
/* 	       convergence to L.C(m==0).. *\/ */
	  
/*     /\* Additional integration, getting attractor. *\/ */
/*     if(j<traj_trans){ */
/*       if(l>0){ */
/* 	fprintf(stdout,"Too few peaks.\n"); */
/* 	fprintf(stdout,"Getting attractor...%d\n",j+1); */
/* 	run(mynum.method,mynum.trans_time,0,1,1); */
/* 	fprintf(stdout,"Done.\n"); */
/*       } */
/*       if(m>0){ */
/* 	fprintf(stdout,"Did not get L.C\n"); */
/* 	fprintf(stdout,"Getting attractor...%d\n",j+1); */
/* 	run(mynum.method,mynum.trans_time,0,1,1); */
/* 	fprintf(stdout,"Done.\n"); */
/*       } */
/*     } */
/*     j++; */
/*     if(j==traj_trans){ */
/*       not_get_flag=TRUE; */
/*       break; */
/*     } */
/*   } */

/*   /\*============================================================================== */
/*     FORMING ERRORS */
/*     ==============================================================================*\/ */
/*   /\* Creating calculation errors for intersections. A bit of sophisticated code.*\/ */
/*   if(per_method == 0){ */
/*     for(k=0; k < n; k++){ //Looping over certain type of var(local dim). */
/*       cmppairs=0; */
/*       for(l=0; l < p-1; l++){//How many vars of certain type(k). */
/* 	for(i=k+l*n; i < DIM-n; i+=n){//Actual stepping over certain type(local dim) of var. */
/* 	  for(j=0; j < min_n_isec[k]; j++)//Looping over intersections. */
/* 	    E_tper[j][cmppairs][k] = tper[j][k+l*n] - tper[j][i+n]; */
/* 	  cmppairs+=1; */
/* 	  if(cmppairs > MAX_N_CMPPAIRS){fprintf(stdout,"Warning(tper):too much values to store.Missing some.\n"); */
/* 	    fprintf(stdout,"Try to descrease total time of integration.\n");break; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   // Creating calculation errors for peaks and time of peaks. Peculiar enough. */
/*   if(per_method == 1){ */
/*     for(k=0; k < n; k++){//Looping over certain type of var(local dim). */
/*       cmppairs=0; */
/*       for(l=0; l < p-1; l++){//How many vars of certain type(k). */
/* 	for(i=k+l*n; i < DIM-n; i+=n){//Actual stepping over certain type(local dim) of var. */
/* 	  for(j=0; j < min_n_peak[k]; j++){//Looping over peaks. */
/* 	    // If You decide to uncomment next line You should add for counting */
/* 	    //+to the following code. It is superfluous. Not recommended. */
/* 	    // Moreover, You should uncomment initializing string at the begining of this */
/* 	    //+function. */
/* 	    //E_peak[j][cmppairs][k] = x_peak[j][k+l*n] - x_peak[j][i+n]; */
/* 	    E_t_peak[j][cmppairs][k] = t_peak[j][k+l*n] - t_peak[j][i+n]; */
/* 	  } */
/* 	  cmppairs+=1; */
/* 	  if(cmppairs > MAX_N_CMPPAIRS){fprintf(stdout,"Warning(peaks):too much values to store.Missing some.\n"); */
/* 	    fprintf(stdout,"Try to descrease total time of integration.\n");break; */
/* 	  } */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   //================================================================================== */
/*   //================================================================================== */
/*   //\****************************************************************************************** */
/*   // FORMATION OF COUNT[K]S */
/*   // Comparing calculation errors with certain level of error defined by the user(D-function). */
/*   // ******************************************************************************************* */

/*   if(per_method == 0){ */
/*     for(k=0; k<n; k++){ */
/*       period=0.0; */
/*       for(l=k;l<DIM;l+=n){ period=period + big_per[l]/p;} */
/*       for(j=0; j<cmppairs; j++){ */
/* 	for(i=0; i<min_n_isec[k];i++){ */
/* 	  if(E_tper[i][j][k] < 0) */
/* 	    E_tper[i][j][k] = -E_tper[i][j][k]; */
/* 	  if(E_tper[i][j][k] > D(eps_abs_tper,eps_rel_tper,period)) */
/* 	    count[k]+=1; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   if(per_method == 1){ */
/*     for(k=0; k<n; k++){ */
/*       period=0.0; */
/*       for(l=k;l<DIM;l+=n){period=period + big_per[l]/p;} */
/*       for(j=0; j<cmppairs; j++){ */
/* 	for(i=0; i<min_n_peak[k];i++){ */
/* 	  if(E_t_peak[i][j][k] < 0) */
/* 	    E_t_peak[i][j][k] = -E_t_peak[i][j][k]; */
/* 	  if(E_t_peak[i][j][k] > D(eps_abs_peak,eps_rel_peak,period)) */
/* 	    count[k]+=1; */
/* 	} */
/*       } */
/*     } */
/*   } */
/*   //COMPARING AMPLITUDES */
/*   for(k=0;k<n; k++){ */
/*     for(i=0;i<p-1;i++){ */
/*       for(j=k+i*n; j<DIM-n; j+=n){ */
/* 	E=ampl[k+i*n] - ampl[j+n]; */
/* 	if( E < 0 ) */
/* 	  E = -E; */
/* 	if(E > D(eps_abs_am,eps_rel_am,(float)ampl[k+i*n])) */
/* 	  count_a[k]+=1; */
/*       } */
/*     } */
/*   } */
/*   //================================================================================================ */
/*   //FORMING REGIME RATIOS FROM COUNTS. */
/*   if(per_method == 0){ */
/*     for(i=0; i<n ; i++) */
/*       regime_ratio[i] = (double)count[i]/((double)cmppairs*(double)min_n_isec[i]); */
/*   } */
/*   if(per_method == 1){ */
/*     for(i=0; i<n ; i++) */
/*       regime_ratio[i] = (double)count[i]/((double)cmppairs*(double)min_n_peak[i]); */
/*   } */
/*   for(i=0; i<n; i++) */
/*     if(regime_ratio[i] > 1.0) */
/*       fprintf(stdout,"RATIO GREATER THAN ONE!"); */
/*   //========================================================================================================== */

/*   //REPORTING ABOUT IN-PHASE, ANTI-PHASE and ALL-ANTI-PHASE("WAVE") REGIMES IF THEY EXIST.ALSO ABOUT IHLC. */
/*   //fprintf(stdout,"Regime ratios:\n"); */

/*   regime=LC; */
/*   for(i=0;i<n;i++){ */
/*     reg[i]=LC; */
/*     if( count_a[i] != 0){ */
/*       regime=LC+IHLC; */
/*       reg[i]+=IHLC; */
/*       continue;     */
/*     } */
/*     if(regime_ratio[i] < eps_inregime_ratio){ */
/*       regime=HLC+INLC; */
/*       reg[i]+=HLC+IHLC; */
/*       continue;     */
/*     } */
	  
/*     if(per_method == 0){ */
/*       E = ((double)cmppairs-1.0)*(double)min_n_isec[i]; */
/*       E1 = (double)cmppairs*(double)min_n_isec[i];} */
/*     if(per_method == 1){ */
/*       E = ((double)cmppairs-1.0)*(double)min_n_peak[i]; */
/*       E1 = (double)cmppairs*(double)min_n_peak[i];} */
	  
/*     if((regime_ratio[i] > eps_inregime_ratio) && (regime_ratio[i] <= E/E1)){ */
/*       regime=HLC+OUTLC;/\*PANTI*\/ */
/*       reg[i]+=HLC+OUTLC;/\*PANTI*\/ */
/*       continue;     */
/*     } */
/*     if(regime_ratio[i] > E/E1){ */
/*       regime=HLC+ANTILC;/\*PANTIALL*\/ */
/*       reg[i]+=HLC+ANTILC;/\*PANTIALL*\/ */
/*       continue;     */
/*     } */

/*   } */
/*   /\******************************************************************************\/ */
/*   for(i=0;i<n;i++){ */
/*     if(reg[i]!=0 && (reg[i]==HLC+ANTILC)){/\* If any of var type is in PANTIALL regime then*\/ */
/*       regime=HLC+ANTILC;            /\* total REGIME is PANTIALL. *\/ */
/*       break; */
/*     } */
/*   } */
/*   for(i=0;i<n;i++){ */
/*     if(reg[i]!=0 && (reg[i]==IHLC)){/\* If any of trajectory is IHLC then total REGIME is IHLC. *\/ */
/*       regime=LC+IHLC; */
/*       break; */
/*     } */
/*   } */
/*   for(i=0;i<n;i++) */
/*     fprintf(stdout,"Var type #%d: %d\n",i+1,reg[i]); */
	
/*   fprintf(stdout,"REGIME: %d\n",regime); */
/*   return 0; */
/* } */

int get_varInd(int init_varInd, int varInds[])
{
  /* The function calculates the variable indices for the multi-dimensional dynamics
     checking. varInds[] must be of size LDIM (number of sub-systems). NOTE: var
     indices start from 0 (C-style). */
  int i,j,k;
  i = (int)DIM/LDIM;/* number of variables in a sub-system */
  /* We need the offset: the number of sub-system where initial var index init_varInd
     is pointing to. */
  j = 1;/* the initial offset is unity */
  /* Calculate the offset */
  while(init_varInd >= j*i)
    j++;
  /* Calculate which (0-based row number) element in a sub-system init_varInd points to. */
  for(k=0;k<LDIM;k++){
    varInds[k] = (init_varInd - (j-1)*i) + k*i;
    //printf("varInds[%d] = %d\n",k,varInds[k]);
  }
  
  return 0;
}

void free_regime_stat(regStat *stat)
{/* Free the regStat structure */
  if(stat != NULL){
    free(stat->os_pd);
    free(stat->os_up);
    free(stat->os_hg_tr);
    free(stat->os_ip_tr);
    free(stat->os_op_tr);
    free(stat);
  }
}

void copy_regime_stat(regStat *dst, regStat *src)
{
  int i;
  dst->N = src->N;
  dst->ss = src->ss;
  dst->os = src->os;
  dst->os_pd = (int *)malloc(dst->os*sizeof(int));
  for(i=0;i<dst->os;i++)
    dst->os_pd[i] = src->os_pd[i];
  dst->os_up_cn = src->os_up_cn;
  dst->os_up = (int *)malloc(dst->os_up_cn*sizeof(int));
  for(i=0;i<dst->os_up_cn;i++)
    dst->os_up[i] = src->os_up[i];
  dst->hg = src->hg;
  dst->ih = src->ih;
  dst->ss_hg = src->ss_hg;
  dst->os_hg = src->os_hg;
  dst->os_hg_tr = (int *)malloc(dst->os*sizeof(int));
  for(i=0;i<dst->os;i++)
    dst->os_hg_tr[i] = src->os_hg_tr[i];
  dst->os_ih = src->os_ih;
  dst->ss_ih = src->ss_ih;
  dst->os_ip = src->os_ip;
  dst->os_ip_tr = (int *)malloc(dst->os*sizeof(int));
  for(i=0;i<dst->os;i++)
    dst->os_ip_tr[i] = src->os_ip_tr[i];
  dst->os_op = src->os_op;
  dst->os_op_tr = (int *)malloc(dst->os*sizeof(int));
  for(i=0;i<dst->os;i++)
    dst->os_op_tr[i] = src->os_op_tr[i];
  dst->ud = src->ud;
  dst->mx = src->mx;
}

regStat *regime_stat(regReport report[],const int report_size, FILE *out)
{/* Calculate the statistics of found regimes, given the array of reports. See
    trajectory.h for the struct members of the statistics structure.
    FILE *out is a stream for the output (out == NULL => no output).
    This function frees the report.
 */
  int i,j;
  /* Statistics struct */
  regStat *stat = (regStat *)malloc(sizeof(regStat));
  /* Init the stat */
  stat->N = report_size;/* total number of regimes analyzed */
  stat->ss = 0;
  stat->os = 0;
  stat->os_pd = (int *)malloc(sizeof(int));
  stat->os_up = (int *)malloc(sizeof(int));
  stat->os_up_cn = 0;
  stat->hg = 0;
  stat->ih = 0;
  stat->os_hg = 0;
  stat->ss_hg = 0;
  stat->os_hg_tr = (int *)malloc(sizeof(int));
  stat->os_ih = 0;
  stat->ss_ih = 0;
  stat->os_ip = 0;
  stat->os_ip_tr = (int *)malloc(sizeof(int));
  stat->os_op = 0;
  stat->os_op_tr = (int *)malloc(sizeof(int));
  stat->ud = 0;
  stat->mx = 0;
  for(i=0;i<report_size;i++){
    if(!report[i].regime){
      (stat->ss)++;
      if(report[i].homog)
	(stat->ss_hg)++;
      if(!report[i].homog)
	(stat->ss_ih)++;
    }
    if(report[i].regime > 0){
      (stat->os)++;
      stat->os_pd = (int *)realloc(stat->os_pd,stat->os * sizeof(int));
      stat->os_hg_tr = (int *)realloc(stat->os_hg_tr,stat->os * sizeof(int));
      stat->os_op_tr = (int *)realloc(stat->os_op_tr,stat->os * sizeof(int));
      stat->os_pd[stat->os-1] = report[i].regime;
      if(report[i].homog){
	(stat->os_hg)++;
	stat->os_hg_tr[stat->os - 1] = 1;/* HOS */
      }
      else if(!report[i].homog){
	(stat->os_ih)++;
	stat->os_hg_tr[stat->os - 1] = 0;/* IHOS */
      }
      if(report[i].ph_shift >= 0.0001){
	(stat->os_op)++;
	stat->os_op_tr[stat->os - 1] = 1;/* OPOS */
	stat->os_ip_tr[stat->os - 1] = 0;
      }
      else if((report[i].ph_shift < 0.0001) && (report[i].ph_shift > 0)){
	(stat->os_ip)++;
	stat->os_ip_tr[stat->os - 1] = 1;/* IPOS */
	stat->os_op_tr[stat->os - 1] = 0;
      }
      else{
	/* NOT IPOS and NOT OPOS (IHOS?) */
	stat->os_ip_tr[stat->os - 1] = 0;
	stat->os_op_tr[stat->os - 1] = 0;
      }
    }
    if(report[i].homog)
      (stat->hg)++;
    if(!report[i].homog)
      (stat->ih)++;
    if(report[i].regime == -1)
      (stat->ud)++;
    if(report[i].regime == -2)
      (stat->mx)++;
  }
  /* Analyze periodicity */
  if(stat->os){
    stat->os_up_cn++;
    stat->os_up[0] = stat->os_pd[0];
    for(i=1;i<stat->os;i++){
      j = 0;
      while((j < stat->os_up_cn) && (stat->os_up[j] != stat->os_pd[i]))
	j++;
      if(j == stat->os_up_cn){ /* unique solution */
	stat->os_up_cn++;
	stat->os_up = (int *)realloc(stat->os_up,stat->os_up_cn*sizeof(int));
	stat->os_up[stat->os_up_cn-1] = stat->os_pd[i];
      }
      if(j < stat->os_up_cn){} /* not unique solution */
    }
  }
  fprintf(stdout,"-------------------\n");
  fprintf(stdout,"Regimes statistics:\n");
  fprintf(stdout,"-------------------\n");
  fprintf(stdout,"Number of regimes: %d\n",stat->N);
  fprintf(stdout,"Steady States: %d (%G%%)\n",stat->ss,
	  100*((double)stat->ss/(double)stat->N));
  fprintf(stdout,"Oscillatory: %d (%G%%)\n",stat->os,
	  100*((double)stat->os/(double)stat->N));
  fprintf(stdout,"\tPeriodicity (unique):");
  for(i=0;i<stat->os_up_cn;i++)
    fprintf(stdout," %d",stat->os_up[i]);
  fprintf(stdout,"\n");
  fprintf(stdout,"Homogeneous steady state: %d (%G%%)\n",stat->ss_hg,
	  100*((double)stat->ss_hg/(double)stat->N));
  fprintf(stdout,"In-homogeneous steady state: %d (%G%%)\n",stat->ss_ih,
	  100*((double)stat->ss_ih/(double)stat->N));
  fprintf(stdout,"Homogeneous oscillatory: %d (%G%%, %G%% of homogeneous)\n",
	  stat->os_hg,100*((double)stat->os_hg/(double)stat->N),
	  100*((double)stat->os_hg/(double)stat->hg));
  fprintf(stdout,"In-homogeneous oscillatory: %d (%G%%, %G%% of in-homogeneous)\n",
	  stat->os_ih,100*((double)stat->os_ih/(double)stat->N),
	  100*((double)stat->os_ih/(double)stat->ih));
  if(stat->os){
    fprintf(stdout,"In-phase oscillatory: %d (%G%%, %G%% of oscillatory)\n",
	    stat->os_ip,100*((double)stat->os_ip/(double)stat->N),
	    100*((double)stat->os_ip/(double)stat->os));
    fprintf(stdout,"Out-of-phase oscillatory: %d (%G%%, %G%% of oscillatory)\n",
	    stat->os_op,100*((double)stat->os_op/(double)stat->N),
	    100*((double)stat->os_op/(double)stat->os));
  }
  fprintf(stdout,"Mixed: %d (%G%%)\n",stat->mx,100*((double)stat->mx/(double)stat->N));
  fprintf(stdout,"Undetermined: %d (%G%%)\n",stat->ud,
	  100*((double)stat->ud/(double)stat->N));
  /* Output to the file if any. Mostly for testing, statistics is shown above. */
  //statND_to_file(out,stat,1);

  free_report(report);
  return stat;
}

void statND_to_file(FILE *out, regStat *stat[], int stat_size)
{
  /*
    The output format will be:
    SS OS HSS  IHSS  HOS(1)...HOS(n)  IHOS(1)...IHOS(n) IPOS(1)...IPOS(n)  OPOS(1)...OPOS(n)
    where n signifies the number of the last periodicity found, e.g. if found 1,2
    and 4 periodicities the output would be: HOS(1)  HOS(2)  HOS(4)  IHOS(1) and
    so on. These are the column headers, whereas the real content is the number of
    the corresponding regimes.
  */
  if(out == NULL)
    return;
  int i,j,k,counter;
  /* Find all unique periodicities in the stat array */
  int up_cn;
  int *up = perdcty_in_stat_array(stat,stat_size,&up_cn);
  if(up == NULL){
    fprintf(stderr,
	    "Error(statND_to_file): unique periodicity could not be resolved.\n");
  }
  /************************** PRINT THE HEADER **************************/
  /* Print out some regimes in any case */
  fprintf(out,"SS OS HSS IHSS ");
  /* Print OS regimes only if there is some periodicity. Print all possible periodicities. */
  if(up_cn){
    for(i=0;i<up_cn;i++){
      fprintf(out,"HOS(%d) ",up[i]);
      printf("HOS(%d) \n",up[i]);
    }
    for(i=0;i<up_cn;i++){
      fprintf(out,"IHOS(%d) ",up[i]);
    }
    for(i=0;i<up_cn;i++){
      fprintf(out,"IPOS(%d) ",up[i]);
    }
    for(i=0;i<up_cn;i++){
      fprintf(out,"OPOS(%d) ",up[i]);
    }
  }
  /* Print MX and UD regimes at the end */
  fprintf(out,"MX UD\n");
  /*************************** PRINT THE VALUES **************************/
  for(j=0;j<stat_size;j++){
    /************ SS vs OS PRINT *************/
    fprintf(out,"%d ",stat[j]->ss);
    printf("SS=%d ",stat[j]->ss);
    fprintf(out,"%d ",stat[j]->os);
    printf("OS=%d \n",stat[j]->os);
    /************** SS PRINT *****************/
    fprintf(out,"%d ",stat[j]->ss_hg);
    fprintf(out,"%d ",stat[j]->ss_ih);
    if(stat[j]->os){/************** OS PRINT ****************/
      if(stat[j]->os_hg){/* HOS */
	counter = 0;
	for(i=0;i<up_cn;i++){/* HOS(up[i]) */
	  for(k=0;k<stat[j]->os;k++){
	    if((stat[j]->os_pd[k] == up[i]) && (stat[j]->os_hg_tr[k])){
	      /* Periodicity and Homogeneity satisfy, count... */
	      counter++;
	    }
	  }
	  fprintf(out,"%d ",counter);
	}
      }
      else
	fprintf(out,"0 ");
      if(stat[j]->os_ih){/* IHOS */
	counter = 0;
	for(i=0;i<up_cn;i++){/* IHOS(up[i]) */
	  for(k=0;k<stat[j]->os;k++){
	    if((stat[j]->os_pd[k] == up[i]) && (!stat[j]->os_hg_tr[k])){
	      /* Periodicity and non-Homogeneity satisfy, count... */
	      counter++;
	    }
	  }
	  fprintf(out,"%d ",counter);
	}
      }
      else
	fprintf(out,"0 ");
      if(stat[j]->os_ip){/* IPOS */
	counter = 0;
	for(i=0;i<up_cn;i++){/* IPOS(up[i]) */
	  for(k=0;k<stat[j]->os;k++){
	    if((stat[j]->os_pd[k] == up[i]) && (stat[j]->os_ip_tr[k])){
	      /* Periodicity and in-phase-ness(not out-of-phase) satisfy, count... */
	      counter++;
	    }
	  }
	  fprintf(out,"%d ",counter);
	}
      }
      else
	fprintf(out,"0 ");
      if(stat[j]->os_op){/* OPOS */
	counter = 0;
	for(i=0;i<up_cn;i++){/* OPOS(up[i]) */
	  for(k=0;k<stat[j]->os;k++){
	    if((stat[j]->os_pd[k] == up[i]) && (stat[j]->os_op_tr[k])){
	      /* Periodicity and out-of-phase-ness satisfy, count... */
	      counter++;
	    }
	  }
	  fprintf(out,"%d ",counter);
	}
      }
      else
	fprintf(out,"0 ");
    }
    if(stat[j]->mx)/****************** MX PRINT *****************/
      fprintf(out,"%d ",stat[j]->mx);
    else
      fprintf(out,"0 ");
    if(stat[j]->ud)/****************** UD PRINT *****************/
      fprintf(out,"%d ",stat[j]->ud);
    else
      fprintf(out,"0 ");
    fprintf(out,"\n");
  }

  free(up);
}

int *perdcty_in_stat_array(regStat *stat[],int stat_size, int *pd_count)
{/* Find all unique periodicities observed in the stat array without repetitions. */
  int i,j,k;
  int *perdcty = (int *)malloc(sizeof(int));
  if(perdcty == NULL){
    fprintf(stderr,"Error: perdcty array was not allocated\n");
    return NULL;
  }
  *pd_count = 0;
  for(i=0;i<stat_size;i++){
    //printf("STAT[%d].os_up_cn = %d\n",i,stat[i]->os_up_cn);
    //printf("STAT[%d].os_pd[0] = %d\n",i,stat[i]->os_pd[0]);
    for(k=0;k<stat[i]->os_up_cn;k++){
      //printf("STAT[%d].os_up[%d] = %d\n",i,k,stat[i]->os_up[k]);
      j = 0;
      while((j < *pd_count) && (perdcty[j] != stat[i]->os_up[k]))
	j++;
      if(j == *pd_count){/* unique periodicity */
	(*pd_count)++;
	perdcty = (int *)realloc(perdcty,*pd_count*sizeof(int));
	if(perdcty == NULL){
	  fprintf(stderr,"Error: perdcty array was not re-allocated\n");
	  return NULL;
	}
	perdcty[*pd_count-1] = stat[i]->os_up[k];
	//printf("PERDCTY[%d]=%d\n",*pd_count-1,stat[i]->os_up[k]);
      }
      if(j < *pd_count){/* Not unique periodicity */}
    }
  }

  return perdcty;
}

char *report_translate(regReport report)
{/* This function translates the dynamics regimes report into the human language */

  /* Output character string */
  char *out;
  int out_len;
  /* Try to estimate the length of the string */
  /* Every 5 is for two letter code + parentheses. Every 6 is for numerical value
     inside the parentheses. Homogeneity report requires two numerical values inside
     parentheses, thus having 15. */
  out_len = 10 + 15 + 10 + 15 + 10 + 30;/* e.g. OS-1(T=13.5)/IP(0T)/H(1,1) */
  out = (char *)malloc(out_len*sizeof(char));
  if(out == NULL){
    fprintf(stderr,"Error: cannot allocate memory (report_translate,out)\n");
    return NULL;
  }
  memset(out,'\0',(size_t)out_len);
  /* Check if the report is NULL, when the check was not successful. */
  if(&report == NULL){
    out = strcat(out,"CANNOT TRANSLATE THE REPORT!");
    return out;
  }
  char *tmp = (char *)malloc(25*sizeof(char));/* 25 is enough */
  if(tmp == NULL){
    fprintf(stderr,"Error: cannot allocate memory (report_translate,tmp)\n");
    return NULL;
  }
  if(!report.regime)
    out = strcat(out,"SS(0)");
  else if(report.regime > 0){
    sprintf(tmp,"OS-%d(T=%G)",report.regime,report.period);
    out = strcat(out,tmp);
    /* Report phase only if OS and LDIM > 1 */
    if(LDIM > 1){
      /* Out-of-phase oscillations */
      if(report.ph_shift > 0.0001){
	out = strcat(out,"/OP");
      }
      else if(report.ph_shift >= 0 && report.ph_shift < 0.0001){
	/* In-phase oscillations, i.e. phase shift = 0 */
	out = strcat(out,"/IP");
      }
      else{/* No phase shift for negative values: IHOS */}
      if(report.ph_shift >= 0)
	sprintf(tmp,"(%GT)",report.ph_shift);
      else
	sprintf(tmp,"");
      out = strcat(out,tmp);
    }
  }
  else if(report.regime == -1){
    fprintf(stdout,"Warning: regime was not determined\n");
    out = strcat(out,"-/");
  }
  else if(report.regime == -2){
    fprintf(stdout,"Warning: both stationary and non-stationary dynamics were detected\n");
    out = strcat(out,"mixed/");
  }
  else
    fprintf(stderr,"Error: regime cannot be equal %d\n",report.regime);
  /* Report homogeneity: only for LDIM > 1 */
  if(LDIM > 1){
    if(report.homog){
      out = strcat(out,"/H");
    }
    else{
      out = strcat(out,"/IH");
    }
    sprintf(tmp,"(bg=%G,ag=%G)",report.b_gain,report.a_gain);
    out = strcat(out,tmp);
  }

  free(tmp);
  return out;
}

regReport *get_sys_dynamics(regReport *report,const int report_size)
{/* Wrapper function for calling get_dynamics_1d() and get_dynamics_nd() altogether with
    all intermediate variable definitions etc. */
  int j;
  /* Report output */
  report = (regReport *)realloc(report,(report_size+1)*sizeof(regReport));
  regReport *tmp;
  /* NOTE:LDIM is the global */
  int varInds[LDIM];/* indices of the same variables from different sub-systems */
  int sl_am_len[LDIM];/* size of the slope amplitudes arrays for each sub-system */
  int regime[LDIM];/* numerical indicator of the regime in each sub-system */
  trajSlope *slope[LDIM];/* structure holding the slope characteristics */

  /* Get indices of the corresponding variables in each sub-system. Note we take the
     first variable used for the graphics output. This can be changed. Also note that
     data_name is the global variable. */
  get_varInd(graph.yInd[0],varInds);
  for(j=0;j<LDIM;j++){
    if(LDIM>1){
      if(traj_verb == 1)
	printf("***Checking sub-system %d(%s):\n",j+1,var_name[varInds[j]]);
    }
    else{
      if(traj_verb == 1)
	printf("***Checking variable %s\n",var_name[varInds[j]]);
    }
    slope[j] = get_dynamics_1d(data_name,varInds[j],slope[j],&sl_am_len[j],&regime[j]);
    /* We escape the further process, since single system dynamics check failed */
    if(slope[j] == NULL){
      printf("Error: CANNOT finish the dynamics test.\n");
      return NULL;
    }
  }
  tmp = get_dynamics_nd(slope,sl_am_len,regime);
  report[report_size] = *tmp;
  free_report(tmp);

  return report;
}

int find_sys_lag(trajSlope *slope[],int sam_len[],int *lag,
		double *b_gain, double *a_gain)
{/* This function finds the LAG between the subsystems. The lag is nonzero in case of
    the out-of-phase regime, i.e. non-stationary regime. Additionally, if the lag is
    not found for an oscillatory regime, the regime is considred in-homogeneous. */
  /* NOTE: LDIM is the number of sub-systems. It is the global constant. */
  /* This function should not be called upon LDIM = 1 case. */
  int i,j,idx=0,k;
  /* We have LDIM physical systems. Each system's last slope amplitude is compared
     with all other systems' slope amplitudes: starting from the last and proceeding
     to the first until the similar (as determined by the relative error eps_traj_rel)
     amplitudes are found. Additionally, the compared slopes must have the same
     direction (ascending or descending). Thus the array of lags is formed (by all
     mutual comparisons). */
  for(i=0;i<LDIM;i++){
    for(j=i+1;j<LDIM;j++){
      for(k=sam_len[j]-1;k>0;k--){
	//printf("%G vs %G\n",slope[i]->ampl[sam_len[i]-1],slope[j]->ampl[k]);
	if((fabs(slope[i]->ampl[sam_len[i]-1] - slope[j]->ampl[k]) <
	    (eps_traj_rel*slope[i]->ampl[sam_len[i]-1])) &&	\
	   (slope[i]->ascend[sam_len[i]-1] == slope[j]->ascend[k])
	   ){
	  /* Amplitudes are similar, form the base gain */
	  if(slope[i]->base[sam_len[i]-1] > slope[j]->base[k])
	    b_gain[idx] = slope[i]->base[sam_len[i]-1] / slope[j]->base[k];
	  else
	    b_gain[idx] = slope[j]->base[k] / slope[i]->base[sam_len[i]-1];
	  if(slope[i]->ampl[sam_len[i]-1] > slope[j]->ampl[k])
	    a_gain[idx] = slope[i]->ampl[sam_len[i]-1] / slope[j]->ampl[k];
	  else
	    a_gain[idx] = slope[j]->ampl[k] / slope[i]->ampl[sam_len[i]-1];
	  break;
	}
      }
      if(k == 0){
	if(traj_verb == 2)
	  printf("Could not determine the lag[%d] => IH?\n",idx);
	lag[idx] = -1;/* could not determine */
      }
      else{
	/* After break k will hold the lag, if not -1 */
	lag[idx] = (sam_len[j]-1) - k;
      }
      if(traj_verb == 2)
	printf("#%d) Systems lag = %d",idx,lag[idx]);
      if(lag[idx] >= 0){
	if(traj_verb == 2)
	  printf(", bg = %lf, ag = %lf\n",b_gain[idx],a_gain[idx]);
      }
      else{
	if(traj_verb == 2)
	  printf("\n");
      }
      idx++;
    }
  }

  return 0;
}

void free_report(regReport *report)
{
  free(report);
}

regReport *get_dynamics_nd(trajSlope *slope[],int sam_len[],int regime[])
{
  /* The function determines the overal multi-dimensional dynamical regime of the
     system under study OR summarizes the functioning of the 1-D system. The input
     slope is gotten freed here. */
  int i,j,k,m;
  int ssize = LDIM;/* trajSlope array size */
  /* We need to check:
     1. the variable levels in different sub-systems (homogeneity test)
     2. and the time moments of the slopes (phase test)
     However, for the SS dynamics we do not need the phase test.
   */
  double tmp = 0.0;/* some temporary variable */
  double b_gain = 0.0;/* maximum base gain */
  double a_gain = 0.0;/* maximum amplitude gain */
  double p_gain = 0.0;/* maximum phase gain */
  double p_shift = 0.0;/* maximum phase shift */
  double period[ssize];/* period of the sub-systems */
  regReport *report;/* report of the collective regime */
  report = (regReport *)malloc(sizeof(regReport));
  /* ***************************** */
  int *lag = (int *)malloc((LDIM*(LDIM-1)/2)*sizeof(int));
  double *bg = (double *)malloc((LDIM*(LDIM-1)/2)*sizeof(double));
  double *ag = (double *)malloc((LDIM*(LDIM-1)/2)*sizeof(double));
  double *psh = (double *)malloc((LDIM*(LDIM-1)/2)*sizeof(double));
  /*** OVERALL DYNAMICS: process the regime identificators ***/
  /* Sum up all numerical identificators of the dynamics in all sub-systems */
  for(i=0;i<ssize;i++)
    tmp += (double)regime[i];
  /* If the sum is negative, regime is undetermined */
  if(tmp < 0){
    report->regime = -1;/* negative regime: undetermined */
  }
  /* If the sum equals 0, it is purely SS */
  else if((double)fabs(tmp/ssize) < 0.00001){/* Practically zero! */
    report->regime = 0;/* steady state regime */
  }
  /* If it is not zero AND equals to the first regime numeric value,
     the regime is pure and OS (since not SS) */
  else if((double)fabs((tmp/ssize) - regime[0]) < 0.00001){
    report->regime = (int)(tmp/ssize);
  }
  /* Otherwise, the regime is mixed, which is basically an indicator of the strange
     behavior. */
  else{
    report->regime = -2;/* mixed regime */
  }

  /* *** Base and Amplitude gains as well as S-LAG *** */
  if(report->regime > 0){/* Non-SS regime */
    /* Find out the System-lag and the corresponding gains */
    find_sys_lag(slope,sam_len,lag,bg,ag);
    /* Phase shift */
    /* Mutual comparisons */
    k = 0;
    for(i=0;i<ssize;i++){
      for(j=i+1;j<ssize;j++){
	if(lag[k] >= 0){/* If the lag exists */
	  /* we use t0 difference for the phase shift */
	  if(slope[i]->t0[sam_len[i]-1] > slope[j]->t0[sam_len[j]-1-lag[k]])
	    psh[k] = slope[i]->t0[sam_len[i]-1] - slope[j]->t0[sam_len[j]-1-lag[k]];
	  else
	    psh[k] = slope[j]->t0[sam_len[j]-1-lag[k]] - slope[i]->t0[sam_len[i]-1];
	}
	else{/* Lag does not exist => IH oscillations. No shift can be determined. */
	  /* Negative number as a marker for IH oscillations, since:
	     slope[i]->t0[sam_len[i]-1] > slope[i]->t0[sam_len[i]-2]
	     Make it approximately equal to the period.
	   */
	  psh[k] = - (slope[i]->t0[sam_len[i]-1] - slope[i]->t0[sam_len[i]-1-regime[i]*2]);
	  /* If the lag does not exist == > IH OS. We must calculate the base and
	     amplitude gains here: the SAME code as for any other regime, including SS. */
	  /* Base gains */
	  /* The gain is always >= 1 for we divide the larger number by smaller one. */
	  if(slope[i]->base[sam_len[i]-1] > slope[j]->base[sam_len[j]-1])
	    bg[k] = slope[i]->base[sam_len[i]-1]/slope[j]->base[sam_len[j]-1];
	  else
	    bg[k] = slope[j]->base[sam_len[j]-1]/slope[i]->base[sam_len[i]-1];
	  /* Maximum Amplitude gain */
	  if(slope[i]->ampl[sam_len[i]-1] > slope[j]->ampl[sam_len[j]-1])
	    ag[k] = slope[i]->ampl[sam_len[i]-1]/slope[j]->ampl[sam_len[j]-1];
	  else
	    ag[k] = slope[j]->ampl[sam_len[j]-1]/slope[i]->ampl[sam_len[i]-1];
	}
	k++;
      }
      /* Period estimation: period is formed by the LAST (#SUBPER)*2 slopes. #SUBPER ==
	 regime. */
      period[i] = slope[i]->t0[sam_len[i]-1] - slope[i]->t0[sam_len[i]-1-regime[i]*2];
      if(traj_verb == 2)
	printf("period %d: %G",i+1,period[i]);
      if(regime[i] > 1){
	for(m=0;m<regime[i];m++){
	  if(m == 0){
	    if(traj_verb == 2)
	      printf(" = %G",slope[i]->t0[sam_len[i]-1] -
		     slope[i]->t0[sam_len[i]-1-2]);
	  }
	  else{
	    if(traj_verb == 2)
	      printf(" + %G",
		     slope[i]->t0[sam_len[i]-1-m*2]-slope[i]->t0[sam_len[i]-1-(m+1)*2]);
	  }
	}
      }
      if(traj_verb == 2)
	printf("\n");
    }
  }
  else{/* for ANY other regime */
    k = 0;
    for(i=0;i<ssize;i++){/* i<ssize-1(for backup) */
      for(j=i+1;j<ssize;j++){/* j=i+1(for backup) */
	/* Base gains */
	/* The gain is always >= 1 for we divide the larger number by smaller one. */
	if(slope[i]->base[sam_len[i]-1] > slope[j]->base[sam_len[j]-1])
	  bg[k] = slope[i]->base[sam_len[i]-1]/slope[j]->base[sam_len[j]-1];
	else
	  bg[k] = slope[j]->base[sam_len[j]-1]/slope[i]->base[sam_len[i]-1];
	/* Maximum Amplitude gain */
	if(slope[i]->ampl[sam_len[i]-1] > slope[j]->ampl[sam_len[j]-1])
	  ag[k] = slope[i]->ampl[sam_len[i]-1]/slope[j]->ampl[sam_len[j]-1];
	else
	  ag[k] = slope[j]->ampl[sam_len[j]-1]/slope[i]->ampl[sam_len[i]-1];
	k++;
      }
      /* Free the slopes: this function is called from get_sys_dynamics() where slope
	 is initialized as array of pointers, i.e. not malloc-ed, thus no need in
	 free-ing. */
      //free_slope(slope[i]);
    }
  }

  /* *********************** */
  /***** FORM THE REPORT *****/
  if(report->regime <= 0){/* SS,mixed or undetermined */
    report->ph_shift = 0.0;/* no phase shift for SS */
    report->period = 0.0;/* no period for SS */
  }
  if(report->regime > 0){/* Non-SS positive regime */
    report->ph_shift =
      (gsl_stats_max(psh,1,LDIM*(LDIM-1)/2))/gsl_stats_max(period,1,LDIM);
    report->period = period[0];
    k = 0;
    for(i=0;i<ssize;i++){
      for(j=i+1;j<ssize;j++){
	if(traj_verb == 2)
	  printf("phase shift %d: %G\n",k+1,psh[k]/period[i]);
	k++;
      }
    }
  }
  /*** HOMOGENEITY ***/
  /* base OR amplitude is different => inhomogeneous */
  if((gsl_stats_max(bg,1,LDIM*(LDIM-1)/2)) > (1+eps_traj_rel) ||	\
     (gsl_stats_max(ag,1,LDIM*(LDIM-1)/2)) > (1+eps_traj_rel)){
    report->homog = 0;
    report->b_gain = gsl_stats_max(bg,1,LDIM*(LDIM-1)/2);
    report->a_gain = gsl_stats_max(ag,1,LDIM*(LDIM-1)/2);;
  }
  else{
    report->homog = 1;
    report->b_gain = gsl_stats_max(bg,1,LDIM*(LDIM-1)/2);
    report->a_gain = gsl_stats_max(ag,1,LDIM*(LDIM-1)/2);
  }

  /* Free the system lag, bg and ag. */
  free(lag);
  free(bg);
  free(ag);
  free(psh);

  return report;
}

trajSlope *get_dynamics_1d(char const *fname, int varInd, trajSlope *slope,
			   int *samp_len,int *regime)
{
  /* The function computes the dynamical regime of the 1-D dynamical system. */
  int i,j,k,m,n;
  int frame_size;
  double **frame;
  trajPeak *peak;
  int sl_ampl_size;
  FILE *in = fopen(fname,"r");
  if(in == NULL){
    fprintf(stderr,"Error: could not open file `%s'\n",fname);
    return NULL;
  }
  else{
    /* We load the first deterministic frame */
    frame = load_frame(in,frame,&frame_size);
    if(frame == NULL){
      fprintf(stderr,"Error: frame was not loaded\n");
      return NULL;
    }
    /* ********************************************* */
    /* regime = -1 ===> not known */
    /* regime = 0 ===> steady state(SS) */
    /* regime = 1...n ===> period-n oscillations(OS) */
    /* ********************************************* */
    /* Make regime undetermined from the very beginning */
    *regime = -1;
    /* Get info on peaks and troughs */
    peak = peak_trough2(frame,frame_size,varInd,peak);
    fclose(in);/* We dont need the file stream as well */
    sl_ampl_size = peak->size-1;/* The size of sl_ampl is known beforehand */
    if(sl_ampl_size == -1){/* There is no single peak/trough in the system */
      /* If the user trusts his/her time scale, then it is SS */
      *regime = 0;
      if(traj_verb == 1)
	printf("SS(%d): no peak/trough was found.\n",*regime);
      if(traj_verb == 2)
	printf("Entering checking frame...\n");
      i = ss_by_frame(frame,frame_size,varInd);
      if(i){
	printf("(Could be not enough data.)\n");
      }
      /* Simple rules to get the regime determined: this case can happen when all the
	 values in the frame are equal and in a SS.*/
      slope = (trajSlope *)malloc(sizeof(trajSlope));
      slope->ampl = (double *)malloc(sizeof(double));
      slope->t0 = (double *)malloc(sizeof(double));
      slope->base = (double *)malloc(sizeof(double));
      slope->ampl[0] = 0.0;
      slope->t0[0] = frame[0][frame_size-1];
      slope->base[0] = frame[varInd+1][frame_size-1];
      *samp_len = 1;
      //return NULL;
      free_frame(frame);
      return slope;
    }
    free_frame(frame);/* We dont need the frame anymore */
    /* Calculate slope characteristics */
    slope = slope_ampl(*peak, slope);
    free_peak(peak);/* We dont need peak struct anymore */
    /* Check the sanity of the amplitudes */
    //slope = check_slope_ampl(slope,&sl_ampl_size);
    /* Do we have enough amplitudes to proceed after the checking? */
    if(sl_ampl_size < 2){
      /* We need at least 2 amplitudes. However, 2 is not enough as well. */
      fprintf(stderr,"Warning: only %d amplitude(s), get more data.\n",sl_ampl_size);
      /* HERE THE TRAJECTORY-ASSISTANT MIGHT COME INTO PLAY */
      if(sl_ampl_size > 0){
	*regime = 0;
	*samp_len = 1;
	return slope;
      }
      return NULL;
    }
    if(sl_ampl_size < 4){
      fprintf(stdout,"Warning: too few amplitudes, analysis might be hampered.\n");
      /* HERE THE TRAJECTORY-ASSISTANT MIGHT COME INTO PLAY */
    }
    /* Copy the size of the amplitudes array to the returned value */
    *samp_len = sl_ampl_size;
    /* ************************* */
    /* DYNAMICAL REGIME CHECKING */
    /* ************************* */
    /* Check the last slope amplitude, compare it with the tolerance, which is
       system-wide and equal to eps_traj_abs. */
    /* double mx_ampl = gsl_stats_max(slope->ampl,1,sl_ampl_size); */
    /* double att_app_rate[sl_ampl_size-1]; */
    /* Form the rate of approaching to an attractor. */
    /* printf("Attractor approaching rate:\n"); */
    /* for(i=0;i<sl_ampl_size-1;i++){ */
    /*   /\* att_app_rate[i] = 100*(slope->ampl[sl_ampl_size-1] - slope->ampl[i]) / *\/ */
    /*   /\* 	slope->ampl[sl_ampl_size-1]; *\/ */
    /*   att_app_rate[i] = 100*(slope->ampl[i+1] - slope->ampl[i])/	\ */
    /* 	(slope->ampl[sl_ampl_size-1]*(slope->t0[i+1] - slope->t0[i]));  */
    /*   printf("%G ",att_app_rate[i]); */
    /* } */
    /* printf("\n"); */
    /* ====================================================================== */
    /********************** START THE DYNAMICS CHECKING *********************/
    if(slope->ampl[sl_ampl_size-1] < eps_traj_abs){
      /********************************* SS 1st PATH *****************************/
      /* The first very simple SS-checking rule. If it fails proceed to the more
	 complicated checking. */
      *regime = 0;
      if(traj_verb == 1)
	fprintf(stdout,"SS(%d): abs.tol = %G.\n",*regime,eps_traj_abs);
      /******************************* END: SS 1st PATH ***************************/
    }
    if(*regime == -1){/* If regime is still undetermined */
      /********************************* OS 1st PATH ******************************/
      /* 1. Check amplitudes and bases of every other slope, since 2 slopes form the
	 complete cycle. Based on the comparison results, propose the number of
	 sub-periods j. */
      m = 0;/* Initial value for the #subperiods */
      /* Iterate by comparing the last slope with the preceding ones. */
      for (k = sl_ampl_size-1; k > 0; k -= 1){
	/* Set the #subperiods to zero */
	j = 0;
	for (i = k-2; i > 0; i -= 2){
	  j++;/* Increase the #subperiods */
	  /* the last slope (ind=sl_ampl_size-1) is compared with every other slope
	     behind it starting with a slope with ind=sl_ampl_size-3. Jump of the
	     i-variable is 2 backwards. */
	  /* Two conditions must be satisfied: relative equality of amplitudes and
	     bases. Relatively equal = equal within the boundaries set by the relative
	     error. Note: the base comparison takes place against the amplitude. */
	  if( (fabs(slope->ampl[i] - slope->ampl[k]) <	\
	       (eps_traj_rel*slope->ampl[k])) &&		\
	      (fabs(slope->base[i] - slope->base[k]) <	\
	       (eps_traj_rel*slope->ampl[k]))){
	    /* Once we have found the condition, escape the loop. So the comparable
	       slopes are equal within the given relative error. Keep the j value
	       before. */
	    if(k == (sl_ampl_size-1))
	      /* For the first comparison hold the j value in m */
	      m = j;
	    break;

	  }
	  /* But if we reached the end of the loop (without breaking), no periodicity
	     is found and j must get zero value to invoke the further analysis of the
	     trajectory (see below OS 2nd PATH). */
	  if(i == 1 || i == 2)/* End of the loop is reached */
	    j = 0;/* Make the j zero to indicate no reasonable #subperiods is found */
	}
	if(j == 0)/* Escape further comparisons since the 1st one failed */
	  break;
	if(m != j)/* If j is not equal to the previous value held in m, then break */
	  break;
      }
      if(m){/* If the #subperiods exists, i.e. some positive ineteger. */
	/* Passed the 1st check: deterministic cyclic trajectory with period j. */
	*regime = m;
	/* The total number of comparisons possible is equal to sl_ampl_size-1 - 2.
	   sl_ampl_size-1 - k comparisons succeeded. */
	if(traj_verb == 1){
	  printf("# subperiods is %d: %d of %d (%G%%) comparisons support the conclusion.\n",
		 m,sl_ampl_size-1-k,sl_ampl_size-2,
		 100*((double)sl_ampl_size-1-k)/(sl_ampl_size-2));
	  printf("OS(%d): period-%d,rel.tol=%G.\n",*regime,m,eps_traj_rel);
	}
      }
      /******************************* END: OS 1st PATH ****************************/
    }
    if(*regime == -1){/* If regime is still undetermined */
      /********************************* SS 2nd PATH *****************************/
      /* Try to find the successful descending slope amplitudes. NOTE: this could
	  indicate the approaching to both SS and OS, however, the numerical entities
	  calculated must give a guess to the user. */
      if((slope->ampl[sl_ampl_size-2] - slope->ampl[sl_ampl_size-1]) > \
	 (eps_traj_rel * slope->ampl[sl_ampl_size-2])){
	j = 1;/* counter for the number of decreasing slope amplitudes */
	for(i=sl_ampl_size-2;i>0;i--){
	  if((slope->ampl[i-1] - slope->ampl[i]) >
	     (eps_traj_rel * slope->ampl[i-1])){
	    j++;
	  }
	  else
	    /* Break here since we are interested in monotonical decrease without
	       interruptions. */
	    break;
	}
	/* min 50% of the slope amplitudes do decrease */
	if(((double)100*(j+1)/(sl_ampl_size) > 50)){
	  *regime = 0;
	  if(traj_verb == 1){
	    printf("Last ampl. (%G) / Abs. tol. (%G) = %G\n",slope->ampl[sl_ampl_size-1],eps_traj_abs,slope->ampl[sl_ampl_size-1]/eps_traj_abs);
	    printf("SS(%d): %d of %d(%G%%) slopes demonstrate monotonical decrease\n",*regime,
		   j+1,sl_ampl_size,(double)100*(j+1)/(sl_ampl_size));
	  }
	  /* Total decrease in amplitudes > 50% ==> SS, if <= 50% warning message */
	  if((100*(1-slope->ampl[sl_ampl_size-1]/slope->ampl[sl_ampl_size-1-j]) <= \
	      50)){
	    if(traj_verb == 1)
	      printf("Warning: total amplitude decay is %G%%\n",
		     100*(1-slope->ampl[sl_ampl_size-1]/slope->ampl[sl_ampl_size-1-j]));
	  }
	  if(traj_verb == 1)
	    printf("SS(%d): total decrease=%G%%,rel.tol=%G\n",*regime,
		   100*(1-slope->ampl[sl_ampl_size-1]/slope->ampl[sl_ampl_size-1-j]),
		   eps_traj_rel);
	}
	else{
	  /* What if this is OS regime, how can we reach the OS-checking from here???
	  */
	  if(traj_verb == 1){
	    printf("SS(%d): failed to pass 50%% criterium for monotonical decrease.\n",*regime);
	    printf("SS(%d): %d of %d(%G%%) slopes demonstrate monotonical decrease\n",*regime,
		   j+1,sl_ampl_size,(double)100*(j+1)/(sl_ampl_size));
	  }
	}
      }
      /******************************* END: SS 2nd PATH ***************************/
    }
    if(*regime == -1){
      /********************************* OS 2nd PATH ******************************/
      /* m = 0 means amplitudes are not equal at all. Try to find ascending
	 amplitudes pointing to the cycle that was not reached by the integration. By
	 the analogy with not reached SS and descending slope amplitudes.*/
      if(m == 0){
	/* Do mutual comparisons before giving up saying that cannot recognize the
	   regime. This is for finding monotonical increase of the amplitudes. */
	if((slope->ampl[sl_ampl_size-1] - slope->ampl[sl_ampl_size-2]) > \
	   (eps_traj_rel*slope->ampl[sl_ampl_size-2])){
	  j++;/* First comparison succeeded. */
	  for(i=sl_ampl_size-2;i>0;i--){
	    if((slope->ampl[i] - slope->ampl[i-1]) > (eps_traj_rel*slope->ampl[i-1])){
	      j++;
	    }
	    else/* Non monotonical increase of the amplitudes. */
	      break;
	  }
	  /* min 50% of slopes increase ==> OS. */
	  if(((double)(100*(j+1)/sl_ampl_size) > 50)){
	    *regime = 1;
	    if(traj_verb == 1){
	      printf("%d of %d(%G%%) slopes demonstrate monotonical increase\n",
		     j+1,sl_ampl_size,(double)100*(j+1)/(sl_ampl_size));
	      printf("OS(%d): total increase=%G%%,rel.tol=%G\n",*regime,
		     100*slope->ampl[sl_ampl_size-1]/slope->ampl[sl_ampl_size-1-j],
		     eps_traj_rel);
	    }
	  }
	  else{
	    if(traj_verb == 1){
	      printf("OS(%d): failed to pass 50%% criterium for monotonical increase.\n",*regime);
	      printf("OS(%d): %d of %d(%G%%) slopes demonstrate monotonical increase\n",
		     *regime,j+1,sl_ampl_size,(double)100*(j+1)/(sl_ampl_size));
	    }
	  }
	}
	/******************************* END: OS 2nd PATH ****************************/
      }
      /* else{ */
      /*   /\* Passed the 1st check: deterministic cyclic trajectory with period j. *\/ */
      /*   *regime = m; */
      /*   /\* The total number of comparisons possible is equal to sl_ampl_size-1 - 2. */
      /*     sl_ampl_size-1 - k comparisons succeeded. *\/ */
      /*   printf("Lag is %d: %d of %d (%G%%) comparisons support the conclusion.\n",m, */
      /* 	 sl_ampl_size-1-k,sl_ampl_size-2,100*((double)sl_ampl_size-1-k)/(sl_ampl_size-2)); */
      /*   printf("OS(%d): period-%d,rel.tol=%G.\n",*regime,m,eps_traj_rel); */
      /* } */
    }
    if(*regime == -1){
      printf("Trajectory System Error: could not recognize regime.\n");
    }
  }
  return slope;
}

int ss_by_frame(double **frame, int frame_size, int varInd)
{
  /* This function checks how long back the trajectory has equal values for
     varInd. This function can be entered only when there are no peaks/troughs found
     in the trajectory. */
  int i;
  i = frame_size - 1;/* last value in the trajectory */
  while((frame[varInd+1][i] == frame[varInd+1][i-1]) && (i > 0)){
    i--;
  }
  /* Remember frame[0][i] is the time/independent-variable of the trajectory. */
  if(traj_verb == 1)
    printf("SS(0): equal for %G%% of time.\n",
	   100*(frame[0][frame_size-1] - frame[0][i])/frame[0][frame_size-1]);
  /* Use 50% criterion over time.*/
  if( (frame[0][frame_size-1] - frame[0][i]) >= (0.5*frame[0][frame_size-1]) ){
    /* SS over more or equal to 50% of the time */
    return 0;
  }
  else{
    /* Otherwise not SS */
    return 1;
  }
}

trajSlope *check_slope_ampl(trajSlope *slope, int *sampl_size)
{
  int i = 0;/* Main looping var, holding the amplitude resumption index */
  int j = 0;/* Helping var, holding the amplitude drop index */
  int k;/* Looping var for cutting the non-relevant values */
  int num_removed = 0;/* Number of removed amplitudes */
  int tot_num_ampl = *sampl_size;/* Total number of amplitudes before checking */
  int num_drop = 0;/* Number of amplitude drops. Needed for single drop trajectories. */
  /* Find the max of all amplitudes and take 1% from it. This will be the global
     error level. Either amplitude of the comparison pair should be larger than
     that. */
  double eps = 0.01*gsl_stats_max(slope->ampl,1,(*sampl_size));
  while(i<(*sampl_size-1)){
    if((slope->ampl[i] > eps) || (slope->ampl[i+1] > eps)){
      if((20*slope->ampl[i]) < slope->ampl[i+1]){/* Amplitude boom first */
	/* First <i+1> values are removed, <i+1> holds the postion of the amplitude
	   resumption. */
	num_removed += i+1;
	/* Shift the array by <i+1> position leftwards(to the beginning). */
	for(k=i+1;k<*sampl_size;k++){
	  slope->ampl[k-i-1] = slope->ampl[k];
	  slope->t0[k-i-1] = slope->t0[k];
	}
	*sampl_size = *sampl_size - i - 1;/* Reduce the size and reallocate */
	slope->ampl = (double *)realloc(slope->ampl,(*sampl_size)*sizeof(double));
	slope->t0 = (double *)realloc(slope->t0,(*sampl_size)*sizeof(double));
	i = 0;/* Reset <i> to the beginning, since we removed first elements. */
	continue;/* Do not need to go further, it implied i++ we don't need. */
      }
      if( (0.05*slope->ampl[i]) > slope->ampl[i+1] ){/* Amplitude drop */
	num_drop++;
	j = i;/* We copy the position before the large drop in the amplitude */
	i++;/* Move to the next value */
	/* Move next until we get a large boom or end of the array */
	while(i<(*sampl_size-1)){
	  if((slope->ampl[i] > eps) || (slope->ampl[i+1] > eps)){
	    if((20*slope->ampl[i]) < slope->ampl[i+1])
	      break;
	  }
	  i++;
	}
	if(i == (*sampl_size - 1)){
	  /* No large boom appeared => No removal is needed. */
	  continue;
	}
	/* After the above <i> holds the index(0-based) of the value which is the last
	   of the amplitude array OR the last before the normal size amplitudes
	   resumption. Everything between <j> and <i>(excluding <j>, including <i>)
	   must be cut and removed. */
	num_removed += i-j;
	*sampl_size = *sampl_size - (i-j);/* Reduce the size of the array */
	/* Cut the non-relevant values: elements are preserved until <j>, after <j> we
	   need to do the shift of the values. */
	for(k=j+1;k<*sampl_size;k++){
	  slope->ampl[k] = slope->ampl[k+(i-j)];
	  slope->t0[k] = slope->t0[k+(i-j)];
	}
	/* Reallocate memory for the reduced array */
	slope->ampl = (double *)realloc(slope->ampl,(*sampl_size)*sizeof(double));
	slope->t0 = (double *)realloc(slope->t0,(*sampl_size)*sizeof(double));
	i = j;/* To reset the starting point after the reallocation */
      }
    }
    i++;
  }
  if(num_removed)
    if(traj_verb == 1)
      printf("Slope check: removed %d/%d slope(s).\n",num_removed,tot_num_ampl);
  //gplot_slope_ampl(slope,*sampl_size);
  
  return slope;
}

void free_slope(trajSlope *slope)
{
  free(slope->ampl);
  free(slope->t0);
  free(slope->base);
  free(slope->ascend);
  free(slope);
}

trajSlope *slope_ampl(const trajPeak peak,trajSlope *slope)
{
  int i;
  slope = (trajSlope *)malloc(sizeof(trajSlope));
  slope->ampl = (double *)malloc((peak.size-1)*sizeof(double));
  slope->t0 = (double *)malloc((peak.size-1)*sizeof(double));
  slope->base = (double *)malloc((peak.size-1)*sizeof(double));
  slope->ascend = (int *)malloc((peak.size-1)*sizeof(int));

  /* FILE *out; */
  /* /\* if((fopen("slopeAmpl.dat","r")) != NULL) *\/ */
  /* /\*   out = fopen("slopeAmpl.dat","a"); *\/ */
  /* /\* else *\/ */
  /*   out = fopen("slopeAmpl.dat","w"); */
  
  if(slope == NULL){
    fprintf(stderr,"Error: could not allocate trajSlope object\n");
    return NULL;
  }
  for(i=0;i<peak.size-1;i++){/* Iterate over the found peaks */
    /* Define the slope as the difference between successive peaks */
    slope->ampl[i] = peak.x[i+1] - peak.x[i];
    /* Define the first point as the slope's base */
    slope->base[i] = peak.x[i];
    /* If the slope amplitude is negative, it means the descending slope */
    if(slope->ampl[i] < 0){
      /* Make the amplitude positive */
      slope->ampl[i] = -slope->ampl[i];
      /* And re-define the slope's base */
      slope->base[i] = peak.x[i+1];
      /* Indicate the descending slope */
      slope->ascend[i] = 0;
    }
    else/* otherwise the slope is ascending */
      slope->ascend[i] = 1;
    /* Take the time at which the slope has started */
    slope->t0[i] = peak.t[i];
    /* For the output */
    if(traj_verb == 2)
      fprintf(stdout,"%G %G %G\n",slope->t0[i],slope->ampl[i],slope->base[i]);
  }
  if(traj_verb == 2)
    fprintf(stdout,"\n\n");
  /* fclose(out); */
  /* gnuplot_cmd(plot_handle,"set term x11 1\n"); */
  /* gnuplot_cmd(plot_handle,"set xlabel \"Time\"\n"); */
  /* gnuplot_cmd(plot_handle,"set ylabel \"Slope Amplitude\"\n"); */
  /* gnuplot_cmd(plot_handle,"plot 'slopeAmpl.dat' ps 3 not\n"); */
  return slope;
}

void gplot_slope_ampl(trajSlope *sampl,int const sampl_size)
{
  int i;
  FILE *out = fopen("slopeAmpl1.dat","w");
  for(i=0;i<sampl_size;i++)
    fprintf(out,"%G %G\n",sampl->t0[i],sampl->ampl[i]);
  fclose(out);
  gnuplot_cmd(plot_handle,"set term x11 1\n");
  gnuplot_cmd(plot_handle,"replot 'slopeAmpl1.dat' ps 3 not\n");
}

int steady_state(double *ss_factor)
{/*This function determines Steady State regime existence being based
   on calculating a difference between variables in the computed
   trajectory. This function returns back how many points are passed
   the criterion of SS.*/

  int i,j;
  double E[DIM];/*It is the difference(error)*/
  double eps=eps_abs_am;/*We use Amplitude absolute error here*/
  for(j=0;j<DIM;j++){
    i=n_steps-1;/*Indexing in C starts from zero=0*/
    while(i>0){/*We moving backward*/
      E[j]=xs[i][j]-xs[i-1][j];if(E[j]<0) E[j]=-E[j];/*Forming the error*/
      if(E[j]>eps)
	{printf("Step of returning %d\n",i);break;}
      i--;
    }
    ss_factor[j]=(ts[n_steps-1]-ts[i])/ts[n_steps-1];
  }
    
  return 0;
}


int attr_conv(double **dat, int const size)
{
  /* Attractor convergence function */
  /* dat[][] is the 2d array of data points and size is the size of the array */
  /* This function takes the last point of the data and compares it with all previous
     points trying to find the closest vicinity in the phase space. If it finds the
     vicinity with a given accuracy it reports success(0) in getting the
     attractor. Dimension of the data is DIM, to be defined elsewhere. */
  FILE *out;
  double dist;/* the Euclidean distance between points*/
  int i,j,k;
  i = size-1;
  while(i >= 0){
    k = size-1;
    while(k >= 0){
      if(k == i) {
	k--;
	continue;
      }
      dist = 0.0;
      for(j=0;j<DIM;j++)/* Compute the distance */
	dist += (dat[k][j]-dat[i][j])*(dat[k][j]-dat[i][j]);/*avoiding pow function*/
      dist = sqrt(dist);
      if(dist < 0.01) {/* How to determine this error??? */
	printf("k/i/size=%d/%d/%d\n",k,i,size);
	return 0;
      }/* we are at the attractor */
      k--;
    }
    i--;
  }
  fprintf(stdout,"Did not get attractor.\n");
  return 1;/* we are NOT at the attractor, needing more integration or chaos(?) */
}

int period_sort(double * per_ratio)
{
  // Criterion of convergence to limit cycle attractor.
  int i,j,k;
  int count[MAXDIM];
  float period;
  double E[MAX_N_PERIODS-1][MAXDIM]; 
  if(per_method == 0){
    for(i=0; i<DIM; i++){
      count[i]=0;
      k=n_subperiods[i];
      if(n_subperiods[i] == 0) k=k+1;
      for(j=0; j<n_isec[i]-1-k; j++){
	E[j][i] = (double)per[j][i] - (double)per[j+k][i];
	if(E[j][i] < 0)
	  E[j][i] = -E[j][i];
      }
      //Set to unity all elements of per_ratio(this means no convergence to L.C.)
      *(per_ratio+i)=1.0;
    }
    for(i=0;i<DIM; i++){
      k=n_subperiods[i];
      if(n_subperiods[i] == 0) k=k+1;
      for(j=0; j<n_isec[i]-1-k; j++){
	period = (per[j][i] + per[j+k][i]) / 2;
	if(E[j][i] > D(eps_abs_per,eps_rel_per,period))
	  count[i]++;
      }
    }
    for(i=0;i<DIM;i++){
      *(per_ratio+i) = (double)count[i]/((double)n_isec[i] - 1.0 - \
    (double)k);
      if( *(per_ratio+i) < 0 )
	*(per_ratio+i) = 1.0;
      if(*(per_ratio+i) > eps_per_ratio)
	fprintf(stdout,"Did not get attractor: %lf for var #%d(%s)\n",*(per_ratio+i),i+1,var_name[i]);
    }
  }
  if(per_method == 1){
    for(i=0; i<DIM; i++){
      count[i]=0;
      k=n_subperiods[i];
      if(n_subperiods[i] == 0) k=k+1;
      for(j=0; j<n_peaks[i]-1-k; j++){
	E[j][i] = (double)per[j][i] - (double)per[j+k][i];
	if(E[j][i] < 0)
	  E[j][i] = -E[j][i];
      }
      //Set to unity all elements of per_ratio(this means no convergence to L.C.
      *(per_ratio+i)=1.0;
    }
    for(i=0;i<DIM; i++){
      k=n_subperiods[i];
      if(n_subperiods[i] == 0) k=k+1;
      for(j=0; j<n_peaks[i]-1-k; j++){
	period = (per[j][i] + per[j+k][i]) / 2;
	if(E[j][i] > D(eps_abs_per,eps_rel_per,period))
	  count[i]+=1;
      }
    }
    for(i=0;i<DIM;i++){
      *(per_ratio+i) = (double)count[i]/((double)n_peaks[i] - 1.0 - \
    (double)k);
      if( *(per_ratio+i) < 0 )
	*(per_ratio+i) = 1.0;
      if(*(per_ratio+i) > eps_per_ratio)
	fprintf(stdout,"Did not get attractor: %lf for var #%d(%s)\n",*(per_ratio+i),i+1,var_name[i]);
    }
  }
  return 0;
}

int subper(int * spectrum_per,int * n_subperiods)
{/*This functions computes subperiods of the oscillations*/
  int i,j;

  int n_false[MAXDIM],n_true[MAXDIM];
  double E[MAX_N_PERIODS-1][MAXDIM]; 
  float period;
  if(per_method == 0){
    for(i=0; i<DIM; i++){
      n_false[i]=0;
      n_true[i]=0;
      for(j=0; j<n_isec[i]-2; j++){
	E[j][i] = (double)per[0][i] - (double)per[j+1][i];
	if(E[j][i] < 0)
	  E[j][i] = -E[j][i];
      }
      // Set to zeros ALL spectrum_per's elements
      for(j=0;j<MAX_N_ISEC-2;j++){
	*(spectrum_per+j*MAXDIM+i)=FALSE;
      }
      // Set to zeros ALL n_subperiods's elements
      *(n_subperiods+i)=0;
    }
    for(i=0; i<DIM; i++){
      period = per[0][i];
      for(j=0; j<n_isec[i]-2; j++){
	if(E[j][i] > D(eps_abs_per,eps_rel_per,period)){
	  *(spectrum_per+j*MAXDIM+i) = FALSE;
	  n_false[i]+=1;}
	if(E[j][i] < D(eps_abs_per,eps_rel_per,period)){
	  *(spectrum_per+j*MAXDIM+i) = TRUE;
	  if(n_true[i] == 0)
	    *(n_subperiods+i) = j+1;
	  n_true[i]+=1;}
      }
    }
  }
  if(per_method == 1){
    for(i=0; i<DIM; i++){
      n_false[i]=0;
      n_true[i]=0;
      for(j=0; j<n_peaks[i]-2; j++){
	E[j][i] = (double)per[0][i] - (double)per[j+1][i];
	if(E[j][i] < 0)
	  E[j][i] = -E[j][i];
      }
      // Set to zeros ALL spectrum_per's elements
      for(j=0;j<MAX_N_ISEC-2;j++){
	*(spectrum_per+j*MAXDIM+i)=FALSE;
      }
      // Set to zeros ALL n_subperiods's elements
      *(n_subperiods+i)=0;
    }
    for(i=0; i<DIM; i++){
      period = per[0][i];
      for(j=0; j<n_peaks[i]-2; j++){
	if(E[j][i] > D(eps_abs_per,eps_rel_per,period)){
	  *(spectrum_per+j*MAXDIM+i) = FALSE;
	  n_false[i]+=1;}
	if(E[j][i] < D(eps_abs_per,eps_rel_per,period)){
	  *(spectrum_per+j*MAXDIM+i) = TRUE;
	  if(n_true[i] == 0)
	    *(n_subperiods+i) = j+1;
	  n_true[i]+=1;}
      }
    }
  }
  for(i=0; i<DIM; i++){
    if(n_false[i] == 0)
      *(n_subperiods+i) = 0;
  }
  return 0;
}

int get_index_trans(const double **frame, const int frame_size, const double trans_time)
{/* This function is intended for finding index in the loaded <frame> close to the
    time threshold <trans_time> defined by the user. This function is useful, for
    example, when one needs to skip some initial amount of data corresponding to the
    transient dynamics until the system reaches the attractor. The found index is
    pointing to the time value >= trans_time. The function returns the index. */
  /* ************************************************************************* */
  /* Time is always stored in the first row of frame. */
  int i;
  for(i=0;i<frame_size;i++)
    if(frame[0][i] >= trans_time)
      return i;
  return -1;/* Not found */
}

int mol_dist(const int bin)
{/* The function creates the molecular/species distribution of the stochastic time
    series. It takes the current data file as a source for downloading the
    frame. Only the first stochastic frame is loaded. Number of bins for the
    histogram is specified by bin.*/
  int data_offset,i,frame_size;
  double max_x[DIM],min_x[DIM],ampl[DIM];
  /* First get info about the time series */
  int *info;
  FILE *in = fopen(data_name,"r");
  if((info=get_info_data(info,in)) == NULL){
    fprintf(stderr,"Error occurred. Exit.");
    return 10;
  }
  fclose(in);
  /* Then we load the frames */
  double **frame;
  in = fopen(data_name,"r");
  if(info[0])
    data_offset = 1;/* We need to skip first frame, since it is determ */
  else/* Only if not complex run then it is possible to have hist for determ */
    data_offset = 0;
  for(i=0;i<data_offset;i++){/* Skip the offset frames, at max 1 */
    frame = load_frame(in,frame,&frame_size);
    if(frame == NULL){
      fprintf(stderr,"Error: frame was not loaded.\n");
      return 14;
    }
    free_frame(frame);
  }
  /* Load only one (first) stoch frame */
  printf("Warning: only 1st frame is analyzed.\n");
  frame = load_frame(in,frame,&frame_size);
  if(frame == NULL){
    fprintf(stderr,"Error: frame was not loaded.\n");
    return 15;
  }
  /* Finally, we compute histogram */
  double **hist;
  if(per_var == -1){
    /* Automatic detection of the variable: var with the largest absolute amplitude*/
    for(i=0;i<DIM;i++){/* Note that first row of frame is Time */
      max_x[i] = gsl_stats_max(*(frame+i+1),1,frame_size);
      min_x[i] = gsl_stats_min(*(frame+i+1),1,frame_size);
      ampl[i] =  max_x[i] - min_x[i];
    }
    /* Take variable with the largest amplitude */
    perVarInd = gsl_stats_max_index(ampl,1,DIM);
    hist = compute_hist_per(hist,frame[perVarInd+1],	\
			    frame_size,bin);
  }
  else{
    perVarInd = per_var;
    hist = compute_hist_per(hist,frame[perVarInd+1],	\
			    frame_size,bin);
  }
  write_histPer_file(hist,bin,"hist.dat","w");
  if(graph_flag){/* print histogram */
    graph_set_labels(graph,"mdist");
    gnuplot_cmd(plot_handle,"set title \"n=%d, for %s var\"\n",
		frame_size,var_name[perVarInd]);
    gnuplot_cmd(plot_handle,"plot 'hist.dat' u 1:2 w boxes ls 1\n");
  }
  free_frame(frame);
  for(i=0;i<2;i++)
    free(hist[i]);
  free(hist);
  return 0;
}

int ss_stab(const double *xss)
{
  /* This function computes the stability of the steady state (SS). It
     returns 0 in case the SS is stable and 1 otherwise. It
     takes as the input values vector xss of the SS. */
  int i,j;
  gsl_complex tmp;
  double dfdx[DIM*DIM];
  double dfdt[DIM];
  double dxdot[DIM];/* Derivatives of the increments, i.e. deviations from SS */
  jac(t,xss,dfdx,dfdt,&mu);
  /* Allocating workspace */
  gsl_eigen_nonsymm_workspace *w=gsl_eigen_nonsymm_alloc(DIM);
  /* Creating the matrix out of Jacobian */
  gsl_matrix_view jacob=gsl_matrix_view_array(dfdx,DIM,DIM);
  /* printf("Reconstructing Jacobian:\n"); */
  /* printf("J=\n"); */
  /* for(i=0;i<DIM;i++){ */
  /*   for(j=0;j<DIM;j++) */
  /*     printf("%G ",gsl_matrix_get(&jacob.matrix,i,j)); */
  /*   printf("\n"); */
  /* } */
  /* eigenvalues are in eval */
  gsl_vector_complex *eval=gsl_vector_complex_alloc(DIM);
  gsl_eigen_nonsymm(&jacob.matrix,eval,w);
  j=0;
  for(i=0;i<DIM;i++){
    tmp=gsl_vector_complex_get(eval,i);
    printf("L_%d=%G + %Gi\n",i+1,GSL_REAL(tmp),GSL_IMAG(tmp));
    if(GSL_REAL(tmp)>1e-7) j++;
  }
  /* Freeing the memory */
  gsl_vector_complex_free(eval);
  gsl_eigen_nonsymm_free(w);
  if(j>0) return 1;
  else return 0;
}
