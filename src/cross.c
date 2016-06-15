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

#include "init.h"
#include "errors.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>
#ifndef TISEC_CHUNK_SIZE
#define TISEC_CHUNK_SIZE 500
#endif /* #ifndef TISEC_CHUNK_SIZE */
#ifndef TRUE
#define TRUE 1
#endif /* #ifndef TRUE */
#ifndef FALSE
#define FALSE 0
#endif /* #ifndef FALSE */


void write_autocorr_file(double autocorr[],const int frame_size,
			 char *out_fname,const char *mode)
{
  int i;
  FILE *out = fopen(out_fname,mode);
  if(out == NULL)
    fprintf(stderr,"Could not open file: %s\n",out_fname);
  else{
    for(i=0;i<frame_size;i++)
      fprintf(out,"%G\n",autocorr[i]);
    fprintf(out,"\n\n");
  }
  fclose(out);
}

double *compute_autocorrelation(double **frame, const int frame_size, \
				double *autocorr, const int var_ind)
{/* This function computes the auto-correlation function for a single frame, i.e. for
    a single stochastic trajectory and for a specified variable(defined with var_ind,
    which must be in C style, i.e. starting with 0; note also that frame contains
    time in the first column). It returns the pointer to 2d array of autocorrelation
    function values for all DIM variables. The memory allocation for autocorr is done
    within the function. */
  int i;
  int eff_frame_size;/*Effective frame size, the one that the frame was allocated for*/
  if(frame_size <= FRAME_N_LINES)/* No realloc-ion took place when frame was read */
    eff_frame_size = FRAME_N_LINES;
  else
    eff_frame_size = frame_size;/* For the reallocated frame's memory */
  autocorr = (double *)malloc(frame_size*sizeof(double));
  for(i=0;i<frame_size-1;i++){/* Ignore the last value to avoid NaN */
    autocorr[i] = 
      gsl_stats_covariance(*(frame+var_ind+1)+i,1,	\
      			   *(frame+var_ind+1),1,frame_size-i);
  }
  for(i=1;i<frame_size-1;i++)/* Normalize all values except the first */
    autocorr[i] = autocorr[i] /	autocorr[0];
  //      gsl_stats_covariance(*(frame+var_ind+1),1,*(frame+var_ind+1),1,frame_size);
  autocorr[0] = autocorr[0] / autocorr[0];/* Formally... */
  return autocorr;
}

void thist_init()
{
  hist = (char *)malloc(((strlen("hist.dat"))+1)*sizeof(char));
  strncpy(hist,"hist.dat",strlen("hist.dat"));
}

void thist_free()
{
  free(hist);
}

void write_histPer_file(double **histPer,const int bin,
			 char *out_fname,const char *mode)
{
  int i;
  FILE *out = fopen(out_fname,mode);
  if(out == NULL)
    fprintf(stderr,"Could not open file: %s\n",out_fname);
  else{
    for(i=0;i<bin+1;i++)/* we always have 1 more than 'bin' bins */
      fprintf(out,"%G %G\n",histPer[1][i],histPer[0][i]);
    fprintf(out,"\n\n");
  }
  fclose(out);
}

double **compute_hist_per(double **histPer,double *period,
			  const int numPer,const int bin)
{/* The function computes the normalized histogram of periods. In the DINAMICA
    conventions that histPer, returned by the function, is 2 rows double array: 1st
    row gives percent of values falls into a bin, center for which is returned in the
    2nd row. period is a double array containing the period values, while bin is the
    number of bins for the histogram to compute. numPer is the number of values in
    period. The actual number of bins = bin+1, since the upper bound of a bin is
    exclusive (while lower bound is inclusive) and max element of period is not
    included, thus additionaly bin is created for this single value.*/
  int i,j,n;
  double perMax,perMin,dT;
  histPer = (double **)malloc(2*sizeof(double *));
  for(i=0;i<2;i++)
    histPer[i] = (double *)malloc((bin+1)*sizeof(double));
  perMax = gsl_stats_max(period,1,numPer);
  perMin = gsl_stats_min(period,1,numPer);
  dT = (perMax - perMin) / bin;
  for(i=0;i<bin+1;i++){/* Additional bin to include perMax value */
    n = 0;
    for(j=0;j<numPer;j++){
      if((period[j] >= (i*dT+perMin)) && (period[j] < ((i+1)*dT+perMin)))
	n++;
    }
    histPer[0][i] = (double)n / (double)numPer;
    histPer[1][i] = perMin + i*dT + dT/2;/* centers of the bins */
  }
  return histPer;
}

size_t get_var_per(char *data_fname, const int data_offset)
{
  /* This function determines which variable to use for the period calculation. For
     this, it gets loads single stoch frame and determines the most abundant species
     in it, by simply getting the difference between max and min of all variables and
     comparing all these values between each other. It returns the index of the
     desired variable (in C style, starting from 0). data_offset is used to ignore
     that number of frames from the beginning of the source file. Note that the total
     number of frames in the source must be at least (data_offset+1), otherwise,
     expect undefined outcome of load_frame(). */
  int i,frame_size;
  double **frame;
  double diff[DIM];/* For max and min difference values */
  FILE *load = fopen(data_fname,"r");
  if(load == NULL){
    fprintf(stderr,"Could not open the file %s\n",data_fname);
    return -1;/* Not valid value */
  }
  else{
    for(i=0;i<data_offset;i++){/* Ignore first data_offset frames of the data */
      frame = load_frame(load,frame,&frame_size);
      free_frame(frame);
    }
    /* Load the next frame */
    frame = load_frame(load,frame,&frame_size);
    if(frame == NULL)
      return -1;
    for(i=0;i<DIM;i++)/* Note that first row of frame is Time */
      diff[i] = gsl_stats_max(*(frame+i+1),1,frame_size) -
	gsl_stats_min(*(frame+i+1),1,frame_size);
    free_frame(frame);
    fclose(load);
    return (gsl_stats_max_index(diff,1,DIM));
  }
}

void write_perStoch_file(double perStoch[],const int numPer,
			 char *out_fname,const char *mode)
{
  int i;
  FILE *out = fopen(out_fname,mode);
  if(out == NULL)
    fprintf(stderr,"Could not open file: %s\n",out_fname);
  else{
    for(i=0;i<numPer;i++)
      fprintf(out,"%G\n",perStoch[i]);
    fprintf(out,"\n\n");
  }
  fclose(out);
}

double *compute_period(double *period, int *nper, double cross_level,
		       const int meth, const char *data_fname,
		       const int n_frame, const int data_offset)
{
  /* This function computes periods for a given variable from a set of
     frames(timeseries) in one data file. n_frame is the number of frames to
     process. data_offset is the number of frames to ignore starting from the
     beginning of the data file, for example, in complex run first frame is the
     deterministic simulation. `period' might not be initialized(assigned with memory
     block) elsewhere. data_fname is the name of the data file to load. meth is the
     type of the algorithm for the period calculation. meth can be: "0" for Poincare
     sectioning and "1" for autocorrelation methods, respectively. *nper is the total
     number of periods computed. NOTE: the function does not distinguish between the
     deterministic and stochastic parts in the source file, given the simulation
     method is complex. Returned value is NULL in the case of memory allocation
     problems (for frames, tisec or periods) and in the case of 0 number of periods,
     i.e. *nper=0. Meanwhile to avoid additional computation max, min and amplitude
     values are computed by the function, given that these entities are global to the
     program. The globals also are cross, perVarInd, per_var, per_thresh...*/
  int i,j,k,frame_size,nisec;
  double **frame,*acorr,*tisec;
  double ampl[DIM],max_x[DIM],min_x[DIM];
  //FILE *outfn = fopen("out","w");
  FILE *load = fopen(data_fname,"r");
  if(period != NULL)
    period = 0;
  if(meth == 1){/* Known number of periods for autocorrelation method */
    period = (double *)realloc(period,n_frame*sizeof(double));
    if(period == NULL){
      fprintf(stderr,"Error:period was NOT initialized\n");
      return NULL;
    }
    (*nper) = n_frame;
  }
  else{
    (*nper) = 0;
  }
  if(load == NULL)
    fprintf(stderr,"Error: could not open file \"%s\"\n",data_fname);
  else{
    for(i=0;i<data_offset;i++){/* Ignore first data_offset frames of the data */
      frame = load_frame(load,frame,&frame_size);
      if(frame == NULL){
	fprintf(stderr,"Error: frame was not loaded\n");
	return NULL;
      }
      free_frame(frame);
    }
    for(i=0;i<n_frame;i++){
      frame = load_frame(load,frame,&frame_size);
      if(frame == NULL){
	fprintf(stderr,"Error: frame was not loaded\n");
	return NULL;
      }
      //write_frame(outfn,frame,frame_size);
      /* Determine the variable for the period, while computing max, min and
	 amplitude. */
      if(per_var == -1){
	/* Only when required for determining the per var we compute the */
	/* amplitudes and max/min values */
	for(j=0;j<DIM;j++){/* Note that first row of frame is Time */
	  max_x[j] = gsl_stats_max(*(frame+j+1),1,frame_size);
	  min_x[j] = gsl_stats_min(*(frame+j+1),1,frame_size);
	  ampl[j] =  max_x[j] - min_x[j];
	}
	/* Take variable with the largest amplitude */
	perVarInd = gsl_stats_max_index(ampl,1,DIM);
      }
      else
	perVarInd = per_var;
      /*** Poincare stuff ***/
      if(meth == 0){
	/* Define max and min. Compute section as (max-min)/cross_level + min. Find
	   the intersections of the section with the trajectory */
	tisec = crossing(&nisec,&cross,cross_level,frame,frame_size,perVarInd);
	if(tisec == NULL)
	  return NULL;
	/* Compute the periods */
	if(nisec > 1){/* Only when there were at least 2 intersections */
	  period = period_cross(period,nper,nisec,tisec,per_thresh);
	  if(period == NULL)
	    return NULL;
	}
	/* Free memory for the intersections */
	free(tisec);
      }
      /*** Autocorrelation stuff ***/
      if(meth == 1){
	/* Compute autocorrelation */
	acorr = compute_autocorrelation(frame,frame_size,acorr,perVarInd);
	/* Write the autocorrelation to a file */
	if((fopen("acorr.dat","r")) == NULL)
	  write_autocorr_file(acorr,frame_size,"acorr.dat","w");
	else
	  write_autocorr_file(acorr,frame_size,"acorr.dat","a");
	/* Find the first peak of the autocorrelation giving the period */
	if((k = compute_peak_autocorr(acorr,frame_size)) >= 0)
	  period[i] = frame[0][k];
      }
      free_frame(frame);
    }
    fclose(load);
    //fclose(outfn);
  }
  return period;
}

double * period_cross(double *period, int *nper, const int nisec,
		      const double *tisec, const double per_thresh)
{/* The function computes period array from the array of intersections of the
    Poincare section with the trajectory for the given variable. period is set to
    NULL in this function, it can be redundant though. nper is the number of periods
    computed so far, i.e. the period array is already of nper size when comes here
    and additionally nper is subject to change inside the function. tisec in the
    current array of intersection time moments and nisec is number of intersections
    in the current tisec. The function will return NULL in period in two cases:
    memory allocation problem or there were not enough intersections (less than 2) to
    compute a single period value, thus, it is preferrable to check the number of
    intersections nisec before this call. per_thresh is used to merge the periods
    less than the value with the subsequent period untill the sum exceeds the
    per_thresh.*/
  int i;
  double tmp;
  if((period != NULL) && ((*nper) == 0))
    period = 0;
  i = 1;
  while(i < nisec){/* Iterate over the intersections */
    tmp = tisec[i] - tisec[i-1];
    /* If less than threshold, iterate untill we get enough by summing */
    while((tmp < per_thresh) && (i < nisec)){
      i++;
      tmp += (tisec[i] - tisec[i-1]);
    }
    if(i >= nisec)/* We met the end of the intersection array */
      break;
    /* We got the desired period, let us allocate the memory for it */
    period = (double *)realloc(period,((*nper)+1)*sizeof(double));
    if(period == NULL){/* Check for bad allocation */
      fprintf(stderr,"Error(period_cross): could not realloc the period\n");
      return NULL;
    }
    period[(*nper)] = tmp;
    (*nper)++;
    i++;
  }
  
  return period;
}

double * crossing(int *nisec, double *cross, const double cross_level,
		  double **frame, const int frame_size, const int varInd)
{/* The function computes the Poincare section (cross) for the variable (varInd) in
    the data frame (frame) of size (frame_size), and the time moments (tisec which is
    returned) of intersections between the trajectory and the section. These
    intersections can be later used to compute periods for the data. The section is
    made on the height being the fraction of the amplitude(max-min) of the
    signal. The fraction is determined by the cross_level. So to have section on the
    half of the amplitude cross_level should be equal to 2, i.e. CROSS =
    (MAX-MIN)/CROSS_LEVEL + MIN. frame is not changed within this call.
 */
  int i;
  double *tisec;
  int alloc_count;
  /*Automatic Poincare section*/
  *cross = (gsl_stats_max(frame[varInd+1],1,frame_size) -
	    gsl_stats_min(frame[varInd+1],1,frame_size)) / cross_level +
    gsl_stats_min(frame[varInd+1],1,frame_size);

  (*nisec) = 0;
  /* Allocate memory for the 1st time for tisec. To avoid multiple realloc
     invocations, we will allocate a big size chunk of memory. Every time the tisec
     container is full we will resize the container with realloc. This would reduce the
     number of realloc invocations. This requres an addition counter on how many times
     we did the allocation of memory. */
  alloc_count = 1;
  tisec = (double *)malloc(alloc_count*TISEC_CHUNK_SIZE*sizeof(double));
  if(tisec == NULL){
    fprintf(stderr,"Error: tisec was NOT initialized.\n");
    return NULL;
  }
  /* Start iterate over the timeseries */
  for(i=1;i<frame_size;i++){
    if((frame[varInd+1][i] >= (*cross)) && (frame[varInd+1][i-1] < (*cross))){
      /* Found the intersection, calculate it */
      tisec[(*nisec)] = ( (*cross)*(frame[0][i-1] - frame[0][i]) -	\
			  frame[0][i-1]*frame[varInd+1][i] +		\
			  frame[0][i]*frame[varInd+1][i-1] ) /		\
	(frame[varInd+1][i-1] - frame[varInd+1][i]);
      /* Increase the number of found intersections */
      (*nisec)++;
      /*Reallocate memory for the tisec, when nisec is larger than TISEC_CHUNK_SIZE*/
      if((*nisec) == (alloc_count*TISEC_CHUNK_SIZE)){
	alloc_count++;/* increase the counter */
	/* Realloc */
	tisec = (double *)realloc(tisec,alloc_count*TISEC_CHUNK_SIZE*sizeof(double));
	/* Check for errors */
	if(tisec == NULL){
	  fprintf(stderr,"Error: tisec was NOT re-initialized: nisec = %d.\n",*nisec);
	  return NULL;
	}
      }
    }
  }
  /* We do not need to resize the final tisec to get away the unused blocks of
     memory, since we also use nisec variable everywhere where tisec is needed. */
  return tisec;
}

int compute_peak_autocorr(const double *autocorr, const int size)
{
  /* This function computes the first peak of the autocorrelation function. The
     autocorrelation function is supplied with double array 'autocorr' of size
     'size'. The peak is returned as an index the array 'autocorr', in order to get
     the period, we have to put this index to the first sector/row (Time moments) of
     the frame, this autocorr was computed from. */
  int i,pstv_intrsc = 0;
  double peak;
  /* We do an assumption that our autocorrelation function crosses the zero to
     negative values at least once. So we determine the first peak as a maximum value
     between intersections of first positive slope of the function and the following
     negative slope of the function. */
  for(i=1;i<size;i++){
    /* We can skip first value of autocorr, being always 1.We need only the first
       peak of autocorr. */
    if(((autocorr[i]-autocorr[i-1]) > 0) &&		\
       (autocorr[i] > 0) && (autocorr[i-1] < 0)){
      /* Positive slope, the function crosses the 0 value. */
      pstv_intrsc = i;/* The first intersection is found */
      continue;/* continue further */
    }
    if(((autocorr[i]-autocorr[i-1]) < 0) &&		\
       (autocorr[i-1] > 0) && (autocorr[i] < 0) &&	\
       (pstv_intrsc > 0)){
      /* Negative slope, the function crosses the 0 value and the preceding positive
	 slope intersection is found */
      peak = pstv_intrsc +
	gsl_stats_max_index(autocorr+pstv_intrsc,1,i-pstv_intrsc);
      break;
    }
  }
  if(i == size)/* We got the end of the loop,i.e. no peak was found */
    return -1;/* Nonscence negative value to return */
  return peak;
}

/* int intersec(int n_isec[],double *tper, int(*func)(double,const double [],double[],void *)) */
/* { */
/*   int i,j,m; */
/*   double static x_prev[MAX_LOCAL_DIM*MAX_TOTAL_DIM];  */
/*   double static t_prev[MAX_LOCAL_DIM*MAX_TOTAL_DIM]; */

/*   func(t,x,f,&mu); */
/*   j=1; */
/*   while((n_steps == 0) && (j <= DIM)){ */
/*     x_prev[j-1] = x[j-1]; */
/*     t_prev[j-1] = t; */
/*     n_isec[j-1]=0; */
/*     for(i=0;i<MAX_N_ISEC;i++){ /\* Set to zero DIM*MAX_N_ISEC elements in tper array *\/ */
/*       *(tper+i*MAXDIM+j-1)=0.0;        /\* (not MAXDIM*MAX_N_ISEC) *\/ */
/*     } */
/*     j++; */
/*   } */
/*   for(i=1; i <= DIM; i++ ) { */
/*     if(*(f+i-1) >= 0){ */
/*       if((x[i-1] < x_cross[i-1]) || (n_steps == 0)) { */
/* 	x_prev[i-1] = x[i-1]; */
/* 	t_prev[i-1] = t; */
/*       } */
/*       if( x[i-1] == x_cross[i-1] ) { */
/* 	*(tper+ n_isec[i-1]*MAXDIM+i-1) = (t - t_prev[i-1])*(x_cross[i-1] - x[i-1])/(x[i-1] - x_prev[i-1]) + t; */
/* 	x_prev[i-1] = x[i-1]; */
/* 	t_prev[i-1] = t; */
/* 	fprintf(stdout,"BINGO\n"); */
/* 	n_isec[i-1]++; */
/*       } */
/*       if( (x[i-1] > x_cross[i-1]) && (x_prev[i-1] < x_cross[i-1]) ) { */
/* 	*(tper+ n_isec[i-1]*MAXDIM+i-1) = (t - t_prev[i-1])*(x_cross[i-1] - x[i-1])/(x[i-1] - x_prev[i-1]) + t; */
/* 	x_prev[i-1] = x[i-1]; */
/* 	t_prev[i-1] = t; */
/* 	n_isec[i-1]++; */
/* 	if(n_isec[i-1] >= MAX_N_ISEC){ */
/* 	  fprintf(stdout,"Too many intersections\n"); */
/* 	  return -1; */
/* 	} */
/* 	else if(((max_x[i-1]-min_x[i-1])<D(eps_abs_am,eps_rel_am,(float)0.0)) \ */
/* 		&& (n_steps != 0)){ */
/* 	  return -2; */
/* 	} */
/*       } */
/*     } */
/*   } */

/*   return 0; */
/* } */

int period(float *per)
{
  int i,k,j;
    for(k=1; k<=DIM; k++){
      if(per_method==0) j=n_isec[k-1]-1;
      if(per_method==1) j=n_peaks[k-1]-1;
      for(i=1; i<=j; i++){
	*(per +(i-1)*MAXDIM+k-1) = 0.0;
	if(per_method==0)
	  *(per+(i-1)*MAXDIM+k-1) = (float)tper[i][k-1] - (float)tper[i-1][k-1];
	if(per_method==1)
	  *(per+(i-1)*MAXDIM+k-1) = (float)t_peak[i][k-1] - (float)t_peak[i-1][k-1];
      }
    }
  return 0;
}

int period_average(float *per_aver,float *big_per)
{
  int i,j,k,l;
  int cnt;
  for(i=0; i<DIM; i++){
    *(big_per+i) = 0;
    k=n_subperiods[i];
    if(n_subperiods[i] == 0) k=k+1;
    for(l=0;l<k;l++){
      cnt=0;
      *(per_aver+l*MAXDIM+i) = 0;
      if(per_method==0){
	for(j=l; j<n_isec[i]-1; j+=k){
	  *(per_aver+l*MAXDIM+i) += per[j][i];
	  cnt+=1;
	}
      }
      if(per_method==1){
	for(j=l; j<n_peaks[i]-1; j+=k){
	  *(per_aver+l*MAXDIM+i) += per[j][i];
	  cnt+=1;
	}
      }
      *(per_aver+l*MAXDIM+i) /= cnt;
      *(big_per+i) += *(per_aver+l*MAXDIM+i);
    }
  }
    
  return 0;
}

double * ampl_multi_frames(double *ampl,int *ampl_cum_size)
{
  /* Function collects amplitudes from several frames and return a single array of
     them. */
  trajPeak *peak;
  double **frame,*ampl_tmp;
  int frame_size,*info,ampl_size,i,j;
  *ampl_cum_size = 0;
  ampl = (double *)malloc(sizeof(double));
  char *fname = (char *)malloc((strlen(data_name)+4)*sizeof(char));
  fname = strcpy(fname,data_name);
  if(ma_span){
    printf("I analyze .ma file since MA span is specified.\n");
    fname = strcat(fname,".ma");
  }
  /* First read the header info */
  FILE *in = fopen(fname,"r");
  if(in == NULL){
    printf("Error openning file %s\n",fname);
    return NULL;
  }
  info = get_info_data(info,in);
  fclose(in);
  /* Second read the frame(s) */
  in = fopen(fname,"r");
  for(i=0;i<info[1];i++){
    frame = load_frame(in,frame,&frame_size);
    peak = peak_trough2(frame,frame_size,graph.yInd[0],peak);
    ampl_tmp = peak_ampl(*peak,ampl_tmp,&ampl_size);
    *ampl_cum_size += ampl_size;
    ampl = (double *)realloc(ampl,*ampl_cum_size*sizeof(double));
    //printf("Amplitudes from frame %d:\n",i+1);
    for(j=0;j<ampl_size;j++){
      ampl[*ampl_cum_size - ampl_size + j] = ampl_tmp[j];
      //printf("%G\n",ampl_tmp[j]);
    }
    free_frame(frame);
    free(peak);
    free(ampl_tmp);
  }
  free(info);
  free(fname);

  return ampl;
}

double * peak_ampl(trajPeak peak,double *ampl,int *ampl_size)
{
  /* Calculate the amplitudes from the trajPeak object (see peak_trough function(s)).
   */
  int i,j,start_index;
  double min,max;
  /* First allocate more than we need: peak.size is more than number of amplitudes. */
  ampl = (double *)malloc(peak.size*sizeof(double));
  j = 0;
  /* Determine the starting index. */
  if(peak.peak_true_table[1])
    start_index = 1;
  if(peak.peak_true_table[2])
    start_index = 2;
  for(i=start_index;i<peak.size-1;i+=2){
    /* Determine the peak (MAX) and take the min from two surrounding troughs
       (MIN). Amplitude is MAX-MIN. Iterate by 2, i.e. every other value is true
       peak.
    */
    /* peak.x[i] MUST be a true peak */
    max = peak.x[i];
    min = peak.x[i-1];/* the preceding trough */
    if(min > peak.x[i+1])
      min = peak.x[i+1];/* the next trough */
    /* Ampl = max - min */
    ampl[j] = max - min;
    //printf("max = %G, min = %G, ampl = %G\n",max,min,ampl[j]);
    j++;
  }
  ampl = (double *)realloc(ampl,j*sizeof(double));
  *ampl_size = j;
  return ampl;
}

void free_peak(trajPeak *peak)
{
  free(peak->t);
  free(peak->x);
  free(peak->peak_true_table);
  free(peak);
}

trajPeak *peak_trough2(double **frame,const int frame_size,const int varInd,
		       trajPeak *peak)
{
  int i,j;
  peak = (trajPeak *)malloc(sizeof(trajPeak));
  if(peak == NULL){
    fprintf(stderr,"Error(peak_trough): could not init peak\n");
    return NULL;
  }
  peak->size = 0;
  peak->t = (double *)malloc((peak->size)*sizeof(double));
  peak->x = (double *)malloc((peak->size)*sizeof(double));
  peak->peak_true_table = (int *)malloc((peak->size)*sizeof(int));
  j = 1;
  while(j<(frame_size-1)){/* Skip the first and the last values */
    i = 0;
    if((frame[varInd+1][j] - frame[varInd+1][j-1]) > 0){
      /* ASCENDING PART */
      if(frame[varInd+1][j] == frame[varInd+1][j+1]){
	/* We have plateau */
	i = j+1;/* Fix the point where the plateau starts */
	while(frame[varInd+1][j] == frame[varInd+1][j+1])
	  j++;/* Skip the interval of plateau */
      }
      if((frame[varInd+1][j] - frame[varInd+1][j+1]) < 0){
	/* Continue of ascending */
	j++;
	continue;
      }
      else if((frame[varInd+1][j] - frame[varInd+1][j+1]) > 0){
	/* We got peak: plateau between <i> and <j>(including <i>, excluding <j>)*/
	/* Find average between <i> and <j>. In case of non-integral we take ceiling
	   of that. <i> will hold the index of the peak. */
	if(i){/* If we were at plateau before */
	  if(fmod(i+j,2) > 0.1)/* Make sure of non-integral value */
	    i = (int)ceil((i+j)/2);
	  else
	    i = (i+j) / 2;
	}
	else/* Oppositely we have peak at <j>. */
	  i = j;
	(peak->size)++;
	peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
	peak->t[peak->size-1] = frame[0][i];/* Time moment of peak */
	peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
	peak->x[peak->size-1] = frame[varInd+1][i];/* Value of var in peak */
	peak->peak_true_table =
	  (int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
	peak->peak_true_table[peak->size-1] = 1;/* It is the peak */
      }
    }
    if((frame[varInd+1][j] - frame[varInd+1][j-1]) < 0){
      /* DESCENDING PART */
      if(frame[varInd+1][j] == frame[varInd+1][j+1]){
	/* We have plateau */
	i = j+1;/* Fix the point where the plateau starts */
	while(frame[varInd+1][j] == frame[varInd+1][j+1])
	  j++;/* Skip the interval of plateau */
      }
      if((frame[varInd+1][j] - frame[varInd+1][j+1]) > 0){
	/* Continue descending */
	j++;
	continue;
      }
      else if((frame[varInd+1][j] - frame[varInd+1][j+1]) < 0){
	/* We got trough: plateau between <i> and <j>(including <i>, excluding <j>)*/
	/* Find average between <i> and <j>. In case of non-integral we take ceiling
	   of that. <i> will hold the index of the peak. */
	if(i){
	  if(fmod(i+j,2) > 0.1)/* Make sure of non-integral value */
	    i = (int)ceil((i+j)/2);
	  else
	    i = (i+j) / 2;
	}
	else
	  i = j;
	(peak->size)++;
	peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
	peak->t[peak->size-1] = frame[0][i];/* Time moment of trough */
	peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
	peak->x[peak->size-1] = frame[varInd+1][i];/* Value of var in trough */
	peak->peak_true_table =
	  (int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
	peak->peak_true_table[peak->size-1] = 0;/* It is not a peak==trough */
      }
    }
    j++;
  }
  return peak;
}

trajPeak *peak_trough1(double **frame,const int frame_size,const int varInd,
		       trajPeak *peak)
{
  int j;
  peak = (trajPeak *)malloc(sizeof(trajPeak));
  if(peak == NULL){
    fprintf(stderr,"Error(peak_trough): could not init peak\n");
    return NULL;
  }
  peak->size = 0;
  peak->t = (double *)malloc((peak->size)*sizeof(double));
  peak->x = (double *)malloc((peak->size)*sizeof(double));;
  peak->peak_true_table = (int *)malloc((peak->size)*sizeof(int));;
  j = 1;
  while(j<(frame_size-1)){/* Skip the first and the last values */
    if((frame[varInd+1][j] > frame[varInd+1][j-1]) &&
       (frame[varInd+1][j] > frame[varInd+1][j+1])){
      (peak->size)++;
      peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
      peak->t[peak->size-1] = frame[0][j];/* Time moment of peak */
      peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
      peak->x[peak->size-1] = frame[varInd+1][j];/* Value of var in peak */
      peak->peak_true_table =
	(int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
      peak->peak_true_table[peak->size-1] = 1;/* It is the peak */
    }
    else if((frame[varInd+1][j] < frame[varInd+1][j-1]) &&
	    (frame[varInd+1][j] < frame[varInd+1][j+1])){
      /* DESCENDING PART */
      (peak->size)++;
      peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
      peak->t[peak->size-1] = frame[0][j];/* Time moment of trough */
      peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
      peak->x[peak->size-1] = frame[varInd+1][j];/* Value of var in trough */
      peak->peak_true_table =
	(int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
      peak->peak_true_table[peak->size-1] = 0;/* It is not a peak==trough */
    }
    else if(frame[varInd+1][j] == frame[varInd+1][j-1]){
      printf("Equal:frame[%d][%d] = %lf!\n",varInd+1,j,frame[varInd+1][j]);
    }
    j++;
  }
  return peak;
}

trajPeak *peak_trough(double **frame,const int frame_size,const int varInd,
		      trajPeak* peak)
{
  int i,j;
  peak = (trajPeak *)malloc(sizeof(trajPeak));
  if(peak == NULL){
    fprintf(stderr,"Error(peak_trough): could not init peak\n");
    return NULL;
  }
  peak->size = 0;
  peak->t = (double *)malloc((peak->size)*sizeof(double));
  peak->x = (double *)malloc((peak->size)*sizeof(double));;
  peak->peak_true_table = (int *)malloc((peak->size)*sizeof(int));;
  j = 1;
  while(j<(frame_size-1)){/* Skip the first and the last values */
    i = 0;
    if((frame[varInd+1][j] - frame[varInd+1][j-1]) > eps_abs_int){
      /* ASCENDING PART */
      if(fabs(frame[varInd+1][j] - frame[varInd+1][j+1]) < eps_abs_int){
	/* We have plateau */
	i = j+1;/* Fix the point where the plateau starts */
	while(fabs(frame[varInd+1][j] - frame[varInd+1][j+1]) < eps_abs_int)
	  j++;/* Skip the interval of plateau */
      }
      if((frame[varInd+1][j] - frame[varInd+1][j+1]) < -eps_abs_int){
	/* Continue of ascending */
	j++;
	continue;
      }
      else if((frame[varInd+1][j] - frame[varInd+1][j+1]) > eps_abs_int){
	/* We got peak: plateau between <i> and <j>(including <i>, excluding <j>)*/
	/* Find average between <i> and <j>. In case of non-integral we take ceiling
	   of that. <i> will hold the index of the peak. */
	if(i){/* If we were at plateau before */
	  if(fmod(i+j,2) > 0.1)/* Make sure of non-integral value */
	    i = (int)ceil((i+j)/2);
	  else
	    i = (i+j) / 2;
	}
	else/* Oppositely we have peak at <j>. */
	  i = j;
	(peak->size)++;
	peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
	peak->t[peak->size-1] = frame[0][i];/* Time moment of peak */
	peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
	peak->x[peak->size-1] = frame[varInd+1][i];/* Value of var in peak */
	peak->peak_true_table =
	  (int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
	peak->peak_true_table[peak->size-1] = 1;/* It is the peak */
      }
    }
    if((frame[varInd+1][j] - frame[varInd+1][j-1]) < -eps_abs_int){
      /* DESCENDING PART */
      if(fabs(frame[varInd+1][j] - frame[varInd+1][j+1]) < eps_abs_int){
	/* We have plateau */
	i = j+1;/* Fix the point where the plateau starts */
	while(fabs(frame[varInd+1][j] - frame[varInd+1][j+1]) < eps_abs_int)
	  j++;/* Skip the interval of plateau */
      }
      if((frame[varInd+1][j] - frame[varInd+1][j+1]) > eps_abs_int){
	/* Continue descending */
	j++;
	continue;
      }
      else if((frame[varInd+1][j] - frame[varInd+1][j+1]) < -eps_abs_int){
	/* We got trough: plateau between <i> and <j>(including <i>, excluding <j>)*/
	/* Find average between <i> and <j>. In case of non-integral we take ceiling
	   of that. <i> will hold the index of the peak. */
	if(i){
	  if(fmod(i+j,2) > 0.1)/* Make sure of non-integral value */
	    i = (int)ceil((i+j)/2);
	  else
	    i = (i+j) / 2;
	}
	else
	  i = j;
	(peak->size)++;
	peak->t = (double *)realloc(peak->t,(peak->size)*sizeof(double));
	peak->t[peak->size-1] = frame[0][i];/* Time moment of trough */
	peak->x = (double *)realloc(peak->x,(peak->size)*sizeof(double));
	peak->x[peak->size-1] = frame[varInd+1][i];/* Value of var in trough */
	peak->peak_true_table =
	  (int *)realloc(peak->peak_true_table,(peak->size)*sizeof(int));
	peak->peak_true_table[peak->size-1] = 0;/* It is not a peak==trough */
      }
    }
    j++;
  }
  return peak;
}

int peak()
{
  int i,j,k;
  double f_prev_peak[DIM];
  double f_prev_cavity[DIM];
  double local_ampl[DIM];
  
  for(i=0;i<DIM;i++){
    n_peaks[i]=0;
    n_cavities[i]=0;
    f_prev_peak[i]=0.0;
    f_prev_cavity[i]=0.0;
    for(j=0;j<MAX_N_PEAKS;j++){
      t_peak[j][i]=0.0;
      x_peak[j][i]=0.0;
      t_cavity[j][i]=0.0;
      x_cavity[j][i]=0.0;
    }
  }
  for(i=1;i<write_count;i++){
    func_odeiv(ts[i],xs[i],f,&mu);
    for(j=0;j<DIM;j++){
      /*Determining peaks...*/
      if(f[j]>0)
	f_prev_peak[j]=f[j];
      if((f[j]<0) && (f_prev_peak[j]>0)){
	if(xs[i-1][j]<xs[i][j]){
	  t_peak[n_peaks[j]][j]=ts[i];
	  x_peak[n_peaks[j]][j]=xs[i][j];
	}
	if(xs[i-1][j]>xs[i][j]){
	  t_peak[n_peaks[j]][j]=ts[i-1];
	  x_peak[n_peaks[j]][j]=xs[i-1][j];
	}
	f_prev_peak[j]=f[j];
	n_peaks[j]++;
      }
      /*Determining cavities...*/
      if(f[j]<0)
	f_prev_cavity[j]=f[j];
      if((f[j]>0) && (f_prev_cavity[j]<0)){
	if(xs[i-1][j]<xs[i][j]){
	  t_cavity[n_cavities[j]][j]=ts[i];
	  x_cavity[n_cavities[j]][j]=xs[i][j];
	}
	if(xs[i-1][j]>xs[i][j]){
	  t_cavity[n_cavities[j]][j]=ts[i-1];
	  x_cavity[n_cavities[j]][j]=xs[i-1][j];
	}
	f_prev_cavity[j]=f[j];
	n_cavities[j]++;
      }
      /*Determine whether we are on the L.C.*/
      /*There is NO possibility of two peaks(or two cavities) in a row.*/
      if((n_peaks[j]==n_cavities[j]) && (n_peaks[j]>0)){
	/*It is a little bit expensive because amplitude(...) again
	  iterates over DIM ...*/
	amplitude(local_ampl,t_cavity[n_cavities[j]-1][j],t_peak[n_peaks[j]-1][j]);
	/*Check if local amplitude for at least one of the variables
	  is less than the error, then return.*/
	if(local_ampl[j] < eps_abs_am){
	  fprintf(stdout,"Stopped computing peaks\n");
	  return 99;
	}
      }
    }
  }
  
  return 0;
}

/* int peak_def(int n_peaks[],double *t_peak, double *x_peak,int (*func)(double,const double [],double [], void *)) */
/* { */
/*   int i,j; */
/*   double static x_prev[MAXDIM]; */
/*   double static f_prev[MAXDIM]; */
/*   double static t_prev[MAXDIM]; */

/*   func(t,x,f,&mu); */

/*   j=1; */
/*   while((n_steps == 0) && (j <= DIM)){ */
/*     n_peaks[j-1]=0; */
/*     f_prev[j-1] = 0; */
/*     x_prev[j-1] = 0; */
/*     t_prev[j-1] = 0; */
/*     for(i=0;i<MAX_N_PEAKS;i++){ /\* We set to zero these values: for MAXDIM*MAX_N_PEAKS elements *\/ */
/*       *(t_peak+i*MAXDIM+j-1)=0.0;       /\* (not for DIM*MAX_N_PEAKS) *\/ */
/*       *(x_peak+i*MAXDIM+j-1)=0.0; */
/*     } */
/*     j++; */
/*   } */
/*   for(i=1; i <= DIM; i++) { */
/*     if(*(f+i-1) > 0) { */
/*       f_prev[i-1] = *(f+i-1); */
/*       x_prev[i-1] = x[i-1]; */
/*       t_prev[i-1] = t; */
/*     } */
/*     if( *(f+i-1) == 0 ) { */
/*       *(t_peak+n_peaks[i-1]*MAXDIM + i-1) = t; */
/*       *(x_peak+n_peaks[i-1]*MAXDIM + i-1) = x[i-1]; */
/*       fprintf(stdout,"BINGO\n"); */
/*       n_peaks[i-1]+=1; */
/*     } */
/*     if((*(f+i-1) < 0) && (f_prev[i-1] > 0))  { */
/*       if( x_prev[i-1] < x[i-1] ) { */
/* 	*(t_peak+n_peaks[i-1]*MAXDIM + i-1) = t; */
/* 	*(x_peak+n_peaks[i-1]*MAXDIM + i-1) = x[i-1]; */
/*       }  */
/*       if( x_prev[i-1] >= x[i-1] ){ */
/* 	*(t_peak+n_peaks[i-1]*MAXDIM + i-1) = t_prev[i-1]; */
/* 	*(x_peak+n_peaks[i-1]*MAXDIM + i-1) = x_prev[i-1]; */
/*       } */
/*       f_prev[i-1] = *(f+i-1); */
/*       n_peaks[i-1]++; */
/*       if(n_peaks[i-1] >= MAX_N_PEAKS &&				\ */
/* 	 (max_x[i-1]-min_x[i-1])>D(eps_abs_am,eps_rel_am,(float)0.0)){ */
/* 	fprintf(stdout,"Too many peaks\n"); */
/* 	return -1; */
/*       } */
/*       else if(((max_x[i-1]-min_x[i-1])<D(eps_abs_am,eps_rel_am,(float)0.0)) \ */
/* 	      && (n_steps != 0)){ */
/* 	return -2; */
/*       } */
		 
/*     } */
/*   } */
/*   return 0; */
/* } */

int amplitude(double * ampl, double t0, double t1)
{
  /*Temporal values of max and min*/
  double maxX[DIM];
  double minX[DIM];
  /*Get max and min from desired interval (t0,t1)*/
  max_min_x(maxX,minX,t0,t1);
  /***********************************************/
  int i;
  for(i=0; i<DIM; i++)
    *(ampl + i) = 0;
  for(i=0; i<DIM; i++)
    *(ampl + i) = maxX[i] - minX[i];

  return 0;
}
