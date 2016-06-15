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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/

/* This file defines/describes the functions and routines used in plotting procedures
   of Dinamica. */

#include "init.h"
//#include <plot.h>
#include <stdio.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>
#define PLOT_NUM_MIN_CHAR 10//Num of char in strings determining numbers in plots
#define USER_COORD_EXTEN 0.2//fraction of the user coord to extend from the ends
#define XY_AXES_EXTEN 0.05//fraction to extend xy-axes from the +\infty ends
/* -\infty ends are left as min of the corresponding x and y vars */

void gplot_frames(const char *data, const int n_frames, const int offset, const int varInd)
{
  /* This function puts the time series(frame(s)) to the graphics output.
   INPUT:
   data - data file name
   n_frames - # frames to put
   offset - # of frames to skip from the beginning of the data file
   varInd - index of the variable to print
  */
  int i;
  /* Set the labels and the terminal */
  gnuplot_cmd(plot_handle,"set term x11 0");
  gnuplot_cmd(plot_handle,"set xlabel 'Time'");
  gnuplot_cmd(plot_handle,"set ylabel '%s'",var_name[varInd]);
  gnuplot_cmd(plot_handle,"plot '%s' i %d u 1:%d w l t 'frame-1'",data,offset,varInd+2);
  for(i=offset+1;i<n_frames;i++){
    gnuplot_cmd(plot_handle,"replot '%s' i %d u 1:%d w l t 'frame-%d'",data,i,varInd+2,i+1);
  }
}

int gnuplot_interp()
{/* This function differs from all _interp function from input_interpreter.c. It
    emulates the gnuplot command line. */
  int Len = 100;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  while(1){
    fprintf(stdout,"%s>",gnuplot_prompt);
    fgets(symbuf,Len,stdin);
    if(symbuf[0] == 'q'){
      fprintf(stdout,"Exit gnuplot prompt...\n");
      break;
    }
    if(strchr(symbuf,'\n') == 0){
      fprintf(stdout,"Limit one line input to %d characters.",Len);
      fprintf(stdout,"Use \"\" for continuation.\n");
      continue;
    }
    gnuplot_cmd(plot_handle,"%s",symbuf);
  }
  return 0;
}

void gplot_results(const char *data_name, const int complex_flag,
		   const char *method)
{
  int i;
  char *ma_name = (char *)malloc((strlen(data_name)+4)*sizeof(char));
  if(ma_name == NULL){
    fprintf(stderr,"Error(gplot_results): malloc failed\n");
    return;
  }
  /* Define the name of MA trajectory */
  ma_name = strcpy(ma_name,data_name);
  ma_name = strcat(ma_name,".ma");

  /* Error check */
  if((graph.xInd > (DIM-1)) || (graph.yInd[0] > (DIM-1)) ||
     (graph.yInd[1] > (DIM-1)) || (graph.yInd[2] > (DIM-1))){
    fprintf(stderr,"Graph: out of var indices. Change to default.\n");
    graph.xInd = -1;
    graph.yInd[0] = 0;
    graph.yInd[1] = -2;
    graph.yInd[2] = -2;
  }
  graph_set_labels(graph,"tsphase");
  /* Grid */
  if(graph.grid_flag)
    gnuplot_cmd(plot_handle,"set grid\n");

  if((complex_flag == 0) || (complex_flag == 1)){
    if((strcmp(method,"discrete") == 0) && (graph.xInd == -1)){
      gnuplot_cmd(plot_handle,"plot '%s' i 0 u %d:%d w steps t \"%s(%s)\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method,
		  var_name[graph.yInd[0]]);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w steps t \"%s(%s)\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method,
		      var_name[graph.yInd[i]]);
      	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w steps t \"%s-MA(%s)\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method,
		    var_name[graph.yInd[0]]);
    }
    else if(graph.xInd != -1){
      gnuplot_cmd(plot_handle,"plot '%s' i 0 u %d:%d w l t \"%s\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w l t \"%s\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method);
      	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w l t \"%s-MA\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method);
    }
    else{
      gnuplot_cmd(plot_handle,"plot '%s' i 0 u %d:%d w l t \"%s(%s)\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method,
		  var_name[graph.yInd[0]]);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w l t \"%s(%s)\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method,
		      var_name[graph.yInd[i]]);
	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 0 u %d:%d w l t \"%s-MA(%s)\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method,
		    var_name[graph.yInd[0]]);
    }
  }
  else{/* Second run fo complex method */
    if((strcmp(method,"discrete") == 0) && (graph.xInd == -1)){
      gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w steps t \"%s(%s)\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method,
		  var_name[graph.yInd[0]]);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w steps t \"%s(%s)\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method,
		      var_name[graph.yInd[i]]);
      	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w steps t \"%s-MA(%s)\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method,
		    var_name[graph.yInd[0]]);
    }
    else if(graph.xInd != -1){/* Phase portrait */
      gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method);
      	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s-MA\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method);
    }
    else{
      gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s(%s)\"\n",
		  data_name,graph.xInd+2,graph.yInd[0]+2,method,
		  var_name[graph.yInd[0]]);
      i = 1;
      while(i<3){
	if(graph.yInd[i] > -1)
	  gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s(%s)\"",
		      data_name,graph.xInd+2,graph.yInd[i]+2,method,
		      var_name[graph.yInd[i]]);
      	i++;
      }
      if(ma_span)
	gnuplot_cmd(plot_handle,"replot '%s' i 1 u %d:%d w l t \"%s-MA(%s)\"",
		    ma_name,graph.xInd+2,graph.yInd[0]+2,method,
		    var_name[graph.yInd[0]]);
    }
  }
}

void graph_set_labels(const struct coordNet graph, const char *type)
{
  char *tmp = (char *)calloc(20,sizeof(char));
  int i = 0;
  if(strcmp(type,"tsphase")==0){
    gnuplot_cmd(plot_handle,"set term x11 0 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    if(graph.xInd == -1)
      gnuplot_cmd(plot_handle,"set xlabel \"Time\"\n");
    else
      gnuplot_cmd(plot_handle,"set xlabel \"%s\"\n",var_name[graph.xInd]);
    while(i<3){
      if(graph.yInd[i] > -1){
	strcat(tmp,var_name[graph.yInd[i]]);
	strcat(tmp,",");
      }
      i++;
    }
    gnuplot_cmd(plot_handle,"set ylabel \"%s\"\n",tmp);
  }
  if(strcmp(type,"tmap")==0){
    gnuplot_cmd(plot_handle,"set term x11 1 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    gnuplot_cmd(plot_handle,"set xlabel \"T(n)\"\n");
    gnuplot_cmd(plot_handle,"set ylabel \"T(n+1)\"\n");
    gnuplot_cmd(plot_handle,"set title \"n=%d\"\n",nPerStoch);
  }
  if(strcmp(type,"ac")==0){
    gnuplot_cmd(plot_handle,"set term x11 1 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    gnuplot_cmd(plot_handle,"set xlabel \"Time points\"\n");
    gnuplot_cmd(plot_handle,"set ylabel \"Autocorrelation\"\n");
    gnuplot_cmd(plot_handle,"set title \"Autocorrelation function\"\n");
  }
  if(strcmp(type,"thist")==0){
    gnuplot_cmd(plot_handle,"set term x11 1 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    gnuplot_cmd(plot_handle,"set boxwidth 0.8 relative\n");
    gnuplot_cmd(plot_handle,"set style fill solid\n");
    gnuplot_cmd(plot_handle,"set style line 1 lw 2 lc rgb \"green\"\n");
    gnuplot_cmd(plot_handle,"set xlabel \"Period\"\n");
    gnuplot_cmd(plot_handle,"set ylabel\n");
    gnuplot_cmd(plot_handle,"set title \"n=%d,m=%G,sd=%G,CV^2=%G,F=%G\"\n",
		nPerStoch,
		gsl_stats_mean(perStoch,1,nPerStoch),
		gsl_stats_sd(perStoch,1,nPerStoch),
		pow(gsl_stats_sd(perStoch,1,nPerStoch)
		    /gsl_stats_mean(perStoch,1,nPerStoch),2),
		pow(gsl_stats_sd(perStoch,1,nPerStoch),2)
		/gsl_stats_mean(perStoch,1,nPerStoch));
  }
  if(strcmp(type,"mdist")==0){
    gnuplot_cmd(plot_handle,"set term x11 1 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    gnuplot_cmd(plot_handle,"set boxwidth 0.8 relative\n");
    gnuplot_cmd(plot_handle,"set style fill solid\n");
    gnuplot_cmd(plot_handle,"set style line 1 lw 2 lc rgb \"green\"\n");
    gnuplot_cmd(plot_handle,"set xlabel \"# %s species\"\n",var_name[perVarInd]);
    gnuplot_cmd(plot_handle,"set ylabel\n");
  }
  if(strcmp(type,"ampldist")==0){
    gnuplot_cmd(plot_handle,"set term x11 1 enh\n");
    gnuplot_cmd(plot_handle,"reset\n");
    gnuplot_cmd(plot_handle,"set boxwidth 0.8 relative\n");
    gnuplot_cmd(plot_handle,"set style fill solid\n");
    gnuplot_cmd(plot_handle,"set style line 1 lw 2 lc rgb \"green\"\n");
    gnuplot_cmd(plot_handle,"set xlabel \"Amplitudes\"\n");
    gnuplot_cmd(plot_handle,"set ylabel\n");
  }
  free(tmp);
}

void init_graph()
{
  /* graph.xInd = -1; Initialized within read_ode.c*/
  /* graph.yInd[0-2] = 0; Initialized within read_ode.c*/
  graph.xUlim = 20;
  graph.xLlim = 0;
  graph.yUlim = 20;
  graph.yLlim = 0;
  graph.nTics = 5;
  graph.grid_flag = 0;
  plot_handle = gnuplot_init();
  gnuplot_cmd(plot_handle,"set pointsize 2\n");
}

void send_to_eps(const char *fname)
{
  gnuplot_cmd(plot_handle,"set term postscript eps color enh\n");
  gnuplot_cmd(plot_handle,"set output \"%s\"\n",fname);
  gnuplot_cmd(plot_handle,"replot\n");
  gnuplot_cmd(plot_handle,"set output\n");
  gnuplot_cmd(plot_handle,"set term X11\n");//Reset to X11 back
  gnuplot_cmd(plot_handle,"set pointsize 2\n");
  fprintf(stdout,"Sent to %s\n",fname);
}
