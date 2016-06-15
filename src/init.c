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
   Moscow, Russia / Tampere, Finland*/
/****************************************************************************************/

#include "init.h"
#include "errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <plot.h>


int init_command_line()
{
  //for initials, at start up or by show_initials().
  init_name = (char *)calloc(FNAME_SIZE,sizeof(char));
  // for writing trajectory, data points.
  data_name = (char *)malloc(FNAME_SIZE*sizeof(char));
  strcpy(data_name,"dat");/* Default data file name */
  /* Special affix for the output file, containing common entities, computed for every
     trajectory. */
  out_name = (char *)malloc(FNAME_SIZE*sizeof(char));
  strcpy(out_name,"_out");/* Default affix for the output file name */
  /*Configuration file.*/
  conf_name = (char *)malloc(FNAME_SIZE*sizeof(char));//Config file(read_config.c)
  /* Input script file: must be init to zeros */
  input_name = (char *)calloc(FNAME_SIZE,sizeof(char));

  return 0;
}

int init()
{
  int i,j,k;
  /******************************************************
   * Variables assignments (ALL INITIAL VALUES)
   * ****************************************************/
  p_mynum = &mynum;
  p_mu = &mu;
  /***********************
   * Menu prompts
   * *********************/
  main_prompt = "";
  numerics_prompt = "numerics";
  graphics_prompt = "graphics";
  cross_prompt = "cross";
  periodics_prompt = "periodics";
  file_prompt = "file";
  cont_prompt = "continue";
  rand_prompt = "random";
  errors_prompt = "errors";
  sing_prompt = "singularity";
  traj_prompt = "trajectory";
  /* *********************** */
  /* Submenu prompts */
  /* *********************** */
  /* Histogram submenu of periodics menu */
  periodics_hist_prompt = (char *)malloc((strlen(periodics_prompt) +		\
				     strlen("/histogram")+2)*sizeof(char));
  periodics_hist_prompt = strcpy(periodics_hist_prompt,periodics_prompt);
  periodics_hist_prompt = strcat(periodics_hist_prompt,"/histogram");
  gnuplot_prompt = "gnuplot('q' to exit)";
  /* End of histogram submenu */
  /* Menu buffers */
  /* History buffers */
  buf_hist = (char ***)malloc(MAX_N_HIST_ENT*sizeof(char **));
  for(i=0;i<MAX_N_HIST_ENT;i++){
    buf_hist[i] = (char **)malloc(MAX_N_ARG*sizeof(char *));
    for(j=0;j<MAX_N_ARG;j++)
      buf_hist[i][j] = (char *)calloc(MAX_ARG_LEN,sizeof(char));
  }
  /* Commands and their arguments. See read_menu(...) function. */
  cmd = (char **)malloc(MAX_N_ARG*sizeof(char *));
  for(i=0;i<MAX_N_ARG;i++)/* We allocate memory for all possible arguments */
    cmd[i] = (char *)calloc(MAX_ARG_LEN,sizeof(char));
  /*Loading constants and parameters from the file --- conf_name. */
  if(strlen(conf_name) != 0) 
    if((read_conf(conf_name)) == -1)/*read_config is the binary conf file~(.bcf)*/
      fprintf(stdout,"Config file: could not open the file.\n");
  /*************************************/
  /*Converting new style entities to old one for compatibility*/
  if((fmod((double)DIM,(double)LDIM))>0.5){
    /*If reminder of DIM/LDIM == 1 or more than it is wrongly defined
      system.*/
    /*DIM comes from number of ODEs.
      LDIM is specified by the user in #system directive of ode
      file.*/
    /*Error in ODE system specification*/
    fprintf(stderr,"Fatal: not symmetric ODE system\n");
    return 151;
  }
  /*New style -> old style converter*/
  mynum.total_dim = LDIM;
  mynum.local_dim = (int)DIM/LDIM;
  mynum.num_par = PARDIM;
  if((strlen(method) >= METH_NAME_LEN) || (strlen(method2) >= METH_NAME_LEN)){
    fprintf(stderr,"Error: method has a very long name.\n");
    fprintf(stderr,"Max allowed length is %d.\n",METH_NAME_LEN-1);
    return 143;
  }
  strncpy(mynum.method,method,(size_t)strlen(method));
  strncpy(mynum.method2,method2,(size_t)strlen(method2));
  /********************************************
  Initializing values
  *********************************************/
  mynum.global_buffer = BUFFER;/*Reading default*/
  xs = (double **)calloc(BUFFER,sizeof(double *));
  for(i=0;i<BUFFER;i++){
    xs[i] = (double *)calloc(DIM,sizeof(double));
  }
  ts = (double *)calloc(BUFFER,sizeof(double));
  //ampl = (double *)calloc(DIM,sizeof(double));
  //max_x = (double *)calloc(DIM,sizeof(double));
  //min_x = (double *)calloc(DIM,sizeof(double));
  x_cross = (double *)calloc(DIM,sizeof(double));
  /* Vars */
  y = (long int *)malloc(DIM*sizeof(long int));
  x = (double *)malloc(DIM*sizeof(double));
  f = (double *)malloc(DIM*sizeof(double));
  for(i=0;i<DIM;i++){
    x[i] = xin[i];
    y[i] = yin[i];
  }
  ksi = (double *)malloc(DIM*sizeof(double));/*Noise is zero*/
  /* Propensities */
  pr = (double *)malloc(NREAC*sizeof(double));
  /* Auxillaries */
  a = (double *)malloc(AUXNUM*sizeof(double));
  /************************************************************/
  /*RANDOM number parameters initialization...*/
  random_init();
  /*CONTINUATION parameters intitialization...*/
  cont_init();
  /*ROOT FINDING parameters initialization...*/
  multiroot_init();
  /*Lyapunov parameters initialization...*/
  lyap_init();
  /* Graphics init */
  init_graph();
  /* periods histogram init */
  thist_init();
  /* T-system init */
  traj_init();
  //******************************************************
  // Loading init. conditions from the file -- init_name.
  if(strlen(init_name) != 0){
    FILE *init;
    init = fopen(init_name,"r");
    if(init == NULL){
      fprintf(stdout,"Cannot open file with initial conditions.\n");}
    else{
      j = 1;
      while(j<=DIM) {
	fscanf(init,"%lf",&x[j-1]);
	y[j-1] = (long int)x[j-1];/*discrete*/
	j += 1;
      }
      fprintf(stdout,"Initial points\n");
      j = 1;
      while(j<=DIM){
	xin[j-1] = x[j-1];
	yin[j-1] = y[j-1];
	fprintf(stdout,"%lf\n",xin[j-1]);
	j += 1;
      }
      fclose(init);
    }
  }

  return 0;
}

void din_close()
{/* The function closes everything that needs to be closed */
  int i,j;
  free(init_name);
  free(data_name);
  free(conf_name);
  free(input_name);
  free(out_name);
  for(j=0;j<MAX_N_HIST_ENT;j++){
    for(i=0;i<MAX_N_ARG;i++)
      free(buf_hist[j][i]);
    free(buf_hist[j]);
  }
  free(buf_hist);
  for(i=0;i<MAX_N_ARG;i++)
    free(cmd[i]);
  free(cmd);
  /* The following is deprecated */
  for(i=0;i<BUFFER;i++)
    free(xs[i]);
  free(xs);
  free(ts);
  /* *************************** */
  //free(ampl);
  //free(max_x);
  //free(min_x);
  free(x_cross);
  free(x);
  free(y);
  free(xin);
  free(xin_par);
  free(yin);
  free(f);
  free(ksi);
  free(pr);
  free(a);
  for(i=0;i<DIM;i++)
    free(var_name[i]);
  free(var_name);
  for(i=0;i<AUXNUM;i++)
    free(aux_name[i]);
  free(aux_name);
  gnuplot_close(plot_handle);
  if(!perDet)
    free(perDet);
  if(!perStoch)
    free(perStoch);
  /* Free Random */
  random_free();
  /*Free Lyapunov*/
  lyap_free();
  /* Free periods histogram */
  thist_free();
  /* Deleting unnecessary files after finishing the program */
  printf("Deleting unnecessary files...\n");
  if((fopen("slopeAmpl.dat","r")) != NULL){
    if(remove("slopeAmpl.dat") != 0)
      printf("Error deleting file slopeAmpl.dat\n");
  }
  if((fopen("slopeAmpl1.dat","r")) != NULL){
    if(remove("slopeAmpl1.dat") != 0)
      printf("Error deleting file slopeAmpl.dat\n");
  }
  if((fopen("hist.dat","r")) != NULL){
    if(remove("hist.dat") != 0)
      printf("Error deleting file hist.dat\n");
  }
  if((fopen("acorr.dat","r")) != NULL){
    if(remove("acorr.dat") != 0)
      printf("Error deleting file acorr.dat\n");
  }
  if((fopen("tmap.dat","r")) != NULL){
    if(remove("tmap.dat") != 0)
      printf("Error deleting file tmap.dat\n");
  }
  if((fopen("rand.throw","r")) != NULL){
    if(remove("rand.throw") != 0)
      printf("Error deleting file rand.throw\n");
  }
  printf("Done.\n");
}

int isVar(char const *name,char **var_list,int var_list_size)
{
  int i=0;
  if(var_list_size <= 0)/* No variables in the system! */
    return -10;
  while(i < var_list_size){
    if((strcmp(name,var_list[i]))==0)
      break;
    i++;
  }
  if(i == var_list_size)
    return -1;
  
  return i;/*We return var index as we found the condition */
}

int isPar(char const *name, const int npar)
{
  int i=0;
  if(npar <= 0)/* no parameters in the system */
    return -10;
  while((strcmp(name,par_name[i]))!=0){
    i++;
    if(i == npar){
      return -1;/* We did not meet the condition */
    }
  }
  return i;/*We return par index as we found the condition */
}
