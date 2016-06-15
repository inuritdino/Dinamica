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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland*/
/****************************************************************************************/
/* It is the core integration algorithm of DINAMICA. It uses four
   methods: my own Runge-Kutta 4-th order without any modifications,
   Runge-Kutta-Fehlberg (4,5), Runge-Kutta Cash-Karp (4,5) and
   Runge-Kutta Prince-Dormand (8,9). For each integrating method there
   is comparising test (if strcmp...) for string in the mynum.method
   struct member. The only function run() uses three parameters: total
   time of integration, write flag and transient flag. 
   Total time is value of upper limit of independent variable till
   what we integrate.
   Write flag(wr) is flag for whether we should write output
   trajectory file.
   Transient flag(tr) is flag for whether our integration is
   transient, e.g. is used for getting to attractor in slow
   systems. Transient integration is always free of writing output
   trajectory file.
   Common integration may be either with writing to disk or without
   it.*/  

#include "init.h"
#include "errors.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>/* calculating the cpu time for the integration */
#include <gsl/gsl_math.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_statistics_double.h>
//#include <plot.h>

int run_lin()
{
  int i,j,k;
  double t1 = mynum.total_time;
  double h = mynum.step;
  double x_lin[DIM];
  int npr;
  t = 0.0;
  for(i=0;i<DIM;i++)
    x_lin[i]=x[i];
  npr = 1;
  while(t <= t1)
    {
      if(t<h){
	printf("Initials:\n");
	for(i=0;i<DIM;i++)
	  printf("%G ",x[i]);
	printf("\n");
	for(i=0;i<DIM;i++)
	  printf("%G ",x_lin[i]);
	printf("\n");
	printf("-----------\n");
      }
      if(t>(t1/10)*npr){
	for(i=0;i<DIM;i++)
	  printf("%G ",x[i]);
	printf("\n");
	for(i=0;i<DIM;i++)
	  printf("%G ",x_lin[i]);
	printf("\n");
	printf("-----------\n");
	npr++;
      }
      solver_lin(x,x_lin,f,h,ode_lin);
      solver(x,h,func_odeiv);
      
      t = t+h;
    }

  return 0;
}

int run (const char *method, const double time, const int short wr,
	 const int short tr, int nruns)
{
  clock_t start, end;/* Start and End ticks of the processes */
  double cpu_time_used;/* The CPU time used */
  int i,j,k,n_wr_step,outPrint = 1;
  /*** If <tr> is set, <nruns> must be 1 ***/
  if(tr)
    nruns = 1;
  /* Static flag showing what run we are in. 0 -- not complex OR before the first run
     and after the second one of complex, 1 -- complex first run and 2 --
     complex second run.*/
  int short static complex_run_flag = 0;
  //  printf("complex flag = %d\n",complex_run_flag);
  /*** Data file ***/
  FILE *data;
  if(strlen(data_name) == 0) strcpy(data_name,"dat");
  if(wr && !tr){/* We need to open file */
    /* if(strcmp(method,"complex") != 0){ */
    /*   /\* simple method: do not open file for complex one. *\/ */
    /*   if(complex_run_flag != 2) */
    /* 	data = fopen(data_name,"w"); */
    /*   else */
    /* 	data = fopen(data_name,"a"); */
    /* } */
    if(!complex_run_flag)
      data = fopen(data_name,"w");
    else
      data = fopen(data_name,"a");
  }
  /*** INFO ARRAY: Technical details of the timeseries to the file ***/
  if(wr && !tr){
    if(!complex_run_flag){
      if(strcmp(method,"complex")==0){
	fprintf(data,"#1,%d,0,%d\n",nruns,lang_flag);
      }
      else{
	fprintf(data,"#0,%d,",nruns);
	if(strcmp(method,"discrete")==0 || lang_flag)
	  fprintf(data,"0,%d\n",lang_flag);/* not deterministic */
	else
	  fprintf(data,"1,%d\n",lang_flag);/* Pure deterministic flag */
      }
      fprintf(data,"#");
      for(i=0;i<DIM;i++)
	fprintf(data,"U(%d)=%s,",i+1,var_name[i]);
      fprintf(data,"\n#");
      for(i=0;i<PARDIM;i++)
	fprintf(data,"%s=%G,",par_name[i],mu.P[i]);
      fprintf(data,"\n");
      if(strcmp(method,"complex")==0)/* Close the stream if complex */
	fclose(data);
    }
  }
  /*** Set initial time, time step etc. ***/
  t = 0.0;
  double t1 = time;
  double h = mynum.step;
  double tau;/*Next reaction time(Gill)*/
  int mur;/*Next reaction(Gill)*/
  double A_sum;/*Sum of all reactions' propensities(Gill)*/
  double y_dbl[DIM];/*Double version of y(Gill)*/
  double y_in[DIM];/*Initial value for each run(swap variable)*/
  n_steps = 0;
  write_count = 0;
  double tmp,tmp1;
  /* Array of random values from standard Gauss distribution (see solver_lang) */
  double Nstd[DIM];
  /*** Random number generation initialization. ***/
  if(strcmp(method,"discrete") == 0 || lang_flag){
    /* Do NOT initialize for complex method */
    if(seed_flag == 2)
      rng_seed = generate_seed();
    rng = init_rng(rng_seed,rng_type,rng);/* rng_type and rng defined in random.h */
  }
  /*** Warning about multiple runs for the deterministic method ***/
  if((strcmp(method,"complex")!=0) && (nruns>1) &&	\
     (strcmp(method,"discrete") != 0) && (!lang_flag)){
    fprintf(stderr,"Method must be <stochastic>\n");
    fprintf(stderr,"If you want to run deterministic simulations multiple");
    fprintf(stderr," times,\nuse transient runs instead (see `rt').\n");
    return 45;
  }
  /*** COMPLEX METHOD: recursive calls to run(...) ***/
  /* Complex method: one deterministic simulation plus nruns of */
  /*  stochastic; requires to specify two methods respectively */
  if((strcmp(method,"complex") == 0)){
    printf("Run deterministic method: %s...\n",mynum.method2);
    complex_run_flag = 1;
    /* First run deterministic simulation */
    if(lang_flag == 0)
      run(method2,time,wr,tr,1);
    else{
      lang_flag = 0;
      run(method2,time,wr,tr,1);
      lang_flag = 1;
    }
    complex_run_flag = 2;
    /* Then, run stochastic method nruns times */
    if(lang_flag == 0){/* we have discrete method, this was checked */
      printf("Run stochastic method: discrete...\n");
      run("discrete",time,wr,tr,nruns);
    }
    else{
      printf("Run stochastic method: langevin...\n");
      run(mynum.method2,time,wr,tr,nruns);
    }
    /* Finished complex run */
    complex_run_flag = 0;
  }
  /*** Own routine for Gillespie algorithm ***/
  else if((strcmp(method,"discrete")==0)){/*Gillespie algorithm*/
    start = clock();/* Start the counter */
    for(i=0;i<DIM;i++)/* Initial value stored in the swap */
      y_in[i] = y[i];
    for(k=0;k<nruns;k++){
      /*Init gillespie*/
      t = 0.0;
      n_wr_step = 0;
      for(i=0;i<DIM;i++){/*Doublify and initialize*/
	y[i] = y_in[i];
	y_dbl[i] = (double)y[i];
      }
      gill_init(y,pr,&A_sum);
      while (t < t1){
	/*Compute auxillary...*/
	auxillary(t,mu.P,y_dbl,a);
	/*Generate two random numbers*/
	tmp1 = gsl_rng_uniform_pos(rng);/*R1*/
	tau = (1/A_sum)*log(1/tmp1);/* Time of the next reaction */
	tmp = gsl_rng_uniform_pos(rng);/*R2*/
	tmp1 = A_sum;
	for(i=NREAC;i>0;i--){/* Choosing the next reaction */
	  tmp1 = tmp1-pr[i-1];
	  if((tmp*A_sum) > tmp1){
	    mur = i;
	    break;
	  }
	}
	if(!tr){/* write the output to file */
	  if(wr){
	    while((t+tau) >= (n_wr_step*mynum.smp_frq)){
	      fprintf(data,"%G ",(double)n_wr_step*mynum.smp_frq);
	      for(j=0;j<DIM;j++)
		fprintf(data,"%ld ",y[j]);
	      fprintf(data,"\n");
	      n_wr_step++;
	    }
	  }
	}
	t = t + tau;/*Forward in time*/
	update(y,mur);/*Update system according to the reaction*/
	/*Recompute sum of the propensities*/
	A_sum = 0.0;
	for(i=0;i<DIM;i++)/*Doublify*/
	  y_dbl[i]=(double)y[i];
	propensity(y_dbl,mu.P,pr);/* Recompute propensities, all of them */
	for(i=0;i<NREAC;i++)
	  A_sum = A_sum + pr[i];
      }
      if(!tr)
	fprintf(data,"\n\n");/*Separating data in the file*/
      if(nruns >= 5){
	i = (int)(outPrint*nruns/5);
	if(k == i){
	  outPrint++;
	  printf("Completed: %d/%d\n",k+1,nruns);
	}
      }
    }
    end = clock();
  }
  /* Fixed step integrators realized in solver and solver_lang functions. These are
     also applied to the stochastic Ito-type differential equations. */
  else if((!strcmp(method,"run-kut4")) || \
	  (!strcmp(method,"eu")) || (!strcmp(method,"milst"))) {
    start = clock();
    for(i=0;i<DIM;i++)/* Initial value stored in the swap */
      y_in[i] = x[i];
    for(k=0;k<nruns;k++){
      t = 0.0;
      n_steps = 0;
      write_count = 0;
      for(i=0;i<DIM;i++)/*initialize: take from the swap*/
	x[i] = y_in[i];
      while(t < t1) {
	/*Compute auxillary...*/
	auxillary(t,mu.P,x,a);
	if(!tr){
	  if(n_steps == write_count*mynum.write_step) {
	    /* if(write_count < BUFFER)//mynum.global_buffer */
	    /*   write_data_array(); */
	    /* if(buffer_overfull_check()) break; */
	    if(wr){
	      //i = write_data_file(data);
	      if(write_data_file(data))
		break;
	    }
	    write_count++;
	  }
	}
	if(!lang_flag){
	  if( (solver(x,h,func_odeiv)) )
	    return 100;//Solver for run-kut4 and Euler, own.
	}
	else{
	  for(i=0;i<DIM;i++)
	    Nstd[i] = gsl_ran_ugaussian(rng);
	  if( (solver_lang(x,h,func_odeiv,Nstd)) )
	    return 100;
	}
	t = t+h;
	n_steps++;
      }
      if(!tr)
	fprintf(data,"\n\n");/*Separating data in the file*/
      if(nruns >= 5){
	i = (int)(outPrint*nruns/5);
	if(k == i){
	  outPrint++;
	  printf("Completed: %d/%d\n",k+1,nruns);
	}
      }
    }
    end = clock();
  }
  else{
    start = clock();
    gsl_odeiv2_step *s;
    if(strcmp(method,"rkf45") == 0)
      s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45,DIM);
    if(strcmp(method,"rk8pd") == 0)
      s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rk8pd,DIM);
    if(strcmp(method,"rkck") == 0)
      s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkck,DIM);
    if(strcmp(method,"bsimp") == 0)
      s = gsl_odeiv2_step_alloc(gsl_odeiv2_step_bsimp,DIM);
    gsl_odeiv2_control *c = gsl_odeiv2_control_standard_new (eps_abs_int,eps_rel_int,a_y,a_dydt);
    gsl_odeiv2_evolve *e = gsl_odeiv2_evolve_alloc (DIM);
    gsl_odeiv2_system sys = {func_odeiv,jac,DIM,&mu};
    
    for(i=0;i<DIM;i++)/* Initial value stored in the swap */
      y_in[i] = x[i];
    for(k=0;k<nruns;k++){
      /* This multiple runs for the deterministic methods might be used only for
      Langevin simulation for now. */
      /* if(k>0) */
      /* 	gsl_odeiv2_step_reset(s);/\* Resetting the stepping function *\/ */
      t = 0.0;
      n_steps = 0;
      write_count = 0;
      for(i=0;i<DIM;i++)/*initialize: take from the swap*/
	x[i] = y_in[i];
      while(t < t1) {
	/*Compute auxillary...*/
	auxillary(t,mu.P,x,a);
	if(!tr){
	  if(n_steps == write_count*mynum.write_step) {
	    /* if(write_count < BUFFER)//mynum.global_buffer */
	    /*   write_data_array(); */
    	    /* if(buffer_overfull_check()) break; */
	    if(wr){
	      //i = write_data_file(data);
	      if(write_data_file(data))
		break;
	    }
	    write_count++;
	  }
	}
	int status = gsl_odeiv2_evolve_apply (e, c, s,
					     &sys, &t, t1, &h, x);
	if(status != GSL_SUCCESS)
	  break;
	n_steps++;
      }//while (t<t1)
      if(!tr)
	fprintf(data,"\n\n");/*Separating data in the file*/
      if(nruns >= 5){
	i = (int)(outPrint*nruns/5);
	if(k == i){
	  outPrint++;
	  printf("Completed: %d/%d\n",k+1,nruns);
	}
      }
      /* Reset the stepping function */
      gsl_odeiv2_step_reset(s);
    }//k over nruns
    //mynum.global_buffer = BUFFER;
    gsl_odeiv2_evolve_free(e);
    gsl_odeiv2_control_free(c);
    gsl_odeiv2_step_free(s);

    end = clock();
  }
  /***Closing data file stream.***/
  if(wr && !tr && (strcmp(method,"complex")!=0))
    fclose(data);
  /***Free random generator***/
  if(strcmp(method,"discrete") == 0 || lang_flag)
    gsl_rng_free(rng);
  /***Compute MA***/
  char tmp_name[strlen(data_name)+4];// = (char *)malloc((strlen(data_name)+4)*sizeof(char));  
  if(!complex_run_flag && !tr && wr){
    if(ma_span){
      moving_average(data_name,ma_span);/* Calculate MA */
      strncpy(tmp_name,data_name,strlen(data_name));
      strncat(tmp_name,".ma",3);
    }
  }
  /***Plotting results***/
  if(wr && !tr && graph_flag && (strcmp(method,"complex")!=0)){
    gplot_results(data_name,complex_run_flag,method);
  }
  /*** Analysis of the trajectory ***/
  if(!complex_run_flag && !tr && wr){
    //if complex method or other method not included inside complex run
    if(ma_span){
      analyze_traj(tmp_name);
    }
    else
      analyze_traj(data_name);
  }
  //free(tmp_name);
  /*** Final printing to the user ***/
  cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
  if(strcmp(method,"complex")!=0){
    if((strcmp(method,"discrete")==0) || (lang_flag))
      fprintf(stdout,"Time: %G. Number of runs: %d. Time used: %G sec.\n",time,nruns,cpu_time_used);
    else{
      if(lang_flag)
	fprintf(stdout,"Time: %G. Number of runs: %d. Time used: %G sec.\n",time,nruns,cpu_time_used);
      else
	fprintf(stdout,"Time: %G. Number of steps: %d. Time used: %G sec.\n",time,n_steps,cpu_time_used);
    }
  }
  return 0;
}

int analyze_traj(const char *fname)
{/* The function is intended to compute the common entities out from the data file */
  int i,j;
  int complex_true,nRuns,pure_det,lf,*info;
  char *out,*tmpstr;
  FILE *tmp;
  /* if(cross_level[1] == 0){/\* The second value is zero *\/ */
  /*   out = (char *)malloc((strlen(fname)+strlen(out_aff)+1)*sizeof(char)); */
  /*   /\* Copy the <fname> to the <out> *\/ */
  /*   out = strcpy(out,fname); */
  /*   out = strcat(out,out_aff); */
  /* } */
  /* else{ */
  /*   tmpstr = (char *)malloc(20*sizeof(char)); */
  /*   sprintf(tmpstr,"%G",cross_level[0]); */
  /*   out = (char *)malloc((strlen(fname)+strlen(tmpstr)+1)*sizeof(char)); */
  /* } */
  out = (char *)malloc((strlen(fname)+4)*sizeof(char));
  for(j=0;j<3;j++){
    if(fabs(cross_level[j]) < 0.00001)/* Zero! */
      continue;
    else{
      /* tmpstr = (char *)malloc(20*sizeof(char)); */
      /* sprintf(tmpstr,".c%G",cross_level[j]);/\*e.g, '.c4' extension for cross_level = 4 *\/ */
      /* if(tmpstr==NULL){ */
      /* 	fprintf(stderr,"Error(sprintf/analyze_traj): could not allocate\n"); */
      /* 	return 10; */
      /* } */
      //out = (char *)malloc((strlen(fname)+strlen(tmpstr)+1)*sizeof(char));
      out = strcpy(out,fname);
      //out = strcat(out,tmpstr);
      out = strcat(out,".do");/* the new format for the (d)inamica (o)utput */

      tmp = fopen(fname,"r");
      if(tmp == NULL){
	fprintf(stderr,"Error opening file `%s'\n",fname);
      }
      if((info=get_info_data(info,tmp)) == NULL){
	fprintf(stderr,"Error occurred. Exit.\n");
	return 1;
      }
      fclose(tmp);
      printf("complex=%d,nRuns=%d,pure_det=%d,lf=%d\n",info[0],info[1],info[2],info[3]);
      complex_true = info[0];
      nRuns = info[1];
      pure_det = info[2];
      lf = info[3];
      if(complex_true){/* Complex run */
	perDet = compute_period(perDet,&nPerDet,cross_level[j],
				per_method,fname,1,0);
	perStoch = compute_period(perStoch,&nPerStoch,cross_level[j],
				  per_method,fname,nRuns,1);
      }
      else{/* Not complex run */
	if(!pure_det){
	  /* Stochastic method: langevin or discrete*/
	  nPerDet = 0;
	  perStoch = compute_period(perStoch,&nPerStoch,cross_level[j],
				    per_method,fname,nRuns,0);
	}
	else{
	  /* Only deterministic method */
	  nPerStoch = 0;
	  perDet = compute_period(perDet,&nPerDet,cross_level[j],
				  per_method,fname,nRuns,0);
	}
      }
      /*** Report about the analysis ***/
      if(j==0)
	tmp = fopen(out,"w");
      else
	tmp = fopen(out,"a");
      /* Print trajectory info array */
      set_info_data(info,tmp);

      /* /\* Print parameter values *\/ */
      /* fprintf(tmp,"#PARAMETERS:\n"); */
      /* for(i=0;i<PARDIM;i++){ */
      /* 	fprintf(tmp,"#%s=%G\n",par_name[i],mu.P[i]); */
      /* } */
      /* /\* Print variable names *\/ */
      /* fprintf(tmp,"#VARIABLES:\n"); */
      /* for(i=0;i<DIM;i++){ */
      /* 	fprintf(tmp,"#U(%d)=%s\n",i+1,var_name[i]); */
      /* } */
      /* fprintf(tmp,"#**********\n"); */
      /* Print period info */
      if(nPerDet){
	fprintf(tmp,"#PerDet(%s,c%G): mean = %G, std = %G, n = %d\n",
		var_name[perVarInd],cross_level[j],
		gsl_stats_mean(perDet,1,nPerDet),
		gsl_stats_sd(perDet,1,nPerDet),nPerDet);
	for(i=0;i<nPerDet;i++)
	  fprintf(tmp,"%G\n",perDet[i]);
	fprintf(tmp,"\n\n\n");
      }
      if(nPerStoch){
	if(!per_method)/* Poincare sections method */
	  fprintf(tmp,"#PerStoch(%s,c%G=%G): mean = %G, std = %G, n = %d\n",
		  var_name[perVarInd],cross_level[j],cross,
		  gsl_stats_mean(perStoch,1,nPerStoch),
		  gsl_stats_sd(perStoch,1,nPerStoch),nPerStoch);
	else/* autocorrelation method */
	  fprintf(tmp,"#PerStoch(%s,ac): mean = %G, std = %G, n = %d\n",
		  var_name[perVarInd],
		  gsl_stats_mean(perStoch,1,nPerStoch),
		  gsl_stats_sd(perStoch,1,nPerStoch),nPerStoch);
	for(i=0;i<nPerStoch;i++)
	  fprintf(tmp,"%G\n",perStoch[i]);
	fprintf(tmp,"\n");
      }
      fprintf(tmp,"\n\n");

      fclose(tmp);
      fprintf(stdout,"Output is %s\n",out);
    }
  }
  //free(info); /* not needed anymore */
  free(out);
  return 0;
}

int moving_average(const char *fname, const int span)
{/* This function computes the moving average (MA) of the simulated trajectory stored
    in fname. The span of the MA is span and must be some positive odd integer. The
    resulting MA-ed trajectory is stored in the output file, which name is
    `fname`.ma. The first and the last values of the trajectory that cannot be
    averaged with the given span are takes as they are. */
  /* Init the iteration variables */
  int i,j,k,m,n,*info;
  /* Init the temporary frame variable and tmp to store the temporary summation */
  double **frame,tmp;
  /* Open the incoming file for reading */
  FILE *in = fopen(fname,"r");
  /* Get the info from the file */
  if((info = get_info_data(info,in)) == NULL){
    fprintf(stderr,"Error occurred (fname = %s). Exit.\n",fname);
    return 1;
  }
  /* Get the output file name */
  char *fname_out = (char *)malloc((strlen(fname) + 4)*sizeof(char));
  fname_out = strcpy(fname_out,fname);
  fname_out = strcat(fname_out,".ma");
  /* Open the output file for writing */
  FILE *out = fopen(fname_out,"w");
  /* Set info data for the output file */
  set_info_data(info,out);
  /* Determine the number of frames in the data */
  if(info[0] == 1)
    info[1] = info[1] + 1;
  /* Reading frame by frame the data, compute the MA and write the output */
  printf("Calculating MA...");
  for(i=0;i<info[1];i++){/* iterate over the total number of frames */
    frame = load_frame(in,frame,&j);/* Load the frame */
    for(m=0;m<DIM;m++){/* over all dimensions, time array left unchanged */
      for(k=0;k<j;k++){/* over all data points in a frame */
	tmp = 0.0;/* init with zero temporary holder of the sum(see below) */
	if(k < (span-1)/2){/* not enough points at the beginning */
	  for(n=0;n<2*k+1;n++)/* form the sum */
	    tmp += frame[m+1][n];
	  frame[m+1][k] = tmp / (2*k+1);
	}
	else if(k >= j-(span-1)/2){/* not enough points at the end */
	  for(n=k-(j-1-k);n<j;n++)
	    tmp += frame[m+1][n];
	  frame[m+1][k] = tmp / (2*(j-1-k)+1);
	}
	else{/* average in the middle of the series, unmodified span */
	  for(n=0;n<span;n++)/* form the sum */
	    tmp += frame[m+1][k-(span-1)/2+n];
	  frame[m+1][k] = tmp / span;/* calculate the average */
	}
      }
    }
    write_frame(out,frame,j);
  }
  printf("Done.\n");
  //free(info);
  /* Close the output and incoming streams */
  fclose(in);
  fclose(out);
  
  return 0;
}

int buffer_overfull_check(){
  int i;
  if(write_count == BUFFER){
    fprintf(stdout,"Buffer(%d) is overfull(t=%lf), continue?[Y/n]\n",n_steps,t);
    i=getchar();
    while(i!='n' && i!='N' && i!='\n' && i!='y' && i!='Y'){
      fprintf(stdout,"Please `Y' or 'N':\n");
      i=getchar();
    }
    if(i=='n' || i=='N') return 1;
    if(i=='\n' || i=='y' || i=='Y'){
      BUFFER += BUF_INCR;
      ts = (double *)realloc(ts,BUFFER*sizeof(double));
      if(ts == NULL)
	printf("Error: couldn't allocate memory for time storage\n");
      xs = (double **)realloc(xs,BUFFER*sizeof(double *));
      if(xs == NULL)
	printf("Error: couldn't allocate memory for variable storage\n");
      for(i=0;i<BUFFER;i++){
	xs[i] = (double *)realloc(xs[i],DIM*sizeof(double));
	if(xs[i] == NULL)
	  printf("Error: couldn't allocate memory for variable points storage\n");
      }
      fprintf(stdout,"Allocated successfully %d points.\n",BUFFER);
    }
  }
  return 0;
}

int write_data_file(FILE *data){
  int i;
  fprintf(data,"%G ",t);
  for(i=0;i<DIM;i++){
    if(gsl_isnan(x[i])){
      printf("Error: NaN !\n");
      return 1;
    }
    fprintf(data,"%.6lf ",x[i]);
  }
  for(i=0;i<AUXNUM;i++)
    fprintf(data,"%.6lf ",a[i]);
  fprintf(data,"\n");
  return 0;
}

void write_data_array(){
  int i;
  ts[write_count] = t;
  for(i=0;i<DIM;i++)
    xs[write_count][i] = x[i];
}

char *app_ext(char const * basename, char const * ext, char * out)
{/* This function appends extension specified to the basename of file
    names*/
  int i,j;
  //char *out;
  out=(char *)calloc((strlen(basename)+strlen(ext)+2),sizeof(char));
  for(i=0;i<strlen(basename);i++)
    out[i]=basename[i];
  out[i]='.';i++;
  for(j=0;j<strlen(ext);j++)
    out[i+j]=ext[j];
  out[i+j]='\0';

  return out;

}
