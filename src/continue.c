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
#include "continue.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_statistics_double.h>
#define TRUE 1
#define FALSE 0



int cont_init()
{
  par1 = 0;
  par2 = 1;
  pvalues = 0;
  npvalues = 0;
  P1in = mu.P[par1];
  P2in = mu.P[par2];
  P1stepin = 0.01;
  P2stepin = 0.01;
  P1max = 2*P1in;
  P1min = 0;
  P2max = 2*P2in;
  P2min = 0;
  
  write_num = 1;
  num_per_to_get = 0;
  dyn_check_flag_cont = 0;

  return 0;
}

int run_extend_rnd(const int pind,const double pvalues[],const int npvalues)
{
  int i,tmp;
  double ttime = mynum.total_time;/* Default total time */
  regStat *stat_all[npvalues];
  for(i=0;i<npvalues;i++){
    /* Change the par value */
    mu.P[pind] = pvalues[i];
    printf("|***** %s = %G *****|\n",par_name[pind],mu.P[pind]);
    /* First, let's run transiently to get an attractor */
    if((run(method,mynum.trans_time,0,1,1)) != 0)
      return 100;
    /* Try to get an estimate of the period */
    if(num_per_to_get && !per_method){
      /* if the method appropriate and user wants it*/
      tmp = graph_flag;/* Remember the graph_flag value */
      if(graph_flag)
	graph_flag = 0;/* Turn off graphic output, no need in this */
      if((run(method,mynum.total_time,1,0,1)) != 0)/* Run */
	return 100;
      graph_flag = tmp;/* Turn on graphic output again*/
      if(nPerDet || nPerStoch){
	if(nPerDet){
	  printf("Mean det per = %G\n",gsl_stats_mean(perDet,1,nPerDet));
	  ttime = num_per_to_get * gsl_stats_mean(perDet,1,nPerDet);
	}
	if(nPerStoch){
	  printf("Max stoch per = %G\n",gsl_stats_max(perStoch,1,nPerStoch));
	  //	  if(ttime < num_per_to_get * gsl_stats_max(perStoch,1,nPerStoch))
	  ttime = num_per_to_get * gsl_stats_max(perStoch,1,nPerStoch);
	}
      }
      else{/* If none period is found, then left default total time */
	ttime = mynum.total_time;
	if(!i)
	  fprintf(stdout,"Total time left unchanged %G.\n",ttime);
      }
    }
    /* Finally, run seriously */
    fprintf(stdout,"Run with time = %G\n",ttime);
    stat_all[i] = rnd_init_cond_run();
    //copy_regime_stat(&(stat_all[i]),stat);
    /* printf("STAT_ALL[%d].os_up_cn = %d\n",i,stat_all[i].os_up_cn); */
    /* printf("STAT_ALL[%d].os_pd[0] = %d\n",i,stat_all[i].os_pd[0]); */
    /* printf("STAT_ALL[%d].os_up[0] = %d\n",i,stat_all[i].os_up[0]); */
  }
  FILE *out = fopen("regstat.dat","w");
  if(out == NULL)
    fprintf(stderr,"Error: failed to open \"regstat.dat\"\n");
  statND_to_file(out,stat_all,npvalues);
  fclose(out);
  for(tmp=0;tmp<npvalues;tmp++)
    free_regime_stat(stat_all[tmp]);

  return 0;
}

int run_extend(const int pind,const double pvalues[],const int npvalues,
	       const int nRuns)
{
  int i,tmp;
  double ttime = mynum.total_time;/* Default total time */
  char *basename = (char *)malloc(FNAME_SIZE*sizeof(char));
  /* <basename> will hold the original data file name */
  strcpy(basename,data_name);
  /* <app> is the modified according to the par value */
  char *app = (char *)malloc(FNAME_SIZE*sizeof(char));
  for(i=0;i<npvalues;i++){
    /* Change the par value */
    mu.P[pind] = pvalues[i];
    printf("|***** %s = %G *****|\n",par_name[pind],mu.P[pind]);
    /* Copy the <basename> to the new data file name */
    strcpy(data_name,basename);
    /* Create the <app> from the par value */
    sprintf(app,"_%s%G",par_name[pind],mu.P[pind]);
    /* Append the <app> to the current file name <data_name> */
    if((strlen(data_name)+strlen(app)) < FNAME_SIZE)
      data_name = strncat(data_name,app,strlen(app));
    else{
      fprintf(stderr,"Error(run_extend): exceeded FNAME_SIZE: %s+%s\n",data_name,app);
      return 10;
    }
    /* Now, we have <data_name> we can run */
    /* First, let's run transiently to get an attractor */
    if((run(method,mynum.trans_time,0,1,1)) != 0)
      return 100;
    /* Try to get an estimate of the period */
    if(num_per_to_get && !per_method){
      /* if the method appropriate and user wants it*/
      tmp = graph_flag;/* Remember the graph_flag value */
      if(graph_flag)
	graph_flag = 0;/* Turn off graphic output, no need in this */
      if((run(method,mynum.total_time,1,0,1)) != 0)/* Run */
	return 100;
      graph_flag = tmp;/* Turn on graphic output again*/
      if(nPerDet || nPerStoch){
	if(nPerDet){
	  printf("Mean det per = %G\n",gsl_stats_mean(perDet,1,nPerDet));
	  ttime = num_per_to_get * gsl_stats_mean(perDet,1,nPerDet);
	}
	if(nPerStoch){
	  printf("Max stoch per = %G\n",gsl_stats_max(perStoch,1,nPerStoch));
	  //	  if(ttime < num_per_to_get * gsl_stats_max(perStoch,1,nPerStoch))
	  ttime = num_per_to_get * gsl_stats_max(perStoch,1,nPerStoch);
	}
      }
      else{/* If none period is found, then left default total time */
	ttime = mynum.total_time;
	if(!i)
	  fprintf(stdout,"Total time left unchanged %G.\n",ttime);
      }
    }
    /* Finally, run seriously */
    fprintf(stdout,"Run with time = %G\n",ttime);
    if((run(method,ttime,1,0,nRuns)) != 0)
      return 100;
  }
  /* Copy back the basename to the data file name */
  strcpy(data_name,basename);

  free(basename);
  free(app);
  return 0;
}




int continuation()
{
  int i,j,l;

  int write_count=0;
  int n=1;

  FILE *param;
  FILE *param_point;
  param=fopen("par","w");
  param_point=fopen("par_p","w");

  P1step=P1stepin;
  P2step=P2stepin;
  P1=P1in;
  P2=P2in;
  // BIG loop over 1-st parameter.
  while((P1 <= P1max) && (P1 >= P1min))
    {
      mu.P[par1]=P1;

      // BIG loop over 2-nd parameter.
      while((P2 <= P2max) && (P2 >= P2min))
	{
	  mu.P[par2]=P2;

	  fprintf(stdout,"***********************\n");
	  //Integrating...
	  run(method,mynum.total_time,write_flag,0,1);
	  //Getting attractor by additional integrating.
	  for(i=0;i<DIM;i++){
	    j=0;
	    while(per_ratio[i] > eps_per_ratio){
	      run(method,mynum.trans_time,0,1,1);
	      j++;
	      if(j>10){
		fprintf(stdout,"!Could not get attractor!\n");
		break;
	      }
	    }
	  }
	  //Determining trajectory.
	  // ???
	  //Reporting to stdout.
	  fprintf(stdout,"%s: %G\t%s: %G\n",par_name[par1],P1,par_name[par2],P2);
	  
	  //Now initial is current, both for pars and vars.
 	  if(cont_sol==regime){
	    for(j=0;j<DIM;j++){xin[j]=x[j];}
	    //P2in=P2;
	  }
	   	  
	  //Point is not successive...
	  if(cont_sol!=regime){
	    P2=P2-P2step;
	    for(j=0;j<DIM;j++){x[j]=xin[j];}
	    fprintf(param,"%Gf %Gf\n",P1,P2);
	    write_count++;
	    if(write_count == n*write_num){ fprintf(param_point,"%lf %lf\n",P1,P2);
	      for(l=0;l<DIM;l++)
		fprintf(param_point,"%lf\n",x[l]);
	      n++;
	    }
	    break;
	  }
	  //Boundary reflection.
	  if((P2 > P2max) || (P2 < P2min)) {P2=P2-P2step; P2step=-P2step;break;}

	  P2 = P2 + P2step;
	}// END of 2-nd par.
      P1 = P1 + P1step;
    }//END of 1-st par.
  fclose(param);
  fclose(param_point);

  return 0;

}
