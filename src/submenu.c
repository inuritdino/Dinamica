/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012,2013 Elias Potapov. */
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
   Lebedev Physical Inst., Dep. of Theoretical Physics.
   Moscow, Russia */
/****************************************************************************************/
/* submenu.c is the file for interpreting commands(such as
   input_interpreter.c), but this interpretation is going on submenu
   commands. It also uses read_menu(...) function from
   input_interpreter.c file. Iterpretation functions are called in way
   of MENU_SUBMENU_interp(...), where MENU is appropriate menu of the
   program and SUBMENU is appropriate submenu of the program.
*/

#include "init.h"
#include "thistogram.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int periodics_hist_interp(char **buffer,const int depth_level)
{
  int static Bin = 10;
  int i,j,k;
  int Len = 20;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  double **histPer;
  if(strcasecmp(buffer[depth_level],"exit") == 0 ||			\
     strcasecmp(buffer[depth_level],"quit") == 0 ||			\
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving HISTOGRAM submenu...\n");
    return 1000;
  }
  else if(strcasecmp(buffer[depth_level],"ls")==0){
    fprintf(stdout,"Periodics/Histogram menu:\n");
    fprintf(stdout,"(Sh)ow\n");
    fprintf(stdout,"(R)un\n");
    fprintf(stdout,"(B)ins\n");
    //fprintf(stdout,"(B)ounds\n");
    //fprintf(stdout,"(D)eltaT\n");
    fprintf(stdout,"(O)utput file name\n");
    //fprintf(stdout,"(W)rite output flag\n");
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"sh")==0 ||	\
	  strcasecmp(buffer[depth_level],"show")==0){
    fprintf(stdout,"Bins:\n");
    fprintf(stdout,"\t%d\n",Bin);
    fprintf(stdout,"Output file name:\n");
    fprintf(stdout,"\t%s\n",hist);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"b")==0 ||		\
	  strcasecmp(buffer[depth_level],"bin")==0 ||		\
	  strcasecmp(buffer[depth_level],"bins")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter number of bins:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0)
      fprintf(stdout,"No change\n");
    else{
      Bin = atoi(symbuf);
      if(Bin < 0){
	fprintf(stdout,"Cannot be negative, inverting\n");
	Bin = -Bin;
      }
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"o")==0 ||		\
	   strcasecmp(buffer[depth_level],"out")==0 ||		\
	   strcasecmp(buffer[depth_level],"output")==0){
    fprintf(stdout,"Enter file name for output:\n");
    symbuf = read_input_line(symbuf,Len);
    if(*symbuf == 0){
      fprintf(stdout,"No changes\n");
    }
    else{
      hist = (char *)realloc(hist,(strlen(symbuf)+1)*sizeof(char));
      strncpy(hist,symbuf,strlen(symbuf));
    }
    fprintf(stdout,"Filename is %s\n",hist);
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"r")==0 ||	\
	   strcasecmp(buffer[depth_level],"run")==0){
    if(!nPerStoch){
      fprintf(stdout,"Cannot compute period distribution.\n");
      return 10;
    }
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      if(Bin != 0)
	histPer = compute_hist_per(histPer,perStoch,nPerStoch,Bin);
      else{
	Bin = 10;
	fprintf(stdout,"Running with %d bins\n",Bin);
	histPer = compute_hist_per(histPer,perStoch,nPerStoch,Bin);
      }
    }
    else{
      Bin = atoi(symbuf);
      if(Bin < 0){
	fprintf(stdout,"Cannot be negative, inverting\n");
	Bin = -Bin;
      }
      histPer = compute_hist_per(histPer,perStoch,nPerStoch,Bin);
    }
    write_histPer_file(histPer,Bin,hist,"w");
    if(graph_flag){/* print histogram */
      graph_set_labels(graph,"thist");
      gnuplot_cmd(plot_handle,"plot '%s' u 1:2 w boxes ls 1\n",hist);
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"")==0)
    return 0;
  else {
    fprintf(stdout,"No such command\n");
    return 0;
  }
}
