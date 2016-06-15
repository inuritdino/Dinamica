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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
/****START OF THE PROGRAM****/
/* Main file of DINAMICA LIBRARY, it reads command line input
   arguments, sets environment, starts reading MAIN menu input*/

#include "init.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
//#include <plot.h>

int main(int argc, char** argv)
{
  int i;/*Looping variable*/
  init_command_line();

  if(argc == 1) {
    fprintf(stderr,"\nNo arguments passed to the program.\n");
  }
  for(i=1;i<argc;i++){
    if(argv[i][0] == '-') {
      switch(argv[i][1]) {
      case 'o': strncpy(out_name,argv[++i],FNAME_SIZE-1);
	break;
      case 'f': strncpy(input_name,argv[++i],FNAME_SIZE-1);//script file
	printf("The script file is: %s\n",input_name);
	break;
      case 'g': graph_flag = 0;//No graphics
	break;
      case 'd': strncpy(data_name,argv[++i],FNAME_SIZE-1);//data file
	break;
      case 'p': strncpy(init_name,argv[++i],FNAME_SIZE-1);
	break;
      case 'c': strncpy(conf_name,argv[++i],FNAME_SIZE-1);
	break;
      default: fprintf(stderr,"\nNo arguments passed to the program or"); 
	fprintf(stderr," I cannot recognize them.\n");
	break;
      }
    }
    else{
      /*Configuration file.*/
      strcpy(conf_name,argv[i]);
    }
  }
/* Copyright notice. All license stuff. */  
  copyleft();
/**********************************************************/
/* Initiate ALL stuff
 * ********************************************************/
  init();
/**********************************************************
 Reading input in main menu. The program rotates inside this. 
 **********************************************************/
  printf("\n");
  /* If input_name does not exist or empty then stream==NULL and reading is from the
    stdin. */
  input_stream = fopen(input_name,"r");
  while(1){
    //read_menu(cmd,main_prompt);
    i = parse_command_line(cmd,main_prompt,input_stream);
    if(i == 100)
      break;
    if(i == 200)
      continue;
    if((main_interp(cmd,0)) == 1000)
      break;
  }
  fclose(input_stream);
  din_close();
  return 0;
}

int max_min_x(double max_x[],double min_x[],double t_start,double t_stop)
{/*This function computes maximum and minimum values of variables in a
   given time interval (t_start,t_stop).*/
  int i,j;
  /*Set to initial values.*/
  int start=0;
  int stop=n_steps-1;
  double th;/*Help variable.*/
  if(t_start>t_stop){/*Interchange t_start<-->t_stop*/
    th=t_stop; t_stop=t_start; t_start=th;
  }
  /*Find start time in storage array.*/
  for(i=0;i<write_count;i++){
    if(t_start > ts[i]) continue;
    if(t_start <= ts[i]) {start=i;break;}/*'=' sign for the case when t_start=ts[0].*/
  }
  /*Find stop time in storage array.*/
  for(i=start;i<write_count;i++){
    if(t_stop >= ts[i]) continue;
    if(t_stop < ts[i]) {stop=i;break;}
  }
  /*Set max and min to first values in the specified interval.*/
  for(i=0;i<DIM;i++){
    max_x[i]=xs[start][i]; min_x[i]=xs[start][i];}
  /*Find max and min values in specified interval.*/
  for(j=start; j<=stop; j++){
    for(i=0; i<DIM; i++) {
      if(max_x[i] <= xs[j][i]) max_x[i] = xs[j][i];
      if(min_x[i] >= xs[j][i]) min_x[i] = xs[j][i];
    }
  }
  return 0;
}

void copyleft()
{
  fprintf(stdout,"\n\nDINAMICA Ver. %s (<%s>)\nCopyright 2008, 2009, 2010, 2011, 2012, 2013, 2014 Elias Potapov\n\n",PACKAGE_VERSION,PACKAGE_URL);
  fprintf(stdout,"Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, \n2005, 2006, 2007, 2008, 2009, 2010, 2011 The GSL Team\n\n");
  fprintf(stdout,"This software uses the gnuplot_i library written by N.Devillard\n(see <http://ndevilla.free.fr/gnuplot/>).\n\n");
  fprintf(stdout,"This program comes with ABSOLUTELY NO WARRANTY;\nfor details type `warranty' or simply `w'\n");
  fprintf(stdout,"This is free software, and you are welcome to redistribute it\n");
  fprintf(stdout,"under certain conditions; see GNU General Public License for details.\n\n");
  fprintf(stdout,"Report bugs to <%s>\n\n\n",PACKAGE_BUGREPORT);
}
