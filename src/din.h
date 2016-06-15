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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
#include "config.h"
#include <stdio.h>

int short Just_Compile;/* flag not to fire up dinamica, just produce .c output */
int short No_Transfer;/* flag not to transfer to .c output */
char *odefname;/*String for .ode file name*/
char *basename;/*Basename for all file names*/
char *cfname;/* C source filename  */
char **odestr;/*ODE specification string*/
char **initstr; /*Array storing initial condition strings: `init x=0.3424,y=3.455' etc*/
char ***jacstr;/*Jacobian specification string*/
char **tjacstr;/*Time derivative of RHS: when not ODE*/
char **propstr;/*Strings containing the propensity declarations(Gill)*/
//char *upstr;/*The current update vector string(Gill)*/
int **upvec;/* The update vector of integers(Gill) */
char **lang_amend_str;/* Strings containing the langevin stoch amendments formula */

struct func {
  char *func_str;
  char *func_name;
  char **arg_list;
  int num_arg;
};
struct func *function;
struct aux {
  char * name;
  char * formula;
};
struct aux *aux;
char *gcc_;/*Compiler command*/
char *run_;/*Execute command*/
char * extr_base(const char *, char *);
char * modify_command(char *,const char *,const int);
char * app_ext(const char *, const char *, char *);
int read_source(const char *);
int process_odefile(const char *);
int process_odestr(int, FILE *,const char *, char **, char **,int);
int isPar(const char *,const int);
int isVar(const char *, char **, int);
int isFunc(char const *);
int isSym(const char *, const char);
int print_heading(FILE *);
int print_tail(FILE *);
int find_eqsign(const char *);
int gen_conf(const char *);
int init();
int process_par(char *);
int get_parameter_name(char *);
int process_inter_par(const char *);
void process_update_vector(char *, int *);
int process_init(char *);
int how_many_eqsigns(char const *);
int buf_strip(char *, const int);
