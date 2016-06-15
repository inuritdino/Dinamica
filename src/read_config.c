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
   Lebedev Physical Inst., Dep. of Theoretical Physics.
   Moscow, Russia */
/****************************************************************************************/
/* read_config.c file contains main algorithms to read configuration
   file of DINAMICA. There are two functions: read_config(...) and
   buf_check(...). First one reads the specified file. Second one
   deletes all blank positions and newlines and tabs from input.
*/

#include "init.h"
#include "errors.h"
#include "continue.h"
#include "random.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
// For reading MAX_N_STR strings from config file. May be increased.
#define MAX_N_STR 100
//#define STRLEN 100
/* Some errors */
#define BAD_FILE -1


int read_conf(char *fname)
{
  int i;
  FILE *file;
  fprintf(stdout,"Reading `%s'...",fname);
  file = fopen(fname,"r");
  if(file==NULL){
    fprintf(stderr,"Error: could not open file `%s'\n",fname);
    return BAD_FILE;
  }
  /**************************************************************/
  fread(&DIM,sizeof(int),(size_t)1,file);
  fread(&LDIM,sizeof(int),(size_t)1,file);
  fread(&PARDIM,sizeof(int),(size_t)1,file);
  fread(&FNUM,sizeof(int),(size_t)1,file);
  fread(&AUXNUM,sizeof(int),(size_t)1,file);
  fread(&NREAC,sizeof(int),(size_t)1,file);
  var_name = (char **)malloc(DIM*sizeof(char *));
  for(i=0;i<DIM;i++){
    var_name[i] = (char *)malloc(VAR_NAME_SIZE*sizeof(char));
    fread(var_name[i],sizeof(char),VAR_NAME_SIZE,file);
  }
  xin = (double *)malloc(DIM*sizeof(double));
  yin = (long int *)malloc(DIM*sizeof(long int));
  xin_par = (int *)malloc(DIM*sizeof(int));
  fread(xin_par,sizeof(int),(size_t)DIM,file);
  fread(xin,sizeof(double),(size_t)DIM,file);
  fread(yin,sizeof(long int),(size_t)DIM,file);
  for(i=0;i<PARDIM;i++)
    fread(par_name[i],sizeof(char),PAR_NAME_SIZE,file);
  fread(mu.P,sizeof(double),(size_t)PARDIM,file);
  int auxlen;/*Length of the aux name string, not permanent(!)*/
  aux_name = (char **)malloc(AUXNUM*sizeof(char *));
  for(i=0;i<AUXNUM;i++){
    fread(&auxlen,sizeof(int),(size_t)1,file);
    aux_name[i] = (char *)malloc(auxlen*sizeof(char));
    fread(aux_name[i],sizeof(char),auxlen,file);
  }
  fread(&mynum.total_time,sizeof(double),(size_t)1,file);
  fread(&mynum.trans_time,sizeof(double),(size_t)1,file);
  fread(&mynum.step,sizeof(double),(size_t)1,file);
  fread(&mynum.write_step,sizeof(int),(size_t)1,file);
  fread(&mynum.smp_frq,sizeof(double),(size_t)1,file);
  fread(method,sizeof(char),(size_t)METH_NAME_LEN,file);
  fread(method2,sizeof(char),(size_t)METH_NAME_LEN,file);
  fread(&traj_trans,sizeof(int short),(size_t)1,file);
  fread(&per_method,sizeof(int short),(size_t)1,file);
  fread(&per_var,sizeof(int short),(size_t)1,file);
  fread(&per_thresh,sizeof(double),(size_t)1,file);
  fread(&cross_level[0],sizeof(double),(size_t)1,file);
  fread(&ma_span,sizeof(int),(size_t)1,file);
  fread(&BUFFER,sizeof(int),(size_t)1,file);
  fread(&BUF_INCR,sizeof(int),(size_t)1,file);
  fread(&write_flag,sizeof(int short),(size_t)1,file);
  fread(&jac_flag,sizeof(int short),(size_t)1,file);
  fread(&lang_flag,sizeof(int short),(size_t)1,file);
  fread(&graph_flag,sizeof(int short),(size_t)1,file);
  fread(&eps_abs_tper,sizeof(double),(size_t)1,file);
  fread(&eps_rel_tper,sizeof(double),(size_t)1,file);
  fread(&eps_abs_int,sizeof(double),(size_t)1,file);
  fread(&eps_rel_int,sizeof(double),(size_t)1,file);
  fread(&a_y,sizeof(double),(size_t)1,file);
  fread(&a_dydt,sizeof(double),(size_t)1,file);
  fread(&eps_abs_per,sizeof(double),(size_t)1,file);
  fread(&eps_rel_per,sizeof(double),(size_t)1,file);
  fread(&eps_abs_peak,sizeof(double),(size_t)1,file);
  fread(&eps_rel_peak,sizeof(double),(size_t)1,file);
  fread(&eps_abs_am,sizeof(double),(size_t)1,file);
  fread(&eps_rel_am,sizeof(double),(size_t)1,file);
  fread(&eps_per_ratio,sizeof(double),(size_t)1,file);
  fread(&eps_inregime_ratio,sizeof(double),(size_t)1,file);
  fread(&graph,sizeof(struct coordNet),(size_t)1,file);
  
  fclose(file);
  fprintf(stdout,"Done.\n");
  return 0;
}


int read_config(char *filename)
{
  int i,j,l;
  int pos,count,count1;
  int eqsign;
  int length=100;
  FILE *file;
  file=fopen(filename,"r");
  if(file == NULL)
    return BAD_FILE;
  char *temp=(char *)calloc(length,sizeof(char));
  char *left=(char *)calloc(length,sizeof(char));
  char *right=(char *)calloc(length,sizeof(char));
  fprintf(stdout,"Starting to read %s configuration file...\n",filename);
  while((fgets(temp,length,file))!=NULL)/*Read line*/
    {
      //buf_check(temp);/*Remove nulls,tabs, spaces*/
      eqsign=0;/*Initially set eqsign to zero*/
      if(temp[0]=='#')/*We deal with comments*/
	continue;
      i=0;
      while(temp[i]!='=' && temp[i]!='\n'){
	left[i]=temp[i];
	i++;
      }
      left[i]='\0';
      if(temp[i]!='\n'){
	eqsign=i;
	i++;
      }
      while(temp[i]!='\n'){
	right[i-eqsign-1]=temp[i];
	i++;
      }
      right[i]='\0';
      if(eqsign==0){/*A keyword*/
	if(strcmp(left,"var")==0)
	  set_var(file);
	if(strcmp(left,"par")==0)
	  set_par(file);
	if(strcmp(left,"aux")==0)
	  set_aux(file);
      }
      else/*A equality*/
	set(left,right);
      memset(temp,'\0',length);
      memset(left,'\0',length);
      memset(right,'\0',length);
    }
  fprintf(stdout,"Done.\n");
  fclose(file);
  return 0;
}

int set(char * left, char *right)
{
  if(strcmp(left,"td")==0)
    DIM=atoi(right);
  if(strcmp(left,"ld")==0)
    LDIM=atoi(right);
  if(strcmp(left,"pd")==0)
    PARDIM=atoi(right);
  if(strcmp(left,"fd")==0)
    FNUM=atoi(right);
  if(strcmp(left,"ad")==0)
    AUXNUM=atoi(right);
  if(strcmp(left,"nr")==0)
    NREAC=atoi(right);
  if(strcmp(left,"tt")==0)
    mynum.total_time=atof(right);
  if(strcmp(left,"tr")==0)
    mynum.trans_time=atof(right);
  if(strcmp(left,"st")==0)
    mynum.step=atof(right);
  if(strcmp(left,"ws")==0)
    mynum.write_step=atoi(right);
  if(strcmp(left,"meth")==0){
    //method=(char *)calloc(strlen(right),sizeof(char));
    if(strlen(right) < METH_NAME_LEN)
      strcpy(method,right);
    else{
      fprintf(stderr,"Could not write the method.\n");
      fprintf(stderr,"Method name is too long.\n");
    }
  }
  if(strcmp(left,"tjt")==0)
    traj_trans=atoi(right);
  if(strcmp(left,"pm")==0)
    per_method=atoi(right);
  if(strcmp(left,"gb")==0)
    BUFFER=atoi(right);
  if(strcmp(left,"bi")==0)
    BUF_INCR=atoi(right);
  if(strcmp(left,"wf")==0)
    write_flag=atoi(right);
  /*Errors*/
  if(strcmp(left,"eat")==0)
    eps_abs_tper=atof(right);
  if(strcmp(left,"ert")==0)
    eps_rel_tper=atof(right);
  if(strcmp(left,"eai")==0)
    eps_abs_int=atof(right);
  if(strcmp(left,"eri")==0)
    eps_rel_int=atof(right);
  if(strcmp(left,"ay")==0)
    a_y=atof(right);
  if(strcmp(left,"ady")==0)
    a_dydt=atof(right);
  if(strcmp(left,"eap")==0)
    eps_abs_per=atof(right);
  if(strcmp(left,"erp")==0)
    eps_rel_per=atof(right);
  if(strcmp(left,"eapk")==0)
    eps_abs_peak=atof(right);
  if(strcmp(left,"erpk")==0)
    eps_rel_peak=atof(right);
  if(strcmp(left,"eaa")==0)
    eps_abs_am=atof(right);
  if(strcmp(left,"era")==0)
    eps_rel_am=atof(right);
  if(strcmp(left,"epr")==0)
    eps_per_ratio=atof(right);
  if(strcmp(left,"eir")==0)
    eps_inregime_ratio=atof(right);
  /*End fo errors*/
  return 0;
}

int set_var(FILE *file)
{
  int i,c;
  int static vindex=0;
  int length=50;
  char *left=calloc(length,sizeof(char));
  char *right=calloc(length,sizeof(char));
  i=0;
  while((c=getc(file))!='='){
    if(c=='\n')
      return 0;
    left[i]=(char)c;
    i++;
  }
  left[i]='\0';
  i=0;
  while((c=getc(file))!='\n'){
    right[i]=(char)c;
    i++;
  }
  right[i]='\0';
  var_name=(char **)realloc(var_name,(vindex+1)*sizeof(char *));
  var_name[vindex]=(char *)realloc(var_name[vindex],VAR_NAME_SIZE*sizeof(char));
  for(i=0;i<strlen(left);i++)
    var_name[vindex][i]=left[i];
  xin[vindex]=atof(right);/*Set initials*/
  x[vindex]=xin[vindex];/*Set default values of x.*/
  vindex++;
  set_var(file);

  return 1;
}

int set_par(FILE *file)
{
  int i,c;
  int static pindex=0;
  int length=50;
  char *left=calloc(length,sizeof(char));
  char *right=calloc(length,sizeof(char));
  i=0;
  while((c=getc(file))!='='){
    if(c=='\n')
      return 0;
    left[i]=(char)c;
    i++;
  }
  left[i]='\0';
  i=0;
  while((c=getc(file))!='\n'){
    right[i]=(char)c;
    i++;
  }
  right[i]='\0';
  for(i=0;i<strlen(left);i++)
    par_name[pindex][i]=left[i];
  mu.P[pindex]=atof(right);
  pindex++;
  set_par(file);

  return 1;
}

int set_aux(FILE *file)
{
  int i,c;
  int static aindex=0;
  int length=50;
  char *left=calloc(length,sizeof(char));

  i=0;
  while((c=getc(file))!='\n'){
    left[i]=(char)c;
    i++;
  }
  left[i]='\0';
  if((strcmp(left,"endaux"))==0)
    return 0;
  aux_name=(char **)realloc(aux_name,(aindex+1)*sizeof(char *));
  aux_name[aindex]=(char *)realloc(aux_name[aindex],VAR_NAME_SIZE*sizeof(char));
  for(i=0;i<strlen(left);i++)
    aux_name[aindex][i]=left[i];

  aindex++;
  set_aux(file);

  return 1;
}

int buf_check(char *buffer)
{/* This function removes newlines, spaces, and tabs and substitutes them with '\0' */
  /* characters with merging the gaps occuring after the substitution */
  int i,j;
  int len=strlen(buffer);
  //  printf("`%s'[%d] =>",buffer,len);
  for(i=0;i<len;i++){
    if(buffer[i]=='\n' || buffer[i]==' ' || buffer[i]=='\t'){
      buffer[i]='\0';
      for(j=i;j<len-1;j++)
	buffer[j]=buffer[j+1];
    }
  }
  //  printf("`%s'\n",buffer);
  return 0;
}
