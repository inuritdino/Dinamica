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
   Tampere University of Technology, Dep. of Signal Processing.
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
/*
  This file defines (mainly) the main function for `dinamica' part of
  DINAMICA. `dinamica' compiles the .c source file, which is generated
  from .ode source file, against DINAMICA library(libdin.*) to produce
  second part of DINAMICA -- executable binary of the program with
  user defined equations(ODE).
 */
#include "din.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define ARG_LEN 50

int
main(int argc, char **argv)
{
  int c,i,j;/*Working variables*/
  char *tmp;/* Temporary string, primarily used for app_ext(...) */
  int length = 10;/*length of a string*/

  /* These values initialized right here, no init function for them:( */
  Just_Compile = 0;
  No_Transfer = 0;
  /* Temporary string for compilation command forming */
  int Len = 200;
  char *tmp_cc = (char *)malloc(Len*sizeof(char));
  /* Arguments for the .din executable */
  int arg_to_din_count = 0;
  char **arg_to_din = (char **)malloc(10*sizeof(char *));
  for(i=0;i<10;i++)
    arg_to_din[i] = (char *)malloc(ARG_LEN*sizeof(char));
  if(argc==1){/*No arguments--> we need .ode file at least.*/
    fprintf(stdout,"You need to specify .ode file to read at least:\n");
    fprintf(stdout,"Enter file name(no auto-filling):\n");
    odefname = (char *)calloc(length,sizeof(char));
    if(odefname==NULL){
      fprintf(stderr,"Error: I cannot allocate memory for your input.\n");
      	return -10;/*Bad allocating*/
    }
    /*Copying input: online file specification*/
    i = 0;
    while((c=getchar())!='\n'){
      i++;
      odefname[i-1] = (char)c;
      if(i==length-1){
	length += 10;
	odefname = (char *)realloc(odefname,length*sizeof(char));
      }
      if(odefname==NULL){
	fprintf(stderr,"Error: I cannot allocate memory for your input.\n");
	return -10;/*Bad allocating*/
      }
    }
    odefname[i] = '\0';/*Terminatin null character*/
    if(strlen(odefname)==0){
      fprintf(stderr,"No file, exit...\n");
      fprintf(stderr,"If you want analysis mode, put '-a' option.\n");
      return 0;
    }
    basename = extr_base(odefname,basename);
    /****************************/
    /**Process ODE file**/
    if((read_source(odefname))!=0){
      fprintf(stderr,"Parsing error occurred. Interrupting procedure.\n");
	  return 23;
    }
    /****************************/
  }
  else {/*At least one command line argument specified*/
    i = 1;
    while(i < argc){
      if(argv[i][0]=='-'){/*This is for `-x' arguments, where
			    x -- letter or digit*/
	switch(argv[i][1]){
	case 'c': Just_Compile = 1;//Just compilation, no start off.
	  break;
	case 'g': strncpy(arg_to_din[arg_to_din_count++],argv[i],ARG_LEN-1);//No graphics output
	  break;
	case 'd': strncpy(arg_to_din[arg_to_din_count++],argv[i],ARG_LEN-1);
	  strncpy(arg_to_din[arg_to_din_count++],argv[++i],ARG_LEN-1);
	  break;
	case 'f': strncpy(arg_to_din[arg_to_din_count++],argv[i],ARG_LEN-1);
	  strncpy(arg_to_din[arg_to_din_count++],argv[++i],ARG_LEN-1);
	  break;
	case 'p': strncpy(arg_to_din[arg_to_din_count++],argv[i],ARG_LEN-1);
	  strncpy(arg_to_din[arg_to_din_count++],argv[++i],ARG_LEN-1);
	  break;
	default: fprintf(stdout,"Error: I don't know the argument `%s'.\n",argv[i]);
	  break;
	}
      }
      else{/*If we have other arguments, not only in unix-like style `-x'
	     notation.*/
	basename = extr_base(argv[i],basename);
	/***********************/
	/*Read .ode file, which name is first argument*/
	if((read_source(argv[i])) != 0){
	  fprintf(stderr,"Parsing error occurred. Interrupting procedure.\n");
	  return 23;
	}
	/***********************/
	/*argv[i] is copied to fname*/
	odefname = (char *)calloc((strlen(argv[i])+1),sizeof(char));
	if(odefname==NULL){
	  fprintf(stderr,"Error: memory allocating: argv[%d]-->fname\n",i);
	  return -10;
	}
	strcpy(odefname,argv[i]);
      }
      i++;
    }
  }
  /****************************************************************/
  /* The system is defined in .ode file and transfered to .c file */
  if(!No_Transfer)
    cfname = app_ext(basename,"c",cfname);
  /*Otherwise the system is defined in external .c file indicated in the
    #include statement of the .ode*/
  /************************************************************************/
  //printf("DIN_LDFLAGS=%s;len=%lu\n",DIN_LDFLAGS,strlen(DIN_LDFLAGS));
  fprintf(stdout,"Preparing for `%s' compiling...\n",cfname);
  gcc_ = modify_command(gcc_,"gcc -Wall",0);
  gcc_ = modify_command(gcc_,cfname,0);
  while((strlen(DIN_CPPFLAGS)+strlen(LIBDIN_PATH)+strlen(DIN_LDFLAGS)+strlen("-L -o -L./")) >= Len-1){
    Len = Len*2;
    tmp_cc = (char *)realloc(tmp_cc,Len*sizeof(char));
  }
  if(strcmp(LIBDIN_PATH,"NONE/lib")!=0){
    //printf("LIBDIN_PATH=%s;len=%lu\n",LIBDIN_PATH,strlen(LIBDIN_PATH));
    /* If the prefix was specified */
    tmp_cc = strcat(tmp_cc,"-L");
    tmp_cc = strcat(tmp_cc,LIBDIN_PATH);
    tmp_cc = strcat(tmp_cc," ");
  }
  tmp_cc = strcat(tmp_cc,DIN_CPPFLAGS);
  tmp_cc = strcat(tmp_cc," ");
  tmp_cc = strcat(tmp_cc,DIN_LDFLAGS);
  tmp_cc = strcat(tmp_cc," ");
  /* Search current directory as well */
  tmp_cc = strcat(tmp_cc,"-L");
  tmp_cc = strcat(tmp_cc,getenv("PWD"));/* current directory */
  tmp_cc = strcat(tmp_cc," -o");
  //printf("tmp_cc=`%s'\n",tmp_cc);
  //gcc_ = modify_command(gcc_,"-L./src -L/usr/local/lib -L/opt/local/lib -o",0);
  gcc_ = modify_command(gcc_,tmp_cc,0);
  gcc_ = modify_command(gcc_,app_ext(basename,"din",tmp),0);
  gcc_ = modify_command(gcc_,"-ldin -lgsl -lgslcblas -lm",0);
  fprintf(stdout,"Compiling with: `%s'\n",gcc_);
  system(gcc_);
  if(!Just_Compile){
    fprintf(stdout,"Starting %s...\n",app_ext(basename,"din",tmp));
    run_ = modify_command(run_,"./",1);
    run_ = modify_command(run_,app_ext(basename,"din",tmp),0);
    run_ = modify_command(run_,app_ext(basename,"bcf",tmp),0);
    for(i=0;i<arg_to_din_count;i++)/* Copy arguments to .din  */
      run_ = modify_command(run_,arg_to_din[i],0);
    fprintf(stdout,"%s\n",run_);
    system(run_);
  }
  else {
    fprintf(stdout,"Compilation completed. Use output: `./%s.din %s.bcf'\n",basename,basename);
  }

  if(!Just_Compile){
    if((fopen(app_ext(basename,"bcf",tmp),"r")) != NULL){
      if(remove(app_ext(basename,"bcf",tmp)) != 0)
	printf("Error deleting file %s\n",app_ext(basename,"bcf",tmp));
    }
    if((fopen(app_ext(basename,"din",tmp),"r")) != NULL){
      if(remove(app_ext(basename,"din",tmp)) != 0)
	printf("Error deleting file %s\n",app_ext(basename,"din",tmp));
    }
    if((fopen(app_ext(basename,"c",tmp),"r")) != NULL){
      if(remove(app_ext(basename,"c",tmp)) != 0)
	printf("Error deleting file %s\n",app_ext(basename,"c",tmp));
    }
  }
  free(run_);
  free(tmp_cc);
  return 0;
}

char *extr_base(const char * filename, char * basename)
{
  int i,j;
  j = strlen(filename);/*Save length of the filename*/
  if(filename[j-4]=='.' && filename[j-3]=='o' && filename[j-2]=='d' &&
     filename[j-1]=='e'){
    basename=(char *)calloc((j-3),sizeof(char));
    if(basename == NULL) fprintf(stderr,"Error while reallocating string `%s'\n",basename);
    /*** Own code for copying ***/
    /* for(i=0;i<j-4;i++) */
    /*   basename[i]=filename[i]; */
    /* basename[i]='\0'; */
    basename = strncpy(basename,filename,j-4);
  }
  else if(filename[j-1]=='c' && filename[j-2]=='.'){
    basename = (char *)calloc((j-1),sizeof(char));
    if(basename==NULL) fprintf(stderr,"Error while reallocating string `%s'\n",basename);
    basename = strncpy(basename,filename,j-2);
  }
  else {
    basename = (char *)calloc(j,sizeof(char));
    if(basename==NULL) fprintf(stderr,"Error while reallocating string `%s'\n",basename);
    strcpy(basename,filename);
  }

  return basename;
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

char * modify_command(char * command, const char *append, const int tight)
{
  int i,j;
  /*This function appends string to existing command string to be
    executed. `append' string is usually new command options and
    arguments. `tight' means whether we do want spaces between
    `command' and `append': For example there is a different commands
    `gcc -g' and `./dinamica', first with tight=0, second with
    tight=1*/
  j = strlen(append);
  if(command != NULL){
    i = strlen(command);
    if(!tight)
      command = (char *)realloc(command,(i+j+2)*sizeof(char));
    else
      command = (char *)realloc(command,(i+j+1)*sizeof(char));
  }
  else{
    i = 0;
    if(!tight)
      command = (char *)calloc((j+2),sizeof(char));
    else
      command = (char *)calloc((j+1),sizeof(char));
  }
  for(j=0;command[j]!='\0';j++);
  for(i=0;i<strlen(append);i++)
    command[j+i]=append[i];
  if(!tight){
    command[j+i]=' ';
    command[j+i+1]='\0';}
  else
    command[j+i]='\0';
  
  return command;
}
