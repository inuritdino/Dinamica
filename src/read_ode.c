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
   Tampere University of Technology, Dep. of Mathematics.
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
/*
  This file defines routines which are responsible for the system (ODE or stoch)
  specification processing along with other information in ode-files. The system
  specification is ASCII file describing the system under study. This system's
  specification is intended to mostly mimic .ode source files for XPPAUT program
  written by B. Ermentrout. read_source(...) routine reads .ode source file and
  converts its input to .c source file output. Additionally, the binary configuration
  file is created, which contains a lot of other information like technical
  parameters for the simulations. After this processing the .c source file is to be
  compiled with DINAMICA library (libdin.*) to produce executable program to
  run. Compilation is done with `dinamica' command (see din_main.c for details).
 */
#include "init.h"
#include "din.h"
#include "errors.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#define STRLEN 1000/* STRLEN is a limit to any string length in ode file*/
#define BAD_FILE -1
int read_source(const char * filename)
{
  char *basename,*p,*p1,*p2,*p3;
  basename = extr_base(filename,basename);
  int i,j,k,m,eqsign,neqsigns;/*i,j,k,m-looping variables, eqsign-equal
		 sign '=' position in specification string,
		 neqsigns-number of '=' signs.*/
  int colon,semicolon;/*colon and semicolon position of the propensity
			specification in the stoch discrete approach*/
  int open,close;/* open and close parenthesis positions */
  int initcount = 0;/*Counter for how many times we had `init x=1.2,y=0.2' etc string*/
  int jvarc;/*How many derivatives in Jacobian for single variable*/
  int jvar;/*Variable index of which Jacobian strings are formed*/
  int jaccount = 0;/* How many jac statements read, counter */
  int jdim = 0;/* dimension detected in the Jacobian statements, is checked */
	   /*   after parsing the whole file, if not equal to DIM, then return */
  FILE *file;/* Our source associated stream */
  file = fopen(filename,"r");/*Open stream*/
  if(file == NULL){/* In case of error */
    if(filename == NULL)
      fprintf(stdout,"Error: no file name.\n");
    else
      fprintf(stdout,"Error while reading source file `%s'",filename);
    return 10;
  }/*Bad file*/
  char *temp;/*Create temporary string for reading `filename' string by
	       string*/
  char *temp1;/*For continued lines and other*/
  char *unstrp;/*Unstripped string, when spaces are crucial*/
  temp = (char *)calloc(STRLEN,sizeof(char));/*Allocate memory*/
  temp1 = (char *)calloc(STRLEN,sizeof(char));/*Allocate memory*/
  unstrp = (char *)calloc(STRLEN,sizeof(char));/*Allocate memory*/
  /*We assume that STRLEN is enough for all type of equations*/
  /***************************************************************/
  /*                   PRINT some info to User                   */
  /***************************************************************/
  fprintf(stdout,"Starting to check system's specification:\n");
  /**************************************************************/
  /* Generate INITIAL values */
  /**************************************************************/
  init();
  /**************************************************************/
  /* Start to READ input */
  /**************************************************************/
  while((fgets(temp,STRLEN,file))!=NULL){/*Get string from source `filename'*/
    /* Remember temp holds a newline '\n' */
    strcpy(unstrp,temp);/*Making duplicate: one is stripped, another
			  is not(for some options this is crucial).*/
    if((isSym(temp,(char)'\\'))>=0){  /*Continuation line*/
      fgets(temp1,STRLEN,file);
      if(temp1==NULL){/*Error*/
	fprintf(stderr,"Error: reading continued line\n");
	return 156;}
      temp[(isSym(temp,(char)'\\'))]='\0';/*Replace '\\' continuation
					    character with null character*/
      strcat(temp,temp1);
      /*If there are more continuation lines*/
      while((isSym(temp1,(char)'\\'))>=0){
	fgets(temp1,STRLEN,file);
	if(temp1==NULL){/*Error*/
	  fprintf(stderr,"Error: reading continued line\n");
	  return 156;}
	temp[(isSym(temp,(char)'\\'))]='\0';
	strcat(temp,temp1);
      }
      memset(temp1,'\0',STRLEN);
    }
    buf_strip(temp,STRLEN);/*Remove Space,Tabs and shrink all \0
		      characters*/
    if((strlen(temp)==1) || (strcmp(temp,"done\n")==0) ||	\
       (strcmp(temp,"d\n")==0))	{
      /*Empty line(with only '\n' in it) OR `done' keyword. Ignore them.*/
      memset(temp,'\0',STRLEN);
      continue;
    }
    /************************************************************************/
    /*Everything that starts from hash sign '#' or '%': comments, directives*/
    /************************************************************************/
    if(temp[0] == '#' || temp[0] == '%'){/*Comment or the directives*/
      if(!strncmp(temp,"#include",strlen("#include")) ||	\
	 !strncmp(temp,"%include",strlen("%include"))){/*Include .c file*/
	if(temp[0] == '#'){/* Print warning message about deprecation */
	  printf("Warning: #include is deprecated. Use %%include instead.\n");
	}
	/* The system is directly from .c source */
        char *p = strchr(temp,'\n');
	*p = 0;/* Remove newline, making the end of the line */
	No_Transfer = 1;/* Set the flag positive */
	cfname = malloc((strlen(temp))*sizeof(char));/* Allocate more */
	strcpy(cfname,temp+8);
	printf("C source specified: `%s'\n",cfname);
      }
      if(!strncmp(temp,"#system",strlen("#system")) ||			\
	 !strncmp(temp,"%system",strlen("%system"))){/*# symm. systems(LDIM)*/
	if(temp[0] == '#'){/* Print warning message about deprecation */
	  printf("Warning: #system is deprecated. Use %%system instead.\n");
	}
	char *p = strchr(temp,'\n');
	*p = 0;/* Remove the newline and make the end of the line */
	LDIM = atoi((temp+7));/* +7 exactly after `#system' */
	if(LDIM < 0)/* Negative LDIM is not allowed */
	  LDIM = -LDIM;
	//printf("System arg: %d\n",LDIM);
      }
      //memset(temp,'\0',STRLEN);//seems useless.
      continue;
    }
    /*Find an equal sign(=) position which is required for ODE,
      functions, auxillary entities etc. Also, how many equal signs
      are in the string.*/
    eqsign = find_eqsign(temp);/*Find equal sign '=' position, in C numbering*/
    neqsigns = how_many_eqsigns(temp);
    /*************************************************************/
    /*CHECK for the ODE specification in the line*/
    if((temp[0]=='d' && temp[eqsign-3]=='/' && temp[eqsign-2]=='d' && 
	temp[eqsign-1]=='t') ||						\
       (temp[0]=='o' && temp[1]=='d' && temp[2]=='e') ||		\
       (temp[eqsign-1]=='\'') ||					\
       (temp[eqsign-3]=='(' && temp[eqsign-2]=='t' && temp[eqsign-1]==')')
       ){
      if(DIM==0)
      	var_name = (char **)malloc(sizeof(char *));
      else
	var_name = (char **)realloc(var_name,(DIM+1)*sizeof(char *));
      if(var_name==NULL){
	fprintf(stderr,"Couldn't allocate memory for var_name\n");
	return 100;
      }
      var_name[DIM] = (char *)malloc(VAR_NAME_SIZE*sizeof(char));
      if(var_name[DIM]==NULL){
	fprintf(stderr,"Error: memory for var_name[%d]\n",DIM);
	return 100;
      }
      /* Getting the variable names */
      if(temp[0] == 'd' && temp[eqsign-3]=='/' && temp[eqsign-2]=='d'){
	if((eqsign-4) > (VAR_NAME_SIZE-1)){
	  fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n",
		  (int)VAR_NAME_SIZE,temp[1],temp[2],temp[3]);
	  return 150;
	}
	for(i=1;i<eqsign-3;i++)
	  var_name[DIM][i-1] = temp[i];
	var_name[DIM][i-1] = '\0';
      }
      else if(temp[0] == 'o' && temp[1]=='d' && temp[2]=='e'){
	if((eqsign-3) > (VAR_NAME_SIZE-1)){
	  fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n",
		  (int)VAR_NAME_SIZE,temp[3],temp[4],temp[5]);
	  return 150;
	}
	for(i=3;i<eqsign;i++)/*We get variable names*/
	  var_name[DIM][i-3] = temp[i];
	var_name[DIM][i-3] = '\0';/*Add terminating null character*/
      }
      else if(temp[eqsign-1] == '\''){
	if((eqsign-1) > (VAR_NAME_SIZE-1)){
	  fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n",
		  (int)VAR_NAME_SIZE,temp[0],temp[1],temp[2]);
	  return 150;
	}
	for(i=0;i<eqsign-1;i++)/*We get variable names*/
	  var_name[DIM][i] = temp[i];
	var_name[DIM][i] = '\0';/*Add terminating null character*/
      }
      else if(temp[eqsign-1] == ')'){
	if((eqsign-3) > (VAR_NAME_SIZE-1)){
	  fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n",
		  (int)VAR_NAME_SIZE,temp[0],temp[1],temp[2]);
	  return 150;
	}
	for(i=0;i<eqsign-3;i++)/*We get variable names*/
	  var_name[DIM][i] = temp[i];
	var_name[DIM][i] = '\0';/*Add terminating null character*/
      }
      /* ***** */
      DIM++; /*Let's increase dimension of the system by one*/
      /* Allocating memory for ODE strings */
      if(DIM==1)
      	odestr = (char **)malloc(sizeof(char *));
      else
	odestr = (char **)realloc(odestr,DIM*sizeof(char *));
      if(odestr==NULL){
	fprintf(stderr,"Error: memory odestr\n");
	return 101;
      }
      odestr[DIM-1] = (char *)malloc((strlen(temp)-eqsign)*sizeof(char));
      if(odestr[DIM-1]==NULL){
	fprintf(stderr,"Error: memory odestr[%d]\n",DIM-1);
	return 101;
      }
      for(i=eqsign+1;i<strlen(temp);i++)/*Now we get ODE string*/
	odestr[DIM-1][i-eqsign-1] = temp[i];
      odestr[DIM-1][i-eqsign-1] = '\0';/*Terminating null character*/
      /*Printing the equation parsed, NOTE: '\n' char is in the odestr[DIM-1]*/
      fprintf(stdout,"%s' = %s",var_name[DIM-1],odestr[DIM-1]);
      continue; /*We don't go further, as we have processed this line,
		  let's read another one*/
    }
    /****************************II: ode X=*********************************/
    /* if(temp[0]=='o' && temp[1]=='d' && temp[2]=='e'){ */
    /*   if(DIM==0) */
    /* 	var_name = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	var_name = (char **)realloc(var_name,(DIM+1)*sizeof(char *)); */
    /*   if(var_name==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */
    /*   var_name[DIM] = (char *)malloc(VAR_NAME_SIZE*sizeof(char)); */
    /*   if(var_name[DIM]==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */

    /*   if((eqsign-3) > (VAR_NAME_SIZE-1)){ */
    /* 	fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n", */
    /* 		(int)VAR_NAME_SIZE,temp[3],temp[4],temp[5]); */
    /* 	return 150; */
    /*   } */
    /*   for(i=3;i<eqsign;i++)/\*We get variable names*\/ */
    /* 	var_name[DIM][i-3] = temp[i]; */
    /*   var_name[DIM][i-3] = '\0';/\*Add terminating null character*\/ */
      
    /*   DIM++; /\*Let's increase dimension of the system by one*\/ */
    /*   /\* Allocating memory for ODE strings *\/ */
    /*   if(DIM==1) */
    /* 	odestr = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	odestr = (char **)realloc(odestr,DIM*sizeof(char *)); */
    /*   if(odestr==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr\n"); */
    /* 	return 101; */
    /*   } */
    /*   odestr[DIM-1] = (char *)malloc((strlen(temp)-eqsign)*sizeof(char)); */
    /*   if(odestr[DIM-1]==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr[DIM-1]\n"); */
    /* 	return 101; */
    /*   } */
    /*   for(i=eqsign+1;i<strlen(temp);i++)/\*Now we get ODE string*\/ */
    /* 	odestr[DIM-1][i-eqsign-1] = temp[i]; */
    /*   odestr[DIM-1][i-eqsign-1] = '\0';/\*Terminating null character*\/ */
    /*   /\*Printing the equation parsed, NOTE: '\n' char is in the odestr[DIM-1]*\/ */
    /*   fprintf(stdout,"%s' = %s",var_name[DIM-1],odestr[DIM-1]); */
    /*   continue; /\*We don't go further, as we have processed this line, */
    /* 		  let's read another one*\/ */
    /* } */
    /****************************III: X'=********************************/
    /* if(temp[eqsign-1]=='\''){/\*We have X'= notation*\/ */
    /*   if(DIM==0) */
    /*   	var_name = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	var_name = (char **)realloc(var_name,(DIM+1)*sizeof(char *)); */
    /*   if(var_name==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */
    /*   var_name[DIM] = (char *)malloc(VAR_NAME_SIZE*sizeof(char)); */
    /*   if(var_name[DIM]==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */

    /*   if((eqsign-1) > (VAR_NAME_SIZE-1)){ */
    /* 	fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n", */
    /* 		(int)VAR_NAME_SIZE,temp[0],temp[1],temp[2]); */
    /* 	return 150; */
    /*   } */
    /*   for(i=0;i<eqsign-1;i++)/\*We get variable names*\/ */
    /* 	var_name[DIM][i] = temp[i]; */
    /*   var_name[DIM][i] = '\0';/\*Add terminating null character*\/ */
      
    /*   DIM++; /\*Let's increase dimension of the system by one*\/ */
    /*   /\* Allocating memory for ODE strings *\/ */
    /*   if(DIM==1) */
    /*   	odestr = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	odestr = (char **)realloc(odestr,DIM*sizeof(char *)); */
    /*   if(odestr==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr\n"); */
    /* 	return 101; */
    /*   } */
    /*   odestr[DIM-1] = (char *)malloc((strlen(temp)-eqsign)*sizeof(char)); */
    /*   if(odestr[DIM-1]==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr[DIM-1]\n"); */
    /* 	return 101; */
    /*   } */
    /*   for(i=eqsign+1;i<strlen(temp);i++)/\*Now we get ODE specification */
    /* 	to the string*\/ */
    /* 	odestr[DIM-1][i-eqsign-1] = temp[i]; */
    /*   odestr[DIM-1][i-eqsign-1] = '\0';/\*Terminating null character*\/ */
    /*   /\*Printing the equation parsed, NOTE: '\n' char is in the odestr[DIM-1]*\/ */
    /*   fprintf(stdout,"%s' = %s",var_name[DIM-1],odestr[DIM-1]); */
    /*   continue; /\*We don't go further, as we have processed this line, */
    /* 		  let's read another one*\/       */
    /* } */
    /****************************IV: X(t)=*******************************/
    /* if(temp[eqsign-3]=='(' && temp[eqsign-2]=='t' */
    /*    && temp[eqsign-1]==')'){ */
    /*   if(DIM==0) */
    /* 	var_name = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	var_name = (char **)realloc(var_name,(DIM+1)*sizeof(char *)); */
    /*   if(var_name==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */
    /*   var_name[DIM] = (char *)malloc(VAR_NAME_SIZE*sizeof(char)); */
    /*   if(var_name[DIM]==NULL){ */
    /* 	fprintf(stderr,"Error: memory\n"); */
    /* 	return 100; */
    /*   } */

    /*   if((eqsign-3) > (VAR_NAME_SIZE-1)){ */
    /* 	fprintf(stderr,"Error: too many symbols(>%d) for var name: %c%c%c\n", */
    /* 		(int)VAR_NAME_SIZE,temp[0],temp[1],temp[2]); */
    /* 	return 150; */
    /*   } */
    /*   for(i=0;i<eqsign-3;i++)/\*We get variable names*\/ */
    /* 	var_name[DIM][i] = temp[i]; */
    /*   var_name[DIM][i] = '\0';/\*Add terminating null character*\/ */
      
    /*   DIM++; /\*Let's increase dimension of the system by one*\/ */
    /*   /\* Allocating memory for ODE strings *\/ */
    /*   if(DIM==1) */
    /* 	odestr = (char **)malloc(sizeof(char *)); */
    /*   else */
    /* 	odestr = (char **)realloc(odestr,DIM*sizeof(char *)); */
    /*   if(odestr==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr\n"); */
    /* 	return 101; */
    /*   } */
    /*   odestr[DIM-1] = (char *)malloc((strlen(temp)-eqsign)*sizeof(char)); */
    /*   if(odestr[DIM-1]==NULL){ */
    /* 	fprintf(stderr,"Error: memory odestr[DIM-1]\n"); */
    /* 	return 101; */
    /*   } */
    /*   for(i=eqsign+1;i<strlen(temp);i++)/\*Now we get ODE specification */
    /* 	to the string*\/ */
    /* 	odestr[DIM-1][i-eqsign-1] = temp[i]; */
    /*   odestr[DIM-1][i-eqsign-1] = '\0';/\*Terminating null character*\/ */
    /*   /\*Printing the equation parsed, NOTE: '\n' char is in the odestr[DIM-1]*\/ */
    /*   fprintf(stdout,"%s' = %s",var_name[DIM-1],odestr[DIM-1]); */
    /*   continue; /\*We don't go further, as we have processed this line, */
    /* 		  let's read another one*\/       */
    /* } */
    /***************************************************************/
    /*               END of ODE specification processing           */
    /***************************************************************/
    /***************************************************************/
    /*                          JACOBIAN                           */
    /***************************************************************/
    if(temp[0]=='j' && temp[1]=='a' && temp[2]=='c'){
      jvarc = 0;
      /*Allocating memory for strings. Assume DIM has its final value,
	which mean jac definition goes AFTER ALL ODE DEFINITIONS. It
	is extremely needed since otherwise it is not possible to
	define fully the jacobian, e.g. over variables that have not
	been defined so far.*/
      /* DIM must have its final value by reaching jacobian parsing. The */
      /* following strings check that */
      if(jdim==0)
	jdim = DIM;
      if((DIM!=jdim) && (jdim!=0)){
	fprintf(stderr,"ERROR: the jacobian cannot be complete! Exit.\n");
	return 456;
      }
      /* ********Initializing jacobian strings********** */
      if(!jaccount){
	jacstr = (char ***)malloc(DIM*sizeof(char **));
	for(i=0;i<DIM;i++)
	  jacstr[i] = (char **)malloc(DIM*sizeof(char *));
      }
      memset(temp1,'\0',STRLEN);
      for(i=3;i<eqsign;i++)/*Copy name of the variable*/
	temp1[i-3] = temp[i];
      if((jvar = isVar(temp1,var_name,DIM))>=0){/*DIM parameter must be
      changed if e.g. time var `t' is under consideration(not ODE
      case)*/
	/*Above IF determines which var Jac is being defined for. Thus
	  the sequence of jac statements in .ode file is not crucial, but
	  all of them MUST be after ALL ode statements.*/
	memset(temp1,'\0',STRLEN);
	j = eqsign+1;
	m = j;
	while(1){
	  if((char)temp[m]==',' || (char)temp[m]==';' || (char)temp[m]=='\n'){
	    if((char)temp[m+1]=='\n')/*If entry is empty then break*/
	      break;
	    temp1[m-j] = '\n';
	    temp1[m-j+1]='\0';
	    /*Allocating memory jacobian strings*/
	    jacstr[jvar][jvarc] = (char *)malloc((strlen(temp1)+1)*sizeof(char));
	    memset(jacstr[jvar][jvarc],'\0',strlen(temp1));
	    /*Copy temp1 to jacstr*/
	    strcpy(jacstr[jvar][jvarc],temp1);
	    /*Print jac strings: NOTE '\n' must be part of the string.*/
	    printf("d(%s')/d%s = %s",var_name[jvar],var_name[jvarc],jacstr[jvar][jvarc]);
	    jvarc++;
	    if((char)temp[m]=='\n')
	      break;
	    m++;
	    j = m;
	    continue;
	  }
	  temp1[m-j] = temp[m];
	  m++;
	}
	if(jvarc>DIM){
	  printf("Warning: you defined something more than DIM in jac\n");
	}
	if(jvarc<DIM){
	  printf("Warning: you have defined jac not completely\n");
	}
      }
      else{
	printf("ERROR: Unknown variable name for Jacobian: `%s'\n",temp1);
	return 365;
      }
      memset(temp1,'\0',STRLEN);
      memset(temp,'\0',STRLEN);
      jaccount++;
      continue;
    }
    if(temp[0]=='t' && temp[1]=='j' && temp[2]=='a' & temp[3]=='c'){
      /*Allocating memory for derivative of RHS over time (not ODE)*/
      tjacstr = (char **)malloc(DIM*sizeof(char *));
      memset(temp1,'\0',STRLEN);
      for(i=4;i<eqsign;i++)/*Copy name of the variable*/
	temp1[i-4] = temp[i];
      if((jvar = isVar(temp1,var_name,DIM))>=0){/*DIM parameter must be
      changed if e.g. time var `t' is under consideration(not ODE case)*/
	memset(temp1,'\0',STRLEN);
	j = eqsign+1;
	while((char)temp[j]!='\n'){
	  temp1[j-eqsign-1] = temp[j];j++;}
	temp1[j-eqsign-1] = '\n';
	temp1[j-eqsign] = '\0';
	tjacstr[jvar] = (char *)malloc((strlen(temp1)+1)*sizeof(char));
	memset(tjacstr[jvar],'\0',strlen(temp1));
	strcpy(tjacstr[jvar],temp1);
      }
      else 
	printf("Unknown variable name for Jacobian: `%s'\n",temp1);
      memset(temp1,'\0',STRLEN);
      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*                 END of JACOBIAN processing                  */
    /***************************************************************/
    /***************************************************************/
    /*                          FUNCTIONS                          */
    /***************************************************************/
    if(temp[eqsign-1] == ')'){
      /*Allocating memory for the next function to process.*/
      if(FNUM==0){
	function = (struct func *)malloc(sizeof(struct func));
	if(function == NULL){
	  fprintf(stderr,"Error: mem alloc for function (FNUM=%d)\n",FNUM);
	  return 12;
	}
      }
      else{
	function = (struct func *)realloc(function,(FNUM+1)*sizeof(struct func));
	if(function == NULL){
	  fprintf(stderr,"Error: mem alloc for function (FNUM=%d)\n",FNUM);
	  return 12;
	}
      }
      function[FNUM].num_arg = 0;
      /*Get the name of the function.*/
      /*Allocating memory for the function's name. (eqsign+1) length
	is enough and even more.*/
      function[FNUM].func_name = (char *)malloc((eqsign+1)*sizeof(char));
      i = 0;
      while(temp[i]!='('){
	function[FNUM].func_name[i] = temp[i];
	i++;
      }
      function[FNUM].func_name[i] = '\0';/*Trailing null character.*/
      if((isFunc(function[FNUM].func_name)!=FNUM) &&
	 (isFunc(function[FNUM].func_name)>=0))
	fprintf(stdout,"Warning: the function with the same name: `%s'\n",function[FNUM].func_name);
      j = i;/*j now stores position of the opening parenthesis '('*/
      /*temp1 is now free, we can use it.*/
      /*Copy list of arguments to temp1(without parentheses).*/
      memset(temp1,'\0',STRLEN);
      for(i=j+1;i<eqsign-1;i++){/*Copying*/
	temp1[i-j-1] = temp[i];
      }
      temp1[i-j-1] = '\0';/*Terminating null character.*/
      /*Get the list of arguments one by one.*/
      i = 0;
      j = 0;
      while(i<=strlen(temp1) && strlen(temp1)!=0){
	if(temp1[i]==',' || i==strlen(temp1)){
	  function[FNUM].num_arg++;/*Increase number of arguments.*/
	  if(function[FNUM].num_arg == 1)
	    function[FNUM].arg_list =
	      (char **)malloc(function[FNUM].num_arg*sizeof(char *));
	  else
	    function[FNUM].arg_list =
	      (char **)realloc(function[FNUM].arg_list,
			       function[FNUM].num_arg*sizeof(char *));
	  if(function[FNUM].arg_list==NULL)
	    fprintf(stderr,"Error: memory allocation (function.arg_list)\n");
	  function[FNUM].arg_list[function[FNUM].num_arg-1] =
	    (char *)malloc((i-j+1)*sizeof(char));
	  if(function[FNUM].arg_list[function[FNUM].num_arg-1]==NULL)
	    fprintf(stderr,"Error: memory allocation (function.arg_list[ArgNum-1])\n");
	  for(k=j;k<i;k++)/*Copying*/
	    function[FNUM].arg_list[function[FNUM].num_arg-1][k-j] = temp1[k];
	  function[FNUM].arg_list[function[FNUM].num_arg-1][k-j] = '\0';
	  j = i+1;/*Mark the start of the next arg right after ',' position)*/	  
	}
	i++;
      }
      /*Get right hand side string representing function.*/
      function[FNUM].func_str = (char *)malloc((strlen(temp)-eqsign-1)*sizeof(char));
      for(i=eqsign+1;i<strlen(temp);i++)/*Now we get RHS of the function*/
	function[FNUM].func_str[i-eqsign-1] = temp[i];
      function[FNUM].func_str[i-eqsign-1] = '\0';/*Terminating null character*/
      /**Rebuilding the functions.**/
      printf("Rebuilding the function:\n");
      printf("%s",function[FNUM].func_name);
      printf("(");
      if(strlen(temp1)==0)
	printf(")");
      for(i=0;i<function[FNUM].num_arg;i++){
	if(i==function[FNUM].num_arg-1){
	  printf("%s)",function[FNUM].arg_list[i]);
	  break;
	}
	printf("%s,",function[FNUM].arg_list[i]);
      }
      printf("=");
      printf("%s",function[FNUM].func_str);
      /**End of rebuilding.**/
      FNUM++;/*Increase number of functions which have been read so far.*/
      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*                  End of FUNCTIONS processing                */
    /***************************************************************/
    /***************************************************************/
    /*                     AUXILLARY entities                      */
    /***************************************************************/
    if(temp[0]=='a' && temp[1]=='u' && temp[2]=='x'){
      /*Allocating memory for a next aux to process.*/
      if(AUXNUM==0)
	aux = malloc(sizeof(struct aux));
      else
	aux = realloc(aux,(AUXNUM+1)*sizeof(struct aux));
      /*Memory for the name of the aux. eqsign chars is more than enough.*/
      aux[AUXNUM].name = (char *)malloc(eqsign*sizeof(char));
      for(i=3;i<eqsign;i++)
	aux[AUXNUM].name[i-3] = temp[i];
      aux[AUXNUM].name[i-3] = '\0';/*Terminating null character.*/
      /*Get formula for the aux.*/
      /*Memory for the string*/
      aux[AUXNUM].formula = (char *)malloc((strlen(temp)-eqsign-1)*sizeof(char));
      for(i=eqsign+1;i<strlen(temp);i++)/*Get formula*/
	aux[AUXNUM].formula[i-eqsign-1] = temp[i];
      aux[AUXNUM].formula[i-eqsign-1] = '\0';/*Terminating null character.*/
      /*Rebuilding the aux*/
      printf("Rebuilding auxillary:\n");
      printf("%s",aux[AUXNUM].name);
      printf("=");
      printf("%s",aux[AUXNUM].formula);
      /*End of rebuilding.*/
      AUXNUM++;
      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*                 End of AUXILLARY entities                   */
    /***************************************************************/
    /*PARAMETERS*/
    /*Structural parameters of the system being defined*/
    if(unstrp[0]=='p' && (unstrp[1]==' ' || unstrp[3]==' ' || unstrp[5]==' ')){
      /*We use unstriped lines. So searching for the spaces.*/
      process_par(temp);
      memset(temp,'\0',STRLEN);
      memset(unstrp,'\0',STRLEN);
      continue;
    }
    /*Techniacal parameters*/
    if(temp[0]=='@'){/*Internal parameters.*/
      if((process_inter_par(temp))!=0)
	return 45;
      memset(temp,'\0',STRLEN);
      continue;
    }
    /*INITIALs*/
    if(unstrp[0]=='i' && (unstrp[1]==' ' || unstrp[4]==' ')){/*Initial conditions*/
      initcount++;
      if(initcount == 1)
	initstr = (char **)malloc(sizeof(char *));
      else
	initstr = (char **)realloc(initstr,initcount*sizeof(char *));
      initstr[initcount-1] = (char *)malloc(STRLEN*sizeof(char));
      strcpy(initstr[initcount-1],temp);/*Copy to buffer of initial strings*/
      memset(temp,'\0',STRLEN);
      memset(unstrp,'\0',STRLEN);
      continue;
    }
    /****************************************************************/
    /*                    Gillespie specifications                  */
    /****************************************************************/
    if((temp[0]=='g' && temp[1]=='i' && temp[2]=='l' && temp[3]=='l')
       || (temp[0]=='g')){
      /*Gillespie specifications MUST be after all ODE specifications,
	since we use here DIM*/
      NREAC++;
      for(i=0;i<strlen(temp);i++){
	if(temp[i]==':') colon = i;
	if(temp[i]==';') semicolon = i;
      }
      /*Propensity strings*/
      /* if(NREAC==1) */
      /* 	propstr=(char **)malloc(sizeof(char *)); */
      /* else */
      propstr=realloc(propstr,NREAC*sizeof(char *));
      if(!propstr){
	fprintf(stderr,"Error: could not allocate memory for propstr\n");
	free(propstr);propstr=0;
	return 1345;
      }
      propstr[NREAC-1]=(char *)malloc((semicolon-colon+1)*sizeof(char));
      if(!propstr[NREAC-1]){
	fprintf(stderr,"Error: couldn't get memory for propstr[NREAC-1]\n");
	return 1344;
      }
      /*Copying propensity string*/
      memcpy(propstr[NREAC-1],(temp+colon+1),(semicolon-colon-1));
      propstr[NREAC-1][semicolon-colon-1] = '\n';
      propstr[NREAC-1][semicolon-colon] = 0;
      if((strlen(temp)-semicolon) == 2){/*No vector is specified*/
	fprintf(stderr,"Error: no update vector found:\n");
	fprintf(stderr,"`%s'",temp);
	return 1;
      }
      /* if(NREAC==1) */
      /* 	upvec=(int **)malloc(sizeof(int *)); */
      /* else */
      upvec=(int **)realloc(upvec,NREAC*(sizeof(int *)));
      if(!upvec){
	fprintf(stderr,"Error: could not allocate memory for upvec\n");
	free(upvec);upvec=0;
	return 1341;
      }
      upvec[NREAC-1]=(int *)calloc(DIM,sizeof(int));/* We want clearance of memory */
      if(!upvec[NREAC-1]){
	fprintf(stderr,"Error: couldn't get memory for upvec[NREAC-1]\n");
	return 1340;
      }
      /*Get update vector*/
      process_update_vector((temp+semicolon+1),upvec[NREAC-1]);
      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*          LANGEVIN amendments, i.e. stochastic terms         */
    /***************************************************************/
    if(!strncmp(temp,"lang",4)){
      if((strchr(temp,'\n')) != NULL)/* Remove the newline */
	*(strchr(temp,'\n')) = '\0';
      if(NLANG == 0){/* Init for DIM strings once */
	lang_amend_str = (char **)malloc(DIM*sizeof(char *));
	if(!lang_amend_str){
	  fprintf(stderr,"Error: no space for lang_amend_str.\n");
	  free(lang_amend_str); return 4252;
	}
	for(i=0;i<DIM;i++){
	  /* Init to zeros to be empty */
	  lang_amend_str[i] = (char *)calloc(STRLEN,sizeof(char));
	  if(lang_amend_str[i] == NULL){
	    fprintf(stderr,"Error: no space for lang_amend_str[%d]\n",i);
	    free(lang_amend_str);
	    return 4253;
	  }
	}
      }
      /* String format is: lang x=0.1, y=0.01 */
      p1 = temp + 4;
      while((p = strsep(&p1,",")) != NULL){
	i = 0;
	while((p2 = strsep(&p,"=")) != NULL){
	  if(i == 0){
	    j = isVar(p2,var_name,DIM);
	    printf("lang(%s) = ",var_name[j]);
	  }
	  if(i == 1){
	    /* Copy the RHS */
	    lang_amend_str[j] = strcpy(lang_amend_str[j],p2);
	    /* Append the newline, needed for parsing */
	    lang_amend_str[j] =	strcat(lang_amend_str[j],"\n");
	    printf("%s",lang_amend_str[j]);
	  }
	  i++;
	}
	NLANG++;
      }
      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*                      FIXED ENTITIES                         */
    /***************************************************************/
    if(neqsigns==1){/*This is very general rule, but is restricted to
      the case of one equal sign in the string. Also, this conditional
      is placed after all others because of its generality.*/
      /*WE TREAT FIXED ENTITIES AS AUXILLARY.*/
      /*Allocating memory for a next aux to process.*/
      if(AUXNUM==0)
	aux=malloc(sizeof(struct aux));
      else
	aux=realloc(aux,(AUXNUM+1)*sizeof(struct aux));
      /*Memory for the name of the aux. eqsign chars is more than enough.*/
      aux[AUXNUM].name=(char *)malloc(eqsign*sizeof(char));
      for(i=0;i<eqsign;i++)
	aux[AUXNUM].name[i]=temp[i];
      aux[AUXNUM].name[i]='\0';/*Terminating null character.*/
      /*Get formula for the aux.*/
      /*Memory for the string*/
      aux[AUXNUM].formula=(char *)malloc((strlen(temp)-eqsign-1)*sizeof(char));
      for(i=eqsign+1;i<strlen(temp);i++)/*Get formula*/
	aux[AUXNUM].formula[i-eqsign-1]=temp[i];
      aux[AUXNUM].formula[i-eqsign-1]='\0';/*Terminating null character.*/

      /*Rebuilding the aux*/
      printf("Rebuilding fixed entity:\n");
      printf("%s",aux[AUXNUM].name);
      printf("=");
      printf("%s",aux[AUXNUM].formula);
      /*End of rebuilding.*/
      AUXNUM++;

      memset(temp,'\0',STRLEN);
      continue;
    }
    /***************************************************************/
    /*                    END of FIXED ENTITIES                    */
    /***************************************************************/
    fprintf(stderr,"Unknown line `%s'\n",temp);
    memset(temp,'\0',STRLEN);
  }
  /****************End of input reading***************************/
  /***********ERROR: no equations*********************************/
  if(DIM==0){
    fprintf(stderr,"ERROR: I have read NO equations.\n");
    return 478;
  }
  if(jaccount<DIM && jaccount>0){
    fprintf(stderr,"ERROR: incomplete Jacobian.\n");
    return 290;
  }
  if(jaccount>0)
    jac_flag=1;/* Setting up the flag, Jac is specified in the source */
  /*****************************************************************/
  free(temp);free(temp1);free(unstrp);/*Freed memory of temporary strings*/
  fclose(file);/*Close stream associated with .ode source file*/
  /************Preparing initial conditions****************/
  /* This is done after DIM and PARDIM are known */
  xin = (double *)malloc(DIM*sizeof(double));
  yin = (long int *)malloc(DIM*sizeof(long int));
  xin_par = (int *)malloc(DIM*sizeof(int));
  /*Make default initial conditions: zeros(xin and yin) and -1(xin_par), which is
      overwritten by init statements.*/
  for(i=0;i<DIM;i++){
    xin[i] = 0.0;
    yin[i] = 0;
    xin_par[i] = -1;
  }
  for(i=0;i<initcount;i++)
    process_init(initstr[i]);
  /*************************************************************/
  /*Processing of the source file .ode to output source file .c*/
  /*************************************************************/
  if(!No_Transfer)/* Transfer only if there is no .c source with system */
    process_odefile(basename);
  /************************/
  /*Let's print par names detected*/
  printf("Active parameters: ");
  for(i=0;i<PARDIM;i++)
    fprintf(stdout,"`%s', ",par_name[i]);
  fprintf(stdout,"\n");
  
  /*************************************************************/
  /*****************************************************************/
  /* Forming configuration file*/
  /*****************************************************************/
  gen_conf(basename);
  return 0;
}

int process_odefile(char const basename[])
{
  int i,j,k;
  /**************************************************************/
  /*Form C names list for variables:
   continuous(var_c_name) and discrete(disc_var_c_name).
   In Dinamica, the former ones are `x[i]' and the latter ones are
   `y[i]' where i is the index of the variable. Range is the same from
   1 to DIM for both types of variables. IMPORTANT: this implies that
   only systems with the same dimension can be implemented.*/
  char **var_c_name,**disc_var_c_name;
  var_c_name=(char **)malloc(DIM*sizeof(char *));
  disc_var_c_name=(char **)malloc(DIM*sizeof(char *));
  for(i=0;i<DIM;i++){
    var_c_name[i]=(char *)malloc(VAR_NAME_SIZE*sizeof(char));
    sprintf(var_c_name[i],"x[%d]",i);
    disc_var_c_name[i]=(char *)malloc(VAR_NAME_SIZE*sizeof(char));
    sprintf(disc_var_c_name[i],"y[%d]",i);
  }
  /**************************************************************/
  
  char *c_src_filename;
  //c_src_filename=(char *)calloc((strlen(basename)+2),sizeof(char));
  c_src_filename=app_ext(basename,"c",c_src_filename);
  if(c_src_filename==NULL)
    fprintf(stderr,"Error: allocating memory for `filename'\n");
  
  FILE * c_src_file;
  /*Open stream*/
  c_src_file=fopen(c_src_filename,"w");
  /*Error check*/
  if(c_src_file==NULL)
    fprintf(stderr,"Error: open file `%s'.\n",c_src_filename);
  /******************************************************************/
  /* TRANSFER ODE system,functions and auxillaries etc. to source .c file */
  /******************************************************************/
  fprintf(stdout,"Starting transfer to `%s'...\n",c_src_filename);
  /*PRINT some heading information to the .c source file, like
    Copyright, Non-Warranty, License*/
  print_heading(c_src_file);
  /* PRINT body of the file */
  fprintf(c_src_file,"#include <math.h>\n");
  fprintf(c_src_file,"/*\n <Space for user comments>\n*/\n\n\n");
  /*PRINT functions...*/
  for(i=0;i<FNUM;i++){
    if((function[i].num_arg)!=0){
      fprintf(c_src_file,"double func_%s(double t, double P[], double %s",function[i].func_name,function[i].arg_list[0]);}
    if((function[i].num_arg)==0){
      fprintf(c_src_file,"double func_%s(double t, double P[])\n",function[i].func_name);
    }
    for(j=1;j<(function[i].num_arg);j++){
      if((function[i].arg_list[j])==NULL) break;
      fprintf(c_src_file,", double %s",(function[i].arg_list[j]));
    }

    if((function[i].num_arg)!=0)
      fprintf(c_src_file,")\n");/*Closing parenthesis*/
    fprintf(c_src_file,"{\n");
    fprintf(c_src_file," return ");
    process_odestr(1,c_src_file,function[i].func_str,function[i].arg_list,function[i].arg_list,function[i].num_arg);
    fprintf(c_src_file,"}\n\n");
  }
  fprintf(c_src_file,"\n");
  /* PRINT LANGEVIN AMENDMENTS */
  fprintf(c_src_file,"int lang_amend(double t, const double x[], ");
  fprintf(c_src_file,"double ksi[], double P[])\n");
  fprintf(c_src_file,"{\n");
  for(i=0;i<DIM;i++){
    fprintf(c_src_file," ksi[%d] = ",i);
    /* if zero is the 1st symbol, nothing is specified OR no lang amendments at all */
    if((NLANG == 0) || (lang_amend_str[i][0] == '\0'))
      fprintf(c_src_file," 0;\n");
    else
      process_odestr(0,c_src_file,lang_amend_str[i],var_name,var_c_name,DIM);
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /* PRINT RHS... */
  fprintf(c_src_file,"int rhs(double t, const double x[], ");
  fprintf(c_src_file,"double f[], double P[])\n");
  fprintf(c_src_file,"{\n");
  for(i=0;i<DIM;i++){/*THE MAINLY MAIN LOOP: Looping over dimension */
    fprintf(c_src_file," f[%d] = ",i);
    //    printf("I am going to process:\n %s",odestr[i]);
    /*Next function provides the MAIN LOOP over characters in odestr[]*/
    process_odestr(0,c_src_file, odestr[i],var_name,var_c_name,DIM);
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /*PRINT auxillary...*/
  fprintf(c_src_file,"int auxillary(double t, double P[], double x[], double * aux)\n");
  fprintf(c_src_file,"{\n");
  for(i=0;i<AUXNUM;i++){
    fprintf(c_src_file," aux[%d] = ",i);
    process_odestr(0,c_src_file,aux[i].formula,var_name,var_c_name,DIM);
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /*PRINT Jacobian... */
  fprintf(c_src_file,"int jacobian(double t, const double x[], ");
  fprintf(c_src_file,"double dfdx[], double dfdt[], double P[])\n");
  fprintf(c_src_file,"{\n");
  if(jac_flag){
    for(i=0;i<DIM;i++){
      for(j=0;j<DIM;j++){
	if(jacstr!=NULL){
	  fprintf(c_src_file," dfdx[%d] = ",(j+i*DIM));
	  process_odestr(0,c_src_file,jacstr[i][j],var_name,var_c_name,DIM);}
	else
	  fprintf(c_src_file," dfdx[%d] = 0;\n",(j+i*DIM));
      }
      if(tjacstr!=NULL){
	fprintf(c_src_file," dfdt[%d] = ",i);
	process_odestr(0,c_src_file,tjacstr[i],var_name,var_c_name,DIM);}
      else
	fprintf(c_src_file," dfdt[%d] = 0;\n",i);
    }
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /*PRINT propensity computing function*/
  fprintf(c_src_file,"int propensity(const double y[], ");
  fprintf(c_src_file,"const double P[], double A[])\n{\n");
  for(i=0;i<NREAC;i++){
    fprintf(c_src_file," A[%d] = ",i);
    process_odestr(0,c_src_file,propstr[i],var_name,disc_var_c_name,DIM);
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /*PRINT the function updating chemical species*/
  fprintf(c_src_file,"int update(long int *y, const int mu)\n{\n");
  for(i=0;i<NREAC;i++){
    fprintf(c_src_file," if(mu==%d){\n",i+1);
    for(j=0;j<DIM;j++){
      if(upvec[i][j]!=0)
	fprintf(c_src_file,"  y[%d] = y[%d] + (%d);\n",j,j,upvec[i][j]);
    }
    fprintf(c_src_file," }\n");
  }
  fprintf(c_src_file,"\n return 0;\n}\n\n");
  /******************************************************************/
  /*                     END of TRANSFER                            */
  /******************************************************************/
  /* Print some helpful tail information to the .c source file for
     improving the user's readability*/
  print_tail(c_src_file);

  fclose(c_src_file);
  fprintf(stdout,"Done.\n");

  return 0;
}

int process_odestr(int ENCLOSE,FILE *c_src_file, const char * odestr,
		   char **var_list, char **var_c_list, int var_list_size)
{
  int i,j,k,m,nick;/*Working variables*/
  char *symbuf;/*Symbols buffer*/
  char *temp;/*Temporary string*/
  size_t length;/*Length of the symbol buffer*/
  int hat=0;/*Position of the '^' symbol,see below this special case*/
  /*String is empty, transfer just zero.*/
  if(odestr==NULL){
    fprintf(c_src_file,"0;\n");
    return 0;
  }
  /*Before any transfer. If you want to surround your expression with
    parentheses just do it with positive value of ENCLOSE.*/
  if(ENCLOSE)
    fprintf(c_src_file,"(");
  /********************************************************/
  if(odestr[0]=='^'){/*Error: we have first '^'*/
    fprintf(stderr,"Error: processing ODE `");
    k=0;/*Error reporting output*/
    while(odestr[k]!='\n'){
      fprintf(stderr,"%c",odestr[k]);k++;}
    fprintf(stderr,"': ");
    for(k=0;k<3;k++){
      if(odestr[k]!='\n')
	fprintf(stderr,"%c",odestr[k]);}
    fprintf(stderr,"...\n");
    fprintf(stderr,"Quit.(1)\n");/*Return code to user*/
    return 1;/* Bad equation */
  }
  for(j=0;j<strlen(odestr);j++){/*THE MAIN LOOP: Looping over characters in the odestr*/
    if(odestr[j]=='^'){/* '^'-- special case as we should replace
			     it with pow(<left>,<right>) function
			     from math.h. <left>, <right> -- left
			     and right operands respectively. */
      if(j==hat){/*This means we are at the '^' again for the second
		   time. It is possible when we dealt <left> complex
		   expressions (i.e.in parenthesis) OR it is a real
		   error*/
	/*****************<LEFT> is complex expression****************/
	if(odestr[j-1]==')'){/*Previous run was on complex
				  expression, Then We need <right>
				  processing to run*/
	  /*Transfer COMMA between <left>(complex expression we came
	    from) and <right> that will be processed soon*/
	  fprintf(c_src_file,",");
	  /*Allocating memory for the buffer*/
	  length=(size_t)40;
	  symbuf=(char *)calloc(length,sizeof(char));
	  if(symbuf==NULL){
	    fprintf(stderr,"Error: memory allocating for temporary buffer: ");
	    fprintf(stderr,"<right> complex, <left> complex.\n");}
	  j=hat+1;/*Step forward*/
	  if(odestr[j]=='('){/*Another complex expression on the <right>*/
	    j++;/*Jump one position forward from '('*/
	    m=1;/*m is counter for matched pairs `('/`)'*/
	    k=0;/*For copying*/
	    while((odestr[j]!=')' && m!=0) || (odestr[j]==')' && m!=0)){
	      if(k==0 && odestr[j]=='(') m++;
	      if(k==0 && odestr[j]==')') m--;
	      symbuf[k]=odestr[j];/*Copying*/
	      if(k==length-1){
		length+=40;
		symbuf=(char *)realloc(symbuf,length*sizeof(char));
		if(symbuf==NULL){
		  fprintf(stderr,"Error: memory reallocating for temporary buffer: ");
		  fprintf(stderr,"<right> complex, <left> complex.\n");}
	      }
	      j++;k++;
	      if(odestr[j]=='(')/*Another open `(' inside other pair
				  `('/`)'*/
		m++;
	      if(odestr[j]==')') m--;/*Counter -- we found one close `)'*/
	      if(odestr[j]=='\n'){ /*Error we got end of the ODE
					  specification and
					  still no closing
					  `)'*/
		fprintf(stderr,"Error: no closing `)'\n");break;}
	    }
	    if(odestr[j]=='\n') break; /*Break from the main loop*/
	    symbuf[k]='\0';
	    /*Transfer opening `('*/
	    fprintf(c_src_file,"(");
	    /*NOW, processing this string...*/
	    process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	    /********************************/
	    /*Transfer closing parenthesis `)' of <right> complex expr*/
	    fprintf(c_src_file,")");
	    /*Transfer closing parenthesis `)' of pow(_,_) function*/
	    fprintf(c_src_file,")");
	    /*j is still on the closing `)', so step back to treat it
	      in the next iteration of the main loop*/
	    continue;
	  }
	  else{/*If after '^' we have not a opening `('*/
	    length=0;
	    while(isalnum(odestr[j]) || odestr[j]=='.'){
	      length++;j++;}
	    /*Allocating memory for the buffer*/
	    symbuf=(char *)calloc(length,sizeof(char));
	    if(symbuf==NULL){
	      fprintf(stderr,"Error: memory allocating for temporary buffer: ");
	      fprintf(stderr,"<left> complex, <right> simple.\n");}
	    /*Copying buffer to symbuf, which we use*/
	    for(k=0;k<length;k++)
	      symbuf[k]=odestr[hat+1+k];
	    symbuf[k]='\0';
	    /*Again processing this string*/
	    process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	    /********************************/
	    /*Transfer closing `)' of pow(_,_) function*/
	    fprintf(c_src_file,")");
	    j--;/*As we have j at one position forward to the end of
		  selected buffer symbuf*/
	    continue;
	  }
	}
	/* Following when error does occur*/
	fprintf(stderr,"Error: I will not process the same '^' case twice.\n");
	break;
      }/*END OF TREATING <left> complex expressions*/
      /*********TRANSFER `pow(' beginning of the pow function************/
      fprintf(c_src_file,"pow(");
      /******************************************************************/
      hat=j;/*Store position of the '^' hat symbol*/
      /*Watching BACK*/
      j--;/*Initially switch back one position*/
      if(odestr[j]==')'){/*If we have complex expression
			      powered*/
	/*Memory*/
	length=(size_t)40; symbuf=(char *)calloc(length,sizeof(char));
	if(symbuf==NULL){
	  fprintf(stderr,"Error: memory allocating for temporary buffer: ");
	  fprintf(stderr,"<left> complex.\n");}
	m=1;/*Matching*/
	while((odestr[j]!='(' && m!=0) || (odestr[j]=='(' && m!=0)){/*Looping until get open '('*/
	  j--;/*Still jumping back*/
	  if(odestr[j]=='(') m--;
	  if(odestr[j]==')') m++;
	  if(j==-1) break;/*Get out of this while loop*/
	}
	if(j==-1){/*Error: let's get out of the main loop*/
	  fprintf(stderr,"Error: no opening `(': ");
	  fprintf(stderr," I searched till the beginning of the `%s'.\n",odestr);
	  break;
	}
	m=1;/*Matching algorithm*/
	k=0;
	j++;/*Stepping forward from opening `('*/
	while((odestr[j]!=')' && m!=0) || (odestr[j]==')' && m!=0)){
	  if(k==0 && odestr[j]=='(') m++;
	  if(k==0 && odestr[j]==')') m--;
	  symbuf[k]=odestr[j];/*Copy*/
	  if(k==length-1){
	    length+=40;
	    symbuf=(char *)realloc(symbuf,length*sizeof(char));
	    if(symbuf==NULL){
	      fprintf(stderr,"Error: memory reallocating: ");
	      fprintf(stderr,"<left> complex.\n");}
	  }
	  k++;j++;
	  if(odestr[j]==')') m--;
	  if(odestr[j]=='(') m++;
	}
	symbuf[k]='\0';/*Terminating null character*/
	/*Transfer opening `(' of complex expression*/
	fprintf(c_src_file,"(");
	/* Process complex expr*/
	process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	/********************************/
	/*Transfer terminating close `)' of <left> complex expr*/
	fprintf(c_src_file,")");
	/*Now j is at closing `)', next is hat '^', that's OK, continue*/
	continue;
      }else {/*********<left> expr is SIMPLE*********/
	while((isalnum(odestr[j]) || odestr[j]=='.')){
	  j--;/*Jump through char or digits BACK; j==-1 condition for
		cases where we've passed through the
		beginning of the odestr buffer*/
	  if(j==-1) break;
	}
	for(k=j+1;k<hat;k++){/*Copy <left> symbol*/
	  if(k-j-1>length-1){/*Reallocating memory in case of exhausting*/
	    length+=40;
	    symbuf=(char *)realloc(symbuf,length*sizeof(char));
	    if(symbuf==NULL){/*Error handling*/
	      fprintf(stderr,"Error: memory reallocating: ");
	      fprintf(stderr,"<left> simple.\n");
	    }
	  }
	  symbuf[k-j-1]=odestr[k];
	}
	symbuf[k]='\0';/*Terminating null character*/
	/**************************FUNCTIONS***************************/
	/*                          F(X)^Y                           */
	if(function!=NULL && (isFunc(symbuf))>=0 && odestr[j]=='('){
	  /*We have function with the same name and must be sure that
	    the next symbol is '(': as all functions should appear.*/
	  k=j;/*Mark the position.*/
	  length=40;/*string length.*/
	  temp=(char *)malloc(length*sizeof(char));
	  fprintf(c_src_file,"(func_%s(t,P",function[isFunc(symbuf)].func_name);
	  j++;/*Step forward*/
	  /*Store in temp arguments' list(without parentheses).*/
	  while(odestr[j]!=')'){
	    /*This while loop also gets j at the correct position for the
	      next main loop iterations.*/
	    temp[j-k-1]=odestr[j];
	    if((j-k-1)==(length-1)){/*Reallocating in case of overbuffering.*/
	      length+=length;
	      temp=(char *)realloc(temp,length*sizeof(char));
	    }
	    j++;
	  }
	  temp[j-k-1]='\0';/*trailing null char.*/
	  if(strlen(temp)==0){
	    /*If no arguments to the function.*/
	    fprintf(c_src_file,"))");
	    continue;
	  }
	  /*If there are some arguments to the function we proceed with
	    the comma "," delimiting first default arguments -- t and
	    P[] from other that follow.*/
	  fprintf(c_src_file,",");
	  /*We do not need symbuf: now it will hold the arguments one by
	    one of the function.*/
	  length=40;
	  symbuf=(char *)realloc(symbuf,length*sizeof(char));
	  /*Get the list of arguments one by one.*/
	  k=0;i=0;
	  while(k<=strlen(temp)){
	    if(temp[k]==',' || k==strlen(temp)){
	      for(m=i;m<k;m++){/*Copying to symbuf.*/
		symbuf[m-i]=temp[m];
		if((m-i)==(length-1)){/*Reallocating*/
		  length+=length;
		  symbuf=(char *)realloc(symbuf,length*sizeof(char));
		}
	      }
	      symbuf[m-i]='\0';
	      i=k+1;/*Hold new argument start position.*/
	      /*Normal processing of the argument.*/
	      process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	      if(temp[k]==',')
		fprintf(c_src_file,",");/*Transfer delimiting comman ","*/
	    }
	    k++;
	  }
	  /*Transfer closing parentheses of the function.*/
	  fprintf(c_src_file,"))");
	  /*Transfer comma `,' of pow(_,_) function*/
	  fprintf(c_src_file,",");
	  /*Without `continue' statement.*/
	}
	/**************************END FUNCTIONS***************************/
	if(isalpha(symbuf[0])){/*Then we are dealing with name...*/
	  if((isVar(symbuf,var_list,var_list_size))>=0)/* If we met var name */
	    fprintf(c_src_file,"%s",var_c_list[isVar(symbuf,var_list,var_list_size)]);/*Transfer*/
	  else if((isVar(symbuf,var_list,var_list_size))<0){/* If we did not meet var name */
	    if(!PARDIM){/*If we didnt find par name yet */
	      for(k=0;k<strlen(symbuf);k++)/* Copy par name */
		par_name[PARDIM][k]=symbuf[k];
	      fprintf(c_src_file,"P[%d]",PARDIM);/*transfer*/
	      PARDIM++;/* Increase number of par */
	    }
	    else if(PARDIM && (isPar(symbuf,PARDIM))>=0){/*If we've already found some
							   par names and newly
							   found name is par which
							   we met before*/
	      fprintf(c_src_file,"P[%d]",isPar(symbuf,PARDIM));
	    }
	    else if(PARDIM && (isPar(symbuf,PARDIM))<0) {/* Newly found name is really NEW par name */
	      for(k=0;k<strlen(symbuf);k++)
		par_name[PARDIM][k]=symbuf[k];
	      fprintf(c_src_file,"P[%d]",PARDIM);
	      PARDIM++;/* Increase number of par */
	    } else fprintf(stderr,"Error: I dont know what else is possible.\n");
	  }
	  else{fprintf(stderr,"Error: could not determine what kind of <left> symbol ");
	    fprintf(stderr,"`%s'", symbuf);break;}
	}
	if(isdigit(symbuf[0]) || symbuf[0]=='.'){/*Then we are dealing with number*/
	  /*Transfer*/
	  for(k=0;k<strlen(symbuf);k++)
	    fprintf(c_src_file,"%c",symbuf[k]);
	}
	if((!isalnum(symbuf[0]))){/*Error: neither char nor digit*/
	  fprintf(stderr,"Error: Neither char nor digit: ");
	  fprintf(stderr,"<left> simple.\n");
	}
      }
      /*Transfer `,'(comma) delimiting arguments of pow(_,_) function
	in case of simple <left> expression*/
      fprintf(c_src_file,",");
      /*Without `continue' statement.*/
      /***************************************************************/
      j=hat;/*Square one: return back from where we started to read BACK*/
      /*Watching FORWARD*/
      memset(symbuf,'\0',length);
      j++;/*Initially switch forward one position*/
      if(odestr[j]=='('){/*If we have something powered by complex
			      expression*/
	m=1;/*Match pair of parenthesis `('/`)'*/
	k=0;
	j++;/*Switch forward from `(' to next position*/
	while((odestr[j]!=')' && m!=0) || (odestr[j]==')' && m!=0)){/*Looping until get close ')'*/
	  if(k==0 && odestr[j]=='(') m++;
	  if(k==0 && odestr[j]==')') m--;
	  symbuf[k]=odestr[j];/*Copy the complex <right> expression*/
	  if(k==length-1){
	    length+=40;
	    symbuf=(char *)realloc(symbuf,length*sizeof(char));
	    if(symbuf==NULL){
	      fprintf(stderr,"Error: memory reallocating: ");
	      fprintf(stderr,"<right> complex, <left> simple\n");
	    }
	  }
	  k++;
	  j++;/*Still jumping forward*/
	  if(odestr[j]=='(') m++;/*Matching algorithm*/
	  if(odestr[j]==')') m--;/*Matching algorithm*/
	  if(odestr[j]=='\n') break;/*Get out of this while loop*/
	}
	symbuf[k]='\0';/*Terminating null character*/
	if(odestr[j]=='\n'){/*Error: let's get out of the main loop*/
	  fprintf(stderr,"Error: no closing `)': ");
	  fprintf(stderr,"I searched till the end of the `%s'.\n",odestr);
	  break;
	}
	/*Transfer opening `(' of complex expression*/
	fprintf(c_src_file,"(");
	/*****Processing the complex expression on the <right>*****/
	process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	/**********************************************************/
	/*Transfer closing parenthesis `)' of <right> complex expression*/
	fprintf(c_src_file,")");
	/*Transfer closing parenthesis `)' of pow(_,_) function*/
	fprintf(c_src_file,")");
	/*Now j is at closing `)', so next position will be treated in the
	  next iteration of the main loop, continue..*/
	continue;
      } else {/***********<right> expr is SIMPLE***************/
	while(isalnum(odestr[j]) || odestr[j]=='.') j++;/*Jump through char or digits FORWARD*/
	for(k=hat+1;k<j;k++)
	  symbuf[k-hat-1]=odestr[k];
	symbuf[k-hat-1]='\0';/*Trailing null charactter.*/

	/**************************FUNCTIONS***************************/
	/*                           Y^F(X)                           */
	if(function!=NULL && (isFunc(symbuf))>=0 && odestr[j]=='('){
	  /*We have function with the same name and must be sure that
	    the next symbol is '(': as all functions should appear.*/
	  k=j;/*Mark the position.*/
	  length=40;/*string length.*/
	  temp=(char *)malloc(length*sizeof(char));
	  fprintf(c_src_file,"(func_%s(t,P",function[isFunc(symbuf)].func_name);
	  j++;/*Step forward*/
	  /*Store in temp arguments' list(without parentheses).*/
	  while(odestr[j]!=')'){
	    /*This while loop also gets j at the correct position for the
	      next main loop iterations.*/
	    temp[j-k-1]=odestr[j];
	    if((j-k-1)==(length-1)){/*Reallocating in case of overbuffering.*/
	      length+=length;
	      temp=(char *)realloc(temp,length*sizeof(char));
	    }
	    j++;
	  }
	  temp[j-k-1]='\0';/*trailing null char.*/
	  if(strlen(temp)==0){
	    /*If no arguments to the function.*/
	    fprintf(c_src_file,"))");
	    continue;
	  }
	  /*If there are some arguments to the function we proceed with
	    the comma "," delimiting first default arguments: t and P[].*/
	  fprintf(c_src_file,",");
	  /*We do not need symbuf: now it will hold the arguments one by
	    one of the function.*/
	  length=40;
	  symbuf=(char *)realloc(symbuf,length*sizeof(char));
	  /*Get the list of arguments one by one.*/
	  k=0;i=0;
	  while(k<=strlen(temp)){
	    if(temp[k]==',' || k==strlen(temp)){
	      for(m=i;m<k;m++){/*Copying to symbuf.*/
		symbuf[m-i]=temp[m];
		if((m-i)==(length-1)){/*Reallocating*/
		  length+=length;
		  symbuf=(char *)realloc(symbuf,length*sizeof(char));
		}
	      }
	      symbuf[m-i]='\0';
	      i=k+1;/*Hold new argument start position.*/
	      /*Normal processing of the argument.*/
	      process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	      if(temp[k]==',')
		fprintf(c_src_file,",");/*Transfer delimiting comman ","*/
	    }
	    k++;
	  }
	  /*Transfer closing parentheses of the function.*/
	  fprintf(c_src_file,"))");
	  /*Transfer closing parenthesis `)' of pow(_,_) function*/
	  fprintf(c_src_file,")");
	
	  continue;
	}
	/**************************END FUNCTIONS***************************/
	
	if(isalpha(symbuf[0])){/*Then we are dealing with name...*/
	  if((isVar(symbuf,var_list,var_list_size))>=0)/* If we met var name */
	    fprintf(c_src_file,"%s",var_c_list[isVar(symbuf,var_list,var_list_size)]);/*Transfer*/
	  else if((isVar(symbuf,var_list,var_list_size))<0){/* If we did not meet var name */
	    if(!PARDIM){/*If we didnt find par name yet */
	      for(k=0;k<strlen(symbuf);k++)/* Copy par name */
		par_name[PARDIM][k]=symbuf[k];
	      fprintf(c_src_file,"P[%d]",PARDIM);/*transfer*/
	      PARDIM++;/* Increase number of par */
	    }
	    else if(PARDIM && (isPar(symbuf,PARDIM))>=0){/*If we've already found some
							   par names and newly
							   found name is par which
							   we met before*/
	      fprintf(c_src_file,"P[%d]",isPar(symbuf,PARDIM));
	    }
	    else if(PARDIM && (isPar(symbuf,PARDIM))<0) {/* Newly found name is really NEW par name */
	      for(k=0;k<strlen(symbuf);k++)
		par_name[PARDIM][k]=symbuf[k];
	      fprintf(c_src_file,"P[%d]",PARDIM);
	      PARDIM++;/* Increase number of par */
	    } else fprintf(stderr,"Error: I dont know what else is possible.\n");
	  }
	  else{fprintf(stderr,"Error: could not determine what kind of <right> symbol ");
	    fprintf(stderr,"`%s'", symbuf);break;}
	}
	if(isdigit(symbuf[0]) || symbuf[0]=='.'){/*Then we are dealing with number*/
	  /*Transfer*/
	  for(k=0;k<strlen(symbuf);k++)
	    fprintf(c_src_file,"%c",symbuf[k]);
	}
	if((!isalnum(symbuf[0]))){/*Error: neither char nor digit*/
	  fprintf(stderr,"Error: Neither char nor digit: ");
	  fprintf(stderr,"<right> simple.\n");
	}
      }
      j--;/*For correct passing through the odestr*/
      /*Transfer final `)'(closing parenthesis) of pow(_,_) function*/
      fprintf(c_src_file,")");
      /**************************************************************/
      continue;/*Yes, Let's continue after processing all '^'
		 related stuff*/
    }
    /*******************************************************************/
    /*************THE END OF '^' SPECIAL CASE INTERPRETATION************/
    /*******************************************************************/
    if(odestr[j]=='\n'){/*NEWLINE*/
      if(ENCLOSE)/*Enclosing your expression.*/
	fprintf(c_src_file,")");
      fprintf(c_src_file,";\n");/*Semicolon and Newline at the end of
				  one equation.*/
      break;/*We have read whole equation, then break*/
    }
    if(isalpha(odestr[j])){/* ALPHABETICAL CHAR */
      length=(size_t)40;/*Initially symbol is 40 char long*/
      symbuf=(char *)calloc(length,sizeof(char));/*Allocate
						   memory*/
      if(symbuf==NULL){/*Error handling*/
	fprintf(stderr,"Error: memory allocating for temporary buffer:	");
	fprintf(stderr,"<alphabetical char>.\n");
      }
      k=0;/*Looping over symbuf*/
      symbuf[k]=odestr[j];/* Get character to symbol buffer*/
      j++;/* Take another character */
      k++;/*Symbol buffer loop also increases*/
      while(isalnum(odestr[j])){/* Following can be alpha or digit*/
	if(k>length-1){
	  length+=40;
	  symbuf=(char *)realloc(symbuf,length*sizeof(char));
	  if(symbuf==NULL){/*Error handling*/
	    fprintf(stderr,"Error: memory reallocating for temporary buffer: ");
	    fprintf(stderr,"<alphabetical char>.\n");
	  }
	}
	symbuf[k]=odestr[j];/*Copy character */
	j++;k++;/*Take another character */
      }
      symbuf[k]='\0';
      /**************************FUNCTIONS****************************/
      if(function!=NULL && (isFunc(symbuf))>=0 && odestr[j]=='('){
	/*We have function with the same name and must be sure that
	  the next symbol is '(': as all functions should appear.*/
	k=j;/*Mark the position.*/
	length=40;/*string length.*/
	temp=(char *)malloc(length*sizeof(char));
	fprintf(c_src_file,"(func_%s(t,P",function[isFunc(symbuf)].func_name);
	j++;/*Step forward*/
	/*Store in temp arguments' list(without parentheses).*/
	while(odestr[j]!=')'){
	  /*This while loop also gets j at the correct position for the
	    next main loop iterations.*/
	  temp[j-k-1]=odestr[j];
	  if((j-k-1)==(length-1)){/*Reallocating in case of overbuffering.*/
	    length+=length;
	    temp=(char *)realloc(temp,length*sizeof(char));
	  }
	  j++;
	}
	temp[j-k-1]='\0';/*trailing null char.*/
	if(strlen(temp)==0){
	  /*If no arguments to the function.*/
	  fprintf(c_src_file,"))");
	  continue;
	}
	/*If there are some arguments to the function we proceed with
	  the comma "," delimiting first default arguments: t and P[].*/
	fprintf(c_src_file,",");
	/*We do not need symbuf: now it will hold the arguments one by
	  one of the function.*/
	length=40;
	symbuf=(char *)realloc(symbuf,length*sizeof(char));
	/*Get the list of arguments one by one.*/
	k=0;i=0;
	while(k<=strlen(temp)){
	  if(temp[k]==',' || k==strlen(temp)){
	    for(m=i;m<k;m++){/*Copying to symbuf.*/
	      symbuf[m-i]=temp[m];
	      if((m-i)==(length-1)){/*Reallocating*/
		length+=length;
		symbuf=(char *)realloc(symbuf,length*sizeof(char));
	      }
	    }
	    symbuf[m-i]='\0';
	    i=k+1;/*Hold new argument start position.*/
	    /*Normal processing of the argument.*/
	    process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	    if(temp[k]==',')
	      fprintf(c_src_file,",");/*Transfer delimiting comman ","*/
	  }
	  k++;
	}
	/*Transfer closing parentheses of the function.*/
	fprintf(c_src_file,"))");

	continue;
      }
      /**************************END FUNCTIONS***************************/
      j--;/*Decrease j as we need current character to be processed by
	    the parser in the next iteration of the main loop and
	    because now it is not processed.*/
      
      /* Compare symbuf we got with var names */
      if((isVar(symbuf,var_list,var_list_size))>=0){/* If we met var name */
	if(odestr[j+1]!='^')/*This is NOT special case,see above*/
	  fprintf(c_src_file,"%s",var_c_list[isVar(symbuf,var_list,var_list_size)]);
	continue;}
      else if((isVar(symbuf,var_list,var_list_size))<0){/* If we did not meet var name */
	if(!PARDIM){/*If we didnt find par name yet */
	  for(k=0;k<strlen(symbuf);k++)/* Copy par name */
	    par_name[PARDIM][k]=symbuf[k];
	  if(odestr[j+1]!='^')
	    fprintf(c_src_file,"P[%d]",PARDIM);/*transfer*/
	  PARDIM++;/* Increase number of par */
	  continue;/*No need to compare further*/
	}
	if(PARDIM && (isPar(symbuf,PARDIM))>=0){/*If we've already found some
					      par names and newly
					      found name is par which
					      we met before*/
	  if(odestr[j+1]!='^')
	    fprintf(c_src_file,"P[%d]",isPar(symbuf,PARDIM));
	  continue;/*No need to compare further*/
	}
	if(PARDIM && (isPar(symbuf,PARDIM))<0) {/* Newly found name is really NEW par name */
	  for(k=0;k<strlen(symbuf);k++)
	    par_name[PARDIM][k]=symbuf[k];
	  if(odestr[j+1]!='^')
	    fprintf(c_src_file,"P[%d]",PARDIM);
	  PARDIM++;/* Increase number of par */
	  continue;/*No need to compare further*/
	}
      }
      /*In principal, the next condition is never met, because any
      unrecognized symbol is treated as a parameter.*/
      else{fprintf(stderr,"Error: could not determine what kind of <alphabetical> symbol: ");
	fprintf(stderr,"`%s'\n", symbuf);break;}
      continue;/* No need to go for further checking */
    }
    if(isdigit(odestr[j]) || odestr[j]=='.'){/*NUMBERS*/
      k=j;/*Storing start position of the number*/
      j++;
      while(isdigit(odestr[j]) ||
	    odestr[j]=='.') j++;/*Jumping through number*/
      if(odestr[j]!='^'){/*Transfer if we are not dealing with special case '^'*/
	while(k<j)
	  fprintf(c_src_file,"%c",odestr[k++]);
      }
      j--;
      continue;
    }
    if(odestr[j]=='('){/*COMPLEX EXPRESSIONS in parenthesis `('
			    and `)'*/
      /*The only thing to do here is to search for closing `)' and
	determine whether after it there is an '^'---special
	case. If there is no '^' after closing `)' then just
	transfer open `(' and correctly return to the main loop.*/
      length=(size_t)40;
      symbuf=(char *)calloc(length,sizeof(char));
      if(symbuf==NULL){
	fprintf(stderr,"Error: memory allocating: ");
	fprintf(stderr,"<complex>.\n");
      }
      k=0;
      j++;/*Jump one symbol further*/
      m=1;/*Count for matched pairs `('/`);*/
      while((odestr[j]!=')' && m!=0) || (odestr[j]==')' && m!=0)){
	if(k==0 && odestr[j]=='(') m++;
	if(k==0 && odestr[j]==')') m--;
	symbuf[k]=odestr[j];
	if(k==length-1){
	  length+=40;
	  symbuf=(char *)realloc(symbuf,length*sizeof(char));
	  if(symbuf==NULL){
	    fprintf(stderr,"Error: memory reallocating: ");
	    fprintf(stderr,"<complex>.\n");
	  }
	}
	if(odestr[j]=='\n'){/*This means the end of ODE
			      specification string*/
	  break;
	}
	k++;j++;/*Further jumping*/
	if(odestr[j]=='(') m++;
	if(odestr[j]==')') m--;
      }
      symbuf[k]='\0';/*Terminating null character*/
      if(odestr[j]=='\n'){
	fprintf(stderr,"Error: there is no closing `)': ");/*Error case*/
	fprintf(stderr,"I searched till the end of the `%s'.\n",odestr);
	break;/*Another break for end of ODE specification to get out
		of the main loop.*/
      }
      j++;/*Yet another jump forward for '^' checking*/
      if(odestr[j]!='^'){/*If not special case*/
	/*Transfer open `('*/
	fprintf(c_src_file,"(");
	/***Processing complex expression***/
	process_odestr(0,c_src_file,symbuf,var_list,var_c_list,var_list_size);
	/***********************************/
	/*Transfer close `)'*/
	fprintf(c_src_file,")");
	j--;
      } else {j--;/*Will jump to '^' symbol at the next iteration of the
		    main loop.*/}
      continue;
    }
    /*Transfer any other characer*/
    fprintf(c_src_file,"%c",odestr[j]);/*Transfer character*/
  }

  return 0;
}

int get_parameter_name(char *symbuf)
{/* This function checks the symbuf if it is already existing name or the new
    parameter name. In both cases it returns the index of the parameter name,
    additionally, in the latter case it copies the new name to par_name. Otherwise,
    the function returns -1 (practically it is impossible). Note par_name and PARDIM
    are considered as global variables here. */
  int i;
  if((i = isPar(symbuf,PARDIM)) >= 0)
    /* We have found the existing name */
    return i;
  else if((i = isPar(symbuf,PARDIM)) < 0){
    /* We couldn't find name in data base OR there are no names at all. */
    if(PAR_NAME_SIZE > strlen(symbuf))
      i = strlen(symbuf);
    else
      i = PAR_NAME_SIZE;
    strncpy(par_name[PARDIM],symbuf,i);
    par_name[PARDIM][PAR_NAME_SIZE-1] = '\0';
    PARDIM++;/* Increase number of found parameters */
    return (PARDIM-1);/* return the index of new parameter */
  }
  else{/* Seems like an impossible option */
    fprintf(stderr,"Something undefined happened.\n");
    return -1;
  }
}

void process_update_vector(char *upstr,int *upvec)
{
  int k = 0;
  char *p = strchr(upstr, ',');
  char *p1;
  p1 = strchr(upstr,'\n');/* Find newline */
  *p1 = 0;/* Delete newline, we do not need it */
  while(p){
    /* Check the variable name */
    *p = 0;/* Nullify comma place in the upstr to get string manipulation */
    /* upstr contains exactly one position of the update vector */
    k = strtol(upstr,&p1,10);
    //printf("Got:k=%d for %s ==> ",k,p1);
    if(k == 0){
      if(*p1 == '-'){/* Single minus '-' mean -1 */
	k = -1;
	p1++;
      }
      else if(*p1 == '+'){/* Single minus '+' mean +1 */
	k = 1;
	p1++;
      }
      else
	k = 1;/* No integer means 1 */
    }
    //printf("changed: k=%d for %s\n",k,p1);
    if(isVar(p1,var_name,DIM)>-1)
      upvec[isVar(p1,var_name,DIM)] = k;
    else
      fprintf(stdout,"No var \"%s\" was found!Check your system.\n",p1);
    /* printf("%s",upstr); */
    /* printf("(%d->",k); */
    /* printf("%s),",p1); */
    upstr = p+1;
    p = strchr(upstr, ',');
  }
  k = strtol(upstr,&p1,10);
  //printf("Got:k=%d for %s ==> ",k,p1);
  if(k == 0){
    if(*p1 == '-'){/* Single minus '-' mean -1 */
      k = -1;
      p1++;
    }
    else if(*p1 == '+'){/* Single minus '+' mean +1 */
      k = 1;
      p1++;
    }
    else
      k = 1;/* No integer means 1 */
  }
  //printf("changed:k=%d for %s\n",k,p1);
  if(isVar(p1,var_name,DIM)>-1)
    upvec[isVar(p1,var_name,DIM)] = k;
  else
    fprintf(stdout,"No var \"%s\" was found!Check your system.\n",p1);
  /* printf("(%d->",k); */
  /* printf("%s)",p1); */
  /* printf("%s\n",upstr); */
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

int isFunc(char const *name)
{
  int i=0;
  while((strcmp(name,function[i].func_name))!=0){
    i++;
    if(i==FNUM){/*FNUM is global.*/
      return -1;/*We did not find function with the same name*/
    }
  }
  return i;/*We return function index as we met the condition.*/
}

int isSym(const char * buf,const char c)
{
  int i;
  for(i=0;i<strlen(buf);i++)
    if(buf[i]==c)/*`c' symbol*/
      return i;/*Return TRUE,index of the symbol*/

  return -1;/*No `c' symbol*/
}

int find_eqsign(char const * string)
{
  int j;
  for(j=0;j<strlen(string);j++)
    if(string[j]=='=')
      return j;/*return position of the '=' sign*/
  return -1;/* Equal sign not found*/
}

int how_many_eqsigns(char const * string)
{
  int j,i=0;
  for(j=0;j<strlen(string);j++)
    if(string[j]=='=')
      i++;
  return i;
}

int buf_strip(char *buffer, const int len)
{/* NOTE: ENTERS(NEWLINES) ARE NOT REMOVED */
  /* This function removes spaces, and tabs and substitutes them with '\0' */
  /* characters with merging the gaps occuring after the substitution */
  int i,j;
  //int len=strlen(buffer);
  for(i=0;i<len;i++){
    if(buffer[i]==' ' || buffer[i]=='\t'){
      buffer[i]='\0';
      for(j=i;j<len-1;j++)
	buffer[j]=buffer[j+1];
    }
  }
  return 0;
}

int print_heading(FILE *file)
{/*This functions prints some general information in the heading of
   the source .c file, like Copyright, Non-Warranty, Description etc*/
  if(file == NULL){
    fprintf(stderr,"Error: wrongly fopened file(print_heading).\n");
    return -1;
  }
  fprintf(file,"/*****************************************");
  fprintf(file,"*****************************************/\n");
  fprintf(file,"/* Copyright 2008,2009,2010,2011,2012,2013,2014 Elias Potapov. */\n");
  fprintf(file,"/* Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, \n");
  fprintf(file,"2003, 2004, 2005, 2006, 2007,");
  fprintf(file,"\n2008,2009,2010,2011\n The GSL Team. */\n");
  fprintf(file,"/*****************************************");
  fprintf(file,"*****************************************/\n");
  fprintf(file,"/* This file is automatically generated by DINAMICA\n");
  fprintf(file,"\tand is part of DINAMICA. */\n");
  fprintf(file,"/* DINAMICA is free software: you can redistribute ");
  fprintf(file,"it and/or modify */\n/* it under the terms of ");
  fprintf(file,"the GNU General Public License as published by */\n");
  fprintf(file,"/* the Free Software Foundation, either version 3 of the License, or */\n");
  fprintf(file,"/* (at your option) any later version. */\n\n");
  fprintf(file,"/* DINAMICA is distributed in the hope");
  fprintf(file," that it will be useful, */\n/* but WITHOUT ANY WARRANTY;");
  fprintf(file," without even the implied warranty of */\n");
  fprintf(file,"/* MERCHANTABILITY or FITNESS FOR A");
  fprintf(file," PARTICULAR PURPOSE.  See the */\n");
  fprintf(file,"/* GNU General Public License for more details. */\n");
  fprintf(file,"/* You should have received a copy");
  fprintf(file," of the GNU General Public License */\n");
  fprintf(file,"/* along with DINAMICA. ");
  fprintf(file," If not, see <http://www.gnu.org/licenses/>. */\n");
  fprintf(file,"/*****************************************");
  fprintf(file,"*****************************************/\n");
  fprintf(file,"/*****************************************");
  fprintf(file,"*****************************************/\n");
  fprintf(file,"/* Original author is Elias Potapov");
  fprintf(file," <elias.potapov@gmail.com>\n");
  fprintf(file,"   Lomonosov Moscow State University,");
  fprintf(file," Biophysics Dep..\n");
  fprintf(file,"   Lebedev Physical Inst.,");
  fprintf(file," Dep. of Theoretical Physics.\n");
  fprintf(file,"   Moscow, Russia */\n");
  fprintf(file,"/*****************************************");
  fprintf(file,"*****************************************/\n");

  return 0;
}

int print_tail(FILE *file)
{/*This functions prints some helpful information in the tail of
   the source .c file, like User parameters as they are seen in the
   program, user varibles etc. This would improve readability
   and help to report bugs.*/
  if(file == NULL){
    fprintf(stderr,"Error: wrongly fopened file(print_tail).\n");
    return -1;
  }
  fprintf(file,"/**************************************************/\n");
  fprintf(file,"/********IMPROVING READABILITY INFORMATION.*******\n");
  fprintf(file,"Variables: x[i], where i -- index; y[i] - discrete.\n");
  fprintf(file,"Parameters: P[i], where i -- index.\n");
  fprintf(file,"NOTE: i starts from 0(zero) in C, not from 1(one).\n");
  fprintf(file,"NOTE: parameters are numbered in the order\n");
  fprintf(file,"they appear in the ODE specification.\n");
  fprintf(file,"**************************************************\n");
  fprintf(file,"Please, verify your system. Report bugs.\n");
  fprintf(file,"It will help to improve the .ode -->.c converter.\n");
  fprintf(file,"**************************************************/\n");
  return 0;
}

int gen_conf(const char *basename)
{
  int i;
  FILE * file;
  char *filename;
  filename = app_ext(basename,"bcf",filename);
  fprintf(stdout,"Writing `%s'...",filename);
  file = fopen(filename,"w");
  if(file == NULL){
    fprintf(stderr,"Error: could not open file `%s'\n",filename);
    return BAD_FILE;
  }
  /**************************************************************/
  fwrite(&DIM,sizeof(int),(size_t)1,file);
  fwrite(&LDIM,sizeof(int),(size_t)1,file);
  fwrite(&PARDIM,sizeof(int),(size_t)1,file);
  fwrite(&FNUM,sizeof(int),(size_t)1,file);
  fwrite(&AUXNUM,sizeof(int),(size_t)1,file);
  fwrite(&NREAC,sizeof(int),(size_t)1,file);
  for(i=0;i<DIM;i++)
    fwrite(var_name[i],sizeof(char),VAR_NAME_SIZE,file);
  fwrite(xin_par,sizeof(int),(size_t)DIM,file);
  fwrite(xin,sizeof(double),(size_t)DIM,file);
  fwrite(yin,sizeof(long int),(size_t)DIM,file);
  for(i=0;i<PARDIM;i++)
    fwrite(par_name[i],sizeof(char),PAR_NAME_SIZE,file);
  fwrite(mu.P,sizeof(double),(size_t)PARDIM,file);
  int auxlen;/*Length of the aux name string, not permanent(!)*/
  for(i=0;i<AUXNUM;i++){
    auxlen=strlen(aux[i].name)+1;
    fwrite(&auxlen,sizeof(int),(size_t)1,file);
    fwrite(aux[i].name,sizeof(char),auxlen,file);
  }
  fwrite(&mynum.total_time,sizeof(double),(size_t)1,file);
  fwrite(&mynum.trans_time,sizeof(double),(size_t)1,file);
  fwrite(&mynum.step,sizeof(double),(size_t)1,file);
  fwrite(&mynum.write_step,sizeof(int),(size_t)1,file);
  fwrite(&mynum.smp_frq,sizeof(double),(size_t)1,file);
  /* if(method==NULL) */
  /*   method=(char *)calloc(METH_NAME_LEN,sizeof(char)); */
  fwrite(method,sizeof(char),(size_t)METH_NAME_LEN,file);
  /* if(method2==NULL) */
  /*   method2=(char *)calloc(METH_NAME_LEN,sizeof(char)); */
  fwrite(method2,sizeof(char),(size_t)METH_NAME_LEN,file);
  fwrite(&traj_trans,sizeof(int short),(size_t)1,file);
  fwrite(&per_method,sizeof(int short),(size_t)1,file);
  fwrite(&per_var,sizeof(int short),(size_t)1,file);
  fwrite(&per_thresh,sizeof(double),(size_t)1,file);
  fwrite(&cross_level[0],sizeof(double),(size_t)1,file);
  fwrite(&ma_span,sizeof(int),(size_t)1,file);
  fwrite(&BUFFER,sizeof(int),(size_t)1,file);
  fwrite(&BUF_INCR,sizeof(int),(size_t)1,file);
  fwrite(&write_flag,sizeof(int short),(size_t)1,file);
  fwrite(&jac_flag,sizeof(int short),(size_t)1,file);
  fwrite(&lang_flag,sizeof(int short),(size_t)1,file);
  fwrite(&graph_flag,sizeof(int short),(size_t)1,file);
  fwrite(&eps_abs_tper,sizeof(double),(size_t)1,file);
  fwrite(&eps_rel_tper,sizeof(double),(size_t)1,file);
  fwrite(&eps_abs_int,sizeof(double),(size_t)1,file);
  fwrite(&eps_rel_int,sizeof(double),(size_t)1,file);
  fwrite(&a_y,sizeof(double),(size_t)1,file);
  fwrite(&a_dydt,sizeof(double),(size_t)1,file);
  fwrite(&eps_abs_per,sizeof(double),(size_t)1,file);
  fwrite(&eps_rel_per,sizeof(double),(size_t)1,file);
  fwrite(&eps_abs_peak,sizeof(double),(size_t)1,file);
  fwrite(&eps_rel_peak,sizeof(double),(size_t)1,file);
  fwrite(&eps_abs_am,sizeof(double),(size_t)1,file);
  fwrite(&eps_rel_am,sizeof(double),(size_t)1,file);
  fwrite(&eps_per_ratio,sizeof(double),(size_t)1,file);
  fwrite(&eps_inregime_ratio,sizeof(double),(size_t)1,file);
  fwrite(&graph,sizeof(struct coordNet),(size_t)1,file);
  
  fclose(file);
  fprintf(stdout,"Done.\n");
  return 0;
}

int init()
{/*Here we set initial value for global working variables. THE
   DEFAULTS.*/
  int i;
  LDIM = 1;DIM = 0;PARDIM = 0;FNUM = 0;AUXNUM = 0;NREAC = 0;
  NLANG = 0;
  mynum.total_time = 20;
  mynum.trans_time = 50;
  mynum.step = 0.02;
  mynum.write_step = 1;
  mynum.smp_frq = 1.0;
  traj_trans = 1;
  //method = (char *)malloc((strlen("rkf45")+1)*sizeof(char));
  strcpy(method,"rkf45");
  per_method = 0;//Poincare sections
  per_var = -1;//Automatic detection of var for periods
  per_thresh = 0;//Period threshold
  cross_level[0] = 4;//Cross level
  cross_level[1] = 0;
  cross_level[2] = 0;
  ma_span = 0;//MA span
  BUFFER = 10000;
  BUF_INCR = 5000;
  mynum.global_buffer = BUFFER;
  write_flag = 1;
  /*Errors*/
  eps_abs_tper = 1e-4;
  eps_rel_tper = 0.1;
  eps_abs_int = 1e-4;
  eps_rel_int = 0.0;
  a_y = 1.0;
  a_dydt = 0.0;
  eps_abs_per = 1e-4;
  eps_rel_per = 0.02;
  eps_abs_peak = 1e-4;
  eps_rel_peak = 0.02;
  eps_abs_am = 1e-3;
  eps_rel_am = 0.01;
  eps_per_ratio = 1e-5;
  eps_inregime_ratio = 1e-5;
  /* Reading source inits */
  var_name = NULL;
  odestr = NULL;
  initstr = NULL;
  jacstr = NULL;
  tjacstr = NULL;
  propstr = NULL;
  upvec = NULL;
  cfname = NULL;
  /* Graph system */
  if(HAVE_GNUPLOT)
    graph_flag = 1;
  else
    graph_flag = 0;
  graph.xInd = -1;
  graph.yInd[0] = 0;
  graph.yInd[1] = -2;
  graph.yInd[2] = -2;
  
  return 0;
}

int process_par(char * parstr)
{/*ODE file format for parameters' values example: "param Q=0.1,b=70,alpha=30".
   Instead of `param' keyword there may be `par' or just `p'.*/
  int i,j;
  char *p1,*p2;

  if(strncmp(parstr,"param",strlen("param"))==0)
    i = strlen("param");
  else if(strncmp(parstr,"par",strlen("par"))==0)
    i = strlen("par");
  else if(strncmp(parstr,"p",strlen("p"))==0)
    i = strlen("p");
  else
    fprintf(stderr,"Error: this is not parameter line.\n");

  p1 = parstr + i;
  while((p2 = strchr(p1,'=')) != NULL){
    *p2 = 0;/* terminate here, so p1 contains par name */
    i = get_parameter_name(p1);
    //printf("Got par=%s, index=%d\t",p1,i);
    p1 = p2 + 1;/* moving one position forward */
    if((p2 = strchr(p1,',')) != NULL)/* next comma */
      *p2 = 0;/* terminate here, so p1 contains par value */
    else if((p2 = strchr(p1,'\n')) != NULL)/* the last entry before newline */
      *p2 = 0;/* terminate here, so p1 contains par value */
    /* Setting the parameter value */
    if(!isalpha(*p1))
      mu.P[i] = atof(p1);
    else
      printf("Error: par %s has wrong value!\n",par_name[i]);
    //printf("mu.P[%d] = %G\n",i,mu.P[i]);
    p1 = p2 + 1;/* moving one position forward */
  }
  return 0;
}

int process_init(char *initstr)
{/*ODE file format for initial values example: "init a1=0.342,a2=7.23403,b1=3.0".
   Instead of `init' keyword there may be `i'.  This function should be invoked after
   all equations are defined, namely, when the var name space is known. Note also,
   that this function accepts initstr without any spaces and tabs between symbols.
 */
  /*Dont forget to allocate memory for *xin and *yin before calling
    this function*/
  int i,j;
  char *p1,*p2;

  if(strncmp(initstr,"init",strlen("init"))==0)
    i = strlen("init");
  else if(strncmp(initstr,"i",strlen("i"))==0)
    i = strlen("i");
  else fprintf(stderr,"Error: this is not initial values line.\n");

  p1 = (initstr + i);
  while((p2 = strchr(p1,'=')) != NULL){/* We found next '=' sign */
    *p2 = 0;/* Terminate at '=' sign, so p1 is var name */
    //printf("Read var=%s\n",p1);
    j = isVar(p1,var_name,DIM);/* j contains index of var or error value */
    if(j < 0)/* Print error message */
      fprintf(stderr,"Init: unknown var=%s\n",p1);
    p1 = p2 + 1;/* Moving one position forward */
    if((p2 = strchr(p1,',')) != NULL)/* We find next ',' sign */
      *p2 = 0;/* terminate here, so p1 contains value/parameter */
    else if ((p2 = strchr(p1,'\n')) != NULL)/* we find the end of line */
      *p2 = 0;/* terminate here, so p1 contains value/parameter */
    //printf("Read value/par=%s\n",p1);
    if(isalpha(*p1)){/* We have non-numeric value on the left hand side */
      i = get_parameter_name(p1);
      if(j >= 0){
	xin[j] = mu.P[i];
	xin_par[j] = i;/* Save index of the parameter */
	//printf("xin[%d] = %G.\t",j,xin[j]);
	yin[j] = (long)mu.P[i];
	//printf("yin[%d] = %G\n",j,yin[j]);
      }
    }
    else{/* Otherwise it is numeric and converted to real value */
      if(j >= 0){/* We have the variable */
	xin[j] = atof(p1);
	//printf("xin[%d] = %G.\t",j,xin[j]);
	yin[j] = atol(p1);
	//printf("yin[%d] = %G\n",j,yin[j]);
      }
    }
    p1 = p2 + 1;
  }
  return 0;
}

int process_inter_par(const char * iparstr)
{/*ODE file format for Dinamica internal parameters' values example:
   "@ total=100,dt=0.2".*/
  int i,j,k,length = 20;
  char *symbuf;
  char *value;
  symbuf = (char *)calloc(length,sizeof(char));
  value = (char *)calloc(length,sizeof(char));
  i = 1;/*As first symbol is `@' in the string*/
  j = 0;
  while(iparstr[i]!='\n'){
    while(iparstr[i]!='='){
      if(j==length){
	length += 20;
	symbuf = (char *)realloc(symbuf,length*sizeof(char));
	value = (char *)realloc(value,length*sizeof(char));
      }
      symbuf[j] = iparstr[i];
      i++;
      j++;
    }
    symbuf[j] = '\0';
    i++;
    j = 0;/*Terminating null char;moving on from '=' sign; nulling counter j*/
    while(iparstr[i]!=',' && iparstr[i]!='\n'){
      if(j==length){
	length += 20;
	symbuf = (char *)realloc(symbuf,length*sizeof(char));
	value = (char *)realloc(value,length*sizeof(char));
      }
      value[j] = iparstr[i];
      i++;
      j++;
    }
    value[j] = '\0';
    if(iparstr[i]!='\n')
      i++;
    j = 0;
    /*Next we distinguish internal parameters and set the options*/
    if(strcmp(symbuf,"total")==0)
      mynum.total_time = atof(value);
    else if(strcmp(symbuf,"dt")==0)
      mynum.step = atof(value);
    else if(strcmp(symbuf,"trans")==0)
      mynum.trans_time = atof(value);
    else if(strcmp(symbuf,"method")==0 || strcmp(symbuf,"m")==0){
      //method = (char *)calloc(METH_NAME_LEN,sizeof(char));
      if(strlen(value)>=METH_NAME_LEN){
	fprintf(stderr,"@Error: meth name must not exceed %d char in length\n",METH_NAME_LEN);
	return 45;
      }
      strcpy(method,value);
    }
    else if(strcmp(symbuf,"method2")==0 || strcmp(symbuf,"m2")==0){
      //method2 = (char *)calloc(METH_NAME_LEN,sizeof(char));
      if(strlen(value)>=METH_NAME_LEN){
	fprintf(stderr,"@Error: method name must not exceed %d char in length\n",METH_NAME_LEN);
	return 45;
      }
      strcpy(method2,value);
    }
    else if(strcmp(symbuf,"buffer")==0)
      BUFFER = atof(value);
    else if(strcmp(symbuf,"bufinc")==0)
      BUF_INCR = atof(value);
    else if(strcmp(symbuf,"ws")==0)
      mynum.write_step = atoi(value);
    else if(strcmp(symbuf,"sf")==0)
      mynum.smp_frq = atof(value);
    else if(strcmp(symbuf,"wf")==0)
      write_flag = atoi(value);
    else if(strcmp(symbuf,"lf")==0 || strcmp(symbuf,"lang")==0)
      lang_flag = atoi(value);
    else if(strcmp(symbuf,"gf")==0 || strcmp(symbuf,"graph")==0){
      if(HAVE_GNUPLOT)
	graph_flag = atoi(value);
      else
	fprintf(stderr,"Warning: Gnuplot disabled during Dinamica compilation.\n");
    }
    else if(strcmp(symbuf,"xaxis")==0 || strcmp(symbuf,"xax")==0){
      if(isalpha(*value))
	graph.xInd = isVar(value,var_name,DIM);
      else
	graph.xInd = atoi(value) - 1;
      if(graph.xInd < 0 || graph.xInd > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: xaxis\n");
	graph.xInd = -1;
      }
    }
    else if(strcmp(symbuf,"yaxis")==0 || strcmp(symbuf,"yax")==0){
      if(isalpha(*value))
	graph.yInd[0] = isVar(value,var_name,DIM);
      else
	graph.yInd[0] = atoi(value) - 1;
      if(graph.yInd[0] < 0 || graph.yInd[0] > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: yaxis\n");
	graph.yInd[0] = 0;
      }
    }
    else if(strcmp(symbuf,"yaxis1")==0 || strcmp(symbuf,"yax1")==0){
      if(isalpha(*value))
	graph.yInd[0] = isVar(value,var_name,DIM);
      else
	graph.yInd[0] = atoi(value) - 1;
      if(graph.yInd[0] < 0 || graph.yInd[0] > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: yaxis1\n");
	graph.yInd[0] = 0;
      }
    }
    else if(strcmp(symbuf,"yaxis2")==0 || strcmp(symbuf,"yax2")==0){
      if(isalpha(*value))
	graph.yInd[1] = isVar(value,var_name,DIM);
      else
	graph.yInd[1] = atoi(value) - 1;
      if(graph.yInd[1] < 0 || graph.yInd[1] > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: yaxis2\n");
	graph.yInd[1] = -2;
      }
    }
    else if(strcmp(symbuf,"yaxis3")==0 || strcmp(symbuf,"yax3")==0){
      if(isalpha(*value))
	graph.yInd[2] = isVar(value,var_name,DIM);
      else
	graph.yInd[2] = atoi(value) - 1;
      if(graph.yInd[2] < 0 || graph.yInd[2] > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: yaxis3\n");
	graph.yInd[2] = -2;
      }
    }
    else if(strcmp(symbuf,"pm")==0 || strcmp(symbuf,"permeth")==0){
      if((atoi(value)) == 1 || (atoi(value)) == 0)
	per_method = atoi(value);
      else
	fprintf(stderr,"@Error: wrong per_method. Left default.\n");
    }
    else if(strcmp(symbuf,"pv")==0 || strcmp(symbuf,"pervar")==0){
      if(isalpha(*value))
	per_var = isVar(value,var_name,DIM);
      else
	per_var = atoi(value) - 1;
      if(per_var < -1 || per_var > (DIM-1)){
	fprintf(stderr,"@Error: wrong input: per_var. Left default.\n");
	per_var = -1;
      }
    }
    else if(strcmp(symbuf,"c")==0 || strcmp(symbuf,"cross")==0){
      if(atof(value) < 0)
	fprintf(stderr,"@Error: cross_value less than 0. Left default\n");
      else
	cross_level[0] = atof(value);
    }
    else{
      fprintf(stderr,"@Error: could not recognize the option: %s=%s\n",symbuf,value);
    }
  }
  return 0;
}
