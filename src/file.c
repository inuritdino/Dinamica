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
/* This file contains routines for operations on data files comprising */
/* time series (usually) like converting to other file formats, split */
/* data sets into different files etc */
#include "init.h"
#include <string.h>
#define MAX_N_DIG 20 /* min 20 */

double ** load_frame(FILE *data, double **frame, int *n_pts)
{
  /*Loads the frame containing the timeseries from data. The 'data' must be
    initialized elsewhere, because this function load only segment from the file,
    separated from other segments by double '\n' indicator. This indicator is the
    convention in DINAMICA. This function should be used inside caller function which
    load the whole data file in some kind of for/while loop. It returns 2d array of
    doubles(frame) and number of lines/points read(n_pts).*/
  int i,k,end_frame_count = 0;
  int realloc_count = 1;
  char *symbuf,*start_symbuf,*tmp;
  char delim = ' ';
  *n_pts = 0;
  frame = (double **)malloc((DIM+1)*sizeof(double *));
  if(frame == NULL){
    fprintf(stderr,"Error: could not malloc the frame.\n");
    return NULL;
  }
  for(i=0;i<(DIM+1);i++){
    frame[i] = (double *)malloc(FRAME_N_LINES*sizeof(double));
    if(frame[i] == NULL){
      fprintf(stderr,"Error: could not malloc the frame[%d].\n",i);
      return NULL;
    }
  }
  symbuf = (char *)malloc(FRAME_LINE_LEN*sizeof(char));
  start_symbuf = symbuf;/* Store the position of symbuf */
  /* Start reading line by line */
  while((symbuf=fgets(symbuf,FRAME_LINE_LEN,data)) != NULL){
    if((*symbuf) == '#')/* Technical line in the beginning of the data */
      continue;
    if((*n_pts) == (realloc_count*FRAME_N_LINES)){
      /* We need reallocation of memory */
      realloc_count++;
      for(i=0;i<(DIM+1);i++){
	frame[i] = (double *)realloc(frame[i],
				     realloc_count*FRAME_N_LINES*sizeof(double));
	if(frame[i] == NULL){
	  fprintf(stderr,"Error: could not realloc the frame[%d].\n",i);
	  return NULL;
	}
      }
    }
    k = 0;/* Number of values read in a line */
    tmp = start_symbuf;/* Difference between symbuf and tmp determines the value */
    if(strchr(symbuf,'\n') == NULL){/* No newline means too many symbols in the line*/
      fprintf(stdout,"Warning: you have more than %d bytes in the line.\n",FRAME_LINE_LEN);
      fprintf(stdout,"Exit.\n");
      break;
    }
    *(strchr(symbuf,'\n')) = 0;/* Remove newline */
    if(*symbuf == 0){
      end_frame_count++;/* one '\n' in the symbuf */
      if(end_frame_count == 2){
	/* End of a frame, i.e. 2 '\n's occured in a row, this is assured by restoring
	   the value of end_frame_count to 0, after this if statement(see below) */
	free(start_symbuf);/* Free the symbuf original line */
	return frame;/* Return success */
      }
      continue;/* No need to go on down the code */
    }
    end_frame_count = 0;/* Restore the value of end_frame_count */
    /* Analyze the read line */
    while(1){
      symbuf = strchr(symbuf,delim);/* Find next delimiter, end of the value */
      /* Even for the very last value in the string we have the delimiter following
	 that value (see integrator.c). This should NOT be changed then. */
      if(symbuf == NULL)/* No delimiter is found, we break */
	break;
      *symbuf = 0;/* Remove the delimiter and put "end of the string" */
      frame[k++][*n_pts] = atof(tmp);/* Get the value into the frame */
      tmp = ++symbuf;/* Put the start of the next value */
    }
    symbuf = start_symbuf;/* Restore the start position of symbuf */
    (*n_pts)++;/* Increase number of points(=lines) read */
  }
  fprintf(stdout,"Error: End of file reached. ");
  fprintf(stdout,"Frame was not terminated.\n");
  return NULL;/* not success */
}

void free_frame(double **frame)
{
  int i;
  for(i=0;i<DIM+1;i++)
    free(frame[i]);
  free(frame);
}

void write_frame(FILE *data, double **frame,const int frame_size)
{
  /* Writes the frame to the output file, the descriptor for which must be allocated
     elsewhere. frame_size is the number of points in the frame, actual points
     containing real values(returned by load_frame()). The number of values
     for each of the point is DIM+1. */
  int i,j;
  if(data == NULL)
    fprintf(stderr,"Error: file descriptor is bad\n");
  else{
    for(i=0;i<frame_size;i++){
      for(j=0;j<DIM+1;j++)
	fprintf(data,"%.5lf ",*(*(frame+j)+i));//frame[j][i]);
      fprintf(data,"\n");
    }
    fprintf(data,"\n\n");/* Separating consequent frames in the file */
  }
}

int read_out_file(const char *fname)
{
  int indicator = 0;
  char *line = (char *)malloc((MAX_ARG_LEN+1)*sizeof(char));
  FILE *fid = fopen(fname,"r");
  if(fid == NULL){
    fprintf(stderr,"Error(read_out_file): could not open file.\n");
    return 10;
  }
  while( (fgets(line,MAX_ARG_LEN,fid)) != NULL){
    if(strlen(line) == 1){/* Only one char, which must be \n */
      /* Skip these lines */
      continue;
    }
    if(indicator == 1){
      nPerStoch++;
      /* Must it be here MALLOC first if perStoch == NULL ??? Tricky place. */
      perStoch = (double *)realloc(perStoch,nPerStoch*sizeof(double));
      perStoch[nPerStoch-1] = atof(line);
    }
    if(strncmp(line,"#PerStoch(",strlen("#PerStoch(")) == 0){
      indicator = 1;
      nPerStoch = 0;
    }
  }
  free(line);
  fclose(fid);
  return 0;
}

int set_info_data(int *info,FILE *output)
{
  /* See the explanation about the info array in get_info_data() function below */
  int i = 0;
  fprintf(output,"#");
  /* The info is 5 numbers, the last one is from ma_span */
  while(i < 4){
    fprintf(output,"%d,",info[i]);
    i++;
  }
  if(ma_span)
    fprintf(output,"%d\n",ma_span);
  else
    fprintf(output,"\n");

  return 0;
}

int *get_info_data(int *info, FILE *input)
{/* The function reads the input from *input and gets the technical information of
    the timeseries stored in the source, like number of simulation runs(nRuns) stored
    and if there was a complex run(complex) or it was a pure deterministic
    traj(pure_det). input must be opend from elsewhere. It returns the array of
    integers, position of each flag/number is pre-defined. */

  char *tmp = malloc((DIM+1)*MAX_N_DIG*sizeof(char));
  char *p = tmp;
  int col = 0;
  /* The info is 4 numbers + ma_span if exists. NOTE: ma_span should not be in the
  original, not MA, data. */
  info = (int *)malloc(4*sizeof(int));
  /* Check if input is opened. */
  if(input == NULL){
    fprintf(stderr,"Error(get_info): perhaps no file was opened.\n");
    return NULL;
  }
  else{/* input is O.K. */
    fgets(tmp,(DIM+1)*MAX_N_DIG-1,input);/* Take the line */
    if((*tmp) != '#'){
      fprintf(stderr,"Error: no configuration line, old format data.\n");
      return NULL;
    }
    else if((strchr(tmp,'\n')) == NULL){
      fprintf(stderr,"Error: I could not read the whole line.\n");
      return NULL;
    }
    else{
      tmp++;/* Move away from the first '#' */
      while((*p) != '\n'){/* Loop until first newline */
	if((*p) == ','){
	  (*p) = 0;
	  info[col++] = atoi(tmp);
	  tmp = p + 1;
	}
	p++;
      }
      if((*p) == '\n'){
	(*p) = 0;
	info[col] = atoi(tmp);
      }
    }
  }
  return info;
}

int convert_data_file(char const *filename, int short ff)
{/* filename -- data file name; ff is file format to convert to. ff can */
  /* be: ff=1 -> csv(comma separated data file*/
  /* ff=2 -> TBD*/
  /* NOTE: default DINAMICA format is space separated data file and 2
     newlines (\n) separated data sets. */
  /* ndata -- number of data sets to read, is to be determined from the data source
     information line. The numbering goes from <filename>-I<.csv>, where I goes from
     0 to (ndata-1)*/
  /* The output is written to hard drive. If method is complex, the deterministic
     output is not numbered, i.e. <filename>.csv. */
  int i,outfn_count = 0,newline_count = 0;
  int *info;
  int ndata;
  char *ext;/* extension of the output */
  FILE *in,*out;
  in = fopen(filename,"r");
  if(in==NULL){
    fprintf(stderr,"Error: no file with name `%s'\n",filename);
    return 100;
  }
  if((info=get_info_data(info,in)) == NULL){
    fprintf(stderr,"Error: could not read the data information.\n");
    return 199;
  }
  else{
    /* complex flag + nRuns is the total number of data sets in the source */
    ndata = info[0] + info[1];
  }
  /*** Separator and extension strings ***/
  char *sep;
  if(ff==1){
    sep = ",";
    ext = ".csv";
  }
  else {
    fprintf(stderr,"Error: dont know the format\n");
    return 10;
  }
  /* Allocating temporary and filename strings */
  /* DIM+1 variables written as "%.5lf" meaning at min 7 digits/symbols, we put
     MAX_N_DIG that is to be 20 at min for certainty */
  char *tmp = (char *)malloc((DIM+1)*MAX_N_DIG*sizeof(char));
  /* Output filename, e.g. original name appended with ".csv"*/
  char *outfn = (char *)malloc((strlen(filename)+10)*sizeof(char));
  if(ndata == 1){
    if(ff==1){
      strcpy(outfn,filename);
      strcat(outfn,ext);
    }
    out = fopen(outfn,"w");
    printf("Converting from `%s' to `%s'\n",filename,outfn);
    /* Printing variables to the header of the output */
    fprintf(out,"Time,");
    for(i=0;i<DIM;i++){
      fprintf(out,"%s,",var_name[i]);
    }
    fprintf(out,"\n");
    /* Start actual transfer */
    while((fgets(tmp,(DIM+1)*MAX_N_DIG,in))!=NULL){
      if((strchr(tmp,'\n')) == NULL){
	fprintf(stderr,"Error: could not get whole line.\n");
	return 101;
      }
      for(i=0;i<strlen(tmp);i++){
	if(tmp[i]==' ')/* Default space separator changed to `sep' */
	  fprintf(out,"%s",sep);
	else
	  fprintf(out,"%c",tmp[i]);
      }
    }
    fclose(out);
  }
  else{/* More than 1 data sets in the source: stoch or complex */
    printf("Converting from `%s' to ",filename);
    outfn = (char *)realloc(outfn,
			    (strlen(filename)+20)*sizeof(char));
    while(outfn_count < ndata){
      /* FORMING OUTPUT FILENAMES */
      if((outfn_count==0) && (info[0])){
	/* Determ data set */
	strcpy(outfn,filename);
	strcat(outfn,ext);
      }
      else{
	/* Stoch sets */
	strcpy(outfn,filename);
	sprintf(tmp,"-%d",(outfn_count-info[0]));
	strcat(outfn,tmp);
	strcat(outfn,ext);
      }
      /* ******************* */
      out = fopen(outfn,"w");
      printf("`%s'\n",outfn);
      /* Printing variables to the header of the output */
      fprintf(out,"Time,");
      for(i=0;i<DIM;i++){
	fprintf(out,"%s,",var_name[i]);
      }
      fprintf(out,"\n");
      while((fgets(tmp,(DIM+1)*MAX_N_DIG,in))!=NULL){
	if((strchr(tmp,'\n')) == NULL){
	  fprintf(stderr,"Error: could not get whole line.\n");
	  return 101;
	}
	if((strlen(tmp)==1) || (tmp[0]=='#')){
	  /* Just newline or comment*/
	  newline_count++;
	  continue;/* Skip the rest */
	}
	if(newline_count==2){
	  /* Final condition for a data set, 2 newlines must be in
	     a row in the data file */
	  newline_count = 0;
	  fclose(out);
	  break;
	}
	newline_count = 0;/* Only 2 newlines in a row count */
	for(i=0;i<strlen(tmp);i++){
	  if(tmp[i]==' ')
	    fprintf(out,"%s",sep);
	  else
	    fprintf(out,"%c",tmp[i]);
	}
      }
      outfn_count++;
    }
  }
  fclose(in);
  free(tmp);free(outfn);
  
  return 0;
}
