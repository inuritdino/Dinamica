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
   Moscow, Russia / Tampere, Finland */
/****************************************************************************************/
/* input_interpreter.c is the file for interpreting any command of DINAMICA typed in
   menu or any submenu. The main function here is read_menu(...) that reads any input
   from user in any menu/submenu. Other functions here are of the form
   MENU_interp(...), where MENU is appropriate menu of the program.  */

#include "init.h"
#include "continue.h"
#include "errors.h"
#include "random.h"
#include "singularity.h"
#include "lyapunov.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <gsl/gsl_statistics_double.h>
#define MAX_SYM_READ 100

/*This file describes interpretation of the commands*/
void copy_history()
{
  int i,j;
  /* Copy the commands if it is not a history command "!" */
  if(strcmp(cmd[0],"!") && strlen(cmd[0])){
    for(i=MAX_N_HIST_ENT-1;i>0;i--){
      if(strlen(buf_hist[i-1][0])){
	for(j=0;j<MAX_N_ARG;j++)
	  strcpy(buf_hist[i][j],buf_hist[i-1][j]);
      }
    }
    for(i=0;i<MAX_N_ARG;i++)
      strcpy(buf_hist[0][i],cmd[i]);
  }
}

int parse_command_line(char **cmd,const char *prompt,FILE *stream)
{/* Read the command line and separate all that was read into the words-commands. The
    array of commands is returned. Reading is from stream or stdin(when
    stream==NULL). */
  int narg = 0,i;
  char *symbuf,*tofree,*token,*p1,*p2;
  /* Print the prompt first */
  fprintf(stdout,"%s>",prompt);
  /* Free the memory of the old calls */
  for(i=0;i<MAX_N_ARG;i++)
    memset(cmd[i],'\0',MAX_ARG_LEN);
  /* Read the line */
  symbuf = read_command_line(stream,&i);
  if(i == 2)/* EOF reached */
    return 100;
  if(i == 3)/* Empty line */
    return 200;
  
  if(strchr(symbuf,'#') != 0){/* Ignore the comments */
    *(strchr(symbuf,'#')) = 0;
  }
  tofree = symbuf;/* for subsequent freeing of the memory. */
  if(strchr(symbuf,'(') == NULL){
    while((token = strsep(&symbuf," \t")) != NULL){
      if(strlen(token) != 0){
	strncpy(cmd[narg++],token,MAX_ARG_LEN-1);
      }
    }
  }
  else{/* Group of args is found */
    p1 = strchr(symbuf,'(');
    p2 = strchr(p1,')');/* must not be NULL */
    *p1 = 0;
    *p2 = 0;
    while((token = strsep(&symbuf," \t")) != NULL){
      /* Copy the part before opening ( */
      if(strlen(token) != 0){
	strncpy(cmd[narg++],token,MAX_ARG_LEN-1);
      }
    }
    /* Copy the group argument */
    strncpy(cmd[narg++],p1+1,MAX_ARG_LEN-1);
    //symbuf = p2 + 1;
    p2 = p2+1;
    while((token = strsep(&p2," \t")) != NULL){
      if(strlen(token) != 0){
      	strncpy(cmd[narg++],token,MAX_ARG_LEN-1);
      }
    }
  }
  free(tofree);
  /* History account */
  copy_history();
  return 0;
}

int read_menu(char **cmd,const char *prompt)
{
  int i,k;
  char c;
  /* Free the memory of the old calls */
  for(i=0;i<MAX_N_ARG;i++)
    memset(cmd[i],'\0',MAX_ARG_LEN);
  fprintf(stdout,"%s>",prompt);
  k = 0;/*Number of commands read*/
  while((c=getchar()) == ' ' || c == '\t');/* Ignore first spaces and tabs */
  i = 0;
  while(c!='\n'){
    if(c==' ' || c=='\t'){/*Separate words*/
      cmd[k][i++] = '\0';
      k++;
      i = 0;
      while((c=getchar()) == ' ' || c == '\t');/* Ignore additional spaces and
      tabs before end '\n' of input */
      if((k == MAX_N_ARG) && (c != '\n')){
	fprintf(stderr,"Error: max # argin exceeded. Use only max # argin=%d\n",MAX_N_ARG);
	break;
      }
      continue;
    }
    cmd[k][i++] = c;
    c = getchar();
    if(c == '\n'){/* Commands DO NOT contain trailing newlines */
        cmd[k][i] = '\0';/* Trailing zero */
	k++;/* Increase number of read arguments */
    }
  }
  /* Copy the command to the history */
  copy_history();
  return 0;
}

int main_interp(char **buffer,const int depth_level)
{
  int i,j,k,*info;
  FILE *fname;
  int Len = 50;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  if(strcasecmp(buffer[depth_level],"exit") == 0 || \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Quiting normally...\n");
    free(symbuf);
    return 1000;
  }
  else if(strcmp(buffer[depth_level],"!") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      printf("History is:\n");
      for(j=MAX_N_HIST_ENT;j>0;j--){
	if(strlen(buf_hist[j-1][0])){
	  printf("%d) ",j);
	  for(i=0;i<MAX_N_ARG;i++)
	    printf("%s ",buf_hist[j-1][i]);
	  printf("\n");
	}
      }
    }
    else{
      if(!isalpha(symbuf[0])){
	main_interp(buf_hist[atoi(symbuf)-1],0);
      }
    }
  }
  else if(strcasecmp(buffer[depth_level],"ls") == 0) {
    fprintf(stdout,"MAIN menu:\n");
    fprintf(stdout,"(R)un*\n");
    fprintf(stdout,"(R)un (t)ransient*\n");
    fprintf(stdout,"(R)un (i)nitial*\n");
    fprintf(stdout,"(C)alculate*\n");
    fprintf(stdout,"(F)ile/\n");
    fprintf(stdout,"(N)umerics/\n");
    fprintf(stdout,"(P)arameters*\n");
    fprintf(stdout,"(V)ariables*\n");
    fprintf(stdout,"P(e)riodics/\n");
    fprintf(stdout,"(G)raphics/\n");
    fprintf(stdout,"(I)nitials*\n");
    fprintf(stdout,"(T)rajectory/\n");
    fprintf(stdout,"C(o)ntinue/\n");
    fprintf(stdout,"R(a)ndom/\n");
    fprintf(stdout,"(Er)rors/\n");
    fprintf(stdout,"(S)ingularity/\n");
    fprintf(stdout,"(L)yapunov/\n");
    fprintf(stdout,"(R)un (l)inear*\n");
  }
  else if((strcasecmp(buffer[depth_level],"v") == 0) ||		\
	  (strcasecmp(buffer[depth_level],"var") == 0) ||	\
	  (strcasecmp(buffer[depth_level],"variables") == 0)){
    fprintf(stdout,"Variables:\n");
    j = 2;
    for(i=0;i<DIM;i++){
      fprintf(stdout,"U(%d) = %s\t",i+1,var_name[i]);
      if(i == j){
	printf("\n");
	j += 3;
      }
    }
    printf("\n");
  }
  else if((strcasecmp(buffer[depth_level],"graphics") == 0) ||	\
	  (strcasecmp(buffer[depth_level],"g") == 0)){
    if(strlen(buffer[depth_level+1])!=0)
      graphics_interp(buffer,depth_level+1);
    else{
      while(1){
	//read_menu(buffer,graphics_prompt);
	i = parse_command_line(cmd,graphics_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	if((graphics_interp(buffer,0)) == 1000)
	  break; 
      }
    }

  }
  else if((strcasecmp(buffer[depth_level],"numerics") == 0) ||	\
	  (strcasecmp(buffer[depth_level],"n") == 0)){
    if(strlen(buffer[depth_level+1])!=0)
      numerics_interp(buffer,depth_level+1);
    else{
      while(1){
	//read_menu(buffer,numerics_prompt);
	i = parse_command_line(cmd,numerics_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	if((numerics_interp(buffer,0)) == 1000)
	  break; 
      }
    }

  }
  else if((strcasecmp(buffer[depth_level],"run") == 0) || \
	  (strcasecmp(buffer[depth_level],"r") == 0)){
    /*Second argument of `run' must be number of simulations for the
    same init conditions, which is useful only for stochastic
    simulations*/
    if(strlen(cmd[1])==0)
      run(method,mynum.total_time,write_flag,0,1);
    else
      run(method,mynum.total_time,write_flag,0,atoi(cmd[1]));

  }
  else if((strcasecmp(buffer[depth_level],"runt") == 0) || \
	  (strcasecmp(buffer[depth_level],"rt") == 0)){
    run(method,mynum.trans_time,0,1,1);
				
  }
  else if((strcasecmp(buffer[depth_level],"runi") == 0) || \
	  (strcasecmp(buffer[depth_level],"ri") == 0)){
    for(i=0;i<DIM;i++){
      if(xin_par[i] >= 0){/* Check the init-param link */
	xin[i] = mu.P[xin_par[i]];
	yin[i] = mu.P[xin_par[i]];
      }
      x[i] = xin[i];
      y[i] = yin[i];
    }
    if(strlen(cmd[1])==0)
      run(method,mynum.total_time,write_flag,0,1);
    else
      run(method,mynum.total_time,write_flag,0,atoi(cmd[1]));

  }
  else if (strcasecmp(buffer[depth_level],"c") == 0 ||		\
	   strcasecmp(buffer[depth_level],"cal") == 0 ||	\
	   strcasecmp(buffer[depth_level],"calc") == 0 ||	\
	   strcasecmp(buffer[depth_level],"calculate") == 0){
    if(ma_span){/* Calculate MA if ma_span is non-zero */
      if((*(data_name+strlen(data_name)-1) == 'a') &&	\
	 (*(data_name+strlen(data_name)-2) == 'm') &&	\
	 (*(data_name+strlen(data_name)-3) == '.')){
	/*The file name contains .ma extension, perhaps already a MA trajectory*/
	printf("Looks like it is already an MA data.\n");
	printf("Do you still want to calculate MA?[Y/n]\n");
	if(getchar() == 'n'){
	  printf("Interrupting...\n");
	  free(symbuf);
	  return 0;
	}
      }
      /* Calculate the MA trajectory */
      moving_average(data_name,ma_span);
      /* Append the .ma extension for analyze_traj to analyze the .ma file. */
      symbuf = (char *)realloc(symbuf,(strlen(data_name)+4)*sizeof(char));
      symbuf = strcpy(symbuf,data_name);/* First copy data_name */
      symbuf = strcat(symbuf,".ma");/* ...then cat .ma extension */
    }
    else{/* Copy the data_name to symbuf to further analyze */
      symbuf = (char *)realloc(symbuf,(strlen(data_name)+1)*sizeof(char));
      symbuf = strcpy(symbuf,data_name);
    }
    if((analyze_traj(symbuf)) != 0){/* Recompute common entities */
      free(symbuf);
      return 10;
    }
    if(graph_flag){/* Plot the data */
      /* Remove .ma extension, since it will be added in gplot_results() */
      if(ma_span)
	*(symbuf+strlen(symbuf)-3) = '\0';/* This must be equal to data_name */
      if(strcmp(symbuf,data_name)!=0){
	printf("ERROR: %s != %s",symbuf,data_name);
      }
      fname = fopen(symbuf,"r");
      if(fname == NULL){
	fprintf(stderr,"Error(f ld): could not open file.\n");
	return 101;
      }
      if((info=get_info_data(info,fname)) == NULL){
	fprintf(stderr,"Error: could not get data info\n");
      }
      fclose(fname);
      if(info[0]){//method is complex
	/* First, plot determ trajectory */
	gplot_results(data_name,1,"rkf45");//any determ method, not "discrete"
	/* Plot stoch trajectory */
	if(!info[3])//lang_flag
	  gplot_results(data_name,2,"discrete");
	else
	  gplot_results(data_name,2,"rkf45");//any determ method, not "discrete"
      }
      else{//method is not complex
	if((!info[3]) && (!info[2]))//lang_flag and pure_det are false
	  gplot_results(data_name,0,"discrete");//...then discrete method
	else
	  gplot_results(data_name,0,"rkf45");//determ or langevin alone...
      }
      fprintf(stdout,"Graph is updated.\n");
    }
  }
  else if( (strcasecmp(buffer[depth_level],"parameters") == 0) || \
	   (strcasecmp(buffer[depth_level],"p") == 0) ) {
    get_parameter(buffer,depth_level+1);
  }
  else if( (strcasecmp(buffer[depth_level],"initials") == 0) || \
	   (strcasecmp(buffer[depth_level],"i") == 0) ) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter filename for initials:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if((*symbuf) != 0){
      if((load_initials(symbuf)) != 0){
	fprintf(stderr,"Error occurred.\n");
	free(symbuf);
	return 10;
      }
    }
    printf("Init conditions:\n");
    printf("Var. = Continuous | Discrete\n");
    for(i=0;i<DIM;i++)
      printf("%s(0) = %G | %ld\n",var_name[i],xin[i],yin[i]);
    printf("\n");
  }
  else if( (strcasecmp(buffer[depth_level],"periodics") == 0) || \
	   (strcasecmp(buffer[depth_level],"e") == 0) ) {
    if(strlen(buffer[depth_level+1])!=0)
      periods_interp(buffer,depth_level+1);
    else{
      while(1){
	//read_menu(cmd,periodics_prompt);
	i = parse_command_line(cmd,periodics_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	if((periods_interp(cmd,0)) == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"trajectory") == 0) || \
	   (strcasecmp(buffer[depth_level],"t") == 0) ) {
    if(strlen(buffer[depth_level+1])!=0)
      traj_interp(buffer,depth_level+1);
    else{
      while(1){
	//read_menu(cmd,traj_prompt);
	i = parse_command_line(cmd,traj_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	if((traj_interp(cmd,0)) == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"file") == 0) || \
	   (strcasecmp(buffer[depth_level],"f") == 0) ) {
    ret_num=0;
    if(strlen(buffer[depth_level+1])!=0)
      file_interp(buffer,depth_level+1);
    else{
      while(!ret_num){
	//read_menu(cmd,file_prompt);
	i = parse_command_line(cmd,file_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num=file_interp(cmd,0);
	if(ret_num == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"continue") == 0) || \
	   (strcasecmp(buffer[depth_level],"o") == 0) ){
    ret_num=0;
    if(strlen(buffer[depth_level+1])!=0)
      cont_interp(buffer,depth_level+1);
    else{
      while(!ret_num){
	//read_menu(cmd,cont_prompt);
	i = parse_command_line(cmd,cont_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num=cont_interp(cmd,0);
	if(ret_num == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"random") == 0) || \
	   (strcasecmp(buffer[depth_level],"a") == 0) ){
    ret_num=0;
    if(strlen(buffer[depth_level+1])!=0)
      rand_interp(buffer,depth_level+1);
    else{
      while(!ret_num){
	//read_menu(cmd,rand_prompt);
	i = parse_command_line(cmd,rand_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num=rand_interp(cmd,0);
	if(ret_num == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"errors")==0) || \
	   (strcasecmp(buffer[depth_level],"errors")==0) || \
	   (strcasecmp(buffer[depth_level],"er")==0)){
    ret_num=0;
    if(strlen(buffer[depth_level+1])!=0)
      errors_interp(buffer,depth_level+1);
    else{
      while(!ret_num){
	//read_menu(cmd,errors_prompt);
	i = parse_command_line(cmd,errors_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num = errors_interp(cmd,0);
	if(ret_num == 1000)
	  break; 
      }
    }
  }
  else if( (strcasecmp(buffer[depth_level],"sing")==0) || \
	   (strcasecmp(buffer[depth_level],"s")==0)){
    ret_num=0;
    if(strlen(cmd[1])!=0)
      sing_interp(cmd[1]);
    else{
      while(!ret_num){
	//read_menu(cmd,sing_prompt);
	i = parse_command_line(cmd,sing_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num=sing_interp(cmd[0]);
	if(ret_num==1000)
	  break;
      }
    }
  }
  else if ( (strcasecmp(buffer[depth_level],"lyap")==0) || \
	    (strcasecmp(buffer[depth_level],"l")==0) ){
    ret_num=0;
    if(strlen(buffer[depth_level+1]) != 0)
      lyap_interp(buffer,depth_level+1);
    else{
      while(!ret_num){
	//read_menu(cmd,"lyapunov");
	i = parse_command_line(cmd,"lyapunov",input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	ret_num = lyap_interp(cmd,0);
	if(ret_num==1000)
	  break;
      }
    }
  }
  else if (strcasecmp(buffer[depth_level],"rl")==0){
    run_lin();
  }
  else if( (strcasecmp(buffer[depth_level],"warranty")==0) || (strcasecmp(buffer[depth_level],"w")==0) ){
    fprintf(stdout,"THERE IS NO WARRANTY FOR THE PROGRAM, TO THE EXTENT PERMITTED BY\n");
    fprintf(stdout,"APPLICABLE LAW.  EXCEPT WHEN OTHERWISE STATED IN WRITING THE COPYRIGHT\n");
    fprintf(stdout,"HOLDERS AND/OR OTHER PARTIES PROVIDE THE PROGRAM \"AS IS\" WITHOUT WARRANTY\n");
    fprintf(stdout,"OF ANY KIND, EITHER EXPRESSED OR IMPLIED, INCLUDING, BUT NOT LIMITED TO,\n");
    fprintf(stdout,"THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR\n");
    fprintf(stdout,"PURPOSE.  THE ENTIRE RISK AS TO THE QUALITY AND PERFORMANCE OF THE PROGRAM\n");
    fprintf(stdout,"IS WITH YOU.  SHOULD THE PROGRAM PROVE DEFECTIVE, YOU ASSUME THE COST OF\n"); 
    fprintf(stdout,"ALL NECESSARY SERVICING, REPAIR OR CORRECTION.\n");
  }
  else if(strcasecmp(buffer[depth_level],"")==0) {}
  else {
    fprintf(stdout,"No such command\n");
  }
  free(symbuf);
  return 0;
}

int traj_interp(char **buffer,const int depth_level)
{
  int i,j,k,m;
  int Len = 20;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  regReport *report;
  if(strcasecmp(buffer[depth_level],"exit") == 0 ||   \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0){
    printf("Leaving TRAJECTORY menu...\n");
    free(symbuf);
    return 1000;
  }
  else if(strcasecmp(buffer[depth_level],"ls") == 0 ||	\
     strcasecmp(buffer[depth_level],"l") == 0){
    fprintf(stdout,"Trajectory menu:\n");
    fprintf(stdout,"(D)ynamical regime\n");
    fprintf(stdout,"(M)olecular distribution\n");
    fprintf(stdout,"(V)erbose level\n");
    fprintf(stdout,"(A)bsolute tolerance\n");
    fprintf(stdout,"(R)elative tolerance\n");
    //fprintf(stdout,"(R)elative (h)omogeneity tolerance\n");
  }
  else if(!strcasecmp(buffer[depth_level],"sh") ||	\
	  !strcasecmp(buffer[depth_level],"show")){
    fprintf(stdout,"Absolute tolerance: %G\n",eps_traj_abs);
    fprintf(stdout,"Relative tolerance: %G\n",eps_traj_rel);
    //fprintf(stdout,"Relative tolerance for homogeneity test: %G\n",eps_traj_rel_hg);
    fprintf(stdout,"Verbose level: %d\n",traj_verb);
  }
  /* else if(strcasecmp(buffer[depth_level],"rh") == 0){ */
  /*   read_next_command_wrapper(buffer,depth_level, */
  /* 			      "Relative tolerance for homogeneity", */
  /* 			      (void *)&eps_traj_rel_hg,1); */
  /* } */
  else if(strcasecmp(buffer[depth_level],"r") == 0){
    read_next_command_wrapper(buffer,depth_level,"Relative tolerance",
			      (void *)&eps_traj_rel,1);
  }
  else if(strcasecmp(buffer[depth_level],"a") == 0){
    read_next_command_wrapper(buffer,depth_level,"Absolute tolerance",
			      (void *)&eps_traj_abs,1);
  }
  else if(strcasecmp(buffer[depth_level],"v") == 0){
    read_next_command_wrapper(buffer,depth_level,"Verbose level (0,1,2)",
			      (void *)&traj_verb,0);
  }
  else if(strcasecmp(buffer[depth_level],"regime") == 0 ||	\
     strcasecmp(buffer[depth_level],"d") == 0){

    report = get_sys_dynamics(report,0);
    if(report != NULL)
      fprintf(stdout,"========\nDYNAMICS:\n%s\n",report_translate(*report));
    free_report(report);
  }
  else if(strcasecmp(buffer[depth_level],"m") == 0 ||	\
     strcasecmp(buffer[depth_level],"mol") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      i = 10;/* Number of bins for the histogram */
    }
    else{
      i = atoi(symbuf);
      if(i<=0){
	printf("Error(# bins): set to default.\n");
	i = 10;
      }
    }
    printf("# bins: %d\n",i);
    mol_dist(i);
  }
  else if(strcasecmp(buffer[depth_level],"") == 0){}
  else {
    fprintf(stdout,"No such command\n");
  }
  free(symbuf);
  return 0;
}

int graphics_interp(char **buffer,const int depth_level)
{
  int i,j;
  int Len = 200;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  char *p,*p1;
  if(strcasecmp(buffer[depth_level], "exit") == 0 || \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving GRAPHICS menu...\n");
    free(symbuf);
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || \
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"Graphics menu:\n");
    fprintf(stdout,"(Sh)ow summary\n");
    fprintf(stdout,"(G)raphics toggle\n");
    fprintf(stdout,"Gnuplot (p)rompt\n");
    fprintf(stdout,"(X)-axis\n");
    fprintf(stdout,"(Y)-axis\n");
    fprintf(stdout,"G(r)id toggle\n");
    fprintf(stdout,"(S)end to eps\n");

  }
  else if(strcasecmp(buffer[depth_level],"send") == 0 ||	\
	  strcasecmp(buffer[depth_level],"s") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Dump the current graph to eps.\n");
      fprintf(stdout,"File name for the .eps:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No changes.\n");
    }else {
      send_to_eps(symbuf);
    }

  }
  else if(strcasecmp(buffer[depth_level],"show") == 0 ||	\
	  strcasecmp(buffer[depth_level],"sh") == 0){
    fprintf(stdout,"Graphics is ");
    if(graph_flag)
      fprintf(stdout,"on\n");
    else
      fprintf(stdout,"off\n");
    fprintf(stdout,"Graphics output: X11\n");
    if(graph.xInd == -1)
      fprintf(stdout,"X-axis: Time, t\n");
    else
      fprintf(stdout,"X-axis: U(%d) = %s\n",graph.xInd+1,var_name[graph.xInd]);
    i = 0;
    fprintf(stdout,"Y-axis: ");
    while(i<3){
      if(graph.yInd[i] > -1)
	fprintf(stdout,"U(%d) = %s; ",graph.yInd[i]+1,var_name[graph.yInd[i]]);
      i++;
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"Grid: ");
    if(graph.grid_flag)
      fprintf(stdout,"on\n");
    else
      fprintf(stdout,"off\n");

  }
  else if(strcasecmp(buffer[depth_level],"g") == 0 ||	\
	  strcasecmp(buffer[depth_level],"graph") == 0){
    if(!HAVE_GNUPLOT)
      fprintf(stdout,"Error: Dinamica was compiled without support of Gnuplot.\n");
    else{
      if((symbuf = read_next_command(buffer,depth_level)) == NULL){
	if(graph_flag){
	  graph_flag = 0;
	  fprintf(stdout,"Graphics is off\n");
	}
	else{
	  graph_flag = 1;
	  fprintf(stdout,"Graphics is on\n");
	}
      }
      else{
	graph_flag = atoi(symbuf);
	if(graph_flag)
	  printf("Graphics is on\n");
	else
	  printf("Graphics is off\n");
      }
    }
  }
  else if(strcasecmp(buffer[depth_level],"grid") == 0 ||	\
	  strcasecmp(buffer[depth_level],"r") == 0){
    if(graph.grid_flag){
      graph.grid_flag = 0;
      fprintf(stdout,"Grid is off\n");
    }
    else{
      graph.grid_flag = 1;
      fprintf(stdout,"Grid is on\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"x") == 0 ||		\
	  strcasecmp(buffer[depth_level],"xaxis") == 0 ||	\
	  strcasecmp(buffer[depth_level],"xax") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Type var name or index (0 = t, 1 = U(1), etc.)\n");
      fprintf(stdout,"X-axis:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No changes.\n");
    }else {
      if(isalpha(symbuf[0])){
	if(isVar(symbuf,var_name,DIM) > -1)
	  graph.xInd = isVar(symbuf,var_name,DIM);
	else
	  fprintf(stdout,"Do not know the variable\n");
      }
      else{
	if((atoi(symbuf) > DIM) || (atoi(symbuf) < 0))
	  fprintf(stdout,"Error: out of the range\n");
	else
	  graph.xInd = atoi(symbuf) - 1;
      }
    }
    if(graph.xInd == -1)
      fprintf(stdout,"X-axis: Time, t\n");
    else
      fprintf(stdout,"X-axis: U(%d) = %s\n",graph.xInd+1,var_name[graph.xInd]);
    graph_set_labels(graph,"tsphase");

  }
  else if(strcasecmp(buffer[depth_level],"y") == 0 ||		\
	  strcasecmp(buffer[depth_level],"yaxis") == 0 ||	\
	  strcasecmp(buffer[depth_level],"yax") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Type var name or index (0 = disable, 1 = U(1), etc.)\n");
      fprintf(stdout,"Multiple choices can be separated by space.\n");
      fprintf(stdout,"Y-axis:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No changes.\n");
    }else {
      i = 0;
      p1 = symbuf;
      while(i<3){
	if((p = strchr(p1,' ')) != 0){/* We can find space */
	  *p = 0;/* then terminate at space */
	}
	if(isalpha(*p1)){
	  if(isVar(p1,var_name,DIM) > -1)
	    graph.yInd[i] = isVar(p1,var_name,DIM);
	  else
	    fprintf(stdout,"Do not know the variable\n");
	}
	else{
	  if((atoi(p1) > DIM) || (atoi(p1) < 0))
	    fprintf(stdout,"Error: out of the range\n");
	  else{
	    /* Y-axis is not allowed to be time */
	    if(atoi(p1) == 0)
	      graph.yInd[i] = -2;/* Disabled */
	    else
	      graph.yInd[i] = atoi(p1) - 1;
	  }
	}
	if(p)/* p is not NULL */
	  p1 = p + 1;/* Move position right after the last found space */
	else /* p is NULL */
	  break;
	i++;
      }
    }
    i = 0;
    fprintf(stdout,"Y-axis: ");
    while(i<3){
      if(graph.yInd[i] > -1)
	fprintf(stdout,"U(%d) = %s; ",graph.yInd[i]+1,var_name[graph.yInd[i]]);
      i++;
    }
    fprintf(stdout,"\n");
    graph_set_labels(graph,"tsphase");

  }
  else if(strcasecmp(buffer[depth_level],"p") == 0 ||		\
	  strcasecmp(buffer[depth_level],"gnuplot") == 0 ||	\
	  strcasecmp(buffer[depth_level],"prompt") == 0){
    gnuplot_interp();

  }
  else if(strcasecmp(buffer[depth_level],"")==0) {}
  else {
    fprintf(stdout,"No such command\n");
  }
  free(symbuf);
  return 0;
}

int numerics_interp(char **buffer,const int depth_level)
{
  int i,j;
  int Len=20;
  char *symbuf=(char *)malloc(Len*sizeof(char));
  /* Store the Jacobian values */
  double dfdx[DIM*DIM];
  double dfdt[DIM];
  if(strcasecmp(buffer[depth_level], "exit") == 0 || \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving NUMERICS menu...\n");
    free(symbuf);
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || \
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"Numerics menu:\n");
    fprintf(stdout,"(Tt)ime\n");
    fprintf(stdout,"(Tr)time\n");
    fprintf(stdout,"(T)ra(j) transient time\n");
    fprintf(stdout,"(S)tep\n");
    fprintf(stdout,"N(p)ar\n");
    fprintf(stdout,"(B)uffer\n");
    fprintf(stdout,"Buffer (I)ncrement\n");
    fprintf(stdout,"(W)riting step\n");
    fprintf(stdout,"(S)ampling (f)requency\n");
    fprintf(stdout,"(M)ethod\n");
    fprintf(stdout,"(J)acobian test\n");
    
  }
  else if(strcasecmp(buffer[depth_level],"j") == 0){
    printf("General form of an equation:\n");
    printf("dX/dt = f(X) + g(X)*<Wiener process>\n");
    printf("Current point:\n");
    printf("(");
    for(i=0;i<DIM;i++){
      if(i==DIM-1)
	printf("%s)",var_name[i]);
      else
	printf("%s,",var_name[i]);
    }
    printf(" = ");
    printf("(");
    for(i=0;i<DIM;i++){
      if(i==DIM-1)
	printf("%.4lf)",x[i]);
      else
	printf("%.4lf,",x[i]);
    }
    printf("\n");
    printf("Jacobian for f(X):\n");
    jac_general(t,x,dfdx,dfdt,&mu,func_odeiv);
    /* Printing */
    for(i=0;i<DIM;i++)
      printf("\t%s",var_name[i]);
    printf("\n");
    for(i=0;i<DIM;i++){
      for(j=0;j<DIM;j++){
	if(j==0)
	  printf("(%s')%.4lf ",var_name[i],dfdx[i*DIM+j]);
	else
	  printf("%.4lf ",dfdx[i*DIM+j]);
      }
      printf("\n");
    }
    printf("Jacobian for g(X):\n");
    jac_general(t,x,dfdx,dfdt,&mu,lang_amend);
        /* Printing */
    for(i=0;i<DIM;i++)
      printf("\t%s",var_name[i]);
    printf("\n");
    for(i=0;i<DIM;i++){
      for(j=0;j<DIM;j++){
	if(j==0)
	  printf("(%s')%.4lf ",var_name[i],dfdx[i*DIM+j]);
	else
	  printf("%.4lf ",dfdx[i*DIM+j]);
      }
      printf("\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"ttime") == 0 || \
	  strcasecmp(buffer[depth_level],"tt") == 0) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter total time:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored. Total time unchanged.\n");
    }else {
      if((atof(symbuf)) < 0)
	fprintf(stdout,"Time must be positive.\n");
      else if ((atof(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero.\n");
      else
	mynum.total_time = atof(symbuf);
      fprintf(stdout,"Total time of integration now is %G\n",mynum.total_time);}
    if((mynum.total_time/(mynum.step*mynum.write_step)>mynum.global_buffer\
	&& mynum.write_step!=0) \
       && (strcasecmp(method,"run-kut4")==0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen total time,increase step\n");
    }

  }
  else if((strcasecmp(buffer[depth_level],"trtime")==0) || \
	  (strcasecmp(buffer[depth_level],"tr") == 0)) {
   if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter transient time:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored.Transient time unchanged.\n");
    }else{
      if((atof(symbuf)) < 0)
	fprintf(stdout,"Trans time must be positive.\n");
      else if ((atof(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero.\n");
      else
	mynum.trans_time = atof(symbuf);
      fprintf(stdout,"Transient time of integration now is %G\n",mynum.trans_time);}
    if((mynum.trans_time/(mynum.step*mynum.write_step) > mynum.global_buffer \
  && mynum.write_step != 0) \
       && (strcasecmp(method,"run-kut4") == 0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen transition time,increase step\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"tj")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter new value:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0)
      fprintf(stdout,"No changes.\n");
    else{
      if((atoi(symbuf)) < 0)
	fprintf(stdout,"Error: must be positive.\n");
      else if ((atoi(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero.\n");
      else
	traj_trans = atoi(symbuf);
      fprintf(stdout,"Trajectory trans is %d\n",traj_trans);
    }

  }
  else if(strcasecmp(buffer[depth_level],"step") == 0 || \
	  strcasecmp(buffer[depth_level],"s") == 0) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter step:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored.Step unchanged.\n");
    }else{
      if((atof(symbuf)) < 0){
	fprintf(stdout,"Step must be positive.Sign changed.\n");
	mynum.step = -atof(symbuf);}
      else if((atof(symbuf)) == 0)
	fprintf(stdout,"Step must not be equal zero.Left unchanged.\n");
      else
	mynum.step = atof(symbuf);
      fprintf(stdout,"Step of integration now is %G\n",mynum.step);}
    if((mynum.total_time/(mynum.step*mynum.write_step) > mynum.global_buffer \
  && mynum.write_step!=0) \
       && (strcasecmp(method,"run-kut4") == 0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen total time,increase step\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"buffer") == 0 || \
	  strcasecmp(buffer[depth_level],"b") == 0) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter number of points to be stored(global buffer):\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored.Global buffer unchanged.\n");
    }else{
      if((atoi(symbuf)) < 0)
	fprintf(stdout,"Error: less than zero.Left unchanged.\n");
      else if ((atoi(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero.\n");
      else{
	BUFFER = atoi(symbuf);
	mynum.global_buffer = BUFFER;}
      fprintf(stdout,"Global buffer now is %d\n",BUFFER);}
    /*Reallocating memory for storage arrays*/
    ts = (double *)realloc(ts,mynum.global_buffer*sizeof(double));
    xs = (double **)realloc(xs,mynum.global_buffer*sizeof(double *));
    for(i=0;i<mynum.global_buffer;i++)
      xs[i] = (double *)realloc(xs[i],DIM*sizeof(double));
    if((mynum.total_time/(mynum.step*mynum.write_step) > mynum.global_buffer \
  && mynum.write_step != 0)\
       && (strcasecmp(method,"run-kut4") == 0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen total time,increase step\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"i")==0 || \
	  strcasecmp(buffer[depth_level],"increment")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter number of points to be incremented:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf)==0){
      fprintf(stdout,"No value was stored. Buffer Increment unchanged.\n");
    }else{
      if((atoi(symbuf)) < 0)
	fprintf(stdout,"Error: less than zero.Left unchanged.\n");
      else if((atoi(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero\n");
      else
	BUF_INCR = atoi(symbuf);
      fprintf(stdout,"Buffer increment now is %d\n",BUF_INCR);
    }

  }
  else if(strcasecmp(buffer[depth_level],"w") == 0 || \
	  strcasecmp(buffer[depth_level],"ws") == 0) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter step for writing to the file:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored.\n");
    }else{
      if((atoi(symbuf)) < 0)
	fprintf(stdout,"Error: cannot be less than zero.Left unchanged.\n");
      else if((atoi(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero\n");
      else
	mynum.write_step = atoi(symbuf);
      fprintf(stdout,"Writing step now is %d\n",mynum.write_step);
    }
    if((mynum.total_time/(mynum.step*mynum.write_step) > mynum.global_buffer) \
  && mynum.write_step != 0 \
       && (strcasecmp(method,"run-kut4") == 0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen total time,increase step\n");
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"sf") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      printf("Enter sampling frequency of discrete simulations:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored. Sampling unchanged.\n");
    }else{
      if((atof(symbuf)) < 0){
	fprintf(stdout,"Cannot be less than zero.Invert.\n");
	mynum.smp_frq = -atof(symbuf);
      }
      else if((atof(symbuf)) == 0)
	fprintf(stdout,"Error: must not be zero\n");
      else
	mynum.smp_frq = atof(symbuf);
      fprintf(stdout,"Sampling frequency now is %.2lf\n",mynum.smp_frq);
    }

  }
  else if(strcasecmp(buffer[depth_level],"method") == 0 || \
	  strcasecmp(buffer[depth_level],"m") == 0) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"(0) Euler 2-nd order\n");
      fprintf(stdout,"(1) Runge-Kutta 4-th order\n");
      fprintf(stdout,"(2) Embedded Runge-Kutta-Fehlberg (4,5) method\n");
      fprintf(stdout,"(3) Embedded Runge-Kutta Prince-Dormand (8,9) method\n");
      fprintf(stdout,"(4) Embedded Runge-Kutta Cash-Karp (4,5) method\n");
      fprintf(stdout,"(5) Discrete\n");
      fprintf(stdout,"(6) Toggle Langevin flag\n");
      fprintf(stdout,"(7) Complex\n");
      fprintf(stdout,"(8) Implicit Bulirsch-Stoer method\n");
      fprintf(stdout,"(9) Milstein method for stochastic ODE\n");
      fprintf(stdout,"Choose one.\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No value was stored. Method unchanged.\n");
    }
    else if(isalpha(*symbuf)){
      fprintf(stdout,"Should be numeric. Method unchanged.\n");
    }
    else{
      j = atoi(symbuf);
      if(j < 0){
	fprintf(stdout,"Error: less than zero. Left unchanged.\n");
      }
      if(j == 0) strcpy(method,"eu");
      if(j == 1) strcpy(method,"run-kut4");
      if(j == 2) strcpy(method,"rkf45");
      if(j == 3) strcpy(method,"rk8pd");
      if(j == 4) strcpy(method,"rkck");
      if(j == 5) strcpy(method,"discrete");
      if(j == 6) {
	if(strcmp(method,"discrete")==0){
	  printf("Langevin is not valid for discrete method. Left unchanged.\n");
	  free(symbuf);
	  return 0;
	}
	if(lang_flag==0){
	  lang_flag=1; printf("Langevin = true\n");}
	else{
	  if(!strcmp(method,"milst"))
	    printf("Cannot switch off the Langevin flag. Method is \"Milstein\".\n");
	  else{
	    lang_flag=0; printf("Langevin = false\n");
	  }
	}
      }
      if(j == 7){
	if(lang_flag==0 && (strcmp(method,"discrete")==0)){
	  /* Ask about deterministic method to use along with discrete
	     one */
	  printf("I will use default rkf45 method along with discrete one in this case\n");
	  printf("For other deterministic method:\n");
	  printf("1) Choose that method.\n");
	  printf("2) Toggle complex (langevin must be false).\n");
	  strcpy(method2,"rkf45");
	  strcpy(method,"complex");
	  /* Complex + lang_flag==0 => rkf45 + discrete */
	}
	else if(strcmp(method,"complex")==0){
	  printf("Error: method is already complex\n");
	}
	else if(lang_flag==0 && (strcmp(method,"discrete")!=0)){
	  /* Deterministic method is chosen, set up the discrete */
	  /* Complex + lang_flag==0 => `current determ. method' + discrete */
	  strcpy(method2,method);
	  strcpy(method,"complex");
	}
	else if(lang_flag!=0 && (strcmp(method,"discrete")==0)){
	  /* Error */
	  printf("Error: no complex option for two stochastic methods\n");
	}
	else if(lang_flag!=0 && (strcmp(method,"discrete")!=0)){
	  /* Complex + lang_flag==1 => Stoch Diff Eq */
	  strcpy(method2,method);
	  strcpy(method,"complex");
	}
      }
      if(j == 8) {
	strcpy(method,"bsimp");
	printf("Warning: this method requires Jacobian calculation.\n");
      }
      if(j == 9){
	if(strcmp(method,"discrete"))
	  strcpy(method,"milst");
	else{
	  printf("The \"discrete\" method is not valid. Turning off.\n");
	  strcpy(method,"milst");
	}
	if(!lang_flag){
	  printf("Langevin flag must be set. Turning on.\n");
	  lang_flag = 1;
	}
      }
      if(j > 9){
	fprintf(stdout,"Error: more than we have.Left unchanged.\n");
      }
      if(strcmp(method,"complex")!=0)
	fprintf(stdout,"Method is %s\n",method);
      else{
	if(lang_flag==0)
	  fprintf(stdout,"Method is %s:%s/discrete\n",method,method2);
	else
	  fprintf(stdout,"Method is %s:%s/Langevin\n",method,method2);
      }
    }
    if((mynum.total_time/(mynum.step*mynum.write_step) > mynum.global_buffer \
	&& mynum.write_step != 0)					\
       && (strcasecmp(method,"run-kut4") == 0)){
      fprintf(stdout,"Warning: You cannot get all points to storage array! Some of data will be lost\n");
      fprintf(stdout,"Lessen total time,increase step or write step,change algorithm.\n");
    }

  }
  else if(strcasecmp(buffer[depth_level],"show") == 0 || \
	  strcasecmp(buffer[depth_level],"sh") == 0){
    printf("* Dimension: %d\n",DIM);
    printf("* Number of systems: %d\n",LDIM);
    printf("* Number of parameters: %d\n",PARDIM);
    if(NREAC)
      printf("* Number of reactions: %d\n",NREAC);
    printf("* Number of user functions: %d\n",FNUM);
    printf("* Number of auxillary entities: %d\n",AUXNUM);
    printf("***************************\n");
    printf("Total time: %G\n",mynum.total_time);
    printf("Transient time: %G\n",mynum.trans_time);
    printf("Step: %G\n",mynum.step);
    printf("Writing step: %d\n",mynum.write_step);
    printf("Sampling frequency: %.2lf\n",mynum.smp_frq);
    printf("Method: %s\n",method);
    if(strcmp(method,"complex")==0)
      printf("|->Method2: %s\n",method2);
    if(lang_flag==0)
      printf("Langevin flag: false\n");
    else
      printf("Langevin flag: true\n");

  }
  else if(strcasecmp(buffer[depth_level],"")==0){}
  else {
    fprintf(stdout,"No such command\n");
  }
  free(symbuf);
  return 0;
}

int periods_interp(char **buffer, const int depth_level) 
{
  int i,j,k,l,m;
  double tmp;
  int Len = 20;
  char *p;
  char * symbuf = (char *)malloc(Len*sizeof(char));
  if(strcasecmp(buffer[depth_level],"exit") == 0 || \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving PERIODICS menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || \
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"PERIODICS menu:\n");
    fprintf(stdout,"(R)e(c)ompute\n");
    fprintf(stdout,"(P)eriod\n");
    fprintf(stdout,"(C)ross level\n");
    fprintf(stdout,"Period (v)ariable\n");
    fprintf(stdout,"Period t(h)reshold\n");
    fprintf(stdout,"(MA) span\n");
    fprintf(stdout,"(S)ubperiod\n");
    fprintf(stdout,"(R)atio\n");
    fprintf(stdout,"A(m)plitude\n");
    fprintf(stdout,"Me(t)hod\n");
    fprintf(stdout,"(A)utocorrelation (p)lot\n");
    fprintf(stdout,"(Tm)ap\n");
    fprintf(stdout,"(Th)istogram\n");
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"ma") == 0) || \
	   (strcasecmp(buffer[depth_level],"span") == 0) ){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"MA span = %d\n",ma_span);
      fprintf(stdout,"Enter moving average span:\n");
      fprintf(stdout,"(0 for disabling MA calculation)\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No change.\n");
    }
    else{
      if(atoi(symbuf) == 0){
	ma_span = atoi(symbuf);
	fprintf(stdout,"MA disabled.\n");
      }
      else if(atoi(symbuf) < 0){
	fprintf(stdout,"Error: Span cannot be negative integer\n");
	fprintf(stdout,"Left unchanged.\n");
      }
      else if(fmod(atoi(symbuf),2) < 0.001){
	printf("Span is even integer. Adding unity.\n");
	ma_span = atoi(symbuf) + 1;
	fprintf(stdout,"MA span = %d\n",ma_span);
      }
      else{
	ma_span = atoi(symbuf);
	fprintf(stdout,"MA span = %d\n",ma_span);
      }
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"h") == 0) || \
	   (strcasecmp(buffer[depth_level],"thresh") == 0) ){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Period threshold = %G\n",per_thresh);
      fprintf(stdout,"Enter period threshold:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No change.\n");
    }
    else{
      per_thresh = atof(symbuf);
      if(per_thresh > 0)
	fprintf(stdout,"Period threshold = %G\n",per_thresh);
      else{
	fprintf(stdout,"Negative value!Set to default.\n");
	per_thresh = 0;
      }
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"rc") == 0) || \
	   (strcasecmp(buffer[depth_level],"recomp") == 0) ) {
    printf("%s will be analyzed...\n",data_name);
    if((analyze_traj(data_name)) != 0)/* Recompute common entities */
      return 10;
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"period") == 0) || \
	   (strcasecmp(buffer[depth_level],"p") == 0) ) {
    fprintf(stdout,"Periods by ");
    if(per_method == 0)
      fprintf(stdout,"Poincare sections (cross=%G) ",cross);
    if(per_method == 1)
      fprintf(stdout,"autocorrelation ");
    if(nPerDet || nPerStoch)
      fprintf(stdout,"(U(%d)=%s):\n",perVarInd+1,var_name[perVarInd]);
    else
      fprintf(stdout,"\n");
    /* Deterministic periods */
    if(nPerDet){
      fprintf(stdout,"Determ.: mean=%G, std=%G, skew=%G, n=%d\n",
	      gsl_stats_mean(perDet,1,nPerDet),
	      gsl_stats_sd(perDet,1,nPerDet),
	      gsl_stats_skew(perDet,1,nPerDet),
	      nPerDet);
      if(nPerDet > 10){/* Truncate the output */
	fprintf(stdout,"(Output is truncated to 10 values)\n");
	for(i=0;i<5;i++)
	  fprintf(stdout,"%G ",perDet[i]);
	fprintf(stdout,"... ");
	for(i=nPerDet-5;i<nPerDet;i++)
	  fprintf(stdout,"%G ",perDet[i]);
	fprintf(stdout,"\n");
      }
      else{
	for(i=0;i<nPerDet;i++)
	  fprintf(stdout,"%G ",perDet[i]);
	fprintf(stdout,"\n");
      }
    }
    /* Stochastic periods */
    if(nPerStoch){
      fprintf(stdout,"Stoch.: mean=%G, std=%G, skew=%G, n=%d\n",
	      gsl_stats_mean(perStoch,1,nPerStoch),
	      gsl_stats_sd(perStoch,1,nPerStoch),
	      gsl_stats_skew(perStoch,1,nPerStoch),
	      nPerStoch);
      if(nPerStoch > 10){/* Truncate the output */
	fprintf(stdout,"(Output is truncated to 10 values)\n");
	for(i=0;i<5;i++)
	  fprintf(stdout,"%G ",perStoch[i]);
	fprintf(stdout,"... ");
	for(i=nPerStoch-5;i<nPerStoch;i++)
	  fprintf(stdout,"%G ",perStoch[i]);
	fprintf(stdout,"\n");
      }
      else{
	for(i=0;i<nPerStoch;i++)
	  fprintf(stdout,"%G ",perStoch[i]);
	fprintf(stdout,"\n");
      }
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"subperiod") == 0) || \
	   (strcasecmp(buffer[depth_level],"s") == 0) ) {
    fprintf(stdout,"Spectra:\n");
    if(per_method == 0){
      for(k=1; k <= DIM; k++){
	fprintf(stdout,"U(%d)=%s\n",k,var_name[k-1]);
	fprintf(stdout,"Subperiods: %d\n",n_subperiods[k-1]);
	for(i=1; i <= n_isec[k-1]-2; i++){
	  fprintf(stdout,"%d ",spectrum_per[i-1][k-1]);
	  if(i == n_isec[k-1]-2) fprintf(stdout,"\n");
	}
      }
    }
    if(per_method == 1){
      for(k=1; k <= DIM; k++){
	fprintf(stdout,"U(%d)=%s\n",k,var_name[k-1]);
	fprintf(stdout,"Subperiods: %d\n",n_subperiods[k-1]);
	for(i=1; i <= n_peaks[k-1]-2; i++){
	  fprintf(stdout,"%d ",spectrum_per[i-1][k-1]);
	  if(i == n_peaks[k-1]-2) fprintf(stdout,"\n");
	}
      }
    }
    return 0;
  }

  else if( (strcasecmp(buffer[depth_level],"cross") == 0) ||	\
	   (strcasecmp(buffer[depth_level],"c") == 0) ) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Cross level(s) = ");
      for(j=0;j<3;j++){
	if(fabs(cross_level[j]) < 0.00001)/* Zero! */
	  break;
	if(cross_level[j] > 0)
	  fprintf(stdout,"%G ",cross_level[j]);
      }
      fprintf(stdout,"\n");
      fprintf(stdout,"[Poincare section = (max - min)/cross_level + min]\n");
      fprintf(stdout,"Enter cross level(2):\n");
      fprintf(stdout,"[Space delimits values (max 3). 0 to disable.]\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No change.\n");
    }
    else{
      j = 0;/* number of values */
      while((*symbuf) != '\0'){
	if(j > 2){
	  fprintf(stderr,"Error: cannot read more than 3 values.\n");
	  break;
	}
	p = 0;
	if(strchr(symbuf,' ') != NULL){
	  p = strchr(symbuf,' ');
	  *p = 0;
	}
	cross_level[j] = atof(symbuf);
	if(cross_level[j] < 0){
	  fprintf(stdout,"Negative value!Set to default.\n");
	  cross_level[j] = 4;
	}
	j++;

	if(p)
	  symbuf = p+1;
	else
	  break;
	while((*symbuf) == ' ')
	  symbuf++;
      }
      fprintf(stdout,"Cross level(s) = ");
      for(j=0;j<3;j++){
	if(fabs(cross_level[j]) < 0.00001)/* Zero! */
	  break;
	if(cross_level[j] > 0)
	  fprintf(stdout,"%G ",cross_level[j]);
      }
      fprintf(stdout,"\n");
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"v") == 0) ||	\
	   (strcasecmp(buffer[depth_level],"pv") == 0) ) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Periods are calculated for the variable.\n");
      fprintf(stdout,"Choose 0 for automatic detection.\n");
      fprintf(stdout,"Period variable = ");
      if(per_var == -1)
	fprintf(stdout,"automatic (max amplitude var)\n");
      else
	fprintf(stdout,"%s\n",var_name[per_var]);
      fprintf(stdout,"Enter period variable:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No change.\n");
    }
    else{
      if(isalpha(*symbuf))
	per_var = isVar(symbuf,var_name,DIM);
      else
	per_var = atoi(symbuf) - 1;
      if(per_var < -1 || per_var > (DIM-1)){
	fprintf(stderr,"Error: wrong input. Left default.\n");
	per_var = -1;
      }
      fprintf(stdout,"Period variable = ");
      if(per_var == -1)
	fprintf(stdout,"automatic(max amplitude chosen)\n");
      else
	fprintf(stdout,"%s\n",var_name[per_var]);
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"ratio") == 0) || \
	   (strcasecmp(buffer[depth_level],"r") == 0) ) {
    fprintf(stdout,"Ratio:\n");
    for(i=0;i<DIM;i++){
      fprintf(stdout,"U(%d)=%s: ",i+1,var_name[i]);
      fprintf(stdout,"%lf\n",per_ratio[i]);
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"amplitude") == 0) || \
	   (strcasecmp(buffer[depth_level],"m") == 0) ) {
    read_next_command_wrapper(buffer,depth_level,"number of bins",
			      (void *)&m,0);
    double *ampl,**hist;
    int ampl_size;
    int bins = m;
    ampl = ampl_multi_frames(ampl,&ampl_size);
    printf("Amplitudes for %s:\n",var_name[graph.yInd[0]]);
    for(i=0;i<ampl_size;i++){
      if(i >= 10){
	printf("...Truncated to first 10 values.\n");
	break;
      }
      printf("%G ",ampl[i]);
    }
    printf("\n==========\n");
    printf("Mean=%G, std=%G, skew=%G (n=%d)\n",gsl_stats_mean(ampl,1,ampl_size),
	   gsl_stats_sd(ampl,1,ampl_size),
	   gsl_stats_skew(ampl,1,ampl_size),ampl_size);
    /* Calculate histogram */
    printf("Using %d bins (m=%d)\n",bins,m);
    hist = compute_hist_per(hist,ampl,ampl_size,bins);
    write_histPer_file(hist,bins,"hist.dat","w");
    /* Print the histogram */
    if(graph_flag){
      graph_set_labels(graph,"ampldist");
      gnuplot_cmd(plot_handle,"set title \"n=%d\"\n",ampl_size);
      gnuplot_cmd(plot_handle,"plot 'hist.dat' u 1:2 w boxes ls 1\n");
    }
    for(i=0;i<2;i++)
      free(hist[i]);
    free(hist);
    free(ampl);
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"method") == 0) || \
	   (strcasecmp(buffer[depth_level],"t") == 0) ){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Choose the method for period calculation:\n");
      if(per_method == 0){
	fprintf(stdout,"(0)By Poincare section(default)*\n");
	fprintf(stdout,"(1)By auto-correlation\n");
      }
      if(per_method == 1){
	fprintf(stdout,"(0)By Poincare section(default)\n");
	fprintf(stdout,"(1)By auto-correlation*\n");
      }
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }else{
      i = per_method;
      per_method = atoi(symbuf);
      if((per_method != 0) && (per_method != 1)){
	fprintf(stdout,"Wrong choice! Unchanged.\n");
	per_method = i;
      }
      printf("Method is: ");
      if(per_method == 0) printf("section\n");
      if(per_method == 1) printf("auto-correlation\n");
    }
    return 0;
  }
  else if ((strcasecmp(buffer[depth_level],"tmap")==0) || \
	   (strcasecmp(buffer[depth_level],"tm")==0)){
    if(nPerStoch){
      fprintf(stdout,"Number of periods to analyze: %d\n",nPerStoch);
      fprintf(stdout,"Output file name is `tmap.dat'\n");
      FILE *ftmap;
      ftmap=fopen("tmap.dat","w");
      for(i=0;i<nPerStoch-1;i++)
	fprintf(ftmap,"%lf %lf\n",perStoch[i],perStoch[i+1]);
      fclose(ftmap);
      graph_set_labels(graph,"tmap");
      gnuplot_cmd(plot_handle,"plot 'tmap.dat' u 1:2 w p ps 2\n");
    }
    else{
      fprintf(stderr,"Error: Could not find the calculated periods.\n");
    }
    return 0;
  }
  else if (!strcasecmp(buffer[depth_level],"ap")){
    /* Plot the autocorrelation function */
    FILE *acf;
    if((acf = fopen("acorr.dat","r")) == NULL){
      fprintf(stderr,"Error: file \"acorr.dat\" is not available.\n");
      free(symbuf);
      return 0;
    }
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   printf(""); */
    /* } */
    graph_set_labels(graph,"ac");
    gnuplot_cmd(plot_handle,"plot 'acorr.dat' i 0 w p\n");
    free(symbuf);
    return 0;
  }
  else if ((strcasecmp(buffer[depth_level],"th")==0) || \
	   (strcasecmp(buffer[depth_level],"thist")==0) || \
	   (strcasecmp(buffer[depth_level],"thistogram")==0)){
    if(strlen(buffer[depth_level+1]) != 0){
      periodics_hist_interp(buffer,depth_level+1);
    }
    else{
      while(1){
	//read_menu(buffer,periodics_hist_prompt);
	i = parse_command_line(cmd,periodics_hist_prompt,input_stream);
	if(i == 100)
	  break;
	if(i == 200)
	  continue;
	if(periodics_hist_interp(buffer,0) == 1000)
	  break;
      }
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

int file_interp(char **buffer, const int depth_level)
{
  int i,j,*info;
  FILE *fname;
  int Len = 100;//filename can be long with directory specification
  char *symbuf = (char *)malloc(Len*sizeof(char));
  if(strcasecmp(buffer[depth_level], "exit") == 0 || strcasecmp(buffer[depth_level],"quit") == 0 || strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving FILE menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || \
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"FILE menu:\n");
    fprintf(stdout,"(S)av(i)n - save current point to initials\n");
    fprintf(stdout,"(S)av(d)at - write current data,trajectory file\n");
    fprintf(stdout,"(L)oad (d)ata - load data file, recompute common entities\n");
    fprintf(stdout,"(L)oad (o)utput data - load output data file with periods\n");
    fprintf(stdout,"(C)onvert - convert data file\n");
    fprintf(stdout,"(D)ata file - change data file to write for next runs\n");
    fprintf(stdout,"(O)utput modifier - appended to the data filename\n");
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"lo")==0){
    char *out = (char *)malloc((strlen(data_name)+30)*sizeof(char));
    char *tmp = (char *)malloc(25*sizeof(char));
    if((check_next_arg(buffer,depth_level,out,Len) != 0)){
      /* If not provided in the command, guess the name */
      sprintf(tmp,"_c%G",cross_level[0]);
      out = strcpy(out,data_name);
      if(ma_span)
	out = strcat(out,".ma");
      out = strcat(out,tmp);
    }
    fprintf(stdout,"File to look: %s\n",out);
    i = read_out_file(out);
    if(!i)
      fprintf(stdout,"# stoch periods found: %d\n",nPerStoch);
    free(out);
    free(tmp);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"lf")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Current file to look into: %s\n",data_name);
      fprintf(stdout,"How many frames to load:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){/* No input */
      fprintf(stdout,"No changes\n");
    }
    else{
      double **frame;
      int size_frame;
      FILE *out = fopen("frame","w");
      fname = fopen(data_name,"r");
      for(i=0;i<atoi(symbuf);i++){
	frame = load_frame(fname,frame,&size_frame);
	if(frame == NULL)
	  fprintf(stderr,"Perhaps you want more frames than there exist\n");
	write_frame(out,frame,size_frame);
      }
      free(frame);
      fclose(fname);
      fclose(out);
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"o")==0 || \
	  strcasecmp(buffer[depth_level],"out")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Current output modifier: %s\n",out_name);
      fprintf(stdout,"Enter new name:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No changes\n");
      fprintf(stdout,"Output modifier: %s\n",out_name);
    }
    else{
      if(Len < FNAME_SIZE)
      	strncpy(out_name,symbuf,Len);
      else
      	strncpy(out_name,symbuf,FNAME_SIZE);
      fprintf(stdout,"New output modifier: %s\n",out_name);
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"d")==0 || \
	  strcasecmp(buffer[depth_level],"data")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Current file name: %s\n",data_name);
      fprintf(stdout,"Enter new name:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"No changes\n");
      fprintf(stdout,"Data file: %s\n",data_name);
    }
    else{
      if(Len < FNAME_SIZE)
      	strncpy(data_name,symbuf,Len);
      else
      	strncpy(data_name,symbuf,FNAME_SIZE);
      fprintf(stdout,"New data file: %s\n",data_name);
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"c")==0 || \
	  strcasecmp(buffer[depth_level],"convert")==0){
    if(strlen(buffer[depth_level+1])!=0)
      convert_data_file(data_name,1);
    else
      convert_data_file(data_name,1);
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"savin")==0) || \
	   (strcasecmp(buffer[depth_level],"si") == 0)) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Filename for initials:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    fname = fopen(symbuf,"w");
    for(i=0;i<DIM;i++)
      fprintf(fname,"%lf\n",x[i]);
    fclose(fname);
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"savdat")==0) || \
	   (strcasecmp(buffer[depth_level],"sd")==0) ){
    fprintf(stdout,"Obsolete for now. Maybe retrieved in the future.\n");
    /* if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){ */
    /*   fprintf(stdout,"Filename for data:\n");
	 symbuf = read_input_line(symbuf,Len);*/
    /* } */
    /* if(strlen(symbuf) == 0){ */
    /*   fprintf(stdout,"Nothing was written\n"); */
    /* }else{ */
    /*   fname = fopen(symbuf,"w"); */
    /*   for(j=0;j<write_count;j++){ */
    /* 	fprintf(fname,"%.5lf ",ts[j]); */
    /* 	for(i=0;i<DIM;i++){ */
    /* 	  fprintf(fname,"%.5lf ",xs[j][i]); */
    /* 	} */
    /* 	fprintf(fname,"\n"); */
    /*   } */
    /*   fclose(fname); */
    /* } */
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"ld")==0) || \
	   (strcasecmp(buffer[depth_level],"load")==0) ) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Data file to load:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read\n");
    }else{
      if(ma_span){/* Calculate MA if ma_span is non-zero */
	if((*(symbuf+strlen(symbuf)-1) == 'a') && (*(symbuf+strlen(symbuf)-2) == 'm')
	   && (*(symbuf+strlen(symbuf)-3) == '.')){
	  /*The file name contanes .ma extension, perhaps already a MA trajectory*/
	  printf("Looks like it is already an MA data.\n");
	  printf("Do you still want to calculate MA?[Y/n]\n");
	  if(getchar() == 'n'){
	    printf("Interrupting...\n");
	    return 0;
	  }
	}
	moving_average(symbuf,ma_span);
	if((strlen(symbuf)+4) < Len)
	  symbuf = strcat(symbuf,".ma");
	else{
	  symbuf = (char *)realloc(symbuf,(strlen(symbuf)+4)*sizeof(char));
	  symbuf = strcat(symbuf,".ma");
	}
      }
      if((analyze_traj(symbuf)) != 0)/* Recompute common entities */
	return 10;
      if(graph_flag){/* Plot the data */
	/* Remove .ma extension since gplot_results will add it */
	if(ma_span)
	  *(symbuf+strlen(symbuf)-3) = '\0';
	fname = fopen(symbuf,"r");
	if(fname == NULL){
	  fprintf(stderr,"Error(f ld): could not open file.\n");
	  return 101;
	}
	if((info=get_info_data(info,fname)) == NULL){
	  fprintf(stderr,"Error: could not get data info\n");
	}
	fclose(fname);
	if(info[0]){//method is complex
	  /* First, plot determ trajectory */
	  gplot_results(symbuf,1,"rkf45");//any determ method, not "discrete"
	  /* Plot stoch trajectory */
	  if(!info[3])//lang_flag
	    gplot_results(symbuf,2,"discrete");
	  else
	    gplot_results(symbuf,2,"rkf45");//any determ method, not "discrete"
	}
	else{//method is not complex
	  if((!info[3]) && (!info[2]))//lang_flag and pure_det are false
	    gplot_results(symbuf,0,"discrete");//...then discrete method
	  else
	    gplot_results(symbuf,0,"rkf45");//determ or langevin alone...
	}
	fprintf(stdout,"Graph is updated.\n");
      }
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"")==0)
    return 0;
  else{
    fprintf(stdout,"No such command\n");
    return 0;
  }
}

int cont_interp(char **buffer,const int depth_level)
{
  int i,j,k;
  int Len = 100;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  char *p;
  
  if(strcasecmp(buffer[depth_level], "exit") == 0 || strcasecmp(buffer[depth_level],"quit") == 0 || strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving CONTINUE menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) ||	\
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"CONTINUE menu:\n");
    fprintf(stdout,"(Sh)ow\n");
    fprintf(stdout,"(P)arameter\n");
    fprintf(stdout,"R(a)nge\n");
    fprintf(stdout,"Parameter (v)alues\n");
    fprintf(stdout,"(R)un\n");
    fprintf(stdout,"(R)un (r)andom init\n");
    fprintf(stdout,"Run on (d)ata\n");
    fprintf(stdout,"Number of periods to (g)et\n");
    fprintf(stdout,"(C)heck dynamics\n");
    
    /* fprintf(stdout,"(P)roblem (p)arameters\n"); */
    /* fprintf(stdout,"(I)nitial\n"); */
    /* fprintf(stdout,"(B)ounds\n"); */
    /* fprintf(stdout,"(S)tep\n"); */
    /* fprintf(stdout,"S(o)lution\n"); */
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"c") == 0){
    if(dyn_check_flag_cont){
      dyn_check_flag_cont = 0;
      fprintf(stdout,"Checking dynamics OFF\n");
    }
    else{
      dyn_check_flag_cont = 1;
      fprintf(stdout,"Checking dynamics ON\n");
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"g") == 0 ||	\
	  strcasecmp(buffer[depth_level],"get") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter number of periods needed for each parameter set:\n");
      fprintf(stdout,
	      "(I will try to estimate the period and adjust the simulation time.\n ");
      fprintf(stdout,"Enter 0 to switch the option off.)\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read.\n");
    }
    else{
      if(!isalpha(*symbuf))
	num_per_to_get = atoi(symbuf);
      else
	printf("Something went wrong.\n");
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"d") == 0 ||	\
	  strcasecmp(buffer[depth_level],"data") == 0){
    char *app = (char *)malloc(FNAME_SIZE*sizeof(char));
    char *base = (char *)malloc(Len*sizeof(char));
    if((check_next_arg(buffer,depth_level,base,Len) != 0)){
      fprintf(stdout,"Files basename:\n");
      fprintf(stdout,"(if none specified, data filename is used)\n");
      base = read_input_line(base,Len);
    }
    if(strlen(base) == 0){
      strncpy(base,data_name,Len-1);
    }
    fprintf(stdout,"Computing common entities...\n");
    fprintf(stdout,"Files basename: %s\n",base);
    for(i=0;i<npvalues;i++){
      mu.P[par1] = pvalues[i];
      fprintf(stdout,"|***** %s = %G *****|\n",par_name[par1],mu.P[par1]);
      strncpy(symbuf,base,Len-1);
      sprintf(app,"_%s%G",par_name[par1],mu.P[par1]);
      strncat(symbuf,app,strlen(app));
      if((analyze_traj(symbuf)) != 0)
	return 10;
    }
    free(app);
    free(base);
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"p") == 0) ||	\
	   (strcasecmp(buffer[depth_level],"par") == 0) ) {
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Parameter name or the index(1-based):\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read.\n");
    }
    else{
      if(isalpha(symbuf[0])){
	for(i=0;i<PARDIM;i++)
	  if(strcasecmp(symbuf,par_name[i])==0)
	    par1 = i;
      }
      else{/* Numeric index */
	par1 = atoi(symbuf) - 1;/* Make it 0-based */
	if(par1 >= PARDIM){
	  fprintf(stderr,"Error: out of range par index, max is %d\n",PARDIM);
	  par1 = 0;
	}
      }
      fprintf(stdout,"Parameter: %s\n",par_name[par1]);
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"range")==0 || \
	  strcasecmp(buffer[depth_level],"a")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter range:\n");
      fprintf(stdout,"(format \"initial:[step:]final\", default step=1)\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read.\n");
    }
    else{
      p = strchr(symbuf,':');
      if(p!=NULL)
	*p = 0;/* Terminate at ":" */
      P1in = atof(symbuf);
      symbuf = p+1;/* Move to the next position */
      if(strchr(symbuf,':') != NULL){/* we find a ':' sign */
	p = strchr(symbuf,':');
	if(p == symbuf){/* Double ':' sign, default step */
	  P1stepin = 1;
	}
	else{/* Step was mentioned */
	  *p = 0;
	  P1stepin = atof(symbuf);
	}
	symbuf = p+1;
      }
      else{/* Single ':' sign, default step */
	P1stepin = 1;
      }
      P1max = atof(symbuf);
      P1min = P1in;
    }
    npvalues = 0;
    while(P1in <= P1max){
      if(npvalues == 0)
	pvalues = (double *)malloc(sizeof(double));
      else
	pvalues = (double *)realloc(pvalues,(npvalues+1)*sizeof(double));
      pvalues[npvalues] = P1in;
      npvalues++;
      P1in = P1in + P1stepin;
    }
    P1in = P1min;
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"v") == 0 || \
	  strcasecmp(buffer[depth_level],"val") == 0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Enter values for the parameter to take:\n");
      fprintf(stdout,"(space is a separator, additional spaces ignored)\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read.\n");
    }
    else{
      npvalues = 0;/* number of pvalues */
      char *tofree, *token;
      tofree = symbuf;
      while((token = strsep(&symbuf," \t")) != NULL){
	if(strlen(token) != 0){
	  if(npvalues == 0)
	    pvalues = (double *)malloc((npvalues+1)*sizeof(double));
	  else
	    pvalues = (double *)realloc(pvalues,(npvalues+1)*sizeof(double));
	  pvalues[npvalues] = atof(token);
	  npvalues++;
	}
      }

      /* while((*symbuf) != '\0'){ */
      /* 	p = 0; */
      /* 	if(strchr(symbuf,' ') != NULL){ */
      /* 	  p = strchr(symbuf,' '); */
      /* 	  *p = 0; */
      /* 	} */
      /* 	if(npvalues == 0) */
      /* 	  pvalues = (double *)malloc((npvalues+1)*sizeof(double)); */
      /* 	else */
      /* 	  pvalues = (double *)realloc(pvalues,(npvalues+1)*sizeof(double)); */
      /* 	pvalues[npvalues] = atof(symbuf); */
      /* 	npvalues++; */
      /* 	if(p) */
      /* 	  symbuf = p+1; */
      /* 	else */
      /* 	  break; */
      /* 	while((*symbuf) == ' ') */
      /* 	  symbuf++; */
      /* } */
      fprintf(stdout,"Read %d values\n",npvalues);
    }
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"run") == 0) ||	\
	   (strcasecmp(buffer[depth_level],"r") == 0) ){
    /* <k> will hold the number of runs for each case */
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Number of runs for each case:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(strlen(symbuf) == 0){
      fprintf(stdout,"Nothing was read. Default is 1.\n");
      k = 1;
    }
    else{
      if(isdigit(symbuf[0]))
	k = atoi(symbuf);
      else{
	fprintf(stderr,"Wrong format! Left default.\n");
	k = 1;
      }
    }
    run_extend(par1,pvalues,npvalues,k);
    return 0;
  }
  else if( (strcasecmp(buffer[depth_level],"rr") == 0)){
    run_extend_rnd(par1,pvalues,npvalues);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"pp") == 0){
    // First problem parameter.
    k=0;
    fprintf(stdout,"Enter 1st problem parameter:\n");
    fgets(buf,10,stdin);
    for(j=0;j<10;j++)
      if((*(buf+j) == (char)'\n') || (*(buf+j) == (char)' ') || (*(buf+j) == (char)'\t'))
	*(buf+j) = (char)'\0';
    for(i=0; i<mynum.num_par; i++){
      if(strcasecmp(buf,par_name[i]) == 0){
	P1=mu.P[i]; par1=i;
	k+=1;}
    }
    if(k == 0)
      fprintf(stdout,"Error: no such parameter\n");
    if(k > 1)
      fprintf(stdout,"Many parameters for the same name?\n");
    // Second problem parameter.
    k=0;
    fprintf(stdout,"Enter 2nd problem parameter:\n");
    fgets(buf2,10,stdin);
    for(j=0;j<10;j++)
      if((*(buf2+j) == (char)'\n') || (*(buf2+j) == (char)' ') || (*(buf2+j) == (char)'\t'))
	*(buf2+j) = (char)'\0';
    for(i=0; i<mynum.num_par; i++){
      if(strcasecmp(buf2,par_name[i]) == 0){
	P2=mu.P[i]; par2=i;
	k+=1;}
    }
    if(k == 0)
      fprintf(stdout,"Error: no such parameter\n");
    if(k > 1)
      fprintf(stdout,"Many parameters for the same name?\n");
		
    return 0;
  }
    else if( strcasecmp(buffer[depth_level],"sh") == 0 ){
    fprintf(stdout,"Continuation parameters:\n");
    fprintf(stdout,"%s\n",par_name[par1]);
    /* fprintf(stdout,"Initials:\n"); */
    /* fprintf(stdout,"%G\t%G\n",P1in,P2in); */

    fprintf(stdout,"Parameter values:\n");
    if(npvalues < 20){
      for(i=0;i<npvalues;i++)
	fprintf(stdout,"%G ",pvalues[i]);
    }
    else{
      for(i=0;i<10;i++)
	fprintf(stdout,"%G ",pvalues[i]);
      fprintf(stdout," ...  ");
      for(i=npvalues-10;i<npvalues;i++)
	fprintf(stdout,"%G ",pvalues[i]);
    }
    fprintf(stdout,"\n");
    fprintf(stdout,"Number of periods to get:\n");
    if(!num_per_to_get)
      fprintf(stdout,"off\n");
    else
      fprintf(stdout,"%d\n",num_per_to_get);
    
    /* fprintf(stdout,"Parameters range:\n"); */
    /* fprintf(stdout,"[%G %G],[%G %G]\n",P1min,P1max,P2min,P2max); */

    /* fprintf(stdout,"Initial step for continuation:\n"); */
    /* fprintf(stdout,"%G  %G\n",P1stepin,P2stepin); */

    /* fprintf(stdout,"Problem solution:\n"); */
    /* if(cont_sol==0) */
    /*   fprintf(stdout,"(Is not specified)\n"); */
    /* else */
    /*   fprintf(stdout,"%d\n",cont_sol); */

    return 0;
  }
  else if( strcasecmp(buffer[depth_level],"i") == 0 ){
    fprintf(stdout,"Enter starting values for parameters:\n");
    fprintf(stdout,"Par. #1(%s):\n",par_name[par1]);
    fgets(buf,20,stdin);
    if(strlen(buf) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P1in=atof(buf);
    }
    //fscanf(stdin,"%lf",&P1);
    fprintf(stdout,"Par. #2(%s):\n",par_name[par2]);
    fgets(buf2,20,stdin);
    if(strlen(buf2) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P2in=atof(buf2);
    }
    //fscanf(stdin,"%lf",&P2);

    return 0;
  }
  else if( strcasecmp(buffer[depth_level],"b") == 0 ) {
    // par #1
    fprintf(stdout,"Bounds for par. #1(%s)(newline are delimiters of min and max values):\n",
	    par_name[par1]);
    fgets(buf,20,stdin);
    if(strlen(buf) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P1min=atof(buf);}
    fgets(buf2,20,stdin);
    if(strlen(buf2) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P1max=atof(buf2);}
    //checking if P1min < P1max
    if(P1min > P1max){
      double tmp=P1max;
      P1max=P1min;
      P1min=tmp;
    }
			
    //par #2
    fprintf(stdout,"Bounds for par. #2(%s)(newline are delimiters of min and max values):\n",
	    par_name[par2]);
    fgets(buf,20,stdin);
    if(strlen(buf) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P2min=atof(buf);}
    fgets(buf2,20,stdin);
    if(strlen(buf2) == 1){
      fprintf(stdout,"No changes made.\n");
    }else{
      P2max=atof(buf2);}
    //checking if P2min < P2max
    if(P2min > P2max){
      double tmp=P2max;
      P2max=P2min;
      P2min=tmp;
    }

    return 0;
  }
  else if( strcasecmp(buffer[depth_level],"s") == 0 ){
    fprintf(stdout,"Initial step for continuation.\n");

    fprintf(stdout,"Par. #1(%s):\n",par_name[par1]);
    fgets(buf,20,stdin);
    if(strlen(buf) == 1){
      fprintf(stdout,"Left unchanged.\n");
    }else{
      P1stepin=atof(buf);
    }

    fprintf(stdout,"Par. #2(%s):\n",par_name[par2]);
    fgets(buf,20,stdin);
    if(strlen(buf) == 1){
      fprintf(stdout,"Leaven unchanged.\n");
    }else{
      P2stepin=atof(buf);
    }

    return 0;
  }
  else if( strcasecmp(buffer[depth_level],"o") == 0 ){
    fprintf(stdout,"Currently not supported, while changing whole trajectory system\n");

    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"")==0)
    return 0;
  else {
    fprintf(stdout,"No such command\n");
    return 0;
  }
}

int rand_interp(char **buffer,const int depth_level)
{
  int i,j,k,l,m;
  FILE *out,*in;
  double **frame;
  double *rand;
  double noise_ampl[DIM];
  char *symbuf;
  char *p,*p1;
  regStat *stat;
  if(strcasecmp(buffer[depth_level], "exit") == 0 || strcasecmp(buffer[depth_level],"quit") == 0 || strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving RANDOM menu...\n");
    return 1000;
  }
  else if(strcasecmp(buffer[depth_level],"ls") == 0) {
    fprintf(stdout,"RANDOM menu:\n");
    fprintf(stdout,"(Sh)ow summary\n");
    fprintf(stdout,"General:\n");
    fprintf(stdout,"\t(D)istribution\n");
    fprintf(stdout,"\t(A)lgorithm\n");
    fprintf(stdout,"\tS(e)ed\n");
    fprintf(stdout,"\t(S)eries\n");
    fprintf(stdout,"Random throwing:\n");
    fprintf(stdout,"\t(R)un\n");
    fprintf(stdout,"\t(P)oints\n");
    fprintf(stdout,"\t(C)heck dynamics\n");
    fprintf(stdout,"\t(Rad)ius\n");
    fprintf(stdout,"\tCo(n)straints\n");
    fprintf(stdout,"\tResults (f)ile\n");
    fprintf(stdout,"\t(F)ull trajectories to output\n");


    return 0;
  }
  else if(!strcasecmp(buffer[depth_level],"f")){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      if(only_init_flag){
	only_init_flag = 0;
      }
      else{
	only_init_flag = 1;
      }
    }
    else{
      only_init_flag = atoi(symbuf);
      /* Invert the flag */
      if(only_init_flag)
	only_init_flag = 0;
      else
	only_init_flag = 1;
    }
    if(only_init_flag)
      printf("Only initial conditions to output.\n");
    else
      printf("Full trajectories output.\n");
    
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"sh") == 0){
    fprintf(stdout,"***** General *****\n");
    fprintf(stdout,"Algorithm: %s\n",rng_type_env_val);
    fprintf(stdout,"Seed: %ld",rng_seed);
    if(seed_flag==1) fprintf(stdout,"(manual)\n");
    if(seed_flag==2) fprintf(stdout,"(automatic)\n");
    fprintf(stdout,"Amplitudes (Langevin amendments):\n");
    lang_amend(t,x,noise_ampl,mu.P);
    for(i=0;i<DIM;i++)
      fprintf(stdout,"%s\t",var_name[i]);
    fprintf(stdout,"\n");
    for(i=0;i<DIM;i++)
      fprintf(stdout,"%G\t",noise_ampl[i]);
    fprintf(stdout,"\n");
    fprintf(stdout,"Distribution: %s: pars: ",distribution_type);
    fprintf(stdout,"[");
    for(j=0;j<distribution_n_par;j++){
      if(j==distribution_n_par-1){
	fprintf(stdout,"%G]\n",distribution_par[j]);break;
      }
      fprintf(stdout,"%G,",distribution_par[j]);
    }
    fprintf(stdout,"***** Random throwing *****\n");
    fprintf(stdout,"Num. points:    %d\n",num_poin);
    fprintf(stdout,"Radius: %G\n",Rnd_Rad);
    fprintf(stdout,"Minimum: %G\n",Rnd_Low);
    fprintf(stdout,"Constraints:\n");
    for(j=0;j<DIM;j++){
      if(Rnd_Rad_Bnd_Idx[j]){
	printf("%s: (%G,%G)\t",var_name[j],Rnd_Rad_Bnd[j],Rnd_Rad_Bnd[j+DIM]);
      }
    }
    printf("\n");
    fprintf(stdout,"Check dynamics: ");
    if(dyn_check_flag)
      fprintf(stdout,"yes\n");
    else
      fprintf(stdout,"no\n");
    fprintf(stdout,"Output file: %s\n",rnd_throw_fname);
    fprintf(stdout,"Full trajectories to file: ");
    if(!only_init_flag)
      fprintf(stdout,"yes\n");
    else
      fprintf(stdout,"no (only initial points)\n");
    
    return 0;
  }
  else if (!strcasecmp(buffer[depth_level],"n") ||	\
	   !strcasecmp(buffer[depth_level],"con")){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"Type var name or index (1 = U(1), etc.) to change its constraints.\n");
      fprintf(stdout,"Multiple choices can be separated by SPACE.\n");
      fprintf(stdout,"Constraints for:\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }else {
      int var_idx_inquired[DIM];
      for(i=0;i<DIM;i++)/* Nullify the index */
	var_idx_inquired[i] = 0;
      p1 = symbuf;
      while(1){
	if((p = strchr(p1,' ')) != 0 || (p = strchr(p1,'\t')) != 0){
	  /* We can find space or tab, but not many of these in a row*/
	  *p = 0;/* then terminate at space */
	}
	if(strlen(p1) == 0){/* empty argument */
	  /* Just empty code: not to go into other conditions */
	}
	else if(isalpha(*p1)){
	  if(isVar(p1,var_name,DIM) > -1){
	    Rnd_Rad_Bnd_Idx[isVar(p1,var_name,DIM)] = 1;
	    var_idx_inquired[isVar(p1,var_name,DIM)] = 1;
	  }
	  else
	    fprintf(stdout,"Do not know the variable\n");
	}
	else{
	  if((atoi(p1) > DIM) || (atoi(p1) <= 0))
	    fprintf(stdout,"Error: out of the range\n");
	  else{
	    Rnd_Rad_Bnd_Idx[atoi(p1) - 1] = 1;
	    var_idx_inquired[atoi(p1) - 1] = 1;
	  }
	}
	if(p)/* p is not NULL */
	  p1 = p + 1;/* Move position right after the last found space */
	else /* p is NULL */
	  break;
      }
      /* For every non-zero Idx define the boudaries */
      free(symbuf);/* Free symbuf from the previous call */
      k = 0;/* the depth level increment */
      for(i=0;i<DIM;i++){
	if(var_idx_inquired[i]){
	  if((symbuf = read_next_command(buffer,depth_level+k+1)) == NULL){
	    fprintf(stdout,"Type the min and max boundaries for %s:\n",var_name[i]);
	    fprintf(stdout,"(Separate MIN and MAX by SPACE.)\n");
	    symbuf = read_command_line(NULL,&j);
	    k++;
	  }
	  if(symbuf == 0){
	    printf("Constraints are NOT set.\n");
	  }else {
	    p1 = symbuf;/* Make a copy of symbuf, to FREE symbuf later */
	    j = 0;
	    while(j < 2){
	      if((p = strchr(p1,' ')) != 0){/* We can find space */
		*p = 0;/* then terminate at space */
	      }
	      if(isalpha(*p1)){
		printf("Error: should be numerical value.\n");
	      }
	      else{
		if(j == 0)
		  Rnd_Rad_Bnd[i] = atof(p1);
		if(j == 1){
		  Rnd_Rad_Bnd[i+DIM] = atof(p1);
		  if(Rnd_Rad_Bnd[i] >= Rnd_Rad_Bnd[i+DIM]){
		    printf("Warning: MIN >= MAX. Commute.\n");
		    Rnd_Rad_Bnd[i+DIM] = Rnd_Rad_Bnd[i];
		    Rnd_Rad_Bnd[i] = atof(p1);
		  }
		}
	      }
	      if(p)
		p1 = p + 1;
	      else
		break;
	      j++;
	    }
	  }
	  printf("Constraints for %s: (%G,%G) \n",var_name[i],Rnd_Rad_Bnd[i],
		 Rnd_Rad_Bnd[i+DIM]);
	  free(symbuf);/* Free symbuf used for this reading */
	  /* fprintf(stdout,"Type the min and max boundaries for %s:\n",var_name[i]); */
	  /* fprintf(stdout,"(Separate MIN and MAX by SPACE.)\n"); */
	  /* symbuf = read_command_line(NULL,&j); */
	  /* p1 = symbuf;/\* Make a copy of symbuf, to FREE symbuf later *\/ */
	  /* if(symbuf == 0)/\* Null *\/ */
	  /*   printf("Constraints are NOT set.\n"); */
	  /* else{ */
	  /*   j = 0; */
	  /*   while(j < 2){ */
	  /*     if((p = strchr(p1,' ')) != 0){/\* We can find space *\/ */
	  /* 	*p = 0;/\* then terminate at space *\/ */
	  /*     } */
	  /*     if(isalpha(*p1)){ */
	  /* 	printf("Error: should be numerical value.\n"); */
	  /*     } */
	  /*     else{ */
	  /* 	if(j == 0) */
	  /* 	  Rnd_Rad_Bnd[i] = atof(p1); */
	  /* 	if(j == 1){ */
	  /* 	  Rnd_Rad_Bnd[i+DIM] = atof(p1); */
	  /* 	  if(Rnd_Rad_Bnd[i] >= Rnd_Rad_Bnd[i+DIM]){ */
	  /* 	    printf("Warning: MIN >= MAX. Commute.\n"); */
	  /* 	    Rnd_Rad_Bnd[i+DIM] = Rnd_Rad_Bnd[i]; */
	  /* 	    Rnd_Rad_Bnd[i] = atof(p1); */
	  /* 	  } */
	  /* 	} */
	  /*     } */
	  /*     if(p) */
	  /* 	p1 = p + 1; */
	  /*     else */
	  /* 	break; */
	  /*     j++; */
	  /*   } */
	  /* } */
	  /* printf("Constraints for %s: (%G,%G) \n",var_name[i],Rnd_Rad_Bnd[i], */
	  /* 	 Rnd_Rad_Bnd[i+DIM]); */
	  /* free(symbuf);/\* Free symbuf used for this reading *\/ */
	}
      }
    }
    /* MUST RETURN FROM FUNCTION HERE, DUE TO free(symbuf) calls cannot be called
  twice on the same pointer */
    return 0;
  }
  else if(!strcasecmp(buffer[depth_level],"radius") || \
	  !strcasecmp(buffer[depth_level],"rad")){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"Radius = %G\n",Rnd_Rad);
      fprintf(stdout,"Enter new radius:\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }
    else{
      Rnd_Rad = atof(symbuf);
      for(i=0;i<DIM;i++){
	if(!Rnd_Rad_Bnd_Idx[i])
	  Rnd_Rad_Bnd[i+DIM] = Rnd_Rad;
      }
      fprintf(stdout,"Current value is %G\n",Rnd_Rad);
    }
    free(symbuf);
    return 0;
  }
  else if(!strcasecmp(buffer[depth_level],"file") ||	\
	  !strcasecmp(buffer[depth_level],"f")){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"Results file is %s\n",rnd_throw_fname);
      fprintf(stdout,"Enter new file name:\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }
    else{
      rnd_throw_fname = (char *)realloc(rnd_throw_fname,(strlen(symbuf)+1)*sizeof(char));
      rnd_throw_fname = strcpy(rnd_throw_fname,symbuf);
      fprintf(stdout,"Current file is %s\n",rnd_throw_fname);
    }
    free(symbuf);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"d") == 0){
    fprintf(stdout,"UNDER RECONSTRUCTION\n");
    return 0;
    fprintf(stdout,"Distribution: %s\n",distribution_type);
    fprintf(stdout,"# of pars: %d\n",distribution_n_par);
    fprintf(stdout,"=========\n");
    fprintf(stdout,"(1) Gaussian\n");
    fprintf(stdout,"(2) Exponential\n");
    fgets(buf,10,stdin);
    while(strlen(buf)==1){
      fprintf(stdout,"Cannot understand.\n");
      fgets(buf,10,stdin);
    }
    m=atoi(buf);
    while(m!=1 && m!=2){
      fprintf(stdout,"Cannot understand.\n");
      fgets(buf,10,stdin);
      m=atoi(buf);
    }
    if(m==1) {
      distribution_type="gauss";
      distribution_n_par=2;
      distribution_par=(double *)realloc(distribution_par,distribution_n_par*sizeof(double));
      distribution_par[0]=0;
      distribution_par[1]=1;
    }
    if(m==2) {
      printf("Sorry: it is not supported yet.\n");
      /* distribution_type="exp"; */
/*       distribution_n_par=1; */
/*       distribution_par=(double *)realloc(distribution_par,sizeof(double)); */
/*       distribution_par[0]=1; */
    }
    
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"a") == 0){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"(1) MT19937\n");
      fprintf(stdout,"(2) Combined multiple recursive generator\n");
      fprintf(stdout,"(3) The fifth-order multiple recursive generator\n");
      fprintf(stdout,"(4) Tausworthe generator 1\n");
      fprintf(stdout,"(5) Tausworthe generator 2\n");
      fprintf(stdout,"(6) RANLXS generator (default)\n");
      fprintf(stdout,"Choose one...\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }else{
      k = atoi(symbuf);
      if(k <= 0){
	fprintf(stdout,"Error: must be positive.\n");
      }
      if(k > 6){
	fprintf(stdout,"Error: more than we have here.\n");
      }
      if(k == 1)
	rng_type_env_val = "mt19937";
      if(k == 2)
	rng_type_env_val = "cmrg";
      if(k == 3)
	rng_type_env_val = "mrg";
      if(k == 4)
	rng_type_env_val = "taus";
      if(k == 5)
	rng_type_env_val = "taus2";
      if(k == 6)
	rng_type_env_val = "ranlxd2";
    }
    free(symbuf);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"p") == 0){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"Points: %d\n",num_poin);
      fprintf(stdout,"How many points?\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }else{
      num_poin = atoi(symbuf);
      if(num_poin < 0){
  	fprintf(stdout,"Error: must be positive. Sign change.\n");
  	num_poin = -num_poin;
      }
    }
    free(symbuf);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"e") == 0){
    if((symbuf = read_next_command(buffer,depth_level)) == NULL){
      fprintf(stdout,"How to get seed:\n");
      fprintf(stdout,"(1) manually\n");
      fprintf(stdout,"(2) automatically\n");
      symbuf = read_command_line(NULL,&i);
    }
    if(symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }
    else{
      m = atoi(symbuf);
      if(m == 1) {
	seed_flag = 1;
	fprintf(stdout,"Enter your seed:\n");
	i = 0;
	while((l = getchar())!='\n'){
	  rng_seed_env_val[i]=(char)l;
	  i++;
	  if(i == SEED_SIZE){
	    fprintf(stdout,"Cannot take more than %d symbols\n",(int)SEED_SIZE);
	    break;
	  }
	}
	rng_seed_env_val[i]='\0';
	rng_seed = atoi(rng_seed_env_val);
	memset(rng_seed_env_val,'\0',SEED_SIZE);
      }
      if(m==2) seed_flag=2;
    }
    free(symbuf);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"s") == 0){
    //Generating SEED for random generator.
    if(seed_flag==2)
      rng_seed = generate_seed();
    rng = init_rng(rng_seed,rng_type,rng);
    //=====
    random_dist(num_poin*DIM,rng_seed);
    for(j=0;j<num_poin*DIM;j++)
      printf("%lf\n",u[j]);
    printf("=====\n");
    printf("Total: %d\nDIM: %d\nPer variable: %d.\n",num_poin*DIM,DIM,num_poin);

    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"c") == 0){
    if(dyn_check_flag){
      dyn_check_flag = 0;
      fprintf(stdout,"Checking dynamics OFF\n");
    }
    else{
      dyn_check_flag = 1;
      fprintf(stdout,"Checking dynamics ON\n");
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"r") == 0  || \
	  strcasecmp(buffer[depth_level],"run")==0) {
    stat = rnd_init_cond_run();/* see random.c for definition */
    if(dyn_check_flag)
      free_regime_stat(stat);

    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"")==0)
    return 0;
  else {
    fprintf(stdout,"No such command\n");
    return 0;
  }
}

int errors_interp(char **buffer, const int depth_level)
{
  int Len = 20;
  char *symbuf = (char *)malloc(Len*sizeof(char));
  int k;
  if(strcasecmp(buffer[depth_level], "exit") == 0 || strcasecmp(buffer[depth_level],"quit") == 0 || strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving ERRORS menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"ERRORS menu:\n");
    fprintf(stdout,"(C)ross\n");
    fprintf(stdout,"(I)ntegrator\n");
    fprintf(stdout,"(P)eriod\n");
    fprintf(stdout,"P(e)ak\n");
    fprintf(stdout,"(A)mplitude\n");
    fprintf(stdout,"Pe(r)iod Ratio\n");
    fprintf(stdout,"Re(g)ime Ratio\n");
    fprintf(stdout,"(Sh)ow\n");
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"show")==0 || strcasecmp(buffer[depth_level],"sh")==0){
    fprintf(stdout,"cross:\n");
    fprintf(stdout,"\tAbs.:%G\n",eps_abs_tper);
    fprintf(stdout,"\tRel.:%G\n",eps_rel_tper);
    fprintf(stdout,"integrator:\n");
    fprintf(stdout,"\tAbs.:%G\n",eps_abs_int);
    fprintf(stdout,"\tRel.:%G\n",eps_rel_int);
    fprintf(stdout,"\tscaling:\n");
    fprintf(stdout,"\t\tvar.:%G\n",a_y);
    fprintf(stdout,"\t\tderiv.:%G\n",a_dydt);
    fprintf(stdout,"period:\n");
    fprintf(stdout,"\tAbs.:%G\n",eps_abs_per);
    fprintf(stdout,"\tRel.:%G\n",eps_rel_per);
    fprintf(stdout,"peak:\n");
    fprintf(stdout,"\tAbs.:%G\n",eps_abs_peak);
    fprintf(stdout,"\tRel.:%G\n",eps_rel_peak);
    fprintf(stdout,"amplitude:\n");
    fprintf(stdout,"\tAbs.:%G\n",eps_abs_am);
    fprintf(stdout,"\tRel.:%G\n",eps_rel_am);
    fprintf(stdout,"period ratio:%G\n",eps_per_ratio);
    fprintf(stdout,"regime ratio:%G\n",eps_inregime_ratio);
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"cross")==0 ||
	  strcasecmp(buffer[depth_level],"c")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Cross errors:\n");
      fprintf(stdout,"(1) Abs.\n");
      fprintf(stdout,"(2) Rel.\n");
      symbuf = read_input_line(symbuf,Len);
      *(symbuf+1) = 0;/* Take just 1 symbol */
    }
    k = atoi(symbuf);
    if(k == 1){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Abs. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_abs_tper = atof(symbuf);
    }
    else if(k == 2){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Rel. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_rel_tper = atof(symbuf);      
    }
    else {
      fprintf(stdout,"No changes.\n");
    }
    return 0;
  }
  else if(strcasecmp(buffer[depth_level],"integrator")==0 || \
	  strcasecmp(buffer[depth_level],"int")==0 || \
	  strcasecmp(buffer[depth_level],"i")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Integrator errors:\n");
      fprintf(stdout,"(1) Abs.\n");
      fprintf(stdout,"(2) Rel.\n");
      fprintf(stdout,"Scaling factors:\n");
      fprintf(stdout,"(3) Variable\n");
      fprintf(stdout,"(4) Derivative\n");
      symbuf = read_input_line(symbuf,Len);
      *(symbuf+1) = 0;/* Take just 1 symbol */
    }
    k = atoi(symbuf);
    if(k == 1){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Abs. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_abs_int=atof(symbuf);
    }
    else if(k == 2){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Rel. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_rel_int=atof(symbuf);
    }
    else if(k == 3){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Var. scaling:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      a_y=atof(symbuf);
    }
    else if(k == 4){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Deriv. scaling:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      a_dydt=atof(symbuf);
    }
    else {
      fprintf(stdout,"No changes.\n");
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"period")==0 || \
	   strcasecmp(buffer[depth_level],"per")==0 || \
	   strcasecmp(buffer[depth_level],"p")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Period errors:\n");
      fprintf(stdout,"(1) Abs.\n");
      fprintf(stdout,"(2) Rel.\n");
      symbuf = read_input_line(symbuf,Len);
      *(symbuf+1) = 0;/* Take just 1 symbol */
    }
    k = atoi(symbuf);
    if(k == 1){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Abs. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_abs_per=atof(symbuf);
    }
    else if(k == 2){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Rel. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_rel_per=atof(symbuf);
    }
    else {
      fprintf(stdout,"No changes.\n");
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"peak")==0 || \
	   strcasecmp(buffer[depth_level],"e")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Peak errors:\n");
      fprintf(stdout,"(1) Abs.\n");
      fprintf(stdout,"(2) Rel.\n");
      symbuf = read_input_line(symbuf,Len);
      *(symbuf+1) = 0;/* Take just 1 symbol */
    }
    k = atoi(symbuf);
    if(k == 1){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Abs. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_abs_peak = atof(symbuf);
    }
    else if(k == 2){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Rel. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_rel_peak = atof(symbuf);
    }
    else {
      fprintf(stdout,"No changes.\n");
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"amplitude")==0 || \
	   strcasecmp(buffer[depth_level],"amp")==0 || \
	   strcasecmp(buffer[depth_level],"am")==0 || \
	   strcasecmp(buffer[depth_level],"a")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Amplitude errors:\n");
      fprintf(stdout,"(1) Abs.\n");
      fprintf(stdout,"(2) Rel.\n");
      symbuf = read_input_line(symbuf,Len);
      *(symbuf+1) = 0;/* Take just 1 symbol */
    }
    k = atoi(symbuf);
    if(k == 1){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Abs. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_abs_am = atof(symbuf);
    }
    else if(k == 2){
      if((check_next_arg(buffer,depth_level+1,symbuf,Len) != 0)){
	fprintf(stdout,"Rel. error:\n");
	symbuf = read_input_line(symbuf,Len);
      }
      eps_rel_am = atof(symbuf);
    }
    else {
      fprintf(stdout,"No changes.\n");
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"perratio")==0 || \
	   strcasecmp(buffer[depth_level],"pratio")==0 || \
	   strcasecmp(buffer[depth_level],"r")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Period ratio:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No changes.\n");
    } else {
      eps_per_ratio = atof(symbuf);
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"regratio")==0 || \
	   strcasecmp(buffer[depth_level],"rratio")==0 || \
	   strcasecmp(buffer[depth_level],"g")==0){
    if((check_next_arg(buffer,depth_level,symbuf,Len) != 0)){
      fprintf(stdout,"Regime ratio:\n");
      symbuf = read_input_line(symbuf,Len);
    }
    if(*symbuf == 0){
      fprintf(stdout,"No changes.\n");
    } else {
      eps_inregime_ratio = atof(symbuf);
    }
    return 0;
  }
  else if (strcasecmp(buffer[depth_level],"")==0)
    return 0;
  else {
    fprintf(stdout,"No such command\n");
    return 0;
  }
}

int sing_interp(char *buffer)
{
  int i;
  if(strcasecmp(buffer, "exit") == 0 || strcasecmp(buffer,"quit") == 0 || strcasecmp(buffer,"q") == 0) {
    printf("Leaving SINGULARITY menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer,"ls") == 0) || (strcasecmp(buffer,"l") == 0) ){
    fprintf(stdout,"SINGULARITY menu:\n");
    fprintf(stdout,"(R)un\n");
    fprintf(stdout,"(F)ile\n");
    fprintf(stdout,"(R)andom (f)lag\n");
    fprintf(stdout,"(F)ile (f)lag\n");
    fprintf(stdout,"(P)oints\n");
    fprintf(stdout,"(Up)per\n");
    fprintf(stdout,"(Lo)wer\n");
    fprintf(stdout,"(Sh)ow\n");
    fprintf(stdout,"***Technical parameters***\n");
    fprintf(stdout,"(I)terations\n");
    fprintf(stdout,"(S)olutions\n");
    fprintf(stdout,"(T)est\n");
    fprintf(stdout,"(E)rror\n");
    return 0;
  }
  else if(strcasecmp(buffer,"sh")==0 || strcasecmp(buffer,"show")==0){
    fprintf(stdout,"File: %s\n",ss_name);
    if(RND_FLAG)
      fprintf(stdout,"Random flag: yes\n");
    else
      fprintf(stdout,"Random flag: no\n");
    if(FILE_FLAG)
      fprintf(stdout,"File flag: yes\n");
    else
      fprintf(stdout,"File flag: no\n");
    fprintf(stdout,"Points: %d\n",PTS);
    fprintf(stdout,"Upper bound: %G\n",MAX_SNG);
    fprintf(stdout,"Lower bound: %G\n",MIN_SNG);
    fprintf(stdout,"***Technical parameters***\n");
    fprintf(stdout,"Iterations: %d\n",MAX_ITER);
    fprintf(stdout,"Roots: %d\n",MAX_N_SOL);
    fprintf(stdout,"Testing error: %G\n",TEST);
    fprintf(stdout,"Solution error: %G\n",SOL_ERROR);
    return 0;
  }
  else if (strcasecmp(buffer,"run")==0 || strcasecmp(buffer,"r")==0){
    multiroot_find();

    return 0;
  }
  else if (strcasecmp(buffer,"file")==0 || strcasecmp(buffer,"f")==0){
    fprintf(stdout,"Current filename: %s\n",ss_name);
    fprintf(stdout,"New filename:\n");
    memset(buf,'\0',100);
    fgets(buf,20,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<20){
      if(ss_name == NULL)
	ss_name=(char *)calloc(strlen(buf),sizeof(char));
      else{
	free(ss_name);
	ss_name=(char *)calloc(strlen(buf),sizeof(char));
      }
    }
    else{
      fprintf(stdout,"No new name OR name is more than 20 characters long.\n");}
    for(i=0;i<strlen(buf);i++)
      *(ss_name+i)=*(buf+i);

    return 0;
  }
  else if (strcasecmp(buffer,"rf")==0){
    if(RND_FLAG)
      fprintf(stdout,"Current random flag: yes\n");
    else
      fprintf(stdout,"Current random flag: no\n");
    fprintf(stdout,"New choice:(yes or no, 0 or not 0)\n");
    memset(buf,'\0',100);
    fgets(buf,5,stdin);
    buf_check(buf);
    if(strlen(buf) && isalpha(buf[0])){
      if(strcasecmp(buf,"yes")==0)
	RND_FLAG=1;
      else if(strcasecmp(buf,"no")==0)
	RND_FLAG=0;
      else
	fprintf(stdout,"Did not get the word. No changes.\n");
    }
    else if(strlen(buf) && !isalpha(buf[0])){
      RND_FLAG=atoi(buf);
    }
    else
      fprintf(stdout,"No changes.\n");
    
    return 0;
  }
  else if (strcasecmp(buffer,"ff")==0){
    if(FILE_FLAG)
      fprintf(stdout,"Current file flag: yes\n");
    else
      fprintf(stdout,"Current file flag: no\n");
    fprintf(stdout,"New choice:(yes or no, 0 or not 0)\n");
    memset(buf,'\0',100);
    fgets(buf,5,stdin);
    buf_check(buf);
    if(strlen(buf) && isalpha(buf[0])){
      if(strcasecmp(buf,"yes")==0){
	if(ss_name == NULL)
	  fprintf(stdout,"Tell me what filename is.\n");
	FILE_FLAG=1;
      }
      else if(strcasecmp(buf,"no")==0)
	FILE_FLAG=0;
      else
	fprintf(stdout,"Did not get the word. No changes.\n");
    }
    else if(strlen(buf) && !isalpha(buf[0])){
      FILE_FLAG=atoi(buf);
      if(FILE_FLAG && ss_name == NULL)
	fprintf(stdout,"Tell me what filename is.\n");
    }
    else
      fprintf(stdout,"No changes.\n");
    
    return 0;
  }
  else if (strcasecmp(buffer,"p")==0 || strcasecmp(buffer,"point")==0\
	   || strcasecmp(buffer,"points")==0){
    fprintf(stdout,"Current number of points: %d\n",PTS);
    fprintf(stdout,"(Remember: it works only if random flag is YES)\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      PTS=atoi(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"up")==0 ||\
	   strcasecmp(buffer,"upper")==0){
    fprintf(stdout,"Current max value: %G\n",MAX_SNG);
    fprintf(stdout,"(Remember: it works only if random flag is YES)\n");
    fprintf(stdout,"New max value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      MAX_SNG=atof(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"lo")==0 ||\
	   strcasecmp(buffer,"lower")==0){
    fprintf(stdout,"Current min value: %G\n",MIN_SNG);
    fprintf(stdout,"(Remember: it works only if random flag is YES)\n");
    fprintf(stdout,"New min value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      MIN_SNG=atof(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"i")==0 || strcasecmp(buffer,"iter")==0){
    fprintf(stdout,"Max number of iterations for each set of \
initials.\n");
    fprintf(stdout,"Current value: %d\n",MAX_ITER);
    fprintf(stdout,"New value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      MAX_ITER=atoi(buf);
    else
      fprintf(stdout,"No changes.\n");
    
    return 0;
  }
  else if (strcasecmp(buffer,"s")==0 || strcasecmp(buffer,"sol")==0){
    fprintf(stdout,"Max number of roots to find.\n");
    fprintf(stdout,"Current value: %d\n",MAX_N_SOL);
    fprintf(stdout,"New value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      MAX_N_SOL=atoi(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"t")==0 || strcasecmp(buffer,"test")==0){
    fprintf(stdout,"User-specified precision for successful iteration\
 process.\n");
    fprintf(stdout,"Current value: %G\n",TEST);
    fprintf(stdout,"New value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      TEST=atof(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"e")==0 || strcasecmp(buffer,"er")==0){
    fprintf(stdout,"Error to distinguish different solutions.\n");
    fprintf(stdout,"Current value: %G\n",SOL_ERROR);
    fprintf(stdout,"New value:\n");
    memset(buf,'\0',100);
    fgets(buf,15,stdin);
    buf_check(buf);
    if(strlen(buf) && strlen(buf)<15)
      SOL_ERROR=atof(buf);
    else
      fprintf(stdout,"No changes.\n");

    return 0;
  }
  else if (strcasecmp(buffer,"")==0)
    return 0;
  else {
    fprintf(stdout,"No such command\n");
    return 0;
  }
}

int lyap_interp(char **buffer,const int depth_level)
{
  int i;
  double *tmp;
  char *symbuf;
  if(strcasecmp(buffer[depth_level], "exit") == 0 || \
     strcasecmp(buffer[depth_level],"quit") == 0 || \
     strcasecmp(buffer[depth_level],"q") == 0) {
    printf("Leaving LYAPUNOV menu...\n");
    return 1000;
  }
  else if( (strcasecmp(buffer[depth_level],"ls") == 0) || \
	   (strcasecmp(buffer[depth_level],"l") == 0) ){
    fprintf(stdout,"LYAPUNOV menu:\n");
    fprintf(stdout,"(R)un\n");
    fprintf(stdout,"(D)istance threshold\n");
    fprintf(stdout,"(I)nitial distance\n");
    fprintf(stdout,"(R)un (c)ontinuation\n");
    fprintf(stdout,"(P)arameter\n");
    fprintf(stdout,"Parameter ma(x)\n");
    fprintf(stdout,"Parameter mi(n)\n");
    fprintf(stdout,"Parameter (s)tep\n");
    fprintf(stdout,"(O)utput file for MLE analysis\n");
    fprintf(stdout,"(R)un (s)pectrum\n");
    fprintf(stdout,"Run (s)pectrum (c)ontinuation\n");
    fprintf(stdout,"Output (f)ile for LCE analysis\n");
  }
  else if ( (strcasecmp(buffer[depth_level],"sh")==0) ){
    fprintf(stdout,"Threshold value (DTR):\t %G\n",DTR);
    fprintf(stdout,"Initial distance (DI):\t %G\n",DI);
    fprintf(stdout,"Parameter:\t %s=%G\n",par_name[LPI],mu.P[LPI]);
    fprintf(stdout,"Parameter max:\t %G\n",LPARMAX);
    fprintf(stdout,"Parameter min:\t %G\n",LPARMIN);
    fprintf(stdout,"Parameter step:\t %G\n",LDP);
    fprintf(stdout,"Compute time:\t %G(change from numerics)\n",mynum.total_time);
    fprintf(stdout,"MLE continuation output: %s\n",mle_file);
    fprintf(stdout,"LCE spectrum continuation output: %s\n",lyap_spec_file);
  }
  else if ( (strcasecmp(buffer[depth_level],"r")==0) || \
	    (strcasecmp(buffer[depth_level],"run")==0) ){
    tmp = (double *)malloc(sizeof(double));
    tmp = mle(tmp);
    free(tmp);
  }
  else if ( (strcasecmp(buffer[depth_level],"d")==0) ){
    read_next_command_wrapper(buffer,depth_level,"threshold (DTR)",(void *)&DTR,1);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Enter threshold (DTR) value:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes.\n"); */
    /* } */
    /* else{ */
    /*   DTR = atof(symbuf); */
    /*   fprintf(stdout,"New DTR value is %G\n",DTR); */
    /*   free(symbuf); */
    /* } */
  }
  else if ( (strcasecmp(buffer[depth_level],"i")==0) ){
    read_next_command_wrapper(buffer,depth_level,"initial distance (DI)",(void *)&DI,1);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Enter initial distance (DI) value:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes.\n"); */
    /* } */
    /* else{ */
    /*   DI = atof(symbuf); */
    /*   fprintf(stdout,"New DI value is %G\n",DI); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"rc")==0){
    mle_cont();
  }
  else if (strcasecmp(buffer[depth_level],"p")==0){
    symbuf = (char *)malloc(100*sizeof(char));
    read_next_command_wrapper(buffer,depth_level,"parameter name",(void *)symbuf,2);
    if(isPar(symbuf,PARDIM) > -1)
      LPI = isPar(symbuf,PARDIM);
    else{
      fprintf(stderr,"Error: no parameter with name %s\n",symbuf);
      fprintf(stderr,"Parameter (%s) unchanged.\n",par_name[LPI]);
    }
    free(symbuf);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Enter parameter name:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes.\n"); */
    /* } */
    /* else{ */
    /*   for(i=0;i<PARDIM;i++){ */
    /* 	if(strcmp(symbuf,par_name[i])==0){ */
    /* 	  LPI = i; */
    /* 	  break; */
    /* 	} */
    /* 	if(i+1==PARDIM){ */
    /* 	  fprintf(stdout,"No parameter with this name (%s).\n",symbuf); */
    /* 	  return 1; */
    /* 	} */
    /*   } */
    /*   fprintf(stdout,"New parameter: %s\n",par_name[LPI]); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"x")==0){
    read_next_command_wrapper(buffer,depth_level,"max parameter",(void *)&LPARMAX,1);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Max parameter value:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes\n"); */
    /* }else{ */
    /*   LPARMAX = atof(symbuf); */
    /*   fprintf(stdout,"New par max value is %G\n",LPARMAX); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"n")==0){
    read_next_command_wrapper(buffer,depth_level,"min parameter",(void *)&LPARMIN,1);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Min parameter value:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes\n"); */
    /* }else{ */
    /*   LPARMIN = atof(symbuf); */
    /*   fprintf(stdout,"New par max value is %G\n",LPARMIN); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"s")==0){
    read_next_command_wrapper(buffer,depth_level,"parameter step",(void *)&LDP,1);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"Parameter step:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes\n"); */
    /* }else{ */
    /*   LDP = atof(symbuf); */
    /*   fprintf(stdout,"New par max value is %G\n",LDP); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"o")==0){
    read_next_command_wrapper(buffer,depth_level,"MLE continuation output",
			      (void *)mle_file,2);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"File name for MLE continuation output:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes\n"); */
    /* }else{ */
    /*   mle_file = strcpy(mle_file,symbuf); */
    /*   fprintf(stdout,"New file name %s\n",mle_file); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"rs")==0){
    tmp = (double *)malloc(DIM*sizeof(double));
    tmp = lyap_spec(tmp);
    free(tmp);
  }
  else if (strcasecmp(buffer[depth_level],"sc")==0){
    lyap_spec_cont();
  }
  else if (strcasecmp(buffer[depth_level],"f")==0){
    read_next_command_wrapper(buffer,depth_level,"LCE spectrum continuation output",
			      (void *)mle_file,2);
    /* if((symbuf = read_next_command(buffer,depth_level)) == NULL){ */
    /*   fprintf(stdout,"File name for LCE spectrum continuation output:\n"); */
    /*   symbuf = read_command_line(NULL); */
    /* } */
    /* if(symbuf == 0){ */
    /*   fprintf(stdout,"No changes\n"); */
    /* }else{ */
    /*   lyap_spec_file = strcpy(lyap_spec_file,symbuf); */
    /*   fprintf(stdout,"New file name %s\n",lyap_spec_file); */
    /*   free(symbuf); */
    /* } */
  }
  else if (strcasecmp(buffer[depth_level],"")==0) {}
  else {
    fprintf(stdout,"No such command\n");
  }
  return 0;
}

int get_parameter(char **buffer, const int depth_level){
  int i,j;
  int Len = 20;
  char * symbuf = (char *)malloc(Len*sizeof(char));
  /* Call with an argument */
  if(*(buffer[depth_level])){
    if(Len > MAX_ARG_LEN)
      strncpy(symbuf,buffer[depth_level],MAX_ARG_LEN);
    else
      strncpy(symbuf,buffer[depth_level],Len);
    for(i=0;i<PARDIM;i++){
      if(strcmp(symbuf,par_name[i]) == 0){
	if(check_next_arg(buffer,depth_level,symbuf,Len)){
	  fprintf(stdout,"value:\n");
	  symbuf = read_input_line(symbuf,Len);
	}
	/*Check for zero input*/
	if(*symbuf == 0){
	  fprintf(stdout,"No changes.\n");
	  break;
	}else{
	  if(isalpha(*symbuf)){
	    printf("Character! Left unchanged.\n");
	    break;
	  }
	  mu.P[i] = atof(symbuf);
	  for(j=0;j<DIM;j++){/* Check par-xin link */
	    if(xin_par[j] == i){
	      x[j] = mu.P[i];/* Set current var to the value */
	      y[j] = (long)mu.P[i];/* Set current var to the value */
	    }
	  }
	  break;
	}
      }
      if(i+1 == PARDIM)
	fprintf(stdout,"\nNo parameter with such a name.\n");
    }
  }
  /* Print parameters */
  fprintf(stdout,"Parameters:\n");
  for(i=0;i<PARDIM;i++)
    fprintf(stdout,"%s = %G\n",par_name[i],p_mu->P[i]);
  /* Call without an argument */
  if(!(*buffer[depth_level])){
    fprintf(stdout,"Enter parameter name:\n");
    symbuf = read_input_line(symbuf,Len);
    if(*symbuf == 0){
      fprintf(stdout,"No changes.\n");
    }else{
      for(i=0; i<PARDIM; i++){
	if( strcmp(symbuf,par_name[i]) == 0 ){
	  fprintf(stdout,"value:\n");
	  symbuf = read_input_line(symbuf,Len);
	  /*Check for zero input*/
	  if(*symbuf == 0){
	    fprintf(stdout,"No changes.\n");
	    break;
	  }else{
	    if(isalpha(*symbuf)){
	      printf("Character! Left unchanged.\n");
	      break;
	    }
	    mu.P[i] = atof(symbuf);
	    for(j=0;j<DIM;j++){/* Check par-xin link */
	      if(xin_par[j] == i){
		x[j] = mu.P[i];/* Set current var to the value */
		y[j] = (long)mu.P[i];/* Set current var to the value */
	      }
	    }
	    break;
	  }
	}
	if(i+1==PARDIM)
	  fprintf(stdout,"\nNo parameter with such a name.\n");
      }
    }
  }
  return 0;
}

int load_initials(char *init_name){
  /* The function loads the init conditions and stores them in the system-wide
     variables xin[]. If there are less than DIM numbers in init_name file then the
     function fills the rest values with zeros and makes a warning message. If there
     are more than DIM values in the source, the function takes only DIM values and
     ignores the rest. The function returns non-zero value when failed to open the
     file, otherwise, the return value is zero, even when not enough and too many
     points were found.*/
  int i,j;
  double tmp;
  FILE *init = fopen(init_name,"r");
  if(init == NULL){
    fprintf(stderr,"Error: cannot read the file.\n");
    return 10;
  }
  else{
    i = 0;
    while(i < DIM){
      j = fscanf(init,"%lf",(xin+i));
      if(j == EOF){
	fprintf(stderr,"Warning: not enough points.\n");
	break;
      }
      i++;
    }
    for(j=i;j<DIM;j++)/* Setting to zeros the rest */
      xin[j] = 0.0;
  }
  /* Copy the xin to yin(discrete initial conditions) */
  for(i=0;i<DIM;i++)
    yin[i] = (long int)ceil(xin[i]);/* take the ceiling */
  
  fclose(init);
  return 0;
}

int check_next_arg(char **buffer,const int depth_level, char *next_arg,
		   const size_t Len){
  /* If the next argument in the cmd exists, the function copies it to the *next_arg
  (a temporary array) and returns 0(true). Otherwise, it returns positive value(1)
  and does not change the *next_arg */
  if(strlen(buffer[depth_level+1]) != 0){
    strncpy(next_arg,buffer[depth_level+1],Len-1);
    *(next_arg+Len-1) = 0;/* Terminate string */
    if(strlen(buffer[depth_level+1])>Len-1){
      printf("Too many symbols in `%s'.",buffer[depth_level+1]);
      printf(" I read: `%s'\n",next_arg);
    }
    return 0;
  }
  else{
    return 1;
  }
}

char *read_input_line(char *symbuf, const int Len)
{
  /* This function reads a single input from stdin, removes newline
     and returns what it has read */
  fgets(symbuf,Len,stdin);
  if(strchr(symbuf,'\n') != NULL){
    *(strchr(symbuf,'\n')) = 0;
  }
  else{
    fprintf(stderr,"Input is too long: >%d char, read: %s\n",Len-1,symbuf);
  }
  return symbuf;
}


char * read_next_command(char **buffer,const int depth_level){
  /* If the next argument(next to the current depth_level) in the cmd exists, the
  function copies it to the output and returns pointer to that. Otherwise, it returns
  NULL. */
  char *buff = (char *)malloc(MAX_SYM_READ*sizeof(char));
  if((strlen(buffer[depth_level+1]) != 0) &&\
     (strlen(buffer[depth_level+1]) < MAX_SYM_READ)){
    buff = strncpy(buff,buffer[depth_level+1],MAX_SYM_READ);
  }
  else{
    free(buff);
    return NULL;
  }
  return buff;
}

char *read_command_line(FILE* stream, int *status)
{
  /* This function reads a single input from stream or stdin, removes newline and
     returns what it has read. The buffer allocation is done within the function.
     The memory should be cleared elsewhere. status is a variable showing the status
     of the fulfilled commands.
     status:
     0 - success: command was read and buffer was allocated.
     1 - memory allocation failed.
     2 - EOF reached (useful when reading from real file stream)
     3 - empty buffer, just newline was there
     4 - too long input, could not copy.
     When status = 1,2,3, or 4 the return pointer is always NULL.
  */
  char *buff = (char *)malloc(MAX_SYM_READ*sizeof(char));
  if(buff == NULL){
    *status = 1;
    return NULL;
  }
  if(stream == NULL)
    buff = fgets(buff,MAX_SYM_READ,stdin);
  else
    buff = fgets(buff,MAX_SYM_READ,stream);
  if(buff == NULL){/* end-of-file is reached */
    *status = 2;
    return NULL;
  }
  if(strchr(buff,'\n') != NULL){/* Clear trailing newline */
    *(strchr(buff,'\n')) = 0;
    if(strlen(buff) == 0){/* empty buff means nothing, return NULL */
      free(buff);
      *status = 3;
      return NULL;
    }
  }
  else{
    fprintf(stderr,"Input is too long: >%d char, read: \n<%s>\n",
	    MAX_SYM_READ-1,buff);
    free(buff);
    *status = 4;
    return NULL;
  }
  *status = 0;
  return buff;
}

int convert_numeric_input(const char *input, void *out, const int type)
{/* This function reads any numerical input, including the conversion of the input
    from variables/parameters to the corresponding numbers. type
    stands for the type of variable the out is. type can be:
    0 - int
    1 - double
 */
  if(isdigit(input[0])){
    /* If the 1st letter is numeric, we do the conversion */
    if(type == 0)
      *(int *)out = (int)strtod(input,(char **)NULL);
    else if(type == 1)
      *(double *)out = strtod(input,(char **)NULL);
  }
  else if(isVar(input,var_name,DIM) != -1){
    /* If the input is a var name */
    if(type == 0)
      *(int *)out = (int)x[isVar(input,var_name,DIM)];
    if(type == 1)
      *(double *)out = x[isVar(input,var_name,DIM)];
    printf("Warning (convert_numeric_input): I took value");
    printf(" %G from the continuous variable x[%d].\n",
	   *(double *)out,isVar(input,var_name,DIM));
  }
  else if(isPar(input,PARDIM) != -1){
    /* If the input is a par name */
    if(type == 0)
      *(int *)out = (int)mu.P[isPar(input,PARDIM)];
    if(type == 1)
      *(double *)out = mu.P[isPar(input,PARDIM)];
  }
  else{
    /* Wrong input */
    fprintf(stderr,"(convert_numeric_input): Wrong input.\n");
    return 1;
  }

  return 0;
}

int read_next_command_wrapper(char **buffer,const int depth_level,
				const char *entity_name,
				void *out, const int entity_type)
{/* Read the next argument and process the value, all-in-one function. entity_type
    stands for the type of variable the entity is. entity_type can be:
    0 - int
    1 - double
    2 - char *
 */
  /* First read the command line in automatic fashion or by inviting with the prompt
   */
  int i;
  char *symbuf = (char *)malloc(MAX_SYM_READ*sizeof(char));
  if(symbuf == NULL){
    fprintf(stderr,"Error(read_next_command_wrapper): ");
    fprintf(stderr,"the symbuf cannot be allocated.\n");
    return 1;
  }
  if((symbuf = read_next_command(buffer,depth_level)) == NULL){
    fprintf(stdout,"Enter new value for %s:\n",entity_name);
    symbuf = read_command_line(NULL,&i);
  }
  /* Then check if the input was actually there */
  if(symbuf == 0){
    fprintf(stdout,"No changes applied.\n");
  }
  else{
    if(entity_type != 2){
      convert_numeric_input(symbuf,out,entity_type);
      if(entity_type == 0)
	fprintf(stdout,"New %s value is %d\n",entity_name,*(int *)out);
      else if(entity_type == 1)
	fprintf(stdout,"New %s value is %G\n",entity_name,*(double *)out);
    }
    else{/* Entity type is string. */
      strncpy((char *)out,symbuf,strlen(symbuf));
      fprintf(stdout,"New %s value is %s\n",entity_name,(char *)out);
    }
  }
  free(symbuf);
  return 0;
}
