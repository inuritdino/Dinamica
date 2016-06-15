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
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/
// max number of systems
#define MAX_TOTAL_DIM 10
// max number of ODEs in each system
#define MAX_LOCAL_DIM 15
// Total dimension
#define MAXDIM MAX_LOCAL_DIM*MAX_TOTAL_DIM
// max number of parameters.
#define	MAX_N_PAR 100
// max number of intersections with Poincare sec. Normally 100.
#define MAX_N_ISEC 100
// max number of peaks. Normally 100.
#define MAX_N_PEAKS 100
/* MAX NUMBER OF PERIODS DEFINES*/
#define MAX_N_PERIODS MAX_N_ISEC-1
#if MAX_N_PEAKS > MAX_N_ISEC
#undef MAX_N_PERIODS
#define MAX_N_PERIODS MAX_N_PEAKS-1
#endif
#if MAX_N_ISEC > MAX_N_PEAKS
#undef MAX_N_PERIODS
#define MAX_N_PERIODS MAX_N_ISEC-1
#endif
/* END OF MAX NUMBER OF PERIODS DEFINES*/
#define MAX_N_SUBPER 100
#define VAR_NAME_SIZE 10
#define PAR_NAME_SIZE 10
#define MAX_N_ARG 10 /* max num of arguments */
#define MAX_ARG_LEN 100 /* max argument length */
#define MAX_N_HIST_ENT 10 /* max num of history entries */
#define FNAME_SIZE 200
#define METH_NAME_LEN 21
/* Frame constants */
#define FRAME_LINE_LEN 2000
#define FRAME_N_LINES 100
//Header files
#include "config.h"
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_rng.h>
#include "trajectory.h"
#include "graph.h"
#include "gnuplot_i.h"
#include "thistogram.h"

int ret_num;
/*buf is a multipurpose char buffer, primarily for the command line
  processing.*/
char ch,buf[MAX_ARG_LEN],buf2[MAX_ARG_LEN];
char ***buf_hist;
char **cmd;
struct num {
  double total_time;
  double trans_time;
  double step;
  int total_dim;
  int local_dim;
  int num_par;
  int global_buffer;/*Current buffer size of storage arrays*/
  int write_step;
  double smp_frq;/*Sample frequency: How frequently to write gillespie outcomes*/
  char method[METH_NAME_LEN];
  char method2[METH_NAME_LEN];/*Deterministic Method for the complex method option*/
} mynum;// Components of NUMERICS menu.
char method[METH_NAME_LEN];/*Integration method. The mynum.method should be removed with time.*/
char method2[METH_NAME_LEN];/* If method=complex, method2 is determ. method to use */
typedef struct {
  double P[MAX_N_PAR];
} MU;
MU mu;//Parameters of the system.
struct num* p_mynum;
MU* p_mu;
int DIM;/*Total dimension of the system--mathematical*/
int LDIM;/*Local dimension of the system--physical(#cells,#Brusselators etc)*/
int PARDIM;/*Parametric dimension*/
//int SDIM;/* Stochastic dimension, i.e. number of variables in the stoch system */
int NREAC;/*Number of chemical reactions for stochastic discrete approach*/
int NLANG;/*Number of Langevin amendments, terms multiplied be Gaussian noise N(0,1)*/
int BUFFER;/*Default buffer size for storage arrays*/
int BUF_INCR;/*it is increment of the buffer, when it becomes large during integration */
int FNUM;/*Number of functions defined by the user*/
int AUXNUM;/*Number of auxillary entities*/
double *x,*f;/*Vars and RHS of ODEs.*/
long int *y;/*Discrete modeling variables.*/
double *pr;/*Propensities of chemical reactions*/
int **upvec;/*Update vector(Gill)*/
double *ksi;/*Amplitude of noise terms.*/
double *a;/*Auxillary values.*/
double **xs,*ts;
int *xin_par;//Array of parameter-initials links
double *xin;//Initials.
long int *yin;/*Discrete initials*/
double t;//TIME
char *conf_name,*ode_name,*init_name, *data_name, *out_name, *input_name;
FILE *input_stream;
char par_name[MAX_N_PAR][PAR_NAME_SIZE];//Names of parameters.
//char var_name[MAXDIM][VAR_NAME_SIZE];//Names of vars.
char **var_name;/*var names.*/
char **aux_name;/*auxillaries' names.*/
int n_steps;//Number of steps of integrator.
int write_count;//How many points are written to data file or storage
//double *max_x,*min_x;//Max and min of vars.
double *x_cross;//Poincare section.
double cross;//Poincare section, new style.
double cross_level[3];//cross = (max-min)/cross_level + min, up to 3 values
int n_isec[MAXDIM];//Number of intersections with Poincare level.
int n_peaks[MAXDIM];//Number of peaks of oscillating trajectory.
int n_cavities[MAXDIM];//Number of cavities of oscillating trajectory.
double t_peak[MAX_N_PEAKS][MAXDIM];//Time when peak arise.
double x_peak[MAX_N_PEAKS][MAXDIM];//Value of var of the peak.
double t_cavity[MAX_N_PEAKS][MAXDIM];//Time when cavity arise.
double x_cavity[MAX_N_PEAKS][MAXDIM];//Value of var of the cavity.
double tper[MAX_N_ISEC][MAXDIM];//Time when intersected with Poincare section.
float per[MAX_N_PERIODS][MAXDIM];//Periods.
int short per_method;//Method periods: Poincare(0) or autocorrelation(1).
double *perDet;//Deterministic periods
double *perStoch;//stochastic periods
int nPerDet;// Number of determ. periods computed out from nPer.
int nPerStoch;//Total number of periods.
int short perVarInd;//Variable which was used for period calculation.
int short per_var;//From which var to get period, it can be detected automatically.
double per_thresh;//Period threshold
//double *ampl;//Amplitude of oscillations.
int ma_span;//Span for the moving average trajectory

float per_aver[MAX_N_SUBPER][MAXDIM];//Averaged periods.
int n_subperiods[MAXDIM];//Subperiods if exist.
int spectrum_per[MAX_N_PERIODS-1][MAXDIM];//Spectra(see manual).
double ss_factor[MAXDIM];//Criterion for SS convergence
double per_ratio[MAXDIM];//Period ratio(criterion for L.C. convergence).
int short regime;//Regime detection.
int short reg[MAX_LOCAL_DIM];//Regime for each var. type.

float big_per[MAXDIM];//If subperiods exist.
int short traj_trans;/* How many times to do transition to attractor in traj() */
/* VARIOUS FLAGS */
int short not_get_flag;
int short write_flag;//If == 0 then no write to file, but storage
		     //arrays.If != 0 write to file.
int short lang_flag;//If Langevin approach is applicable(=1)
int short jac_flag;/* If jacobian is provided in the source file =1,
		      if numerical =0 */
/****************/
char *main_prompt;
char *numerics_prompt;
char *graphics_prompt;
char *cross_prompt;
char *periodics_prompt;
char *traj_prompt;
char *file_prompt;
char *cont_prompt;
char *rand_prompt;
char *errors_prompt;
char *periodics_hist_prompt;
char *gnuplot_prompt;
int short cont_sol;
char *sing_prompt;
/*
  Prototyping the functions
*/
// main.c
int max_min_x(double [],double [],double,double);
void copyleft();
// init.c
int init_command_line();
int init();
void din_close();
int isVar(char const *,char **,int);
int isPar(char const *, const int);
// input_interpreter.c
void copy_history();
int parse_command_line(char **,const char*,FILE *);
int read_menu(char **, const char *);
int main_interp(char **, const int);
int graphics_interp(char **, const int);
int numerics_interp(char **, const int);
int periods_interp(char **, const int);
int traj_interp(char **, const int);
int file_interp(char **, const int);
int cont_interp(char **,const int);
int rand_interp(char **,const int);
int lyap_interp(char **,const int);
int sing_interp(char *);
int errors_interp(char **, const int);
int get_parameter(char **, const int);
int load_initials(char *);
int check_next_arg(char **,const int, char *, const size_t);
char *read_next_command(char **,const int);
char *read_input_line(char *, const int);
char *read_command_line(FILE *, int *);
int convert_numeric_input(const char *input,void *out,const int);
int read_next_command_wrapper(char **,const int,const char *,void *,const int);
//ode.c
int func_odeiv(double ,const double *, double *, void *);
int jac(double, const double *, double *, double *, void *);
int jac_general(double, const double *, double *, double *, void *,int(*)());
int ode_lin(double,const double *,double *,const double *,double *,void *);
// integrator.c
int run(const char *,const double,const int short,const int short, int);
int run_lin();
int analyze_traj(const char *);
int moving_average(const char *,const int);
int write_data_file(FILE *);
void write_data_array();
int buffer_overfull_check();
char *app_ext(char const *,char const *,char *);
// cross.c
void write_autocorr_file(double *,const int,char *, const char *);
double * compute_autocorrelation(double **, const int, double *, const int);
void thist_init();
void thist_free();
void write_histPer_file(double **,const int,char *,const char *);
double **compute_hist_per(double **,double *,const int,const int);
size_t get_var_per(char *,const int);
void write_perStoch_file(double *,const int, char *,const char *);
double *compute_period(double *,int *, double, const int,
		       const char *,const int, const int);
int compute_peak_autocorr(const double *, const int);
int amplitude(double *,double,double);
trajPeak *peak_trough(double **,const int,const int,trajPeak *);
trajPeak *peak_trough1(double **,const int,const int,trajPeak *);
void free_peak(trajPeak *);
trajPeak *peak_trough2(double **,const int,const int,trajPeak *);
double *peak_ampl(trajPeak, double *,int *);
double *ampl_multi_frames(double *ampl, int *);
int peak_def(int *,double *, double *, int (*)());
int intersec(int *, double *,int (*)());
double *period_cross(double *, int *, const int, const double *, const double);
double *crossing(int *, double *, const double, double **,
	     const int, const int);
//int cross(double [], double *, double *, int short);
int period(float *);
int period_average(float *, float *);
// continue.c
int run_extend_rnd(const int,const double [],const int);
int run_extend(const int,const double [],const int,const int);
int cont_init();
int continuation();
// random.c
regStat *rnd_init_cond_run();
gsl_rng * init_rng(const long int, const gsl_rng_type *, gsl_rng *);
double * rgen_uni(const int,const long int);
int random_dist(const int,const long int);
unsigned long int generate_seed();
int random_init();
void random_free();
// read_config.c
int read_config(char *);/* deprecated */
int buf_check(char *);
//int buf_strip(char *,const int);
int read_conf(char *);
int set(char *,char *);
int set_var(FILE *);
int set_aux(FILE *);
int set_par(FILE *);
// trajectory.c
int traj_init();
int ss_stab(const double *);
double D(double, double, double );
//int traj();
//int st_st_detect();
int get_varInd(int, int []);
void free_regime_stat(regStat *);
void copy_regime_stat(regStat *,regStat *);
regStat *regime_stat(regReport [],const int,FILE *);
void statND_to_file(FILE*, regStat **, int);
int *perdcty_in_stat_array(regStat **, int, int *);
char *report_translate(regReport);
regReport *get_sys_dynamics(regReport *,const int);
int find_sys_lag(trajSlope *[], int [],int *,double *,double *);
void free_report(regReport *);
regReport *get_dynamics_nd(trajSlope *[],int [],int []);
trajSlope *get_dynamics_1d(char const *,int,trajSlope *,int *,int *);
int ss_by_frame(double **,int,int);
trajSlope *check_slope_ampl(trajSlope *,int *);
void free_slope(trajSlope *);
trajSlope *slope_ampl(const trajPeak,trajSlope *);
void gplot_slope_ampl(trajSlope *,int const);
int steady_state(double *);
int period_sort(double *);
int subper(int *, int *);
int attr_conv(double **, int const);
int get_index_trans(const double **,const int,const double);
int mol_dist(const int);
// solver.c
int solver(double *,double, int(*)());
int solver_lang(double *,double, int(*)(), double []);
int solver_lin(const double *,double *,double *, double, int(*)());
int gill_init(long int *, double *, double *);
// submenu.c
int periodics_hist_interp(char **, const int);
// singularity.c
int func_multiroot(const gsl_vector *, void *, gsl_vector *);
int multiroot_find();
int print_state(size_t, gsl_multiroot_fsolver *);
int multiroot_init();
// lyapunov.c
double * mle (double *);
int mle_cont();
int lyap_init();
void lyap_free();
double * lyap_spec(double *);
int lyap_spec_cont();
// file.c
double ** load_frame(FILE *, double **, int *);
void free_frame(double **);
void write_frame(FILE *, double **, const int);
int read_out_file(const char *);
int set_info_data(int *,FILE *);
int *get_info_data(int *,FILE *);
int convert_data_file(char const *, int short);
// graph.c
void gplot_frames(const char *,const int,const int,const int);
int gnuplot_interp();
void gplot_results(const char *, const int, const char *);
void graph_set_labels(const struct coordNet,const char *);
void send_to_eps(const char *);
void init_graph();
int init_plotter(const char *);
void close_plotter(const int);
void plot_coord(const int, struct coordNet *);
char *get_label(char *,const int);

// User defined .c file
int lang_amend(double,const double [],double [],double []);
int rhs(double,const double [],double [],double []);
int auxillary(double,double [],double [],double []);
int propensity(const double [],const double [], double []);
int jacobian(double, const double [],double [],double [],double []);
int update(long int *, const int);
  
