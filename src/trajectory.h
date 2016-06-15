/*****************************************************************************************/
/* Copyright 2008, 2009,2010,2011,2012,2013,2014 Elias Potapov. */
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
   Tampere University of Technology, Department of Signal Processing
   Moscow, Russia / Tampere, Finland */
/****************************************************************************************/

typedef struct {
  double *t;/* time of the peak/trough */
  double *x;/* the variable value at the peak */
  int size;/* number of peaks */
  int *peak_true_table;/* true table for the peaks(vs troughs) */
} trajPeak;

typedef struct {
  double *ampl;/* the amplitude of the slope */
  double *t0;/* the beginning of the slope in time */
  double *base;/* the level from which the ampl is calculated, the lowest point of
		  the slope */
  int *ascend;/* boolean array determining whether slope is ascending(1) or descending(0) */
} trajSlope;

typedef struct {/* Report structure for the regimes */
  int regime;/* numerical value of the regime */
  double ph_shift;/* max phase shift for the system, as a fraction of period */
  double period;/* period */
  int homog;/* homogeneity flag: true(>0) if homogeneous solution */
  double b_gain;/* base gain */
  double a_gain;/* amplitude gain */
} regReport;

typedef struct{
  int N;/* total number of regimes analyzed */
  int ss;/* number of stable steady states */
  int os;/* number of oscillatory regimes */
  int *os_pd;/* periodicity array for OS regimes */
  int *os_up;/* unique periodicity values */
  int os_up_cn;/* counter on how many unique periodicity values */
  int hg;/* number of homogeneous regimes */
  int ih;/* number of inhomogeneous regimes */
  int ss_hg;/* number of homogeneous SS */
  int os_hg;/* number of homogeneous oscillatory */
  int *os_hg_tr;/* true array of homogeneous OS */
  int os_ih;/* number of inhomogeneous oscillatory */
  int ss_ih;/* number of inhomogeneous SS */
  int os_ip;/* number of in-phase OS */
  int *os_ip_tr;/* true array of in-phase OS */
  int os_op;/* number of out-of-phase oscillatory regimes */
  int *os_op_tr;/* true array of out-of-phase OS */
  int ud;/* number of undetermined regimes */
  int mx;/* number of mixed regimes */
} regStat;

int traj_verb;/* verbose level */
double eps_traj_abs, eps_traj_rel;/* Main absolute and relative tolerances for
				     T-system */
double eps_traj_rel_hg;/* The relative error for homogeneity test */
