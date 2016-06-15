/*****************************************************************************************/
/* Copyright 2008,2009,2010 Elias Potapov. */
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

// This file describes errors for integrating, comparison of periods and so on.

/**************************
 * For trajectory identification(also peaks comparison errors used for
 * that, see below)
 **************************/
double eps_abs_tper,eps_rel_tper;

/*******************************
 * Errors for integration
 *******************************/
double eps_abs_int, eps_rel_int;
/* Next are the scaling factors for Gnu Scientific Library ODE control
   functions(read GSL reference)*/
double a_y, a_dydt;

/*****************************
 * Period comparison errors
 * **************************/
double eps_abs_per, eps_rel_per;

/*****************************
 * Peaks comparison
 * **************************/
double eps_abs_peak, eps_rel_peak;

/*****************************
 * Amplitudes comparison. Detection of IHLC and SS(special function
 * for SS is st_st_detect()).
 **************************/
double eps_abs_am, eps_rel_am;

/****************************
 * Period ratio comparison
 * *************************/
double eps_per_ratio;

/*****************************
 * Regime ratio comparison. Error for inphase regime.
 * **************************************************/
double eps_inregime_ratio;

