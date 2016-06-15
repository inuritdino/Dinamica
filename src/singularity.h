/*****************************************************************************************/
/* Copyright 2008, 2009 Elias Potapov. */
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

short int RND_FLAG; /* Whether we use random generator */
int PTS; /* how many points for random generator to generate */
char *ss_name; /* name of the output file */
int FILE_FLAG; /* whether we write output */
double SOL_ERROR; /* Error for comparing solutions */
double TEST;
double MAX_SNG; /* Upper bound for random dispersion */
double MIN_SNG; /* Lower bound for random dispersion */
int MAX_ITER; /*max num of iterations of each initial point provided*/
int MAX_N_SOL;  /* max num of solutions(roots) to find */


