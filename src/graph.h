/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011,2012 Elias Potapov. */
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

#include "gnuplot_i.h"

int short graph_flag;//Flag whether we need graphics(TRUE if !=0)
struct coordNet {
  int xInd;//Index of a var(or other if not var) for x-axis
  int yInd[3];//Index of a vars(or other if not var) for y-axis
  double xUlim;//Limits for the x-axis
  double xLlim;
  double yUlim;//Limits for the y-axis
  double yLlim;
  int short nTics;//Simple implementation of tics: number of major ones.
  int short grid_flag;//Grid yes/no non-0/0.
};
struct coordNet graph;
gnuplot_ctrl * plot_handle;
//int short openpl_flag;//Flag whether we have a page of graphics open.

