/*****************************************************************************************/
/* Copyright 2008,2009,2010,2011 Elias Potapov. */
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
   Tampere University of Technology, Dep. of Signal Processing
   Moscow, Russia / Tampere, Finland
*/
/****************************************************************************************/

double DI;/*Initial distance from reference orbit*/
double DTR;/*Threshold distance value*/
double LPARMAX;/*Max value of parameter when varying for MLE*/
double LPARMIN;/*Min value of parameter when varying for MLE*/
double LDP;/*Step for varying parameter*/
int short LPI;/*Index of the parameter of interest*/
char *mle_file;/*Filename for output parameter-MLE vector*/
char *lyap_spec_file;/*Filename for output parameter-LCE spectrum*/
