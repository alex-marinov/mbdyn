/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */			   

/******************************************************************************

	Set di typedef da includere sempre

******************************************************************************/


#ifndef _MYSTDDEF
#define _MYSTDDEF

#include <iostream.h>
#include <iomanip.h>

// Definizioni generali
#ifdef __alpha
typedef int flag; // Usato come variabile booleana di ritorno
#else
typedef long int flag; // Usato come variabile booleana di ritorno
#endif 
typedef char     byte; // Usato come mattone fondamentale di tipi



// Definizioni matematiche
// typedef unsigned int  Index;
typedef int           Int;
typedef long int      Lint;
typedef double        Real;


#define OUTSET(prec, width) \
          setprecision(prec) << setw(width)

#endif
