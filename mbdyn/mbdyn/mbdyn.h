/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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

#ifndef MBDYN_H
#define MBDYN_H

/* Global typedefs (unused yet) */
typedef double mbReal;
typedef int mbInt;
typedef long int mbLong;
typedef long int mbFlag;

/* Global variables */
extern int fSilent;

/* Global macros */
#define silent_cout(arg) \
    do { \
        if (!::fSilent) { \
            cout << arg; \
        } \
    } while (0)

/* Debug levels (from 0x0001 to 0x0080 are reserved) */
static const long int MYDEBUG_INPUT               = 0x00000100;
static const long int MYDEBUG_ASSEMBLY            = 0x00000200;
static const long int MYDEBUG_DERIVATIVES         = 0x00000400;
static const long int MYDEBUG_FSTEPS              = 0x00000800;
static const long int MYDEBUG_MEM                 = 0x00001000;
static const long int MYDEBUG_MPI                 = 0x00002000;
static const long int MYDEBUG_PRED                = 0x00004000;
static const long int MYDEBUG_RESIDUAL            = 0x00008000;
static const long int MYDEBUG_SOL                 = 0x00010000;
static const long int MYDEBUG_INIT                = 0x00020000;
static const long int MYDEBUG_OUTPUT              = 0x00040000;

#endif /* MBDYN_H */

