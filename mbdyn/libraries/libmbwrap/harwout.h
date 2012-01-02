/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 * the Free Software Foundation (version 2 of the License).
 * 
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

#ifndef HARWOUT_H
#define HARWOUT_H

#ifdef USE_HARWELL

enum {
   MA28AD_99999,
   MA28AD_99998,
   MA28AD_99997,
   MA28AD_99996,
   MA28AD_99995,
   MA28AD_99994,
   MA28AD_99993,
   MA28AD_99992,
   MA28AD_99991,
   
   MA28BD_99999,
   MA28BD_99998,
   MA28BD_99997,
   MA28BD_99996,
   MA28BD_99995,
   MA28BD_99994,
   
   MA28CD_99999,
   MA28CD_99998,
   MA28CD_99997,

   MA28ID_99999,
   MA28ID_99998,

   MA30AD_99999,
   MA30AD_99998,
   MA30AD_99997,
   MA30AD_99996,
   MA30AD_99995,
   MA30AD_99994,
   MA30AD_99993,
   MA30AD_99992,
   
   MA30BD_99999,
   
   MC23A_180,
   MC23A_200
};

extern int harwell_error(int, const char* const v[]);

#endif /* USE_HARWELL */

#endif /* HARWOUT_H */
