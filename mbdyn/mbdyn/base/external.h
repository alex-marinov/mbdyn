/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* external Communications */


#ifdef USE_MPI 

#ifndef EXTERNAL_HH
#define EXTERNAL_HH

#include "ac/mpi.h"
#include <list>
#define INTERF_COMM_LABEL 1000
extern std::list<MPI::Intercomm>  InterfaceComms; 

namespace External {

	enum ExtMessage{
		EMPTY,
		INITIAL,
		REGULAR,
		CLOSE,
		ERROR
	};
	
	void SendNull(void);

	void SendInitial(void);

	void SendFreeze(void);
	
	void SendRegular(void);

	void SendClose(void);

	void SendError(void);

}


#endif /* EXTERNAL_HH */

#endif /* USE_MPI */
