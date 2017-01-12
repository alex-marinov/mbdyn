/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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


#include "mbconfig.h" 		/* This goes first in every *.c,*.cc file */


#include "external.h"

#ifdef USE_EXTERNAL

namespace External {

/* setta il sistema per una connessione con /dev/null */
void SendNull(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'N';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

/* setta il sistema per una comunicare che gli ultimi dati inviati 
   sono le condizioni iniziali post assemblaggio e passi fittizi */
void SendInitial(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'I';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

/* setta il sistema per una ricevere le stesse forze del passo precedente */
void SendFreeze(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'F';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

/* setta il sistema per una connessione regolare con il codice interfaciato */
void SendRegular(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'R';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

/* setta il sistema per la chiusura del calcolo */
void SendClose(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'Q';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

/* setta il sistema per la chiusura del calcolo in seguito ad un errore */
void SendError(void) {
	std::list<MPI::Intercomm>::iterator iComm; 
	std::list<MPI::Intercomm>::const_iterator iComm_end = InterfaceComms.end(); 
	char buff = 'X';
	for(iComm = InterfaceComms.begin(); iComm != iComm_end; ++iComm) {
		iComm->Send(&buff, 1, MPI::CHAR, 0, 0);
	}
}

}
#endif /* USE_EXTERNAL */
