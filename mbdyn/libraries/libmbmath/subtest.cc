/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include "mbconfig.h"

#include "submat.h"
#include "spmapmh.h"

int
main(void)
{
	FullMatrixHandler FMH(4);
	SpMapMatrixHandler SMH(4);

	FullSubMatrixHandler FSMH(2, 2);
	FSMH.ResizeReset(2, 2);

	FSMH.PutRowIndex(1, 3);
	FSMH.PutRowIndex(2, 4);

	FSMH.PutColIndex(1, 1);
	FSMH.PutColIndex(2, 2);

	FSMH(1, 1) = 11.;
	FSMH(1, 2) = 12.;
	FSMH(2, 1) = 21.;
	FSMH(2, 2) = 22.;

	std::cout << "FSMH: " << std::endl << FSMH << std::endl;

	FMH.Reset();
	FSMH.AddTo(FMH);

	std::cout << "FMH = FSMH: " << std::endl << FMH << std::endl;

	FMH.Reset();
	FSMH.SubFrom(FMH);

	std::cout << "FMH = -FSMH: " << std::endl << FMH << std::endl;

	FMH.Reset();
	FSMH.AddToT(FMH);

	std::cout << "FMH = FSMH^T: " << std::endl << FMH << std::endl;

	FMH.Reset();
	FSMH.SubFromT(FMH);

	std::cout << "FMH = -FSMH^T: " << std::endl << FMH << std::endl;

	SMH.Reset();
	FSMH.AddTo(SMH);

	std::cout << "SMH = FSMH: " << std::endl << SMH << std::endl;

	SMH.Reset();
	FSMH.SubFrom(SMH);

	std::cout << "SMH = -FSMH: " << std::endl << SMH << std::endl;

	SMH.Reset();
	FSMH.AddToT(SMH);

	std::cout << "SMH = FSMH^T: " << std::endl << SMH << std::endl;

	SMH.Reset();
	FSMH.SubFromT(SMH);

	std::cout << "SMH = -FSMH^T: " << std::endl << SMH << std::endl;

	return 0;
}
