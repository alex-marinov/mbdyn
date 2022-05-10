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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "netcdf_output_file.h"
#include "output.h"

virtual void NetCDFOutput::Open(const int format) {
	if (!IsOpen()) {
		m_pBinFile = new netCDF::NcFile(_sPutExt((char*)(OutputHandler::psExt[OutputHandler::NETCDF])), netCDF::NcFile::replace, format); // using the default (nc4) mode was seen to drasticly reduce the writing speed, thus using classic format
		//~ NC_FILL only applies top variables, not files or groups in netcdf-cxx4
		// also: error messages (throw) are part of the netcdf-cxx4 interface by default...

		// Let's define some dimensions which could be useful
		DimTime = CreateDim("time");
		DimV1 = CreateDim("Vec1", 1);
		DimV3 = CreateDim("Vec3", 3);
	}

	return;
};
