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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

/* libreria per il calcolo delle partizioni */
#include <metis.h>
#undef ASSERT /* kill annoying redefiniton message */
#include "mynewmem.h"

int
mbdyn_METIS_PartGraph(int iTotVertices,
		int *pXadj,
		int *pAdjncy,
		int *pVertexWgts,
		int *pCommWgts,
		int *pEdgeWgts,
		int DataCommSize,
		int *pParAmgProcs)
{
	/* required by METIS_PartGraphVKway API, but otherwise ignored */
	/* 0: C-style numbering [0..n-1]; 1: F77-style numbering (1..n) */
	int	numflag = 0;
	/* if options[0] == 0, the rest is ignored */
	int	options[5] = { 0 };
	/* total communication volume */
	int	volume = 0;
	/* weiht flags */
	int	wgtflag;

	idxtype	*xadj = 0,
		*adjncy = 0,
		*vwgt = 0,
		*adjwgt = 0,
		*part = 0;

	if (sizeof(idxtype) == sizeof(int)) {
		xadj = pXadj;
		adjncy = pAdjncy;
		vwgt = pVertexWgts;
		adjwgt = pCommWgts;
		part = pParAmgProcs;

	} else {
		SAFENEWARR(xadj, idxtype, iTotVertices);
		SAFENEWARR(adjncy, idxtype, pXadj[iTotVertices]);
		if (pVertexWgts) {
			SAFENEWARR(vwgt, idxtype, iTotVertices);
		}
		if (pCommWgts) {
			SAFENEWARR(adjwgt, idxtype, pXadj[iTotVertices]);
		}
		SAFENEWARR(part, idxtype, iTotVertices);

		for (int i = 0; i < iTotVertices; i++) {
			xadj[i] = pXadj[i];
			part[i] = pParAmgProcs[i];
		}

		if (pVertexWgts) {
			for (int i = 0; i < iTotVertices; i++) {
				vwgt[i] = pVertexWgts[i];
			}
		}

		for (int i = 0; i < pXadj[iTotVertices]; i++) {
			adjncy[i] = pAdjncy[i];
		}

		if (pCommWgts) {
			for (int i = 0; i < pXadj[iTotVertices]; i++) {
				adjwgt[i] = pCommWgts[i];
			}
		}
	}

	if (pVertexWgts && pCommWgts) {
		wgtflag = 3;
	} else if (pVertexWgts) {
		wgtflag = 2;
	} else if (pCommWgts) {
		wgtflag = 1;
	} else {
		wgtflag = 0;
	}

	METIS_PartGraphVKway(&iTotVertices,
			xadj,
			adjncy,
			vwgt,
			adjwgt,
			&wgtflag,
			&numflag,
			&DataCommSize,
			options,
			&volume,
			part);


	if (sizeof(idxtype) != sizeof(int)) {
		for (int i = 0; i < iTotVertices; i++) {
			pParAmgProcs[i] = part[i];
		}

		SAFEDELETEARR(xadj);
		SAFEDELETEARR(adjncy);

		if (pVertexWgts) {
			SAFEDELETEARR(vwgt);
		}

		if (pCommWgts) {
			SAFEDELETEARR(adjwgt);
		}

		SAFEDELETEARR(part);
	}

	/* NOTE: the manual suggests to use
	 * METIS_PartGraphRecursive if DataCommSize < 8 */

	return 0;
}

