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
 
#ifndef SOLVERDIAGNOSTICS_H
#define SOLVERDIAGNOSTICS_H

#include "drive.h"

class SolverDiagnostics {
protected:
 	unsigned OutputFlags;
	DriveCaller *pOutputMeter;

	enum {
		OUTPUT_NONE		= 0x0000,

		OUTPUT_ITERS		= 0x0001,
		OUTPUT_RES		= 0x0002,
		OUTPUT_SOL		= 0x0004,
		OUTPUT_JAC		= 0x0008,

		OUTPUT_STEP		= (OUTPUT_ITERS | OUTPUT_RES | OUTPUT_SOL | OUTPUT_JAC),

		OUTPUT_BAILOUT		= 0x0010,
		OUTPUT_MSG		= 0x0020,
		OUTPUT_COUNTER		= 0x0040,

		OUTPUT_DEFAULT		= OUTPUT_MSG,

		OUTPUT_MASK		= 0x00FF
	};
public:

	SolverDiagnostics(unsigned OF = OUTPUT_DEFAULT, DriveCaller *pOM = 0);
	virtual ~SolverDiagnostics(void);
	
	void SetNoOutput(void);

	void SetOutputMeter(DriveCaller *pOM);

	void SetOutputFlags(unsigned OF);
	void AddOutputFlags(unsigned OF);
	void DelOutputFlags(unsigned OF);

	inline bool outputMeter(void) const {
		return (!pOutputMeter || pOutputMeter->dGet());
	};

	inline bool outputIters(void) const {
		return outputMeter() && (OutputFlags & OUTPUT_ITERS);
	};
 
	inline bool outputRes(void) const {
		return outputMeter() && (OutputFlags & OUTPUT_RES);
	};
 
	inline bool outputSol(void) const {
		return outputMeter() && (OutputFlags & OUTPUT_SOL);
	};
 
	inline bool outputJac(void) const {
		return outputMeter() && (OutputFlags & OUTPUT_JAC);
	};

	inline bool outputStep(void) const {
		return outputMeter() && (OutputFlags & OUTPUT_STEP);
	};

	inline bool outputBailout(void) const {
		// NOTE: this is NOT conditioned by output meter
		return (OutputFlags & OUTPUT_BAILOUT);
	};

	inline bool outputCounter(void) const {
		return (OutputFlags & OUTPUT_COUNTER);
	};

	// all messages not protected behind any other condition
	// must be protected by a "if (outputMsg())" condition
	inline bool outputMsg(void) const {
		return (OutputFlags & OUTPUT_MSG);
	};
};

#endif /* SOLVERDIAGNOSTICS_H */

