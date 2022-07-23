/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

/*
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2022(-2022) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

// This code provides interfaces to INRIA's Siconos library
// https://nonsmooth.gricad-pages.univ-grenoble-alpes.fr/siconos/index.html
// https://github.com/siconos/siconos

#ifdef USE_SICONOS
#include "ls.h"
#include "linsol.h"
#include "vh.h"
#include "siconosmh.h"

class SiconosSolutionManager: public SolutionManager {
protected:
     explicit SiconosSolutionManager(NumericsMatrix_types eMatType, integer iDim, integer iNumNz);

public:
     virtual ~SiconosSolutionManager();

#ifdef DEBUG
     virtual void IsValid() const override;
#endif

     virtual SiconosMatrixHandler* pMatHdl() const override final;

     virtual SiconosVectorHandler* pResHdl() const override final;

     virtual SiconosVectorHandler* pSolHdl() const override final;

     virtual bool bGetConditionNumber(doublereal& dCond) const override;

     virtual void MatrReset() override final;
     virtual void MatrInitialize() override final;
     virtual void Solve() override final;

     SiconosIndexMap* pGetIndexMap() { return &oRowMap; }
protected:
     SiconosIndexMap oRowMap;
     mutable SiconosVectorHandler x;
     mutable SiconosMatrixHandler A;
};

class SiconosDenseSolutionManager: public SiconosSolutionManager {
public:
     explicit SiconosDenseSolutionManager(integer iDim);
     virtual ~SiconosDenseSolutionManager();
};

class SiconosSparseSolutionManager: public SiconosSolutionManager {
public:
     SiconosSparseSolutionManager(integer iDim, integer iNumNz);
     virtual ~SiconosSparseSolutionManager();
};

#endif
