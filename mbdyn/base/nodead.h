/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2022
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2022) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#ifndef NODEAD_H
#define NODEAD_H

#include "myassert.h"
#include "node.h"
#include "sp_gradient.h"

class ScalarNodeAd: virtual public ScalarNode {
public:
     ScalarNodeAd(unsigned int uL, const DofOwner* pDO, flag fOut);
     virtual ~ScalarNodeAd();

     inline void GetX(doublereal& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetX(sp_grad::SpGradient& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetX(sp_grad::GpGradProd& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     virtual void UpdateJac(const VectorHandler& Y, doublereal dCoef) override;

protected:
     doublereal XY;
};

class ScalarDifferentialNodeAd: virtual public ScalarDifferentialNode, public ScalarNodeAd {
public:
     ScalarDifferentialNodeAd(unsigned int uL, const DofOwner* pDO,
                              const doublereal& dx, const doublereal& dxp, flag fOut);
     virtual ~ScalarDifferentialNodeAd();

     inline void GetXPrime(doublereal& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetXPrime(sp_grad::SpGradient& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const;
     inline void GetXPrime(sp_grad::GpGradProd& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const;
};

class ScalarAlgebraicNodeAd: virtual public ScalarAlgebraicNode, public ScalarNodeAd {
public:
     ScalarAlgebraicNodeAd(unsigned int uL, const DofOwner* pDO,
                           doublereal dx, flag fOut);

     virtual ~ScalarAlgebraicNodeAd();
};

inline void ScalarNodeAd::GetX(doublereal& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     dX = dGetX();
}

inline void ScalarNodeAd::GetX(sp_grad::SpGradient& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     const doublereal dDeriv = GetDofType(0) == DofOrder::DIFFERENTIAL ? -dCoef : -1;
     dX.Reset(dGetX(), iGetFirstColIndex() + 1, dDeriv);
}

inline void ScalarNodeAd::GetX(sp_grad::GpGradProd& dX, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     const doublereal dDeriv = GetDofType(0) == DofOrder::DIFFERENTIAL ? -dCoef : -1;

     dX.Reset(dGetX(), XY * dDeriv);
}

inline void ScalarDifferentialNodeAd::GetXPrime(doublereal& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     dXPrime = dGetXPrime();
}

inline void ScalarDifferentialNodeAd::GetXPrime(sp_grad::SpGradient& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     dXPrime.Reset(dGetXPrime(), iGetFirstColIndex() + 1, -1.);
}

inline void ScalarDifferentialNodeAd::GetXPrime(sp_grad::GpGradProd& dXPrime, doublereal dCoef, sp_grad::SpFunctionCall func) const
{
     dXPrime.Reset(dGetXPrime(), -XY);
}

#endif /* NODEAD_H */
