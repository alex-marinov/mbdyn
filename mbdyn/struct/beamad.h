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

#ifndef BEAMAD_H
#define BEAMAD_H

#include <array>

#include "beam.h"
#include "strnodead.h"
#include "sp_matvecass.h"

class BeamAd: virtual public Beam {
public:
     BeamAd(unsigned int uL,
            const StructNodeAd* pN1, const StructNodeAd* pN2, const StructNodeAd* pN3,
            const Vec3& F1, const Vec3& F2, const Vec3& F3,
            const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3,
            const Mat3x3& r_I, const Mat3x3& rII,
            const ConstitutiveLaw6D* pD_I, const ConstitutiveLaw6D* pDII,
            OrientationDescription ood,
            flag fOut);

     virtual ~BeamAd();

     virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const override;

     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    enum sp_grad::SpFunctionCall func);

protected:
     virtual void
     AddInternalForces(sp_grad::SpColVector<doublereal, 6>& AzLoc, unsigned int iSez);
     
     virtual void
     AddInternalForces(sp_grad::SpColVector<sp_grad::SpGradient, 6>& AzLoc, unsigned int iSez);

     virtual void
     AddInternalForces(sp_grad::SpColVector<sp_grad::GpGradProd, 6>& AzLoc, unsigned int iSez);
        
     template <typename T>
     inline void
     AssReactionForce(sp_grad::SpGradientAssVec<T>& WorkVec,
                      const std::array<sp_grad::SpColVectorA<T, 3>, NUMSEZ>& p,
                      const std::array<sp_grad::SpColVectorA<T, 6>, NUMSEZ>& Az,
                      const std::array<sp_grad::SpColVectorA<T, 3>, NUMNODES>& X) const;
     
     template <typename T>
     static sp_grad::SpColVector<T, 3>
     InterpState(const sp_grad::SpColVector<T, 3>& v1,
                 const sp_grad::SpColVector<T, 3>& v2,
                 const sp_grad::SpColVector<T, 3>& v3,
                 Section Sec);

     template <typename T>
     sp_grad::SpColVector<T, 3>
     InterpDeriv(const sp_grad::SpColVector<T, 3>& v1,
                 const sp_grad::SpColVector<T, 3>& v2,
                 const sp_grad::SpColVector<T, 3>& v3,
                 Section Sec);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLoc);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& AzLoc) {}

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::GpGradProd, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& AzLoc) {}
     
protected:
     const std::array<const StructNodeAd*, NUMNODES> pNode;
};

class ViscoElasticBeamAd: public ViscoElasticBeam, public BeamAd {
public:
     ViscoElasticBeamAd(unsigned int uL,
                        const StructNodeAd* pN1,
                        const StructNodeAd* pN2,
                        const StructNodeAd* pN3,
                        const Vec3& F1,
                        const Vec3& F2,
                        const Vec3& F3,
                        const Mat3x3& R1,
                        const Mat3x3& R2,
                        const Mat3x3& R3,
                        const Mat3x3& r_I,
                        const Mat3x3& rII,
                        const ConstitutiveLaw6D* pD_I,
                        const ConstitutiveLaw6D* pDII,
                        OrientationDescription ood,
                        flag fOut);
     
     virtual ~ViscoElasticBeamAd();
     
     virtual SubVectorHandler&
     AssRes(SubVectorHandler& WorkVec,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;
     
     virtual VariableSubMatrixHandler&
     AssJac(VariableSubMatrixHandler& WorkMat,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr) override;

     virtual void
     AssJac(VectorHandler& JacY,
            const VectorHandler& Y,
            doublereal dCoef,
            const VectorHandler& XCurr,
            const VectorHandler& XPrimeCurr,
            VariableSubMatrixHandler& WorkMat) override;

     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    enum sp_grad::SpFunctionCall func);
     
protected:
     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gPrime,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& Omega,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LPrime,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefPrimeLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLoc);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& gPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& Omega,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& LPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefPrimeLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& AzLoc) {}

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::GpGradProd, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& gPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& Omega,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 3>, NUMSEZ>& LPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& DefPrimeLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::GpGradProd, 6>, NUMSEZ>& AzLoc) {}
};

#endif /* BEAMAD_H */
