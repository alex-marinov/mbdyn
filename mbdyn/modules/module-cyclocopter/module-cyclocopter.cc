/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2009
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
 * Author: Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>
 *
 * Reworked as runtime loadable module by
 * Pierangelo Masarati <pierangelo.masarati@polimi.it>
 */

/* Elementi di rotore */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <limits>
#include <cmath>

#include "dataman.h"
#include "userelem.h"
#include "indvel.h"

/* CyclocopterInflow - begin */

/* Classe base per i modelli di influsso del ciclocottero:
 * tutti gli altri modelli sono ottenuti ereditando questa
 * classe
*/

class CyclocopterInflow
: virtual public Elem, public UserDefinedElem, public InducedVelocity {
protected:
	const StructNode* pRotor;

	bool bFlagAverage;
	// == true: usa la media delle forze sul giro
	//	per calcolare la velocità indotta per il giro successivo
	// == false: usa il valore istantaneo (eventualmente filtrato)


	doublereal dRadius;		// Rotor radius
	doublereal dSpan;		// Blade length
	doublereal dArea;		// Cylinder longitudinal area

	doublereal dKappa;		// Hover correction coefficient

	DriveOwner Weight;
	doublereal dWeight;

	Mat3x3 RRot;
	/* dati che servono per il calcolo dell'output */
	Mat3x3 RRotorTranspose;
	doublereal dUindMean;

	// Coefficienti del filtro di butterworth del second'ordine
	doublereal a1, a2, b0, b1, b2;

	// Coefficienti del filtro di butterworth del second'ordine
	void SetFilterCoefficients(doublereal dOmegaFilter, doublereal dDeltaT);

public:
	CyclocopterInflow(unsigned int uL, const DofOwner* pDO);
	virtual ~CyclocopterInflow(void);

	virtual Elem::Type GetElemType(void) const;
	virtual InducedVelocity::Type GetInducedVelocityType(void) const;

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void SetValue(DataManager *pDM,
		VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);

	// Relativo ai ...WithDofs
	virtual void SetInitialValue(VectorHandler& X);

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void 
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat, 
		      const VectorHandler& XCurr);
   	SubVectorHandler& 
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

CyclocopterInflow::CyclocopterInflow(unsigned int uL, const DofOwner* pDO)
: Elem(uL, flag(0)),
UserDefinedElem(uL, pDO),
InducedVelocity(uL, 0, 0, flag(0)),
pRotor(0),
bFlagAverage(false),
dRadius(0.),
dSpan(0.),
dArea(0.),
Weight(0),
dWeight(0.),
RRot(::Eye3),
RRotorTranspose(::Eye3),
dUindMean(0.),
a1(0.), a2(0.), b0(1.), b1(0.), b2(0.)
{
	NO_OP;
}

CyclocopterInflow::~CyclocopterInflow(void)
{
	NO_OP;
}

Elem::Type
CyclocopterInflow::GetElemType(void) const
{
	return Elem::LOADABLE;
}

InducedVelocity::Type
CyclocopterInflow::GetInducedVelocityType(void) const
{
	return InducedVelocity::USER_DEFINED;
}

void
CyclocopterInflow::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterInflow::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << "0."                	 /* 9 */
                        << " " << "0."                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << "0."                	 /* 14 */
                        << " " << "0."   		 /* 15 */
                        << " " << "0."           	 /* 16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

std::ostream&
CyclocopterInflow::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

void
CyclocopterInflow::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

void
CyclocopterInflow::SetInitialValue(VectorHandler& X)
{
	NO_OP;
}

void
CyclocopterInflow::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(2);
	connectedNodes[0] = pCraft;
	connectedNodes[1] = pRotor;
}

unsigned int
CyclocopterInflow::iGetInitialNumDof(void) const
{
	return 0;
}

void 
CyclocopterInflow::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
CyclocopterInflow::InitialAssJac(
	VariableSubMatrixHandler& WorkMat, 
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler& 
CyclocopterInflow::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

void
CyclocopterInflow::SetFilterCoefficients(doublereal dOmegaFilter,
	doublereal dDeltaT)
{
	/* Butterworth discrete low-pass filter coefficients */
	if (dDeltaT > 0. && dOmegaFilter > 0.) {
		doublereal dTmp = 4. + 2.*sqrt(2.)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter;
		a1 = (-8. + 2.*dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;
		a2 = (4. - 2.*sqrt(2.)*dOmegaFilter*dDeltaT + dDeltaT*dDeltaT*dOmegaFilter*dOmegaFilter)/dTmp;

		dTmp = dOmegaFilter*dOmegaFilter*dDeltaT*dDeltaT/dTmp;
		b0 = dTmp;
		b1 = 2.*dTmp;
		b2 = dTmp;
	}
}

static bool
ReadRotorData(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel,
	const StructNode *& pCraft,
	Mat3x3& rrot,
	const StructNode *& pRotor)
{
     	/* aircraft node */
	pCraft = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	/* rotor orientation with respect to aircraft */
     	rrot = ::Eye3;
     	if (HP.IsKeyWord("orientation")) {
     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);

     	} else if (HP.IsKeyWord("hinge")) {
		silent_cerr("InducedVelocity(" << uLabel << "): deprecated keyword \"hinge\"; use \"orientation\" instead at line " << HP.GetLineData() << std::endl);

     		ReferenceFrame RF(pCraft);
     		rrot = HP.GetRotRel(RF);
     	}

     	/* rotor node */
     	pRotor = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	return true;
}

static bool
ReadUniform(DataManager* pDM,
	MBDynParser& HP,
	unsigned int uLabel,
	bool& bFlagAve,
	doublereal& dR,
	doublereal& dL,
	DriveCaller *& pdW,
	doublereal& dOmegaFilter,
	doublereal& dKappa,
	doublereal& dDeltaT)
{
	bFlagAve = HP.GetYesNoOrBool();

	dR = HP.GetReal();
	if (dR <= 0.) {
		silent_cerr("ReadUniform(" << uLabel << "): "
			"illegal null or negative radius"
			"for rotor" << uLabel << " at line " << HP.GetLineData()
			<< std::endl);
		return false;
	}

	dL = HP.GetReal();
	if (dL <= 0.) {
		silent_cerr("ReadUniform(" << uLabel << "): "
			"illegal null or negative blade"
			"length for rotor" << uLabel
			<< " at line " << HP.GetLineData()
			<< std::endl);
		return false;
	}

	pdW = 0;
	if (HP.IsKeyWord("delay")) {
		pdW = HP.GetDriveCaller();

	} else {
		SAFENEW(pdW, NullDriveCaller);
	}

	dOmegaFilter = 0.;
	if (HP.IsKeyWord("omegacut")) {
		dOmegaFilter = HP.GetReal();
		if (dOmegaFilter <= 0) {
			silent_cerr("Illegal null or negative filter"
				"cut frequency for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}
	} else {
		dOmegaFilter = 0.;
	}

	dKappa = 1.;
	if (HP.IsKeyWord("kappa")) {
		dKappa = HP.GetReal();
		if (dKappa <= 0) {
			silent_cerr("Illegal null or negative hover"
				"correction coefficient" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}
	}

	dDeltaT = 0.;
	if (HP.IsKeyWord("timestep")) {
		dDeltaT = HP.GetReal();
		if (dDeltaT <= 0) {
			silent_cerr("Illegal null or negative time"
				"step for rotor" << uLabel
				<< " at line " << HP.GetLineData()
				<< std::endl);
			return false;
		}

	} else {
		dDeltaT = 0.;
	}

	return true;
}

/* CyclocopterInflow - end */


/* CyclocopterNoInflow - begin */

class CyclocopterNoInflow
: virtual public Elem, public CyclocopterInflow {
public:
	CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterNoInflow(void);

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	}
	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	}
	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	}

};

CyclocopterNoInflow::CyclocopterNoInflow(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal no inflow ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));
}

CyclocopterNoInflow::~CyclocopterNoInflow(void)
{
	NO_OP;
}

SubVectorHandler&
CyclocopterNoInflow::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	if (fToBeOutput()) {
		RRotorTranspose = pCraft->GetRCurr()*RRot;
		RRotorTranspose = RRotorTranspose.Transpose();
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterNoInflow::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);
	}
}

Vec3
CyclocopterNoInflow::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{

	return Zero3;
}

/* CyclocopterNoInflow - end */

/* CyclocopterUniform1D - begin */

/*

From Moble: the induced velocity is opposite to the
the force generated by the rotor in the direction 3
of the rotor reference (the direction 1 of this
reference must be align with the rotation axis of
the rotor). It should make sense just for a multiblade
rotor (not for the one-blade rotor) when the generated
force is mainly in one direction.

*/

class CyclocopterUniform1D
: virtual public Elem, public CyclocopterInflow {
protected:
	Vec3 RRot3;
	Mat3x3 RRotor;

	mutable doublereal dUindMeanPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	doublereal dTz, dTzMean;
	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	doublereal Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterUniform1D(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterUniform1D::CyclocopterUniform1D(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRot3(::Zero3),
RRotor(::Eye3),
dUindMeanPrev(0.),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
dTz(0.), dTzMean(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(0.), Uk_1(0.), Uk_2(0.), Yk(0.), Yk_1(0.), Yk_2(0.)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal uniform 1D ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterUniform1D::~CyclocopterUniform1D(void)
{
	NO_OP;
}

void
CyclocopterUniform1D::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << "0."                	 /* 11 */
                        << " " << "0."                	 /* 12 */
                        << " " << "0."                	 /* 13 */
                        << " " << FMeanOut                	 /* 14 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

void
CyclocopterUniform1D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
	/* calcolo la forza media sul giro generata dal rotore */
	dTzMean += dTz;
	FMean += F;
	iStepCounter++;
	// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean /= iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			dTzMean = dTzMean/iStepCounter;
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*copysign(std::sqrt(std::abs(dTzMean)/(2*dRho*dArea)), dTzMean);
			dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
			dTzMean = 0.;
		}
		FMean = ::Zero3;
		iStepCounter = 0;
	}
	dAzimuthPrev = dAzimuth;

	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dUindMeanPrev = dUindMean;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;
	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

SubVectorHandler&
CyclocopterUniform1D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity (Moble version)*/
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRot3 = RRotor.GetVec(3);
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	dTz = RRot3*Res.Force();
	if (!bFlagAverage) {
		/* filtro le forze */
		Uk = dTz;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		dTz = Yk;
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*copysign(std::sqrt(std::abs(dTz)/(2*dRho*dArea)), dTz);

		dUindMean = (1 - dWeight)*dUindMean + dWeight*dUindMeanPrev;
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterUniform1D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{

	/* colcolo la posizione azimutale della prima pala */
	if (bFlagIsFirstBlade == true) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}


	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform1D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return RRot3*dUindMean;
}

/* CyclocopterUniform1D - end */

/* CyclocopterUniform2D - begin */

/*

Uniform inflow: the induced velocity is opposite to
the force generated by the rotor in the plane
perdendicular to the rotor rotation axis. The
rotor reference must have direction 1 aligned
with the rotor rotation axis!

*/

class CyclocopterUniform2D
: virtual public Elem, public CyclocopterInflow {
protected:
	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

public:
	CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterUniform2D(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterUniform2D::CyclocopterUniform2D(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRotor(::Eye3),
dUind(::Zero3), dUindPrev(::Zero3),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(::Zero3), Uk_1(::Zero3), Uk_2(::Zero3), Yk(::Zero3), Yk_1(::Zero3), Yk_2(::Zero3)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal uniform 2D ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterUniform2D::~CyclocopterUniform2D(void)
{
	NO_OP;
}

void
CyclocopterUniform2D::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean                	 /* 8 */
                        << " " << dAzimuth                	 /* 9 */
                        << " " << iStepCounter                	 /* 10 */
                        << " " << dUind                	 /* 11-13 */
                        << " " << FMeanOut                	 /* 14-16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}

}

void
CyclocopterUniform2D::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
#if 0
	if (bFlagAverage) {
		/* calcolo la forza media sul giro generata dal rotore */
		FMean = FMean + F;
		iStepCounter++;
		// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
		if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
			FMean = FMean/iStepCounter;
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema
	 		* rotore */
			dUind = 0.;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

			FMean = 0.;
			iStepCounter = 0;
		}
	}
#endif

	/* calcolo la forza media sul giro generata dal rotore */
	FMean = FMean + F;
	iStepCounter++;
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema
	 		* rotore */
			dUind = ::Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
		}

		FMean = ::Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;

	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

SubVectorHandler&
CyclocopterUniform2D::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	if (!bFlagAverage) {
		/* filtro le forze */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;
		/* Forza nel piano normale all'asse di rotazione */
		doublereal dT = sqrt(F(2)*F(2) + F(3)*F(3));
		/* Velocità indotta: calcolata in base alla dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
		/* Componenti della velocità indotta nel sistema
	 	* rotore */
		dUind = ::Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterUniform2D::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* colcolo la posizione azimutale della prima pala */
	// if (bFlagIsFirstBlade && bFlagAverage) {
	if (bFlagIsFirstBlade) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}

	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F, M, X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterUniform2D::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	//printf("%f %f %f\n",dUind(1),dUind(2),dUind(3));
	return RRotor*dUind;
}

/* CyclocopterUniform2D - end */


/* CyclocopterPolimi - begin */

/*

The induced velocity is opposite to
the force generated by the rotor in the plane
perdendicular to the rotor rotation axis. The
rotor reference must have direction 1 aligned
with the rotor rotation axis!
The inflow velocity distribution is constant along
the cylinder span and is function of r:
Vi(r) = Kc*cos((pi/2)*(r/R)) + Ks*sin(pi*(r/R))
where the cofficients Kc and Ks are based on the
momentum theory

*/

class CyclocopterPolimi
: virtual public Elem, public CyclocopterInflow {
protected:
	Mat3x3 RRotor;
	Vec3 dUind;
	mutable Vec3 dUindPrev;

	doublereal dXi;

	bool bFlagIsFirstBlade;

	doublereal dAzimuth, dAzimuthPrev;

	Vec3 F, FMean, FMeanOut;

	unsigned int iStepCounter;

	/* dati per il filtraggio delle forze */
	Vec3 Uk, Uk_1, Uk_2, Yk, Yk_1, Yk_2;

	unsigned int iCounter;
	unsigned int iRotationCounter;

	doublereal dpPrev, dp;

	bool flag_print;

public:
	CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterPolimi(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// Restituisce ad un elemento la velocità indotta
	// in base alla posizione azimuthale (usata per
	// iterare tra calcolo della valocità indotta e
	// calcolo delle forze aerodinamiche - solo
	// per il modello di influsso KARI per il ciclocottero)
#if 0
	virtual void GetInducedVelocityIter(const Vec3& X, const Vec3& T, doublereal *UindM, doublereal dTn0, doublereal dTn_dUindM) {
		NO_OP;
	};
#endif

	// Restituisce la velocità indotta dalla metà superiore del rotore
	// ( solo per KARI ciclocottero)
	virtual doublereal GetW(const Vec3& X) const {
		return 0.;
	};

	virtual doublereal GetPsi(const Vec3& X) const {
		return 0.;
	};

	virtual Mat3x3 GetRRotor(const Vec3& X) const {
		return ::Zero3x3;
	};
};

CyclocopterPolimi::CyclocopterPolimi(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, flag(0)),
CyclocopterInflow(uL, pDO),
RRotor(::Eye3),
dUind(::Zero3), dUindPrev(::Zero3),
dXi(0.),
bFlagIsFirstBlade(true),
dAzimuth(0.), dAzimuthPrev(0.),
F(::Zero3), FMean(::Zero3), FMeanOut(::Zero3),
iStepCounter(0),
Uk(::Zero3), Uk_1(::Zero3), Uk_2(::Zero3), Yk(::Zero3), Yk_1(::Zero3), Yk_2(::Zero3),
iCounter(0), iRotationCounter(0),
dpPrev(0.), dp(0.),
flag_print(true)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
"									\n"
"Module: 	Cyclocopter						\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
"based on work by							\n"
"		Mattia Mattaboni <mattia.mattaboni@mail.polimi.it>	\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
"		Politecnico di Milano					\n"
"		http://www.aero.polimi.it/				\n"
" Description:	This module implements induced velocity models		\n"
"		for cycloidal rotors.					\n"
"									\n"
"	All rights reserved.						\n"
"\n"
" Usage:\n"
"	user element: <label> , cycloidal Polimi ,\n"
"		<aircraft_node_label> ,\n"
"		[ orientation , (OrientationMatrix) <orientation> , ]\n"
"		<rotor_node_label>\n"
"		(bool) <average> ,\n"
"		<rotor_radius> ,\n"
"		<blade_span>\n"
"		[ , delay , (DriveCaller) <delay> ]\n"
"		[ , omegacut , <cut_frequency> ]\n"
"		[ , kappa , <hover_correction_coefficient> ]\n"
"		[ , timestep , <time_step> ]\n"
"		[ , <output_data> ]\n"
"	;\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	DriveCaller *pdW = 0;
	doublereal dOmegaFilter;
	doublereal dDeltaT;
	if (!ReadUniform(pDM, HP, uLabel, bFlagAverage, dRadius, dSpan, pdW, dOmegaFilter, dKappa, dDeltaT)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));

	dArea = 2*dRadius*dSpan;
	Weight.Set(pdW);

	SetFilterCoefficients(dOmegaFilter, dDeltaT);
}

CyclocopterPolimi::~CyclocopterPolimi(void)
{
	NO_OP;
}

void
CyclocopterPolimi::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	bFlagIsFirstBlade = true;
	/* calcolo la forza media sul giro generata dal rotore */
	FMean = FMean + F;;
	iStepCounter++;
	// if ((dAzimuth > 0. && dAzimuthPrev < 0.) || (dAzimuth < 0. && dAzimuthPrev > 0.)) {
	if ((dAzimuth > 0. && dAzimuthPrev < 0.)) {
		FMean = FMean/iStepCounter;
		FMeanOut = FMean;
		if (bFlagAverage) {
			/* Forza nel piano normale all'asse di rotazione */
			doublereal dT = sqrt(FMean(2)*FMean(2) + FMean(3)*FMean(3));
			/* Velocità indotta: calcolata in base alla dT */
			doublereal dRho = dGetAirDensity(GetXCurr());
			dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
			/* Componenti della velocità indotta nel sistema
	 		* rotore */
			dUind = Zero3;
			if (dT > std::numeric_limits<doublereal>::epsilon()) {
				dUind(2) = dUindMean*FMean(2)/dT;
				dUind(3) = dUindMean*FMean(3)/dT;
			}
			dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
			dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
			dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);
			/* angolo di cui è ruotata la trazione */
			dXi = atan2(FMean(3), FMean(2)) - M_PI/2.;
		}

		FMean = Zero3;
		iStepCounter = 0;
	}

	dAzimuthPrev = dAzimuth;

	dUindPrev = dUind;

	/* aggiorno ingressi e uscite del filtro */
	Yk_2 = Yk_1;
	Yk_1 = Yk;
	Uk_2 = Uk_1;
	Uk_1 = Uk;

	dWeight = Weight.dGet();
	if (dWeight < 0.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay < 0.0; using 0.0" << std::endl);
		dWeight = 0.;

	} else if (dWeight > 1.) {
		silent_cout("Rotor(" << GetLabel() << "): "
			"delay > 1.0; using 1.0" << std::endl);
		dWeight = 1.;
	}

	InducedVelocity::AfterConvergence(X, XP);
}

void
CyclocopterPolimi::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
                OH.Loadable()
                        << std::setw(8) << GetLabel()   /* 1 */
                        << " " << RRotorTranspose*Res.Force()     /* 2-4 */
                        << " " << RRotorTranspose*Res.Moment()    /* 5-7 */
                        << " " << dUindMean              /* 8 */
                        << " " << dUind                	 /* 9 -11*/
                        << " " << dXi                	 /* 12 */
                        << " " << dAzimuth               /* 13 */
                        << " " << FMeanOut                	 /* 14-16 */
                        << std::endl;

                /* FIXME: check for parallel stuff ... */
                for (int i = 0; ppRes && ppRes[i]; i++) {
                        OH.Loadable()
                                << std::setw(8) << GetLabel()
                                << ":" << ppRes[i]->GetLabel()
                                << " " << ppRes[i]->pRes->Force()
                                << " " << ppRes[i]->pRes->Moment()
                                << std::endl;
                }
	}
}

SubVectorHandler&
CyclocopterPolimi::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	/* UNIFORM induced velocity */
	/* Trasporta della matrice di rotazione del rotore */
	RRotor = pCraft->GetRCurr()*RRot;
	RRotorTranspose = RRotor.Transpose();
	/* Forze nel sistema rotore */
	F = RRotorTranspose*Res.Force();
	if (!bFlagAverage) {
		/* filtro le forze */
		Uk = F;
		Yk = -Yk_1*a1 - Yk_2*a2 + Uk*b0 + Uk_1*b1 + Uk_2*b2;
		F = Yk;
		/* Forza nel piano normale all'asse di rotazione */
		doublereal dT = sqrt(F(2)*F(2) + F(3)*F(3));
		/* Velocità indotta: calcolata in base alla dT */
		doublereal dRho = dGetAirDensity(GetXCurr());
		dUindMean = dKappa*sqrt(dT/(2*dRho*dArea));
		/* Componenti della velocità indotta nel sistema
	 	* rotore */
		dUind = Zero3;
		if (dT > std::numeric_limits<doublereal>::epsilon()) {
			dUind(2) = dUindMean*F(2)/dT;
			dUind(3) = dUindMean*F(3)/dT;
		}
		dUind(1) = (1 - dWeight)*dUind(1) + dWeight*dUindPrev(1);
		dUind(2) = (1 - dWeight)*dUind(2) + dWeight*dUindPrev(2);
		dUind(3) = (1 - dWeight)*dUind(3) + dWeight*dUindPrev(3);

		dUindMean = sqrt(dUind(1)*dUind(1) + dUind(2)*dUind(2) + dUind(3)*dUind(3));
		/* angolo di cui è ruotata la trazione */
		dXi = atan2(F(3), F(2)) - M_PI/2.;
	}

	ResetForce();
	WorkVec.Resize(0);

	return WorkVec;
}

void
CyclocopterPolimi::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	/* colcolo la posizione azimutale della prima pala */
	if (bFlagIsFirstBlade) {
		Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));
		doublereal d1 = XRel(2);
		doublereal d2 = XRel(3);
		dAzimuth = atan2(d2, d1);
		bFlagIsFirstBlade = false;
	}

	/* Sole se deve fare l'output calcola anche il momento */
	if (fToBeOutput()) {
		Res.AddForces(F,M,X);
		InducedVelocity::AddForce(pEl, pNode, F, M, X);

	} else {
		Res.AddForce(F);
	}
}

Vec3
CyclocopterPolimi::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	Vec3 XRel(RRotorTranspose*(X - pRotor->GetXCurr()));

	doublereal d1 = XRel(2);
	doublereal d2 = XRel(3);

	/* dPsi0 non serve a nulla perchè uso l'angolo
	 * relativo: (dp-dXi)!!! */
	doublereal dpp = atan2(d2, d1);

	doublereal r = sqrt(d1*d1 + d2*d2)*cos(dpp - dXi);

	return RRotor*((dUind*(M_PI/2.))*cos((M_PI/2.)*(r/dRadius)));
}

/* CyclocopterPolimi - end */

#if 0 // KARI NOT IMPLEMENTED YET!

/* CyclocopterKARI - begin */

/*

From:

"A New VTOL UAV Cyclocopter with Cycloidal Blades System",

Chul Yong Yun, Illkyung Park,
Ho Yong Lee, Jai Sang Jung, In Seong Hwang,
Seung Jo Kim,
Sung Nam Jung

Presented at the American Helicopter Society 60th Annual Forum,
Baltimore, MD, June 7-10, 2004.
*/

class CyclocopterKARI
: virtual public Elem, public CyclocopterInflow {
public:
	CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~CyclocopterKARI(void);

	// Elaborazione stato interno dopo la convergenza
	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	// output; si assume che ogni tipo di elemento sappia,
	// attraverso l'OutputHandler, dove scrivere il proprio output
	virtual void Output(OutputHandler& OH) const;

	// Contributo al file di Restart
	virtual std::ostream& Restart(std::ostream& out) const;

	// assemblaggio residuo
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

#if 0
	// Somma alla trazione il contributo di un elemento
	virtual void
	AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X);
#endif

	// Restituisce ad un elemento la velocita' indotta
	// in base alla posizione azimuthale
	virtual Vec3 GetInducedVelocity(Elem::Type type,
		unsigned uLabel, unsigned uPnt, const Vec3& X) const;

	// *******PER IL SOLUTORE PARALLELO********
	// Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	// utile per l'assemblaggio della matrice di connessione fra i dofs
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	// ************************************************
};

CyclocopterKARI::CyclocopterKARI(unsigned int uL, const DofOwner* pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uL, 0),
CyclocopterInflow(uL, pDO),
pRotor(0),
RRot(::Eye3)
{
	if (!ReadRotorData(pDM, HP, uLabel, pCraft, RRot, pRotor)) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ppRes = ReadResSets(pDM, HP);

	SetOutputFlag(pDM->fReadOutput(HP, Elem::INDUCEDVELOCITY));
}

CyclocopterKARI::~CyclocopterKARI(void)
{
	NO_OP;
}

void
CyclocopterKARI::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
{
	NO_OP;
}

void
CyclocopterKARI::Output(OutputHandler& OH) const
{
	NO_OP;
}

std::ostream&
CyclocopterKARI::Restart(std::ostream& out) const
{
	return out << "# cyclocopter: not implemented yet" << std::endl;
}

SubVectorHandler&
CyclocopterKARI::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	return WorkVec;
}

#if 0
void
CyclocopterKARI::AddForce(const Elem *pEl, const StructNode *pNode, const Vec3& F, const Vec3& M, const Vec3& X)
{
	NO_OP;
}
#endif

Vec3
CyclocopterKARI::GetInducedVelocity(Elem::Type type,
	unsigned uLabel, unsigned uPnt, const Vec3& X) const
{
	return Zero3;
}

void
CyclocopterKARI::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(1);
	connectedNodes[0] = pCraft;
}

/* CyclocopterKARI - end */

#endif // KARI NOT IMPLEMENTED YET!

bool
mbdyn_cyclocopter_set(void)
{
	UserDefinedElemRead *rf;
	
	rf = new UDERead<CyclocopterNoInflow>;
	if (!SetUDE("cyclocopter" "no" "inflow", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter no inflow\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterUniform1D>;
	if (!SetUDE("cyclocopter" "uniform" "1D", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter uniform 1D\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterUniform2D>;
	if (!SetUDE("cyclocopter" "uniform" "2D", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter uniform 2D\""
			<< std::endl);

		return false;
	}

	rf = new UDERead<CyclocopterPolimi>;
	if (!SetUDE("cyclocopter" "Polimi", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter Polimi\""
			<< std::endl);

		return false;
	}

#if 0
	rf = new UDERead<CyclocopterKARI>;
	if (!SetUDE("cyclocopter" "KARI", rf)) {
		delete rf;

		silent_cerr("module-cyclocopter: "
			"unable to register \"cyclocopter KARI\""
			<< std::endl);

		return false;
	}
#endif

	return true;
}

#ifdef MBDYN_MODULE

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!mbdyn_cyclocopter_set()) {
		silent_cerr("cyclocopter: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // MBDYN_MODULE
