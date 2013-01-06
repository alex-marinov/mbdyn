/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* Discrete Control:
 * writes the output of selected dofs to a subprocess that generates the
 * input; then applies the input to selected dofs as an abstract force
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <cstdlib>

#include "discctrl.h"

/* DiscreteControlProcess - begin */

DiscreteControlProcess::~DiscreteControlProcess(void)
{
	NO_OP;
}

/* DiscreteControlProcess - end */

/* DiscreteControlARXProcess_Debug - begin */

int
DiscreteControlARXProcess_Debug::ReadMatrix(std::istream& In,
	doublereal* pd,
	unsigned int iNumRows,
	unsigned int iNumCols,
	unsigned int iNumSubMats,
	const char* sMatName)
{
	// Matrices alpha_i

	// for every k + 1 = 1->p (every matrix)
	for (unsigned int k = 0; k < iNumSubMats; k++) {
		int iShift = iNumCols*k;

		// for every row
		for (unsigned int i = 0; i < iNumRows; i++) {
			doublereal* pdTmp = pd + iShift + iNumSubMats*iNumCols*i;

			// for every column
			for (unsigned int j = 0; j < iNumCols; j++) {
				In >> pdTmp[j];
				if (!In) {
					silent_cerr("Error: unexpected end of stream while reading "
						<< sMatName << '_' << k + 1
						<< '(' << i + 1 << ',' << j + 1 << ')'
						<< std::endl);

					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
				}

				DEBUGLCOUT(MYDEBUG_INIT, sMatName << '_' << k + 1
					<< "(" << i << "," << j << ") = " << pdTmp[j]
					<< std::endl);
			}
		}
	}

	return 0;
}

DiscreteControlARXProcess_Debug::DiscreteControlARXProcess_Debug(integer iNumOut,
	integer iNumIn,
	integer iOrdA,
	integer iOrdB,
	const std::string& infile)
: iNumOutputs(iNumOut),
iNumInputs(iNumIn),
iOrderA(iOrdA),
iOrderB(iOrdB),
pdA(0),
pdY(0),
pdB(0),
pdU(0),
pdU0(0),
iRefA(0),
iRefB(0)
{
	/* The control model is handled as follows:
	 * the outputs and inputs are stored in two vector:
	 * Y = [y(k-1), ..., y(k-pa)], U = [u(k-1), ..., u(k-pb)]
	 * as soon as new y, u are available, the oldes values are replaced by
	 * the newest and the vector is used as a queue, i.e. it starts from
	 * a floating reference and goes back to the beginning when it's
	 * over.
	 * Mathices alpha, beta are stacked in two big matrtices:
	 * A = [alpha_1, ..., alpha_pa], B = [beta_1, ..., beta_pb]
	 * the input is calculated as A*Y + B*U
	 */

	// At present at least one input and one output are required
	ASSERT(iNumInputs > 0);
	ASSERT(iNumOutputs > 0);
	ASSERT(iOrderA > 0);
	ASSERT(iOrderB > 0);

	// Size of working area
	int iSize = iNumInputs*iNumOutputs*iOrderA	// Matrix A
		+ iNumInputs*iNumInputs*iOrderB		// Matrix B
		+ iNumOutputs*iOrderA			// Vector Y
		+ iNumInputs*iOrderB			// Vector U
		+ iNumInputs;				// Vector u(k)

	doublereal *pd = 0;
	SAFENEWARR(pd, doublereal, iSize);

	pdA = pd;
	pd += iNumInputs*iNumOutputs*iOrderA;

	pdB = pd;
	pd += iNumInputs*iNumInputs*iOrderB;

	pdY = pd;
	pd += iNumOutputs*iOrderA;

	pdU = pd;
	pd += iNumInputs*iOrderB;

	pdU0 = pd;

	for (int i = iSize; i-- > 0; ) {
		pdA[i] = 0.;
	}

	/* In the input file, for human readability matrices are written as:
	 * alpha_1,
	 * alpha_2,
	 * ...
	 * alpha_p,
	 * beta_1,
	 * ...
	 * beta_p
	 */

	std::ifstream In(infile.c_str());
	if (!In) {
		silent_cerr("DiscreteControlARXProcess_Debug: "
			"unable to open control data file \"" << infile << "\""
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	ReadMatrix(In, pdA, iNumInputs, iNumOutputs, iOrderA, "alpha");
	ReadMatrix(In, pdB, iNumInputs, iNumInputs, iOrderB, "beta");
}

DiscreteControlARXProcess_Debug::~DiscreteControlARXProcess_Debug(void)
{
	/* It must be non-null;
	 * matrices and arrays come from one allocation only */
	if (pdA != 0) {
		SAFEDELETEARR(pdA);
	}
}

void
DiscreteControlARXProcess_Debug::GetInput(std::vector<doublereal>& dIn)
{
	for (int i = iNumInputs; i-- > 0; ) {
		dIn[i] = pdU0[i];
	}
}

void
DiscreteControlARXProcess_Debug::PutOutput(const std::vector<doublereal>& dOut,
	const std::vector<doublereal>& dIn,
	const std::vector<doublereal>& /* dDesiredOut */ )
{
	/* At present dDesiredOutput is not used;
	 * it will be when the system is controlled
	 * and ID'd at the same time */

	// Moves reference backwards
	if (iRefA == 0) {
		iRefA = iOrderA;
	}
	iRefA--;

	doublereal* pdOff = pdY + iNumOutputs*iRefA;
	for (int i = iNumOutputs; i-- > 0; ) {
		pdOff[i] = dOut[i];
	}

	// Moves reference backwards
	if (iRefB == 0) {
		iRefB = iOrderB;
	}
	iRefB--;

	pdOff = pdU + iNumInputs*iRefB;
	for (int i = iNumInputs; i-- > 0; ) {
		/* dIn contiene la somma dell'ingresso esogeno e dell'ingresso
		 * di controllo, ovvero l'ingresso totale applicato al sistema;
		 * pdU0 contiene il solo ingresso generato dal controllore */
		pdOff[i] = dIn[i];
		pdU0[i] = 0.;
	}

	// Matrices are row-oriented
	doublereal* pdMatOff = pdA + iOrderA*iNumOutputs*iNumInputs;
	// Y shiftato di iRefA
	doublereal* pdVecOff = pdY + iRefA*iNumOutputs;
	for (int i = iNumInputs; i-- > 0; ) {
		for (int j = iRefA*iNumOutputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdY[j];
		}
		for (int j = (iOrderA-iRefA)*iNumOutputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdVecOff[j];
		}
	}

	pdMatOff = pdB + iOrderB*iNumInputs*iNumInputs;
	pdVecOff = pdU + iRefB*iNumInputs;
	for (int i = iNumInputs; i-- > 0; ) {
		for (int j = iRefB*iNumInputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdU[j];
		}
		for (int j = (iOrderB-iRefB)*iNumInputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdVecOff[j];
		}
	}
}

/* DiscreteControlProcess_Debug - end */


/* DiscreteIdentProcess_Debug - begin */

DiscreteIdentProcess_Debug::DiscreteIdentProcess_Debug(integer iNumOut,
	integer iNumIn,
	integer iOrdA,
	integer iOrdB,
	ForgettingFactor* pf,
	PersistentExcitation* px,
	unsigned f_proc,
	const std::string& sf)
: iNumOutputs(iNumOut), iNumInputs(iNumIn),
iOrderA(iOrdA), iOrderB(iOrdB), pId(0), pPx(px),
outfile(sf)
{
	ASSERT(pf != 0);
	ASSERT(px != 0);

	switch (f_proc) {
	case DISCPROC_ARX:
		SAFENEWWITHCONSTRUCTOR(pId, IdentARXProcess,
			IdentARXProcess(iNumOut, iNumIn, iOrdA, iOrdB, pf));
		break;

	case DISCPROC_ARMAX:
		SAFENEWWITHCONSTRUCTOR(pId, IdentARMAXProcess,
			IdentARMAXProcess(iNumOut, iNumIn, iOrdA, iOrdB, pf));
		break;

	default:
		silent_cerr("DiscreteIdentProcess_Debug: "
			"unknown model type " << f_proc << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// temporaneo ?!?
	if (!outfile.empty()) {
		out.open(outfile.c_str());
	}
}

DiscreteIdentProcess_Debug::~DiscreteIdentProcess_Debug(void)
{
	if (pId != 0) {
		SAFEDELETE(pId);
	}

	if (pPx != 0) {
		SAFEDELETE(pPx);
	}

	// temporaneo ?!?
	if (!outfile.empty()) {
		out.close();
	}
}

/* Returns the new control input values in array dIn */
void
DiscreteIdentProcess_Debug::GetInput(std::vector<doublereal>& dIn)
{
	/* azzera per sicurezza */
	for (integer i = iNumInputs; i-- > 0; ) {
		dIn[i] = 0.;
	}

	pPx->AddInput(&dIn[0]);
}

const unsigned long BUFSIZE = 102400;
static doublereal buf[BUFSIZE];

/* Sets the new measures (and eventually the input) */
void
DiscreteIdentProcess_Debug::PutOutput(const std::vector<doublereal>& dOut,
	const std::vector<doublereal>& dIn,
	const std::vector<doublereal>& /* dDesiredOut */ )
{
	pId->Update(&dOut[0], &dIn[0]);

	if (!outfile.empty()) {
		integer size = pId->iGetSize()*pId->iGetNumOutput();
		if (size*sizeof(doublereal) > BUFSIZE) {
			silent_cerr("buffer is too small" << std::endl);

		} else {
			pId->GetTheta(buf);
			for (integer i = 0; i < pId->iGetSize()*pId->iGetNumOutput(); i++) {
				out << std::setw(16) << buf[i];
			}
			out << std::setw(16) << pId->dGetForgettingFactor() << std::endl;
		}
	}
}

/* DiscreteIdentProcess_Debug - end */


/* DAC_Process_Debug - begin */

DAC_Process_Debug::DAC_Process_Debug(integer iNumOut,
	integer iNumIn,
	integer iOrdA,
	integer iOrdB,
	ForgettingFactor* pf,
	GPCDesigner* pd,
	PersistentExcitation* px,
	DriveCaller* pTrig,
	std::vector<DriveCaller *>& vDesOut,
	const std::string& sf,
	unsigned f_proc)
: iNumOutputs(iNumOut),
iNumInputs(iNumIn),
iOrderA(iOrdA),
iOrderB(iOrdB),
iOrderMd(0),
pdBase(0),
pdTheta(0),
pdA(0),
pdY(0),
pdB(0),
pdU(0),
f_proc(f_proc),
pdC(0),
pdE(0),
pdMd(0),
pdYd(0),
pdU0(0),
iRefA(0),
iRefB(0),
iRefMd(0),
pId(0), pCD(pd), pPx(px), Trigger(pTrig),
outfile(sf)
{
	/* The control model is handled as follows:
	 * the outputs and inputs are stored in two vector:
	 * Y = [y(k-1), ..., y(k-pa)], U = [u(k-1), ..., u(k-pb)]
	 * as soon as new y, u are available, the oldes values are replaced by
	 * the newest and the vector is used as a queue, i.e. it starts from
	 * a floating reference and goes back to the beginning when it's
	 * over.
	 * Mathices alpha, beta are stacked in two big matrtices:
	 * A = [alpha_1, ..., alpha_pa], B = [beta_1, ..., beta_pb]
	 * the input is calculated as A*Y + B*U
	 */

	/* At present at least one input and one output are required */
	ASSERT(iNumInputs > 0);
	ASSERT(iNumOutputs > 0);
	ASSERT(iOrderA > 0);
	ASSERT(iOrderB > 0);

	ASSERT(pCD != 0);
	ASSERT(pPx != 0);

	iOrderMd = pCD->iGetPredStep()-pCD->iGetPredHor();
	ASSERT(iOrderMd > 0);

	/* Size of working area */
	int iSize =
		iNumOutputs*iOrderA		// Vector Y
		+iNumInputs*iOrderB		// Vector U
		+iNumOutputs*iOrderMd		// Vector Yd
		+iNumInputs			// Vector u(k)
		+iNumOutputs*(iNumOutputs*iOrderA + iNumInputs*(iOrderB + 1));	// Theta

	if ((f_proc & DISCPROC_MA) == DISCPROC_MA) {
		iSize += iNumOutputs*iOrderA			// Vector E
			+ iNumOutputs*iNumOutputs*iOrderA;	// Theta e' piu' grande!
	}

	SAFENEWARR(pdBase, doublereal, iSize);

	pdY = pdBase;
	pdU = pdY + iNumOutputs*iOrderA;
	pdYd = pdU + iNumInputs*iOrderB;
	pdU0 = pdYd + iNumOutputs*iOrderMd;
	pdTheta = pdU0 + iNumInputs;

	if ((f_proc & DISCPROC_MA) == DISCPROC_MA) {
		pdE = pdU0 + iNumInputs;
		pdTheta = pdE + iNumOutputs*iOrderA;
	}

	for (int i = iSize; i-- > 0; ) {
		pdBase[i] = 0.;
	}

	/* In the input file, for human readability matrices are written as:
	 * alpha_1,
	 * alpha_2,
	 * ...
	 * alpha_p,
	 * beta_1,
	 * ...
	 * beta_p
	 */

	/* Si costruisce l'identificatore ARX */
	switch (f_proc) {
	case DISCPROC_ARX:
		SAFENEWWITHCONSTRUCTOR(pId, IdentARXProcess,
			IdentARXProcess(iNumOut, iNumIn, iOrdA, iOrdB, pf));
		break;

	case DISCPROC_ARMAX:
		SAFENEWWITHCONSTRUCTOR(pId, IdentARMAXProcess,
			IdentARMAXProcess(iNumOut, iNumIn, iOrdA, iOrdB, pf));
		break;

	default:
		silent_cerr("DAC_Process_Debug: "
			"unknown identification type " << f_proc
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!vDesOut.empty()) {
		vDesiredOut.resize(iNumOutputs);

		for (integer i = 0; i < iNumOutputs; i++) {
			vDesiredOut[i] = 0;
			SAFENEWWITHCONSTRUCTOR(vDesiredOut[i],
				DriveOwner,
				DriveOwner(vDesOut[i]));
		}
	}

	// temporaneo ?!?
	if (!outfile.empty()) {
		out.open(outfile.c_str());
		// mettere un test?
	}
}


DAC_Process_Debug::~DAC_Process_Debug(void)
{
	// matrices and arrays come from one allocation only

	if (pdBase != 0) {
		SAFEDELETEARR(pdBase);
	}

	if (pId != 0) {
		SAFEDELETE(pId);
	}

	if (pCD != 0) {
		SAFEDELETE(pCD);
	}

	if (pPx) {
		SAFEDELETE(pPx);
	}

	// can be null (no desired output
	for (std::vector<DriveOwner *>::iterator i = vDesiredOut.begin();
		i != vDesiredOut.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	// temporaneo ?!?
	if (!outfile.empty()) {
		out.close();
	}
}

void
DAC_Process_Debug::GetInput(std::vector<doublereal>& dIn)
{
	// control input
	for (int i = iNumInputs; i-- > 0; ) {
		dIn[i] = pdU0[i];
	}

	// persistent excitation
	pPx->AddInput(&dIn[0]);
}

void
DAC_Process_Debug::PutOutput(const std::vector<doublereal>& dOut,
	const std::vector<doublereal>& dIn,
	const std::vector<doublereal>& dDesiredOut)
{
	/* At present dDesiredOutput is not used;
	 * it will be when the system is controlled
	 * and ID'd at the same time */

	ASSERT(pId != 0);
	ASSERT(pCD != 0);

	pId->Update(&dOut[0], &dIn[0]);

	if (!outfile.empty()) {
		pId->GetTheta(pdTheta);

		for (integer i = 0; i < pId->iGetSize()*pId->iGetNumOutput(); i++) {
			out << std::setw(16) << pdTheta[i];
		}
		out << std::setw(16) << pId->dGetForgettingFactor() << std::endl;
	}

	if (Trigger.dGet()) {
		pCD->DesignControl(pdTheta, &pdA, &pdB, &pdMd, &pdC);
	}

	// Moves reference backwards
	if (iRefA == 0) {
		iRefA = iOrderA;
	}
	iRefA--;

	doublereal* pdOff = pdY + iNumOutputs*iRefA;
	for (int i = iNumOutputs; i-- > 0; ) {
		pdOff[i] = dOut[i];
	}

	if ((f_proc & DISCPROC_MA) == DISCPROC_MA) {
		pdOff = pdE + iNumOutputs*iRefA;
		pId->GetErr(pdOff);
	}

	/* Moves reference backwards */
	if (iRefB == 0) {
		iRefB = iOrderB;
	}
	iRefB--;

	pdOff = pdU + iNumInputs*iRefB;
	for (int i = iNumInputs; i-- > 0; ) {
		/* dIn contiene la somma dell'ingresso esogeno e dell'ingresso
		 * di controllo, ovvero l'ingresso totale applicato al sistema;
		 * pdU0 contiene il solo ingresso generato dal controllore */
		pdOff[i] = dIn[i];
		pdU0[i] = 0.;
	}

	// Moves reference backwards
	if (!vDesiredOut.empty()) {
		if (iRefMd == 0) {
			iRefMd = iOrderMd;
		}
		iRefMd--;

		doublereal* pdOff = pdYd + iNumOutputs*iRefMd;
		for (int i = iNumOutputs; i-- > 0; ) {
			pdOff[i] = vDesiredOut[i]->dGet();
		}

		if (!dDesiredOut.empty()) {
			for (int i = iNumOutputs; i-- > 0; ) {
				pdOff[i] += dDesiredOut[i];
			}
		}
	}

	if (pdB == 0) {
		return;
	}

	// Matrices are row-oriented
	doublereal* pdMatOff = pdA + iOrderA*iNumOutputs*iNumInputs;
	// Y shiftato di iRefA
	doublereal* pdVecOff = pdY + iRefA*iNumOutputs;
	for (int i = iNumInputs; i-- > 0; ) {
		for (int j = iRefA*iNumOutputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdY[j];
		}
		for (int j = (iOrderA-iRefA)*iNumOutputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdVecOff[j];
		}
	}

	pdMatOff = pdB + iOrderB*iNumInputs*iNumInputs;
	pdVecOff = pdU + iRefB*iNumInputs;
	for (int i = iNumInputs; i-- > 0; ) {
		for (int j = iRefB*iNumInputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdU[j];
		}
		for (int j = (iOrderB - iRefB)*iNumInputs; j-- > 0; ) {
			pdMatOff--;
			pdU0[i] += pdMatOff[0]*pdVecOff[j];
		}
	}

	if ((f_proc & DISCPROC_MA) == DISCPROC_MA) {
		pdMatOff = pdC + iOrderA*iNumOutputs*iNumInputs;
		// E shiftato di iRefA
		pdVecOff = pdE + iRefA*iNumOutputs;
		for (int i = iNumInputs; i-- > 0; ) {
			for (int j = iRefA*iNumOutputs; j-- > 0; ) {
				pdMatOff--;
				pdU0[i] += pdMatOff[0]*pdE[j];
			}
			for (int j = (iOrderA - iRefA)*iNumOutputs; j-- > 0; ) {
				pdMatOff--;
				pdU0[i] += pdMatOff[0]*pdVecOff[j];
			}
		}
	}

	if (!vDesiredOut.empty()) {
		pdMatOff = pdMd + iOrderMd*iNumOutputs*iNumInputs;
		// Yd shiftato di iRefMd
		pdVecOff = pdYd + iRefMd*iNumOutputs;
		for (int i = iNumInputs; i-- > 0; ) {
			for (int j = iRefMd*iNumOutputs; j-- > 0; ) {
				pdMatOff--;
				pdU0[i] += pdMatOff[0]*pdYd[j];
			}
			for (int j = (iOrderMd - iRefMd)*iNumOutputs; j-- > 0; ) {
				pdMatOff--;
				pdU0[i] += pdMatOff[0]*pdVecOff[j];
			}
		}
	}
}

/* DAC_Process_Debug - end */

/* DiscreteControlElem - begin */

DiscreteControlElem::DiscreteControlElem(unsigned int uL,
	const DofOwner* pDO,
	integer iNumOut,
	std::vector<ScalarValue *>& vOut,
	std::vector<DriveCaller *>& vOutSF,
	integer iNumIn,
	ScalarDof* pIn,
	DiscreteControlProcess* p,
	integer iNIt,
	flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pDCP(p),
bNewStep(true),
iNumIter(iNIt),
iCurrIter(0),
iNumOutputs(iNumOut),
vOutputs(vOut),
dOut(iNumOutputs),
iNumInputs(iNumIn),
pInputs(pIn),
dIn(iNumInputs)
{
	/* The inputs should be all scalar dofs; an abtract force is
	 * applied to every node;
	 * the outputs can be of every kind */

	/* Allocs and resets memory for input and output current values */
	ASSERT(iNumInputs > 0);
	ASSERT(iNumOutputs > 0);

	vOutScaleFact.resize(iNumOutputs);

	for (int i = iNumOutputs; i-- > 0; ) {
		vOutScaleFact[i] = 0;
		SAFENEWWITHCONSTRUCTOR(vOutScaleFact[i],
			DriveOwner,
			DriveOwner(vOutSF[i]));
	}

	for (std::vector<doublereal>::iterator i = dIn.begin();
		i != dIn.end(); ++i)
	{
		*i = 0.;
	}

	for (std::vector<doublereal>::iterator i = dOut.begin();
		i != dOut.end(); ++i)
	{
		*i = 0.;
	}
}

DiscreteControlElem::~DiscreteControlElem(void)
{
	for (std::vector<DriveOwner *>::iterator i = vOutScaleFact.begin();
		i != vOutScaleFact.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	for (std::vector<ScalarValue *>::iterator i = vOutputs.begin();
		i != vOutputs.end(); ++i)
	{
		SAFEDELETE(*i);
	}

	// Allocated by DataManager
	SAFEDELETEARR(pInputs);

	// Must destroy the other processes too
	SAFEDELETE(pDCP);
}

Electric::Type
DiscreteControlElem::GetElectricType(void) const
{
	return Electric::DISCRETECONTROL;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
DiscreteControlElem::Restart(std::ostream& out) const
{
	out << "  electric: " << GetLabel() << ", discrete control; "
		"# to be implemented" << std::endl;
	return out;
}

static const std::vector<doublereal> dEmptyDesiredOut;

void
DiscreteControlElem::AfterConvergence(const VectorHandler& X,
	const VectorHandler& XP)
{
	// Sets the flag for a new step
	if (++iCurrIter == iNumIter) {
		iCurrIter = 0;

		ASSERT(bNewStep == false);
		bNewStep = true;

		// Gets output from solution
		for (int iCnt = iNumOutputs; iCnt-- > 0; ) {
			dOut[iCnt] = vOutScaleFact[iCnt]->dGet()*vOutputs[iCnt]->dGetValue();
		}

		// Gets input from solution
		for (int iCnt = iNumInputs; iCnt-- > 0; ) {
			dIn[iCnt] = pInputs[iCnt].dGetValue();
		}

		// Puts measures to control subprocess
		// DCP knows the length of dOut and dIn
		pDCP->PutOutput(dOut, dIn, dEmptyDesiredOut);
	}

	// Add extra output as needed
}

/* assemblaggio residuo */
SubVectorHandler&
DiscreteControlElem::AssRes(SubVectorHandler& WorkVec,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	WorkVec.ResizeReset(iNumInputs);

	// At first iteration gets the input from the control subprocess
	if (bNewStep) {
		// Gets the inputs from subprocess
		// DCP knows the length of dIn
		pDCP->GetInput(dIn);

		// resets the flag
		bNewStep = false;
	}

	// Sets the parameters
	for (int iCnt = iNumInputs; iCnt-- > 0; ) {
		// Equivalent to an abstract force
		WorkVec.PutItem(iCnt + 1, pInputs[iCnt].pNode->iGetFirstRowIndex() + 1, dIn[iCnt]);
	}

	return WorkVec;
}

/* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
unsigned int
DiscreteControlElem::iGetNumDof(void) const
{
      return 0;
}

/* esegue operazioni sui dof di proprieta' dell'elemento */
DofOrder::Order
DiscreteControlElem::GetDofType(unsigned int /* i */ ) const
{
	ASSERTMSG(0, "You shouldn't have called this function");
	return DofOrder::UNKNOWN;
}

/* Dimensioni del workspace */
void
DiscreteControlElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = iNumInputs;
	*piNumCols = 1;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler&
DiscreteControlElem::AssJac(VariableSubMatrixHandler& WorkMat,
	doublereal /* dCoef */ ,
	const VectorHandler& /* XCurr */ ,
	const VectorHandler& /* XPrimeCurr */ )
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}

unsigned int
DiscreteControlElem::iGetNumPrivData(void) const
{
	return iNumInputs;
}

unsigned int
DiscreteControlElem::iGetPrivDataIdx(const char *s) const
{
	/*
	 * legal values:
	 *
	 * "u[<idx>]"	with <idx> ::= positive integer between 1 and iNumInputs
	 */
	switch (s[0]) {
	case 'u':
		s++;
		break;

	default:
		return 0;
	}

	if (s[0] != '[') {
		return 0;
	}
	s++;

	if (s[0] == '-') {
		return 0;
	}

	char *next;
	errno = 0;
	long idx = std::strtol(s, &next, 10);
	int save_errno = errno;
	if (next == s || next[0] != ']') {
		return 0;
	}

	if (save_errno == ERANGE) {
		silent_cerr("DiscreteControlElem(" << GetLabel() << "): "
			"warning, private data index " << std::string(s, next - s)
			<< " overflows" << std::endl);
		return 0;
	}

	if (idx > iNumInputs) {
		return 0;
	}

	return idx;
}

doublereal
DiscreteControlElem::dGetPrivData(unsigned int i) const
{
	if (i < 1 || i > (unsigned int)iNumInputs) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return dIn[i - 1];
}

/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
 * utile per l'assemblaggio della matrice di connessione fra i dofs */
void
DiscreteControlElem::GetConnectedNodes(
	std::vector<const Node *>& connectedNodes) const {
	connectedNodes.resize(iNumInputs + iNumOutputs);
	for (int i = 0; i < iNumInputs; i++) {
		connectedNodes[i] = pInputs[i].pNode;
	}

	unsigned cnt = unsigned(iNumInputs);
	for (int i = 0; i < iNumOutputs; i++) {
		ScalarDofValue *psdv = dynamic_cast<ScalarDofValue *>(vOutputs[i]);
		if (psdv) {
			connectedNodes[cnt] = psdv->pNode;
			cnt++;
		}
	}

	ASSERT(cnt <= connectedNodes.size());

	if (cnt != connectedNodes.size()) {
		connectedNodes.resize(cnt);
	}
}

/* DiscreteControlElem - end */

