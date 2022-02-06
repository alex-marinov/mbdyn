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

/* Discrete Control:
 * writes the output of selected dofs to a subprocess that generates the
 * input; then applies the input to selected dofs as an abstract force
 */

#ifndef DISCCTRL_H
#define DISCCTRL_H

#include <iostream>
#include <fstream>

#include "elec.h"
#include "scalarvalue.h"

#include "id.h"
#include "gpc.h"
#include "px.h"

/* DiscreteControlProcess - begin */
/* This class provides input and output writing to and from the discrete
 * control MBDyn element and the external process that actually computes
 * the control inputs from the control outputs.
 *
 * At present, for developing and debugging purposes, the class itself
 * generates the output.
 *
 * The discrete system is thought in the form:
 *
 * y(k) = alpha_1*y(k-1)+...+alpha_p*y(k-p)
 *        +beta_0*u(k)+...+beta_p*u(k-p)
 *
 * where y(i) is the m x 1 vector of the measures at step i,
 *       alpha_v is the v-th m x m system matrix,
 *       beta_w is the w-th m x r input matrix, and
 *       u(j) is the r x 1 input vector at step j.
 *
 * The controller computes the control input at step k based on the inputs
 * and outputs at previous time steps [and on the desired value
 * for the outputs] by means of the control matrices a_v, b_w
 * [and m_d] :
 *
 * u(k) = a_1*y(k-1)+...+a_p*y(k-p)+b_1*u(k-1)+...+b_p*u(k-p)
 *
 * [ ... + m_s*yd(k+s)+...+m_0*y_d(k)
 * (Since usually the future values of y_d are unknown, the actual oves are
 * considered, that is y_d is shifted backwards by s places) [B]
 *
 * The matrices are organised as follows:
 *
 * A = [a_1,...,a_p]
 * B = [b_1,...,b_p]
 *
 * [ M = [m_s,...,m_0] ]
 *
 * The vectors are stacked as follows:
 *
 * Y = [y(k-1),...,y(k-p)]
 * U = [u(k-1),...,u(k-p)]
 *
 * as soon as a new measure set and a new control set is available,
 * it replaces the oldest one, that is not required any more, and the
 * computation of the input uses the old inputs and outputs vectors
 * as queues, maintaining a reference to the latest ones. This minimizes the
 * operations required for the computation of control inputs for large systems
 */

class DiscreteControlProcess {
public:
	enum {
		DISCPROC_UNKNOWN	= -1,

		DISCPROC_AR		= 0x1U,
		DISCPROC_MA		= 0x2U,
		DISCPROC_X			= 0x4U,

		DISCPROC_ARX		= (DISCPROC_AR | DISCPROC_X),
		DISCPROC_ARMA		= (DISCPROC_AR | DISCPROC_MA),
		DISCPROC_ARMAX		= (DISCPROC_AR | DISCPROC_MA | DISCPROC_X),

		DISCPROC_LAST
	};

public:
	virtual ~DiscreteControlProcess(void);

	// Returns the new control input values in array dIn
	virtual void GetInput(std::vector<doublereal>& dIn) = 0;

	// Sets the new measures (and the input)
	virtual void PutOutput(const std::vector<doublereal>& dOut,
		const std::vector<doublereal>& dIn,
		const std::vector<doublereal>& dDesiredOut) = 0;
};

/* DiscreteControlProcess - end */


/* DiscreteControlARXProcess_Debug - begin */

class DiscreteControlARXProcess_Debug : public DiscreteControlProcess {
protected:
	integer iNumOutputs;	// Number of outputs (measures)
	integer iNumInputs;	// Number of inputs (forces)
	integer iOrderA;	// Order of the system (p)
	integer iOrderB;	// Order of the input (p)
	doublereal* pdA;	// Stack of matrices a_c (r x (p*m))
	doublereal* pdY;	// Stack of output vectors at previous times (p*m)
	doublereal* pdB;	// Stack of matrices b_c (r x (p*r) )
	doublereal* pdU;	// Stack of input vectors at previous times (p*r)
	doublereal* pdU0;	// Current input vector (r)
	integer iRefA;		// Current position of most recent output
	integer iRefB;		// Current position of most recent input

	const std::string infile;

	// Reads the control matrices
	int
	ReadMatrix(std::istream& In, doublereal* pd,
		unsigned int iRows, unsigned int iCols,
		unsigned int iNumSubMats,
		const char* sMatName);

public:
	DiscreteControlARXProcess_Debug(integer iNumOut, integer iNumIn,
		integer iOrdA, integer iOrdB, const std::string& infile);
	virtual ~DiscreteControlARXProcess_Debug(void);

	// Returns the new control input values in array dIn
	void GetInput(std::vector<doublereal>& dIn);

	// Sets the new measures (and the input)
	void
	PutOutput(const std::vector<doublereal>& dOut,
		const std::vector<doublereal>& dIn,
		const std::vector<doublereal>& dDesiredOut);
};

/* DiscreteControlARXProcess_Debug - end */


/* DiscreteIdentProcess_Debug - begin */

class DiscreteIdentProcess_Debug : public DiscreteControlProcess {
protected:
	integer iNumOutputs;		/* Number of outputs (measures) */
	integer iNumInputs;		/* Number of inputs (forces) */
	integer iOrderA;		/* Order of the system (p) */
	integer iOrderB;		/* Order of the input (p) */

	IdentProcess*  pId;		/* Identifier */
	PersistentExcitation* pPx;	/* Excitation */

	// provvisorio?!?
	const std::string outfile;
	std::ofstream out;

public:
	DiscreteIdentProcess_Debug(integer iNumOut, integer iNumIn,
		integer iOrdA, integer iOrdB,
		ForgettingFactor* pf,
		PersistentExcitation* px,
		unsigned f_proc,
		const std::string& sf);
	virtual ~DiscreteIdentProcess_Debug(void);

	// Returns the new control input values in array dIn
	void GetInput(std::vector<doublereal>& dIn);

	// Sets the new measures (and the input)
	void PutOutput(const std::vector<doublereal>& dOut,
		const std::vector<doublereal>& dIn,
		const std::vector<doublereal>& dDesiredOut);
};

/* DiscreteIdentProcess_Debug - end */


/* DAC_Process_Debug - begin */

class DAC_Process_Debug : public DiscreteControlProcess {
protected:
	integer iNumOutputs;	// Number of outputs (measures)
	integer iNumInputs;	// Number of inputs (forces)
	integer iOrderA;	// Order of the system (pa)
	integer iOrderB;	// Order of the input (pb)
	integer iOrderMd;	// Order of the desired output
	doublereal* pdBase;
	doublereal* pdTheta;	// Predicted matrix
	doublereal* pdA;	// Stack of matrices a_c (r x (pa*m))
	doublereal* pdY;	// Stack of output vectors at previous times (p*m)
	doublereal* pdB;	// Stack of matrices b_c (r x (pb*r))
	doublereal* pdU;	// Stack of input vectors at previous times (p*r)
	unsigned f_proc;
	doublereal* pdC;	// Stack of matrices c_c (r x (pa*m))
	doublereal* pdE;	// Stack of error vectors at previous times
	doublereal* pdMd;	// Stack of matrices m_c (r x (pa*?))
	doublereal* pdYd;	// Stack of desired output vectors at following times
	doublereal* pdU0;	// Current input vector (r)
	integer iRefA;		// Current position of most recent output
	integer iRefB;		// Current position of most recent input
	integer iRefMd;		// Current position of most recent desired output

	IdentProcess*  pId;	// Identifier
	GPCDesigner* pCD;	// control designer
	PersistentExcitation* pPx;	// excitation
	DriveOwner Trigger;	// control trigger

	std::vector<DriveOwner*> vDesiredOut;

	// provvisorio?!?
	const std::string outfile;
	std::ofstream out;

public:
	DAC_Process_Debug(integer iNumOut, integer iNumIn,
		integer iOrdA, integer iOrdB,
		ForgettingFactor* pf,
		GPCDesigner* pd,
		PersistentExcitation* px,
		DriveCaller* PTrig,
		std::vector<DriveCaller *>& vDesOut,
		const std::string& sf,
		unsigned f_proc);
	virtual ~DAC_Process_Debug(void);

	// Returns the new control input values in array dIn
	void GetInput(std::vector<doublereal>& dIn);

	// Sets the new measures (and the input)
	void PutOutput(const std::vector<doublereal>& dOut,
		const std::vector<doublereal>& dIn,
		const std::vector<doublereal>& dDesiredOut);
};

/* DAC_Process_Debug - end */


/* DiscreteControlElem - begin */

class DiscreteControlElem : virtual public Elem, public Electric {
protected:
	DiscreteControlProcess* pDCP;
	bool bNewStep;		// Decides whether to read data from control process
	integer iNumIter;	// The control acts at iNumIter*dt time steps
	integer iCurrIter;

	integer iNumOutputs;
	std::vector<ScalarValue *> vOutputs;
	std::vector<DriveOwner*> vOutScaleFact;
	std::vector<doublereal> dOut;

	integer iNumInputs;
	ScalarDof* pInputs;
	std::vector<doublereal> dIn;

public:
	DiscreteControlElem(unsigned int uL, const DofOwner* pDO,
		integer iNumOut,
		std::vector<ScalarValue *>& vOut,
		std::vector<DriveCaller *>& vOutSF,
		integer iNumIn,
		ScalarDof* ppIn,
		DiscreteControlProcess* p,
		integer iNIt,
		flag fOut);
	virtual ~DiscreteControlElem(void);

	virtual Electric::Type GetElectricType(void) const;

	/* Scrive il contributo dell'elemento al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	virtual void
	AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

	/* ritorna il numero di Dofs per gli elementi che sono anche DofOwners */
	virtual unsigned int iGetNumDof(void) const;

	/* esegue operazioni sui dof di proprieta' dell'elemento */
	virtual DofOrder::Order GetDofType(unsigned int /* i */ ) const;

	/* Dimensioni del workspace */
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;

	/* assemblaggio jacobiano */
	virtual VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ );

	/* assemblaggio residuo */
	virtual SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);

	/* Dati privati */
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;   

	/* *******PER IL SOLUTORE PARALLELO******** */
	/* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
	 * utile per l'assemblaggio della matrice di connessione fra i dofs */
	virtual void
	GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	/* ************************************************ */

	/* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

/* DiscreteControlElem - end */

#endif /* DISCCTRL_H */

