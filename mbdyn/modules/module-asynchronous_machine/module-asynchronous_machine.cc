/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 AUTHOR: Reinhard Resch <r.resch@secop.com>
        Copyright (C) 2011(-2014) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

/**
 * \mainpage
 * This code is a simple example on how to implement a user defined element for the free multibody software <A HREF="http://www.aero.polimi.it/mbdyn/">MBDyn</A>.
 * Since no documentation exists for user defined elements in MBDyn at the time of writing I have written this document in the hope that it will be helpful for others who want to write their own user defined elements for MBDyn.<BR>
 *
 * The purpose of this code is the simulation of an asynchronous electrical machine as a part of a multibody model.
 * With this software the coupling between the mechanical and the electrical differential equations can be considered.
 * The theoretical background is based on the book <A HREF="http://books.google.at/books?id=z6B3lcxDgkkC&printsec=frontcover&dq=Maschinendynamik++%26+dresig&hl=de&ei=SwvETYi8NJSw4Ab_luX3BA&sa=X&oi=book_result&ct=result&resnum=1&ved=0CDIQ6AEwAA#v=onepage&q&f=false">
 * Maschinendynamik Springer 2009, Dresig, Holzweißig.</A><BR>
 * In the formulas in this book the stator resistance is neglected and the free run slip is assumed to be zero.
 * Also small perturbations around the mean angular velocity are assumed.<BR><BR>
 * The input parameters for the simulation are:<BR>
 *  Breakdown torque \f$M_K\f$.<BR>
 *  Slip at the breakdown torque \f$s_K\f$.<BR>
 *  Synchronous angular velocity \f$\Omega_S=2\,\pi\,\frac{f}{p}\f$ (might be a function of the time in case of an frequency inverter).<BR>
 *  A flag that determines whether the power supply is turned on or off \f$\gamma(t)\f$.<BR>
 *  Optional the initial value of the motor torque \f$M\f$.<BR>
 *  Optional the initial value of the first derivative of the motor torque \f$\dot{M}\f$.<BR>
 *
 *  The sense of rotation is determined by the sign of \f$\Omega_S\f$.
 *  If the sign is positive, the sense of rotation is also positive in the reference frame of the stator node.
 *
 *  This code is implemented as an user defined element for <A HREF="http://www.aero.polimi.it/mbdyn/">MBDyn</A>.
 *  The element is connected to two structural nodes. The rotor node and the stator node.
 *  The axis of rotation is assumed to be axis three in the reference frame of the stator node.
 *  It is assumed that the axis of the rotor node is parallel to the axis of the stator node.
 *  The torque is imposed to the structural nodes in direction of axis three in the reference frame of the stator node.
 *  The synchronous angular velocity can be specified by means of a drive caller. This allows to simulate a frequency inverter.
 */

#include "mbconfig.h"           // This goes first in every *.c,*.cc file

#include <cassert>
#include <cstdio>
#include <cmath>
#include <cfloat>
#include <iostream>
#include <iomanip>
#include <limits>

#include "dataman.h"
#include "userelem.h"
#include "module-asynchronous_machine.h"

class asynchronous_machine
: virtual public Elem, public UserDefinedElem
{
public:
	asynchronous_machine(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~asynchronous_machine(void);
	virtual void Output(OutputHandler& OH) const;
	virtual unsigned int iGetNumDof(void) const;
	virtual DofOrder::Order GetDofType(unsigned int i) const;
	virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
	virtual unsigned int iGetNumPrivData(void) const;
	virtual unsigned int iGetPrivDataIdx(const char *s) const;
	virtual doublereal dGetPrivData(unsigned int i) const;
	virtual void SetInitialValue(VectorHandler& X);
	virtual void Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
	virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	virtual int iGetNumConnectedNodes(void) const;
	virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	virtual void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
			SimulationEntity::Hints *ph);
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	virtual VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	virtual SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

private:
	bool IsMotorOn(void) const;

private:
	const StructNode* m_pRotorNode;	/**< node of the rotor where the torque \f$\boldsymbol{f}_3\f$ is imposed and the angular velocity \f$\dot{\boldsymbol{\varphi}}_1\f$ is determined */
	const StructNode* m_pStatorNode;	/**< node of the stator where the torque \f$\boldsymbol{f}_4\f$ is applied and the angular velocity \f$\dot{\boldsymbol{\varphi}}_2\f$ is determined (axis 3 of the stator node is assumed to be the axis of rotation) */
	doublereal m_MK;				/**< breakdown torque \f$M_K\f$	*/
	doublereal m_sK;				/**< slip at the breakdown torque \f$s_K\f$	*/
	DriveOwner m_OmegaS;		/**< drive that returns the synchronous angular velocity \f$\Omega_S=2\,\pi\,\frac{f}{p}\f$ (might be a function of the time in case of an frequency inverter)	*/
	DriveOwner m_MotorOn;		/**< drive that returns whether the power supply is turned on or off	\f$\gamma(t)\f$ */
	doublereal m_M;					/**< holds the motor torque \f$M\f$ after convergence for convenience	*/
	doublereal m_dM_dt;				/**< holds the first derivative of the motor torque \f$\dot{M}\f$ after convergence	*/
	doublereal m_dM_dt2;			/**< holds the second derivative of the motor torque \f$\ddot{M}\f$ after convergence	*/
	doublereal m_omega;				/**< holds the angular velocity \f$\dot{\varphi}\f$ of the rotor relative to the stator around axis 3 of the stator node after convergence	*/
	doublereal m_domega_dt;			/**< holds the angular acceleration \f$\ddot{\varphi}\f$ of the rotor node relative to the stator node after convergence	*/
    static const doublereal sm_SingTol;
};

const doublereal asynchronous_machine::sm_SingTol = std::pow(std::numeric_limits<doublereal>::epsilon(), 0.9);

/**
 * \param uLabel the label assigned to this element in the input file
 * \param pDO
 * \param pDM pointer to the data manager (needed to read structural nodes for example)
 * \param HP reference to the math parser (needed to read from the input file)
 * \brief Constructs the element from the information in the input file.<BR>
 *
 * A description of the exact input file syntax can be retrieved by adding the following statement to the MBDyn input file:<BR>
 * user defined: 1, asynchronous_machine, help;
 */
asynchronous_machine::asynchronous_machine(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: 	Elem(uLabel, flag(0)),
	UserDefinedElem(uLabel, pDO),
	m_pRotorNode(0),
	m_pStatorNode(0),
	m_MK(0.),
	m_sK(0.),
	m_OmegaS(0),
	m_MotorOn(0),
	m_M(0.),
	m_dM_dt(0.),
	m_dM_dt2(0.),
	m_omega(0.),
	m_domega_dt(0.)
{
	if (HP.IsKeyWord("help")) {
		silent_cout(
			"\n"
			"Module: 	asynchronous_machine\n"
			"\n"
			"	This module implements an asynchronous electric machine\n"
			"   according to\n"
			"\n"
			"	Maschinendynamik\n"
			"	Hans Dresig, Franz Holzweißig\n"
			"	Springer, 2009\n"
			"	page 58 equation 1.122 and 1.123\n"
			"\n"
			"Syntax:\n"
			"	asynchronous_machine,\n"
			"		rotor, (label) <rotor_node>,\n"
			"		stator, (label) <stator_node>,\n"
			"		MK, (real) <MK>,\n"
			"		sK, (real) <sK>,\n"
			"		OmegaS, (DriveCaller) <omega_s>\n"
			"		[ , motor on , (DriveCaller) <motor_on> ]\n"
			"		[ , M0 , (real) <M0> ]\n"
			"		[ , MP0 , (real) <MP0> ]\n"
			"\n"
			"       MK ... breakdown torque ( MK > 0 )\n"
			"       sK ... breakdown slip ( sK > 0 )\n"
			"       OmegaS = 2 * pi * f / p * sense_of_rotation ... synchronous angular velocity\n"
			"       motor_on ... power supply is turned off when 0, on otherwise\n"
			"       M0 ... initial torque ( default 0 )\n"
			"       MP0 ... initial torque derivative ( default 0 )\n"
			"\n"
			"       The axis of rotation is assumed to be axis 3 of the reference frame of the stator node.\n"
			<< std::endl);

		if (!HP.IsArg()) {
			// Exit quietly if nothing else is provided
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	if (!HP.IsKeyWord("rotor")) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): keyword \"rotor\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pRotorNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	if (!m_pRotorNode) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HP.IsKeyWord("stator")) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): keyword \"stator\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	m_pStatorNode = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

	if (!m_pStatorNode) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): structural node expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!HP.IsKeyWord("MK")) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): keyword \"MK\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_MK = HP.GetReal();

	if (!HP.IsKeyWord("sK")) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): keyword \"sK\" expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_sK = HP.GetReal();

	if (!HP.IsKeyWord("OmegaS")) {
		silent_cerr("asynchronous_machine(" << GetLabel() << "): keyword \"OmegaS\" expected" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	m_OmegaS.Set(HP.GetDriveCaller());

	if (HP.IsKeyWord("motor" "on") || HP.IsKeyWord("motor_on")) {
		m_MotorOn.Set(HP.GetDriveCaller());
	} else {
		m_MotorOn.Set(new OneDriveCaller());
	}

	if (HP.IsKeyWord("M0")) {
		m_M = HP.GetReal();
	}

	if (HP.IsKeyWord("MP0")) {
		m_dM_dt = HP.GetReal();
	}

	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

	Vec3 e3(m_pStatorNode->GetRCurr().GetCol(3));
	const Vec3& omega1 = m_pRotorNode->GetWCurr();
	const Vec3& omega2 = m_pStatorNode->GetWCurr();

	m_omega = e3.Dot(omega1 - omega2);

	const doublereal OmegaS = m_OmegaS.dGet();

	pDM->GetLogFile() << "asynchronous machine: "
		<< uLabel << " "
		<< m_pRotorNode->GetLabel() << " "
		<< m_pStatorNode->GetLabel() << " "
		<< m_MK << " "
		<< m_sK << " "
		<< OmegaS << " "
		<< IsMotorOn() << " "
		<< m_M << " "
		<< m_dM_dt << " "
		<< m_omega
		<< std::endl;

	// Note: The ODE for the asynchronous machine is defined for positive sign of OmegaS only!
	// At user level the signs of M, dM/dt and d^2M/dt^2 are defined to be positive
	// in positive coordinate system direction in the reference frame of the stator.
	m_M *= copysign(1., OmegaS);
	m_dM_dt *= copysign(1., OmegaS);
	m_dM_dt2 *= copysign(1., OmegaS);
}

/**
 * @param X vector of global degrees of freedom \f$\boldsymbol{y}\f$ after convergence at the current time step.
 * @param XP derivative of the global degrees of freedom \f$\dot{\boldsymbol{y}}\f$ after convergence at the current time step.<BR>
 * \brief This member function is called after each iteration in the nonlinear solution phase.
 *
 * In our element this member function saves the current state of the private degrees of freedom.
 * This is needed for the implementation of the Output() and the dGetPrivData() member functions.
 * If the private degrees of freedom are needed just for Output(), the code that saves it's state should be moved to AfterPredict() since this member function is called only once per time step.
 */
void
asynchronous_machine::Update(const VectorHandler& X,const VectorHandler& XP)
{
	const integer iFirstIndex = iGetFirstIndex();

	const integer intTorqueDerivativeRowIndex 	 = iFirstIndex + 1;
	const integer intTorqueRowIndex 			 = iFirstIndex + 2;
	const integer intOmegaRowIndex 				 = iFirstIndex + 3;

	// save the current state needed by dGetPrivData() and Output()
	m_M 		= X(intTorqueRowIndex);
	m_dM_dt		= X(intTorqueDerivativeRowIndex);
	m_dM_dt2	= XP(intTorqueDerivativeRowIndex);
	m_omega 	= X(intOmegaRowIndex);
	m_domega_dt = XP(intOmegaRowIndex);
}

void asynchronous_machine::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	Update(X, XP);
}

/**
 * @return Returns true if the motor is powered on.<BR>
 * \brief If a drive caller has been specified in the input file the value returned by the drive caller is used to determine if the motor is powered on.
 * If no drive caller has been specified in the input file it is assumed that the motor is always powered on.
 */
bool
asynchronous_machine::IsMotorOn(void) const
{
	return (m_MotorOn.dGet() != 0.);
}

asynchronous_machine::~asynchronous_machine(void)
{
	// destroy private data
#ifdef DEBUG
	std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << std::endl;
#endif
}

/**
 * @param OH OH.Loadable() returns a reference to the stream intended for the output of loadable elements.<BR>
 * \brief Writes private data to a file at each time step.
 *
 * The output format is as follows:<BR>
 * <TABLE>
 * <TR><TD>Column</TD>          <TD>Value</TD>                 <TD>Description</TD></TR>
 * <TR><TD>1</TD>               <TD>GetLabel()</TD>            <TD>The label of the element.</TD></TR>
 * <TR><TD>2</TD>               <TD>\f$M\f$</TD>               <TD>motor torque</TD></TR>
 * <TR><TD>3</TD>               <TD>\f$\dot{M}\f$</TD>         <TD>derivative of the motor torque versus time</TD></TR>
 * <TR><TD>4</TD>               <TD>\f$\ddot{M}\f$</TD>        <TD>second derivative of the motor torque versus time</TD></TR>
 * <TR><TD>5</TD>               <TD>\f$\dot{\varphi}\f$</TD>         <TD>the angular velocity of the rotor node with respect to the stator node</TD></TR>
 * <TR><TD>6</TD>               <TD>\f$\ddot{\varphi}\f$</TD>   <TD>the angular acceleration of the rotor node with respect to the stator node</TD></TR>
 * <TR><TD>7</TD>               <TD>\f$\Omega_S\f$</TD>        <TD>the synchronous angular velocity specified by the a drive</TD></TR>
 * <TR><TD>8</TD>               <TD>\f$\gamma\f$</TD>          <TD>a flag which specifies if the motor is powered on</TD></TR>
 * </TABLE>
 */

void
asynchronous_machine::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		// Note: The ODE for the asynchronous machine is defined for positive sign of OmegaS only!
		// At user level the signs of M, dM/dt and d^2M/dt^2 are defined to be positive
		// in positive coordinate system direction in the reference frame of the stator.
		doublereal OmegaS = m_OmegaS.dGet();
		doublereal M = m_M*copysign(1., OmegaS);
		doublereal dM_dt = m_dM_dt*copysign(1., OmegaS);
		doublereal dM_dt2 = m_dM_dt2*copysign(1., OmegaS);

		if (OH.UseText(OutputHandler::LOADABLE)) {
			// output the current state: the column layout is as follows

			// 1 2     3      4     5         6      7
			// M dM_dt dM_dt2 omega domega_dt OmegaS gamma
			//
			// 0	label of the element
			// 1	motor torque
			// 2	derivative of the motor torque versus time
			// 3	second derivative of the motor torque versus time
			// 4	angular velocity of the rotor node relative to the stator node around axis 3
			// 5	angular acceleration of the rotor node relative to the stator node around axis 3
			// 6	synchronous angular velocity (might be a function of the time in case of an frequency inverter)
			// 7	power supply is turned on or off

			OH.Loadable() << GetLabel()
				<< " " << M
				<< " " << dM_dt
				<< " " << dM_dt2
				<< " " << m_omega
				<< " " << m_domega_dt
				<< " " << OmegaS
				<< " " << IsMotorOn()
				<< std::endl;
		}
	}
}

/**
 * @param piNumRows pointer to a variable which receives the number of rows of the contribution to the residual vector and the jacobian matrix
 * @param piNumCols pointer to a variable which receives the number of columns of the contribution to the jacobian matrix
 * \brief The size of the contribution to the residual and the jacobian matrix is determined by this member function.
 */
void
asynchronous_machine::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 9;
	*piNumCols = 9;
}

/**
 * @return Returns the number of private degrees of freedom. In our case \f$M\f$, \f$\dot{M}\f$ and \f$\dot{\varphi}\f$ are private degrees of freedom only accessible to this element.
 * \brief The number of private degrees of freedom is determined by this member function.
 */
unsigned int
asynchronous_machine::iGetNumDof(void) const
{
	return 3;
}

/**
 * @param i zero based index of the private degree of freedom (zero for the first private degree of freedom, one for second ...)
 * @return Returns if the private degrees of freedom specified by i refer to differential or algebraic variables. In our case all private degrees of freedom are differential variables.
 * \brief The type of the equation of the private degrees of freedom is determined by this member function.
 */
DofOrder::Order
asynchronous_machine::GetDofType(unsigned int i) const
{

	return DofOrder::DIFFERENTIAL;
}

/**
 * @param out output stream that receives the formated output
 * @param prefix should be output in the first column
 *
 * \brief This member function is called if the statement "print: dof description;" is specified in the input file.
 *
 * It prints a short description of the private degrees of freedom \f$M\f$, \f$\dot{M}\f$ and \f$\dot{\varphi}\f$.
 */
std::ostream&
asynchronous_machine::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": motor torque derivative [MP]" << std::endl
		<< prefix << iIndex + 2 << ": motor torque [M]" << std::endl
		<< prefix << iIndex + 3 << ": relative angular velocity [omega]" << std::endl;

	return out;
}

/**
 * @param out output stream that receives the formated output
 * @param prefix should be output in the first column
 *
 * \brief This member function is called if the statement "print: equation description;" is specified in the input file.
 *
 * It prints a short description of the residual of the private degrees of freedom \f$f_1\f$, \f$f_2\f$ and \f$f_9\f$.
 */
std::ostream&
asynchronous_machine::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	integer iIndex = iGetFirstIndex();

	out
		<< prefix << iIndex + 1 << ": motor DAE [f1]" << std::endl
		<< prefix << iIndex + 2 << ": motor torque derivative [f2]" << std::endl
		<< prefix << iIndex + 3 << ": angular velocity derivative [f9]" << std::endl;

	return out;
}

/**
 * @return Returns the number of private data available for this element.
 * \brief This member function is called when a bind statement or a element drive is used to access private data of this element.
 *
 * In our case \f$M\f$, \f$\dot{M}\f$, \f$\ddot{M}\f$, \f$\dot{\varphi}\f$ and \f$\ddot{\varphi}\f$ are available.
 */
unsigned int
asynchronous_machine::iGetNumPrivData(void) const
{
	return 5;
}

/**
 * @param s specifies the name of the private data to be accessed by a bind statement or an element drive caller.
 * @return returns the one based index of the private data.
 * \brief The following string values can be specified in a bind statement or in an element drive caller.
 *
 * <TABLE>
 * <TR><TD>index</TD>   <TD>string</TD>  <TD>variable</TD></TR>
 * <TR><TD>1</TD>       <TD>M</TD>       <TD>\f$M\f$</TD></TR>
 * <TR><TD>2</TD>       <TD>MP</TD>      <TD>\f$\dot{M}\f$</TD></TR>
 * <TR><TD>3</TD>       <TD>MPP</TD>     <TD>\f$\ddot{M}\f$</TD></TR>
 * <TR><TD>4</TD>       <TD>omega</TD>   <TD>\f$\dot{\varphi}\f$</TD></TR>
 * <TR><TD>5</TD>       <TD>omegaP</TD>  <TD>\f$\ddot{\varphi}\f$</TD></TR>
 * </TABLE>
 */
unsigned int
asynchronous_machine::iGetPrivDataIdx(const char *s) const
{
	static const struct {
		int index;
		char name[7];
	}

	data[] = {
			{ 1, "M" },		// motor torque
			{ 2, "MP"},		// derivative of the motor torque versus time diff(M,t) (named MP instead of dM_dt for compatibility with other elements)
			{ 3, "MPP"},	// second derivative of the motor torque  versus time diff(M,t,2) (named MPP instead of dM_dt2 for compatibility with other elements)
			{ 4, "omega"},	// angular velocity of the rotor node relative to the stator node
			{ 5, "omegaP"}	// angular acceleration of the rotor node relative to the stator node (named omegaP instead of domega_dt for compatibility with other elements)
	};

	for (unsigned i = 0; i < sizeof(data) / sizeof(data[0]); ++i ) {
		if (0 == strcmp(data[i].name,s)) {
			return data[i].index;
		}
	}

	silent_cerr("asynchronous_machine(" << GetLabel() << "): no private data \"" << s << "\"" << std::endl);

	return 0;
}

/**
 * @param i the one based index of the private data
 * @return returns the value of the private data addressed by i.
 * \brief This member function is called whenever a bind statement or a element drive is used to access private data.
 */
doublereal
asynchronous_machine::dGetPrivData(unsigned int i) const
{
	// Note: The ODE for the asynchronous machine is defined for positive sign of OmegaS only!
	// At user level the signs of M, dM/dt and d^2M/dt^2 are defined to be positive
	// in positive coordinate system direction in the reference frame of the stator.

	const doublereal OmegaS = m_OmegaS.dGet();

	switch (i) {
	case 1:
		return m_M*copysign(1., OmegaS);
	case 2:
		return m_dM_dt*copysign(1., OmegaS);
	case 3:
		return m_dM_dt2*copysign(1., OmegaS);
	case 4:
		return m_omega;
	case 5:
		return m_domega_dt;
	}

	throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
}

/**
 * See iGetInitialNumDof().
 */
void
asynchronous_machine::SetInitialValue(VectorHandler& X)
{
	return;
}

/**
 * @param WorkMat sparse or full submatrix which receives the contribution to the jacobian matrix \f$Jac\f$
 * @param dCoef the derivative coefficient is defined as \f$\Delta\boldsymbol{y} = dCoef \, \Delta\dot{\boldsymbol{y}}\f$
 * @param XCurr the vector of the global degrees of freedom \f$\boldsymbol{y}\f$ at the current iteration step
 * @param XPrimeCurr the vector of the derivative of the global degrees of freedom \f$\dot{\boldsymbol{y}}\f$<BR>
 * \brief Computes the contribution to the jacobian matrix \f$\boldsymbol{Jac}\f$.<BR>
 *
 * \f[\boldsymbol{Jac}=-\frac{\partial\boldsymbol{f}}{\partial\dot{\boldsymbol{y}}}-dCoef\,\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{y}}\f]<BR>
 * For the definition of \f$\boldsymbol{f}\f$ see AssRes().<BR>
 * \f[\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{y}}=\left(\begin{array}{ccccc} \frac{\partial f_{1}}{\partial y_{1}} & \frac{\partial f_{1}}{\partial y_{2}} & \frac{\partial f_{1}}{\partial\boldsymbol{y}_{3}} & \frac{\partial f_{1}}{\partial\boldsymbol{y}_{4}} & \frac{\partial f_{1}}{\partial y_{5}}\\ \frac{\partial f_{2}}{\partial y_{1}} & \frac{\partial f_{2}}{\partial y_{2}} & \frac{\partial f_{2}}{\partial\boldsymbol{y}_{3}} & \frac{\partial f_{2}}{\partial\boldsymbol{y}_{4}} & \frac{\partial f_{2}}{\partial y_{5}}\\ \frac{\partial\boldsymbol{f}_{3}}{\partial y_{1}} & \frac{\partial\boldsymbol{f}_{3}}{\partial y_{2}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{y}_{3}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{y}_{4}} & \frac{\partial\boldsymbol{f}_{3}}{\partial y_{5}}\\ \frac{\partial\boldsymbol{f}_{4}}{\partial y_{1}} & \frac{\partial\boldsymbol{f}_{4}}{\partial y_{2}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{y}_{3}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{y}_{4}} & \frac{\partial\boldsymbol{f}_{4}}{\partial y_{5}}\\ \frac{\partial f_{5}}{\partial y_{1}} & \frac{\partial f_{5}}{\partial y_{2}} & \frac{\partial f_{5}}{\partial\boldsymbol{y}_{3}} & \frac{\partial f_{5}}{\partial\boldsymbol{y}_{4}} & \frac{\partial f_{5}}{\partial y_{5}}\end{array}\right)\f]<BR>
 * \f[\frac{\partial\boldsymbol{f}}{\partial\boldsymbol{\dot{y}}}=\left(\begin{array}{ccccc} \frac{\partial f_{1}}{\partial\dot{y}_{1}} & \frac{\partial f_{1}}{\partial\dot{y}_{2}} & \frac{\partial f_{1}}{\partial\boldsymbol{\dot{y}}_{3}} & \frac{\partial f_{1}}{\partial\boldsymbol{\dot{y}}_{4}} & \frac{\partial f_{1}}{\partial\dot{y}_{5}}\\ \frac{\partial f_{2}}{\partial\dot{y}_{1}} & \frac{\partial f_{2}}{\partial\dot{y}_{2}} & \frac{\partial f_{2}}{\partial\boldsymbol{\dot{y}}_{3}} & \frac{\partial f_{2}}{\partial\boldsymbol{\dot{y}}_{4}} & \frac{\partial f_{2}}{\partial\dot{y}_{5}}\\ \frac{\partial\boldsymbol{f}_{3}}{\partial\dot{y}_{1}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\dot{y}_{2}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{\dot{y}}_{3}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{\dot{y}}_{4}} & \frac{\partial\boldsymbol{f}_{3}}{\partial\dot{y}_{5}}\\ \frac{\partial\boldsymbol{f}_{4}}{\partial\dot{y}_{1}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\dot{y}_{2}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{\dot{y}}_{3}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{\dot{y}}_{4}} & \frac{\partial\boldsymbol{f}_{4}}{\partial\dot{y}_{5}}\\ \frac{\partial f_{5}}{\partial\dot{y}_{1}} & \frac{\partial f_{5}}{\partial\dot{y}_{2}} & \frac{\partial f_{5}}{\partial\dot{\boldsymbol{y}}_{3}} & \frac{\partial f_{5}}{\partial\boldsymbol{\dot{y}}_{4}} & \frac{\partial f_{5}}{\partial\dot{y}_{5}}\end{array}\right)\f]<BR>
 * \f[\frac{\partial f_{1}}{\partial y_{1}}=\left(2\, s_{K}\,+\frac{sign\left(\Omega_{S}\right)\,\dot{y}_{5}}{s\,\Omega_{S}^{2}}\right)\,\left|\Omega_{S}\right|\f]
 * \f[\frac{\partial f_{1}}{\partial\dot{y}_{1}}=1\f]
 * \f[\frac{\partial f_{1}}{\partial y_{2}}=\left(s_{K}^{2}+s^{2}\right)\,\Omega_{S}^{2}+\frac{sign\left(\Omega_{S}\right)\,\dot{y}_{5}\, sK}{s}\f]
 * \f[\frac{\partial f_{1}}{\partial y_{5}}=\frac{\partial f_{1}}{\partial s}\,\frac{\partial s}{\partial y_{5}}\f]
 * \f[\frac{\partial f_{1}}{\partial s}=-\frac{y_{1}\,\dot{y}_{5}}{s^{2}\,\Omega_{S}}+\left(2\, s\,\Omega_{S}^{2}-\frac{sign\left(\Omega_{S}\right)\:\dot{y}_{5}\, s_{K}}{s^{2}}\right)\, y_{2}-2\, M_{K}\, s_{K}\,\Omega_{S}^{2}\f]
 * \f[\frac{\partial s}{\partial y_{5}}=-\frac{1}{\Omega_{S}}\f]
 * \f[\frac{\partial f_{1}}{\partial\dot{y}_{5}}=\frac{y_{1}}{s\,\Omega_{S}}+sign\left(\Omega_{S}\right)\,\frac{s_{K}}{s}\, y_{2}\f]
 * \f[\frac{\partial f_{2}}{\partial y_{1}}=-1\f]
 * \f[\frac{\partial f_{2}}{\partial\dot{y}_{2}}=1\f]
 * \f[\frac{\partial f_{3}}{\partial y_{2}}=\boldsymbol{R}_{2}\,\boldsymbol{e}_{3}\, sign\left(\Omega_{S}\right)\f]
 * \f[\frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{y}_{4}}=-\left\langle \boldsymbol{R}_{2}^{\left(0\right)}\,\boldsymbol{e}_{3}\, y_{2}\, sign\left(\Omega_{S}\right)\right\rangle \f]
 * \f[\frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{y}_{2}}=-\frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{y}_{2}}\f]
 * \f[\frac{\partial\boldsymbol{f}_{4}}{\partial\boldsymbol{y}_{4}}=-\frac{\partial\boldsymbol{f}_{3}}{\partial\boldsymbol{y}_{4}}\f]
 * \f[\frac{\partial f_{5}}{\partial\boldsymbol{y}_{3}}=\boldsymbol{e}_{3}^{T}\,\boldsymbol{R}_{2}^{T}\,\left\langle \boldsymbol{\omega}_{ref_{1}}\right\rangle =-\left(\left\langle \boldsymbol{\omega}_{ref_{1}}\right\rangle \,\boldsymbol{R}_{2}\,\boldsymbol{e}_{3}\right)^{T}\f]
 * \f[\frac{\partial f_{5}}{\partial\boldsymbol{y}_{4}}=-\boldsymbol{e}_{3}^{T}\,\boldsymbol{R}_{2}^{\left(0\right)^{T}}\,\left\langle \boldsymbol{\omega}_{1}-\boldsymbol{\omega}_{2}\right\rangle -\boldsymbol{e}_{3}^{T}\,\boldsymbol{R}_{2}^{T}\,\left\langle \boldsymbol{\omega}_{ref_{2}}\right\rangle =\left[\left\langle \boldsymbol{\omega}_{1}-\boldsymbol{\omega}_{2}\right\rangle \,\boldsymbol{R}_{2}^{\left(0\right)}\,\boldsymbol{e}_{3}+\left\langle \boldsymbol{\omega}_{ref_{2}}\right\rangle \,\boldsymbol{R}_{2}\,\boldsymbol{e}_{3}\right]^{T}\f]
 * In order to compute the derivatives versus the global degrees of freedom of structural nodes the following rules can be applied:
 * \f[\frac{\partial \left(\boldsymbol{R}_1\,\boldsymbol{v}_1\right)}{\partial \boldsymbol{g}_1}\approx-\left\langle \boldsymbol{R}^{\left(0\right)}\,\boldsymbol{v} \right\rangle \f]<BR>
 * \f[\frac{\partial \left(\boldsymbol{R}_1^{T}\,\boldsymbol{v}\right)}{\partial \boldsymbol{g}_1}\approx\boldsymbol{R}_1^{\left(0\right)T}\,\left\langle\boldsymbol{v}\right\rangle\f]
 * \f[\frac{\partial\boldsymbol{\omega}_1}{\partial\dot{\boldsymbol{g}}_1}\approx\boldsymbol{I}\f]
 * \f[\frac{\partial\boldsymbol{\omega}_1}{\partial\boldsymbol{g}_1}\approx-\left\langle\boldsymbol{\omega}_{1_{ref}}\right\rangle\f]
 * \f$\boldsymbol{R}_1\f$ is the corrected rotation matrix of node 1 at the current iteration which can be obtained by GetRCurr().<BR>
 * \f$\boldsymbol{R}^{(0)}_1\f$ is the predicted rotation matrix of node 1 at the current time step which can be obtained by GetRRef().<BR>
 * \f$\boldsymbol{g}_1\f$ is the vector of Gibbs Rodriguez rotation parameters of node 1.<BR>
 * The rotation matrix \f$\boldsymbol{R}_1\f$ is internally computed as follows:<BR>
 * \f[\boldsymbol{R}_1=\boldsymbol{R}_{1_{\Delta}}\,\boldsymbol{R}^{(0)}_1\f]
 * \f$\boldsymbol{R}_{1_{\Delta}}\f$ is the incremental rotation matrix which can be computed by means of the Gibbs Rodriguez rotation parameters \f$\boldsymbol{g}_1\f$.<BR>
 * \f[\boldsymbol{R}_{1_{\Delta}}=\boldsymbol{I}+\frac{4}{4+\boldsymbol{g}_1^T\,\boldsymbol{g}_1}\,\left(\left\langle\boldsymbol{g}\right\rangle+\frac{1}{2}\,\left\langle\boldsymbol{g}_1\right\rangle\,\left\langle\boldsymbol{g}_1\right\rangle\right)\f]
 * \f$\boldsymbol{\omega}_1\f$ is the corrected angular velocity of node 1 at the current iteration which can be obtained by GetWCurr().<BR>
 * \f[\left\langle\boldsymbol{\omega}_1\right\rangle=\dot{\boldsymbol{R}}_1\,\boldsymbol{R}_1^T\f]
 * \f[\boldsymbol{\omega}_1\approx\dot{\boldsymbol{g}}_1+\boldsymbol{R}_{1_{\Delta}}\,\boldsymbol{\omega}_{1_{ref}}\f]
 * \f$\boldsymbol{\omega}_{1_{ref}}\f$ is the predicted angular velocity of node 1 at the current time step which can be obtained by GetWRef().<BR>
 * \f[\left\langle\boldsymbol{\omega}_{1_{ref}}\right\rangle=\dot{\boldsymbol{R}}_1^{(0)}\,\boldsymbol{R}_1^{(0)T}\f]
 * \f$\boldsymbol{v}\f$ is a const vector in the reference frame \f$\boldsymbol{R}_1\f$.<BR>
 * \f[\boldsymbol{v}=\left(\begin{array}{c} v_{1}\\ v_{2}\\ v_{3}\end{array}\right)\f]
 * \f$\left\langle \boldsymbol{v}\right\rangle\f$ is a matrix that rearranges the components of \f$\boldsymbol{v}\f$ according to the cross product rule:<BR>
 * \f$\left\langle\boldsymbol{v}\right\rangle\,\boldsymbol{a}=-\left\langle\boldsymbol{a}\right\rangle\,\boldsymbol{v}=\boldsymbol{v}\times\boldsymbol{a}=-\boldsymbol{a}\times\boldsymbol{v}\f$<BR>
 * \f[\left\langle \boldsymbol{v} \right\rangle=\left(\begin{array}{ccc} 0 & -v_{3} & v_{2}\\ v_{3} & 0 & -v_{1}\\ -v_{2} & v_{1} & 0\end{array}\right)\f]
 */
VariableSubMatrixHandler&
asynchronous_machine::AssJac(VariableSubMatrixHandler& WorkMat_,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
#ifdef DEBUG
	std::cerr << __FILE__ << ':' << __LINE__ << ':' << __PRETTY_FUNCTION__ << std::endl;
#endif

	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMat_.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	const integer iFirstIndex = iGetFirstIndex();

	const integer intTorqueDerivativeColumnIndex = iFirstIndex + 1;
	const integer intTorqueColumnIndex 			 = iFirstIndex + 2;
	const integer intOmegaColumnIndex 			 = iFirstIndex + 3;

	const integer intTorqueDerivativeRowIndex 	 = iFirstIndex + 1;
	const integer intTorqueRowIndex 			 = iFirstIndex + 2;
	const integer intOmegaRowIndex 				 = iFirstIndex + 3;

	const integer intRotorPositionIndex  = m_pRotorNode->iGetFirstPositionIndex();
	const integer intStatorPositionIndex = m_pStatorNode->iGetFirstPositionIndex();

	const integer intRotorMomentumIndex  = m_pRotorNode->iGetFirstMomentumIndex();
	const integer intStatorMomentumIndex = m_pStatorNode->iGetFirstMomentumIndex();

	/*
	 * M 		 ... motor torque
	 * diff(M,t) ... derivative of the motor torque versus time
	 * g1		 ... rotation parameter of the rotor
	 * g2		 ... rotation parameter of the stator
	 * R2 		 ... rotation matrix of the stator node
	 * omega1    ... rotor angular velocity
	 * omega2 	 ... stator angular velocity
	 *
	 * e3 = ( 0, 0, 1 )^T ... axis of rotation in the reference frame of the stator node
	 *
	 * omega = e3^T * R2^T * ( omega1 - omega2 )
	 *
	 * s = 1 - omega / OmegaS
	 *
	 * 		| y1 |	 | diff(M,t) |
	 * 		| y2 |	 | M 		 |
	 * 	y =	| y3 | = | g1		 |
	 * 		| y4 |	 | g2 		 |
	 * 		| y5 |   | omega     |
	 *
	 * 		| f1 |   |  diff(y1,t) + ( 2 * sK + sign(OmegaS) * diff(y5,t) / ( s * OmegaS^2 ) ) * y1 * abs(OmegaS) + ( ( sK^2 + s^2 ) * OmegaS^2 + sign(OmegaS) * diff(y5,t) * sK / s ) * y2 - 2 * MK * sK * s * OmegaS^2 |
	 * 		| f2 |   |  diff(y2,t) - y1																																				    		   	   					 |																																				    		   |
	 * 	f =	| f3 | = |  R2 * e3 * y2																																				    		   	   					 |
	 * 		| f4 |   | -R2 * e3 * y2																																				    		       					 |
	 * 		| f5 |   |  y5 - e3^T * R2^T * ( omega1 - omega2 )																																	       					 |
	 *
	 * 					      1,       2,       3,       6,    	  9
	 * 		          | df1/dy1, df1/dy2, df1/dy3, df1/dy4, df1/dy5 |  1
	 * 				  | df2/dy1, df2/dy2, df2/dy3, df2/dy4, df2/dy5 |  2
	 * 	diff(f,y) =   | df3/dy1, df3/dy2, df3/dy3, df3/dy4, df3/dy5 |  3
	 * 		          | df4/dy1, df4/dy2, df4/dy3, df4/dy4, df4/dy5 |  6
	 * 		          | df5/dy1, df5/dy2, df5/dy3, df5/dy4, df5/dy5 |  9
	 *
	 * WorkMat = -diff(f,diff(y,t)) - dCoef * diff(f,y)
	 */
	WorkMat.PutRowIndex(1, intTorqueDerivativeRowIndex);
	WorkMat.PutColIndex(1, intTorqueDerivativeColumnIndex);

	WorkMat.PutRowIndex(2, intTorqueRowIndex);
	WorkMat.PutColIndex(2, intTorqueColumnIndex);

	for (int iCnt = 1; iCnt <= 3; ++iCnt) {
		WorkMat.PutRowIndex(iCnt + 2, intRotorMomentumIndex + iCnt + 3);
		WorkMat.PutColIndex(iCnt + 2, intRotorPositionIndex + iCnt + 3);

		WorkMat.PutRowIndex(iCnt + 5, intStatorMomentumIndex + iCnt + 3);
		WorkMat.PutColIndex(iCnt + 5, intStatorPositionIndex + iCnt + 3);
	}

	WorkMat.PutRowIndex(9, intOmegaRowIndex);
	WorkMat.PutColIndex(9, intOmegaColumnIndex);

	Vec3 e3 = m_pStatorNode->GetRCurr().GetCol(3);	// corrected axis of rotation of the stator node
	Vec3 e3_0 = m_pStatorNode->GetRRef().GetCol(3);	// predicted axis of rotation of the stator node
	const Vec3& omega1 = m_pRotorNode->GetWCurr();	// corrected angular velocity of the rotor node
	const Vec3& omega2 = m_pStatorNode->GetWCurr(); // corrected angular velocity of the rotor node
	const Vec3& omega1_ref = m_pRotorNode->GetWRef(); // predicted angular velocity of the rotor node
	const Vec3& omega2_ref = m_pStatorNode->GetWRef(); // predicted angular velocity of the stator node

	const doublereal OmegaS = m_OmegaS.dGet(); // synchronous angular velocity

	const doublereal y1 		= XCurr(intTorqueDerivativeRowIndex);
	const doublereal y2 		= XCurr(intTorqueRowIndex);
	const doublereal y5 		= XCurr(intOmegaRowIndex);
#if 0 // unused
	const doublereal y1_dot 	= XPrimeCurr(intTorqueDerivativeRowIndex);
	const doublereal y2_dot 	= XPrimeCurr(intTorqueRowIndex);
#endif
	const doublereal y5_dot 	= XPrimeCurr(intOmegaRowIndex);

	doublereal s = 1 - y5 / OmegaS; // slip of the rotor

    if ( std::abs(s) < sm_SingTol )
    {
        silent_cerr("\nasynchronous_machine(" << GetLabel() << "): warning slip rate s = " << s << " is close to machine precision!\n");
        //FIXME: avoid division by zero
        //FIXME: results might be inaccurate
        s = copysign(sm_SingTol, s);
    }

	doublereal df1_dy1, df1_dy2, df1_dy5;
	doublereal df1_dy1_dot, df1_dy2_dot, df1_dy5_dot;
	doublereal df2_dy1, df2_dy2_dot;

	if (IsMotorOn()) {
		// power supply is turned on

		// df1_dy1 = diff(f1,y1)
		df1_dy1 = ( 2 * m_sK + copysign(1., OmegaS) * y5_dot / ( s * std::pow(OmegaS,2) ) ) * abs(OmegaS);
		// df1_dy2 = diff(f1,y2)
		df1_dy2 = ( std::pow(m_sK,2) + std::pow(s,2) ) * std::pow(OmegaS,2) + copysign(1., OmegaS) * y5_dot * m_sK / s;
		// df1_dy1_dot = diff(f1,diff(y1,t))
		df1_dy1_dot = 1.;
		// df2_dy1 = diff(f2,dy1)
		df2_dy1 = -1.;
		// df2_dy2_dot = diff(f2,diff(y2,t))
		df2_dy2_dot = 1.;
		// df1_ds = diff(f1,s)
		const doublereal df1_ds = -y1 * y5_dot / ( std::pow(s,2) * OmegaS ) + ( 2 * s * std::pow(OmegaS,2) - copysign(1., OmegaS) * y5_dot * m_sK / std::pow(s,2) ) * y2 - 2 * m_MK * m_sK * std::pow(OmegaS,2);
		// df1_dy5 = diff(f1,y5) = diff(f1,s) * diff(s,y5)
		df1_dy5 = df1_ds * ( -1. / OmegaS );
		// df1_dy5_dot = diff(f1,diff(y5,t))
		df1_dy5_dot = y1 / ( s * OmegaS ) + copysign(1., OmegaS) * y2 * m_sK / s;

	} else {
		// power supply is turned off
		df1_dy1 = 0.;
		df1_dy2 = 1.;
		df1_dy5 = 0.;
		df1_dy1_dot = 0.;
		df1_dy2_dot = 0.;
		df1_dy5_dot = 0.;

		df2_dy1 = 1.;
		df2_dy2_dot = 0;
	}

	// df3_dy2 = diff(f3,y2)
	const Vec3 df3_dy2 = e3 * copysign(1., OmegaS); // df3_dy2 = R2 * e3 * sign(OmegaS)
	// df4_dy2 = diff(f4,y2)
	const Vec3 df4_dy2 = -df3_dy2; // df4_dy2 = -R2 * e3
	// df3_dy4 = diff(f3,y4)
	const Mat3x3 df3_dy4 = -Mat3x3(MatCross, e3_0 * (y2 * copysign(1., OmegaS)) ); // df3_dy4 = -skew( R2^(0) * e3 * y2 * sign(OmegaS) )
	// df4_dy4 = diff(f4,y4)
	const Mat3x3 df4_dy4 = -df3_dy4;
	// df5_dy5 = diff(f5,y5)
	const doublereal df5_dy5 = 1.;
	// df5_dy3_dot_T = transpose(diff(f5,diff(y3,t)))
	const Vec3 df5_dy3_dot_T = -e3; // diff(y5,diff(y3,t)) = -e3^T * R2^T
	// df5_dy4_dot_T = transpose(diff(f5,diff(y4,t)))
	const Vec3 df5_dy4_dot_T = -df5_dy3_dot_T;  // diff(y5,diff(y4,t)) = e3^T * R2^T
	// df5_dy3_T = transpose(diff(y5,y3))
	const Vec3 df5_dy3_T = -omega1_ref.Cross(e3);
	// df5_dy4_T = transpose(diff(y5,y4))
	const Vec3 df5_dy4_T = ( omega1 - omega2 ).Cross( e3_0 ) + omega2_ref.Cross( e3 );

	WorkMat.PutCoef( 1, 1,  -df1_dy1_dot   - dCoef * df1_dy1 );
	WorkMat.PutCoef( 1, 2, 				   - dCoef * df1_dy2 );
	WorkMat.PutCoef( 1, 9,  -df1_dy5_dot   - dCoef * df1_dy5 );
	WorkMat.PutCoef( 2, 1, 				   - dCoef * df2_dy1 );
	WorkMat.PutCoef( 2, 2,  -df2_dy2_dot 						    );
	WorkMat.Put(     3, 2,				   	 df3_dy2   * ( -dCoef ) );
	WorkMat.Put(     3, 6,				   	 df3_dy4   * ( -dCoef ) );
	WorkMat.Put(     6, 2,                 	 df4_dy2   * ( -dCoef ) );
	WorkMat.Put(     6, 6,                 	 df4_dy4   * ( -dCoef ) );
	WorkMat.PutT(    9, 3,	-df5_dy3_dot_T + df5_dy3_T * ( -dCoef ) );
	WorkMat.PutT(    9, 6,  -df5_dy4_dot_T + df5_dy4_T * ( -dCoef ) );
	WorkMat.PutCoef( 9, 9,				     df5_dy5   * ( -dCoef ) );

#ifdef DEBUG
	std::cerr << __FILE__ ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << std::endl;
	std::cerr << "s = " << s << std::endl;
	std::cerr << "WorkMat=" << std::endl
		 << WorkMat << std::endl << std::endl;
#endif

	return WorkMat_;
}

/**
 * @param WorkVec subvector which receives the contribution to the residual \f$\boldsymbol{f}\f$
 * @param dCoef the derivative coefficient is defined as \f$\Delta\boldsymbol{y} = dCoef \, \Delta\dot{\boldsymbol{y}}\f$
 * @param XCurr the vector of the global degrees of freedom \f$\boldsymbol{y}\f$ at the current iteration step
 * @param XPrimeCurr the vector of the derivative of the global degrees of freedom \f$\dot{\boldsymbol{y}}\f$<BR>
 * \brief This member function implements the equation of an asynchronous machine according to
 * <A HREF="http://books.google.at/books?id=z6B3lcxDgkkC&printsec=frontcover&dq=Maschinendynamik++%26+dresig&hl=de&ei=SwvETYi8NJSw4Ab_luX3BA&sa=X&oi=book_result&ct=result&resnum=1&ved=0CDIQ6AEwAA#v=onepage&q&f=false">
 * Maschinendynamik Springer 2009, Dresig, Holzweißig.
 * </A><BR>
 *
 * \f[\ddot{M}+\left(2\, s_{K}+\frac{sign\left(\Omega_{S}\right)\,\ddot{\varphi}}{s\,\Omega_{S}^{2}}\right)\,\dot{M\,}\left|\Omega_{s}\right|+\left[\left(s_{K}^{2}+s^{2}\right)\,\Omega_{S}^{2}+\frac{sign\left(\Omega_{S}\right)\,\ddot{\varphi}\, s_{K}}{s}\right]\, M=2\, M_{K}\, s_{K}\, s\,\Omega_{S}^{2}\f]<BR>
 * The term \f$sign\left(\Omega_{S}\right)\f$ and the absolute operator in \f$\left|\Omega_{s}\right|\f$ are modifications of the original formula since the formula in the literature does not handle negative synchronous angular velocities \f$\Omega_{S}\f$.
 * If the synchronous angular velocity \f$\Omega_{S}\f$ has a negative value, the sense of rotation is assumed to be negative.<BR>
 * \f$M\f$ is motor torque<BR>
 * \f$M_{K}\f$ is the breakdown torque<BR>
 * \f$s_{K}\f$ is the breakdown slip<BR>
 * \f$\Omega_{S}\f$ is the synchronous angular velocity of the machine<BR>
 * \f$\Omega_{S}=2\,\pi\frac{f}{p}\f$<BR>
 * \f$f\f$ is the power frequency<BR>
 * \f$p\f$ is the number of terminal pairs<BR>
 * \f$s\f$ is the actual slip<BR>
 * \f$s=1-\frac{\dot{\varphi}}{\Omega_{S}}\f$<BR>
 * \f$\dot{\varphi}\f$ is the relative angular velocity of the rotor with respect to the stator<BR>
 * The axis of rotation is assumed to be axis 3 in the reference frame of the stator node.<BR>
 * \f$\dot{\varphi}=\boldsymbol{e}_{3}^{T}\cdot\boldsymbol{R}_{2}^{T}\cdot\left(\boldsymbol{\omega}_{1}-\boldsymbol{\omega}_{2}\right)\f$<BR>
 * \f$\boldsymbol{\omega}_{1}\f$ is the angular velocity of the rotor node in the global reference frame which can be obtained by GetWCurr()<BR>
 * \f$\boldsymbol{\omega}_{2}\f$ is the angular velocity of the stator node in the global reference frame<BR>
 * \f$\boldsymbol{R}_{2}\f$ is the rotation matrix of the stator node which can be obtained by GetRCurr()<BR>
 * The subvector of the global degrees of freedom XCurr used by this element is<BR>
 * \f[\boldsymbol{y}=\left(\begin{array}{c} y_{1}\\ y_{2}\\ \boldsymbol{y}_{3}\\ \boldsymbol{y}_{4}\\ y_{5}\end{array}\right)=\left(\begin{array}{c}\dot{M}\\M\\ \boldsymbol{g}_{1}\\ \boldsymbol{g}_{2}\\ \dot{\varphi}\end{array}\right)\f]<BR>
 * \f$\boldsymbol{g}_{1}\f$ is the vector of rotation parameters of the rotor node.<BR>
 * \f$\boldsymbol{g}_{2}\f$ is the vector of the rotation parameters of the stator node.<BR>
 * The index of \f$\boldsymbol{g}_1\f$ and \f$\boldsymbol{g}_2\f$ can be obtained by iGetFirstMomentumIndex() + 4 ... 6.<BR>
 * The index of \f$y_1\f$, \f$y_2\f$ and \f$y_5\f$ in the vector of the global degrees of freedom XCurr can be obtained by iGetFirstIndex() + 1 ... 3.<BR>
 * The subvector of the derivatives of the global degrees of freedom used by this element is<BR>
 * \f[\dot{\boldsymbol{y}}=\left(\begin{array}{c} \dot{y_{1}}\\ \dot{y_{2}}\\ \dot{\boldsymbol{y}_{3}}\\ \dot{\boldsymbol{y}_{4}}\\ \dot{y_{5}}\end{array}\right)=\left(\begin{array}{c} \ddot{M}\\ \dot{M}\\ \dot{\boldsymbol{g}_{1}}\\ \dot{\boldsymbol{g}_{2}}\\ \ddot{\varphi}\end{array}\right)\f]<BR>
 * The additional degree of freedom \f$y_{5}=\dot{\varphi}\f$ is needed since MBDyn does not provide the angular acceleration of a node during the nonlinear solution phase.
 * The value returned by GetWPCurr() is updated only after convergence and can not be used for this reason.
 * The contribution to the residual of this element is<BR>
 * \f[\boldsymbol{f}=\left(\begin{array}{c} f_{1}\\ f_{2}\\ \boldsymbol{f}_{3}\\ \boldsymbol{f}_{4}\\ f_{5}\end{array}\right)\f]<BR>
 * \f[f_{1}=\dot{y}_{1}+\left(2\, s_{K}+\frac{sign\left(\Omega_{S}\right)\,\dot{y}_{5}}{s\,\Omega_{S}^{2}}\right)\, y_{1}\,\left|\Omega_{S}\right|+\left[\left(s_{K}^{2}+s^{2}\right)\,\Omega_{S}^{2}+\frac{sign\left(\Omega_{S}\right)\,\dot{y}_{5}\, s_{K}}{s}\right]\, y_{2}-2\, M_{K}\, s_{K}\, s\,\Omega_{S}^{2}\f]<BR>
 * \f[f_{2}=\dot{y}_{2}-y_{1}\f]<BR>
 * \f$\boldsymbol{f}_3\f$ is the torque imposed to the rotor node.<BR>
 * \f[\boldsymbol{f}_{3}=\boldsymbol{R}_{2}\,\boldsymbol{e}_{3}\, y_{2}\, sign\left(\Omega_{S}\right)\f]<BR>
 * \f$\boldsymbol{f}_{4}\f$ is the torque imposed to the stator node.<BR>
 * \f[\boldsymbol{f}_{4}=-\boldsymbol{f}_{3}\f]<BR>
 * \f[f_{5}=y_{5}-\boldsymbol{e}_{3}^{T}\,\boldsymbol{R}_{2}^{T}\,\left(\boldsymbol{\omega}_{1}-\boldsymbol{\omega}_{2}\right)\f]<BR>
 */

SubVectorHandler&
asynchronous_machine::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
#ifdef DEBUG
	std::cerr << __FILE__ << ':' << __LINE__ << ':' << __PRETTY_FUNCTION__ << std::endl;
#endif
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	const integer iFirstIndex = iGetFirstIndex();

	const integer intTorqueDerivativeRowIndex 	 = iFirstIndex + 1;
	const integer intTorqueRowIndex 			 = iFirstIndex + 2;
	const integer intOmegaRowIndex 				 = iFirstIndex + 3;

	const integer intRotorMomentumIndex  = m_pRotorNode->iGetFirstMomentumIndex();
	const integer intStatorMomentumIndex = m_pStatorNode->iGetFirstMomentumIndex();

	WorkVec.PutRowIndex(1, intTorqueDerivativeRowIndex);
	WorkVec.PutRowIndex(2, intTorqueRowIndex);

	for (int iCnt = 1; iCnt <= 3; ++iCnt) {
		WorkVec.PutRowIndex(iCnt + 2, intRotorMomentumIndex + iCnt + 3);
		WorkVec.PutRowIndex(iCnt + 5, intStatorMomentumIndex + iCnt + 3);
	}

	WorkVec.PutRowIndex(9, intOmegaRowIndex);

	Vec3 e3 = m_pStatorNode->GetRCurr().GetCol(3);
	const Vec3& omega1 = m_pRotorNode->GetWCurr();
	const Vec3& omega2 = m_pStatorNode->GetWCurr();

	const doublereal OmegaS = m_OmegaS.dGet();

	const doublereal y1 		= XCurr(intTorqueDerivativeRowIndex); 				// y1 = diff(M,t)
	const doublereal y2 		= XCurr(intTorqueRowIndex);			  				// y2 = M
	const doublereal y5 		= XCurr(intOmegaRowIndex);							// y5 = omega
	const doublereal y1_dot 	= XPrimeCurr(intTorqueDerivativeRowIndex); 			// diff(y1,t) = diff(M,t,2)
	const doublereal y2_dot 	= XPrimeCurr(intTorqueRowIndex);					// diff(y2,t) = diff(M,t)
	const doublereal y5_dot 	= XPrimeCurr(intOmegaRowIndex); 					// diff(y5,t) = diff(omega,t)

	doublereal s = 1 - y5 / OmegaS;

    if ( std::abs(s) < sm_SingTol )
    {
        silent_cerr("\nasynchronous_machine(" << GetLabel() << "): warning slip rate s = " << s << " is close to machine precision!\n");
        //FIXME: avoid division by zero
        //FIXME: results might be inaccurate
        s = copysign(sm_SingTol, s);
    }

	doublereal f1, f2;

	if (IsMotorOn()) {
		// power supply is switched on
		f1 = y1_dot + ( 2 * m_sK + copysign(1., OmegaS) * y5_dot / ( s * std::pow(OmegaS,2) ) ) * y1 * abs(OmegaS)
			+ ( ( std::pow(m_sK,2) + std::pow(s,2) ) * std::pow(OmegaS,2) + copysign(1., OmegaS) * y5_dot * m_sK / s ) * y2 - 2 * m_MK * m_sK * s * std::pow(OmegaS,2);
		f2 = y2_dot - y1;

	} else {
		// power supply is switched off: torque must be zero
		f1 = y2;
		f2 = y1;
	}

	const Vec3 f3 = e3 * (y2 * copysign(1., OmegaS));
	// const Vec3 f4 = -f3;

	const doublereal f5 = y5 - e3.Dot( omega1 - omega2 );

	WorkVec.PutCoef( 1, f1 );
	WorkVec.PutCoef( 2, f2 );
	WorkVec.Put(     3,	f3 );
	WorkVec.Sub(     6, f3 );
	WorkVec.PutCoef( 9, f5 );

#ifdef DEBUG
	std::cerr
		<< __FILE__ ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << std::endl
		<< "s = " << s << std::endl
		<< "y1 = " << y1 << std::endl
		<< "y1_dot = " << y1_dot << std::endl
		<< "y2 = " << y2 << std::endl
		<< "y2_dot = " << y2_dot << std::endl
		<< "f1 = " << f1 << std::endl
		<< "f2 = " << f2 << std::endl
		<< "f3 = " << f3 << std::endl
		<< "f4 = " << -f3 << std::endl
		<< "f5 = " << f5 << std::endl
		<< "WorkVec=" << std::endl
		<< WorkVec << std::endl << std::endl;
#endif

	return WorkVec;
}
/**
 * @return returns the number of connected nodes.
 * \brief This member function is called if the statement "print: node connection;" is specified in the input file.
 * */
int
asynchronous_machine::iGetNumConnectedNodes(void) const
{
	return 2;
}

/**
 * @param connectedNodes vector which receives pointers to the connected nodes
 * \brief This member function is called if the statement "print: node connection;" is specified in the input file.
 */
void
asynchronous_machine::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(iGetNumConnectedNodes());
	connectedNodes[0] = m_pRotorNode;
	connectedNodes[1] = m_pStatorNode;
}

/**
 * @param X vector of global degrees of freedom \f$\boldsymbol{y}\f$ that receives the initial values provided by this element.
 * @param XP vector of the derivative of the global degrees of freedom \f$\dot{\boldsymbol{y}}\f$
 * \brief This member function is called before the integration starts in order to set the initial values for the private degrees of freedom.
 *
 * In our case the initial values for \f$M\f$, \f$\dot{M}\f$ and \f$\dot{\varphi}\f$ are set in this routine to the values specified in the input file.
 */
void
asynchronous_machine::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	const integer intFirstIndex = iGetFirstIndex();

	//                             1        2        3
	X.Put( intFirstIndex + 1, Vec3(m_dM_dt, m_M,     m_omega) );
	XP.Put(intFirstIndex + 1, Vec3(     0., m_dM_dt,      0.) );
}

/**
 * @param out stream of the restart file
 * \brief This member function is called if "make restart file;" is specified in the input file.
 *
 * It should reproduce the syntax in the input file used to create this element.
 */
std::ostream&
asynchronous_machine::Restart(std::ostream& out) const
{
	// FIXME: mbdyn crashes when "make restart file;" is specified in the input file before this member function is invoked!
	out << "asynchronous_machine, "
		"rotor, " << m_pRotorNode->GetLabel() << ", "
		"stator, " << m_pStatorNode->GetLabel() << ", "
		"MK, " << m_MK << ", "
		"sK, " << m_sK << ", "
		"OmegaS, " << m_OmegaS.pGetDriveCaller()->Restart(out) << ", "
		"motor on, " << m_MotorOn.pGetDriveCaller()->Restart(out) << ", "
		"M0, " << m_M * copysign(1., m_OmegaS.dGet()) << ", "
		"MP0, " << m_dM_dt * copysign(1., m_OmegaS.dGet())  << ";" << std::endl;
	
	return out;
}

/**
 * \brief This member function must be implemented only if the initial assembly feature is requested.
 *
 * The initial assembly phase is needed only if initial values are provided which are not consistent with the algebraic constraints.
 * Since this element does not use algebraic constraints we do not have to implement
 * iGetInitialNumDof(), InitialWorkSpaceDim(), InitialAssJac(), InitialAssRes() and SetInitialValue().
 */
unsigned int
asynchronous_machine::iGetInitialNumDof(void) const
{
	return 0;
}

/**
 * See iGetInitialNumDof().
 */
void
asynchronous_machine::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

/**
 * See iGetInitialNumDof().
 */
VariableSubMatrixHandler&
asynchronous_machine::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	WorkMat.SetNullMatrix();

	return WorkMat;
}

/**
 * See iGetInitialNumDof().
 */
SubVectorHandler&
asynchronous_machine::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	WorkVec.ResizeReset(0);

	return WorkVec;
}

bool
asynchronous_machine_set(void)
{
#ifdef DEBUG
	std::cerr << __FILE__ <<":"<< __LINE__ << ":"<< __PRETTY_FUNCTION__ << std::endl;
#endif

	UserDefinedElemRead *rf = new UDERead<asynchronous_machine>;

	if (!SetUDE("asynchronous_machine", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" 
{

/**
 * \brief This function registers our user defined element for the math parser.
 *
 * It is called when the "module load" statement appears in the input file.
 */
int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!asynchronous_machine_set()) {
		silent_cerr("asynchronous_machine: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

} // extern "C"

#endif // ! STATIC_MODULES

