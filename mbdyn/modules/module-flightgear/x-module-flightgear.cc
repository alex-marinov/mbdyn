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

#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include "userelem.h"
#include "socketstream_out_elem.h"

// -I flightgear-$(FGVERSION)/src/ -I simgear-$(FGVERSION)/

// controls structure
#include "Network/net_ctrls.hxx"

// flight data model structure
#include "Network/net_fdm.hxx"

#if (FGVER != 20400) && (FGVER != 30200) && (FGVER != 30400)
#error "unsupported version"
#endif

class FGOutStreamContentModifier : public StreamContent::Modifier
{
private:
	size_t m_size;
	const char *m_buf;
	FGNetFDM m_fdm;

public:
	FGOutStreamContentModifier(void);
	virtual ~FGOutStreamContentModifier(void);

	virtual void Set(size_t size, const char *buf);
	virtual void Modify(void);

	virtual const void *GetOutBuf(void) const;
	virtual int GetOutSize(void) const;
};

FGOutStreamContentModifier::FGOutStreamContentModifier(void)
{
	NO_OP;
}

FGOutStreamContentModifier::~FGOutStreamContentModifier(void)
{
	NO_OP;
}

void
FGOutStreamContentModifier::Set(size_t size, const char *buf)
{
	m_size = size;
	m_buf = buf;
}

void
FGOutStreamContentModifier::Modify(void)
{
	m_fdm.version = FG_NET_FDM_VERSION;
}

const void *
FGOutStreamContentModifier::GetOutBuf(void) const
{
	return (void *)&m_fdm;
}

int
FGOutStreamContentModifier::GetOutSize(void) const
{
	return sizeof(m_fdm);
}

class FGOutElem : virtual public Elem, public UserDefinedElem {
public:
	friend class FGOutStreamContentModifier;

private:
	SocketStreamElem *m_pEl;

public: 
	FGOutElem(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~FGOutElem(void);

	virtual void SetValue(DataManager*, VectorHandler&, VectorHandler&, SimulationEntity::Hints*);

	virtual std::ostream& Restart(std::ostream&) const;
	virtual void WorkSpaceDim(integer*, integer*) const;
	virtual SubVectorHandler& AssRes(SubVectorHandler&, doublereal, const VectorHandler&, const VectorHandler&);
	virtual VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler&, doublereal, const VectorHandler&, const VectorHandler&);

	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void InitialWorkSpaceDim(integer*, integer*) const;
	virtual VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler&, const VectorHandler&);
	virtual SubVectorHandler& InitialAssRes(SubVectorHandler&, const VectorHandler&);
};


FGOutElem::FGOutElem(unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
UserDefinedElem(uLabel, pDO),
m_pEl(0)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
"\n"
"Module: 	flightgear\n"
"Author: 	Pierangelo Masarati <pierangelo.masarati@polimi.it>\n"
"Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali\n"
"		Politecnico di Milano\n"
"		http://www.aero.polimi.it\n"
"\n"
"	All rights reserved\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	// so far...
	NO_OP;

#if 0
   	SocketStreamElem(unsigned int uL, const std::string& name,
		unsigned int oe,
		UseSocket *pUS, StreamContent *pSC,
		int flags, bool bSendFirst, bool bAbortIfBroken,
		const std::string& sOutFileName, int iPrecision,
		doublereal dShift);
#endif

}

FGOutElem::~FGOutElem(void)
{
	if (m_pEl) {
		delete m_pEl;
		m_pEl = 0;
	}
}

void
FGOutElem::SetValue(DataManager*,
	VectorHandler&, VectorHandler&,
	SimulationEntity::Hints*)
{
	return;
}

std::ostream&
FGOutElem::Restart(std::ostream& out) const
{
	return out << "# FGOutElem(" << GetLabel() << "):"
		" restart not implement yet" << std::endl;
}

void
FGOutElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

SubVectorHandler&
FGOutElem::AssRes(SubVectorHandler& VH, doublereal, const VectorHandler&, const VectorHandler&)
{
	VH.Resize(0);
	return VH;
}

VariableSubMatrixHandler&
AssJac(VariableSubMatrixHandler& MH, doublereal, const VectorHandler&, const VectorHandler&)
{
	MH.SetNullMatrix();
	return MH;
}

unsigned int
FGOutElem::iGetInitialNumDof(void) const
{
	return 0;
}

void
FGOutElem::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
FGOutElem::InitialAssJac(VariableSubMatrixHandler& WM, const VectorHandler&)
{
	WM.SetNullMatrix();
	return WM;
}

SubVectorHandler&
FGOutElem::InitialAssRes(SubVectorHandler& VH, const VectorHandler&)
{
	VH.Resize(0);
	return VH;
}

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	UserDefinedElemRead *rf = new UDERead<FGOutElem>;

	if (!SetUDE("flightgear", rf)) {
		delete rf;
		return false;
	}

	return 0;
}

