/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <dataman.h>
#include "strext.h"

#include <fstream>

/* ExtForce - begin */

/* Costruttore */
ExtForce::ExtForce(unsigned int uL,
	std::vector<StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	std::string& fin,
	std::string& fout,
	flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
Force(uL, EXTERNALFORCE, 0, fOut), 
fin(fin.c_str()), fout(fout.c_str())
{
	ASSERT(nodes.size() == offsets.size());
	Nodes.resize(nodes.size());
	Offsets.resize(nodes.size());
	F.resize(nodes.size());
	M.resize(nodes.size());

	for (unsigned int i = 0; i < nodes.size(); i++) {
		Nodes[i] = nodes[i];
		Offsets[i] = offsets[i];
	}
}

ExtForce::~ExtForce(void)
{
	NO_OP;
}


void
ExtForce::Update(const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	std::ofstream	outf(fout.c_str());

	if (!outf) {
		silent_cerr("unable to open file " << fout.c_str() << std::endl);
		throw ErrGeneric();
	}

	for (unsigned int i = 0; i < Nodes.size(); i++) {
		Vec3 f = Nodes[i]->GetRCurr()*Offsets[i];
		Vec3 x = Nodes[i]->GetXCurr() + f;
		Vec3 v = Nodes[i]->GetVCurr() + Nodes[i]->GetWCurr().Cross(f);
		outf << Nodes[i]->GetLabel()
			<< " " << x
			<< " " << Nodes[i]->GetRCurr()
			<< " " << v
			<< " " << Nodes[i]->GetWCurr()
			<< std::endl;
	}

	/* send */
	outf.close();
}
	
/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	Update(X, XP);
}


SubVectorHandler&
ExtForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	std::ifstream inf(fin.c_str());
	std::vector<bool> done(Nodes.size());

	for (unsigned int i = 0; i < Nodes.size(); i++) {
		done[i] = false;
	}

	while (!inf) {
		sleep(1);
		inf.open(fin.c_str());
	}

	WorkVec.ResizeReset(6*Nodes.size());

	for (int cnt = 0; inf; cnt++) {
		/* assumo che la label sia un intero */
		unsigned l, i;
		doublereal f[3], m[3];

		inf >> l >> f[0] >> f[1] >> f[2] >> m[0] >> m[1] >> m[2];

		if (!inf) {
			break;
		}

		for (i = 0; i < Nodes.size(); i++) {
			if (Nodes[i]->GetLabel() == l) {
				break;
			}
		}

		if (i == Nodes.size()) {
			silent_cerr("ExtForce(" << GetLabel() << "): unknown label " << l << " as " << cnt << "-th node" << std::endl);
			throw ErrGeneric();
		}

		if (done[i]) {
			silent_cerr("ExtForce(" << GetLabel() << "): label " << l << " already done" << std::endl);
			throw ErrGeneric();
		}

		done[i] = true;
		F[i] = Vec3(f);
		M[i] = Vec3(m);

		integer iFirstIndex = Nodes[i]->iGetFirstMomentumIndex();
		for (int r = 1; r <= 6; r++) {
			WorkVec.PutRowIndex(i*6 + r, iFirstIndex + r);
		}

		WorkVec.Add(i*6 + 1, F[i]);
		WorkVec.Add(i*6 + 4, M[i]);
	}

	for (unsigned int i = 0; i < Nodes.size(); i++) {
		if (!done[i]) {
			silent_cerr("ExtForce(" << GetLabel() << "): node " << Nodes[i]->GetLabel() << " not done" << std::endl);
			throw ErrGeneric();
		}
	}

	return WorkVec;
}

void
ExtForce::Output(OutputHandler& OH) const
{
	std::ostream& out = OH.Forces();

	for (unsigned int i = 0; i < Nodes.size(); i++) {
		out << GetLabel() << "." << Nodes[i]->GetLabel()
			<< " " << F[i]
			<< " " << M[i]
			<< std::endl;
	}
}
   
Elem*
ReadExtForce(DataManager* pDM, 
	MBDynParser& HP, 
	unsigned int uLabel)
{
	const char	*s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): unable to get input file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	std::string fin(s);

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): unable to get output file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	std::string fout(s);

	int n = HP.GetInt();
	if (n <= 0) {
		silent_cerr("ExtForce(" << uLabel << "): illegal node number " << n <<
			" at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}

	std::vector<StructNode *> Nodes(n);
	std::vector<Vec3> Offsets(n);

	for (int i = 0; i < n; i++ ) {
		Nodes[i] = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);
		
		ReferenceFrame RF(Nodes[i]);

		if (HP.IsKeyWord("offset")) {
			Offsets[i] = HP.GetPosRel(RF);
		} else {
			Offsets[i] = Vec3(0.);
		}
	}

	flag fOut = pDM->fReadOutput(HP, Elem::FORCE);
	Elem *pEl = 0;
	SAFENEWWITHCONSTRUCTOR(pEl, ExtForce,
		ExtForce(uLabel, Nodes, Offsets, fin, fout, fOut));

	return pEl;
}

