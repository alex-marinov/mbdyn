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

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <errno.h>

/* ExtForce - begin */

/* Costruttore */
ExtForce::ExtForce(unsigned int uL,
	std::vector<StructNode *>& nodes,
	std::vector<Vec3>& offsets,
	std::string& fin,
	bool bRemoveIn,
        std::string& fout,
	bool bNoClobberOut,
	int iSleepTime,
	bool bTightCoupling,
	flag fOut)
: Elem(uL, Elem::FORCE, fOut), 
Force(uL, EXTERNALFORCE, 0, fOut), 
fin(fin.c_str()),
fout(fout.c_str()),
bRemoveIn(bRemoveIn),
bNoClobberOut(bNoClobberOut),
bTightCoupling(bTightCoupling),
bFirstRes(true),
iSleepTime(iSleepTime)
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
	/* If running tight coupling, send kinematics every iteration */
	/* NOTE: tight coupling may need relaxation */
	if (bTightCoupling) {
		Send();
	}
}
	
/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
	/* After prediction, mark next residual as first */
	bFirstRes = true;
}

/*
 * Elaborazione stato interno dopo la convergenza
 */
void
ExtForce::AfterConvergence(const VectorHandler& X, 
	const VectorHandler& XP)
{
	/* If not running tight coupling, send kinematics only at convergence */
	if (!bTightCoupling) {
		Send();
		if (bRemoveIn) {
			Unlink();
		}
	}
}

/*
 * Unlink input file when no longer required
 * used to inform * companion software that a new input file can be written
 */
void
ExtForce::Unlink(void)
{
	if (unlink(fin.c_str()) != 0) {
		int save_errno = errno;

		switch (save_errno) {
		case ENOENT:
			break;

		default:
			silent_cerr("ExtForce(" << GetLabel() << "): "
				<< "unable to delete input file \"" << fin.c_str() 
				<< "\": " << strerror(save_errno) << std::endl);
			throw ErrGeneric();
		}
	}
}

/*
 * Send output to companion software
 */
void
ExtForce::Send(void)
{
	if (bNoClobberOut) {
		bool	bKeepGoing(true);

		for (int cnt = 0; bKeepGoing; cnt++) {
			struct stat	s;

			if (stat(fout.c_str(), &s) != 0) {
				int save_errno = errno;

				switch (save_errno) {
				case ENOENT:
					bKeepGoing = false;
					break;

				default:
					silent_cerr("ExtForce(" << GetLabel() << "): "
						"unable to stat output file \"" << fout.c_str() << "\": "
						<< strerror(save_errno) << std::endl);
					throw ErrGeneric();
				}

			} else {
				silent_cout("ExtForce(" << GetLabel() << "): "
					"output file \"" << fout.c_str() << "\" still present, "
					"try #" << cnt << "; "
					"sleeping " << iSleepTime << " s" << std::endl);
				sleep(iSleepTime);
			}
		}
	}

	std::string tmpout(fout + ".tmp");
	std::ofstream	outf(tmpout.c_str());

	if (!outf) {
		silent_cerr("ExtForce(" << GetLabel() << "): "
			"unable to open file \"" << fout.c_str() << "\"" << std::endl);
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
	rename(tmpout.c_str(), fout.c_str());
}

SubVectorHandler&
ExtForce::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	if (bTightCoupling || bFirstRes) {
	        std::ifstream inf(fin.c_str());
		std::vector<bool> done(Nodes.size());

		for (unsigned int i = 0; i < Nodes.size(); i++) {
			done[i] = false;
		}

		for (int cnt = 0; !inf; cnt++) {
			silent_cout("ExtForce(" << GetLabel() << "): "
				"input file \"" << fin.c_str() << "\" missing, "
				"try #" << cnt << "; "
				"sleeping " << iSleepTime << " s" << std::endl); 
               
			sleep(iSleepTime);
			inf.clear();
			inf.open(fin.c_str());
		}

		WorkVec.ResizeReset(6*Nodes.size());

		for (int cnt = 0; inf; cnt++) {
			/* assume unsigned int label */
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
		}

		for (unsigned int i = 0; i < Nodes.size(); i++) {
			if (!done[i]) {
				silent_cerr("ExtForce(" << GetLabel() << "): node " << Nodes[i]->GetLabel() << " not done" << std::endl);
				throw ErrGeneric();
			}
		}

		if (bRemoveIn) {
			Unlink();
		}
	}

	bFirstRes = false;

	for (unsigned int i = 0; i < Nodes.size(); i++) {
		integer iFirstIndex = Nodes[i]->iGetFirstMomentumIndex();
		for (int r = 1; r <= 6; r++) {
			WorkVec.PutRowIndex(i*6 + r, iFirstIndex + r);
		}

		WorkVec.Add(i*6 + 1, F[i]);
		WorkVec.Add(i*6 + 4, M[i]);
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

	bool bUnlinkIn(false);
	if (HP.IsKeyWord("unlink")) {
		bUnlinkIn = true;
	}

	s = HP.GetFileName();
	if (s == 0) {
		silent_cerr("ExtForce(" << uLabel << "): unable to get output file name "
			"at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric();
	}
	std::string fout(s);

	bool bNoClobberOut(false);
	if (HP.IsKeyWord("no" "clobber")) {
		bNoClobberOut = true;
	}

	int iSleepTime = 1;
	if (HP.IsKeyWord("sleep" "time")) {
		iSleepTime = HP.GetInt();
		if (iSleepTime <= 0 ) {
			silent_cerr("ExtForce(" << uLabel << "): "
				"invalid sleep time " << iSleepTime <<std::endl);
			throw ErrGeneric();
		}
	}

	/* TODO: read loose/tight coupling flag */

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
		ExtForce(uLabel, Nodes, Offsets, fin, bUnlinkIn, fout, bNoClobberOut,
			iSleepTime, false, fOut));

	return pEl;
}

