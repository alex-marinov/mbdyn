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

/* 
 * These routines are part of DataManager, but are especially meant
 * to prepare textual output for ADAMS/View
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_ADAMS

#include "dataman.h"
#include "dataman_.h"

bool
DataManager::bAdamsOutput(void) const
{
	return (ResMode & RES_ADAMS) ? true : false;
}

void 
DataManager::AdamsResOutputInit(void)
{
	OutHdl.AdamsResOpen();
	std::ostream& out = OutHdl.AdamsRes();
	time_t t = time(NULL);

	/*
	 * conto i nodi non dummy reference frame
	 */
	unsigned int nStrNodes = 0;
	for (unsigned int i = 0; i < NodeData[Node::STRUCTURAL].iNum; i++) {
		StructNode *pStr = (StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

		/* Skip relative frame dummy nodes */
		if (pStr->GetStructNodeType() == StructNode::DUMMY) {
			DummyStructNode *pDmy = (DummyStructNode *)pStr;

			if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
				continue;
			}
		}

		nStrNodes++;
	}
   
	iAdamsOutputNodes = nStrNodes;
	iAdamsOutputParts = nStrNodes;

	Elem* p = NULL;
	if (ElemIter.bGetFirst(p)) {
		do {
			ASSERT(p != NULL);
			if (p != NULL) {
				iAdamsOutputParts += p->iGetNumDummyParts();
			}
		} while (ElemIter.bGetNext(p));
	}
   
	/* Cmd Init */
	OutHdl.AdamsCmdOpen();
	std::ostream& cmd = OutHdl.AdamsCmd();

	cmd << sAdamsModelName << std::endl;

#if 0
	cmd 
		<< "!" << std::endl
		<< "model create  &" << std::endl
		<< "   model_name = " << sAdamsModelName << std::endl
		<< "!" << std::endl
		<< "defaults model  &" << std::endl
		<< "   part_name = ground" << std::endl
		<< "!" << std::endl
		<< "defaults coordinate_system  &" << std::endl
		<< "   default_coordinate_system = ." << sAdamsModelName << ".ground" << std::endl;   
#endif

	out
		/* Header Block: File Type and Source Sub-Block */
		<< "FILE TYPE AND SOURCE                           6" << std::endl
		<< "ADAMS RESULTS FILE       M_KGS_N_RACA 0.100000000000000E+01" << std::endl
		<< "ADAMS VERSION - 10.0                    A.RE7" << std::endl
		<< "Mechanical Dynamics, Inc." << std::endl
		<< "Unknown" << std::endl
		<< ctime(&t)
		<< "!" << std::endl
     
		/* Header Block: Data Set Title Sub-Block */
		<< "DATA SET TITLE                                 3" << std::endl
		<< "ADAMS/View model name: " << (sSimulationTitle ? sSimulationTitle : "") << std::endl
		<< "" << std::endl
		<< "!" << std::endl

		/* Header Block: Results Set Sub-Block */
		<< "RESULTS TITLE                                  3" << std::endl
		<< "" << std::endl
		<< "" << std::endl
		<< "!" << std::endl
 
		/* Header Block: Number of Analysis Output Blocks Sub-Block 
		 * FIXME: replace 10000 with the number of steps
		 * actually performed
		 */
		<< "NUMBER OF ANALYSIS BLOCKS                      2" << std::endl;
	adamsNoab = out.tellp();

	const char *sVel = bAdamsVelocity ? "YES" : "NO ";
	const char *sAcc = bAdamsAcceleration ? "YES" : "NO ";
	out
		<< std::setw(8) << 2 << std::setw(8) << 123+3*iAdamsOutputParts << std::endl
		<< "!" << std::endl
     
		/* Rigid Parts Map Block */
		<< "RIGID PARTS MAP                         " << std::setw(8) << 2+iAdamsOutputParts << std::endl
		<< std::setw(8) << iAdamsOutputParts << "YES" << sVel << sAcc << std::endl;

	unsigned int i;
	unsigned int iVelOffset = iAdamsOutputParts;
	unsigned int iAccOffset = bAdamsVelocity ? 2*iAdamsOutputParts : iAdamsOutputParts;
	for (i = 0; i < nStrNodes; i++) {
		out 
			<< std::setw(20) << 2+i  /* 2 perche' 1 e' il "ground" in Adams! */
			<< std::setw(8) << 5+i
			<< std::setw(8) << 5+i+iVelOffset
			<< std::setw(8) << 5+i+iAccOffset << std::endl;
	}
   
	if (ElemIter.bGetFirst(p)) {
		do {
			ASSERT(p != NULL);
			if (p != NULL) {
				for (unsigned int part = 1; part <= p->iGetNumDummyParts(); part++, i++) {
					out
						<< std::setw(20) << 2+i
						<< std::setw(8) << 5+i
						<< std::setw(8) << 5+i+iVelOffset
						<< std::setw(8) << 5+i+iAccOffset << std::endl;
				}
			}
		} while (ElemIter.bGetNext(p));
	}

	out
		<< "!" << std::endl

		/* Point Mass Map Block */
		<< "POINT MASS MAP                                 2" << std::endl
		<< std::setw(8) << 0 << "YESNO NO" << std::endl
		<< "!" << std::endl

		/* Flexible Body Map Block */
		<< "FLEXIBLE BODY MAP                              2" << std::endl
		<< std::setw(8) << 0 << "YESNO NO" << std::endl
		<< "!" << std::endl

		/* User Differential Equations Map Block */
		<< "USER DIFFERENTIAL EQUATIONS MAP                2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Joints Map Block */
		<< "JOINTS MAP                                     2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Joint Primitives Map Block */
		<< "JOINT PRIMITIVES MAP                           2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Motion Generators Map Block */
		<< "MOTION GENERATORS MAP                          2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Gears Map Block */
		<< "GEARS MAP                                      2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl     
		<< "!" << std::endl

		/* Couplers Map Block */
		<< "COUPLERS MAP                                   2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* User Constraints Map Block */
		<< "USER CONSTRAINTS MAP                           2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Sforces Map Block */
		<< "SFORCES MAP                                    2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Spring Damper Forces Map Block */
		<< "SPRINGDAMPER FORCES MAP                        2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Bushings Map Block */
		<< "BUSHINGS MAP                                   2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Beams Map Block */
		<< "BEAMS MAP                                      2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Fields Map Block */
		<< "FIELDS MAP                                     2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* User Define Results Map Block */
		<< "USER DEFINED RESULTS MAP                       2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Variable map Block */
		<< "VARIABLE MAP                                   2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Pinput Map Block */
		<< "PINPUT MAP                                     2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Poutput Map Block */
		<< "POUTPUT MAP                                    2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Transfer Function SISO Map Block */
		<< "TFSISO MAP                                     2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* LSE Map Block */
		<< "LSE MAP                                        2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* GSE Map Block */
		<< "GSE MAP                                        2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* PTCV Map Block */
		<< "PTCV MAP                                       2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* CVCV Map Block */
		<< "CVCV MAP                                       2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Nforces Map Block */
		<< "NFORCES MAP                                    2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Vforces Map Block */
		<< "VFORCES MAP                                    2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Vtorques Map Block */
		<< "VTORQUES MAP                                   2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Gforces map Block */
		<< "GFORCES MAP                                    2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Frictions Map Block */
		<< "FRICTIONS MAP                                  2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Floating Markers Map Block */
		<< "FLOATING MARKERS MAP                           2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Mode Shapes Map Block */
		<< "MODE SHAPES MAP                                8" << std::endl
		<< std::setw(8) << 12 << std::setw(8) << 1 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl
		<< std::setw(8) << 0 << std::setw(8) << 13 << std::endl     
		<< "!" << std::endl

		/* Tires Map Block */
		<< "TIRES MAP                                      2" << std::endl
		<< std::setw(8) << 0 << "NO" << std::endl
		<< "!" << std::endl

		/* Name Map Block */
		<< "NAME MAP                                " << std::setw(8) << 2+2*iAdamsOutputParts << std::endl
		<< std::setw(8) << iAdamsOutputParts << "YES" << std::endl;

	for (i = 0; i < nStrNodes; i++) {
		StructNode *p = (StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];
		unsigned int l = p->GetLabel();
		out 
			<< std::setw(20) << 2+i << "RIGID PART                 1" << std::endl
			<< "." << sAdamsModelName << ".PART_" << l << std::endl;

		Vec3 x(p->GetXCurr());
		Vec3 e(MatR2EulerAngles(p->GetRCurr())*dRaDegr);

		/* Commento con il nome del nodo */
		const char *sName = p->GetName();
		if (sName != NULL) {
			cmd << "! --- " << sName << std::endl;
		} else {
			cmd << "! --- Structural Node " << l << std::endl;
		}
		cmd 
			<< "PART_" << l << std::endl
			<< 2+i << " " << x << " " << e << std::endl;

#if 0 
		cmd
			<< "!" << std::endl
			<< "part create rigid_body name_and_position  &" << std::endl
			<< "   part_name = ." << sAdamsModelName << ".PART_" << 1+i << "  &" << std::endl
			<< "   adams_id = " << 2+i << "  &" << std::endl
			<< "   location = ", x.Write(cmd, ", ") << "  &" << std::endl
			<< "   orientation = ", e.Write(cmd, ", ") << std::endl
			<< "!" << std::endl
			<< "marker create &" << std::endl
			<< "   marker_name = ." << sAdamsModelName << ".PART_" << 1+i << ".MAR_1  &" << std::endl
			<< "   adams_id = " << 2+i << "  &" << std::endl
			<< "   location = ", x.Write(cmd, ", ") << "  &" << std::endl
			<< "   orientation = ", e.Write(cmd, ", ") << std::endl
			<< "!" << std::endl
			<< "geometry create shape ellipsoid  &" << std::endl
			<< "   ellipsoid_name = ." << sAdamsModelName << ".PART_" << 1+i << ".SPHERE_1  &" << std::endl
			<< "   center_marker = ." << sAdamsModelName << ".PART_" << 1+i << ".MAR_1  &" << std::endl
			<< "   x_scale_factor = 0.1414213562  &" << std::endl
			<< "   y_scale_factor = 0.1414213562  &" << std::endl
			<< "   z_scale_factor = 0.1414213562" << std::endl
			<< "!" << std::endl
			<< "part attributes  &" << std::endl
			<< "   part_name = ." << sAdamsModelName << ".PART_" << 1+i << "  &" << std::endl
			<< "   color = MAIZE  &" << std::endl
			<< "   name_visibility = off" << std::endl;
#endif

	}

	if (ElemIter.bGetFirst(p)) {
		do {
			ASSERT(p != NULL);
			
			if (p != NULL) {
				unsigned int nParts = p->iGetNumDummyParts();
				
				if (nParts > 0) {
					unsigned int l = p->GetLabel();
					const char *s = psAdamsElemCode[p->GetElemType()];	     

					/* Nomi e dati delle parti */
					for (unsigned int part = 1; part <= p->iGetNumDummyParts(); part++, i++) {
						/* file intermedio per generazione .cmd */
						const char *sName = p->GetName();
						if (sName != NULL) {
							cmd
								<< "! --- " << sName 
								<< " (part " << part << ")" << std::endl;
						} else {
							cmd
								<< "! --- " << psElemNames[p->GetElemType()] 
								<< " " << l << " (part " << part << ")" << std::endl;
						}
						p->WriteAdamsDummyPartCmd(cmd, part, 2+i);

						/* nome parte */
						out 
							<< std::setw(20) << 2+i << "RIGID PART                 1" << std::endl
							<< "." << sAdamsModelName << "." << s << "_" << l << "_" << part << std::endl;
					}
				}
			}
		} while (ElemIter.bGetNext(p));
	}

	out 
		<< "!" << std::endl;
}

void 
DataManager::AdamsResOutput(integer iBlock, const char *type, const char *id) const
{
	std::ostream& out = OutHdl.AdamsRes();

	unsigned int nStrNodes = iAdamsOutputNodes;
	std::ios::fmtflags oldflags;

	oldflags = out.flags(std::ios::scientific);

	out 
		<< "ANALYSIS OUTPUT BLOCK                   " << std::setw(8) << 5+iAdamsOutputParts << std::endl
		<< std::setw(8) << iBlock;
	std::ios::fmtflags tmpflags;
	tmpflags = out.flags(std::ios::left);
	out << std::setw(20) << type;
	out.flags(tmpflags);
	out << id << std::endl;
	out << std::setw(12) << std::setprecision(5)
		<< dGetTime() << std::endl
		<< std::endl
		<< std::endl;

	for (unsigned int i = 0; i < nStrNodes; i++) {
		Vec3 x;
		doublereal e0;
		Vec3 e;

		StructNode *pNode =
			(StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

		/* Skip relative frame dummy nodes */
		if (pNode->GetStructNodeType() == StructNode::DUMMY) {
			DummyStructNode *pDmy = (DummyStructNode *)pNode;

			if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
				continue;
			}
		}
 
		x = pNode->GetXCurr();
		MatR2EulerParams(pNode->GetRCurr(), e0, e);
		out
			<< std::setw(12) << std::setprecision(5) << x(1)
			<< std::setw(12) << std::setprecision(5) << x(2)
			<< std::setw(12) << std::setprecision(5) << x(3)
			<< std::setw(12) << std::setprecision(5) << e0
			<< std::setw(12) << std::setprecision(5) << e(1)
			<< std::setw(12) << std::setprecision(5) << e(2)
			<< std::setw(12) << std::setprecision(5) << e(3)
			<< std::endl;
	}

	Elem *p = NULL;
	if (ElemIter.bGetFirst(p)) {
		do {
			ASSERT(p != NULL);
			
			if (p != NULL) {
				for (unsigned int part = 1; part <= p->iGetNumDummyParts(); part++) {
					Vec3 x;
					Mat3x3 R;
					doublereal e0;
					Vec3 e;
					
					p->GetDummyPartPos(part, x, R);
					MatR2EulerParams(R, e0, e);
					out 
						<< std::setw(12) << std::setprecision(5) << x(1)
						<< std::setw(12) << std::setprecision(5) << x(2)
						<< std::setw(12) << std::setprecision(5) << x(3)
						<< std::setw(12) << std::setprecision(5) << e0
						<< std::setw(12) << std::setprecision(5) << e(1)
						<< std::setw(12) << std::setprecision(5) << e(2)
						<< std::setw(12) << std::setprecision(5) << e(3)
						<< std::endl;
				}
			}
		} while (ElemIter.bGetNext(p));
	}

	if (bAdamsVelocity) {
		for (unsigned int i = 0; i < nStrNodes; i++) {
			Vec3 v;
			Vec3 w;

			StructNode *pNode =
				(StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

			/* Skip relative frame dummy nodes */
			if (pNode->GetStructNodeType() == StructNode::DUMMY) {
				DummyStructNode *pDmy = (DummyStructNode *)pNode;

				if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
					continue;
				}
			}
 
			v = pNode->GetVCurr();
			w = pNode->GetWCurr();
			
			out
				<< std::setw(12) << std::setprecision(5) << v(1)
				<< std::setw(12) << std::setprecision(5) << v(2)
				<< std::setw(12) << std::setprecision(5) << v(3)
				<< std::setw(12) << std::setprecision(5) << w(1)
				<< std::setw(12) << std::setprecision(5) << w(2)
				<< std::setw(12) << std::setprecision(5) << w(3)
				<< std::endl;
		}

		Elem *p = NULL;
		if (ElemIter.bGetFirst(p)) {
			do {
				ASSERT(p != NULL);
				
				if (p != NULL) {
					for (unsigned int part = 1; part <= p->iGetNumDummyParts(); part++) {
						Vec3 v, w;
						
						p->GetDummyPartVel(part, v, w);
						out 
							<< std::setw(12) << std::setprecision(5) << v.dGet(1)
							<< std::setw(12) << std::setprecision(5) << v.dGet(2)
							<< std::setw(12) << std::setprecision(5) << v.dGet(3)
							<< std::setw(12) << std::setprecision(5) << w.dGet(1)
							<< std::setw(12) << std::setprecision(5) << w.dGet(2)
							<< std::setw(12) << std::setprecision(5) << w.dGet(3)
							<< std::endl;
					}
				}
			} while (ElemIter.bGetNext(p));
		}
	}

	if (bAdamsAcceleration) {
		for (unsigned int i = 0; i < nStrNodes; i++) {
			Vec3 a(0.);
			Vec3 wp(0.);

			StructNode *pNode =
				(StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

			/* Skip relative frame dummy nodes */
			if (pNode->GetStructNodeType() == StructNode::DUMMY) {
				DummyStructNode *pDmy = (DummyStructNode *)pNode;

				if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
					continue;
				}
			} else if (pNode->GetStructNodeType() == StructNode::DYNAMIC) {
				DynamicStructNode *pDSN = dynamic_cast<DynamicStructNode *>(pNode);
				a = pDSN->GetXPPCurr();
				wp = pDSN->GetWPCurr();
			}
			
			out
				<< std::setw(12) << std::setprecision(5) << a(1)
				<< std::setw(12) << std::setprecision(5) << a(2)
				<< std::setw(12) << std::setprecision(5) << a(3)
				<< std::setw(12) << std::setprecision(5) << wp(1)
				<< std::setw(12) << std::setprecision(5) << wp(2)
				<< std::setw(12) << std::setprecision(5) << wp(3)
				<< std::endl;
		}

		Elem *p = NULL;
		if (ElemIter.bGetFirst(p)) {
			do {
				ASSERT(p != NULL);
				
				if (p != NULL) {
					for (unsigned int part = 1; part <= p->iGetNumDummyParts(); part++) {
						Vec3 a(0.), wp(0.);
						
						out 
							<< std::setw(12) << std::setprecision(5) << a.dGet(1)
							<< std::setw(12) << std::setprecision(5) << a.dGet(2)
							<< std::setw(12) << std::setprecision(5) << a.dGet(3)
							<< std::setw(12) << std::setprecision(5) << wp.dGet(1)
							<< std::setw(12) << std::setprecision(5) << wp.dGet(2)
							<< std::setw(12) << std::setprecision(5) << wp.dGet(3)
							<< std::endl;
					}
				}
			} while (ElemIter.bGetNext(p));
		}
	}

	out.flags(oldflags);

	out << "!" << std::endl;
}

void
DataManager::AdamsResOutputFini(void) const
{
	std::ostream& out = OutHdl.AdamsRes();
	out.seekp(adamsNoab, std::ios_base::beg);
	out << std::setw(8) << iOutputBlock;
}
#endif /* USE_ADAMS */

