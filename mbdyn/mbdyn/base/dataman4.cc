/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2008
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

/* Lettura elementi */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "mbdefs.h"
#include "modules.h"

#include <float.h>
#include <vector>
#include <set>

#include "dataman.h"
#include "dataman_.h"

/* Elementi */
#include "autostr.h"   /* Elementi automatici associati ai nodi dinamici */
#include "gravity.h"   /* Elemento accelerazione di gravita' */
#include "body.h"
#include "joint.h"
#include "jointreg.h"
#include "force.h"
#include "beam.h"
#include "beam2.h"
#include "hbeam.h"
#include "aerodyn.h"   /* Classe di base degli elementi aerodinamici */
#include "rotor.h"
#include "aeroelem.h"
#include "aeromodal.h"
#include "instruments.h"
#include "genfm.h"
#ifdef USE_EXTERNAL
#include "aeroext.h"
#endif /* USE_EXTERNAL */

#include "driven.h"    /* driven elements */
#include "genel.h"
#include "j2p.h"       /* bind di elementi a nodi parametrici */

#include "preselem.h"
#include "elec.h"
#include "therm.h"
#include "bulk.h"

#include "drive.h"
#include "tpldrive.h"

#include "loadable.h"

#ifdef USE_RTAI
#include "rtai_out_elem.h"
#endif /* USE_RTAI */

#ifdef USE_SOCKET
#include "socketstream_out_elem.h"
#include "socketstreammotionelem.h"
#endif // USE_SOCKET

static int iNumTypes[Elem::LASTELEMTYPE];

/* enum delle parole chiave */
enum KeyWords {
	UNKNOWNKEYWORD = -1,
	END = 0,
	ELEMENTS,

	GRAVITY,
	BODY,
	AUTOMATICSTRUCTURAL,
	JOINT,
	JOINT_REGULARIZATION,
	COUPLE,
	BEAM3,
	BEAM,		/* same as BEAM3 */
	BEAM2,
	HBEAM,

	AIRPROPERTIES,
	ROTOR,
	AERODYNAMICBODY,
	AERODYNAMICBEAM3,
	AERODYNAMICBEAM,	/* same as AERODYNAMICBEAM3 */
	AERODYNAMICBEAM2,
	AERODYNAMICEXTERNAL,
	AERODYNAMICEXTERNALMODAL,
	AEROMODAL,
	AIRCRAFTINSTRUMENTS,
	GENERICAERODYNAMICFORCE,

	FORCE,

	GENEL,
	ELECTRIC,

	THERMAL,

	HYDRAULIC,

	BULK,
	LOADABLE,
	DRIVEN,

	SOCKETSTREAM_OUTPUT,
	SOCKETSTREAM_MOTION_OUTPUT,
	RTAI_OUTPUT,

	INERTIA,

	EXISTING,
	OUTPUT,
	BIND,

	LASTKEYWORD
};

void
DataManager::ReadElems(MBDynParser& HP)
{
	DEBUGCOUTFNAME("DataManager::ReadElems");

	/* parole chiave del blocco degli elementi */
	const char* sKeyWords[] = {
		"end",

		"elements",

		"gravity",
		"body",
		"automatic" "structural",
		"joint",
			"joint" "regularization",
		"couple",
		"beam3",
		"beam",
		"beam2",
		"hbeam",

		"air" "properties",
		"rotor",
		"aerodynamic" "body",
		"aerodynamic" "beam3",
		"aerodynamic" "beam",
		"aerodynamic" "beam2",
		"aerodynamic" "external",
		"aerodynamic" "external" "modal",
		"aero" "modal",
		"aircraft" "instruments",
		"generic" "aerodynamic" "force",

		"force",

		"genel",
		"electric",

		"thermal",

		"hydraulic",

		"bulk",
		"loadable",
		"driven",

		"stream" "output",
		"stream" "motion" "output",
		"RTAI" "output",

		"inertia",

		"existing",
		"output",
		"bind",

		NULL
	};

	/* tabella delle parole chiave */
	KeyTable K(HP, sKeyWords);

	/* strutture di conteggio degli elementi letti */
	for (int i = 0; i < Elem::LASTELEMTYPE; iNumTypes[i] = ElemData[i].iExpectedNum, i++) {
		NO_OP;
	}

	int iMissingElems = Elems.size();
	DEBUGLCOUT(MYDEBUG_INPUT, "Expected elements: " << iMissingElems << std::endl);

	/* Aggiunta degli elementi strutturali automatici legati ai nodi dinamici */
	if (ElemData[Elem::AUTOMATICSTRUCTURAL].iExpectedNum > 0) {
		StructNode** ppTmpNod = (StructNode**)NodeData[Node::STRUCTURAL].ppFirstNode;
		int iTotNod = NodeData[Node::STRUCTURAL].iNum;

		ElemMapType &ElemMap = ElemData[Elem::AUTOMATICSTRUCTURAL].ElemMap;
		for (StructNode** ppTmp = ppTmpNod; ppTmp < ppTmpNod+iTotNod; ppTmp++) {
			if ((*ppTmp)->GetStructNodeType() == StructNode::DYNAMIC) {
				ASSERT(dynamic_cast<DynamicStructNode*>(*ppTmp) != 0);

				Elem *pTmpEl = 0;
				SAFENEWWITHCONSTRUCTOR(pTmpEl, AutomaticStructElem,
					AutomaticStructElem(dynamic_cast<DynamicStructNode*>(*ppTmp)));
				ElemMap.insert(ElemMapType::value_type((*ppTmp)->GetLabel(), pTmpEl));

				iMissingElems--;
				iNumTypes[Elem::AUTOMATICSTRUCTURAL]--;
				DEBUGLCOUT(MYDEBUG_INPUT,
					"Initializing automatic structural element linked to node "
					<< (*ppTmp)->GetLabel() << std::endl);
			}
		}
	}

	KeyWords CurrDesc;
	while ((CurrDesc = KeyWords(HP.GetDescription())) != END) {
		if (CurrDesc == OUTPUT) {
			DEBUGLCOUT(MYDEBUG_INPUT, "Elements to be output: ");
			Elem::Type Typ;
			switch (KeyWords(HP.GetWord())) {
			case BODY: {
				DEBUGLCOUT(MYDEBUG_INPUT, "bodies" << std::endl);
				Typ = Elem::BODY;
				break;
			}

			case AUTOMATICSTRUCTURAL: {
				DEBUGLCOUT(MYDEBUG_INPUT, "automatic structural" << std::endl);
				Typ = Elem::AUTOMATICSTRUCTURAL;
				break;
			}

			case JOINT: {
				DEBUGLCOUT(MYDEBUG_INPUT, "joints" << std::endl);
				Typ = Elem::JOINT;
				break;
			}

			case BEAM:
			case BEAM3:			/* same as BEAM */
			case BEAM2:
			case HBEAM: {
				DEBUGLCOUT(MYDEBUG_INPUT, "beams" << std::endl);
				Typ = Elem::BEAM;
				break;
			}

			case ROTOR: {
				DEBUGLCOUT(MYDEBUG_INPUT, "rotors" << std::endl);
				Typ = Elem::ROTOR;
				break;
			}

			case AEROMODAL: {
				DEBUGLCOUT(MYDEBUG_INPUT, "aeromodal" << std::endl);
				Typ = Elem::AEROMODAL;
				break;
			}

#ifdef USE_EXTERNAL
			case AERODYNAMICEXTERNAL:
			case AERODYNAMICEXTERNALMODAL: {
				DEBUGLCOUT(MYDEBUG_INPUT, "aerodynamic external" << std::endl);
				Typ = Elem::EXTERNAL;
				break;
			}
#endif /* USE_EXTERNAL */

			case AERODYNAMICBODY:
			case AERODYNAMICBEAM:
			case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
			case AERODYNAMICBEAM2:
			case AIRCRAFTINSTRUMENTS:
			case GENERICAERODYNAMICFORCE: {
				DEBUGLCOUT(MYDEBUG_INPUT, "aerodynamic" << std::endl);
				Typ = Elem::AERODYNAMIC;
				break;
			}

			case FORCE:
			case COUPLE:
			{
				DEBUGLCOUT(MYDEBUG_INPUT, "forces" << std::endl);
				Typ = Elem::FORCE;
				break;
				}

			case GENEL: {
				DEBUGLCOUT(MYDEBUG_INPUT, "genels" << std::endl);
				Typ = Elem::GENEL;
				break;
			}

			case ELECTRIC: {
				DEBUGLCOUT(MYDEBUG_INPUT, "electric" << std::endl);
				Typ = Elem::ELECTRIC;
				break;
			}

			case THERMAL: {
				DEBUGLCOUT(MYDEBUG_INPUT, "thermal" << std::endl);
				Typ = Elem::THERMAL;
				break;
			}

			case HYDRAULIC: {
				DEBUGLCOUT(MYDEBUG_INPUT, "hydraulic elements" << std::endl);
				Typ = Elem::HYDRAULIC;
				break;
			}

			case BULK: {
				DEBUGLCOUT(MYDEBUG_INPUT, "bulk" << std::endl);
				Typ = Elem::BULK;
				break;
			}

			case LOADABLE: {
				DEBUGLCOUT(MYDEBUG_INPUT, "loadable" << std::endl);
				Typ = Elem::LOADABLE;
				break;
			}

			case UNKNOWNKEYWORD: {
				silent_cerr("Error: unknown element type, cannot modify output"
					<< std::endl);
				throw DataManager::ErrGeneric();
			}

			default: {
				silent_cerr("Error: element type " << sKeyWords[CurrDesc]
					<< " at line " << HP.GetLineData() << " is not allowed"
					<< std::endl);
				throw DataManager::ErrGeneric();
			}
			} /* end switch (KeyWords(HP.GetWord()))  */

			/* Elements list */
			while (HP.IsArg()) {
				if (HP.IsKeyWord("range")) {
					unsigned int uL = (unsigned int)HP.GetInt();
					unsigned int uEndL = (unsigned int)HP.GetInt();
					if (uEndL < uL) {
						silent_cerr("End label " << uEndL
							<< " must be larger than "
							"or equal to start label "
							<< uL << std::endl);
						throw ErrGeneric();
					}
					for ( ; uL <= uEndL; uL++) {
						Elem* pE = (Elem*)pFindElem(Typ, uL);
						if (pE != 0) {
							DEBUGLCOUT(MYDEBUG_INPUT, "element " << uL << std::endl);
							pE->SetOutputFlag(flag(1));
						}
					}

				} else {
					unsigned int uL = (unsigned int)HP.GetInt();
					Elem* pE = (Elem*)pFindElem(Typ, uL);
					if (pE == 0) {
						silent_cerr(psElemNames[Typ] << "(" << uL << ") "
							"is not defined; output cannot be modified"
							<< std::endl);
					} else {
						DEBUGLCOUT(MYDEBUG_INPUT, "element " << uL << std::endl);
						pE->SetOutputFlag(flag(1));
					}
				}
			} /* end while (HP.IsArg()) */

		} else if (CurrDesc == INERTIA) {
			unsigned int uIn = (unsigned int)HP.GetInt();

			/* Nome dell'elemento */
			const char *sName = NULL;
			if (HP.IsKeyWord("name")) {
				const char *sTmp = HP.GetStringWithDelims();
				SAFESTRDUP(sName, sTmp);
			}

			Vec3 x(Zero3);
			if (HP.IsKeyWord("position")) {
				x = HP.GetPosAbs(AbsRefFrame);
			}

			Mat3x3 R(Eye3);
			Mat3x3 RT(Eye3);
			if (HP.IsKeyWord("orientation")) {
				R = HP.GetRotAbs(AbsRefFrame);
				RT = R.Transpose();
			}

			doublereal dM(0.);
			Vec3 S(0.);
			Mat3x3 J(0.);

			std::set<unsigned int> Body_labels;
			Elem::Type Type = Elem::UNKNOWN;
			bool bOut(false);
			bool bLog(true);
			while (HP.IsArg()) {
				if (HP.IsKeyWord("output")) {
					if (HP.IsKeyWord("no")) {
						bLog = false;

					} else if (HP.IsKeyWord("yes")) {
						bLog = false;
						bOut = true;

					} else if (HP.IsKeyWord("log")) {
						/* bLog = true */ ;

					} else if (HP.IsKeyWord("both")) {
						bOut = true;

					} else {
						silent_cerr("unknown output mode for inertia "
							<< uIn << " at line "
							<< HP.GetLineData()
							<< std::endl);
						throw ErrGeneric();
					}
					break;

				} else if (HP.IsKeyWord(sKeyWords[BODY])) {
					Type = Elem::BODY;

				} else if (HP.IsKeyWord(sKeyWords[JOINT])) {
					Type = Elem::JOINT;

				} else if (HP.IsKeyWord(sKeyWords[LOADABLE])) {
					Type = Elem::LOADABLE;

#if 0
				} else if (HP.IsKeyWord("...")) {
					/* other types with inertia */
#endif
				}

				if (Type == Elem::UNKNOWN) {
					silent_cerr("inertia " << uIn);
					if (sName) {
						silent_cerr(" (" << sName << ")");
					}
					silent_cerr(" at line " << HP.GetLineData()
						<< ": unknown or undefined element type"
						<< std::endl);
						throw ErrGeneric();
				}

				/*
				 * FIXME: duplicate check?
				 */

				if (HP.IsKeyWord("all")) {
					for (ElemMapType::const_iterator i = ElemData[Type].ElemMap.begin();
						i != ElemData[Type].ElemMap.end();
						i++)
					{
						unsigned int uL = i->second->GetLabel();
						std::set<unsigned int>::const_iterator BL_end = Body_labels.end();

						if (Body_labels.find(uL) == BL_end) {
							Body_labels.insert(uL);
							ElemGravityOwner *pEl = dynamic_cast<ElemGravityOwner *>(i->second);

							dM += pEl->dGetM();
							S += pEl->GetS();
							J += pEl->GetJ();

						} else {
							silent_cerr(psElemNames[Type]
								<< "(" << uL
								<< ") duplicate label at line "
								<< HP.GetLineData()
								<< " (ignored)" << std::endl);
						}
					}

				} else {
					unsigned int uL = (unsigned int)HP.GetInt();
					std::set<unsigned int>::const_iterator BL_end = Body_labels.end();
					if (Body_labels.find(uL) == BL_end) {
						Elem **ppTmpEl = (Elem **)ppFindElem(Type, uL);
						if (ppTmpEl == NULL || ppTmpEl[0] == NULL) {
							silent_cerr("inertia " << uIn
								<< " at line " << HP.GetLineData()
								<< ": unable to find " << psElemNames[Type]
								<< "(" << uL << ")" << std::endl);
							throw ErrGeneric();
						}

						Body_labels.insert(uL);
						ElemGravityOwner *pEl = dynamic_cast<ElemGravityOwner *>(ppTmpEl[0]);
						ASSERT(pEl != 0);

						/* FIXME: temporary... */
						if (Type == Elem::JOINT) {
							Joint *pJ = dynamic_cast<Joint *>(pEl);
							if (pJ == 0 || pJ->GetJointType() != Joint::MODAL) {
								silent_cerr("warning, Joint(" << uL << ") "
									"is not modal" << std::endl);
							}
						}

						dM += pEl->dGetM();
						S += pEl->GetS();
						J += pEl->GetJ();

					} else {
						silent_cerr(psElemNames[Type] << "(" << uL
							<< "): duplicate label at line "
							<< HP.GetLineData() << " (ignored)"
							<< std::endl);
					}
				}
			}  /* end while (HP.IsArg()) */ 

			Vec3 Xcg(0.);
			Mat3x3 Jcg(J);
			if (dM < DBL_EPSILON) {
				silent_cerr("inertia " << uIn
					<< " at line " << HP.GetLineData()
					<< ": mass is null" << std::endl);
			} else {
				Xcg = S/dM;

				/*
				 * FIXME: should also rotate it in the principal
				 * reference frame, and log the angles
				 */
				Jcg += Mat3x3(S, Xcg);
			}

			if (!x.IsNull()) {
				Vec3 Dx = Xcg - x;
				J += Mat3x3(x, x*dM) + Mat3x3(Dx, x*dM) + Mat3x3(x, Dx*dM);
			}

			if (bLog) {
				OutHdl.Log()
					<< "inertia: " << uIn
					<< " (" << ( sName ? sName : "unnamed" ) << ")"
					<< std::endl
					<< "    mass:        " << dM << std::endl
					<< "    Xcg:         " << Xcg << std::endl
					<< "    Xcg-X:       " << (Xcg - x) << std::endl
					<< "    R^T*(Xcg-X): " << RT*(Xcg - x) << std::endl
					<< "    J:           " << RT*J*R << std::endl
					<< "    Jcg:         " << RT*Jcg*R << std::endl;
			}
			if (bOut) {
				silent_cout("inertia: " << uIn
					<< " (" << ( sName ? sName : "unnamed" ) << ")"
					<< std::endl
					<< "    mass:        " << dM << std::endl
					<< "    Xcg:         " << Xcg << std::endl
					<< "    Xcg-X:       " << (Xcg - x) << std::endl
					<< "    R^T*(Xcg-X): " << RT*(Xcg - x) << std::endl
					<< "    J:           " << RT*J*R << std::endl
					<< "    Jcg:         " << RT*Jcg*R << std::endl);
			}
			if (sName) {
				SAFEDELETEARR(sName);
			}

		} else if (CurrDesc == BIND) {
			/* Label dell'elemento */
			unsigned int uL = HP.GetInt();

			/* Tipo dell'elemento */
			Elem::Type t = Elem::UNKNOWN;
			switch (KeyWords(HP.GetWord())) {
			case BODY:
				t = Elem::BODY;
				break;

			case AUTOMATICSTRUCTURAL:
				t = Elem::AUTOMATICSTRUCTURAL;
				break;

			case JOINT:
				t = Elem::JOINT;
				break;

			case FORCE:
			case COUPLE:
				t = Elem::FORCE;
				break;

			case BEAM:
			case BEAM3:			/* same as BEAM */
			case BEAM2:
			case HBEAM:
				t = Elem::BEAM;
				break;

			case ROTOR:
				t = Elem::ROTOR;
				break;

			case AEROMODAL:
				t = Elem::AEROMODAL;
				break;

#ifdef USE_EXTERNAL
			case AERODYNAMICEXTERNAL:
			case AERODYNAMICEXTERNALMODAL:
				t = Elem::EXTERNAL;
				break;
#endif /* USE_EXTERNAL */

			case AERODYNAMICBODY:
			case AERODYNAMICBEAM:
			case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
			case AERODYNAMICBEAM2:
			case AIRCRAFTINSTRUMENTS:
			case GENERICAERODYNAMICFORCE:
				t = Elem::AERODYNAMIC;
				break;

			case GENEL:
				t = Elem::GENEL;
				break;

			case ELECTRIC:
				t = Elem::ELECTRIC;
				break;

			case THERMAL:
				t = Elem::THERMAL;
				break;

			case HYDRAULIC:
				t = Elem::HYDRAULIC;
				break;

			case BULK:
				t = Elem::BULK;
				break;

			case LOADABLE:
				t = Elem::LOADABLE;
				break;

			case RTAI_OUTPUT:
			case SOCKETSTREAM_OUTPUT:
			case SOCKETSTREAM_MOTION_OUTPUT:
				silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
					<< " does not support bind" << std::endl);
			default:
				throw ErrGeneric();
			} /* end switch (KeyWords(HP.GetWord())) */

			Elem* pEl = dynamic_cast<Elem*>(pFindElem(t, uL));
			if (pEl == 0) {
				silent_cerr("can't find " << psElemNames[t] << " (" << uL << ") "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}

			/* Label del nodo parameter */
			uL = HP.GetInt();

			Elem2Param* pNd = dynamic_cast<Elem2Param*>(pFindNode(Node::PARAMETER, uL));
			if (pNd == 0) {
				silent_cerr("can't find parameter node (" << uL << ") "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}

			/* Numero d'ordine del dato privato a cui fare il binding */
			unsigned int i = 0;
			if (HP.IsKeyWord("name") || HP.IsKeyWord("string" /* deprecated */ )) {
				const char *s = HP.GetStringWithDelims();

				ASSERT(s != NULL);

				DEBUGCOUT("binding to " << psElemNames[pEl->GetElemType()]
					<< "(" << pEl->GetLabel() << ") private data "
					"\"" << s << "\"" << std::endl);

				i = pEl->iGetPrivDataIdx(s);

			} else {
				i = HP.GetInt();
			}

			/* indice del dato a cui il parametro e' bound */
			if (i <= 0 || i > pEl->iGetNumPrivData()) {
				silent_cerr("error in private data number " << i << " for element "
					<< psElemNames[t] << " (" << pEl->GetLabel() << ") "
					"at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric();
			}

			/* fa il binding del ParameterNode all'elemento */
			DEBUGLCOUT(MYDEBUG_INPUT, "Binding " << psElemNames[t]
				<< " (" << pEl->GetLabel() << ") "
				"to Parameter " << pNd->GetLabel() << std::endl);
			pNd->Bind(pEl, i);


		/* gestisco a parte gli elementi automatici strutturali, perche'
		 * sono gia' stati costruiti altrove e li devo solo inizializzare;
		 * eventualmente si puo' fare altrimenti */
		} else if (CurrDesc == AUTOMATICSTRUCTURAL) {
			unsigned int uLabel = HP.GetInt();
			Elem* pEl = pFindElem(Elem::AUTOMATICSTRUCTURAL, uLabel);
			if (pEl == 0) {
				silent_cerr("line " << HP.GetLineData()
					<< ": unable to find automatic structural element "
					<< uLabel << std::endl);
				throw ErrGeneric();
			}

			DEBUGCOUT("reading automatic structural element " << uLabel << std::endl);

			/* forse e' il caso di usare il riferimento del nodo? */

			/* NOTE: i primi due sono gestiti direttamente
			 * dagli elementi con inerzia, e quindi non sono usati */
			Vec3 q(HP.GetVecAbs(AbsRefFrame));
			Vec3 g(HP.GetVecAbs(AbsRefFrame));
			Vec3 qp(HP.GetVecAbs(AbsRefFrame));
			Vec3 gp(HP.GetVecAbs(AbsRefFrame));

			DEBUGCOUT("Q  = " << q << std::endl
				<< "G  = " << g << std::endl
				<< "Qp = " << qp << std::endl
				<< "Gp = " << gp << std::endl);

			AutomaticStructElem* pAuto = dynamic_cast<AutomaticStructElem *>(pEl);
			if (pAuto == 0) {
				silent_cerr("line " << HP.GetLineData()
					<< ": Elem(" << uLabel << ") is not AutomaticStructural"
					<< std::endl);
				throw ErrGeneric();
			}
			pAuto->Init(q, g, qp, gp);

        /*  <<<<  D E F A U L T  >>>>  :  Read one element and create it */ 
		/* default: leggo un elemento e lo creo */
		} else {              			
			/* puntatore all'elemento */
			Elem* pE = 0;

			unsigned int uLabel;
			switch (CurrDesc) {
			/* Qui vengono elencati gli elementi unici, che non richiedono label
			 * (per ora: accelerazione di gravita' e proprieta' dell'aria */

			/* Accelerazione di gravita' */
			case GRAVITY: {
				silent_cout("Reading gravity acceleration" << std::endl);

				if (iNumTypes[Elem::GRAVITY]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData()
						<< ": gravity acceleration element is not defined"
						<< std::endl);

					throw DataManager::ErrGeneric();
				}

				uLabel = 1;

				TplDriveCaller<Vec3>* pDC = ReadDC3D(this, HP);

				flag fOut = fReadOutput(HP, Elem::GRAVITY);

				SAFENEWWITHCONSTRUCTOR(pE, Gravity, Gravity(pDC, fOut));

				ElemData[Elem::GRAVITY].ElemMap.insert(ElemMapType::value_type(uLabel, pE));

				break;
			}

			/* Elementi aerodinamici: proprieta' dell'aria */
			case AIRPROPERTIES: {
				silent_cout("Reading air properties" << std::endl);

				if (iNumTypes[Elem::AIRPROPERTIES]-- <= 0) {
					DEBUGCERR("");
					silent_cerr("line " << HP.GetLineData()
						<< ": air properties element is not defined"
						<< std::endl);

					throw DataManager::ErrGeneric();
				}

				uLabel = 1;

				pE = ReadAirProperties(this, HP);
				if (pE != 0) {
					ElemData[Elem::AIRPROPERTIES].ElemMap.insert(ElemMapType::value_type(uLabel, pE));
				}

				break;
			}

			/* Elemento generico: legge la label e fa un nuovo switch */
			default: {
				/* legge la label */
				uLabel = unsigned(HP.GetInt());

				/* in base al tipo, avviene l'allocazione */
				switch (CurrDesc) {
				case DRIVEN: {
					/* Reads the driver */
					DriveCaller* pDC = HP.GetDriveCaller();
					std::vector<std::string> hints;

					for (unsigned i = 1; HP.IsKeyWord("hint"); i++) {
						const char *hint = HP.GetStringWithDelims(HighParser::DEFAULTDELIM, false);
						if (hint == 0) {
							silent_cerr("Driven(" << uLabel << "): "
								"unable to read hint #" << i
								<< " at line " << HP.GetLineData()
								<< std::endl);
							throw ErrGeneric();
						}
						hints.push_back(hint);
					}

					HP.ExpectDescription();
					KeyWords CurrDriven = KeyWords(HP.GetDescription());

#ifdef DEBUG
					switch (CurrDriven) {
					case FORCE:
					case BODY:
					case JOINT:
					case JOINT_REGULARIZATION:
					case COUPLE:
					case BEAM:
					case BEAM3:			/* same as BEAM */
					case BEAM2:
					case HBEAM:

					case ROTOR:
					case AERODYNAMICBODY:
					case AERODYNAMICBEAM:
					case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
					case AERODYNAMICBEAM2:
					case AIRCRAFTINSTRUMENTS:
					case GENERICAERODYNAMICFORCE:

					case GENEL:
					case ELECTRIC:

					case THERMAL:

					case HYDRAULIC:

					case BULK:
					case LOADABLE:
					case EXISTING:
						DEBUGLCOUT(MYDEBUG_INPUT, "OK, this element can be driven" << std::endl);
						break;

					case RTAI_OUTPUT:
					case SOCKETSTREAM_OUTPUT:
					case SOCKETSTREAM_MOTION_OUTPUT:
						silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
							<< " cannot be driven" << std::endl);
						throw ErrGeneric();

					default:
						DEBUGCERR("warning, this element can't be driven" << std::endl);
						break;
					}
#endif /* DEBUG */

					Elem **ppE = 0;
					if (CurrDriven == EXISTING) {
						iMissingElems++;
						CurrDriven = KeyWords(HP.GetWord());
						unsigned int uL = (unsigned int)HP.GetInt();
						if (uL != uLabel) {
							silent_cerr("Error: the driving element must have "
								"the same label of the driven"
								<< std::endl);

							throw DataManager::ErrGeneric();
						}

						/* FIXME: use pFindElem() instead? */
						switch (CurrDriven) {
						case FORCE:
							ppE = ppFindElem(Elem::FORCE, uLabel);
							break;

						case BODY:
							ppE = ppFindElem(Elem::BODY, uLabel);
							break;

						case JOINT:
							ppE = ppFindElem(Elem::JOINT, uLabel);
							break;

						case JOINT_REGULARIZATION:
							ppE = ppFindElem(Elem::JOINT_REGULARIZATION, uLabel);
							break;

						case COUPLE:
							ppE = ppFindElem(Elem::FORCE, uLabel);
							break;

						case BEAM:
						case BEAM3:		/* same as BEAM */
						case BEAM2:
						case HBEAM:
							ppE = ppFindElem(Elem::BEAM, uLabel);
							break;

						case ROTOR:
							ppE = ppFindElem(Elem::ROTOR, uLabel);
							break;

						case AERODYNAMICBODY:
						case AERODYNAMICBEAM:
						case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
						case AERODYNAMICBEAM2:
						case AIRCRAFTINSTRUMENTS:
						case GENERICAERODYNAMICFORCE:
							ppE = ppFindElem(Elem::AERODYNAMIC, uLabel);
							break;

						case GENEL:
							ppE = ppFindElem(Elem::GENEL, uLabel);
							break;

						case ELECTRIC:
							ppE = ppFindElem(Elem::ELECTRIC, uLabel);
							break;

						case THERMAL:
							ppE = ppFindElem(Elem::THERMAL, uLabel);
							break;

						case HYDRAULIC:
							ppE = ppFindElem(Elem::HYDRAULIC, uLabel);
							break;

						case BULK:
							ppE = ppFindElem(Elem::BULK, uLabel);
							break;

						case LOADABLE:
							ppE = ppFindElem(Elem::LOADABLE, uLabel);
							break;

						default:
							DEBUGCERR("warning, this element can't be driven" << std::endl);
							break;
							
						}  /*switch (CurrDriven) */

						if (ppE == NULL) {
							silent_cerr("Error: element " << uLabel
								<< "cannot be driven since it doesn't exist"
								<< std::endl);

							throw DataManager::ErrGeneric();
						}

						pE = *ppE;

						flag fOut = fReadOutput(HP, pE->GetElemType());
						pE->SetOutputFlag(fOut);

					} else {
						unsigned int uDummy = (unsigned int)HP.GetInt();
						if (uDummy != uLabel) {
							silent_cerr("Error: the element label "
								"must be the same of the driving element"
								<< std::endl);

							throw DataManager::ErrGeneric();
						}

						/* Reads the true element */
						ppE = ReadOneElem(HP, uLabel, CurrDriven);
						if (ppE == NULL) {
							DEBUGCERR("");
							silent_cerr("error in allocation of element "
								<< uLabel << std::endl);

							throw ErrMemory();
						}

						pE = *ppE;
					}  /*  if (CurrDriven == EXISTING) {..} else {..}  */

					SimulationEntity::Hints *pHints = 0;
					if (!hints.empty()) {
						for (std::vector<std::string>::const_iterator i = hints.begin();
								i != hints.end();
								i++)
						{
							Hint *ph = pE->ParseHint(this, i->c_str());
							if (ph != 0) {
								if (pHints == 0) {
									pHints = new SimulationEntity::Hints;
								}
								pHints->push_back(ph);
							}
						}
					}

					/* Creates the driver for the element */
					Elem* pEl = 0;
					SAFENEWWITHCONSTRUCTOR(pEl,
						DrivenElem,
						DrivenElem(this, pDC, pE, pHints));

					/* Substitutes the element with the driver */
					pE = *ppE = pEl;

					break;
				}  /* end case DRIVEN: */

                /*  <<<<  N O R M A L   E L E M E N T S  >>>>>  */
				/* Normal element */
				case FORCE:

				case BODY:
				case JOINT:
				case JOINT_REGULARIZATION:
				case COUPLE:
				case BEAM:
				case BEAM3:		/* same as BEAM */
				case BEAM2:
				case HBEAM:

				case ROTOR:
				case AERODYNAMICBODY:
				case AERODYNAMICBEAM:
				case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
				case AERODYNAMICBEAM2:
				case AIRCRAFTINSTRUMENTS:
				case GENERICAERODYNAMICFORCE:
#ifdef USE_EXTERNAL
				case AERODYNAMICEXTERNAL:
				case AERODYNAMICEXTERNALMODAL:
#endif /* USE_EXTERNAL */
				case AEROMODAL:

				case GENEL:
				case ELECTRIC:

				case THERMAL:

				case HYDRAULIC:

				case BULK:
				case LOADABLE:
				case RTAI_OUTPUT:
				case SOCKETSTREAM_OUTPUT:
				case SOCKETSTREAM_MOTION_OUTPUT:
				{
					Elem **ppE = 0;

					/* Nome dell'elemento */
					const char *sName = NULL;
					if (HP.IsKeyWord("name")) {
						const char *sTmp = HP.GetStringWithDelims();
						SAFESTRDUP(sName, sTmp);
					}

#ifdef USE_RUNTIME_LOADING
					if (CurrDesc == LOADABLE && !loadableElemInitialized) {
						module_initialize();
						loadableElemInitialized = true;
					}
#endif // USE_RUNTIME_LOADING

					ppE = ReadOneElem(HP, uLabel, CurrDesc);

					if (ppE != 0 ) {
						pE = *ppE;
						if (sName != 0) {
							pE->PutName(sName);
							SAFEDELETEARR(sName);
						}
					}

					break;
				}  /* end case 'Normal elements'*/

				/* in caso di tipo sconosciuto */
				case UNKNOWNKEYWORD: {
					DEBUGCERR("");
					silent_cerr("error - unknown element type at line "
						<< HP.GetLineData() << std::endl);

					throw DataManager::ErrGeneric();
				}

				default: {
					DEBUGCERR("");
					silent_cerr("error - element type " << sKeyWords[CurrDesc]
						<< " at line " << HP.GetLineData()
						<< " is not allowed " << std::endl);

					throw DataManager::ErrGeneric();
				}
				}  /* end switch (CurrDesc) 'in base al tipo'  */
				
			}  /* end case default: 'Elemento generico' */
			}  /* end switch (CurrDesc) 'Elemento generico' */

			/* verifica dell'allocazione */
			ASSERT(pE != 0);

			/* Aggiorna le dimensioni massime degli spazi di lavoro
			 * (qui va bene perche' il puntatore e' gia' stato verificato) */
			integer iNumRows = 0;
			integer iNumCols = 0;
			pE->WorkSpaceDim(&iNumRows, &iNumCols);
			if (iNumRows > iMaxWorkNumRows) {
				iMaxWorkNumRows = iNumRows;
				DEBUGLCOUT(MYDEBUG_INIT, "Current max work rows number: "
					<< iMaxWorkNumRows << std::endl);
			}
			if (iNumCols > iMaxWorkNumCols) {
				iMaxWorkNumCols = iNumCols;
				DEBUGLCOUT(MYDEBUG_INIT, "Current max work cols number: "
					<< iMaxWorkNumCols << std::endl);
			}

			/* decrementa il totale degli elementi mancanti */
			iMissingElems--;
			
		}  /* end <<<<  D E F A U L T  >>>>  :  Read one element and create it */ 
		   /* end default: leggo un elemento e lo creo */
		 
	}  /* while ((CurrDesc = KeyWords(HP.GetDescription())) != END) */

	if (KeyWords(HP.GetWord()) != ELEMENTS) {
		DEBUGCERR("");
		silent_cerr("<end: elements;> expected at line"
			<< HP.GetLineData() << std::endl);

		throw DataManager::ErrGeneric();
	}

	/* Se non c'e' il punto e virgola finale */
	if (HP.IsArg()) {
		DEBUGCERR("");
		silent_cerr("semicolon expected at line " << HP.GetLineData()
			<< std::endl);

		throw DataManager::ErrGeneric();
	}

	if (iMissingElems > 0) {
		DEBUGCERR("");
		silent_cerr("warning: " << iMissingElems
			<< " elements are missing;" << std::endl);
		for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
			if (iNumTypes[iCnt] > 0) {
				silent_cerr("  " << iNumTypes[iCnt]
					<< ' ' << psElemNames[iCnt] << std::endl);
			}
		}

		throw DataManager::ErrGeneric();
	}

	/* Linka gli elementi che generano forze d'inerzia all'elemento
	 * accelerazione di gravita' */
	if (!ElemData[Elem::GRAVITY].ElemMap.empty()) {
		Gravity* pGrav = dynamic_cast<Gravity *>(ElemData[Elem::GRAVITY].ElemMap.begin()->second);
		ASSERT(pGrav != 0);

		for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
			if (ElemData[iCnt].bGeneratesInertiaForces()
				&& !ElemData[iCnt].ElemMap.empty())
			{
				for (ElemMapType::const_iterator p = ElemData[iCnt].ElemMap.begin();
					p != ElemData[iCnt].ElemMap.end();
					p++)
				{
					ElemGravityOwner *pGO = CastElemGravityOwner(p->second);

					ASSERT(pGO != 0);
					pGO->PutGravity(pGrav);
				}
			}
		}
	}

	/* Linka gli elementi che usano le proprieta' dell'aria all'elemento
	 * proprieta' dell'aria */
	if (!ElemData[Elem::AIRPROPERTIES].ElemMap.empty()) {
		AirProperties* pProp = dynamic_cast<AirProperties *>(ElemData[Elem::AIRPROPERTIES].ElemMap.begin()->second);
		ASSERT(pProp != 0);

		for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
			if (ElemData[iCnt].bUsesAirProperties()
				&& !ElemData[iCnt].ElemMap.empty())
			{
				for (ElemMapType::const_iterator p = ElemData[iCnt].ElemMap.begin();
					p != ElemData[iCnt].ElemMap.end();
					p++)
				{
					AerodynamicElem *pAE = CastAerodynamicElem(p->second);

					ASSERT(pAE != 0);
					pAE->PutAirProperties(pProp);
				}
			}
		}

	} else { /* '' */
		/* Esegue un controllo per vedere se esistono elementi aerodinamici
		 * ma non sono definite le proprieta' dell'aria, nel qual caso
		 * il calcolo deve essere arrestato */
		bool bStop(false);

		for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
			if (ElemData[iCnt].bUsesAirProperties()
				&& !ElemData[iCnt].ElemMap.empty())
			{
				for (ElemMapType::const_iterator p = ElemData[iCnt].ElemMap.begin();
					p != ElemData[iCnt].ElemMap.end();
					p++)
				{
					if (dynamic_cast<AerodynamicElem *>(p->second)->NeedsAirProperties()) {
						if (bStop == 0) {
							silent_cerr("The following aerodynamic elements "
								"are defined: " << std::endl);
							bStop = true;
						}
					}
				}

				silent_cerr(ElemData[iCnt].ElemMap.size() << " "
					<< psElemNames[iCnt] << std::endl);
			}
		}

		if (bStop) {
			silent_cerr("while no air properties are defined; aborting ..."
				<< std::endl);

			throw DataManager::ErrGeneric();
		}
	}

	/* count & initialize element array */
	unsigned iNumElems = 0;
	for (int iCnt = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		iNumElems += ElemData[iCnt].ElemMap.size();
	}
	Elems.resize(iNumElems);
	for (int iCnt = 0, iElem = 0; iCnt < Elem::LASTELEMTYPE; iCnt++) {
		for (ElemMapType::const_iterator p = ElemData[iCnt].ElemMap.begin();
			p != ElemData[iCnt].ElemMap.end();
			p++, iElem++)
		{
			Elems[iElem] = p->second;
		}
	}
	
	DEBUGLCOUT(MYDEBUG_INPUT, "End of elements data" << std::endl);
} /*  End of DataManager::ReadElems()  */

Elem**
DataManager::ReadOneElem(MBDynParser& HP, unsigned int uLabel, int CurrType)
{
	Elem* pE = 0;
	Elem** ppE = 0;

	switch (KeyWords(CurrType)) {
	/* forza */
	case FORCE:
	case COUPLE:
	{
		int iForceType;
		if (KeyWords(CurrType) == FORCE) {
			iForceType = 0;
			silent_cout("Reading force " << uLabel << std::endl);

		} else /* if(KeyWords(CurrType) == COUPLE) */ {
			iForceType = 1;
			silent_cout("Reading couple " << uLabel << std::endl);
		}

		if (iNumTypes[Elem::FORCE]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": force " << uLabel
				<< " exceeds force elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::FORCE, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": force " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		pE = ReadForce(this, HP, uLabel, iForceType);
		if (pE != 0) {
			ppE = &ElemData[Elem::FORCE].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	case BODY: {
		silent_cout("Reading rigid body " << uLabel << std::endl);

		if (iNumTypes[Elem::BODY]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": rigid body " << uLabel
				<< " exceeds rigid body elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::BODY, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": rigid body " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		pE = ReadBody(this, HP, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::BODY].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* vincoli */
	case JOINT: {
		silent_cout("Reading joint " << uLabel << std::endl);

		if (iNumTypes[Elem::JOINT]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": joint " << uLabel
				<< " exceeds joint elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::JOINT, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": joint " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::JOINT].iExpectedNum
			- iNumTypes[Elem::JOINT] - 1;
		DofOwner* pDO = DofData[DofOwner::JOINT].pFirstDofOwner + i;

		pE = ReadJoint(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::JOINT].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		/* attenzione: i Joint aggiungono DofOwner e quindi devono
		 * completare la relativa struttura */
		break;
	}

	/* regolarizzazione vincoli */
	case JOINT_REGULARIZATION: {
		silent_cout("Reading joint regularization " << uLabel << std::endl);

		if (iNumTypes[Elem::JOINT_REGULARIZATION]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": " << psElemNames[Elem::JOINT_REGULARIZATION]
				<< "(" << uLabel << ")"
				" exceeds " << psElemNames[Elem::JOINT_REGULARIZATION]
				<< " elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::JOINT_REGULARIZATION, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": " << psElemNames[Elem::JOINT_REGULARIZATION]
				<< "(" << uLabel << ")"
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		pE = ReadJointRegularization(this, HP, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::JOINT_REGULARIZATION].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		/* attenzione: i Joint aggiungono DofOwner e quindi devono
		 * completare la relativa struttura */
		break;
	}

	/* trave */
	case BEAM:
	case BEAM3:		/* same as BEAM */
	case BEAM2:
	case HBEAM: {
		silent_cout("Reading beam " << uLabel << std::endl);

		if (iNumTypes[Elem::BEAM]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": beam " << uLabel
				<< " exceeds beam elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::BEAM, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": beam " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		switch (KeyWords(CurrType)) {
		case BEAM3:
		case BEAM:	/* same as BEAM3 */
			pE = ReadBeam(this, HP, uLabel);
			break;

		case BEAM2:
			pE = ReadBeam2(this, HP, uLabel);
			break;

		case HBEAM:
			pE = ReadHBeam(this, HP, uLabel);
			break;

		default:
			throw DataManager::ErrGeneric();
		}

		if (pE != 0) {
			ppE = &ElemData[Elem::BEAM].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* Elementi aerodinamici: rotori */
	case ROTOR: {
		silent_cout("Reading rotor " << uLabel << std::endl);

		if (iNumTypes[Elem::ROTOR]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": rotor " << uLabel
				<< " exceeds rotor elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::ROTOR, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": rotor " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::ROTOR].iExpectedNum
			- iNumTypes[Elem::ROTOR] - 1;
		DofOwner* pDO = DofData[DofOwner::ROTOR].pFirstDofOwner + i;

		pE = ReadRotor(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::ROTOR].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* Elementi aerodinamici: modale */
	case AEROMODAL: {
		silent_cout("Reading aero modal " << uLabel << std::endl);

		if (iNumTypes[Elem::AEROMODAL]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": aeromodal " << uLabel
				<< " exceeds aeromodal elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::AEROMODAL, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": aeromodal " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::AEROMODAL].iExpectedNum
			- iNumTypes[Elem::AEROMODAL] - 1;
		DofOwner* pDO = DofData[DofOwner::AEROMODAL].pFirstDofOwner + i;

		pE = ReadAerodynamicModal(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::AEROMODAL].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* Elementi aerodinamici: aeromodal */
	case AERODYNAMICEXTERNAL:
	case AERODYNAMICEXTERNALMODAL: {
#ifdef USE_EXTERNAL
		silent_cout("Reading external element " << uLabel << std::endl);

		if (iNumTypes[Elem::EXTERNAL]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": external element " << uLabel
				<< " exceeds external elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::EXTERNAL, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": external element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::EXTERNAL].iExpectedNum
			- iNumTypes[Elem::EXTERNAL] - 1;

		switch (KeyWords(CurrType)) {
		case AERODYNAMICEXTERNAL:
			pE = ReadAerodynamicExternal(this, HP, uLabel);
			break;

		case AERODYNAMICEXTERNALMODAL:
			pE = ReadAerodynamicExternalModal(this, HP, uLabel);
			break;

		default:
			ASSERTMSG(0, "You shouldn't have reached this point");
			break;
		}
		if (pE != 0) {
			ppE = &ElemData[Elem::EXTERNAL].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}
#else /* !USE_EXTERNAL */
		silent_cerr("You need mpi and -DUSE_AERODYNAMIC_EXTERNAL "
			<< "to use this type of elements." << std::endl);
		throw ErrGeneric();
#endif /* !USE_EXTERNAL */
		break;
	}

	case AERODYNAMICBODY:
	case AERODYNAMICBEAM:
	case AERODYNAMICBEAM3:	/* same as AERODYNAMICBEAM */
	case AERODYNAMICBEAM2:
	case AIRCRAFTINSTRUMENTS:
	case GENERICAERODYNAMICFORCE: {
		silent_cout("Reading aerodynamic element " << uLabel << std::endl);

		if (iNumTypes[Elem::AERODYNAMIC]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": aerodynamic element " << uLabel
				<< " exceeds aerodynamic elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::AERODYNAMIC, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": aerodynamic element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		switch (KeyWords(CurrType)) {
		case AERODYNAMICBODY:
			pE = ReadAerodynamicBody(this, HP, uLabel);
			break;

		case AERODYNAMICBEAM3:
		case AERODYNAMICBEAM:	/* same as BEAM3 */
			pE = ReadAerodynamicBeam(this, HP, uLabel);
			break;

		case AERODYNAMICBEAM2:
			pE = ReadAerodynamicBeam2(this, HP, uLabel);
			break;

		case AIRCRAFTINSTRUMENTS:
			pE = ReadAircraftInstruments(this, HP, uLabel);
			break;

		case GENERICAERODYNAMICFORCE:
			pE = ReadGenericAerodynamicForce(this, HP, uLabel);
			break;

		default:
			ASSERTMSG(0, "You shouldn't have reached this point");
			break;
		}
		if (pE != 0) {
			ppE = &ElemData[Elem::AERODYNAMIC].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* genel */
	case GENEL: {
		silent_cout("Reading genel " << uLabel << std::endl);

		if (iNumTypes[Elem::GENEL]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": genel " << uLabel
				<< " exceeds genel elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::GENEL, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": genel " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::GENEL].iExpectedNum
			- iNumTypes[Elem::GENEL] - 1;
		DofOwner* pDO = DofData[DofOwner::GENEL].pFirstDofOwner + i;

		pE = ReadGenel(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::GENEL].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* elementi idraulici */
	case HYDRAULIC: {
		silent_cout("Reading hydraulic element " << uLabel << std::endl);

		if (iNumTypes[Elem::HYDRAULIC]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": hydraulic element " << uLabel
				<< " exceeds hydraulic elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::HYDRAULIC, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": hydraulic element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::HYDRAULIC].iExpectedNum
			- iNumTypes[Elem::HYDRAULIC] - 1;
		DofOwner* pDO = DofData[DofOwner::HYDRAULIC].pFirstDofOwner + i;

		pE = ReadHydraulicElem(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::HYDRAULIC].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		/* attenzione: gli elementi elettrici aggiungono DofOwner
		 * e quindi devono completare la relativa struttura */

		break;
	}

	/* elementi elettrici */
	case ELECTRIC: {
		silent_cout("Reading electric element " << uLabel << std::endl);

		if (iNumTypes[Elem::ELECTRIC]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": electric element " << uLabel
				<< " exceeds electric elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::ELECTRIC, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": electric element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::ELECTRIC].iExpectedNum
			- iNumTypes[Elem::ELECTRIC] - 1;
		DofOwner* pDO = DofData[DofOwner::ELECTRIC].pFirstDofOwner + i;

		pE = ReadElectric(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::ELECTRIC].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		/* attenzione: gli elementi elettrici aggiungono DofOwner
		 * e quindi devono completare la relativa struttura */

		break;
	}

	/* elementi termici */
	case THERMAL: {
		silent_cout("Reading thermal element " << uLabel << std::endl);

		if (iNumTypes[Elem::THERMAL]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": thermal element " << uLabel
				<< " exceeds thermal elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::THERMAL, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": thermal element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::THERMAL].iExpectedNum
			- iNumTypes[Elem::THERMAL] - 1;
		DofOwner* pDO = DofData[DofOwner::THERMAL].pFirstDofOwner + i;

		pE = ReadThermal(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::THERMAL].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		/* attenzione: gli elementi termici aggiungono DofOwner
		 * e quindi devono completare la relativa struttura */

		break;
	}

	/* elementi bulk */
	case BULK: {
		silent_cout("Reading bulk element " << uLabel << std::endl);

		if (iNumTypes[Elem::BULK]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": bulk element " << uLabel
				<< " exceeds bulk elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::BULK, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": bulk element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		pE = ReadBulk(this, HP, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::BULK].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	/* elementi loadable */
	case LOADABLE: {
		silent_cout("Reading loadable element " << uLabel << std::endl);

		if (iNumTypes[Elem::LOADABLE]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": loadable element " << uLabel
				<< " exceeds loadable elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::LOADABLE, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": loadable element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
		int i = ElemData[Elem::LOADABLE].iExpectedNum
			- iNumTypes[Elem::LOADABLE] - 1;
		DofOwner* pDO = DofData[DofOwner::LOADABLE].pFirstDofOwner + i;

		pE = ReadLoadable(this, HP, pDO, uLabel);
		if (pE != 0) {
			ppE = &ElemData[Elem::LOADABLE].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}

		break;
	}

	case RTAI_OUTPUT:
#ifndef USE_RTAI
		silent_cout("need --with-rtai to allow RTAI mailboxes; "
			"using stream output..." << std::endl);
#endif /* ! USE_RTAI */
		/* fall thru... */

	case SOCKETSTREAM_OUTPUT:
	case SOCKETSTREAM_MOTION_OUTPUT: {
		silent_cout("Reading stream output element " << uLabel
			<< std::endl);

		if (iNumTypes[Elem::SOCKETSTREAM_OUTPUT]-- <= 0) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": stream output element " << uLabel
				<< " exceeds stream output elements number" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* verifica che non sia gia' definito */
		if (pFindElem(Elem::SOCKETSTREAM_OUTPUT, uLabel) != NULL) {
			DEBUGCERR("");
			silent_cerr("line " << HP.GetLineData()
				<< ": stream output element " << uLabel
				<< " already defined" << std::endl);

			throw DataManager::ErrGeneric();
		}

		/* allocazione e creazione */
#ifdef USE_RTAI
		if (mbdyn_rtai_task != 0) {
			switch (KeyWords(CurrType)) {
			case SOCKETSTREAM_OUTPUT:
			case RTAI_OUTPUT:
				silent_cerr("starting RTAI out element" << std::endl);
				pE = ReadRTAIOutElem(this, HP, uLabel);
				break;

			case SOCKETSTREAM_MOTION_OUTPUT:
			default:
				silent_cerr("line " << HP.GetLineData()
					<< ": stream motion element " << uLabel
					<< " not allowed with RTAI" << std::endl);

				throw DataManager::ErrGeneric();
			}
				
		} else
#endif /* USE_RTAI */
		{
			switch (KeyWords(CurrType)) {
			case SOCKETSTREAM_OUTPUT:
			case RTAI_OUTPUT:
#ifdef USE_SOCKET
				pedantic_cout("starting stream element" << std::endl);
				pE = ReadSocketStreamElem(this, HP, uLabel);
#else // ! USE_SOCKET
				silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
					<< "(" << uLabel << ") "
					" not allowed at line " << HP.GetLineData()
					<< " because apparently the current "
					"architecture does not support sockets"
					<< std::endl);
				throw ErrGeneric();
#endif // ! USE_SOCKET
				break;

			case SOCKETSTREAM_MOTION_OUTPUT:
#ifdef USE_SOCKET
				pedantic_cout("starting stream motion element" << std::endl);
				pE = ReadSocketStreamMotionElem(this, HP, uLabel);
#else // ! USE_SOCKET
				silent_cerr(psElemNames[Elem::SOCKETSTREAM_OUTPUT]
					<< "(" << uLabel << ") "
					" not allowed at line " << HP.GetLineData()
					<< " because apparently the current "
					"architecture does not support sockets"
					<< std::endl);
				throw ErrGeneric();
#endif // ! USE_SOCKET
				break;
			default:
				silent_cerr("You shouldn't be here: " __FILE__ << ":" << __LINE__ << std::endl);
				throw DataManager::ErrGeneric();
				break;
			}
		}
		if (pE != 0) {
			ppE = &ElemData[Elem::SOCKETSTREAM_OUTPUT].ElemMap.insert(ElemMapType::value_type(uLabel, pE)).first->second;
		}
		break;
	}

	/* In case the element type is not correct */
	default:
		silent_cerr("You shouldn't be here" << std::endl);

		throw DataManager::ErrGeneric();
	}

	/* Ritorna il puntatore al puntatore all'elemento appena costruito */
	return ppE;
}
