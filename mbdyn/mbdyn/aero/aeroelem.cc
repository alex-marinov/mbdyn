/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

/* Elementi aerodinamici */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <aeroelem.h>
#include <aerodata.h>
#include <dataman.h>
extern "C" {
#include <aerodc81.h>
#include <c81data.h>
}

#if 0
static doublereal dDA = 1.; /* diventera' parente del Delta_t (memoria) */
#endif

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
Aero_output *
Aero_output_alloc(unsigned int iNumPoints)
{
   	DEBUGCOUTFNAME("Aero_output_alloc");
	ASSERT(iNumPoints > 0);
   	Aero_output* p = NULL;
   	SAFENEWARR(p, Aero_output, iNumPoints, DMmm);  
   	return p;
}

inline void
set_alpha(Aero_output* p, Vec3& v)
{
   	DEBUGCOUTFNAME("set_alpha");
	ASSERT(p != NULL);
   	p->alpha = 180./M_PI*atan2(-v.dGet(2), v.dGet(1));
}

inline void
set_f(Aero_output* p, doublereal* pd)
{
   	DEBUGCOUTFNAME("set_f");
	ASSERT(p != NULL);
	ASSERT(pd != NULL);
   	p->f = Vec3(pd[1], pd[0], pd[5]);
}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */


/* AerodynamicBody - begin */

AerodynamicBody::AerodynamicBody(unsigned int uLabel, 
				 const StructNode* pN, Rotor* pR,
				 const Vec3& fTmp, doublereal dS,
				 const Mat3x3& RaTmp,
				 const Shape* pC, const Shape* pF, 
				 const Shape* pV, const Shape* pT,
				 integer iN, AeroData* a,
				 const DriveCaller* pDC, 
				 flag fOut)
: Elem(uLabel, Elem::AERODYNAMIC, fOut), 
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICBODY, fOut), 
InitialAssemblyElem(uLabel, Elem::AERODYNAMIC, fOut),
DriveOwner(pDC),
aerodata(a),
pNode(pN),
pRotor(pR),
fPassiveRotor(0),
f(fTmp),
dHalfSpan(dS/2.),
Ra(RaTmp),
Ra3(RaTmp.GetVec(3)), 
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
GDI(iN),
pdOuta(NULL),
pvdOuta(NULL),
F(0.),
M(0.)
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
, pOutput(NULL)
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
{
   	DEBUGCOUTFNAME("AerodynamicBody::AerodynamicBody");
   
   	ASSERT(pNode != NULL);
   	ASSERT(pNode->GetNodeType() == Node::STRUCTURAL);
   	ASSERT(aerodata != NULL);
   
#ifdef DEBUG
   	if(pRotor != NULL) {
      		ASSERT(pRotor->GetElemType() == Elem::ROTOR);
   	}
#endif /* DEBUG */
   
   	SAFENEWARR(pdOuta, doublereal, 20*iN, DMmm);
   	SAFENEWARR(pvdOuta, doublereal*, iN, DMmm);
   	for (integer i = 20*iN; i-- > 0; ) {
      		pdOuta[i] = 0.;
   	}
   	for (integer i = iN; i-- > 0; ) {
      		pvdOuta[i] = pdOuta+20*i;
   	}
   
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   	if (fToBeOutput()) {
#ifdef USE_EXCEPTIONS
      		try {
#endif /* USE_EXCEPTIONS */
	 		pOutput = Aero_output_alloc(iN);
#ifdef USE_EXCEPTIONS
      		}
      		catch (ErrMemory) {
	 		SetOutputFlag(flag(0));
	 		cerr << "Unable to alloc memory for output"
				" of AerodynamicBody " << GetLabel() << endl;
      		}
#endif /* USE_EXCEPTIONS */
   	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
}

AerodynamicBody::~AerodynamicBody(void)
{
   	DEBUGCOUTFNAME("AerodynamicBody::~AerodynamicBody");
   
   	SAFEDELETEARR(pvdOuta, DMmm);
   	SAFEDELETEARR(pdOuta, DMmm);

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   	if (pOutput != NULL) {
      		SAFEDELETEARR(pOutput, DMmm);
   	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
   
   	SAFEDELETE(aerodata, DMmm);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
void
AerodynamicBody::SetOutputFlag(flag f) 
{
   	DEBUGCOUTFNAME("AerodynamicBody::SetOutputFlag");
   	ToBeOutput::SetOutputFlag(f);
   
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   	if (fToBeOutput()) {
      		/* se non e' gia' stato allocato ... */
      		if (pOutput == NULL) {
#ifdef USE_EXCEPTIONS
	 		try {
#endif /* USE_EXCEPTIONS */
	    			pOutput = Aero_output_alloc(GDI.iGetNum());
#ifdef USE_EXCEPTIONS
	 		}
	 		catch (ErrMemory) {
	    			SetOutputFlag(flag(0));
	    			cerr << "Unable to alloc memory for output"
					" of AerodynamicBody " 
					<< GetLabel() << endl;
	 		}
#endif /* USE_EXCEPTIONS */
      		}
   	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
}

/* Scrive il contributo dell'elemento al file di restart */
ostream&
AerodynamicBody::Restart(ostream& out) const
{
   	DEBUGCOUTFNAME("AerodynamicBody::Restart");
   
   	out << "  aerodynamic body: " << GetLabel() << ", " 
		<< pNode->GetLabel();
   	if (pRotor != NULL) {
      		out << ", rotor, " << pRotor->GetLabel();
   	}
   	out << ", reference, node, ", f.Write(out, ", ")
     		<< ", " << dHalfSpan*2. 
     		<< ", reference, node, 1, ", (Ra.GetVec(1)).Write(out, ", ")
     		<< ", 2, ", (Ra.GetVec(2)).Write(out, ", ")
     		<< ", ";
   	Chord.pGetShape()->Restart(out) << ", ";
   	ForcePoint.pGetShape()->Restart(out) << ", ";
   	VelocityPoint.pGetShape()->Restart(out) << ", ";
   	Twist.pGetShape()->Restart(out) << ", " 
     		<< ", " << GDI.iGetNum() << ", control, ";
   	pGetDriveCaller()->Restart(out) << ", ";
   	aerodata->Restart(out);
   	return out << ';' << endl;
}

SubVectorHandler&
AerodynamicBody::AssRes(SubVectorHandler& WorkVec,
			doublereal /* dCoef */ ,
			const VectorHandler& /* XCurr */ ,
			const VectorHandler& /* XPrimeCurr */ )
{
   	DEBUGCOUTFNAME("AerodynamicBody::AssRes");
   	WorkVec.Resize(6);
   	WorkVec.Reset(0.);
   
   	integer iFirstIndex = pNode->iGetFirstMomentumIndex();
   	for (int iCnt = 1; iCnt <= 6; iCnt++) {
      		WorkVec.fPutRowIndex(iCnt, iFirstIndex+iCnt);
   	}
   
   	AssVec(WorkVec);
   
   	return WorkVec;
}


SubVectorHandler&
AerodynamicBody::InitialAssRes(SubVectorHandler& WorkVec,
			       const VectorHandler& /* XCurr */ )
{
   	DEBUGCOUTFNAME("AerodynamicBody::InitialAssRes");
   
   	WorkVec.Resize(6);
   	WorkVec.Reset(0.);
   
   	integer iFirstIndex = pNode->iGetFirstPositionIndex();
   	for (int iCnt = 1; iCnt <= 6; iCnt++) {
      		WorkVec.fPutRowIndex(iCnt, iFirstIndex+iCnt);
   	}
   
   	AssVec(WorkVec);
   
   	return WorkVec;
}


/* assemblaggio residuo */
void
AerodynamicBody::AssVec(SubVectorHandler& WorkVec)
{
   	DEBUGCOUTFNAME("AerodynamicBody::AssVec");
  
   	doublereal dTng[6];
   	doublereal dW[6];
  
   	/* Dati del nodo */
   	Vec3 Xn(pNode->GetXCurr());
   	Mat3x3 Rn(pNode->GetRCurr());
   	Vec3 Vn(pNode->GetVCurr());
   	Vec3 Wn(pNode->GetWCurr());

   	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
   	Mat3x3 RR(Rn*Ra);
   	Mat3x3 RRT(RR.Transpose());
      
   	/*
	 * Se l'elemento e' collegato ad un rotore,
    	 * si fa dare la velocita' di rotazione
	 */
   	doublereal dOmega = 0.;
   	if (pRotor != NULL) {
      		dOmega = pRotor->dGetOmega();
   	}

   	/* Resetta i dati */
   	F = Vec3(0.);
   	M = Vec3(0.);

   	/* Dati "permanenti" */
   	aerodata->SetAirData(dGetAirDensity(), dGetSoundSpeed());

   	doublereal** pvd = pvdOuta;
   
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   	/* per output in modalita' punti di Gauss */
   	Aero_output* pTmpOutput = pOutput;
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
   
   	/* Ciclo sui punti di Gauss */
   	PntWght PW = GDI.GetFirst();
   	do {
      		doublereal dCsi = PW.dGetPnt();
      		Vec3 Xr(Rn*(f+Ra3*(dHalfSpan*dCsi)));
      		Vec3 Vr(Vn+Wn.Cross(Xr));

      		/* Contributo di velocita' del vento */
      		Vec3 VTmp(0.);
      		if (fGetAirVelocity(VTmp, Xr)) {
	 		Vr -= VTmp;
      		}
      
      		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
      		if (pRotor != NULL) {
	 		Vr += pRotor->GetInducedVelocity(Xr+Xn);	 
      		}

      		/* Copia i dati nel vettore di lavoro dVAM */
      		doublereal dTw = Twist.dGet(dCsi)+dGet();
      		aerodata->SetSectionData(Chord.dGet(dCsi),
			       		 ForcePoint.dGet(dCsi),
					 VelocityPoint.dGet(dCsi),
					 dTw,
					 dOmega);
      
      		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
      		doublereal dCosT = cos(dTw);
		doublereal dSinT = sin(dTw);
		/* Assumo lo svergolamento positivo a cabrare */
		Mat3x3 RTw( dCosT, dSinT, 0.,
			   -dSinT, dCosT, 0.,
			    0.,    0.,    1.);
      		Mat3x3 RRloc(RR*RTw);
		Mat3x3 RRlocT(RRloc.Transpose());
		
		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
      		Vec3 Tmp(RRlocT*Vr);
      		Tmp.PutTo(dW);

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
      		if (fToBeOutput()) {
	 		ASSERT(pOutput != NULL);
	 		set_alpha(pTmpOutput, Tmp);
      		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
      
      		Tmp = RRlocT*Wn;
      		Tmp.PutTo(dW+3);
      
      		/* Funzione di calcolo delle forze aerodinamiche */
      		aerodata->GetForces(dW, dTng, *pvd);

      		/* OUTA */
      		pvd++;
      
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
      		if (fToBeOutput()) {
	 		ASSERT(pOutput != NULL);
	 		set_f(pTmpOutput, dTng);
	 		pTmpOutput++;
      		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
      
      		/* Dimensionalizza le forze */
      		doublereal dWght = PW.dGetWght();
      		Vec3 FTmp(RRloc*(Vec3(dTng)*(dHalfSpan*dWght)));
      		F += FTmp;
      		M += RRloc*(Vec3(dTng+3)*(dHalfSpan*dWght));
      		M += Xr.Cross(FTmp);
      
   	} while(GDI.fGetNext(PW));
   
   	/* Se e' definito il rotore, aggiungere il contributo alla trazione */
   	if (pRotor != NULL && !fPassiveRotor) {
      		pRotor->AddForce(F, M, Xn);
   	}
   
   	/* Sommare il termine al residuo */
   	WorkVec.Add(1, F);
   	WorkVec.Add(4, M);
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBody::Output(OutputHandler& OH) const
{
   	/* Memoria in caso di forze instazionarie */
#if 0 /* E' molto bacato */
   	if (iInst == 1) {
      		doublereal d = 0.;
      		if (pRotor != NULL) {
	 		d = dDA*pRotor->dGetOmega();
      		}
      		if (fabs(d) < 1.e-6) {
	 		d = 1.e-6;
      		}
#if 0
      		for (integer i = 0; i < GDI.iGetNum(); i++) {
	 		__FC_DECL__(coeprd)(&d, pvdOuta[i]);
      		}
#endif /* 0 */
   	} else if (iInst == 2) {
#if 0
      		for (integer i = 0; i < GDI.iGetNum(); i++) {
	 		__FC_DECL__(coeprd)(&dDA, pvdOuta[i]);
      		}
#endif /* 0 */
   	}
#endif /* 0 */

   	/* Output delle forze aerodinamiche F, M su apposito file */
   	if (fToBeOutput()) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
      		ASSERT(pOutput != NULL);
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
      		ostream& out = OH.Aerodynamic() << setw(8) << GetLabel();
      
#if AEROD_OUTPUT == AEROD_OUT_NODE
      		out << " " << setw(8) << pNode->GetLabel()
			<< " ", F.Write(out, " ") << " ", M.Write(out, " ");
#else /* AEROD_OUTPUT != AEROD_OUT_NODE */ 
      		for (int i = 0; i < GDI.iGetNum(); i++) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	 		out << ' ' << pOutput[i]->alpha
				<< ' ' << pOutput[i]->f;
#elif AEROD_OUTPUT == AEROD_OUT_STD
	 		for (int j = 1; j <= 6; j++) {
	    			out << ' ' << pvdOuta[i][j];
	 		}
#endif /* AEROD_OUTPUT == AEROD_OUT_STD */
      		}
#endif /* AEROD_OUTPUT != AEROD_OUT_NODE */
      		out << endl;      
   	}
}

/* AerodynamicBody - end */


/* Legge dati aerodinamici */

static integer
ReadUnsteadyFlag(MBDynParser& HP)
{
   	if (HP.fIsArg()) {
      		integer iInst = HP.GetInt();
      		if (iInst != 0
#if 0
	  	    && iInst != 1 && iInst != 2
#endif /* 0 */
	  	    ) {
	 		if (iInst < 0 || iInst > 2) {
	    			cerr << "illegal unsteady flag";
	 		} else {
	    			cerr << "unsteady aerodynamics"
					" are not tested yet";
	 		}
	 		cerr << " at line " << HP.GetLineData() << endl;
	 		THROW(ErrGeneric());
      		}
      		return iInst;
   	}
   	return 0;
}

static void 
ReadAeroData(DataManager* pDM,
	     MBDynParser& HP,
	     Shape** ppChord,
	     Shape** ppForce,
	     Shape** ppVelocity,
	     Shape** ppTwist,
	     integer* piNumber,
	     DriveCaller** ppDC,
	     AeroData** aerodata)
{
   	DEBUGCOUTFNAME("ReadAeroData");
      
   	/* Keywords */
   	const char* sKeyWords[] = {    
		"naca0012",
		"rae9671",
		"c81"
   	};
   
   	/* enum delle parole chiave */
   	enum KeyWords {
      		UNKNOWN = -1,     
		NACA0012 = 0,
		RAE9671,
		C81,
	
		LASTKEYWORD 
   	};
   
   	/* tabella delle parole chiave */
   	KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   	*ppChord = ReadShape(HP);
   	*ppForce = ReadShape(HP);
   	*ppVelocity = ReadShape(HP);
   	*ppTwist = ReadShape(HP);
   
   	HP.PutKeyTable(K);

   	*piNumber = HP.GetInt();
   	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Gauss points number: " << *piNumber << endl);
  
   	if (HP.IsKeyWord("control")) {      
      		/* Driver di un'eventuale controllo */
      		*ppDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
      		HP.PutKeyTable(K);
   	} else {
      		SAFENEWWITHCONSTRUCTOR(*ppDC,
			NullDriveCaller,
			NullDriveCaller(pDM->pGetDrvHdl()), DMmm);
   	}
 
   	if (HP.fIsArg()) {            
      		switch (HP.IsKeyWord()) {
       		default:
	  		cerr << "unknown profile type at line "
	    			<< HP.GetLineData()
				<< "; using default (NACA0012)" << endl;

       		case NACA0012: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				   "profile is NACA0012" << endl);
	  
	  		integer iInst = ReadUnsteadyFlag(HP);
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(1, iInst),
				DMmm);
	  		break;
       		}
	 
       		case RAE9671: {
	  		DEBUGLCOUT(MYDEBUG_INPUT,
				"profile is RAE9671" << endl);
	  
	  		integer iInst = ReadUnsteadyFlag(HP);
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				STAHRAeroData,
				STAHRAeroData(2, iInst),
				DMmm);
	  		break;
       		}
	 
	 	/*
		 * uso tabelle standard, che vengono cercate in base all'indice
		 * (sistemare) */
       		case C81: {
	  		integer iProfile = HP.GetInt();
	  		const c81_data* data = HP.GetC81Data(iProfile);

	  		DEBUGLCOUT(MYDEBUG_INPUT, 
				"profile data is from file c81 "
				<< iProfile << endl);
	  
	  		SAFENEWWITHCONSTRUCTOR(*aerodata,
				C81AeroData,
				C81AeroData(iProfile, data),
				DMmm);
			break;
       		}
      		}
   	} else {
      		cerr << "unknown profile type at line "
			<< HP.GetLineData()
			<< "; using default (NACA0012)" << endl;
      
      		integer iInst = ReadUnsteadyFlag(HP);
      		SAFENEWWITHCONSTRUCTOR(*aerodata,
			STAHRAeroData,
			STAHRAeroData(1, iInst),
			DMmm);
   	}
}


Elem *
ReadAerodynamicBody(DataManager* pDM,
		    MBDynParser& HP,
		    unsigned int uLabel)
{
   	DEBUGCOUTFNAME("ReadAerodynamicBody");
      
   	/* Nodo */
   	unsigned int uNode = (unsigned int)HP.GetInt();
   
   	DEBUGLCOUT(MYDEBUG_INPUT, "Node: " << uNode << endl);
   
   	/* verifica di esistenza del nodo */
   	StructNode* pNode = pDM->pFindStructNode(uNode);
   	if (pNode == NULL) {
      		cerr << endl
			<< " at line " << HP.GetLineData() 
			<< ": structural node " << uNode
			<< " not defined" << endl;      
		THROW(DataManager::ErrGeneric());
	}
	
	/* Eventuale rotore */
	Rotor* pRotor = NULL;
	if (HP.IsKeyWord("rotor")) {
		uNode = (unsigned int)HP.GetInt();  
		DEBUGLCOUT(MYDEBUG_INPUT,
			   "Linked to Rotor: " << uNode << endl);
		
		/*
		 * verifica di esistenza del rotore       
		 * NOTA: ovviamente il rotore deve essere definito 
		 * prima dell'elemento aerodinamico
		 */
		 Elem* p = (Elem*)(pDM->pFindElem(Elem::ROTOR, uNode));
		 if (p  == NULL) {
		 	cerr << endl
				<< " at line " << HP.GetLineData() 
				<< ": rotor " << uNode
				<< " not defined" << endl;	
				THROW(DataManager::ErrGeneric());
		}
		pRotor = (Rotor*)p->pGet();
	}
	
	ReferenceFrame RF(pNode);
	Vec3 f(HP.GetPosRel(RF));
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Offset: " << f << endl);
	
	Mat3x3 Ra(HP.GetRotRel(RF));
	
	doublereal dSpan = HP.GetReal();
	DEBUGLCOUT(MYDEBUG_INPUT, "Span: " << dSpan << endl);
	
	Shape* pChord = NULL;
	Shape* pForce = NULL;
	Shape* pVelocity = NULL;
	Shape* pTwist = NULL;
	
	integer iNumber = 0;
	DriveCaller* pDC = NULL;
	AeroData* aerodata = NULL;
	
	ReadAeroData(pDM, HP, 
		     &pChord, &pForce, &pVelocity, &pTwist,
		     &iNumber, &pDC, &aerodata);
	
	flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
	
	Elem* pEl = NULL;
	SAFENEWWITHCONSTRUCTOR(pEl, 
		AerodynamicBody,
		AerodynamicBody(uLabel, pNode, pRotor, f, dSpan, Ra,
				pChord, pForce, pVelocity, pTwist,
				iNumber, aerodata, pDC, fOut),
		DMmm);
	
	/* Se non c'e' il punto e virgola finale */
	if (HP.fIsArg()) {
		cerr << ": semicolon expected at line "
			<< HP.GetLineData() << endl;      
		THROW(DataManager::ErrGeneric());
	}
	
	return pEl;
} /* End of DataManager::ReadAerodynamicBody() */


/* AerodynamicBeam - begin */

AerodynamicBeam::AerodynamicBeam(unsigned int uLabel,
				 const Beam* pB, Rotor* pR,
				 const Vec3& fTmp1,
				 const Vec3& fTmp2,
				 const Vec3& fTmp3,
				 const Mat3x3& Ra1Tmp,
				 const Mat3x3& Ra2Tmp,
				 const Mat3x3& Ra3Tmp,
				 const Shape* pC, const Shape* pF, 
				 const Shape* pV, const Shape* pT,
				 integer iN, AeroData* a,
				 const DriveCaller* pDC,
				 flag fOut)
: Elem(uLabel, Elem::AERODYNAMIC, fOut),
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICBEAM, fOut),
InitialAssemblyElem(uLabel, Elem::AERODYNAMIC, fOut),
DriveOwner(pDC),
aerodata(a),
pBeam(pB),
pRotor(pR),
fPassiveRotor(0),
f1(fTmp1),
f2(fTmp2),
f3(fTmp3), 
Ra1(Ra1Tmp),
Ra2(Ra2Tmp),
Ra3(Ra3Tmp),
Ra1_3(Ra1Tmp.GetVec(3)),
Ra2_3(Ra2Tmp.GetVec(3)),
Ra3_3(Ra3Tmp.GetVec(3)), 
Chord(pC),
ForcePoint(pF),
VelocityPoint(pV),
Twist(pT),
GDI(iN),
pdOuta(NULL),
pvdOuta(NULL),
F1(0.),
M1(0.),
F2(0.),
M2(0.),
F3(0.),
M3(0.)
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
, pOutput(NULL)
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
{
   	DEBUGCOUTFNAME("AerodynamicBeam::AerodynamicBeam");
	
	ASSERT(pBeam != NULL);
	ASSERT(pBeam->GetElemType() == Elem::BEAM);
	
	pNode1 = pBeam->pGetNode(1); 
	pNode2 = pBeam->pGetNode(2);
	pNode3 = pBeam->pGetNode(3);
	
	ASSERT(pNode1 != NULL);
	ASSERT(pNode1->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode2 != NULL);
	ASSERT(pNode2->GetNodeType() == Node::STRUCTURAL);
	ASSERT(pNode3 != NULL);
	ASSERT(pNode3->GetNodeType() == Node::STRUCTURAL);
	
#ifdef DEBUG
	if(pRotor != NULL) {      
		ASSERT(pRotor->GetElemType() == Elem::ROTOR);
	}
#endif /* DEBUG */

	SAFENEWARR(pdOuta, doublereal, 20*3*iN, DMmm);
	SAFENEWARR(pvdOuta, doublereal*, 3*iN, DMmm);
	for (integer i = 20*3*iN; i-- > 0; ) {
		pdOuta[i] = 0.;
	}
	for (integer i = 3*iN; i-- > 0; ) {
		pvdOuta[i] = pdOuta+20*i;
	}
	
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	if (fToBeOutput()) {
#ifdef USE_EXCEPTIONS
		try {
#endif /* USE_EXCEPTIONS */
			pOutput = Aero_output_alloc(3*iN);
#ifdef USE_EXCEPTIONS
      		}
      		catch (ErrMemory) {
	 		SetOutputFlag(flag(0));
	 		cerr << "Unable to alloc memory for output"
				" of AerodynamicBeam " << GetLabel() << endl;
      		}
#endif /* USE_EXCEPTIONS */
   	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
}

AerodynamicBeam::~AerodynamicBeam(void)
{
   	DEBUGCOUTFNAME("AerodynamicBeam::~AerodynamicBeam");
	
	SAFEDELETEARR(pvdOuta, DMmm);
	SAFEDELETEARR(pdOuta, DMmm);
	
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	if (pOutput != NULL) {
		SAFEDELETEARR(pOutput, DMmm);
	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

	SAFEDELETE(aerodata, DMmm);
}

/*
 * overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione
 */
void
AerodynamicBeam::SetOutputFlag(flag f)
{
	DEBUGCOUTFNAME("AerodynamicBeam::SetOutputFlag");
	
	ToBeOutput::SetOutputFlag(f);
	
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	if (fToBeOutput()) {
		/* se non e' gia' stato allocato ... */
		if (pOutput == NULL) {
#ifdef USE_EXCEPTIONS
			try {
#endif /* USE_EXCEPTIONS */
				pOutput = Aero_output_alloc(3*GDI.iGetNum());
#ifdef USE_EXCEPTIONS
	 		}
	 		catch (ErrMemory) {
	    			SetOutputFlag(flag(0));
	    			cerr << "Unable to alloc memory for output"
					" of AerodynamicBeam "
					<< GetLabel() << endl;
	 		}
#endif /* USE_EXCEPTIONS */
      		}
   	}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
}

/* Scrive il contributo dell'elemento al file di restart */
ostream&
AerodynamicBeam::Restart(ostream& out) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Restart");
	out << "  aerodynamic beam: " << GetLabel()
		<< ", " << pBeam->GetLabel();
	if (pRotor != NULL) {
		out << ", rotor, " << pRotor->GetLabel();
	}
	out << ", reference, node, ", f1.Write(out, ", ") 
		<< ", reference, node, 1, ", (Ra1.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra1.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f2.Write(out, ", ")      
		<< ", reference, node, 1, ", (Ra2.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra2.GetVec(2)).Write(out, ", ")
		<< ", reference, node, ", f3.Write(out, ", ") 
		<< ", reference, node, 1, ", (Ra3.GetVec(1)).Write(out, ", ")
		<< ", 2, ", (Ra3.GetVec(2)).Write(out, ", ")
		<< ", ";
	Chord.pGetShape()->Restart(out) << ", ";
	ForcePoint.pGetShape()->Restart(out) << ", ";
	VelocityPoint.pGetShape()->Restart(out) << ", ";
	Twist.pGetShape()->Restart(out) << ", " 
		<< GDI.iGetNum() << ", control, ";
	pGetDriveCaller()->Restart(out) << ", ";
	aerodata->Restart(out);
	return out << ';' << endl;
}

/* assemblaggio residuo */
SubVectorHandler&
AerodynamicBeam::AssRes(SubVectorHandler& WorkVec,
			doublereal /* dCoef */ ,
			const VectorHandler& /* XCurr */ ,
			const VectorHandler& /* XPrimeCurr */ )
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssRes");
	WorkVec.Resize(18);
	WorkVec.Reset(0.);
	
	integer iNode1FirstIndex = pNode1->iGetFirstMomentumIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstMomentumIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstMomentumIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.fPutRowIndex(iCnt, iNode1FirstIndex+iCnt);
		WorkVec.fPutRowIndex(6+iCnt, iNode2FirstIndex+iCnt);
		WorkVec.fPutRowIndex(12+iCnt, iNode3FirstIndex+iCnt);
	}
	
	AssVec(WorkVec);
	
	return WorkVec;
}

/* assemblaggio iniziale residuo */
SubVectorHandler&
AerodynamicBeam::InitialAssRes(SubVectorHandler& WorkVec,
			       const VectorHandler& /* XCurr */ )
{
	DEBUGCOUTFNAME("AerodynamicBeam::InitialAssRes");
	WorkVec.Resize(18);
	WorkVec.Reset(0.);
	
	integer iNode1FirstIndex = pNode1->iGetFirstPositionIndex();
	integer iNode2FirstIndex = pNode2->iGetFirstPositionIndex();
	integer iNode3FirstIndex = pNode3->iGetFirstPositionIndex();
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.fPutRowIndex(iCnt, iNode1FirstIndex+iCnt);
		WorkVec.fPutRowIndex(6+iCnt, iNode2FirstIndex+iCnt);
		WorkVec.fPutRowIndex(12+iCnt, iNode3FirstIndex+iCnt);
	}
	
	AssVec(WorkVec);
	
	return WorkVec;
}

static const doublereal d13 = 1./sqrt(3.);

/* assemblaggio residuo */
void
AerodynamicBeam::AssVec(SubVectorHandler& WorkVec)
{
	DEBUGCOUTFNAME("AerodynamicBeam::AssVec");
	
	doublereal dTng[6];
	doublereal dW[6];
	
	/* Dati dei nodi */
	Vec3 Xn1(pNode1->GetXCurr());
	Mat3x3 Rn1(pNode1->GetRCurr());
	Vec3 Vn1(pNode1->GetVCurr());
	Vec3 Wn1(pNode1->GetWCurr());
	
	Vec3 Xn2(pNode2->GetXCurr());
	Mat3x3 Rn2(pNode2->GetRCurr());
	Vec3 Vn2(pNode2->GetVCurr());
	Vec3 Wn2(pNode2->GetWCurr());
	
	Vec3 Xn3(pNode3->GetXCurr());
	Mat3x3 Rn3(pNode3->GetRCurr());
	Vec3 Vn3(pNode3->GetVCurr());
	Vec3 Wn3(pNode3->GetWCurr());
	
	Vec3 f1Tmp(Rn1*f1);
	Vec3 f2Tmp(Rn2*f2);
	Vec3 f3Tmp(Rn3*f3);
	
	Vec3 X1Tmp(Xn1+f1Tmp);
	Vec3 X2Tmp(Xn2+f2Tmp);
	Vec3 X3Tmp(Xn3+f3Tmp);
	
	Vec3 V1Tmp(Vn1+Wn1.Cross(f1Tmp));
	Vec3 V2Tmp(Vn2+Wn2.Cross(f2Tmp));
	Vec3 V3Tmp(Vn3+Wn3.Cross(f3Tmp));
	
	/*
	 * Matrice di trasformazione dal sistema globale a quello aerodinamico
	 */
	Mat3x3 RR1(Rn1*Ra1);
	/* Mat3x3 RR1T(RR1.Transpose()); */
	
	Mat3x3 RR2(Rn2*Ra2);
	Mat3x3 RR2T(RR2.Transpose());
	
	Mat3x3 RR3(Rn3*Ra3);
	/* Mat3x3 RR3T(RR3.Transpose()); */
	
	/*
	 * Paramteri di rotazione dai nodi 1 e 3 al nodo 2 (nell'ipotesi 
	 * che tale trasformazione non dia luogo ad una singolarita')
	 */
	Vec3 g1(gparam(RR2T*RR1));
	Vec3 g3(gparam(RR2T*RR3));   
	
	/*
	 * Se l'elemento e' collegato ad un rotore,
	 * si fa dare la velocita' di rotazione
	 */
	doublereal dOmega = 0.;
	if (pRotor != NULL) {
		dOmega = pRotor->dGetOmega();
   	}

	/* Dati "permanenti" */
	aerodata->SetAirData(dGetAirDensity(), dGetSoundSpeed());
	
	/* Tratto relativo al primo nodo */
	
	/* Resetta i dati */
	F1 = Vec3(0.);
	M1 = Vec3(0.);
	
	doublereal dsi = -1.;
	doublereal dsf = -d13; /* -1./sqrt(3.) */

	doublereal dsm = (dsf+dsi)/2.;
	doublereal dsdCsi = (dsf-dsi)/2.;
	
	/* OUTA */
	doublereal** pvd = pvdOuta;
	
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	/* per output */
	Aero_output* pTmpOutput = pOutput;
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

	/* Ciclo sui punti di Gauss */
	PntWght PW = GDI.GetFirst();
	do {
		doublereal dCsi = PW.dGetPnt();
		doublereal ds = dsm+dsdCsi*dCsi;
		doublereal dXds = DxDcsi3N(ds, Xn1, Xn2, Xn3);
		
		doublereal dN1 = ShapeFunc3N(ds, 1);
		doublereal dN2 = ShapeFunc3N(ds, 2);
		doublereal dN3 = ShapeFunc3N(ds, 3);
		
		Vec3 Xr(X1Tmp*dN1+X2Tmp*dN2+X3Tmp*dN3);
		Vec3 Vr(V1Tmp*dN1+V2Tmp*dN2+V3Tmp*dN3);
		Vec3 Wr(Wn1*dN1+Wn2*dN2+Wn3*dN3);
		
		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xr)) {
			Vr -= VTmp;
		}
		
		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pRotor != NULL) {
			Vr += pRotor->GetInducedVelocity(Xr);
		}
      
      		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(ds);
		dTw += dGet(); /* Contributo dell'eventuale sup. mobile */
		aerodata->SetSectionData(Chord.dGet(ds), 
					 ForcePoint.dGet(ds),
					 VelocityPoint.dGet(ds),
					 dTw,
					 dOmega);
		
		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		doublereal dCosT = cos(dTw);
		doublereal dSinT = sin(dTw);
		/* Assumo lo svergolamento positivo a cabrare */
		Mat3x3 RTw( dCosT, dSinT, 0.,
			   -dSinT, dCosT, 0.,
			    0.,    0.,    1.);
		/*
		 * Allo stesso tempo interpola le g e aggiunge lo svergolamento
		 */
		Mat3x3 RRloc(RR2*Mat3x3(MatR, g1*dN1+g3*dN3)*RTw);
		Mat3x3 RRlocT(RRloc.Transpose());
		
		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		Vec3 Tmp(RRlocT*Vr);
		Tmp.PutTo(dW);
		
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_alpha(pTmpOutput, Tmp);
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

		Tmp = RRlocT*Wr;
		Tmp.PutTo(dW+3);
		
		/* Funzione di calcolo delle forze aerodinamiche */
		aerodata->GetForces(dW, dTng, *pvd);
		
		/* OUTA */
		pvd++;
		
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_f(pTmpOutput, dTng);
			pTmpOutput++;
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
	    
		/* Dimensionalizza le forze */
		doublereal dWght = dXds*dsdCsi*PW.dGetWght();
		Vec3 FTmp(RRloc*(Vec3(dTng)*dWght));
		F1 += FTmp;
		M1 += RRloc*(Vec3(dTng+3)*dWght);
		M1 += (Xr-Xn1).Cross(FTmp);
		
	} while (GDI.fGetNext(PW));
	
	/* Se e' definito il rotore, aggiungere il contributo alla trazione */
	if (pRotor != NULL && !fPassiveRotor) {
		pRotor->AddForce(F1, M1, Xn1);
	}
	
	/* Somma il termine al residuo */
	WorkVec.Add(1, F1);
	WorkVec.Add(4, M1);
	
	
	/* Tratto relativo al secondo nodo */
	
	/* Resetta i dati */
	F2 = Vec3(0.);
	M2 = Vec3(0.);
	
	dsi = -d13; /* -1./sqrt(3.) */
	dsf = d13;  /* 1./sqrt(3.) */

	dsm = (dsf+dsi)/2.;
	dsdCsi = (dsf-dsi)/2.;
	
	/* Ciclo sui punti di Gauss */
	PW = GDI.GetFirst();
	do {
		doublereal dCsi = PW.dGetPnt();
		doublereal ds = dsm+dsdCsi*dCsi;
		doublereal dXds = DxDcsi3N(ds, Xn1, Xn2, Xn3);
		
		doublereal dN1 = ShapeFunc3N(ds, 1);
		doublereal dN2 = ShapeFunc3N(ds, 2);
		doublereal dN3 = ShapeFunc3N(ds, 3);
		
		Vec3 Xr(X1Tmp*dN1+X2Tmp*dN2+X3Tmp*dN3);
		Vec3 Vr(V1Tmp*dN1+V2Tmp*dN2+V3Tmp*dN3);
		Vec3 Wr(Wn1*dN1+Wn2*dN2+Wn3*dN3);
		
		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xr)) {
			Vr -= VTmp;
		}
		
		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		if (pRotor != NULL) {
			Vr += pRotor->GetInducedVelocity(Xr);
		}
		
		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(ds);
		dTw += dGet(); /* Contributo dell'eventuale sup. mobile */
		aerodata->SetSectionData(Chord.dGet(ds), 
					 ForcePoint.dGet(ds),
					 VelocityPoint.dGet(ds),
					 dTw,
					 dOmega);
		
		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		doublereal dCosT = cos(dTw);
		doublereal dSinT = sin(dTw);
		/* Assumo lo svergolamento positivo a cabrare */
		Mat3x3 RTw( dCosT, dSinT, 0.,
			   -dSinT, dCosT, 0.,
			    0.,    0.,    1.);
		/*
		 * Allo stesso tempo interpola le g e aggiunge lo svergolamento
		 */
		Mat3x3 RRloc(RR2*Mat3x3(MatR, g1*dN1+g3*dN3)*RTw);
		Mat3x3 RRlocT(RRloc.Transpose());
		
		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		Vec3 Tmp(RRlocT*Vr);
		Tmp.PutTo(dW);

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_alpha(pTmpOutput, Tmp);
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

		Tmp = RRlocT*Wr;
		Tmp.PutTo(dW+3);
		
		/* Funzione di calcolo delle forze aerodinamiche */
		aerodata->GetForces(dW, dTng, *pvd);
		
		/* OUTA */
		pvd++;

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_f(pTmpOutput, dTng);
			pTmpOutput++;
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

		/* Dimensionalizza le forze */
		doublereal dWght = dXds*dsdCsi*PW.dGetWght();
		Vec3 FTmp(RRloc*(Vec3(dTng)*dWght));
		F2 += FTmp;
		M2 += RRloc*(Vec3(dTng+3)*dWght);
		M2 += (Xr-Xn2).Cross(FTmp);
		
	} while (GDI.fGetNext(PW));
	
	/* Se e' definito il rotore, aggiungere il contributo alla trazione */
	if (pRotor != NULL && !fPassiveRotor) {
		pRotor->AddForce(F2, M2, Xn2);
	}
	
	/* Somma il termine al residuo */
	WorkVec.Add(7, F2);
	WorkVec.Add(10, M2);
		
	/* Tratto relativo al terzo nodo */
	
	/* Resetta i dati */
	F3 = Vec3(0.);
	M3 = Vec3(0.);
	
	dsi = d13; /* 1./sqrt(3.) */
	dsf = 1.;

	dsm = (dsf+dsi)/2.;
	dsdCsi = (dsf-dsi)/2.;
	
	/* Ciclo sui punti di Gauss */
	PW = GDI.GetFirst();
	do {
		doublereal dCsi = PW.dGetPnt();
		doublereal ds = dsm+dsdCsi*dCsi;
		doublereal dXds = DxDcsi3N(ds, Xn1, Xn2, Xn3);
		
		doublereal dN1 = ShapeFunc3N(ds, 1);
		doublereal dN2 = ShapeFunc3N(ds, 2);
		doublereal dN3 = ShapeFunc3N(ds, 3);
		
		Vec3 Xr(X1Tmp*dN1+X2Tmp*dN2+X3Tmp*dN3);
		Vec3 Vr(V1Tmp*dN1+V2Tmp*dN2+V3Tmp*dN3);
		Vec3 Wr(Wn1*dN1+Wn2*dN2+Wn3*dN3);
		
		/* Contributo di velocita' del vento */
		Vec3 VTmp(0.);
		if (fGetAirVelocity(VTmp, Xr)) {
			Vr -= VTmp;
		}
		
		/*
		 * Se l'elemento e' collegato ad un rotore,
		 * aggiunge alla velocita' la velocita' indotta
		 */
		 if (pRotor != NULL) {
		 	Vr += pRotor->GetInducedVelocity(Xr);
		}
		
		/* Copia i dati nel vettore di lavoro dVAM */
		doublereal dTw = Twist.dGet(ds);
		dTw += dGet(); /* Contributo dell'eventuale sup. mobile */
		aerodata->SetSectionData(Chord.dGet(ds), 
					 ForcePoint.dGet(ds),
					 VelocityPoint.dGet(ds),
					 dTw,
					 dOmega);
		
		/*
		 * Lo svergolamento non viene piu' trattato in aerod2_; quindi
		 * lo uso per correggere la matrice di rotazione
		 * dal sistema aerodinamico a quello globale
		 */
		doublereal dCosT = cos(dTw);
		doublereal dSinT = sin(dTw);
		/* Assumo lo svergolamento positivo a cabrare */
		Mat3x3 RTw( dCosT, dSinT, 0.,
			   -dSinT, dCosT, 0.,
			    0.,    0.,    1.);
		/*
		 * Allo stesso tempo interpola le g e aggiunge lo svergolamento
		 */
		Mat3x3 RRloc(RR2*Mat3x3(MatR, g1*dN1+g3*dN3)*RTw);
		Mat3x3 RRlocT(RRloc.Transpose());
		
		/*
		 * Ruota velocita' e velocita' angolare nel sistema
		 * aerodinamico e li copia nel vettore di lavoro dW
		 */
		Vec3 Tmp(RRlocT*Vr);
		Tmp.PutTo(dW);
		
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_alpha(pTmpOutput, Tmp);
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

		Tmp = RRlocT*Wr;
		Tmp.PutTo(dW+3);
		
		/* Funzione di calcolo delle forze aerodinamiche */
		aerodata->GetForces(dW, dTng, *pvd);
		
		/* OUTA */
		pvd++;
		
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		if (fToBeOutput()) {
			ASSERT(pOutput != NULL);
			set_f(pTmpOutput, dTng);
			pTmpOutput++;
		}
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */

		/* Dimensionalizza le forze */
		doublereal dWght = dXds*dsdCsi*PW.dGetWght();
		Vec3 FTmp(RRloc*(Vec3(dTng)*dWght));
		F3 += FTmp;
		M3 += RRloc*(Vec3(dTng+3)*dWght);
		M3 += (Xr-Xn3).Cross(FTmp);
		
	} while (GDI.fGetNext(PW));
	
	/* Se e' definito il rotore, aggiungere il contributo alla trazione */
	if (pRotor != NULL && !fPassiveRotor) {
		pRotor->AddForce(F3, M3, Xn3);
	}
	
	/* Somma il termine al residuo */
	WorkVec.Add(13, F3);
	WorkVec.Add(16, M3);
}

/*
 * output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output
 */
void
AerodynamicBeam::Output(OutputHandler& OH ) const
{
	DEBUGCOUTFNAME("AerodynamicBeam::Output");
	
#if 0 /* Ha seri bugs */
	/* Memoria in caso di forze instazionarie */
	if (iInst == 1) {
		doublereal d = 0.;
		if (pRotor != NULL) {
			d = dDA*pRotor->dGetOmega();
		}
		if (fabs(d) < 1.e-6) {
			d = 1.e-6;
		}
#if 0
		for (integer i = 0; i < 3*GDI.iGetNum(); i++) {
			__FC_DECL__(coeprd)(&d, pvdOuta[i]);
		}
#endif /* 0 */
	} else if (iInst == 2) {
#if 0
		for (integer i = 0; i < 3*GDI.iGetNum(); i++) {
			__FC_DECL__(coeprd)(&dDA, pvdOuta[i]);
		}
#endif /* 0 */
	}
#endif 

	if (fToBeOutput()) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
		ASSERT(pOutput != NULL);
#endif /* AEROD_OUTPUT == AEROD_OUT_PGAUSS */
		ostream& out = OH.Aerodynamic() << setw(8) << GetLabel();
		
#if AEROD_OUTPUT == AEROD_OUT_NODE
		out << " " << setw(8) << pBeam->GetLabel();
		out << " ", F1.Write(out, " ") << " ", M1.Write(out, " ");
		out << " ", F2.Write(out, " ") << " ", M2.Write(out, " ");
		out << " ", F3.Write(out, " ") << " ", M3.Write(out, " ");
#else /* AEROD_OUTPUT != AEROD_OUT_NODE */
		for (int i = 0; i < 3*GDI.iGetNum(); i++) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
			out << ' ' << pOutput[i]->alpha
				<< ' ' << pOutput[i]->f;
#elif AEROD_OUTPUT == AEROD_OUT_STD
			for (int j = 1; j <= 6; j++) {
				out << ' ' << pvdOuta[i][j];
			}
#endif /* AEROD_OUTPUT == AEROD_OUT_STD */
		}
#endif /* AEROD_OUTPUT != AEROD_OUT_NODE */
		out << endl;
	}
}

/* AerodynamicBeam - end */


/* Legge un elemento aerodinamico di trave */

Elem *
ReadAerodynamicBeam(DataManager* pDM,
		    MBDynParser& HP,
		    unsigned int uLabel)
{
	DEBUGCOUTFNAME("ReadAerodynamicBeam");
	
	/* Trave */
	unsigned int uBeam = (unsigned int)HP.GetInt();
	
	DEBUGLCOUT(MYDEBUG_INPUT, "Linked to beam: " << uBeam << endl);
	
	/* verifica di esistenza della trave */
	Beam* pBeam;
	Elem* p = (Elem*)(pDM->pFindElem(Elem::BEAM, uBeam));
	if (p == NULL) {
		cerr << " at line " << HP.GetLineData() 
			<< ": beam " << uBeam
			<< " not defined" << endl;      
		THROW(DataManager::ErrGeneric());
	}
	pBeam = (Beam*)p->pGet();
	
	/* Eventuale rotore */
	Rotor* pRotor = NULL;
	if (HP.IsKeyWord("rotor")) {
		unsigned int uRotor = (unsigned int)HP.GetInt();
		DEBUGLCOUT(MYDEBUG_INPUT,
			"Linked to Rotor: " << uRotor << endl);
		
		/*
		 * verifica di esistenza del rotore
		 * NOTA: ovviamente il rotore deve essere definito 
		 * prima dell'elemento aerodinamico
		 */
		Elem* p = (Elem*)(pDM->pFindElem(Elem::ROTOR, uRotor));
		if (p == NULL) {
			cerr << " at line " << HP.GetLineData() 
				<< ": rotor " << uRotor
				<< " not defined" << endl;	 
			THROW(DataManager::ErrGeneric());
		}
		pRotor = (Rotor*)p->pGet();
	}
	
	/* Nodo 1: */
	
	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode1 = pBeam->pGetNode(1);
	
	ReferenceFrame RF(pNode1);
	Vec3 f1(HP.GetPosRel(RF));      
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 1 offset: " << f1 << endl);
	
	Mat3x3 Ra1(HP.GetRotRel(RF));
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 1 rotation matrix: " << endl << Ra1 << endl);

	/* Nodo 2: */
	
	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode2 = pBeam->pGetNode(2);
	
	RF = ReferenceFrame(pNode2);
	Vec3 f2(HP.GetPosRel(RF));    
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 2 offset: " << f2 << endl);
	
	Mat3x3 Ra2(HP.GetRotRel(RF));   
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 2 rotation matrix: " << endl << Ra2 << endl);

	/* Nodo 3: */
	
	/* Offset del corpo aerodinamico rispetto al nodo */
	const StructNode* pNode3 = pBeam->pGetNode(3);
	
	RF = ReferenceFrame(pNode3);
	Vec3 f3(HP.GetPosRel(RF));   
	DEBUGLCOUT(MYDEBUG_INPUT, "Node 3 offset: " << f3 << endl);
	
	Mat3x3 Ra3(HP.GetRotRel(RF));   
	DEBUGLCOUT(MYDEBUG_INPUT,
		   "Node 3 rotation matrix: " << endl << Ra3 << endl);

	Shape* pChord = NULL;
	Shape* pForce = NULL;
	Shape* pVelocity = NULL;
	Shape* pTwist = NULL;
	
	integer iNumber = 0;
	DriveCaller* pDC = NULL;
	AeroData* aerodata = NULL;
	
	ReadAeroData(pDM, HP, 
		     &pChord, &pForce, &pVelocity, &pTwist,
		     &iNumber, &pDC, &aerodata);
	
	flag fOut = pDM->fReadOutput(HP, Elem::BEAM);
	
	Elem* pEl = NULL;
	
	SAFENEWWITHCONSTRUCTOR(pEl,
		AerodynamicBeam,
		AerodynamicBeam(uLabel, pBeam, pRotor, 
				f1, f2, f3, Ra1, Ra2, Ra3,
				pChord, pForce, pVelocity, pTwist,
				iNumber, aerodata, pDC, fOut),
		DMmm);
	
	/* Se non c'e' il punto e virgola finale */
	if (HP.fIsArg()) {
		cerr << ": semicolon expected at line "
			<< HP.GetLineData() << endl;      
		THROW(DataManager::ErrGeneric());
	}
	
	return pEl;
} /* End of DataManager::ReadAerodynamicBeam() */

