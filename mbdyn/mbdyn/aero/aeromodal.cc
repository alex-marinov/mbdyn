/* Elemento aerodinamico modale */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <aeromodal.h>
#include <dataman.h>


/* Workspace statici usati nel calcolo delle forze */
static doublereal dTng[6];
static doublereal dW[6];
static doublereal dVAM[6];
static doublereal dDA = 1.; /* diventera' parente del Delta_t (memoria) */


#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
Aero_output* Aero_output_alloc(unsigned int iNumPoints)
{
   // DEBUGCOUTFNAME("Aero_output_alloc");
   Aero_output* p = NULL;
   SAFENEWARR(p, Aero_output, iNumPoints, EMmm);  
   return p;
}

inline void set_alpha(Aero_output* p, Vec3& v)
{
   // DEBUGCOUTFNAME("set_alpha");
   p->alpha = 180./M_PI*atan2(-v.dGet(2), v.dGet(1));
}


inline void set_f(Aero_output* p, doublereal* pd)
{
   // DEBUGCOUTFNAME("set_f");
   p->f = Vec3(*(pd+1), *(pd), *(pd+5));
}
#endif /* AEROD_OUTPUT */


/* AerodynamicModal - begin */

AerodynamicModal::AerodynamicModal(unsigned int uLabel, 
				 const StructNode* pN, const Modal* pMJ, Rotor* pR,
				 const Mat3x3& RaTmp,
				 const Shape* pC, const Shape* pF, 
				 const Shape* pV, const Shape* pT,
				 integer iI, integer iN, integer iP,
                                 integer iM, integer NFN, 
				 integer iAP,
                                 Mat3xN* pModeShapest, Mat3xN* pModeShapesr,
                                 FullMatrixHandler* pH,
                                 FullMatrixHandler* pGTKG,
                                 double dL,
                                 Mat3xN *pFNP,
                                 //integer Zaxis,
				 const DriveCaller* pDC, 
				 flag fOut)
: Elem(uLabel, Elem::AERODYNAMIC, fOut), 
AerodynamicElem(uLabel, AerodynamicElem::AERODYNAMICMODAL, fOut), 
InitialAssemblyElem(uLabel, Elem::AERODYNAMIC, fOut),
DriveOwner(pDC),
pModalNode(pN), pModalJoint(pMJ), pRotor(pR),NAeroElems(iAP),
NModes(iM),NFemNodes(NFN),
pPHIt(pModeShapest), pPHIr(pModeShapesr),
Ra(RaTmp), Ra3(RaTmp.GetVec(3)), 
Chord(pC), ForcePoint(pF), VelocityPoint(pV), Twist(pT),
iInst(iI), iProfile(iP), GDI(iN), pdOuta(NULL), pvdOuta(NULL),
F(0.), M(0.),ppR1tot(NULL),
a(NModes, 0.), aPrime(NModes, 0.),
pH(pH), pGTKG(pGTKG), dL(dL),
pFemNodesPosition(pFNP)
//Zaxis(Zaxis)
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
, pOutput(NULL)
#endif /* AEROD_OUTPUT */
{
   // DEBUGCOUTFNAME("AerodynamicModal::AerodynamicModal");
   
   ASSERT(pModalNode != NULL);
   ASSERT(pModalNode->GetNodeType() == Node::STRUCTURAL);
   
#ifdef DEBUG
   if(pRotor != NULL) {      
      ASSERT(pRotor->GetElemType() == Elem::ROTOR);
   }
#endif
   
   SAFENEWARR(pdOuta, doublereal, 20*iN*NAeroElems, EMmm);
   SAFENEWARR(pvdOuta, doublereal*, iN*NAeroElems, EMmm);
   for (integer i = 20*iN*NAeroElems; i-- > 0; ) {
      pdOuta[i] = 0.;
   }
   for (integer i = iN*NAeroElems; i-- > 0; ) {
      pvdOuta[i] = pdOuta+20*i;
   }
   SAFENEWARR(ppR1tot, Mat3x3*, NAeroElems, EMmm);
   for ( int i = 0; i < NAeroElems; i++ ) {
    ppR1tot[i] = NULL;
    SAFENEW(ppR1tot[i], Mat3x3, EMmm);
  }
  
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   if (fToBeOutput()) {
#ifdef USE_EXCEPTIONS
      try {
#endif // USE_EXCEPTIONS	 
	 pOutput = Aero_output_alloc(iN);
#ifdef USE_EXCEPTIONS
      }
      catch (ErrMemory) {
	 SetOutputFlag(flag(0));
	 cerr << "Unable to alloc memory for output of AerodynamicModal " 
	   << GetLabel() << endl;
      }      
#endif // USE_EXCEPTIONS	 
   }
#endif /* AEROD_OUTPUT */

}


AerodynamicModal::~AerodynamicModal(void)
{
   // DEBUGCOUTFNAME("AerodynamicModal::~AerodynamicModal");
   
   SAFEDELETEARR(pvdOuta, EMmm);   
   SAFEDELETEARR(pdOuta, EMmm);
   if (ppR1tot != NULL) {
    for ( int i = 0; i < NAeroElems; i++ ) {
      if (ppR1tot[i] != NULL) {
	SAFEDELETE(ppR1tot[i], EMmm);
      }
    }
    SAFEDELETEARR(ppR1tot, EMmm);
  } 
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   if (pOutput != NULL) {
      SAFEDELETEARR(pOutput, EMmm);
   }
#endif /* AEROD_OUTPUT */
}


/* overload della funzione di ToBeOutput();
 * serve per allocare il vettore dei dati di output se il flag
 * viene settato dopo la costruzione */
void AerodynamicModal::SetOutputFlag(flag f) 
{
   // DEBUGCOUTFNAME("AerodynamicModal::SetOutputFlag");   
   ToBeOutput::SetOutputFlag(f);
   
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   if (fToBeOutput()) {
      /* se non e' gia' stato allocato ... */
      if (pOutput == NULL) {
#ifdef USE_EXCEPTIONS
	 try {
#endif // USE_EXCEPTIONS	 
	    pOutput = Aero_output_alloc(GDI.iGetNum());
#ifdef USE_EXCEPTIONS
	 }
	 catch (ErrMemory) {
	    SetOutputFlag(flag(0));
	    cerr << "Unable to alloc memory for output of AerodynamicModal " 
	      << GetLabel() << endl;
	 }      
#endif // USE_EXCEPTIONS	 
      }	 
   }
#endif /* AEROD_OUTPUT */
}


/* Scrive il contributo dell'elemento al file di restart */
ostream& AerodynamicModal::Restart(ostream& out) const
{

#if 0
   // DEBUGCOUTFNAME("AerodynamicModal::Restart");
   
   out
     << "  aerodynamic modal: " << GetLabel() << ", "
     << pModalNode->GetLabel();
   if(pRotor != NULL) {
      out << ", rotor, " << pRotor->GetLabel();
   }
   out
     << ", reference, node, 1, ", (Ra.GetVec(1)).Write(out, ", ")
     << ", 2, ", (Ra.GetVec(2)).Write(out, ", ")
     << ", ";
   Chord.pGetShape()->Restart(out) << ", ";
   ForcePoint.pGetShape()->Restart(out) << ", ";
   VelocityPoint.pGetShape()->Restart(out) << ", ";
   Twist.pGetShape()->Restart(out) << ", " 
     << iInst << ", " << GDI.iGetNum() << ", control, ";
   pGetDriveCaller()->Restart(out) << ", ";
   if(iProfile == 1) {
      out << "NACA0012";
   } else if(iProfile == 2) {
      out << "RAE9671";
   }
   return out << ';' << endl;
#else
   return out << "  /* aerodynamic modal: not implemented yet */" << endl;
#endif
}


SubVectorHandler& AerodynamicModal::AssRes(SubVectorHandler& WorkVec,
					  doublereal /* dCoef */ ,
					  const VectorHandler& XCurr,
					  const VectorHandler&  XPrimeCurr )
{
   DEBUGCOUTFNAME("AerodynamicModal::AssRes");
   WorkVec.Resize(6+NModes);
   WorkVec.Reset(0.);
   
   integer iFirstIndex = pModalNode->iGetFirstIndex();
   integer iModalIndex = pModalJoint->iGetModalIndex();

   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, 6+iFirstIndex+iCnt);
   }   
   
    for(int iCnt = 1; iCnt <= NModes; iCnt++) {
      WorkVec.fPutRowIndex(6+iCnt, iModalIndex+NModes+iCnt);
    } 
   
   /* Recupera i vettori {a} e {aP}(deformate modali) */  
   for(integer  iCnt=1; iCnt<=NModes; iCnt++) {
     a.Put(iCnt, XCurr.dGetCoef(iModalIndex+iCnt));
     aPrime.Put(iCnt, XPrimeCurr.dGetCoef(iModalIndex+iCnt)); 
   }
  
   AssVec(WorkVec);
   
   return WorkVec;
}


SubVectorHandler& AerodynamicModal::InitialAssRes(SubVectorHandler& WorkVec,
						 const VectorHandler&  XCurr )
{
   DEBUGCOUTFNAME("AerodynamicModal::InitialAssRes");
   
   WorkVec.Resize(6+NModes);
   WorkVec.Reset(0.);
   
   integer iFirstIndex = pModalNode->iGetFirstIndex();
   integer iModalIndex = pModalJoint->iGetModalIndex();

   /* modificare */

   for(int iCnt = 1; iCnt <= 6; iCnt++) {
      WorkVec.fPutRowIndex(iCnt, iFirstIndex+iCnt);
   }   

   for(int iCnt = 1; iCnt <= NModes; iCnt++) {
      WorkVec.fPutRowIndex(6+iCnt, iModalIndex+iCnt);
   }
   
   /* Recupera i vettori {a} e {aP}(deformate modali) */  
   for(integer  iCnt=1; iCnt<=NModes; iCnt++) {
     a.Put(iCnt, XCurr.dGetCoef(iModalIndex+iCnt));
     aPrime.Put(iCnt, XCurr.dGetCoef(iModalIndex+NModes+iCnt)); 
   }
    
   AssVec(WorkVec);
   
   return WorkVec;
}


/* assemblaggio residuo */
void AerodynamicModal::AssVec(SubVectorHandler& WorkVec)
{
   DEBUGCOUTFNAME("AerodynamicModal::AssVec");

   /* Dati del nodo rigido */
   Vec3 X0(pModalNode->GetXCurr());
   Mat3x3 Rn(pModalNode->GetRCurr());
   Vec3 V0(pModalNode->GetVCurr());
   Vec3 W0(pModalNode->GetWCurr());

   Mat3x3 RnT = Rn.Transpose();
   Mat3x3 RR(Rn*Ra);
   Mat3x3 RRT(RR.Transpose());
  
   Mat3x3 RnTot(0.), RnTotT(0.), RRTot(0.), RRTotT(0.);

   /* Se l'elemento e' collegato ad un rotore,
    * si fa dare la velocita' di rotazione */
   
   doublereal dOmega = 0.;
   if(pRotor != NULL) {
      dOmega = pRotor->dGetOmega();
   }   

   /* Resetta i dati */

   Vec3 FTot(0.);
   Vec3 MTot(0.);

   /*
    * Dati "permanenti" (uso la posizione del corpo perche'
    * non dovrebbero cambiare "molto")
    */
   dVAM[0] = dGetAirDensity(X0);
   dVAM[1] = dGetSoundSpeed(X0);

   doublereal** pvd = pvdOuta;
   
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   /* per output */
   Aero_output* pTmpOutput = pOutput;      
#endif /* AEROD_OUTPUT */

   Mat3xN PHIt(NModes,0.), PHIr(NModes,0.);
   MatNx3 PHItT(NModes,0.),PHIrT(NModes,0.);

   VecN AeroDef(NAeroElems*2, 0.);    // deformate dei nodi aerodinamici
   VecN AeroVel(NAeroElems*2, 0.);    // vel. di def. dei nodi aerodinamici
   VecN Fa(NAeroElems*2,0.);

    /* moltiplica la matrice H (naero x nmodes) per il vettore delle deformate
       modali: recupero il vettore delle deformazioni e vel. di deformazioni aerodinamiche */ 

   for(int iCnt = 1; iCnt <= 2*NAeroElems; iCnt ++) {
     double temp1 = 0., temp2 = 0.;
     for(int jCnt = 1; jCnt <= NModes; jCnt ++) {
       temp1 += pH->dGetCoef(iCnt, jCnt)*a.dGet(jCnt);
       temp2 += pH->dGetCoef(iCnt, jCnt)*aPrime.dGet(jCnt);
     }
     AeroDef.Put(iCnt, temp1);
     AeroVel.Put(iCnt, temp2); 
   }

    /* recupero l'offset rigido dei nodi aerodinamici: {xa} = [G]*{xs} */
    
    VecN AeroPointsOffset(2*NAeroElems, 0.);
    for(int iCnt = 1; iCnt <= 2*NAeroElems; iCnt ++) {
      double temp = 0.;
      for(int jCnt = 1; jCnt <= NFemNodes; jCnt ++) {
         temp += pGTKG->dGetCoef(iCnt, jCnt)*pFemNodesPosition->dGet(1,jCnt);
      }
     AeroPointsOffset.Put(iCnt, temp);
    } 

    for(int iAero = 1; iAero <= NAeroElems; iAero ++) { // ciclo esteso agli elementi aerodiamici

      /* Resetta i dati (F e M sono le forze e i momenti per ogni elemento aerodinamico)*/
      F = Vec3(0.);
      M = Vec3(0.);

      Vec3 Xrig(AeroPointsOffset.dGet((iAero-1)*2+1),0.,0.);
      Vec3 gdef(AeroDef.dGet((iAero-1)*2+2),0.,0.);
      Vec3 Wdef(AeroVel.dGet((iAero-1)*2+2),0.,0.);
      Vec3 Xdef(0.,0.,AeroDef.dGet((iAero-1)*2+1));
      Vec3 Vdef(0.,0.,AeroVel.dGet((iAero-1)*2+1));
        
      /* Ciclo sui punti di Gauss */
      PntWght PW = GDI.GetFirst();
      do { 
	 doublereal dCsi = PW.dGetPnt();
         Vec3 Xr(Rn*(Xrig+Ra3*(dL/2.*dCsi)+Xdef));
	 Vec3 Vr(V0+W0.Cross(Xr)+Rn*Vdef);
	 Vec3 Wr(W0+Rn*Wdef);

         DEBUGCOUT("X0 :"<<X0<<endl);
         DEBUGCOUT("V0 :"<<V0<<endl);
         DEBUGCOUT("Xr :"<<Xr<<endl);
         DEBUGCOUT("Vr :"<<Vr<<endl);
         DEBUGCOUT("W0 :"<<W0<<endl);
         DEBUGCOUT("Rn :"<<Rn<<endl);
     
	 
	 /* Contributo di velocita' del vento */
	 /* per adesso e' costante, poi variera' */
	 Vec3 VTmp(0.);
	 if(fGetAirVelocity(VTmp, Xr)) {
	    Vr -= VTmp;
	 }
	 
	 
	 /* Se l'elemento e' collegato ad un rotore,
	  * aggiunge alla velocita' la velocita' indotta */
	 if(pRotor != NULL) {
	    Vr += pRotor->GetInducedVelocity(Xr+X0);  	 
	 }        
	 
	 /* Copia i dati nel vettore di lavoro dVAM */
	 dVAM[2] = Chord.dGet(dCsi);
	 dVAM[3] = ForcePoint.dGet(dCsi);
	 dVAM[4] = VelocityPoint.dGet(dCsi);
	 doublereal dTw = Twist.dGet(dCsi);
	 dTw += dGet(); /* Contributo dell'eventuale superficie mobile */
	 dVAM[5] = dTw;
	 
	 /* Matrici di trasformazione da sdr globale a nodale */
       
	 RnTot  = Rn + Rn*Mat3x3(gdef);
         RnTotT = RnTot.Transpose();
         RRTot  = RnTot*Ra;
         RRTotT = RRTot.Transpose();

	 /* Lo svergolamento non viene piu' trattato in aerod2_; quindi
	  * lo uso per correggere la matrice di rotazione
	  * dal sistema aerodinamico a quello globale */
	 doublereal dCosT = cos(dTw);
	 doublereal dSinT = sin(dTw);
	 /* Assumo lo svergolamento positivo a cabrare */
	 Mat3x3 RTw( dCosT, dSinT, 0.,
		    -dSinT, dCosT, 0.,
		    0.,    0.,    1.);
	 /* aggiungo lo svergolamento */
	 Mat3x3 RRloc(RR*RTw);       
	 Mat3x3 RRlocT(RRloc.Transpose());
	 Mat3x3 RRlocTot(RRTot*RTw);       
	 Mat3x3 RRlocTotT(RRlocTot.Transpose());

	 /* Ruota velocita' e velocita' angolare nel sistema aerodinamico
	  * e li copia nel vettore di lavoro dW */
	 Vec3 Tmp(RRlocT*Vr);
	 Tmp.PutTo(dW);
	 
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	 if (fToBeOutput()) {
	    ASSERT(pOutput != NULL);
	    set_alpha(pTmpOutput, Tmp);
	 }
#endif /* AEROD_OUTPUT */
	 
	 Tmp = RRlocT*Wr;
	 Tmp.PutTo(dW+3);

	 /* Funzione di calcolo delle forze aerodinamiche */
	 __FC_DECL__(aerod2)(dW, dVAM, dTng, *pvd, &iInst,  &dOmega, &iProfile);
       
	 
	 /* OUTA */
	 
	 pvd++;  
	 
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	 if (fToBeOutput()) {
	    ASSERT(pOutput != NULL);
	    set_f(pTmpOutput, dTng);
	    pTmpOutput++;
	 }
#endif /* AEROD_OUTPUT */
	 
	 /* Dimensionalizza le forze */
	  doublereal dWght = PW.dGetWght();
	  Vec3 FTmp(RRloc*(Vec3(dTng)*(dL/2.*dWght)));
	  F += FTmp;
	  M += RRlocTot*(Vec3(dTng+3)*(dL/2.*dWght));
	  M += Xr.Cross(FTmp);
	 
          DEBUGCOUT("FsdrAer: "<<Vec3(dTng)<<endl);
	 
      } while(GDI.fGetNext(PW));   // fine ciclo sui punti di Gauss
      
      FTot += F;
      MTot += M;

      DEBUGCOUT("F: "<<F<<endl);
      DEBUGCOUT("M: "<<M<<endl);

    /* costruisco il vettore delle forze aerodinamiche (serve per il calcolo delle forze 
       aerodinamiche modali) nel sdr locale */
       Fa.Add((iAero-1)*2+1, (RnT*F).dGet(3));
       Fa.Add((iAero-1)*2+2, (RnTotT*M).dGet(1));  
   } // fine ciclo sugli elementi aerodinamici

   DEBUGCOUT("FTOT: "<<FTot<<endl);
   DEBUGCOUT("MTOT: "<<MTot<<endl);
   DEBUGCOUT("FModali: "<<endl);

    /* Forze aerodinamiche modali : {Fm} = [H ]T*{Fa} */   
    for(int iCnt = 1; iCnt <= NModes; iCnt ++) {
       double temp = 0.;
         for(int jCnt = 1; jCnt <= 2*NAeroElems; jCnt ++) {
            temp += pH->dGetCoef(jCnt, iCnt)*Fa.dGet(jCnt);
	 }
     WorkVec.fPutCoef(6+iCnt, temp); 
    }

   /* Se e' definito il rotore, aggiungere il contributo alla trazione */
   if(pRotor != NULL) {
      pRotor->AddForce(FTot, MTot, X0);
   }


   /* Sommare il termine al residuo */
   WorkVec.Add(1, FTot);
   WorkVec.Add(4, MTot);

}


/* output; si assume che ogni tipo di elemento sappia, attraverso
 * l'OutputHandler, dove scrivere il proprio output */
void AerodynamicModal::Output(OutputHandler& OH) const
{
   // DEBUGCOUTFNAME("AerodynamicModal::Output");
   
   /* Memoria in caso di forze instazionarie */
   if (iInst == 1) {
      doublereal d = 0.;
      if(pRotor != NULL) {
	 d = dDA*pRotor->dGetOmega();
      }
      if (fabs(d) < 1.e-6) {
	 d = 1.e-6;
      }
      for (integer i = 0; i < GDI.iGetNum(); i++) {
	 __FC_DECL__(coeprd)(&d, pvdOuta[i]);
      }
   } else if (iInst == 2) {
      for (integer i = 0; i < GDI.iGetNum(); i++) {
	 __FC_DECL__(coeprd)(&dDA, pvdOuta[i]);
      }
   }

   /* Output delle forze aerodinamiche F, M su apposito file */
   if (fToBeOutput()) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
      ASSERT(pOutput != NULL);
#endif /* AEROD_OUTPUT */
      ostream& out = OH.Aerodynamic() << setw(8) << GetLabel();
      
#if AEROD_OUTPUT == AEROD_OUT_NODE
      out << " " << setw(8) << pModalNode->GetLabel()
	<< " ", F.Write(out, " ") << " ", M.Write(out, " ");
#else /* AEROD_OUTPUT */      
      for (int i = 0; i < GDI.iGetNum(); i++) {
#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
	 out << ' ' << pOutput[i]->alpha << ' ' << pOutput[i]->f;
#elif AEROD_OUTPUT == AEROD_OUT_STD /* AEROD_OUTPUT */
	 for (int j = 1; j <= 6; j++) {
	    out << ' ' << pvdOuta[i][j];
	 }
#endif /* AEROD_OUTPUT */
      }
#endif /* AEROD_OUTPUT */
      out << endl;      
   }   
}

/* AerodynamicModal - end */


/* Funzioni di interpolazione */
doublereal ShapeFunc2N(doublereal d, integer iNode)
{
   ASSERT(iNode == 1 || iNode == 2);
   
   switch (iNode) {
	case 1:		
	  return .5*(1.-d);
				
	case 2:		
	  return .5*(1.+d);
		  			      	   
       break;
   }
   
   /* Per evitare warnings */
   return 0.;
}

doublereal DxDcsi2N(const Vec3& X1, const Vec3& X2)
{
   Vec3 DXDcsi(-X1*.5+X2*.5);
   doublereal dd = DXDcsi.Dot();
   if (dd > DBL_EPSILON) {
      return sqrt(dd);
   }
   return 0.;
}

/* Legge un elemento aerodinamico modale */

void ReadAeroModalData(DataManager* pDM,
		       MBDynParser& HP,
		       Shape** ppChord,
		       Shape** ppForce,
		       Shape** ppVelocity,
		       Shape** ppTwist,
		       integer* piInst,
		       integer* piNumber,
		       DriveCaller** ppDC,
		       integer* piProfile)
{
   DEBUGCOUTFNAME("ReadAeroModalData");
      
   /* Keywords */
   const char* sKeyWords[] = {    
	"naca0012",
	"rae9671"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,     
	NACA0012 = 0,
	RAE9671,
	LASTKEYWORD 
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   *ppChord = ReadShape(HP);
   *ppForce = ReadShape(HP);
   *ppVelocity = ReadShape(HP);
   *ppTwist = ReadShape(HP);
   
   HP.PutKeyTable(K);

   *piInst = HP.GetInt();
   if (*piInst != 0
#if 0
       && *piInst != 1 && *piInst != 2
#endif
       ) {      
      if (*piInst < 0 || *piInst > 2) {
	 cerr << "illegal unsteady flag";
      } else {
	 cerr << "unsteady aerodynamics are not tested yet";
      }
      cerr << " at line " << HP.GetLineData() << endl;
      THROW(ErrGeneric());
   }
   *piNumber = HP.GetInt();
   DEBUGLCOUT(MYDEBUG_INPUT, "PK flag flag: " << *piInst << endl
	      << "Gauss points number: " << *piNumber << endl);
  
   if (HP.IsKeyWord("control")) {      
      /* Driver di un'eventuale controllo */
      *ppDC = ReadDriveData(pDM, HP, pDM->pGetDrvHdl());
      HP.PutKeyTable(K);
   } else {      
      SAFENEWWITHCONSTRUCTOR(*ppDC,
			     NullDriveCaller,
			     NullDriveCaller(pDM->pGetDrvHdl()), EMmm);
   }
   
   *piProfile = 1;
   if (HP.fIsArg()) {            
      switch (HP.IsKeyWord()) {
       case NACA0012: {
	  DEBUGLCOUT(MYDEBUG_INPUT, "profile is NACA0012" << endl);
	  *piProfile = 1;
	  break;
       }
       case RAE9671: {
	  DEBUGLCOUT(MYDEBUG_INPUT, "profile is RAE9671" << endl);
	  *piProfile = 2;
	  break;
       }
       default: {
	  cerr << endl << "Unknown profile type at line "
	    << HP.GetLineData() << "; using default (NACA0012)" << endl;
	  *piProfile = 1;
	  break;
       }
      }
   }
}

Elem* ReadAerodynamicModal(DataManager* pDM,
			  MBDynParser& HP,
			  unsigned int uLabel)
{
   DEBUGCOUTFNAME("ReadAerodynamicModal");

   /* formato dell'input: 
    *
    *  label, 
    *  modal node,
    *  reference modal joint,   
    *  number of aerodynamic elements,
    *  lenght of aerodynamic elements,
    *  file gtkg,   
    *  (same as aerodynamic body)
    *  ...
    *  ;
    */             

      
   /* Nodo */
   StructNode* pModalNode = (StructNode*)pDM->ReadNode(HP, Node::STRUCTURAL);

   /* giunto modale collegato */		     
   Elem* pM = pDM->ReadElem(HP, Elem::JOINT);
   Modal* pModalJoint = (Modal*)pM->pGet();
   if (pModalJoint->GetJointType() != Joint::MODAL) {
      cerr << "element " << pModalJoint->GetLabel()
	      << " is required to be a modal joint" << endl;
      THROW(DataManager::ErrGeneric());
   }

   /* numero modi e nodi FEM */
   integer NModes = pModalJoint->uGetNModes();
   integer NFemNodes = pModalJoint->uGetNFemNodes();

   /* recupera gli autovettori */
   Mat3xN *pPHItTmp = NULL;
   Mat3xN *pPHIrTmp = NULL;
   
   SAFENEWWITHCONSTRUCTOR(pPHItTmp, Mat3xN, Mat3xN(NFemNodes*NModes, 0.), EMmm);
   SAFENEWWITHCONSTRUCTOR(pPHIrTmp, Mat3xN, Mat3xN(NFemNodes*NModes, 0.), EMmm);
   
   pPHItTmp = pModalJoint->pGetPHIt();
   pPHIrTmp = pModalJoint->pGetPHIr();
   
   /* Eventuale rotore */
   Rotor* pRotor = NULL;
   if (HP.IsKeyWord("rotor")) {
   
      /* verifica di esistenza del rotore       
       * NOTA: ovviamente il rotore deve essere definito 
       * prima dell'elemento aerodinamico */
      Elem* p = pDM->ReadElem(HP, Elem::ROTOR);
      DEBUGLCOUT(MYDEBUG_INPUT, "Linked to Rotor: " << p->GetLabel() << endl);
      pRotor = (Rotor*)p->pGet();
   }   

   /* legge il numero di elementi aerodinamici */
   integer NAeroElems = HP.GetInt();
   DEBUGLCOUT(MYDEBUG_INPUT, "Number of Aerodynamics Elements: " 
              << NAeroElems << endl);
	
   /* lunghezza degli elementi aerodinamici (per adesso
      e' costante, poi variera' */

   double dL = HP.GetReal();
   DEBUGLCOUT(MYDEBUG_INPUT, "Lenght : " 
              << dL << endl);
      
   
   /* legge la direzione dell'asse z (temoraneo)! */
   /*integer Zaxis = HP.GetInt();
     if (Zaxis <= 0 || Zaxis > 3) {
	  cerr << endl
	    << " at line " << HP.GetLineData() 
	    << " z axis is'nt valid" << endl;	
	  THROW(DataManager::ErrGeneric());
       }
   */
   /* apre il file contenente la matrice GTKG */
   const char *sFileGTKG = HP.GetFileName();
   ifstream fgtkg(sFileGTKG);       
   DEBUGCOUT("Reading Interpolation Matrix Data from file " << sFileGTKG << endl);
   if(!fgtkg) {
     cerr << endl << "Unable to open file " << sFileGTKG << endl;
     THROW(DataManager::ErrGeneric());
   }	
	      
   ReferenceFrame RF(pModalNode);

   Mat3x3 Ra(HP.GetRotRel(RF));

   Shape* pChord = NULL;
   Shape* pForce = NULL;
   Shape* pVelocity = NULL;
   Shape* pTwist = NULL;

   integer iInst = 0;
   integer iNumber = 0;
   DriveCaller* pDC = NULL;
   integer iProfile = 0;
   
   ReadAeroModalData(pDM, HP, 
		     &pChord, &pForce, &pVelocity, &pTwist,
		     &iInst, &iNumber, &pDC, &iProfile);


   /* la matrice 'H' e' quella che lega gli spostamenti dei punti aerodinamici
    * agli spostamenti modali, ed e' data dal prodotto della matrice che
    * lega gli spostamenti dei nodi aerodinamici con quelli dei nodi FEM per
    * la matrice degli autovettori :
    *  u = GTKG*q = GTKG*PHI*a = H*a, dove:
    *  u sono gli spostamenti dei punti aerodinamici
    *  q sono gli spostamenti dei nodi FEM
    *  a sono gli spostamenti modali
    * Questa matrice viene calcolata una sola volta e quindi utilizzata
    * continuamente per aggiornare la posizione dei nodi aerodinamici
    */  

   /* La matrice GTKG ha dimensioni [2 x NAeroElems] x [1 x NFemNodes] in quanto ogni strip
      ha due gradi di liberta' (traslazione verticale e torsione) e ogni nodo Fem ha solo la traslazione.
      La matrice PHI ha quindi dimensioni [NFemNodes x NModes]
    */

   FullMatrixHandler *pGTKG = NULL; 
   FullMatrixHandler PHI(NFemNodes, NModes, 0.);
   FullMatrixHandler *pH = NULL;

   SAFENEWWITHCONSTRUCTOR(pGTKG, FullMatrixHandler, FullMatrixHandler(2*NAeroElems, NFemNodes, 0.), EMmm);
   SAFENEWWITHCONSTRUCTOR(pH, FullMatrixHandler, FullMatrixHandler(2*NAeroElems, NModes, 0.), EMmm);

   doublereal d; 
   for(int iCnt = 1; iCnt <= 2*NAeroElems; iCnt ++) { 
        for(int jCnt = 1; jCnt <= NFemNodes; jCnt ++)    { 
          fgtkg >> d;
          pGTKG->fPutCoef(iCnt, jCnt, d);
        } 
    }
    fgtkg.close();  

    /* estraggo la matrice che contiene i modi di traslazione dei nodi FEM */
    for(int iMode=1; iMode<=NModes; iMode++) {     
          for(int iNode=1; iNode<=NFemNodes; iNode++) {     
	     //PHI.fPutCoef(iNode,iMode,pPHItTmp->dGet(Zaxis, (iMode-1)*NFemNodes+iNode));  
             PHI.fPutCoef(iNode,iMode,pPHItTmp->dGet(3, (iMode-1)*NFemNodes+iNode));   
      } 
    }

   pH->MatMul(*pGTKG, PHI);

   /* recupero gli offset dei nodi FEM */
   Mat3xN *pFemNodesPosition = NULL;
   SAFENEWWITHCONSTRUCTOR(pFemNodesPosition, Mat3xN, Mat3xN(NFemNodes, 0.), EMmm);
   pFemNodesPosition = pModalJoint->pGetFemNodesPosition();

   flag fOut = pDM->fReadOutput(HP, Elem::AERODYNAMIC);
   
   Elem* pEl = NULL;
   SAFENEWWITHCONSTRUCTOR(pEl, 
			  AerodynamicModal,
			  AerodynamicModal(uLabel, pModalNode, pModalJoint, pRotor, Ra,
					  pChord, pForce, pVelocity, pTwist,
					  iInst, iNumber, iProfile, NModes, NFemNodes,
                                          NAeroElems, pPHItTmp, pPHIrTmp, pH, pGTKG, dL,
                                          pFemNodesPosition, //Zaxis,
                                          pDC, fOut),
			  EMmm);
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << endl
	<< ": semicolon expected at line " << HP.GetLineData() << endl;      
      THROW(DataManager::ErrGeneric());
   }   
   
   return pEl;
} /* End of DataManager::ReadAerodynamicModal() */

