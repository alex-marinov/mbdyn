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

/* Rotore */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <rotor.h>
#include <dataman.h>

#ifdef USE_MPI
#ifdef USE_MYSLEEP
#include <mysleep.h>
#endif /* USE_MYSLEEP */

#ifdef MPI_PROFILING
extern "C" {
#include <mpe.h>
#include <stdio.h>
}
#endif /* MPI_PROFILING */
#endif /* USE_MPI */

extern "C" {
#include <mymath.h>
}

/* Rotor - begin */

Rotor::Rotor(unsigned int uL, const DofOwner* pDO,
	     const StructNode* pC, const StructNode* pR, flag fOut)
: Elem(uL, Elem::ROTOR, fOut), 
AerodynamicElem(uL, AerodynamicElem::ROTOR, fOut), 
ElemWithDofs(uL, Elem::ROTOR, pDO, fOut),
#ifdef USE_MPI
pBlockLenght(NULL),
pDispl(NULL),
ReqV(MPI::REQUEST_NULL),
pRotDataType(NULL),
#endif /* USE_MPI */
pCraft(pC), pRotor(pR), 
dOmegaRef(0.), dRadius(0.), dArea(0.), dUMean(0.), dUMeanPrev(0.), 
dWeight(0.), dCorrection(1.),
FTraction(0.), MTraction(0.),
RRotTranspose(0.), RRot3(0.), XCraft(0.), VCraft(0.),
dPsi0(0.), dSinAlphad(1.), dCosAlphad(0.),
dMu(0.), dLambda(1.), dChi(0.),
dVelocity(0.), dOmega(0.),
iNumSteps(0)
{
   ASSERT(pC != NULL);
   ASSERT(pC->GetNodeType() == Node::STRUCTURAL);
   ASSERT(pR != NULL);
   ASSERT(pR->GetNodeType() == Node::STRUCTURAL);
      
   Vec3 R3C((pCraft->GetRCurr()).GetVec(3));
   Vec3 R3R((pRotor->GetRCurr()).GetVec(3));
   if (R3C.Dot(R3R) < 1.-DBL_EPSILON) {
      cerr << "warning, possible misalignment of rotor and craft axes "
	      "for Rotor(" << GetLabel() << ")" << endl;
   }
#ifdef USE_MPI
   RotorComm = MPI::COMM_WORLD.Dup();
#endif /* USE_MPI */
}

Rotor::~Rotor(void)
{
   NO_OP;
}

/* Tipo dell'elemento (usato per debug ecc.) */
Elem::Type Rotor::GetElemType(void) const
{
    return Elem::ROTOR;
}

void Rotor::Output(OutputHandler& OH) const
{
    /* non mi ricordo a cosa serve! */
    ((int&)iNumSteps)++;

    /* updates the umean at the previous step */
    (doublereal&)dUMeanPrev = dUMean; 


    if (fToBeOutput()) {
#ifdef USE_MPI
        if (RotorComm.Get_size() > 1) { 
	    if (RotorComm.Get_rank() == 0) {
		Vec3 TmpF(TmpVecR), TmpM(TmpVecR+3);
		Mat3x3 RT((pCraft->GetRCurr()).Transpose());
		OH.Rotors() << setw(8) << GetLabel() << " " 
			<< (RT*TmpF) << " " << (RT*TmpM) << " " 
			<< dUMean 
#if 1
			<< " "
			<< dVelocity << " " 
			<< dSinAlphad << " " << dCosAlphad << " "
			<< dMu << " " << dLambda << " " << dChi << " "
			<< dPsi0
#endif
			<< endl;
	    }
	} else {
	    Mat3x3 RT((pCraft->GetRCurr()).Transpose());
	    OH.Rotors() << setw(8) << GetLabel() << " "
		    << (RT*FTraction) << " " << (RT*MTraction) << " " 
		    << dUMean 
#if 1
		    << " "
		    << dVelocity << " " 
		    << dSinAlphad << " " << dCosAlphad << " "
		    << dMu << " " << dLambda << " " << dChi << " "
		    << dPsi0
#endif 
		    << endl;
	}
#else /* !USE_MPI */     
	Mat3x3 RT((pCraft->GetRCurr()).Transpose());
	OH.Rotors() << setw(8) << GetLabel() << " " 
		<< (RT*FTraction) << " " << (RT*MTraction) << " " 
		<< dUMean 
#if 1
		<< " "
		<< dVelocity << " " 
		<< dSinAlphad << " " << dCosAlphad << " "
		<< dMu << " " << dLambda << " " << dChi << " "
		<< dPsi0
#endif
		<< endl;
#endif /* !USE_MPI */
    }
}

/* Calcola la posizione azimuthale di un punto generico.
 * X e' il punto di cui e' chiesta la posizione azimuthale, 
 *   nel sistema assoluto.
 * Gli viene sottratta la posizione del rotore (corpo attorno al quale avviene
 * la rotazione). Quindi il vettore risultante viene trasformato
 * nel riferimento del mozzo. 
 * Si trascura l'eventuale componente fuori del piano, dalle altre due si 
 * ricava l'angolo di rotazione relativa, a cui viene sommato l'angolo
 * che il corpo rotore forma con la direzione del vento relativo dPsi0, 
 * calcolata in precedenza.
 */
doublereal Rotor::dGetPsi(const Vec3& X) const
{
    Vec3 XRel(RRotTranspose*(X-XCraft));
    return dPsi0+atan2(XRel.dGet(2), XRel.dGet(1));
}

/* Calcola la distanza di un punto dall'asse di rotazione in coordinate 
 * adimensionali */
doublereal Rotor::dGetPos(const Vec3& X) const
{
   ASSERT(dRadius > 0.);
   Vec3 XRel(RRotTranspose*(X-XCraft));
   doublereal d1 = XRel.dGet(1);
   doublereal d2 = XRel.dGet(2);
   doublereal d = sqrt(d1*d1+d2*d2);

   ASSERT(dRadius > DBL_EPSILON);
   return d/dRadius;   
}

/* Calcola vari parametri geometrici
 * A partire dai corpi che identificano il velivolo ed il rotore
 */
void Rotor::InitParam(void)
{
   ASSERT(pCraft != NULL);
   ASSERT(pRotor != NULL);
   
   /* Trasposta della matrice di rotazione del rotore */
   RRotTranspose = pCraft->GetRCurr();
   RRot3 = RRotTranspose.GetVec(3);
   RRotTranspose = RRotTranspose.Transpose();
   
   /* Posizione del rotore */
   XCraft = pRotor->GetXCurr();

   /* Velocita' angolare del rotore */
   dOmega = (pRotor->GetWCurr()-pCraft->GetWCurr()).Norm();
   
   /* Velocita' di traslazione del velivolo */
   VCraft = -(pRotor->GetVCurr());
   Vec3 VTmp(0.);
   if (fGetAirVelocity(VTmp, pRotor->GetXCurr())) {
      VCraft += VTmp;
   }
  
   /* Velocita' nel sistema del velivolo (del disco?) decomposta */
   VTmp = RRotTranspose*VCraft;
   doublereal dV1 = VTmp.dGet(1);
   doublereal dV2 = VTmp.dGet(2);
   doublereal dV3 = VTmp.dGet(3);
   doublereal dVV = dV1*dV1+dV2*dV2;
   doublereal dV = sqrt(dVV);
       
   /* Angolo di azimuth del rotore */
   dPsi0 = atan2(-dV2, dV1);
   
   /* Angolo di influsso */
   dVelocity = sqrt(dV3*dV3+dVV);
   if (dVelocity > DBL_EPSILON) {
      dSinAlphad = -dV3/dVelocity;
      dCosAlphad = dV/dVelocity;
   } else {
      dSinAlphad = 1.;
      dCosAlphad = 0.;
   }

   /* Parametri di influsso (usano il valore di dUMean al passo precedente) */
   doublereal dVTip = 0.;
   dMu = 0.;
   dLambda = 0.;
   dVTip = dOmega*dRadius;
   if (dVTip > DBL_EPSILON) {
      dMu = (dVelocity*dCosAlphad)/dVTip;
      dLambda = (dVelocity*dSinAlphad+dUMeanPrev/dCorrection)/dVTip;
   }
   
   if (dMu == 0. && dLambda == 0.) {
      dChi = 0.;
   } else {      
      dChi = atan2(dMu, dLambda);
   }
}

/* Velocita' indotta media (uniforme) */
void Rotor::MeanInducedVelocity(void)
{
   /* Trazione nel sistema rotore */
   doublereal dT = RRot3*FTraction;

   /* Velocita' indotta media */
   doublereal dVRef = dOmega*dRadius*sqrt(dMu*dMu+dLambda*dLambda);
   doublereal d = 2.*dGetAirDensity(GetXCurr())*dArea*dVRef;
   doublereal dUMeanTmp = dCorrection*dT/(d+1.);

   dUMean = dUMeanTmp*(1.-dWeight)+dUMeanPrev*dWeight;
}

/* Azzera la trazione */
void Rotor::ResetTraction(void)
{
   FTraction = Vec3(0.);
   MTraction = Vec3(0.);
}

/* assemblaggio jacobiano (nullo per tutti tranne che per il DynamicInflow) */
VariableSubMatrixHandler& Rotor::AssJac(VariableSubMatrixHandler& WorkMat,
					doublereal /* dCoef */ ,
					const VectorHandler& /* XCurr */ ,
					const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering Rotor::AssJac()" << endl);
   WorkMat.SetNullMatrix();
   return WorkMat;	
}

ostream& Rotor::Restart(ostream& out) const
{
   return out << "  rotor: " << GetLabel() << ", " 
     << pCraft->GetLabel() << ", " << pRotor->GetLabel() 
     << ", induced velocity: ";
}

#ifdef USE_MPI
void Rotor::ExchangeTraction(flag fWhat)
{
#ifdef MPI_PROFILING
  MPE_Log_event(33, 0, "start RotorTrust Exchange ");
#endif /* MPI_PROFILING */
  /* Se il rotore è connesso ad una sola macchina non è necessario scambiare messaggi */
  if (RotorComm.Get_size() > 1){
    if (fWhat) {
      /* Scambia F e M */
      FTraction.PutTo(TmpVecS);
      MTraction.PutTo(TmpVecS+3);
      RotorComm.Allreduce(TmpVecS, TmpVecR, 6, MPI::DOUBLE, MPI::SUM);
      for (int i=0; i <= 2; i++) { 
        FTraction.pGetVec()[i] = TmpVecR[i];
        MTraction.pGetVec()[i] = TmpVecR[i+3];
      }
    } else {
      RotorComm.Allreduce(FTraction.pGetVec(), TmpVecR, 3, MPI::DOUBLE, MPI::SUM);
      for (int i=0; i <= 2; i++) { 
        FTraction.pGetVec()[i] = TmpVecR[i];
      }
    }
  }
#ifdef MPI_PROFILING
  MPE_Log_event(34, 0, "end RotorTrust Exchange ");
#endif /* MPI_PROFILING */
}

void Rotor::InitializeRotorComm(MPI::Intracomm* Rot)
{ 
  RotorComm = *Rot;
} 

void Rotor::ExchangeVelocity(void) 
{
  if (RotorComm.Get_size() > 1){
    if (RotorComm.Get_rank() == 0) {
      for (int i=1; i < RotorComm.Get_size(); i++) {
        RotorComm.Send(MPI::BOTTOM, 1, *pRotDataType, i, 100);
      }
    }
    else {
      ReqV = RotorComm.Irecv((void *)MPI::BOTTOM, 1, *pRotDataType, 0, 100);
    }
  }
}
#endif /* USE_MPI */
   
/* Rotor - end */


/* NoRotor - begin */

NoRotor::NoRotor(unsigned int uLabel,
		 const DofOwner* pDO,
		 const StructNode* pCraft, 
		 const StructNode* pRotor,
		 doublereal dR,
		 flag fOut)
: Elem(uLabel, Elem::ROTOR, fOut), 
Rotor(uLabel, pDO, pCraft, pRotor, fOut)
{
  dRadius = dR; /* puo' essere richiesto dal trim */
#ifdef USE_MPI
  if (fToBeOutput()) {
    SAFENEWARR(pBlockLenght, int, 3, DMmm);
    SAFENEWARR(pDispl, MPI::Aint, 3, DMmm);
    for (int i=0; i < 3; i++) {
      pBlockLenght[i] = 1;
    }
    for (int i=0; i < 3; i++) {	
      pDispl[i] = MPI::Get_address(&(XCraft.pGetVec()[i]));	  	
    }
    SAFENEWWITHCONSTRUCTOR(pRotDataType, MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(3, pBlockLenght, pDispl)), SMmm);
    pRotDataType->Commit();
  }
#endif /* USE_MPI */
}

NoRotor::~NoRotor(void)
{
   NO_OP;
}
   
/* assemblaggio residuo */
SubVectorHandler& NoRotor::AssRes(SubVectorHandler& WorkVec,
				  doublereal /* dCoef */ ,
				  const VectorHandler& /* XCurr */ ,
				  const VectorHandler& /* XPrimeCurr */ )
{
  DEBUGCOUT("Entering NoRotor::AssRes()" << endl);
   
#ifdef USE_MPI
  if (fToBeOutput()) {
    ExchangeTraction(fToBeOutput());
  }
  if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */   
    
     /* Velocita' angolare del rotore */
    Vec3 Omega(pRotor->GetWCurr()-pCraft->GetWCurr());
    dOmega = Omega.Dot();
    if (dOmega > DBL_EPSILON) {
      dOmega = sqrt(dOmega);
    }      
    
    if (fToBeOutput()) {      
      XCraft = pRotor->GetXCurr();
    }   
#ifdef USE_MPI 
    if (fToBeOutput()) {      
      ExchangeVelocity();
    }
  }
#endif /* USE_MPI */
  ResetTraction();
  WorkVec.Resize(0);
  return WorkVec;
}

ostream& NoRotor::Restart(ostream& out) const
{
  return Rotor::Restart(out) << "no;" << endl;
}

/* Somma alla trazione il contributo di forza di un elemento generico */
void NoRotor::AddForce(const Vec3& F, const Vec3& M, const Vec3& X)
{
  /*
   * Non gli serve in quanto non calcola velocita' indotta.
   * Solo se deve fare l'output lo calcola
   */
#ifdef USE_MPI
  if (ReqV != MPI::REQUEST_NULL) {
    flag RecvFlag;
    while (1) {
      RecvFlag = ReqV.Test();
      if (RecvFlag) { 
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(300);
#endif /* USE_MYSLEEP */
      }
    }
  }
#endif /* USE_MPI */

  if (fToBeOutput()) {
    FTraction += F;
    MTraction += M + (X - XCraft).Cross(F);
  }   
}

/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3 NoRotor::GetInducedVelocity(const Vec3& /* X */ ) const
{
  return Zero3;
}
  
/* NoRotor - end */


/* UniformRotor - begin */

UniformRotor::UniformRotor(unsigned int uLabel,
			   const DofOwner* pDO,
			   const StructNode* pCraft, 
			   const StructNode* pRotor,
			   doublereal dOR,
			   doublereal dR, 
			   doublereal dW,
			   doublereal dC,
			   flag fOut)
: Elem(uLabel, Elem::ROTOR, fOut), 
Rotor(uLabel, pDO, pCraft, pRotor, fOut)
{
   ASSERT(dOR > 0.);
   ASSERT(dR > 0.);
   ASSERT((dW <= 1.) && (dW >= 0.));
   
   dOmegaRef = dOR;
   dRadius = dR;
   dArea = M_PI*dRadius*dRadius;
   dWeight = dW;
   dCorrection = dC;

#ifdef USE_MPI
  SAFENEWARR(pBlockLenght, int, 7, DMmm);
  SAFENEWARR(pDispl, MPI::Aint, 7, DMmm);
  for (int i=0; i < 7; i++) {
    pBlockLenght[i] = 1;
  }
  for (int i=0; i < 3; i++) {	
    pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));	  	
  }
  pDispl[3] = MPI::Get_address(&dUMeanPrev);
  for (int i=4; i <= 6; i++) {	
    pDispl[i] = MPI::Get_address(&(XCraft.pGetVec()[i-4]));	  	
  }
  SAFENEWWITHCONSTRUCTOR(pRotDataType, MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(7, pBlockLenght, pDispl)), SMmm);
  pRotDataType->Commit();
  
#endif /* USE_MPI */
}

UniformRotor::~UniformRotor(void)
{
   NO_OP;
}

/* assemblaggio residuo */
SubVectorHandler& UniformRotor::AssRes(SubVectorHandler& WorkVec,
				       doublereal /* dCoef */ ,
				       const VectorHandler& /* XCurr */ ,
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering UniformRotor::AssRes()" << endl);  

#ifdef USE_MPI
   ExchangeTraction(fToBeOutput());
   if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */  


     /* Calcola parametri vari */
     Rotor::InitParam();   
     Rotor::MeanInducedVelocity();
   
     
#ifdef DEBUG   
     /* Prova:
	Vec3 XTmp(2.,2.,0.);
	doublereal dPsiTmp = dGetPsi(XTmp);
	doublereal dXTmp = dGetPos(XTmp);
	cout 
	<< "X rotore:  " << pRotor->GetXCurr() << endl
	<< "V rotore:  " << VCraft << endl
	<< "X punto:   " << XTmp << endl
	<< "Omega:     " << dOmega << endl
	<< "Velocita': " << dVelocity << endl
	<< "Psi0:      " << dPsi0 << endl
	<< "Psi punto: " << dPsiTmp << endl
	<< "Raggio:    " << dRadius << endl
	<< "r punto:   " << dXTmp << endl
	<< "mu:        " << dMu << endl
	<< "lambda:    " << dLambda << endl
	<< "cos(ad):   " << dCosAlphad << endl
	<< "sin(ad):   " << dSinAlphad << endl
	<< "UMean:     " << dUMean << endl;
     */        
#endif /* DEBUG */
   
#ifdef USE_MPI 
   }
   ExchangeVelocity();
#endif /* USE_MPI */   
   ResetTraction();
   
     /* Non tocca il residuo */
   WorkVec.Resize(0);
   return WorkVec;  
}

ostream& UniformRotor::Restart(ostream& out) const
{
  return Rotor::Restart(out) << "uniform, " << dRadius << ", " 
	  << dWeight << ", correction, " << dCorrection << ';' << endl;
}

/* Somma alla trazione il contributo di forza di un elemento generico */
void UniformRotor::AddForce(const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
  if (ReqV != MPI::REQUEST_NULL) {
    flag RecvFlag;
    while (1) {
      RecvFlag = ReqV.Test();
      if (RecvFlag) { 
	break;
#ifdef USE_MYSLEEP	
      } else {
	mysleep(300);
#endif /* USE_MYSLEEP */
      }
    }
  }
#endif /* USE_MPI */

   /* Gli serve solo la trazione */
  FTraction += F;
  
  /* Solo se deve fare l'output calcola anche il momento */
   if (fToBeOutput()) {      
     MTraction += M + (X - XCraft).Cross(F);
   }   
}

/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3 UniformRotor::GetInducedVelocity(const Vec3& /* X */ ) const
{
  return RRot3*dUMeanPrev;
};

/* UniformRotor - end */


/* GlauertRotor - begin */

GlauertRotor::GlauertRotor(unsigned int uLabel,
			   const DofOwner* pDO,
			   const StructNode* pCraft,
			   const StructNode* pRotor,
			   doublereal dOR,
			   doublereal dR, 
			   doublereal dW,
			   doublereal dC,
			   flag fOut)
: Elem(uLabel, Elem::ROTOR, fOut),
Rotor(uLabel, pDO, pCraft, pRotor, fOut)
{
   ASSERT(dOR > 0.);
   ASSERT(dR > 0.);
   ASSERT((dW <= 1.) && (dW >= 0.));
   
   dOmegaRef = dOR;
   dRadius = dR;
   dArea = M_PI*dRadius*dRadius;
   dWeight = dW;
   dCorrection = dC;

#ifdef USE_MPI
  SAFENEWARR(pBlockLenght, int, 20, DMmm);
  SAFENEWARR(pDispl, MPI::Aint, 20, DMmm);
  for (int i=0; i < 20; i++) {
    pBlockLenght[i] = 1;
  }
  for (int i=0; i < 3; i++) {	
    pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));	  	
  }
  pDispl[3] = MPI::Get_address(&dUMeanPrev);
  pDispl[4] = MPI::Get_address(&dLambda);
  pDispl[5] = MPI::Get_address(&dMu);
  pDispl[6] = MPI::Get_address(&dChi);
  pDispl[7] = MPI::Get_address(&dPsi0);
  for (int i=8; i <= 10; i++) {	
    pDispl[i] = MPI::Get_address(&(XCraft.pGetVec()[i-8]));	  	
  }
  for (int i=11; i < 20; i++) {	
    pDispl[i] = MPI::Get_address(&(RRotTranspose.pGetMat()[i-11]));	  	
  }
  SAFENEWWITHCONSTRUCTOR(pRotDataType, MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(20, pBlockLenght, pDispl)), SMmm);
  pRotDataType->Commit();
  
#endif /* USE_MPI */
}


GlauertRotor::~GlauertRotor(void)
{
   NO_OP;
}


/* assemblaggio residuo */
SubVectorHandler& GlauertRotor::AssRes(SubVectorHandler& WorkVec,
				       doublereal /* dCoef */ ,
				       const VectorHandler& /* XCurr */ , 
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering GlauertRotor::AssRes()" << endl);  

#ifdef USE_MPI
   ExchangeTraction(fToBeOutput());
   if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */

     /* Calcola parametri vari */
     Rotor::InitParam();   
     Rotor::MeanInducedVelocity();
#ifdef USE_MPI 
   }
   
   ExchangeVelocity();
#endif /* USE_MPI */

   /* Non tocca il residuo */
   ResetTraction();
   WorkVec.Resize(0);
   return WorkVec;  
}


ostream& GlauertRotor::Restart(ostream& out) const
{
   return Rotor::Restart(out) << "Glauert, " << dRadius << ", " 
     << dWeight << ';' << endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void GlauertRotor::AddForce(const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
  if (ReqV != MPI::REQUEST_NULL) {
    flag RecvFlag;
    while (1) {
      RecvFlag = ReqV.Test();
      if (RecvFlag) { 
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(300);
#endif /* USE_MYSLEEP */
      }
    }
  }
#endif /* USE_MPI */

   /* Gli serve solo la trazione */
   FTraction += F;
      
   /* Solo se deve fare l'output calcola anche il momento */
   if (fToBeOutput()) {      
      MTraction += M + (X - XCraft).Cross(F);
   }   
}


/*
 * Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale
 */
Vec3 GlauertRotor::GetInducedVelocity(const Vec3& X) const
{   
   if (dUMeanPrev == 0.) {
      return Vec3(0.);
   }
   
   if (fabs(dLambda) < 1.e-9) {
      return RRot3*dUMeanPrev;
   }
   
   doublereal dr = dGetPos(X);
   doublereal dp = dGetPsi(X);
   doublereal dd = 1.+4./3.*(1.-1.8*dMu*dMu)*tan(dChi/2.)*dr*cos(dp);
   
   return RRot3*(dd*dUMeanPrev);
};



/* GlauertRotor - end */


/* ManglerRotor - begin */

ManglerRotor::ManglerRotor(unsigned int uLabel,
			   const DofOwner* pDO,
			   const StructNode* pCraft, 
			   const StructNode* pRotor,
			   doublereal dOR,
			   doublereal dR, 
			   doublereal dW,
			   doublereal dC,
			   flag fOut)
: Elem(uLabel, Elem::ROTOR, fOut), 
Rotor(uLabel, pDO, pCraft, pRotor, fOut)
{
   ASSERT(dOR > 0.);
   ASSERT(dR > 0.);
   ASSERT((dW <= 1.) && (dW >= 0.));
   
   dOmegaRef = dOR;
   dRadius = dR;
   dArea = M_PI*dRadius*dRadius;
   dWeight = dW;
   dCorrection = dC;

#ifdef USE_MPI
  SAFENEWARR(pBlockLenght, int, 18, DMmm);
  SAFENEWARR(pDispl, MPI::Aint, 18, DMmm);
  for (int i=0; i < 18; i++) {
    pBlockLenght[i] = 1;
  }
  for (int i=0; i < 3; i++) {	
    pDispl[i] = MPI::Get_address(&(RRot3.pGetVec()[i]));	  	
  }
  pDispl[3] = MPI::Get_address(&dUMeanPrev);
  pDispl[4] = MPI::Get_address(&dSinAlphad);
  pDispl[5] = MPI::Get_address(&dPsi0);
  for (int i=6; i <= 8; i++) {	
    pDispl[i] = MPI::Get_address(&(XCraft.pGetVec()[i-6]));	  	
  }
  for (int i=9; i < 18; i++) {	
    pDispl[i] = MPI::Get_address(&(RRotTranspose.pGetMat()[i-9]));	  	
  }
  SAFENEWWITHCONSTRUCTOR(pRotDataType, MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(18, pBlockLenght, pDispl)), SMmm);
  pRotDataType->Commit();
  
#endif /* USE_MPI */
}


ManglerRotor::~ManglerRotor(void)
{
   NO_OP;
}


/* assemblaggio residuo */
SubVectorHandler& ManglerRotor::AssRes(SubVectorHandler& WorkVec,
				       doublereal /* dCoef */ ,
				       const VectorHandler& /* XCurr */ ,
				       const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering ManglerRotor::AssRes()" << endl);  

#ifdef USE_MPI
   ExchangeTraction(fToBeOutput());
   if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */

     /* Calcola parametri vari */
     Rotor::InitParam();   
     Rotor::MeanInducedVelocity();
   
     
#ifdef DEBUG   
     /* Prova: */
     Vec3 XTmp(2.,2.,0.);
     doublereal dPsiTmp = dGetPsi(XTmp);
     doublereal dXTmp = dGetPos(XTmp);
     Vec3 IndV = GetInducedVelocity(XTmp);
     cout 
       << "X rotore:  " << pRotor->GetXCurr() << endl
       << "V rotore:  " << VCraft << endl
       << "X punto:   " << XTmp << endl
       << "Omega:     " << dOmega << endl
       << "Velocita': " << dVelocity << endl
       << "Psi0:      " << dPsi0 << endl
       << "Psi punto: " << dPsiTmp << endl
       << "Raggio:    " << dRadius << endl
       << "r punto:   " << dXTmp << endl
       << "mu:        " << dMu << endl
       << "lambda:    " << dLambda << endl
       << "cos(ad):   " << dCosAlphad << endl
       << "sin(ad):   " << dSinAlphad << endl
       << "UMean:     " << dUMean << endl
       << "iv punto:  " << IndV << endl;
#endif /* DEBUG */

#ifdef USE_MPI 
   }
   ExchangeVelocity();
#endif /* USE_MPI */

   /* Non tocca il residuo */
   ResetTraction();
   WorkVec.Resize(0);
   return WorkVec;  
}
   

ostream& ManglerRotor::Restart(ostream& out) const
{
  return Rotor::Restart(out) << "Mangler, " << dRadius << ", " 
			     << dWeight << ';' << endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void ManglerRotor::AddForce(const Vec3& F, const Vec3& M, const Vec3& X)
{
#ifdef USE_MPI
  if (ReqV != MPI::REQUEST_NULL) {
    flag RecvFlag;
    while (1) {
      RecvFlag = ReqV.Test();
      if (RecvFlag) { 
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(300);
#endif /* USE_MYSLEEP */
      }
    }
  }
#endif /* USE_MPI */

  /* Gli serve solo la trazione */
  FTraction += F;
  
  /* Solo se deve fare l'output calcola anche il momento */
  if (fToBeOutput()) {      
    MTraction += M + (X - XCraft).Cross(F);
  }   
}


/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3 ManglerRotor::GetInducedVelocity(const Vec3& X) const
{
   if (dUMeanPrev == 0.) {
      return Vec3(0.);
   }
   
   doublereal dr = dGetPos(X);
   doublereal dp = dGetPsi(X);
   
   doublereal dr2 = dr*dr;
   doublereal dm2 = 1.-dr2;   
   doublereal dm = 0.;
   if (dm2 > 0.) {
      dm = sqrt(dm2);
   }
   doublereal da = 1.+dSinAlphad;
   if (fabs(da) > 1e-9) {
      da = (1.-dSinAlphad)/da;
   }
   if (fabs(da) > 0.) {
      da = sqrt(da);
   }
   
   doublereal dd = 15./4.*dm*dr2;
   
   /* Primo coefficiente */
   doublereal dc = -15./256.*M_PI*(9.*dr2-4.)*dr*da;
   dd -= 4.*dc*cos(dp);
   
   /* Secondo coefficiente */
   dc = 45./256.*M_PI*pow(da*dr, 3);
   dd -= 4.*dc*cos(3.*dp);
   
   /* Coefficienti pari, da 2 a 10: */
   for (int i = 2; i <= 10; i += 2) {
      dc = pow(-1., i/2-1)*15./8.
	*((dm+i)/(i*i-1.)*(3.-9.*dr2+i*i)+3.*dm)/(i*i-9.)
	*pow(((1.-dm)/(1.+dm))*da, i/2.);
      dd -= 4.*dc*cos(i*dp);
   }
         
   return RRot3*(dd*dUMeanPrev);
};



/* ManglerRotor - end */




const doublereal dM11 = 8./(3.*M_PI);
const doublereal dM22 = -16./(45.*M_PI);
const doublereal dM33 = -16./(45.*M_PI);

/* DynamicInflowRotor - begin */

DynamicInflowRotor::DynamicInflowRotor(unsigned int uLabel,
				       const DofOwner* pDO,
				       const StructNode* pCraft, 
				       const StructNode* pRotor,
				       doublereal dOR,
				       doublereal dR,
				       doublereal dVConstTmp,
				       doublereal dVCosineTmp,
				       doublereal dVSineTmp,
				       flag fOut)
: Elem(uLabel, Elem::ROTOR, fOut),
Rotor(uLabel, pDO, pCraft, pRotor, fOut),
dVConst(dVConstTmp), dVCosine(dVCosineTmp), dVSine(dVSineTmp), 
dL11(0.), dL13(0.), dL22(0.), dL31(0.), dL33(0.)
{
   ASSERT(dOR > 0.);
   ASSERT(dR > 0.);
   
   dOmegaRef = dOR;
   dRadius = dR;
   dArea = M_PI*dRadius*dRadius;
   
   /* Significa che valuta la velocita' indotta media al passo corrente */
   dWeight = 0.;

#ifdef USE_MPI
  SAFENEWARR(pBlockLenght, int, 20, DMmm);
  SAFENEWARR(pDispl, MPI::Aint, 20, DMmm);
  for (int i=0; i < 20; i++) {
    pBlockLenght[i] = 1;
  }
  for (int i=0; i < 3; i++) {	
    pDispl[i] = MPI::Get_address(RRot3.pGetVec()+i);	  	
  }
  pDispl[3] = MPI::Get_address(&dVConst);
  pDispl[4] = MPI::Get_address(&dVCosine);
  pDispl[5] = MPI::Get_address(&dVSine);
  pDispl[6] = MPI::Get_address(&dOmega);
  pDispl[7] = MPI::Get_address(&dPsi0);
  for (int i=8; i <= 10; i++) {	
    pDispl[i] = MPI::Get_address(XCraft.pGetVec()+i-8);	  	
  }
  for (int i=11; i < 20; i++) {	
    pDispl[i] = MPI::Get_address(RRotTranspose.pGetMat()+i-11);	  	
  }
  SAFENEWWITHCONSTRUCTOR(pRotDataType, MPI::Datatype, MPI::Datatype(MPI::DOUBLE.Create_hindexed(20, pBlockLenght, pDispl)), SMmm);
  pRotDataType->Commit();
#endif /* USE_MPI */
}


DynamicInflowRotor::~DynamicInflowRotor(void)
{
   NO_OP;
}


void DynamicInflowRotor::Output(OutputHandler& OH) const
{
   /* non mi ricordo a cosa serve! */
   ((int&)iNumSteps)++;

   (doublereal&)dUMeanPrev = dUMean; /* updates the umean at the previous step */

   /* FIXME: posso usare dei temporanei per il calcolo della trazione
    * totale per l'output, cosi' evito il giro dei cast */
#ifdef USE_MPI
 if (RotorComm.Get_size() > 1) {
   if (RotorComm.Get_rank() == 0) {
     if (fToBeOutput()) {
       Vec3 TmpF(TmpVecR), TmpM(TmpVecR+3);
       Mat3x3 RT((pCraft->GetRCurr()).Transpose());
       OH.Rotors() << setw(8) << GetLabel() << " " 
		   << (RT*TmpF) << " " << (RT*TmpM) << " " << dUMean << " "
		   << dVConst << " " << dVCosine << " " << dVSine  << endl; 
     }
   }
  } else {
    Mat3x3 RT((pCraft->GetRCurr()).Transpose());
	 OH.Rotors() << setw(8) << GetLabel() << " "
		     << (RT*FTraction) << " " << (RT*MTraction) << " " << dUMean << " "
		     << dVConst << " " << dVCosine << " " << dVSine  << endl;
  }
#else /* !USE_MPI */     
   if (fToBeOutput()) {
     Mat3x3 RT((pCraft->GetRCurr()).Transpose());
     OH.Rotors() << setw(8) << GetLabel() << " " 
		 << (RT*FTraction) << " " << (RT*MTraction) << " " 
		 << dUMean << " "
		 << dVConst << " " << dVCosine << " " << dVSine  << endl;
   } 
#endif /* !USE_MPI */     
}
     
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
DynamicInflowRotor::AssJac(VariableSubMatrixHandler& WorkMat,
			   doublereal dCoef,
			   const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   DEBUGCOUT("Entering DynamicInflowRotor::AssJac()" << endl);
#ifdef USE_MPI 
   if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */     
      SparseSubMatrixHandler& WM = WorkMat.SetSparse();
      WM.ResizeInit(5, 0, 0.);
   
      integer iFirstIndex = iGetFirstIndex();
      WM.fPutItem(1, iFirstIndex+1, iFirstIndex+1, dM11+dCoef*dL11);
      WM.fPutItem(2, iFirstIndex+3, iFirstIndex+1, dCoef*dL31);
      WM.fPutItem(3, iFirstIndex+2, iFirstIndex+2, dM22+dCoef*dL22);
      WM.fPutItem(4, iFirstIndex+1, iFirstIndex+3, dCoef*dL13);
      WM.fPutItem(5, iFirstIndex+3, iFirstIndex+3, dM33+dCoef*dL33);

#ifdef USE_MPI
   } else {
      WorkMat.SetNullMatrix();
   }
#endif /* USE_MPI */

   return WorkMat;
}


/* assemblaggio residuo */
SubVectorHandler& DynamicInflowRotor::AssRes(SubVectorHandler& WorkVec,
					     doublereal /* dCoef */ ,
					     const VectorHandler& XCurr, 
					     const VectorHandler& XPrimeCurr)
{
   DEBUGCOUT("Entering DynamicInflowRotor::AssRes()" << endl);
   
#ifdef USE_MPI
   ExchangeTraction(flag(1));
   
   if (RotorComm.Get_rank() == 0) {
#endif /* USE_MPI */

     /* Calcola parametri vari */
     Rotor::InitParam();   
     Rotor::MeanInducedVelocity();
     
#ifdef DEBUG
     /* Prova: */
     Vec3 XTmp(2.,2.,0.);
     doublereal dPsiTmp = dGetPsi(XTmp);
     doublereal dXTmp = dGetPos(XTmp);
     Vec3 IndV = GetInducedVelocity(XTmp);
     cout 
       << "X rotore:  " << pRotor->GetXCurr() << endl
       << "V rotore:  " << VCraft << endl
       << "X punto:   " << XTmp << endl
       << "Omega:     " << dOmega << endl
       << "Velocita': " << dVelocity << endl
       << "Psi0:      " << dPsi0 << endl
       << "Psi punto: " << dPsiTmp << endl
       << "Raggio:    " << dRadius << endl
       << "r punto:   " << dXTmp << endl
       << "mu:        " << dMu << endl
       << "lambda:    " << dLambda << endl
       << "cos(ad):   " << dCosAlphad << endl
       << "sin(ad):   " << dSinAlphad << endl
       << "UMean:     " << dUMean << endl
       << "iv punto:  " << IndV << endl;
#endif /* DEBUG */
     
     WorkVec.Resize(3);
     integer iFirstIndex = iGetFirstIndex();
   
     WorkVec.fPutRowIndex(1, iFirstIndex+1);
     WorkVec.fPutRowIndex(2, iFirstIndex+2);
     WorkVec.fPutRowIndex(3, iFirstIndex+3);
     
     dVConst = XCurr.dGetCoef(iFirstIndex+1);
     
     /* PROVA */
     dUMeanPrev = dUMean = dVConst*dOmega*dRadius;
     
     dVCosine = XCurr.dGetCoef(iFirstIndex+2);
     dVSine = XCurr.dGetCoef(iFirstIndex+3);
     
     doublereal dVConstPrime = XPrimeCurr.dGetCoef(iFirstIndex+1);
     doublereal dVCosinePrime = XPrimeCurr.dGetCoef(iFirstIndex+2);
     doublereal dVSinePrime = XPrimeCurr.dGetCoef(iFirstIndex+3);
     
     /* Trazione nel sistema rotore */
     Vec3 T(RRotTranspose*FTraction);
     
     /* Momento nel sistema rotore-vento */
     doublereal dCosP = cos(dPsi0);
     doublereal dSinP = sin(dPsi0);
     Mat3x3 RTmp( dCosP,  dSinP, 0., 
		  -dSinP,  dCosP, 0.,
		  0.,     0.,    1.);
     Vec3 M(RTmp*(RRotTranspose*MTraction));
       

     /* Ora la trazione non serve piu' */
     ResetTraction();
     
     /* Attenzione: moltiplico tutte le equazioni per dOmega */
     doublereal d = dGetAirDensity(GetXCurr())*dArea*dOmega*(dRadius*dRadius);
     if (d > DBL_EPSILON) {
       
       /* Coefficienti di trazione e momento */
       doublereal dCt = T.dGet(3)/d;
       d *= dRadius;
       doublereal dCl = M.dGet(1)/d;
       doublereal dCm = M.dGet(2)/d;
   
       /* Coefficienti della matrice */
       doublereal dUt = sqrt(dLambda*dLambda+dMu*dMu);
       doublereal dUm = 0.;
       if (dUt > DBL_EPSILON) { 
	 dUm = (dMu*dMu+dLambda*(dLambda+dVConst))/dUt;
       }
       
       d = dChi/2.;      
       doublereal dSinChi2 = sin(d);
       doublereal dCosChi2 = cos(d);
       
       d = 15./64.*M_PI*dSinChi2;
       doublereal dDen;
       dDen = 1.+d*d;
       
       dL11 = dOmega*(2*dUt/dDen);
       d = dOmega*(15./64.*M_PI*dSinChi2*dCosChi2/dDen);
       dL13 = d*dUt;
       dL31 = d*dUm;
       d = dOmega*(dCosChi2*dCosChi2/2.);
       dL22 = -d*dUm;
       dL33 = -d*dUm/dDen;
       
       WorkVec.fPutCoef(1, dCt-dM11*dVConstPrime-dL11*dVConst-dL13*dVSine);
       WorkVec.fPutCoef(2, dCl-dM22*dVCosinePrime-dL22*dVCosine);
       WorkVec.fPutCoef(3, dCm-dM33*dVSinePrime-dL31*dVConst-dL33*dVSine);
     } else {   
       WorkVec.fPutCoef(1, 0.);
       WorkVec.fPutCoef(2, 0.);
       WorkVec.fPutCoef(3, 0.);
     }
     
#ifdef USE_MPI 
     ExchangeVelocity();

   } else {

     /* Ora la trazione non serve piu' */
     ExchangeVelocity();
     ResetTraction();
     WorkVec.Resize(0);
   }
#endif /* USE_MPI */

   return WorkVec;
}
   

/* Relativo ai ...WithDofs */
void DynamicInflowRotor::SetInitialValue(VectorHandler& /* X */ ) const
{
    NO_OP;
}


/* Relativo ai ...WithDofs */
void DynamicInflowRotor::SetValue(VectorHandler& X, VectorHandler& XP) const
{
  integer iFirstIndex = iGetFirstIndex();
  
  for (int iCnt = 1; iCnt <= 3; iCnt++) {
    XP.fPutCoef(iFirstIndex+iCnt, 0.);
  }   
  
  X.fPutCoef(iFirstIndex+1, dVConst);
  X.fPutCoef(iFirstIndex+2, dVCosine);
  X.fPutCoef(iFirstIndex+3, dVSine);
}


/* Restart */
ostream& DynamicInflowRotor::Restart(ostream& out) const
{
   return Rotor::Restart(out) << "dynamic inflow, " << dRadius << ';' << endl;
}


/* Somma alla trazione il contributo di forza di un elemento generico */
void DynamicInflowRotor::AddForce(const Vec3& F, const Vec3& M, const Vec3& X)
{
   /*
    * Gli serve la trazione ed il momento rispetto al rotore, 
    * che si calcola da se'
    */
#ifdef USE_MPI
  if (ReqV != MPI::REQUEST_NULL) {
    flag RecvFlag;
    while (1) {
      RecvFlag = ReqV.Test();
      if (RecvFlag) { 
	break;
#ifdef USE_MYSLEEP
      } else {
	mysleep(300);
#endif /* USE_MYSLEEP */
      }
    }
  }
#endif /* USE_MPI */

   FTraction += F;
   MTraction += M + (X - XCraft).Cross(F);
}


/* Restituisce ad un elemento la velocita' indotta in base alla posizione
 * azimuthale */
Vec3 DynamicInflowRotor::GetInducedVelocity(const Vec3& X) const
{
   if (fabs(dVConst) <= DBL_EPSILON) {
      return Vec3(0.);
   }
   
   doublereal dr = dGetPos(X);
   doublereal dp = dGetPsi(X);
   
   return RRot3*((dRadius*dVConst+dr*(dVCosine*cos(dp)+dVSine*sin(dp)))*dOmega);
};

/* DynamicInflowRotor - end */


/* Legge un rotore */

Elem* ReadRotor(DataManager* pDM,
		MBDynParser& HP, 
		const DofOwner* pDO, 
		unsigned int uLabel)
{
   const char sFuncName[] = "ReadRotor()";
   DEBUGCOUT("Entering " << sFuncName << endl);
   
   const char* sKeyWords[] = {
      "inducedvelocity",
	"no",
	"uniform",
	"glauert",
	"mangler",
	"dynamicinflow"
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	INDUCEDVELOCITY = 0,
	NO,
	UNIFORM,
	GLAUERT,
	MANGLER,
	DYNAMICINFLOW,
	LASTKEYWORD 
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);
      
   /* Per ogni nodo: */
   
   /*     Velivolo */
   unsigned int uNode = (unsigned int)HP.GetInt();
   
   DEBUGCOUT("Craft Node: " << uNode << endl);
   
   /* verifica di esistenza del nodo */
   StructNode* pCraft;
   if ((pCraft = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << endl << sFuncName
	<< " at line " << HP.GetLineData() 
	<< ": craft structural node " << uNode
	<< " not defined" << endl;     
      THROW(DataManager::ErrGeneric());
   }		     
   
   uNode = (unsigned int)HP.GetInt();
   
   DEBUGCOUT("Rotor Node: " << uNode << endl);
   
   /* verifica di esistenza del nodo */
   StructNode* pRotor;
   if ((pRotor = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << endl << sFuncName
	<< " at line " << HP.GetLineData() 
	<< ": rotor structural node " << uNode
	<< " not defined" << endl;     
      THROW(DataManager::ErrGeneric());
   }		     
   
   /* Tipo di velocita' indotta */
   HP.ExpectDescription();
   if (KeyWords(HP.GetDescription()) != INDUCEDVELOCITY) {
      cerr << endl << "\"induced velocity\" expected at line " 
	<< HP.GetLineData() << endl;      
      THROW(DataManager::ErrGeneric());
   }

   Elem* pEl = NULL;
   
   KeyWords InducedType;
   switch (InducedType = KeyWords(HP.GetWord())) {
      
    case NO: {
       DEBUGCOUT("No induced velocity is considered" << endl);
       
       doublereal dR = 0.;
       if (HP.IsKeyWord("radius")) {
	  dR = HP.GetReal();
       }
       
       flag fOut = pDM->fReadOutput(HP, Elem::ROTOR);
       
       SAFENEWWITHCONSTRUCTOR(pEl,
			      NoRotor,
			      NoRotor(uLabel, pDO, pCraft, pRotor, dR, fOut),
			      DMmm);
       break;
    }
      
    case UNIFORM:
    case GLAUERT:
    case MANGLER:
    case DYNAMICINFLOW: {
       doublereal dOR = HP.GetReal();
       DEBUGCOUT("Reference rotation speed: " << dOR << endl);
       if (dOR <= 0.) {
	  cerr << endl << "Illegal null or negative reference speed for rotor "
	    << uLabel << " at line " << HP.GetLineData() << endl;
	  
	  THROW(DataManager::ErrGeneric());
       }
       
       doublereal dR = HP.GetReal();
       DEBUGCOUT("Radius: " << dR << endl);
       if (dR <= 0.) {
	  cerr << endl << "Illegal null or negative radius for rotor " 
	    << uLabel << " at line " << HP.GetLineData() << endl;	  
	  THROW(DataManager::ErrGeneric());
       }	   
       
       if (InducedType == DYNAMICINFLOW) {
	  DEBUGCOUT("Dynamic inflow induced velocity is considered" << endl);
	  
	  doublereal dVConst = 0.;
	  doublereal dVCosine = 0.;
	  doublereal dVSine = 0.;
	  if (HP.IsKeyWord("initial" "value")) {
	     dVConst = HP.GetReal();
	     dVCosine = HP.GetReal();
	     dVSine = HP.GetReal();
	  }	  
	  
	  flag fOut = pDM->fReadOutput(HP, Elem::ROTOR);
	  
	  /* Mettere qui il caso di Dynamic Inflow */
	  SAFENEWWITHCONSTRUCTOR(pEl, 
				 DynamicInflowRotor,
				 DynamicInflowRotor(uLabel, pDO, 
						    pCraft, pRotor, 
						    dOR, dR,
						    dVConst, dVCosine, dVSine,
						    fOut),
				 DMmm);
       } else {   	      	   	      	  
	  /* Legge il coefficiente di peso della velocita' indotta 
	   * ("weight" e' deprecato, si preferisce "delay") 
	   *
	   * nota:
	   * 
	   * U = U_n * ( 1 - dW ) + U_n-1 * dW
	   * 
	   * quindi dW rappresenta il peso che si da' al valore
	   * al passo precedente; in questo modo si introduce un
	   * ritardo euristico (attenzione: il ritardo vero dipende
	   * dal passo temporale) che aiuta ad evitare problemi 
	   * di convergenza. Se si desidera un ritardo "fisico",
	   * conviene forse provare il "Dynamic Inflow".
	   */
	  doublereal dW = 0.;
	  if (HP.IsKeyWord("weight") || HP.IsKeyWord("delay")) {
	     dW = HP.GetReal();
	     DEBUGCOUT("Weight: " << dW << endl);
	     if (dW < 0.) {
		cerr 
		  << "warning, illegal negative weight for uniform rotor "
		  << uLabel << ", switching to 0." 
		  << endl;
		dW = 0.;
	     } else if(dW > 1.) {
		cerr 
		  << "warning, illegal weight greater than 1. for uniform rotor "
		  << uLabel << ", switching to 1." 
		  << endl;
		dW = 1.;
	     }
	  }

	  /* Legge la correzione della velocita' indotta */
	  doublereal dC = 1.;
	  if (HP.IsKeyWord("correction")) {
	     dC = HP.GetReal();
	     DEBUGCOUT("Correction: " << dC << endl);
	     if (dC <= 0.) {
		cerr 
		  << "warning, illegal null or negative correction"
		  " for uniform rotor " << uLabel << ", switching to 1" 
		  << endl;
		dW = 1.;
	     }
	  }
	  
	  flag fOut = pDM->fReadOutput(HP, Elem::ROTOR);
	  
	  switch (InducedType) {
	   case UNIFORM: {		      
	      DEBUGCOUT("Uniform induced velocity is considered" << endl);
	      SAFENEWWITHCONSTRUCTOR(pEl, 
				     UniformRotor,
				     UniformRotor(uLabel, pDO, pCraft, pRotor,
						  dOR, dR, dW, dC, fOut), 
				     DMmm);
	      break;
	   }
	     
	   case GLAUERT: {		      
	      DEBUGCOUT("Glauert induced velocity is considered" << endl);
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     GlauertRotor,
				     GlauertRotor(uLabel, pDO, pCraft, pRotor,
						  dOR, dR, dW, dC, fOut), 
				     DMmm);
	      break;
	   }
	     
	   case MANGLER: {		      
	      DEBUGCOUT("Mangler induced velocity is considered" << endl);
	      
	      SAFENEWWITHCONSTRUCTOR(pEl,
				     ManglerRotor,
				     ManglerRotor(uLabel, pDO, pCraft, pRotor, 
						  dOR, dR, dW, dC, fOut), 
				     DMmm);
	      break;
	   }
	     
	   default: {
	      ASSERTMSG(0, "You shouldn't have reached this point");
	      THROW(DataManager::ErrGeneric());	    
	   }		 
	  }
       }
       
       break;
    }
      
    default: {
       cerr << endl << "unknown induced velocity type at line " 
	 << HP.GetLineData() << endl;       
       THROW(DataManager::ErrGeneric());
    }
   }
   
   /* Se non c'e' il punto e virgola finale */
   if (HP.fIsArg()) {
      cerr << endl << sFuncName
	<< ": semicolon expected at line " << HP.GetLineData() << endl;      
      THROW(DataManager::ErrGeneric());
   }   

   ASSERT(pEl != NULL);
   return pEl;
} /* End of DataManager::ReadRotor() */
