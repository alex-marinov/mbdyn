/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#ifdef USE_ELECTRIC_NODES

#include <discctrl.h>

/* DiscreteControlARXProcess_Debug - begin */

int DiscreteControlARXProcess_Debug::ReadMatrix(std::istream& In, 
						doublereal* pd,
						unsigned int iNumRows, 
						unsigned int iNumCols,
						unsigned int iNumSubMats,
						const char* sMatName)
{   
   /* Matrices alpha_i */
   for (unsigned int k = 0; k < iNumSubMats; k++) {  // for every k+1 = 1->p (every matrix)
      int iShift = iNumCols*k;
      for (unsigned int i = 0; i < iNumRows; i++) {  // for every row
	 doublereal* pdTmp = pd+iShift+iNumSubMats*iNumCols*i;
	 for (unsigned int j = 0; j < iNumCols; j++) {  // for every column
	    In >> pdTmp[j];
	    if (!In) {
	       silent_cerr("Error: unexpected end of stream while reading " 
		 << sMatName << '_'
		 << k+1 << '(' << i+1 << ',' << j+1 << ')' << std::endl);
	       
	       throw ErrGeneric();
	    }
	    
	    DEBUGLCOUT(MYDEBUG_INIT, sMatName << '_' << k+1 << "(" << i 
		       << "," << j << ") = " << pdTmp[j] << std::endl);
	 }
      }
   }
   
   return 0;
}


DiscreteControlARXProcess_Debug::DiscreteControlARXProcess_Debug(integer iNumOut,
								 integer iNumIn,
								 integer iOrdA,
								 integer iOrdB,
								 std::istream& In)
: iNumOutputs(iNumOut),
iNumInputs(iNumIn),
iOrderA(iOrdA),
iOrderB(iOrdB),
pdA(NULL),
pdY(NULL),
pdB(NULL),
pdU(NULL),
pdU0(NULL),
iRefA(0),
iRefB(0)
{
   /* The control model is handled as follows:
    * the outputs and inputs are stored in two vector:
    * Y = [y(k-1), ..., y(k-pa)], U = [u(k-1), ..., u(k-pb)]
    * as soon as new y, u are available, the oldes values are replaced by
    * the newest and the vector is used as a queue, i.e. it starts from
    * a floating reference and goes back to the beginning when it's
    * over.
    * Mathices alpha, beta are stacked in two big matrtices:
    * A = [alpha_1, ..., alpha_pa], B = [beta_1, ..., beta_pb]
    * the input is calculated as A*Y+B*U
    */
   
   
   /* At present at least one input and one output are required */
   ASSERT(iNumInputs > 0);
   ASSERT(iNumOutputs > 0);
   ASSERT(iOrderA > 0);
   ASSERT(iOrderB > 0);
   
   /* Size of working area */
   int iSize = iNumInputs*iNumOutputs*iOrderA     // Matrix A
     +iNumInputs*iNumInputs*iOrderB               // Matrix B
     +iNumOutputs*iOrderA                         // Vector Y
     +iNumInputs*iOrderB                          // Vector U
     +iNumInputs;                                 // Vector u(k)
     
   SAFENEWARR(pdA, doublereal, iSize);
   
   pdB = pdA+iNumInputs*iNumOutputs*iOrderA;
   pdY = pdB+iNumInputs*iNumInputs*iOrderB;
   pdU = pdY+iNumOutputs*iOrderA;
   pdU0 = pdU+iNumInputs*iOrderB;
   
   for (int i = iSize; i-- > 0; ) {
      pdA[i] = 0.;
   }
   
   /* In the input file, for human readability matrices are written as:
    * alpha_1,
    * alpha_2,
    * ...
    * alpha_p,
    * beta_1,
    * ...
    * beta_p
    */

   ReadMatrix(In, pdA, iNumInputs, iNumOutputs, iOrderA, "alpha");
   ReadMatrix(In, pdB, iNumInputs, iNumInputs, iOrderB, "beta");   
}


DiscreteControlARXProcess_Debug::~DiscreteControlARXProcess_Debug(void)
{
   /* It must be non-null;
    * matrices and arrays come from one allocation only */
   SAFEDELETEARR(pdA);
}


void DiscreteControlARXProcess_Debug::GetInput(doublereal* pdIn)
{
   for (int i = iNumInputs; i-- > 0; ) {
      pdIn[i] = pdU0[i];
   }
}


void DiscreteControlARXProcess_Debug::PutOutput(doublereal* pdOut, 
						doublereal* pdIn,
						doublereal* /* pdDesiredOut */ )
{
   /* At present pdDesiredOutput is not used;
    * it will be when the system is controlled
    * and ID'd at the same time */
   
   /* Moves reference backwards */
   if (iRefA == 0) {
      iRefA = iOrderA;
   }
   iRefA--;
   
   doublereal* pdOff = pdY+iNumOutputs*iRefA;
   for (int i = iNumOutputs; i-- > 0; ) {
      pdOff[i] = pdOut[i];
   }

   /* Moves reference backwards */
   if (iRefB == 0) {
      iRefB = iOrderB;
   }
   iRefB--;
   
   pdOff = pdU+iNumInputs*iRefB;
   for (int i = iNumInputs; i-- > 0; ) {
      /* pdIn contiene la somma dell'ingresso esogeno e dell'ingresso 
       * di controllo, ovvero l'ingresso totale applicato al sistema;
       * pdU0 contiene il solo ingresso generato dal controllore */
      pdOff[i] = pdIn[i];
      pdU0[i] = 0.; 
   }
   
   
   /* Matrices are sorted by rows */
   doublereal* pdMatOff = pdA+iOrderA*iNumOutputs*iNumInputs;
   doublereal* pdVecOff = pdY+iRefA*iNumOutputs; // Y shiftato di iRefA
   for (int i = iNumInputs; i-- > 0; ) {
      for (int j = iRefA*iNumOutputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdY[j];
      }
      for (int j = (iOrderA-iRefA)*iNumOutputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdVecOff[j];
      }
   }

   pdMatOff = pdB+iOrderB*iNumInputs*iNumInputs;
   pdVecOff = pdU+iRefB*iNumInputs; 
   for (int i = iNumInputs; i-- > 0; ) {
      for (int j = iRefB*iNumInputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdU[j];
      }
      for (int j = (iOrderB-iRefB)*iNumInputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdVecOff[j];
      }
   }
}

/* DiscreteControlProcess_Debug - end */


/* DiscreteIdentProcess_Debug - begin */

DiscreteIdentProcess_Debug::DiscreteIdentProcess_Debug(integer iNumOut,
						       integer iNumIn,
						       integer iOrdA, 
						       integer iOrdB,
						       ForgettingFactor* pf,
						       PersistentExcitation* px,
						       flag f_armax,
						       const char* sf)
: iNumOutputs(iNumOut), iNumInputs(iNumIn),
iOrderA(iOrdA), iOrderB(iOrdB), pId(NULL), pPx(px),
fout(sf != NULL ? 1 : 0)
{
   ASSERT(pf != NULL);
   ASSERT(px != NULL);

   if (f_armax) {     
      SAFENEWWITHCONSTRUCTOR(pId, 
			     IdentARMAXProcess,
			     IdentARMAXProcess(iNumOut, iNumIn, 
					       iOrdA, iOrdB, pf));
   } else {
      SAFENEWWITHCONSTRUCTOR(pId,
			     IdentARXProcess,
			     IdentARXProcess(iNumOut, iNumIn,
					     iOrdA, iOrdB, pf));
   }
   
   // temporaneo ?!?
   if (sf != NULL) {
      out.open(sf);
   }
}

DiscreteIdentProcess_Debug::~DiscreteIdentProcess_Debug(void)
{
   SAFEDELETE(pId);
   SAFEDELETE(pPx);
   
   // temporaneo ?!?
   if (fout) {      
      out.close();
   }   
}


/* Returns the new control input values in array pdIn */
void DiscreteIdentProcess_Debug::GetInput(doublereal* pdIn) 
{
   /* azzera per sicurezza */
   for (integer i = iNumInputs; i-- > 0; ) {
      pdIn[i] = 0.;
   }
   pPx->AddInput(pdIn);
}


const unsigned long BUFSIZE = 102400;
static doublereal buf[BUFSIZE];

/* Sets the new measures (and eventually the input) */
void DiscreteIdentProcess_Debug::PutOutput(doublereal* pdOut,
					      doublereal* pdIn,
					      doublereal* /* pdDesiredOut */ )
{
   pId->Update(pdOut, pdIn);

   if (fout) {
      integer size = pId->iGetSize()*pId->iGetNumOutput();
      if (size*sizeof(doublereal) > BUFSIZE) {
	 silent_cerr("buffer is too small" << std::endl);
      } else {      	 
	 pId->GetTheta(buf);
	 for (integer i = 0; i < pId->iGetSize()*pId->iGetNumOutput(); i++) {
	    out << std::setw(16) << buf[i];
	 }
	 out << std::setw(16) << pId->dGetForgettingFactor() << std::endl;
      }      
   }  
}

/* DiscreteIdentProcess_Debug - end */


/* DAC_Process_Debug - begin */

DAC_Process_Debug::DAC_Process_Debug(integer iNumOut,
				     integer iNumIn,
				     integer iOrdA,
				     integer iOrdB,
				     ForgettingFactor* pf,
				     GPCDesigner* pd,
				     PersistentExcitation* px,
				     DriveCaller* pTrig,
				     DriveCaller** pvDesOut,
				     const char* sf,
				     flag f)
: iNumOutputs(iNumOut),
iNumInputs(iNumIn),
iOrderA(iOrdA),
iOrderB(iOrdB),
iOrderMd(0),
pdBase(NULL),
pdTheta(NULL),
pdA(NULL),
pdY(NULL),
pdB(NULL),
pdU(NULL),
f_ma(f),
pdC(NULL),
pdE(NULL),
pdMd(NULL),
pdYd(NULL),
pdU0(NULL),
iRefA(0),
iRefB(0),
iRefMd(0),
pId(NULL), pCD(pd), pPx(px), Trigger(pTrig), 
f_md(pvDesOut != NULL ? 1 : 0), pvDesiredOut(NULL),
fout(sf != NULL ? 1 : 0)
{
   /* The control model is handled as follows:
    * the outputs and inputs are stored in two vector:
    * Y = [y(k-1), ..., y(k-pa)], U = [u(k-1), ..., u(k-pb)]
    * as soon as new y, u are available, the oldes values are replaced by
    * the newest and the vector is used as a queue, i.e. it starts from
    * a floating reference and goes back to the beginning when it's
    * over.
    * Mathices alpha, beta are stacked in two big matrtices:
    * A = [alpha_1, ..., alpha_pa], B = [beta_1, ..., beta_pb]
    * the input is calculated as A*Y+B*U
    */
   
   
   /* At present at least one input and one output are required */
   ASSERT(iNumInputs > 0);
   ASSERT(iNumOutputs > 0);
   ASSERT(iOrderA > 0);
   ASSERT(iOrderB > 0);
      
   ASSERT(pCD != NULL);
   ASSERT(pPx != NULL);
   
   iOrderMd = pCD->iGetPredStep()-pCD->iGetPredHor();
   ASSERT(iOrderMd > 0);
   
   /* Size of working area */
   int iSize = 
     iNumOutputs*iOrderA                                        // Vector Y
     +iNumInputs*iOrderB                                        // Vector U
     +iNumOutputs*iOrderMd                                      // Vector Yd
     +iNumInputs                                                // Vector u(k)
     +iNumOutputs*(iNumOutputs*iOrderA+iNumInputs*(iOrderB+1)); // Theta
     
   if (f_ma) {      
      iSize += iNumOutputs*iOrderA                        // Vector E
	+iNumOutputs*iNumOutputs*iOrderA;                 // Theta e' piu' grande!
   }
   
     
   SAFENEWARR(pdBase, doublereal, iSize);
   
   pdY = pdBase;
   pdU = pdY+iNumOutputs*iOrderA;
   pdYd = pdU+iNumInputs*iOrderB;
   pdU0 = pdYd+iNumOutputs*iOrderMd;
   pdTheta = pdU0+iNumInputs;
   
   if (f_ma) {   
      pdE = pdU0+iNumInputs;
      pdTheta = pdE+iNumOutputs*iOrderA;
   }
   
   for (int i = iSize; i-- > 0; ) {
      pdBase[i] = 0.;
   }
   
   /* In the input file, for human readability matrices are written as:
    * alpha_1,
    * alpha_2,
    * ...
    * alpha_p,
    * beta_1,
    * ...
    * beta_p
    */

   /* Si costruisce l'identificatore ARX */
   switch (f_ma) {
    case 0:
      SAFENEWWITHCONSTRUCTOR(pId,
			     IdentARXProcess,
			     IdentARXProcess(iNumOut, iNumIn,
					     iOrdA, iOrdB, pf));
      break;
    case 1:
      SAFENEWWITHCONSTRUCTOR(pId,
			     IdentARMAXProcess,
			     IdentARMAXProcess(iNumOut, iNumIn,
					       iOrdA, iOrdB, pf));
      break;
    default:
      silent_cerr("Unknown type of identification!" << std::endl);
      throw ErrGeneric();
   }
   
   if (f_md) {
      SAFENEWARR(pvDesiredOut, DriveOwner*, iNumOutputs);
      
      for (integer i = 0; i < iNumOutputs; i++) {
	 pvDesiredOut[i] = NULL;
	 SAFENEWWITHCONSTRUCTOR(pvDesiredOut[i],
				DriveOwner, 
				DriveOwner(pvDesOut[i]));
      }
      SAFEDELETEARR(pvDesOut);
   }
   
   // temporaneo ?!?
   if (fout) {
      out.open(sf);
      // mettere un test?
   }
}


DAC_Process_Debug::~DAC_Process_Debug(void)
{
   
   /* matrices and arrays come from one allocation only */
   
   SAFEDELETEARR(pdBase);          /* must be non-null */
   SAFEDELETE(pId);                /* must be non-null */
   SAFEDELETE(pCD);                /* must be non-null */
   SAFEDELETE(pPx);                /* must be non-null */
   if (pvDesiredOut != NULL) {           /* can be null (no desired output */
      for (integer i = iNumOutputs; i-- > 0; ) {
	 SAFEDELETE(pvDesiredOut[i]);
      }
      SAFEDELETEARR(pvDesiredOut);
   }

   // temporaneo ?!?
   if (fout) {      
      out.close();
   }   
}


void DAC_Process_Debug::GetInput(doublereal* pdIn)
{
   /* control input */
   for (int i = iNumInputs; i-- > 0; ) {
      pdIn[i] = pdU0[i];
   }
   
   /* persistent excitation */
   pPx->AddInput(pdIn);
}


void DAC_Process_Debug::PutOutput(doublereal* pdOut, 
				  doublereal* pdIn,
				  doublereal* pdDesiredOut)
{
   /* At present pdDesiredOutput is not used;
    * it will be when the system is controlled
    * and ID'd at the same time */
   
   ASSERT(pId != NULL);
   ASSERT(pCD != NULL);
   
   pId->Update(pdOut, pdIn);
   
   if (fout) {
      pId->GetTheta(pdTheta);
           
      for (integer i = 0; i < pId->iGetSize()*pId->iGetNumOutput(); i++) {
	 out << std::setw(16) << pdTheta[i];
      }
      out << std::setw(16) << pId->dGetForgettingFactor() << std::endl;
   }
      
   if (Trigger.dGet()) {
      pCD->DesignControl(pdTheta, &pdA, &pdB, &pdMd, &pdC);
   }
   
   /* Moves reference backwards */
   if (iRefA == 0) {
      iRefA = iOrderA;
   }
   iRefA--;
   
   doublereal* pdOff = pdY+iNumOutputs*iRefA;
   for (int i = iNumOutputs; i-- > 0; ) {
      pdOff[i] = pdOut[i];
   }    
   
   if (f_ma) {
      pdOff = pdE+iNumOutputs*iRefA;
      pId->GetErr(pdOff);     
   }   
   
   /* Moves reference backwards */
   if (iRefB == 0) {
      iRefB = iOrderB;
   }
   iRefB--;
   
   pdOff = pdU+iNumInputs*iRefB;
   for (int i = iNumInputs; i-- > 0; ) {
      /* pdIn contiene la somma dell'ingresso esogeno e dell'ingresso 
       * di controllo, ovvero l'ingresso totale applicato al sistema;
       * pdU0 contiene il solo ingresso generato dal controllore */
      pdOff[i] = pdIn[i];
      pdU0[i] = 0.; 
   }

   /* Moves reference backwards */
   if (f_md) {      
      if (iRefMd == 0) {
	 iRefMd = iOrderMd;
      }
      iRefMd--;
      
      doublereal* pdOff = pdYd+iNumOutputs*iRefMd;
      for (int i = iNumOutputs; i-- > 0; ) {       
	 pdOff[i] = pvDesiredOut[i]->dGet();
      }
      
      if (pdDesiredOut != NULL) {	
	 for (int i = iNumOutputs; i-- > 0; ) {
	    pdOff[i] += pdDesiredOut[i];
	 }
      }
   }

   if (pdB == NULL) {
      return;
   }
   
   /* Matrices are sorted by rows */
   doublereal* pdMatOff = pdA+iOrderA*iNumOutputs*iNumInputs;
   doublereal* pdVecOff = pdY+iRefA*iNumOutputs; // Y shiftato di iRefA
   for (int i = iNumInputs; i-- > 0; ) {
      for (int j = iRefA*iNumOutputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdY[j];
      }
      for (int j = (iOrderA-iRefA)*iNumOutputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdVecOff[j];
      }
   }
   
   pdMatOff = pdB+iOrderB*iNumInputs*iNumInputs;
   pdVecOff = pdU+iRefB*iNumInputs; 
   for (int i = iNumInputs; i-- > 0; ) {
      for (int j = iRefB*iNumInputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdU[j];
      }
      for (int j = (iOrderB-iRefB)*iNumInputs; j-- > 0; ) {
	 pdMatOff--;
	 pdU0[i] += pdMatOff[0]*pdVecOff[j];
      }
   }
   
   if (f_ma) {
      pdMatOff = pdC+iOrderA*iNumOutputs*iNumInputs;
      pdVecOff = pdE+iRefA*iNumOutputs; // E shiftato di iRefA
      for (int i = iNumInputs; i-- > 0; ) {
	 for (int j = iRefA*iNumOutputs; j-- > 0; ) {
	    pdMatOff--;
	    pdU0[i] += pdMatOff[0]*pdE[j];
	 }
	 for (int j = (iOrderA-iRefA)*iNumOutputs; j-- > 0; ) {
	    pdMatOff--;
	    pdU0[i] += pdMatOff[0]*pdVecOff[j];
	 }
      }
   }
   
   if (f_md) {
      pdMatOff = pdMd+iOrderMd*iNumOutputs*iNumInputs;
      pdVecOff = pdYd+iRefMd*iNumOutputs; // Yd shiftato di iRefMd
      for (int i = iNumInputs; i-- > 0; ) {
	 for (int j = iRefMd*iNumOutputs; j-- > 0; ) {
	    pdMatOff--;
	    pdU0[i] += pdMatOff[0]*pdYd[j];
	 }
	 for (int j = (iOrderMd-iRefMd)*iNumOutputs; j-- > 0; ) {
	    pdMatOff--;
	    pdU0[i] += pdMatOff[0]*pdVecOff[j];
	 }
      }
   }
}

/* DAC_Process_Debug - end */
   
   

















/* DiscreteControlElem - begin */

DiscreteControlElem::DiscreteControlElem(unsigned int uL, 
					 const DofOwner* pDO, 
					 integer iNumOut,
					 ScalarDof* pOut,
					 DriveCaller** ppOutSF,
					 integer iNumIn,
					 ScalarDof* pIn,
					 DiscreteControlProcess* p,
					 integer iNIt,
					 flag fOut)
: Elem(uL, fOut),
Electric(uL, pDO, fOut),
pDCP(p),
fNewStep(1),
iNumIter(iNIt),
iCurrIter(0),
iNumOutputs(iNumOut),
pOutputs(pOut),
pvOutScaleFact(NULL),
pdOut(NULL),
iNumInputs(iNumIn),
pInputs(pIn),
pdIn(NULL)
{
   /* The inputs should be all scalar dofs; an abtract force is 
    * applied to every node;
    * the outputs can be of every kind */
   
   
   /* Allocs and resets memory for input and output current values */
   ASSERT(iNumInputs > 0);
   ASSERT(iNumOutputs > 0);
   
   ASSERT(ppOutSF != NULL);
   SAFENEWARR(pvOutScaleFact, DriveOwner*, iNumOutputs);
   
   for (int i = iNumOutputs; i-- > 0; ) {
      ASSERT(ppOutSF[i] != NULL);
      pvOutScaleFact[i] = NULL;
      SAFENEWWITHCONSTRUCTOR(pvOutScaleFact[i],
			     DriveOwner, 
			     DriveOwner(ppOutSF[i]));
   }
   
   SAFENEWARR(pdIn, doublereal, iNumInputs+iNumOutputs);   
   pdOut = pdIn+iNumInputs;
   
   for (int i = iNumInputs+iNumOutputs; i-- > 0; ) {
      pdIn[i] = 0.;
   }
}


DiscreteControlElem::~DiscreteControlElem(void)
{
   if (pvOutScaleFact != NULL) {
      for (integer i = iNumOutputs; i-- > 0; ) {
	 SAFEDELETE(pvOutScaleFact[i]);
      }
      SAFEDELETEARR(pvOutScaleFact);
   }
   
   /* They must be non-null */  
   SAFEDELETEARR(pdIn);
   
   /* Allocated by DataManager */
   SAFEDELETEARR(pInputs);
   SAFEDELETEARR(pOutputs);
   SAFEDELETE(pDCP); /* Must destroy the other processes too */
}


/* Scrive il contributo dell'elemento al file di restart */
std::ostream& DiscreteControlElem::Restart(std::ostream& out) const
{
   out << "  electric: " << GetLabel() << ", discrete control;" << std::endl;
   return out;
}


void 
DiscreteControlElem::AfterConvergence(const VectorHandler& X,
		const VectorHandler& XP)
{
   /* Sets the flag for a new step */
   if (++iCurrIter == iNumIter) {      
      iCurrIter = 0;
      
      ASSERT(fNewStep == flag(0));
      fNewStep = flag(1);
      
      
      /* Gets output from solution */   
      for (int iCnt = iNumOutputs; iCnt-- > 0; ) {
	 pdOut[iCnt] = pvOutScaleFact[iCnt]->dGet()*pOutputs[iCnt].dGetValue();
      }
      
      /* Gets input from solution */     
      for (int iCnt = iNumInputs; iCnt-- > 0; ) {
	 pdIn[iCnt] = pInputs[iCnt].dGetValue();       
      }
      
      /* Puts measures to control subprocess */
      pDCP->PutOutput(pdOut, pdIn); /* DCP knows how long pdOut and pdIn are */
   }
      
   /* Add extra output as needed */
}


/* assemblaggio residuo */
SubVectorHandler& 
DiscreteControlElem::AssRes(SubVectorHandler& WorkVec,
			    doublereal /* dCoef */ ,
			    const VectorHandler& /* XCurr */ ,
			    const VectorHandler& /* XPrimeCurr */ )
{
   WorkVec.ResizeReset(iNumInputs);
 
   /* At first iteration gets the input from the control subprocess */
   if (fNewStep) {      
      /* Gets the inputs from subprocess */
      pDCP->GetInput(pdIn); /* DCP knows how long pdIn is*/
      
      /* Sets back the flag to zero */
      fNewStep = flag(0);
   }

   /* Sets the parameters */   
   for (int iCnt = iNumInputs; iCnt-- > 0; ) {
      /* Equivalent to an abstract force */
      WorkVec.PutItem(iCnt+1, pInputs[iCnt].pNode->iGetFirstRowIndex()+1, pdIn[iCnt]);     
   }
   
   return WorkVec;
}

/* DiscreteControlElem - end */

#endif /* USE_ELECTRIC_NODES */

