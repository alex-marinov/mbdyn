/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#include <mynewmem.h>
#include <strnode.h>
#include <body.h>
#include <autostr.h>
#include <dataman.h>

/*
 * StructNodeOutput - begin
 */
StructNodeOutput::~StructNodeOutput(void)
{
	NO_OP;
}

BasicStructNodeOutput::~BasicStructNodeOutput(void)
{
	NO_OP;
}

std::ostream&
BasicStructNodeOutput::Output(std::ostream& out, const StructNode *pN) const
{
	return out << pN->GetXCurr()
		<< " " << MatR2EulerAngles(pN->GetRCurr())*dRaDegr
		<< " " << pN->GetVCurr()
		<< " " << pN->GetWCurr()
		<< std::endl;
}

	StructNode *pBaseNode;

RelativeStructNodeOutput::RelativeStructNodeOutput(StructNode *pN)
: pBaseNode(pN)
{
	ASSERT(pBaseNode != NULL);
}

RelativeStructNodeOutput::~RelativeStructNodeOutput(void)
{
	NO_OP;
}

std::ostream&
RelativeStructNodeOutput::Output(std::ostream& out, const StructNode *pN) const
{
	Vec3 Xr = pN->GetXCurr() - pBaseNode->GetXCurr();
	Mat3x3 RT = pBaseNode->GetRCurr().Transpose();

	return out << RT*Xr
		<< " " << MatR2EulerAngles(RT*pN->GetRCurr())*dRaDegr
		<< " " << RT*(pN->GetVCurr() - pBaseNode->GetVCurr() - pBaseNode->GetWCurr().Cross(Xr))
		<< " " << RT*(pN->GetWCurr() - pBaseNode->GetWCurr())
		<< std::endl;
}

/*
 * StructNodeOutput - end
 */

/* StructNode - begin */

/* Costruttore definitivo */
StructNode::StructNode(unsigned int uL, 
		       const DofOwner* pDO,
		       const Vec3& X0, 
		       const Mat3x3& R0,
		       const Vec3& V0, 
		       const Vec3& W0,
		       doublereal dPosStiff,
		       doublereal dVelStiff, 
		       flag fOmRot, 
		       flag fOut)
: Node(uL, pDO, fOut),
RPrev(R0),
RRef(R0),
RCurr(R0),
gRef(0.),
gCurr(0.),
gPRef(0.),
gPCurr(0.),
XPrev(X0), 
XCurr(X0),
VPrev(V0), 
VCurr(V0), 
WPrev(W0),
WRef(W0),
WCurr(W0),
dPositionStiffness(dPosStiff), 
dVelocityStiffness(dVelStiff),
fOmegaRot(fOmRot) 
{
   NO_OP;
}


/* Distruttore (per ora e' banale) */
StructNode::~StructNode(void)
{
   NO_OP;
}


/* Tipo di nodo */
Node::Type 
StructNode::GetNodeType(void) const
{
   return Node::STRUCTURAL;
}


/* Contributo del nodo strutturale al file di restart */
std::ostream& 
StructNode::Restart(std::ostream& out) const
{
   out << "  structural: " << GetLabel() << ", ";
   if (GetStructNodeType() == StructNode::DYNAMIC) {	
      out << "dynamic";
   } else if (GetStructNodeType() == StructNode::STATIC) {	
	out << "static";
   }
   out << ", reference, global, ";
   XCurr.Write(out, ", ")
     << ", reference, global, 1, ", (RCurr.GetVec(1)).Write(out, ", ")
     << ", 2, ", (RCurr.GetVec(2)).Write(out, ", ") 
     << ", reference, global, ",
     VCurr.Write(out, ", ") 
     << ", reference, global, ",
     WCurr.Write(out, ", ") << ", assembly, "
     << dPositionStiffness << ", "
     << dVelocityStiffness << ", " 
     << fOmegaRot << ';' << std::endl;
   
   return out;
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal& 
StructNode::dGetDofValue(int iDof, int iOrder) const
{
   ASSERT(iDof >= 1 && iDof <= 6);
   ASSERT(iOrder == 0 || iOrder == 1);
   if (iDof >= 1 && iDof <= 3) {
      if (iOrder == 0) {
	 return XCurr.dGet(iDof);	 
      } else if (iOrder == 1) {	 
	 return VCurr.dGet(iDof);
      } 
   } else if (iDof >= 4 && iDof <= 6) {	        
      if (iOrder == 1) {
	 return WCurr.dGet(iDof-3);
      } else if (iOrder == 0) {
	 std::cerr << "Node " << GetLabel() 
	   << ": unable to return the angle." << std::endl;
	 THROW(StructNode::ErrGeneric());
      }	 
   } else {
      std::cerr << "Required dof is not available." << std::endl;
      THROW(StructNode::ErrGeneric());
   }

   /* dummy return value to workaround compiler complains */
   static doublereal dmy = 0.;
   return dmy;
}

/* Restituisce il valore del dof iDof al passo precedente;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal& 
StructNode::dGetDofValuePrev(int iDof, int iOrder) const
{
   ASSERT(iDof >= 1 && iDof <= 6);
   ASSERT(iOrder == 0 || iOrder == 1);
   if (iDof >= 1 && iDof <= 3) {
      if (iOrder == 0) {
	 return XPrev.dGet(iDof);	 
      } else if (iOrder == 1) {	 
	 return VPrev.dGet(iDof);
      } 
   } else if (iDof >= 4 && iDof <= 6) {	        
      if (iOrder == 1) {
	 return WPrev.dGet(iDof-3);
      } else if (iOrder == 0) {
	 std::cerr << "Node " << GetLabel() 
	   << ": unable to return the angle." << std::endl;
	 THROW(StructNode::ErrGeneric());
      }	 
   } else {
      std::cerr << "Required dof is not available." << std::endl;
      THROW(StructNode::ErrGeneric());
   }

   /* dummy return value to workaround compiler complains */
   static doublereal dmy = 0.;
   return dmy;
}

/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void 
StructNode::SetDofValue(const doublereal& dValue,
			unsigned int iDof, 
			unsigned int iOrder /* = 0 */ )
{
   ASSERT(iDof >= 1 && iDof <= 6);
   ASSERT(iOrder == 0 || iOrder == 1);
   if (iDof >= 1 && iDof <= 3) {
      if (iOrder == 0) {
	 XCurr.Put(iDof, dValue);
	 return;
      } else if (iOrder == 1) {	 
	 VCurr.Put(iDof, dValue);
	 return;
      } 
   } else if (iDof >= 4 && iDof <= 6) {	        
      if (iOrder == 1) {
	 WCurr.Put(iDof, dValue);
	 return;
      } else if (iOrder == 0) {
	 std::cerr << "Node " << GetLabel() 
	   << ": unable to set the angle." << std::endl;
	 THROW(StructNode::ErrGeneric());	 
      }	 
   } else {
      std::cerr << "Required dof is not available." << std::endl;
      THROW(StructNode::ErrGeneric());
   }   
}


DofOrder::Order 
StructNode::GetDofType(unsigned int i) const
{
	ASSERT(i >= 0 && i < iGetNumDof());
	return DofOrder::DIFFERENTIAL;
}


/* Output del nodo strutturale (da mettere a punto) */
void 
StructNode::Output(OutputHandler& OH) const
{
   if (fToBeOutput()) {      
      OH.StrNodes() << std::setw(8) << GetLabel() << " "
	<< XCurr << " " << MatR2EulerAngles(RCurr)*dRaDegr << " "
	<< VCurr << " " << WCurr << " " << std::endl;
   }   
}


/* Output della soluzione perturbata (modi ...) */
void
StructNode::Output(
		OutputHandler& OH, 
		const VectorHandler& X, 
		const VectorHandler& XP
		) const
{
   if (fToBeOutput()) {
      integer iFirstIndex = iGetFirstIndex();
      Vec3 DX(X, iFirstIndex+1);
      Vec3 Dg(X, iFirstIndex+4);
      Mat3x3 DR(MatR, Dg);

      OH.StrNodes() << std::setw(8) << GetLabel() << " "
	<< (XCurr+DX) << " " << MatR2EulerAngles(DR*RCurr)*dRaDegr << " "
	<< "#" << std::endl;
   }
}

/* Output di un modello NASTRAN equivalente nella configurazione corrente */
void 
StructNode::Output_pch(
		std::ostream& out
		) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		const char *name = GetName();

		out << "$ Node " << GetLabel();
		if (name) {
			out << " (" << name << ")";
		}
		
#define __NASTRAN_FORMAT__ __HACK_NASTRAN_MODES__

		Vec3 eZ = XCurr+RCurr.GetVec(3);
		Vec3 eX = XCurr+RCurr.GetVec(1);

#if __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED__
		out << std::endl
			/* CORD2R with node position and orientation */
			<< "CORD2R  " 
			<< std::setw(8) << GetLabel()
			<< std::setw(8) << 0
			<< std::setw(8) << XCurr.dGet(1)
			<< std::setw(8) << XCurr.dGet(2)
			<< std::setw(8) << XCurr.dGet(3)
			<< std::setw(8) << eZ.dGet(1)
			<< std::setw(8) << eZ.dGet(2)
			<< std::setw(8) << eZ.dGet(3)
			<< "+" << std::setw(1) << 1
			<< std::endl
			<< "+" << std::setw(7) << 1
			<< std::setw(8) << eX.dGet(1)
			<< std::setw(8) << eX.dGet(2)
			<< std::setw(8) << eX.dGet(3)
			<< std::endl
			<< "GRID    "
			<< std::setw(8) << GetLabel()
			<< std::setw(8) << GetLabel()
			<< std::setw(8) << 0.
			<< std::setw(8) << 0.
			<< std::setw(8) << 0.
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FIXED16__
		out << std::endl
			/* CORD2R with node position and orientation */
			<< "CORD2R* " 
			<< std::setw(16) << GetLabel()
			<< std::setw(16) << 0
			<< std::setw(16) << XCurr.dGet(1)
			<< std::setw(16) << XCurr.dGet(2)
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(16) << XCurr.dGet(3)
			<< std::setw(16) << eZ.dGet(1)
			<< std::setw(16) << eZ.dGet(2)
			<< std::setw(16) << eZ.dGet(3)
			<< "*" << std::setw(7) << 2
			<< std::endl
			<< "*" << std::setw(7) << 2
			<< std::setw(16) << eX.dGet(1)
			<< std::setw(16) << eX.dGet(2)
			<< std::setw(16) << eX.dGet(3)
			<< std::endl
			<< "GRID*   "
			<< std::setw(16) << GetLabel()
			<< std::setw(16) << GetLabel()
			<< std::setw(16) << 0.
			<< std::setw(16) << 0.
			<< "*" << std::setw(7) << 1
			<< std::endl
			<< "*" << std::setw(7) << 1
			<< std::setw(16) << 0.
			<< std::endl;
#elif __NASTRAN_FORMAT__ == __NASTRAN_FORMAT_FREE__
		out << std::endl
			/* CORD2R with node position and orientation */
			<< "CORD2R," << GetLabel()
			<< ",0,", XCurr.Write(out, ",")
			<< ",", eZ.Write(out, ",")
#if 0
			<< "," 
#endif
			<< std::endl
#if 1
			<< "," 
#endif
			<< " ", eX.Write(out, ",")
			<< std::endl
			/* grid in CORD2R */
			<< "GRID," << GetLabel() 
			<< "," << GetLabel()
			<< ",0.,0.,0." << std::endl;
#else
#error "unknown NASTRAN format"
#endif
	}
#endif /* __HACK_NASTRAN_MODES__ */
}


void 
StructNode::Output_f06(
		std::ostream& out,
		const VectorHandler& X
		) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		integer iFirstIndex = iGetFirstIndex();
		
		out 
			<< std::setw(13) << GetLabel() << "      G" 
			<< std::setw(18) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+1)
			<< std::setw(15) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+2)
			<< std::setw(15) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+3)

			<< std::setw(15) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+4)
			<< std::setw(15) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+5)
			<< std::setw(15) << std::setprecision(6) 
			<< X.dGetCoef(iFirstIndex+6)
			<< std::endl;
	}
#endif /* __HACK_NASTRAN_MODES__ */
}


void 
StructNode::Output_f06(
		std::ostream& out,
		const VectorHandler& Xr,
		const VectorHandler& Xi
		) const
{
#if defined(__HACK_NASTRAN_MODES__)
	if (fToBeOutput()) {
		integer iFirstIndex = iGetFirstIndex();
		
		out 
			<< "0" << std::setw(12) << GetLabel() << "      G" 
			<< std::setw(18) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+1)
			<< std::setw(15) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+2)
			<< std::setw(15) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+3)

			<< std::setw(15) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+4)
			<< std::setw(15) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+5)
			<< std::setw(15) << std::setprecision(6) 
			<< Xr.dGetCoef(iFirstIndex+6)
			<< std::endl

			<< std::setw(38) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+1)
			<< std::setw(15) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+2)
			<< std::setw(15) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+3)

			<< std::setw(15) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+4)
			<< std::setw(15) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+5)
			<< std::setw(15) << std::setprecision(6) 
			<< Xi.dGetCoef(iFirstIndex+6)
			<< std::endl;
	}
#endif /* __HACK_NASTRAN_MODES__ */
}


/* Aggiorna dati in base alla soluzione */
void 
StructNode::Update(const VectorHandler& X, const VectorHandler& XP)
{   
   integer iFirstIndex = iGetFirstIndex();
   
   XCurr = Vec3(X, iFirstIndex+1);
   VCurr = Vec3(XP, iFirstIndex+1);
    
   /* Nota: i g, gP non vengono incrementati */
   gCurr = Vec3(X, iFirstIndex+4);
   gPCurr = Vec3(XP, iFirstIndex+4);
   
#if 1
   /* Matrice RDelta, incremento di rotazione da predetto a corrente;
    Questo e' piu' efficiente */
   Mat3x3 RDelta(MatR, gCurr);
#else
   /* Questo e' meno efficiente anche se sembra piu' elegante.
    * Il problema e' che per scrivere il manipolatore in forma
    * elegante bisogna aggiungere alla matrice le informazioni
    * di memorizzazione della funzione di manipolazione.
    * Oppure occorre un operatore ternario */
    RDelta = MatR << gCurr;
#endif
  
#if 1
   /* La matrice di rotazione corrente e' data dalla matrice predetta 
    * (costante) moltiplicata per l'incremento totale occorso;
    * la velocita' angolare e' data dalla parte incrementale totale
    * piu' il contributo della velocita' di riferimento (costante) */
   RCurr = RDelta*RRef;   
   WCurr = Mat3x3(MatG, gCurr)*gPCurr+RDelta*WRef;
#else   
   /* Nuovo manipolatore (forse e' meno efficiente) */
   WCurr = (MatG << gCurr)*gPCurr+RDelta*WRef;
#endif
}


/* Aggiorna dati in base alla soluzione */
void 
StructNode::DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP)
{
   integer iFirstIndex = iGetFirstIndex();

   /* Forza configurazione e velocita' al valore iniziale */
   ((VectorHandler&)X).Put(iFirstIndex+1, XCurr);
   ((VectorHandler&)X).Put(iFirstIndex+4, Zero3);
   ((VectorHandler&)XP).Put(iFirstIndex+1, VCurr);
   ((VectorHandler&)XP).Put(iFirstIndex+4, Zero3);
      
   /* Questa parte e' spuria, va messa a punto (?) */
   XPrev = XCurr;
   RPrev = RCurr;
   VPrev = VCurr;
   WPrev = WCurr;   
}


/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
void 
StructNode::InitialUpdate(const VectorHandler& X)
{
   integer iFirstIndex = iGetFirstIndex();
   
   XCurr = Vec3(X, iFirstIndex+1);
   VCurr = Vec3(X, iFirstIndex+7);
   
   /* Nota: g viene incrementato */
   gCurr = Vec3(X, iFirstIndex+4);
   
#if 1
   /* Questo manipolatore e' piu' efficiente */
   Mat3x3 RDelta(MatR, gCurr);
#else
   /* Nuovo manipolatore (e' meno efficiente) */
   Mat3x3 RDelta(MatR << gCurr);
#endif
   
   RCurr = RDelta*RRef;
   WCurr = Vec3(X, iFirstIndex+10);
}

   
/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
void 
StructNode::SetInitialValue(VectorHandler& X) const
{
   integer iFirstIndex = iGetFirstIndex();
   X.Put(iFirstIndex+1, XCurr);
   X.Put(iFirstIndex+4, Vec3(0.));
   X.Put(iFirstIndex+7, VCurr);
   X.Put(iFirstIndex+10, WCurr);		
}


void 
StructNode::BeforePredict(VectorHandler& X, 
			  VectorHandler& XP,
			  VectorHandler& XPr, 
			  VectorHandler& XPPr) const
{
   /* Questa e' la predizione "consistente", ovvero usa come gdl di rotazione
    * i parametri di rotazione "totali" per predire la configurazione al nuovo
    * passo, quindi ritorna in forma incrementale */
   
   /* Calcolo la matrice RDelta riferita a tutto il passo trascorso
    * all'indietro */
   Mat3x3 RDelta(RPrev*RCurr.Transpose());
   
   /* Pongo la R al passo precedente uguale a quella corrente
    * mi servira' se devo ripetere il passo con un diverso Delta t 
    * e per la rettifica dopo la predizione */
   RPrev = RCurr;
   
   integer iFirstPos = iGetFirstIndex()+3;
   
   /* Mi assicuro che g, al passo corrente sia nullo */
   X.Put(iFirstPos+1, Vec3(0.));
   
   /* Calcolo g al passo precedente attraverso la matrice RDelta riferita
    * a tutto il passo. Siccome RDelta e' calcolata all'indietro,
    * i parametri soo gia' con il segno corretto */
   Vec3 gPrev = gparam(RDelta);
   XPr.Put(iFirstPos+1, gPrev);
   
   /* Calcolo gP al passo precedente attraverso la definizione 
    * mediante le Omega. Siccome i parametri sono con il segno meno 
    * e la matrice RDelta e' gia' calcolata all'indietro, l'insieme
    * e' consistente */   
   XPPr.Put(iFirstPos+1, Mat3x3(MatGm1, gPrev)*WPrev);

   /* Metto Omega al passo corrente come gP */
   XP.Put(iFirstPos+1, WCurr);
   
   /* Pongo le Omega al passo precedente uguali alle Omega al passo corrente
    * mi servira' per la correzione dopo la predizione */
   WPrev = WCurr;
   
   XPrev = XCurr;   
   VPrev = VCurr;
}

void 
StructNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
   integer iFirstIndex = iGetFirstIndex();
   
   /* Spostamento e velocita' aggiornati */
   XCurr = Vec3(X, iFirstIndex+1);
   VCurr = Vec3(XP, iFirstIndex+1);
   
   /* Ottengo il g predetto */
   gRef = Vec3(X, iFirstIndex+4);
   
   /* Calcolo la matrice RDelta derivante dalla predizione */
   Mat3x3 RDelta(MatR, gRef);
   
   /* Calcolo la R corrente in base alla predizione */ 
   RCurr = RRef = RDelta*RPrev; 
   
   /* Calcolo la Omega corrente in base alla predizione (gP "totale") */
   gPRef = Vec3(XP, iFirstIndex+4);
      
   /* Calcolo il nuovo Omega */
   WCurr = WRef = Mat3x3(MatG, gRef)*gPRef;  

   /* Resetto i parametri di rotazione e le derivate, g e gP */
   X.Put(iFirstIndex+4, Vec3(0.));
   XP.Put(iFirstIndex+4, Vec3(0.));   

   gCurr = gPCurr = Vec3(0.);
}

/* StructNode - end */


/* DynamicStructNode - begin */

DynamicStructNode::DynamicStructNode(unsigned int uL, 
				     const DofOwner* pDO,
				     const Vec3& X0, 
				     const Mat3x3& R0,
				     const Vec3& V0, 
				     const Vec3& W0,
				     doublereal dPosStiff,
				     doublereal dVelStiff, 
				     flag fOmRot,
				     flag fOut)
: StructNode(uL, pDO, X0, R0, V0, W0, dPosStiff, dVelStiff, fOmRot, fOut)
{
   NO_OP;
}

   
/* Distruttore (per ora e' banale) */
DynamicStructNode::~DynamicStructNode(void)
{
   NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type 
DynamicStructNode::GetStructNodeType(void) const
{
   return StructNode::DYNAMIC;
}

   
/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
 * al posto giusto */
integer 
DynamicStructNode::iGetFirstRowIndex(void) const
{
   return iGetFirstMomentumIndex(); 
}


void 
DynamicStructNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   /* Mettere a punto */
   integer iFirstIndex = iGetFirstIndex();
   X.Put(iFirstIndex+1, XCurr);
   X.Put(iFirstIndex+4, Vec3(0.));
   gRef = gCurr = gPRef = gPCurr = Vec3(0.);
   XP.Put(iFirstIndex+1, VCurr);
   XP.Put(iFirstIndex+4, WCurr);
}

/* DynamicStructNode - end */


/* StaticStructNode - begin */

/* Costruttore definitivo */
StaticStructNode::StaticStructNode(unsigned int uL, 
				   const DofOwner* pDO,
				   const Vec3& X0, 
				   const Mat3x3& R0,
				   const Vec3& V0, 
				   const Vec3& W0,
				   doublereal dPosStiff, 
				   doublereal dVelStiff, 
				   flag fOmRot,
				   flag fOut)
: StructNode(uL, pDO, X0, R0, V0, W0, dPosStiff, dVelStiff, fOmRot, fOut)
{
   NO_OP;
}


/* Distruttore (per ora e' banale) */
StaticStructNode::~StaticStructNode(void)
{
   NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type 
StaticStructNode::GetStructNodeType(void) const
{
   return StructNode::STATIC; 
}


void 
StaticStructNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   /* Mettere a punto */
   integer iFirstIndex = iGetFirstIndex();
   X.Put(iFirstIndex+1, XCurr);
   X.Put(iFirstIndex+4, Vec3(0.));   
   gRef = gCurr = gPRef = gPCurr = Vec3(0.);
   XP.Put(iFirstIndex+1, VCurr);
   XP.Put(iFirstIndex+4, WCurr);   
}

/* StaticStructNode - end */


/* ModalNode - begin */

ModalNode::ModalNode(unsigned int uL, 
				     const DofOwner* pDO,
				     const Vec3& X0, 
				     const Mat3x3& R0,
				     const Vec3& V0, 
				     const Vec3& W0,
				     doublereal dPosStiff,
				     doublereal dVelStiff, 
				     flag fOmRot,
				     flag fOut)
: DynamicStructNode(uL, pDO, X0, R0, V0, W0, dPosStiff, dVelStiff, fOmRot, fOut),
XPPCurr(0.),
WPCurr(0.)
{
   NO_OP;
}

   
/* Distruttore (per ora e' banale) */
ModalNode::~ModalNode(void)
{
   NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type 
ModalNode::GetStructNodeType(void) const
{
   return StructNode::MODAL;
}

   
/* Usato dalle forze astratte, dai bulk ecc., per assemblare le forze
 * al posto giusto */
integer 
ModalNode::iGetFirstRowIndex(void) const
{
   return iGetFirstMomentumIndex(); 
}


void 
ModalNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   /* Mettere a punto */
   integer iFirstIndex = iGetFirstIndex();
   X.Put(iFirstIndex+1, XCurr);
   X.Put(iFirstIndex+4, Vec3(0.));
   gRef = gCurr = gPRef = gPCurr = Vec3(0.);
   XP.Put(iFirstIndex+1, VCurr);
   XP.Put(iFirstIndex+4, WCurr);

  
}

/* Aggiorna dati in base alla soluzione */
void 
ModalNode::Update(const VectorHandler& X, const VectorHandler& XP)
{   
   StructNode::Update(X, XP);
   
   integer iFirstIndex = iGetFirstIndex();
   
   /* aggiorno XPP e WP (servono solo a modal.cc) */
   XPPCurr = Vec3(XP, iFirstIndex+7);
   WPCurr  = Vec3(XP, iFirstIndex+10);      
}

/* ModalNode - end */


/* DummyStructNode - begin */

/* Costruttore definitivo */
DummyStructNode::DummyStructNode(unsigned int uL,
				 const DofOwner* pDO,
				 const StructNode* pN)
: StructNode(uL, pDO, 0., 0., 0., 0., 0., 0., 0, flag(1)), pNode(pN)
{
   ASSERT(pNode != NULL);
}


/* Distruttore (per ora e' banale) */
DummyStructNode::~DummyStructNode(void)
{
   NO_OP;
}


/* Tipo di nodo strutturale */
StructNode::Type 
DummyStructNode::GetStructNodeType(void) const
{
   return StructNode::DUMMY;
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const 
doublereal& DummyStructNode::dGetDofValue(int iDof, int iOrder) const
{
   std::cerr << "DummyStructNode(" << GetLabel() << ") has no dofs" << std::endl;
   THROW(ErrGeneric());   
}


/* Restituisce il valore del dof iDof al passo precedente;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const 
doublereal& DummyStructNode::dGetDofValuePrev(int iDof, int iOrder) const
{
   std::cerr << "DummyStructNode(" << GetLabel() << ") has no dofs" << std::endl;
   THROW(ErrGeneric());   
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void 
DummyStructNode::SetDofValue(const doublereal& dValue,
			     unsigned int iDof, unsigned int iOrder)
{
   std::cerr << "DummyStructNode(" << GetLabel() << ") has no dofs" << std::endl;
   THROW(ErrGeneric());   
}

   
/* Aggiorna dati durante l'iterazione fittizia iniziale */
void 
DummyStructNode::DerivativesUpdate(const VectorHandler& X, const VectorHandler& XP)
{
   /* posso farlo perche' in genere i dummy nodes si limitano
    * a copiare i valori di altri nodi, quindi non alterano 
    * le variabili cinematiche */
   Update(X, XP);
}


/* Aggiorna dati in base alla soluzione durante l'assemblaggio iniziale */
void 
DummyStructNode::InitialUpdate(const VectorHandler& /* X */ )
{
   NO_OP;
}


/* Funzioni di inizializzazione, ereditate da DofOwnerOwner */
void 
DummyStructNode::SetInitialValue(VectorHandler& /* X */ ) const
{
   NO_OP;
}


void 
DummyStructNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   ((DummyStructNode*)this)->Update(X, XP);
}


/* Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator */
void 
DummyStructNode::BeforePredict(VectorHandler& /* X */ ,
			       VectorHandler& /* XP */ ,
			       VectorHandler& /* XPrev */ , 
			       VectorHandler& /* XPPrev */ ) const
{
   NO_OP;
}


void 
DummyStructNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
{
   Update(X, XP);
}

/* DummyStructNode - end */


/* OffsetDummyStructNode - begin */

/* Costruttore definitivo */
OffsetDummyStructNode::OffsetDummyStructNode(unsigned int uL,
					     const DofOwner* pDO,
					     const StructNode* pN,
					     const Vec3& f,
					     const Mat3x3& R)
: DummyStructNode(uL, pDO, pN), f(f), R(R)
{
   /* forzo la ricostruzione del nodo strutturale sottostante */
   __Update();
}


/* Distruttore (per ora e' banale) */
OffsetDummyStructNode::~OffsetDummyStructNode(void)
{
   NO_OP;
}


/* update - interno */
void OffsetDummyStructNode::__Update(void)
{
   RCurr = pNode->GetRCurr();
   XCurr = pNode->GetXCurr()+RCurr*f;
   WCurr = pNode->GetWCurr();
   VCurr = pNode->GetVCurr()+WCurr.Cross(RCurr*f);	
}


/* Tipo di nodo dummy */
DummyStructNode::Type 
OffsetDummyStructNode::GetDummyType(void) const
{
   return DummyStructNode::OFFSET;
}


/* Aggiorna dati in base alla soluzione */
void
OffsetDummyStructNode::Update(const VectorHandler& /* X */ ,
			      const VectorHandler& /* XP */ )
{
   __Update();
}

/* OffsetDummyStructNode - end */


/* RelFrameDummyStructNode - begin */

/* Costruttore definitivo */
RelFrameDummyStructNode::RelFrameDummyStructNode(unsigned int uL,
						 const DofOwner* pDO,
						 const StructNode* pN,
						 const StructNode* pNR,
						 const Vec3& fh,
						 const Mat3x3& Rh)
: DummyStructNode(uL, pDO, pN), pNodeRef(pNR), RhT(Rh.Transpose()), fhT(RhT*fh)
{
   ASSERT(pNodeRef != NULL);

   /*
    * Note: Rh is transposed from the beginning because it is 
    *       never used directly;
    *       fh is premultiplied by Rh.Transpose() for the same reason
    *
    * Formulas:
    *
    * R = RhT * RrT * Rn
    * X = RhT * RrT * (Xn - Xr)
    * W = RhT * RrT * (Wn - Wr)
    * V = RhT * RrT * (Vn - Vr - Wr x (Xn - Xr))
    *
    * by defining
    *
    * Rn = Rr * Rh * R
    * Xn = Xr + Rr * (fh + Rh * X)
    *
    * and differentiating with respect to time
    */
   
   /* forzo la ricostruzione del nodo strutturale sottostante */
   __Update();
}


/* Distruttore (per ora e' banale) */
RelFrameDummyStructNode::~RelFrameDummyStructNode(void)
{
   NO_OP;
}


/* update - interno */
void RelFrameDummyStructNode::__Update(void)
{
   Mat3x3 RrT(pNodeRef->GetRCurr().Transpose());
   Mat3x3 RT(RhT*RrT);
   Vec3 XRel(pNode->GetXCurr()-pNodeRef->GetXCurr());
   
   RCurr = RT*pNode->GetRCurr();
   XCurr = RT*XRel - fhT;
   WCurr = RT*(pNode->GetWCurr()-pNodeRef->GetWCurr());
     
   VCurr = RT*(pNode->GetVCurr()
	       -pNodeRef->GetVCurr()
	       -pNodeRef->GetWCurr().Cross(XRel));
}


/* Tipo di nodo dummy */
DummyStructNode::Type 
RelFrameDummyStructNode::GetDummyType(void) const
{
   return DummyStructNode::RELATIVEFRAME;
}


/* Aggiorna dati in base alla soluzione */
void 
RelFrameDummyStructNode::Update(const VectorHandler& /* X */ ,
				const VectorHandler& /* XP */ )
{
   __Update();
}

/* RelFrameDummyStructNode - end */


/* Legge un nodo strutturale */

Node* 
ReadStructNode(DataManager* pDM, 
	       MBDynParser& HP,
	       DofOwner* pDO,
	       unsigned int uLabel)
{
   const char sFuncName[] = "ReadStructNode()";
   DEBUGCOUT("Entering " << sFuncName << std::endl);
   
   const char* sKeyWords[] = {
      "static",
      "dynamic",
      "modal",
      "dummy",
      
      "offset",
      "relativeframe"   /* temporary */
   };
   
   /* enum delle parole chiave */
   enum KeyWords { 
      UNKNOWN = -1,
      
      STATIC = 0,
      DYNAMIC,
      MODAL,
      DUMMY,
      
      OFFSET,
      RELATIVEFRAME,
      
      LASTKEYWORD
   };
   
   /* tabella delle parole chiave */
   KeyTable K((int)LASTKEYWORD, sKeyWords);
   
   /* parser del blocco di controllo */
   HP.PutKeyTable(K);

   /* lettura dati specifici */
   KeyWords CurrType((KeyWords)HP.IsKeyWord());

   /*
    * explicit node type required; default is no longer "DYNAMIC"
    */
   if (CurrType == UNKNOWN) {
      std::cerr << "missing node type at line " << HP.GetLineData()
      	      << std::endl;
      THROW(ErrGeneric());
   }   
   
#ifdef DEBUG
   if (CurrType == STATIC) {      
      std::cout << "Static structural node" << std::endl;
   } else if (CurrType == DYNAMIC) {
      std::cout << "Dynamic structural node" << std::endl;
   } else if (CurrType == DUMMY) {
      std::cout << "Dummy structural node" << std::endl;
   } else if (CurrType == MODAL) {
      std::cout << "Modal node" << std::endl;
   } else {
      std::cout << "Unknown structural node" << std::endl;
   }   
#endif /* DEBUG */
   
   Node* pNd = NULL;
   if (CurrType == DUMMY) {
      unsigned int uNode = (unsigned int)HP.GetInt();
      
      DEBUGCOUT("Linked to Node " << uNode << std::endl);
      
      /* verifica di esistenza del nodo */
      StructNode* pNode;
      if ((pNode = pDM->pFindStructNode(uNode)) == NULL) {
	 std::cerr << std::endl << sFuncName
	   << " at line " << HP.GetLineData() 
	     << ": structural node " << uNode
	   << " not defined" << std::endl;    
	 THROW(DataManager::ErrGeneric());
      }
      
      KeyWords DummyType = KeyWords(HP.GetWord());
      switch (DummyType) {
       case OFFSET: {
	  ReferenceFrame RF(pNode);
	  Vec3 f(HP.GetPosRel(RF));
	  Mat3x3 R(HP.GetRotRel(RF));
	  
	  SAFENEWWITHCONSTRUCTOR(pNd, 
				 OffsetDummyStructNode,
				 OffsetDummyStructNode(uLabel, pDO, pNode, f, R));
	  break;
       }
	 
       case RELATIVEFRAME: {
	  unsigned int uNodeRef = (unsigned int)HP.GetInt();
	  
	  DEBUGCOUT("Linked to Node " << uNodeRef << std::endl);
	  
	  /* verifica di esistenza del nodo */
	  StructNode* pNodeRef;
	  if ((pNodeRef = pDM->pFindStructNode(uNodeRef)) == NULL) {
	     std::cerr << std::endl << sFuncName
	       << " at line " << HP.GetLineData() 
		 << ": structural node " << uNodeRef
	       << " not defined" << std::endl;    
	     THROW(DataManager::ErrGeneric());
	  }

	  ReferenceFrame RF(pNodeRef);

	  Vec3 fh(Zero3);
	  if (HP.IsKeyWord("position")) {
		  fh = HP.GetPosRel(RF);
	  }

	  Mat3x3 Rh(Eye3);
	  if (HP.IsKeyWord("orientation")) {
		  Rh = HP.GetRotRel(RF);
	  }

	  SAFENEWWITHCONSTRUCTOR(pNd, 
				 RelFrameDummyStructNode,
				 RelFrameDummyStructNode(uLabel, pDO,
					 pNode, pNodeRef, fh, Rh));
	  break;	  
       }
	 
       default: {
	  std::cerr << "unknown dummy node type at line " << HP.GetLineData() << std::endl;
	  THROW(ErrGeneric());
       }
      }
   } else {
      /* posizione (vettore di 3 elementi) */
      if (!HP.IsKeyWord("position")) {
	 pedantic_cerr("StructNode(" << uLabel 
		 << "): missing keyword \"position\" at line " 
		 << HP.GetLineData() << std::endl);
      }
      Vec3 X0(HP.GetPosAbs(AbsRefFrame));   
      DEBUGCOUT("X0 =" << std::endl << X0 << std::endl);
      
      /* sistema di riferimento (trucco dei due vettori) */  
      if (!HP.IsKeyWord("orientation")) {
	 pedantic_cerr("StructNode(" << uLabel 
		 << "): missing keyword \"orientation\" at line " 
		 << HP.GetLineData() << std::endl);
      }
      Mat3x3 R0(HP.GetRotAbs(AbsRefFrame));
      DEBUGCOUT("R0 =" << std::endl << R0 << std::endl);
      
      /* Velocita' iniziali (due vettori di 3 elementi, con la possibilita' 
       * di usare "null" per porli uguali a zero) */
      if (!HP.IsKeyWord("velocity")) {
	 pedantic_cerr("StructNode(" << uLabel 
		 << "): missing keyword \"velocity\" at line " 
		 << HP.GetLineData() << std::endl);
      }
      Vec3 XPrime0(HP.GetVelAbs(AbsRefFrame, X0));

      if (!HP.IsKeyWord("angular" "velocity")) {
	 pedantic_cerr("StructNode(" << uLabel 
		 << "): missing keyword \"angular velocity\" at line " 
		 << HP.GetLineData() << std::endl);
      }
      Vec3 Omega0(HP.GetOmeAbs(AbsRefFrame));
      DEBUGCOUT("Xprime0 =" << std::endl << XPrime0 << std::endl 
		<< "Omega0 =" << std::endl << Omega0 << std::endl);
      
      /* Rigidezza in assemblaggio diversa da quella di default
       * e flag di output */
      doublereal dPosStiff = pDM->dInitialPositionStiffness;
      doublereal dVelStiff = pDM->dInitialVelocityStiffness;
      flag fOmRot = pDM->fOmegaRotates;
      
      if (HP.fIsArg()) {
	 if (HP.IsKeyWord("assembly")) {
	    dPosStiff = HP.GetReal(dPosStiff);
	    dVelStiff = HP.GetReal(dVelStiff);
	    if (HP.IsKeyWord("yes")) {
	       fOmRot = 1;
	    } else if (HP.IsKeyWord("no")) {
	       fOmRot = 0;
	    } else {
	       silent_cerr("use keywords \"yes\" or \"no\"" << std::endl);
	       fOmRot = HP.GetInt(fOmRot);
	    }
	    
	    DEBUGCOUT("Initial position stiffness: " << dPosStiff << std::endl);
	    DEBUGCOUT("Initial velocity stiffness: " << dVelStiff << std::endl);
	    DEBUGCOUT("Omega rotates? : " << (fOmRot ? "yes" : "no") << std::endl);
	 }
      }
      
      pDO->SetScale(pDM->dReadScale(HP, DofOwner::STRUCTURALNODE));
      flag fOut = pDM->fReadOutput(HP, Node::STRUCTURAL);   
      
      /* Se non c'e' il punto e virgola finale */
      if (HP.fIsArg()) {
	 std::cerr << std::endl << sFuncName
	   << ": semicolon expected at line " << HP.GetLineData() << std::endl;      
	 THROW(DataManager::ErrGeneric());
      }   
      
      /* costruzione del nodo */
      if (CurrType == STATIC) {	
	 SAFENEWWITHCONSTRUCTOR(pNd, StaticStructNode,
				StaticStructNode(uLabel, pDO, 
						 X0, R0, 
						 XPrime0, Omega0,
						 dPosStiff, dVelStiff, 
						 fOmRot, fOut));
	 
      } else if(CurrType == DYNAMIC) {
	 SAFENEWWITHCONSTRUCTOR(pNd, DynamicStructNode,
				DynamicStructNode(uLabel, pDO, 
						  X0, R0, 
						  XPrime0, Omega0,
						  dPosStiff, dVelStiff, 
						  fOmRot, fOut));
	 
	 /* Incrementa il numero di elementi automatici dei nodi dinamici */
	 pDM->ElemData[Elem::AUTOMATICSTRUCTURAL].iNum++;
	 
      } else if(CurrType == MODAL) {
	 SAFENEWWITHCONSTRUCTOR(pNd, ModalNode,
				ModalNode(uLabel, pDO, 
					  X0, R0, 
					  XPrime0, Omega0,
					  dPosStiff, dVelStiff, 
					  fOmRot, fOut));
      }
   }
   
   ASSERT(pNd != NULL);

   return pNd;
} /* End of ReadStructNode() */

