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
#include <node.h>
#include <solman.h>


/* Node - begin */

/* Costruttore */
Node::Node(unsigned int uL, const DofOwner* pDO, flag fOut)
: WithLabel(uL), DofOwnerOwner(pDO), ToBeOutput(fOut) 
{
   NO_OP; 
}
 

/* Distruttore banale */
Node::~Node(void)
{
   NO_OP; 
}

   
flag Node::fIsValidIndex(unsigned int i) const
{
   if (i >= 1 && i <= this->iGetNumDof()) {
      return flag(1);
   }
   return flag(0);
}
 

/* esegue operazioni sui dof di proprieta' dell'elemento 
 * in particolare ritorna il tipo di Dof in base all'indice i. Di default
 * i Dof dei nodi sono assunti differenziali */   
#ifdef DEBUG
DofOrder::Order Node::SetDof(unsigned int i) const 
#else
DofOrder::Order Node::SetDof(unsigned int /* i */ ) const 
#endif     
{
   ASSERT(i < this->iGetNumDof());
   return DofOrder::DIFFERENTIAL; 
}


/* Ritorna gli indici di riga e colonna. Tipicamente sono gli stessi */
integer Node::iGetFirstRowIndex(void) const 
{
   return iGetFirstIndex(); 
}


integer Node::iGetFirstColIndex(void) const
{
   return iGetFirstIndex();
}


/* Output di default per nodi di cui non si desidera output */
void Node::Output(OutputHandler& /* OH */ ) const
{
   NO_OP;
}


/* Output di default per nodi di cui non si desidera output */
void
Node::Output(
		OutputHandler& /* OH */,
		const VectorHandler& /* X */ ,
		const VectorHandler& /* XP */
		) const
{
   NO_OP;
}


/* Setta i valori iniziali delle variabili (e fa altre cose) 
 * prima di iniziare l'integrazione */
void Node::SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const 
{
   NO_OP;
}


/* Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator */
void Node::BeforePredict(VectorHandler& /* X */ ,
			 VectorHandler& /* XP */ ,
			 VectorHandler& /* XPrev */ ,
			 VectorHandler& /* XPPrev */ ) const 
{
   NO_OP; 
}

/* Node - end */


/* ScalarNode - begin */

ScalarNode::ScalarNode(unsigned int uL, const DofOwner* pDO, flag fOut)
: Node(uL, pDO, fOut) 
{ 
   NO_OP; 
}
 

ScalarNode::~ScalarNode(void)
{
   NO_OP;
}

/* ScalarNode - end */


/* ScalarDifferentialNode - begin */

ScalarDifferentialNode::ScalarDifferentialNode(unsigned int uL, 
					       const DofOwner* pDO, 
					       const doublereal& dx, 
					       const doublereal& dxp, 
					       flag fOut)
: ScalarNode(uL, pDO, fOut), dX(dx), dXP(dxp)
{
   NO_OP;
}


ScalarDifferentialNode::~ScalarDifferentialNode(void)
{
   NO_OP;
}


/* esegue operazioni sui dof di proprieta' dell'elemento 
 * in particolare ritorna il tipo di Dof in base all'indice i. Di default
 * i Dof dei nodi sono assunti differenziali */   
#ifdef DEBUG
DofOrder::Order ScalarDifferentialNode::SetDof(unsigned int i) const
#else
DofOrder::Order ScalarDifferentialNode::SetDof(unsigned int /* i */ ) const
#endif     
{ 
   ASSERT(i < this->iGetNumDof());
   return DofOrder::DIFFERENTIAL; 
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal& 
#ifdef DEBUG
ScalarDifferentialNode::dGetDofValue(int iDof, int iOrder) const
#else
ScalarDifferentialNode::dGetDofValue(int /* iDof */ , int iOrder) const
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0 || iOrder == 1);
   if(iOrder == 0) {
      return dX;
   }
   return dXP;
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void 
#ifdef DEBUG
ScalarDifferentialNode::SetDofValue(const doublereal& dValue,
				    unsigned int iDof, 
				    unsigned int iOrder)
#else
ScalarDifferentialNode::SetDofValue(const doublereal& dValue,
				    unsigned int /* iDof */ ,
				    unsigned int iOrder)
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0 || iOrder == 1);
   if(iOrder == 0) {
      dX = dValue;
   } else if(iOrder == 1) {
      dXP = dValue;
   }
}

   
/* Funzioni "spurie": consentono l'accesso ai dati privati;
 * sono state definite perche' i nodi astratti sono usati nei
 * modi piu' strani e quindi puo' essere necessario l'accesso */
void ScalarDifferentialNode::SetX(const doublereal& d) 
{
   dX = d;
}

/* only for differential nodes!?! */
void ScalarDifferentialNode::SetXPrime(const doublereal& d)
{
   dXP = d;
};
 

void ScalarDifferentialNode::SetValue(VectorHandler& X, VectorHandler& XP) const
{
   integer iIndex = iGetFirstIndex();
   X.fPutCoef(iIndex+1, dX);
   XP.fPutCoef(iIndex+1, dXP);
}

/*
 * Each node should prepend its type
 */
std::ostream& 
ScalarDifferentialNode::Restart(std::ostream& out) const
{
	return out << ", value, " << dX 
		<< ", derivative, " << dXP << ";" << std::endl;
}

/* ScalarDifferentialNode - end */


/* ScalarAlgebraicNode - begin */

ScalarAlgebraicNode::ScalarAlgebraicNode(unsigned int uL, 
					 const DofOwner* pDO, 
					 doublereal dx, 
					 flag fOut)
: ScalarNode(uL, pDO, fOut), dX(dx) 
{
   NO_OP;
}


ScalarAlgebraicNode::~ScalarAlgebraicNode(void)
{
   NO_OP;
}


/* esegue operazioni sui dof di proprieta' dell'elemento 
 * in particolare ritorna il tipo di Dof in base all'indice i. Di default
 * i Dof dei nodi sono assunti differenziali */   
#ifdef DEBUG
DofOrder::Order ScalarAlgebraicNode::SetDof(unsigned int i) const
#else
DofOrder::Order ScalarAlgebraicNode::SetDof(unsigned int /* i */ ) const
#endif     
{ 
   ASSERT(i < this->iGetNumDof());
   return DofOrder::ALGEBRAIC; 
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal&
#ifdef DEBUG
ScalarAlgebraicNode::dGetDofValue(int iDof, int iOrder) const
#else
ScalarAlgebraicNode::dGetDofValue(int /* iDof */ , int /* iOrder */ ) const
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0);      
   return dX;	       
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void 
#ifdef DEBUG
ScalarAlgebraicNode::SetDofValue(const doublereal& dValue, 
				 unsigned int iDof, 
				 unsigned int iOrder)
#else
ScalarAlgebraicNode::SetDofValue(const doublereal& dValue, 
				 unsigned int /* iDof */ ,
				 unsigned int /* iOrder */ )
#endif
{
   ASSERT(iDof == 1);
   ASSERT(iOrder == 0);      
   dX = dValue;
}

   
/* Funzioni "spurie": consentono l'accesso ai dati privati;
 * sono state definite perche' i nodi astratti sono usati nei
 * modi piu' strani e quindi puo' essere necessario l'accesso */
void ScalarAlgebraicNode::SetX(const doublereal& d)
{
   dX = d;
}


/* only for differential nodes!?! */
void ScalarAlgebraicNode::SetXPrime(const doublereal& /* d */ ) 
{
   DEBUGCERR("Error, setting derivative from algebraic dof" << std::endl);      
   THROW(Node::ErrGeneric());
}


void ScalarAlgebraicNode::SetValue(VectorHandler& X, VectorHandler& /* XP */ ) const
{
   integer iIndex = iGetFirstIndex();
   X.fPutCoef(iIndex+1, dX);
}

/* ScalarAlgebraicNode - end */


/* ParameterNode - begin */

ParameterNode::ParameterNode(unsigned int uL, 
			     const DofOwner* pDO,
			     doublereal dx, 
			     flag fOut)
: ScalarAlgebraicNode(uL, pDO, dx, fOut) 
{
   NO_OP;
}


ParameterNode::~ParameterNode(void)
{
   NO_OP;
}


/* Tipo del nodo (usato solo per debug ecc.) */
Node::Type ParameterNode::GetNodeType(void) const
{
   return Node::PARAMETER;
}


/* Contributo del nodo al file di restart */
std::ostream& ParameterNode::Restart(std::ostream& out) const
{
   return out << "parameter: " << GetLabel() << ", " << dX << ';' << std::endl;
}


flag ParameterNode::fIsValidIndex(unsigned int i) const
{
   if(i == 0) {
      return flag(1);
   }
   return flag(0);
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal& 
ParameterNode::dGetDofValue(int /* iDof */, int /* iOrder */ ) const
{
   return dX;
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void ParameterNode::SetDofValue(const doublereal& dValue, 
				unsigned int /* iDof */ , 
				unsigned int /* iOrder */) {
   dX = dValue;
}


/* Output di default per nodi di cui non si desidera output */
void ParameterNode::Output(OutputHandler& /* OH */ ) const
{
   NO_OP; 
}


void ParameterNode::SetValue(VectorHandler& /* X */ , 
			     VectorHandler& /* XP */ ) const
{
   NO_OP; 
}


/* Aggiorna dati in base alla soluzione */
void ParameterNode::Update(const VectorHandler& /* XCurr */ ,
			   const VectorHandler& /* XPrimeCurr */ )
{
   NO_OP; 
}


void ParameterNode::AfterPredict(VectorHandler& /* X */ , 
				 VectorHandler& /* XP */ )
{
   NO_OP; 
}

/* ParameterNode - end */


/* Node2Scalar - begin */

NodeDof::NodeDof(void) { 
   NO_OP;
}


NodeDof::NodeDof(unsigned int u, int id, Node* p)
: uNode(u), iDofNumber(id), pNode(p) {
   NO_OP;
}


NodeDof::~NodeDof(void) 
{
   NO_OP;
}


Node2Scalar::Node2Scalar(const NodeDof& nd)
: ScalarNode(nd.uNode, nd.pNode->pGetDofOwner(), 0), ND(nd)
{
   NO_OP;
}


Node2Scalar::~Node2Scalar(void)
{
   NO_OP;
}


/* Tipo del nodo (usato solo per debug ecc.) */
Node::Type Node2Scalar::GetNodeType(void) const
{
   return ND.pNode->GetNodeType();
}


/* Contributo del nodo al file di restart */
std::ostream& Node2Scalar::Restart(std::ostream& out) const
{
   out << "# Warning not implemented yet " << std::endl; // adds a remark!
     return out;
}


flag Node2Scalar::fIsValidIndex(unsigned int i) const
{
   if (i == 1) {
      return flag(1);
   }
   return flag(0);
}


/* esegue operazioni sui dof di proprieta' dell'elemento 
 * in particolare ritorna il tipo di Dof in base all'indice i. Di default
 * i Dof dei nodi sono assunti differenziali */   
#ifdef DEBUG
DofOrder::Order Node2Scalar::SetDof(unsigned int i) const 
#else
DofOrder::Order Node2Scalar::SetDof(unsigned int /* i */ ) const 
#endif     
{
   ASSERT(i < this->iGetNumDof());
   return DofOrder::DIFFERENTIAL; 
}


/* Ritorna gli indici di riga e colonna. Tipicamente sono gli stessi */
integer Node2Scalar::iGetFirstRowIndex(void) const
{
   return ND.pNode->iGetFirstRowIndex()+ND.iDofNumber;
}


integer Node2Scalar::iGetFirstColIndex(void) const
{
   return ND.pNode->iGetFirstColIndex()+ND.iDofNumber;
}


/* Restituisce il valore del dof iDof;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
const doublereal& Node2Scalar::dGetDofValue(int iDof, int iOrder) const
{
   if (iDof != 1) {
      THROW(ErrGeneric());
   }
   return ND.pNode->dGetDofValue(ND.iDofNumber+1, iOrder);
}


/* Setta il valore del dof iDof a dValue;
 * se differenziale, iOrder puo' essere = 1 per la derivata */
void Node2Scalar::SetDofValue(const doublereal& dValue, 
			      unsigned int iDof, 
			      unsigned int iOrder)
{
   ASSERT(iDof == 1);
   if (iDof == 1) {
      ((Node*)ND.pNode)->SetDofValue(dValue, ND.iDofNumber+1, iOrder);
   } 
   THROW(ErrGeneric());    
}


/* Output di default per nodi di cui non si desidera output */
void Node2Scalar::Output(OutputHandler& /* OH */ ) const
{
#ifdef DEBUG
   std::cout << "Default Node::Output(), the " << psNodeNames[GetNodeType()] 
     << " doesn't have any output" << std::endl;
#endif      
}


void Node2Scalar::SetValue(VectorHandler& /* X */ , 
			   VectorHandler& /* XP */ ) const
{
   NO_OP;
}

/* Aggiorna dati in base alla soluzione */
void Node2Scalar::Update(const VectorHandler& /* XCurr */ , 
			 const VectorHandler& /* XPrimeCurr */ )
{
   NO_OP;
}


/* Elaborazione vettori e dati prima e dopo la predizione
 * per MultiStepIntegrator */
void Node2Scalar::BeforePredict(VectorHandler& /* X */ ,
				VectorHandler& /* XP */ ,
				VectorHandler& /* XPrev */ ,
				VectorHandler& /* XPPrev */ ) const
{
   NO_OP; 
}


void Node2Scalar::AfterPredict(VectorHandler& /* X */ , 
			       VectorHandler& /* XP */ )
{
   NO_OP;
}

   
/* Funzioni "spurie": consentono l'accesso ai dati privati;
 * sono state definite perche' i nodi astratti sono usati nei
 * modi piu' strani e quindi puo' essere necessario l'accesso */
void Node2Scalar::SetX(const doublereal& d)
{
   SetDofValue(d, 1, 0);
}


/* only for differential nodes!?! */
void Node2Scalar::SetXPrime(const doublereal& d)
{
   SetDofValue(d, 1, 1);
}

/* Node2Scalar - end */


/* ScalarDof - begin */

ScalarDof::ScalarDof(void) 
: pNode(NULL), iOrder(0) 
{ 
   NO_OP; 
}


ScalarDof::ScalarDof(ScalarNode* p, int i) 
: pNode(p), iOrder(i)
{
   NO_OP;
}


ScalarDof::~ScalarDof(void)
{
   
   NO_OP;
}


doublereal ScalarDof::dGetValue(void) const
{
   return pNode->dGetDofValue(1, iOrder);
}

/* ScalarDof - end */
