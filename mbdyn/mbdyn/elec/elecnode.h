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

#ifndef ELECNODE_H
#define ELECNODE_H

#include "node.h"

/* AbstractNode - begin */

/* Nodo astratto, consente l'accesso ad un grado di liberta' messo in comune
 * e privo di significato fisico proprio. Gli elementi di tipo ElectricBulk
 * usano i nodi elettrici. Numerosi elementi elettrici elementari sono 
 * associati a nodi elettrici (impedenze, amplificatori, generatori di 
 * tensione e di corrente, ecc.). Altri elementi, di tipo Electric, possono 
 * essere associati sia a nodi elettrici (tipicamente a coppie di nodi, in 
 * quanto usano differenze di tensione), oppure a nodi astratti,
 * ovvero oggetti che consentono un riferimento esplicito a Dof messi in
 * comune che non hanno un significato fisico preciso */

/* Numero di dof del tipo di nodo - usato anche dal DofManager (?) */
const int iAbstractNodeDofNumber = 1;

class AbstractNode : public ScalarDifferentialNode {
 public:   
   /* Costruttore definitivo (da mettere a punto) */
   AbstractNode(unsigned int uL, const DofOwner* pDO,
		doublereal dx, doublereal dxp, flag fOut)
     : ScalarDifferentialNode(uL, pDO, dx, dxp, fOut) { 
	NO_OP; 
     };
   
   /* Distruttore (per ora e' banale) */
   virtual ~AbstractNode(void) {
      NO_OP;
   };
   
   /* Tipo di nodo */
   virtual Node::Type GetNodeType(void) const {
      return Node::ABSTRACT;
   };
   
   /* Output del nodo */
   virtual void Output(OutputHandler& OH) const;
   
   virtual void AfterPredict(VectorHandler& X, VectorHandler& XP) {
      Update(X, XP);
   };
};

/* AbstractNode - end */


/* ElectricNode - begin */

/* Nodo elettrico, descrive fisicamente un nodo di una rete elettrica.
 * Viene usato con gli elementi ElectricBulk. Numerosi elementi elettrici 
 * elementari sono associati a questi nodi (impedenze, amplificatori, 
 * generatori di tensione e di corrente, ecc.). Altri elementi, di tipo
 * Electric, possono essere associati sia a nodi elettrici (tipicamente a
 * coppie, in quanto usano differenze di tensione), oppure a nodi astratti,
 * ovvero oggetti che consentono un riferimento esplicito a Dof messi in
 * comune che non hanno un significato fisico preciso */

/* Numero di dof del tipo di nodo - usato anche dal DofManager (?) */
const int iElectricNodeDofNumber = iAbstractNodeDofNumber;

class ElectricNode : public AbstractNode {
 public:
   /* Costruttore */
   ElectricNode(unsigned int uL, const DofOwner* pDO, 
		doublereal dx, doublereal dxp, flag fOut)
     : AbstractNode(uL, pDO, dx, dxp, fOut) { 
	NO_OP; 
     };
   
   /* Distruttore (per ora e' banale) */
   virtual ~ElectricNode(void) { 
      NO_OP; 
   };

   /* Tipo di nodo */
   virtual Node::Type GetNodeType(void) const { 
      return Node::ELECTRIC;
   };
};

/* ElectricNode - end */

#endif /* ELECNODE_H */

