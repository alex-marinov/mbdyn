/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2005
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

#ifndef THERMALNODE_H
#define THERMALNODE_H

#include "node.h"

/* ThermalNode - begin */

/* Nodo termico, descrive fisicamente un nodo di una rete struttura termica.
 */

class ThermalNode : public ScalarDifferentialNode {
public:
	/* Costruttore */
	ThermalNode(unsigned int uL, const DofOwner* pDO, 
		doublereal dx, doublereal dxp, flag fOut);
   
	/* Distruttore (per ora e' banale) */
	virtual ~ThermalNode(void);

	/* Tipo di nodo */
	virtual Node::Type GetNodeType(void) const;
	
	/* Output */
	void Output(OutputHandler&OH) const;
};

/* ThermalNode - end */

#endif /* THERMALNODE_H */

