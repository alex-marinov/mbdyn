/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "withlab.h"


/* WithLabel - begin */

WithLabel::WithLabel(unsigned int uL, const std::string& sN) 
: uLabel(uL), sName(sN)
{
	NO_OP;
}
	
WithLabel::~WithLabel(void)
{
	NO_OP;
}

void
WithLabel::PutLabel(unsigned int uL)
{
	uLabel = uL;
}
	
void
WithLabel::PutName(const std::string& sN)
{
	sName = sN;
}
   
unsigned int
WithLabel::GetLabel(void) const
{
	return uLabel; 
}

const std::string&
WithLabel::GetName(void) const
{
	return sName;
}
 
/* WithLabel - end */
