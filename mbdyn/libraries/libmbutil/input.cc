/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

/* Input */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "input.h"

/* InputStream - begin */

/* Costruttore - inizializza il filtro con un reference ad un istream */
InputStream::InputStream(std::istream& in) 
: iStrm(in), uLineNumber(1) 
{
   NO_OP;
}

/* Distruttore banale */
InputStream::~InputStream(void) { 
   NO_OP;
}  

InputStream& operator >> (InputStream& in, int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, long int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, short int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, unsigned int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, unsigned long int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, unsigned short int& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, char& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, float& i)
{
   in.iStrm >> i;
   return in;
}

InputStream& operator >> (InputStream& in, double& i)
{
   in.iStrm >> i;
   return in;
}

#if 0
InputStream& operator >> (InputStream& in, long double& i)
{
   in.iStrm >> i;
   return in;
}
#endif

/* InputStream - end */

