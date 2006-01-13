/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

/* Tabella dei simboli; richiamata da table.cpp, parser.h */

#ifndef TABLE_H
#define TABLE_H

extern "C" {
#include <strings.h>
}

#include "except.h"
#include "mathtyp.h"

const int DEF_TABLE_SIZE = 127;

class Table {      
   friend std::ostream& operator << (std::ostream& out, Table& T);
   
 public:
   class ErrNameAlreadyDefined {};   
   
 private:
   int size;        /* dimensioni dell'array */
   VarList** v;     /* array di VarList */
   Int FindRow(const char* const name) const;
   
 public:
   Table(Int s = DEF_TABLE_SIZE, Int f = 0);
   virtual ~Table(void);
   Var* Put(const char* const name, const TypedValue& v);
   NamedValue* Put(NamedValue* p);
   NamedValue* Get(const char* const name) const;
};

extern std::ostream& operator << (std::ostream& out, Table& T);

#endif /* TABLE_H */

