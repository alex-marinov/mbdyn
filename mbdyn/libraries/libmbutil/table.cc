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

/* Codice della tabella dei simboli 
 * dell'analizzatore di espressioni alfanumeriche */

#include <mbconfig.h>

#include <table.h>

#include <myassert.h>
#include <mynewmem.h>

extern "C" {
#include <mymath.h>
}

enum {
     TABLE_INT,
     TABLE_REAL
};

#define GCC_BUG 1

struct tmp_sym {
   char* name;
#if GCC_BUG != 1
   TypedValue val;
#else /* GCC_BUG == 1 */
   int type;
   Int i_val;
   Real r_val;
#endif /* GCC_BUG == 1 */
};

tmp_sym consts[] = {
#if GCC_BUG != 1
     { "e",         Real(M_E)        }, /* */
     { "pi",        Real(M_PI)       }, /* */
     { "RAND_MAX",  Int(RAND_MAX)    }, /* 2147483647 */
     { "in2m",      Real(.0254)      }, /* conversion inches -> meters */
     { "m2in",      Real(1./.0254)   }, /* conversion meters -> inches */
     { "in2mm",     Real(25.4)       }, /* conversion inches -> millimeters */
     { "mm2in",     Real(1./25.4)    }, /* conversion millimeters -> inches */
     { "ft2m",      Real(.3048)      }, /* conversion feet -> meters */
     { "m2ft",      Real(1./.3048)   }, /* conversion meters -> feet */
     { "lb2kg",     Real(.4535)      }, /* conversion pounds -> kg */
     { "kg2lb",     Real(1./.4535)   }, /* conversion kg -> pounds */
   
   /* add as needed ... */

     { NULL,        Int(0)           }  /* terminale */
#else /* GCC_BUG == 1 */
     { "e",         TABLE_REAL, Int(0),        Real(M_E)      }, /* */
     { "pi",        TABLE_REAL, Int(0),        Real(M_PI)     }, /* */
     { "RAND_MAX",  TABLE_INT,  Int(RAND_MAX), Real(0.)       }, /* 2147483647 */
     { "in2m",      TABLE_REAL, Int(0),        Real(.0254)    }, /* conversion inches -> meters */
     { "m2in",      TABLE_REAL, Int(0),        Real(1./.0254) }, /* conversion meters -> inches */
     { "in2mm",     TABLE_REAL, Int(0),        Real(25.4)     }, /* conversion inches -> millimeters */
     { "mm2in",     TABLE_REAL, Int(0),        Real(1./25.4)  }, /* conversion millimeters -> inches */
     { "ft2m",      TABLE_REAL, Int(0),        Real(.3048)    }, /* conversion feet -> meters */
     { "m2ft",      TABLE_REAL, Int(0),        Real(1./.3048) }, /* conversion meters -> feet */
     { "lb2kg",     TABLE_REAL, Int(0),        Real(.4535)    }, /* conversion pounds -> kg */
     { "kg2lb",     TABLE_REAL, Int(0),        Real(1./.4535) }, /* conversion kg -> pounds */
   
   /* add as needed ... */

     { NULL,        TABLE_INT,  Int(0),        Real(0.)       }  /* terminale */
#endif /* GCC_BUG == 1 */
};


Table::Table(Int s, Int f) : size(s), v(NULL)
{
   static char func_name[] = "Table::Table";

   ASSERT(size > 0);
   DEBUGCOUT("Table::Table" << endl);
   
   SAFENEWARR(v, VarList*, size, MPmm);
   
   for (Int i = size; i-- > 0; ) {
      v[i] = NULL; 
   }
   
   if (!f) {
      return;
   }
   
   tmp_sym* p = consts;
   while (p->name != NULL) {
#if GCC_BUG != 1      
      NamedValue* n = Put(p->name, p->val);
#else /* GCC_BUG == 1 */      
      NamedValue* n = Put(p->name, (p->type == TABLE_REAL ? p->r_val : p->i_val));
#endif /* GCC_BUG == 1 */
      if (n == NULL) { 
	 cerr << func_name << ": unable to insert " << p->name << endl;
	 THROW(ErrGeneric());
      }
      p++;
   }
}

Table::~Table(void)
{
   if (v != NULL) {
      for (Int i = size; i-- > 0; ) {
	 VarList* pn = v[i];
	 while (pn != NULL) {
	    VarList* n = pn;
	    pn = n->next;
	    SAFEDELETE(n->var, MPmm);
	    SAFEDELETE(n, MPmm);
	 }
      }
      SAFEDELETEARR(v, MPmm);
   }
}

Int 
Table::FindRow(const char* const name) const
{
   int ii = 0;
   int i = 0; 
   while (name[i] != '\0') {
      ii = (ii << 1) ^ name[i++]; 
   }
   if (ii < 0) { 
      ii = -ii; 
   }
   
   ASSERT(size > 0);
   return (ii%size);
}

Var * 
Table::Put(const char* const name, const TypedValue& x)
{
   static char func_name[] = "Table::Put()";
   NamedValue* pVar = Get(name);
   if (pVar != NULL) {
      cerr << endl << func_name << ": name \"" << name 
	<< "\" already defined" << endl;
      THROW(Table::ErrNameAlreadyDefined());
   }
   
   int ii = FindRow(name);
   pVar = NULL;
   SAFENEWWITHCONSTRUCTOR(pVar, Var, Var(name, x), MPmm);
   VarList* pList = NULL;
   SAFENEW(pList, VarList, MPmm);
   pList->var = pVar;
   pList->next = v[ii];
   v[ii] = pList;
   return (Var *)pVar;
}

NamedValue *
Table::Put(NamedValue *p)
{
   static char func_name[] = "Table::Put()";
   NamedValue* pVar = Get(p->GetName());
   if (pVar != NULL) {
      cerr << endl << func_name << ": name \"" << p->GetName()
	<< "\" already defined" << endl;
      THROW(Table::ErrNameAlreadyDefined());
   }

   int ii = FindRow(p->GetName());
   pVar = NULL;
   VarList* pList = NULL;
   SAFENEW(pList, VarList, MPmm);
   pList->var = p;
   pList->next = v[ii];
   v[ii] = pList;
   return pList->var;
}

NamedValue* 
Table::Get(const char* const name) const
{
   int ii = FindRow(name);
   VarList* n = v[ii];
   
   
   while (n != NULL) {
      if (strcmp(n->var->GetName(), name) == 0) {
	 return n->var;
      }
      n = n->next;
   }
   return NULL;
}


ostream& 
operator << (ostream& out, Table& T)
{  
   for (Int i = T.size; i-- > 0; ) {      
      VarList* pn = T.v[i];
      while (pn != NULL) {
	 out << "  " << pn->var->GetName() << '<';
	 switch (pn->var->GetType()) {
	  case TypedValue::VAR_INT:
	    out << "int> = " << pn->var->GetVal().GetInt() << endl;
	    break;
	  case TypedValue::VAR_REAL:
	    out << "real> = " << pn->var->GetVal().GetReal() << endl;
	    break;
	  default:
	    out << "unknown type!>" << endl;
	    break;
	 }	 
	     
	 pn = pn->next;  
      }   
   }
   return out;
}
