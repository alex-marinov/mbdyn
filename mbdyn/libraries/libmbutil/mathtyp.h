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

#ifndef MATHTYP_H
#define MATHTYP_H

/* definizione dei tipi */
typedef double Real;
typedef int Int;

/* valori con tipo */
class TypedValue {
 public:
   class ErrUnknownType {};
   class ErrWrongType {};
   class ErrUnknownValue {};

   enum Type {
      VAR_UNKNOWN = -1,
	VAR_INT,
	VAR_REAL,
	VAR_LAST
   };
   
 protected:
   TypedValue::Type type;
   union {
      Int i;
      Real r;
   } v;
   
 public:
   TypedValue(const Int& i);
   TypedValue(const Real& r);
   TypedValue(const TypedValue::Type t);
   TypedValue(const TypedValue& var);
   
   const TypedValue& operator = (const TypedValue& var);

   TypedValue::Type GetType(void) const;
   Int GetInt(void) const;
   Real GetReal(void) const;

   void SetType(TypedValue::Type t);
   const TypedValue& Set(const Int& i);
   const TypedValue& Set(const Real& r);
      
#if _G_HAVE_BOOL 
   bool operator && (const TypedValue& v) const;
   bool operator || (const TypedValue& v) const;
   bool operator > (const TypedValue& v) const;
   bool operator >= (const TypedValue& v) const;
   bool operator == (const TypedValue& v) const;
   bool operator <= (const TypedValue& v) const;
   bool operator < (const TypedValue& v) const;
   bool operator != (const TypedValue& v) const;
#else // _G_HAVE_BOOL
   int operator && (const TypedValue& v) const;
   int operator || (const TypedValue& v) const;
   int operator > (const TypedValue& v) const;
   int operator >= (const TypedValue& v) const;
   int operator == (const TypedValue& v) const;
   int operator <= (const TypedValue& v) const;
   int operator < (const TypedValue& v) const;
   int operator != (const TypedValue& v) const;
#endif // _G_HAVE_BOOL 

   TypedValue operator + (const TypedValue& v) const;
   TypedValue operator - (const TypedValue& v) const;
   TypedValue operator * (const TypedValue& v) const;
   TypedValue operator / (const TypedValue& v) const;

   const TypedValue& operator += (const TypedValue& v);
   const TypedValue& operator -= (const TypedValue& v);
   const TypedValue& operator *= (const TypedValue& v);
   const TypedValue& operator /= (const TypedValue& v);
};

#if _G_HAVE_BOOL
extern bool operator ! (const TypedValue& v);
#else // _G_HAVE_BOOL
extern int operator ! (const TypedValue& v);
#endif // _G_HAVE_BOOL

extern TypedValue operator - (const TypedValue& v);
extern TypedValue operator + (const TypedValue& v);
extern ostream& operator << (ostream& out, const TypedValue& v);


/* classe per la memorizzazione delle variabili */
class NamedValue {
 private:
   char *name;

   void AllocName(const char *const s);

 public:
   NamedValue(const char *const s);
   virtual ~NamedValue(void);

   virtual int IsVar(void) const;

   const char *GetName(void) const;
   virtual TypedValue::Type GetType(void) const = 0;
   virtual TypedValue GetVal(void) const = 0;
};

class Var : public NamedValue {
 private:
   TypedValue value;

 public:
   Var(const char* const s, const TypedValue& v);
   Var(const char* const s, const Real& v);
   Var(const char* const s, const Int& v);
   ~Var(void);

   int IsVar(void) const;
   
   TypedValue::Type GetType(void) const;
   TypedValue GetVal(void) const;
   
   void SetVal(const Real& v);
   void SetVal(const Int& v);
   void SetVal(const TypedValue& v);
};

/* lista delle variabili, sara' sostituita da Table (che la usa per le sue liste) */
struct VarList {
   NamedValue* var;
   
   VarList* next;
};

#endif // MATHTYP_H
