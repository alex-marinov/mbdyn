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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <float.h>
#include <mathp.h>


/* helper per le funzioni built-in*/
Int 
asin_t(Real* d)
{
   if (d[0] > 1. || d[0] < -1.) {
      return 1;
   }
   return 0;
}

Int 
tan_t(Real* d)
{
   Real a = d[0]-int(d[0]/M_PI)*M_PI;
   if (fabs(fabs(a)-M_PI_2) < DBL_EPSILON) {
      return 1;
   }
   return 0;
}

Int 
acosh_t(Real* d)
{
   if (d[0] < 0.) {
      return 1;
   }
   return 0;
}

Int 
atanh_t(Real* d) {
   if (d[0] >= 1. || d[0] <= -1.) {
      return 1;
   }
   return 0;
}

Int 
log_t(Real* d) {
   if (d[0] <= 0.) {
      return 1;
   }
   return 0;
}

Real 
mp_rand(void)
{
   return Real(rand());
}

Real 
mp_rndm(void)
{
   return -1.+2.*(Real(rand())/Real(RAND_MAX));
}

Real 
mp_srnd(Real s)
{
   srand((unsigned int)s);
   return Real(0);
}

Real 
mp_prnt(Real s)
{
   cout << s << endl;
   return s;
}

Real 
mp_ctg(Real s)
{
   return 1./tan(s);
}

Int 
mp_ctg_t(Real* d)
{
   Real a = d[0]-int(d[0]/M_PI)*M_PI;
   if (fabs(a) < DBL_EPSILON) {
      return 1;
   }
   return 0;
}

Real 
mp_actg(Real s)
{
   return atan2(1., s);
}

Real 
mp_actg2(Real s1, Real s2)
{
   return atan2(s2, s1);
}

Real 
mp_ctgh(Real s) 
{
   return 1./tanh(s);
}

Int 
mp_ctgh_t(Real* d)
{
   if (fabs(d[0]) < DBL_EPSILON) {
      return 1;
   }
   return 0;
}

#ifdef __USE_XOPEN
Real 
mp_actgh(Real s)
{
   return atanh(1./s);
}
#endif 

Real 
mp_step(Real s)
{
   if (s > 0.) {
      return 1.;
   }
   return 0.;
}

Real 
mp_ramp(Real s)
{
   if (s > 0.) {
      return s;
   }
   return 0.;
}

Real 
mp_par(Real s)
{
   if (s > 0.) {
      return s*s;
   }
   return 0.;
}


/* tabella statica delle funzioni built-in */
static mathfuncs FuncTable[] = {
     { "asin",       1, { (Real (*)(void))((Real (*)(Real))asin)      }, asin_t,     ""   },
     { "acos",       1, { (Real (*)(void))((Real (*)(Real))acos)      }, asin_t,     ""   },
     { "atan",       1, { (Real (*)(void))((Real (*)(Real))atan)      }, NULL,       ""   },
     { "actan",      1, { (Real (*)(void))((Real (*)(Real))mp_actg)   }, NULL,       ""   },
     { "atan2",      2, { (Real (*)(void))((Real (*)(Real, Real))atan2)     }, NULL,       ""   },
     { "actan2",     2, { (Real (*)(void))((Real (*)(Real, Real))mp_actg2)  }, NULL,       ""   },
     
     { "cos",        1, { (Real (*)(void))((Real (*)(Real))cos)       }, NULL,       ""   },
     { "sin",        1, { (Real (*)(void))((Real (*)(Real))sin)       }, NULL,       ""   },
     { "tan",        1, { (Real (*)(void))((Real (*)(Real))tan)       }, tan_t,      ""   },
     { "ctan",       1, { (Real (*)(void))((Real (*)(Real))mp_ctg)    }, mp_ctg_t,   ""   },
     
     { "cosh",       1, { (Real (*)(void))((Real (*)(Real))cosh)      }, NULL,       ""   },
     { "sinh",       1, { (Real (*)(void))((Real (*)(Real))sinh)      }, NULL,       ""   },
     { "tanh",       1, { (Real (*)(void))((Real (*)(Real))tanh)      }, NULL,       ""   },
     { "ctanh",      1, { (Real (*)(void))((Real (*)(Real))mp_ctgh)   }, mp_ctgh_t,  ""   },
     
#ifdef __USE_XOPEN
     { "acosh",      1, { (Real (*)(void))((Real (*)(Real))acosh)     }, acosh_t,    ""   },
     { "asinh",      1, { (Real (*)(void))((Real (*)(Real))asinh)     }, NULL,       ""   },
     { "atanh",      1, { (Real (*)(void))((Real (*)(Real))atanh)     }, atanh_t,    ""   },
     { "actanh",     1, { (Real (*)(void))((Real (*)(Real))mp_actgh)  }, NULL,       ""   },     
#endif
   
     { "exp",        1, { (Real (*)(void))((Real (*)(Real))exp)       }, NULL,       ""   },
     { "log",        1, { (Real (*)(void))((Real (*)(Real))log)       }, log_t,      ""   },
     
     { "sqrt",       1, { (Real (*)(void))((Real (*)(Real))sqrt)      }, log_t,      ""   },
     
     { "abs",        1, { (Real (*)(void))((Real (*)(Real))fabs)      }, NULL,       ""   },
     { "floor",      1, { (Real (*)(void))((Real (*)(Real))floor)     }, NULL,       ""   },
     { "ceil",       1, { (Real (*)(void))((Real (*)(Real))ceil)      }, NULL,       ""   },
#ifdef __USE_XOPEN
     { "round",      1, { (Real (*)(void))((Real (*)(Real))rint)      }, NULL,       ""   },
#endif

     { "rand",       0, { (Real (*)(void))((Real (*)(void))mp_rand)   }, NULL,       ""   },
     { "random",     0, { (Real (*)(void))((Real (*)(void))mp_rndm)   }, NULL,       ""   },
     { "seed",       1, { (Real (*)(void))((Real (*)(Real))mp_srnd)   }, NULL,       ""   },

     { "step",       1, { (Real (*)(void))((Real (*)(Real))mp_step)   }, NULL,       ""   },
     { "ramp",       1, { (Real (*)(void))((Real (*)(Real))mp_ramp)   }, NULL,       ""   },
     { "par",        1, { (Real (*)(void))((Real (*)(Real))mp_par)    }, NULL,       ""   },
     
     { "print",      1, { (Real (*)(void))((Real (*)(Real))mp_prnt)   }, NULL,       ""   },
     
     /* add more as needed */
     
     { NULL,         0, { (Real (*)(void))NULL  }, NULL,     NULL }
};


/* tipi delle variabili */
struct typenames {
   const char* name;
   TypedValue::Type type;
};

static typenames TypeNames[] = {
     { "integer",  TypedValue::VAR_INT },
     { "real",     TypedValue::VAR_REAL },
     
     { NULL,   TypedValue::VAR_UNKNOWN }
};


TypedValue::TypedValue(const Int& i) : type(TypedValue::VAR_INT)
{
   v.i = i;
}

TypedValue::TypedValue(const Real& r) : type(TypedValue::VAR_REAL)
{
   v.r = r;
}

TypedValue::TypedValue(const TypedValue::Type t) : type(t)
{
   switch (type) {
    case TypedValue::VAR_INT:
      v.i = 0;
      break;
    case TypedValue::VAR_REAL:
      v.r = 0.;
      break;
    default:
      THROW(ErrUnknownType());
   }
}

TypedValue::TypedValue(const TypedValue& var) : type(var.type)
{
   switch (type) {
    case TypedValue::VAR_INT:
      v.i = var.v.i;
      break;
    case TypedValue::VAR_REAL:
      v.r = var.v.r;
      break;
    default:
      THROW(ErrUnknownType());
   }
}

const TypedValue& 
TypedValue::operator = (const TypedValue& var)
{
   switch ((type = var.type)) {
    case TypedValue::VAR_INT:
      Set(var.GetInt());
      break;
    case TypedValue::VAR_REAL:
      Set(var.GetReal());
      break;
    default:
      THROW(ErrUnknownType());
   }
   return *this;
}

TypedValue::Type 
TypedValue::GetType(void) const
{
   return type;
}

Int 
TypedValue::GetInt(void) const 
{
   switch (type) {
    case TypedValue::VAR_INT:
      return v.i;
    case TypedValue::VAR_REAL:
      return Int(v.r);
    default:
      THROW(ErrUnknownType());
   }
   return 0;
}

Real 
TypedValue::GetReal(void) const
{
   switch (type) {
    case TypedValue::VAR_INT:
      return Real(v.i);
    case TypedValue::VAR_REAL:
      return v.r;
    default:
      THROW(ErrUnknownType());
   }
   return 0.;
}

void 
TypedValue::SetType(TypedValue::Type t)
{
   type = t;
}  

const TypedValue& 
TypedValue::Set(const Int& i) 
{
   switch(GetType()) {
    case TypedValue::VAR_INT:
      v.i = i;
      break;
    case TypedValue::VAR_REAL:
      v.r = Real(i);
      break;
    default:
      THROW(ErrUnknownType());
   }   
   return *this;
}

const TypedValue& 
TypedValue::Set(const Real& r)
{
   switch(GetType()) {
    case TypedValue::VAR_INT:
      v.i = Int(r);
      break;
    case TypedValue::VAR_REAL:
      v.r = r;
      break;
    default:
      THROW(ErrUnknownType());
   }   
   return *this;
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator && (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator && (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() && v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator || (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator || (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() || v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator > (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator > (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() > v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator >= (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator >= (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() >= v.GetReal());
}

#if _G_HAVE_BOOL
bool 
TypedValue::operator == (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator == (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() == v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator <= (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator <= (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() <= v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator < (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator < (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() < v.GetReal());
}

#if _G_HAVE_BOOL 
bool 
TypedValue::operator != (const TypedValue& v) const
#else /* _G_HAVE_BOOL */
int 
TypedValue::operator != (const TypedValue& v) const
#endif /* _G_HAVE_BOOL */
{
   return (GetReal() != v.GetReal());
}

TypedValue 
TypedValue::operator + (const TypedValue& v) const
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return TypedValue(GetInt()+v.GetInt());
   }
   return TypedValue(GetReal()+v.GetReal());
}

TypedValue 
TypedValue::operator - (const TypedValue& v) const
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return TypedValue(GetInt()-v.GetInt());
   }
   return TypedValue(GetReal()-v.GetReal());
}

TypedValue 
TypedValue::operator * (const TypedValue& v) const
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return TypedValue(GetInt()*v.GetInt());
   }
   return TypedValue(GetReal()*v.GetReal());
}

TypedValue 
TypedValue::operator / (const TypedValue& v) const
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return TypedValue(GetInt()/v.GetInt());
   }
   return TypedValue(GetReal()/v.GetReal());
}

const TypedValue& 
TypedValue::operator += (const TypedValue& v)
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return Set(GetInt()+v.GetInt());
   }
   type = TypedValue::VAR_REAL;
   return Set(GetReal()+v.GetReal());
}

const TypedValue& 
TypedValue::operator -= (const TypedValue& v)
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return Set(GetInt()-v.GetInt());
   }
   Real d = GetReal()-v.GetReal();
   type = TypedValue::VAR_REAL;
   return Set(d);
}

const TypedValue& 
TypedValue::operator *= (const TypedValue& v)
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return Set(GetInt()*v.GetInt());
   }
   Real d = GetReal()*v.GetReal();
   type = TypedValue::VAR_REAL;
   return Set(d);
}

const TypedValue& 
TypedValue::operator /= (const TypedValue& v)
{
   if ((GetType() == TypedValue::VAR_INT)
       && (v.GetType() == TypedValue::VAR_INT)) {
      return Set(GetInt()/v.GetInt());
   }
   Real d = GetReal()/v.GetReal();
   type = TypedValue::VAR_REAL;
   return Set(d);
}

#if _G_HAVE_BOOL 
bool 
operator ! (const TypedValue& v)
#else /* _G_HAVE_BOOL */
int 
operator ! (const TypedValue& v)
#endif /* _G_HAVE_BOOL */
{
   return (!v.GetReal());
}

TypedValue 
operator - (const TypedValue& v) 
{
   switch (v.GetType()) {
    case TypedValue::VAR_INT:
      return TypedValue(-v.GetInt());
    case TypedValue::VAR_REAL:
      return TypedValue(-v.GetReal());
    default:
      THROW(TypedValue::ErrUnknownType());
   }
   return 0;
}

TypedValue 
operator + (const TypedValue& v) 
{
   return v;
}

ostream& 
operator << (ostream& out, const TypedValue& v) 
{
   switch (v.GetType()) {
    case TypedValue::VAR_INT:
      return out << v.GetInt();
    case TypedValue::VAR_REAL:
      return out << v.GetReal();
    default:
      THROW(TypedValue::ErrUnknownType());
   }
   return out;
}

/* classe per la memorizzazione delle variabili (modificare per usare Table) */

NamedValue::NamedValue(const char *const s)
: name(NULL) {
	AllocName(s);
};

NamedValue::~NamedValue(void)
{
	ASSERT(name != NULL);
	SAFEDELETEARR(name, MPmm);
}

void 
NamedValue::AllocName(const char* const s)
{
   ASSERT(s != NULL);
   SAFESTRDUP(name, s, MPmm);  
}

int 
NamedValue::IsVar(void) const
{
   return 0;
}

const char*
NamedValue::GetName(void) const
{
   ASSERT(name != NULL);
   return name;
}

Var::Var(const char* const s, const TypedValue& v)
: NamedValue(s), value(v) 
{
   NO_OP;
}

Var::Var(const char* const s, const Real& v)
: NamedValue(s), value(v) 
{
   NO_OP;
}

Var::Var(const char* const s, const Int& v)
: NamedValue(s), value(v) 
{
   NO_OP;
}

Var::~Var(void)
{
   NO_OP;
}

int
Var::IsVar(void) const
{
   return 1;
}

TypedValue::Type 
Var::GetType(void) const 
{
   return value.GetType();
}

TypedValue 
Var::GetVal(void) const 
{
   return value;
}

void 
Var::SetVal(const Real& v)
{
   value.Set(v);
}

void 
Var::SetVal(const Int& v)
{
   value.Set(v);
}

void 
Var::SetVal(const TypedValue& v)
{
   switch (GetType()) {
    case TypedValue::VAR_INT:
      value.Set(v.GetInt());
      break;    
    case TypedValue::VAR_REAL:
      value.Set(v.GetReal());
      break;
    default:
      THROW(TypedValue::ErrUnknownValue());
   }
}


MathParser::ErrGeneric::ErrGeneric(void) 
{
   NO_OP;
}

MathParser::ErrGeneric::ErrGeneric(MathParser* p, const char* const s) 
{
   cerr << s << " at line " << p->GetLineNumber() << endl;
}

MathParser::ErrGeneric::ErrGeneric(MathParser* p,
				   const char* const s1,
				   const char* const s2,
				   const char* const s3) 
{
   cerr << "MathParser - " << s1 << s2 << s3
     << " at line " << p->GetLineNumber() << endl;
}

/* gioca con table e stream di ingresso */
Table& 
MathParser::GetSymbolTable(void) const 
{
   return table;
}
   
void 
MathParser::PutSymbolTable(Table& T)
{
   (Table&)table = T;
}

int 
MathParser::GetLineNumber(void) const 
{
   ASSERT(in != NULL);
   return in->GetLineNumber();
}

MathParser::TokenList::TokenList(Token t) 
: t(t), value(Real(0)), name(NULL), next(NULL) 
{
   NO_OP;
}
      
MathParser::TokenList::TokenList(const char* const s) 
: t(NAME), value(Real(0)), name(NULL), next(NULL) 
{
   SAFESTRDUP(name, s, MPmm);
}
      
MathParser::TokenList::TokenList(const TypedValue& v) 
: t(NUM), value(v), name(NULL), next(NULL) 
{
   NO_OP;
}

MathParser::TokenList::~TokenList(void) 
{
   if (t == NAME) {
      SAFEDELETEARR(name, MPmm);
   }
}
   
void 
MathParser::TokenPush(enum Token t) 
{
   TokenList* p = NULL;
   if (t == NUM) {	 
      SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(value), MPmm);
   } else if (t == NAME) {      
      SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(namebuf), MPmm);
   } else {      
      SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(t), MPmm);
   }      
   p->next = tokenlist;
   tokenlist = p;
}

int 
MathParser::TokenPop(void) 
{
   if (tokenlist == NULL) {
      /* stack is empty */
      return 0;
   }
   currtoken = tokenlist->t;
   if (currtoken == NUM) {      
      value = tokenlist->value;
   } else if (currtoken == NAME) {
      ASSERT(tokenlist->name != NULL);
      strcpy(namebuf, tokenlist->name);
   }
   TokenList* p = tokenlist;
   tokenlist = tokenlist->next;
   SAFEDELETE(p, MPmm);
   
   return 1;
}

NamedValue* 
MathParser::GetVar(const char* const s) 
{
   for (VarList* p = varlist; p != NULL; p = p->next) {
      if (p->var != NULL && strcmp(s, p->var->GetName()) == 0) {
	 return p->var;
      }
   }
   return NULL;
}

Var* 
MathParser::NewVar(const char* const s, TypedValue::Type t, const Real& d) 
{
   if (GetVar(s) != NULL) {
      DEBUGCERR("var '" << s << "' already defined!" << endl);
      THROW(ErrGeneric(this, "var '", s, "' already defined!" ));
   }
   Var* v = NULL;
   switch (t) {
    case TypedValue::VAR_INT:
      SAFENEWWITHCONSTRUCTOR(v, Var, Var(s, Int(d)), MPmm);
      break;
    case TypedValue::VAR_REAL:
      SAFENEWWITHCONSTRUCTOR(v, Var, Var(s, d), MPmm);
      break;
    default: 
      THROW(TypedValue::ErrUnknownType());
   }
   VarList* p = NULL;
   SAFENEW(p, VarList, MPmm);
   p->var = v;
   p->next = varlist;
   varlist = p;
   return v;
}

mathfuncs* 
MathParser::GetFunc(const char* const s) const 
{
   for (int i = 0; FuncTable[i].fname != NULL; i++) {
      if (strcmp(s, FuncTable[i].fname) == 0) {
	 return &FuncTable[i];
      }
   }
   return NULL;
}

TypedValue::Type 
MathParser::GetType(const char* const s) const
{
   for (Int i = 0; TypeNames[i].name != NULL; i++) { 
      if (strcmp(s, TypeNames[i].name) == 0) {
	 return TypeNames[i].type;
      }
   }
   return TypedValue::VAR_UNKNOWN;
}

Int 
MathParser::IsType(const char* const s) const 
{
   for (int i = 0; TypeNames[i].name != NULL; i++) {
      if (strcmp(s, TypeNames[i].name) == 0) {
	 return 1;
      }	 
   }
   return 0;
}

Int 
MathParser::IsFunc(const char* const s) const 
{
   for (int i = 0; FuncTable[i].fname != NULL; i++) {
      if (strcmp(s, FuncTable[i].fname) == 0) {
	 return 1;
      }
   }
   return 0;
}

Int 
MathParser::IsKeyWord(const char* const s) const 
{
   if (IsType(s)) {
      return 1;
   }
   if (IsFunc(s)) {
      return 1;
   }
   return 0;
}

enum MathParser::Token 
MathParser::GetToken(void) {
   ASSERT(in != NULL);
   
   if (TokenPop()) { /* se lo trova! */
      return currtoken;
   }
      
   int c = 0;
   
start_parsing:
   
   /* skip spaces */
   while ((c = in->get()), isspace(c)) {
      NO_OP;
   };
   
   if (c == EOF || in->eof()) {
      return (currtoken = ENDOFFILE);
   }
   
   if (c == REMARK) {
      while ((c = in->get()) != '\n') {
	 NO_OP;
      }
      goto start_parsing;
   }
      
   /* number? */
   if (c == '.' || isdigit(c)) {
      static char s[256];
      int f = 0;
      int i = 0;
      
      s[i++] = c;
      
      if (c == '.') {
	 f = 1;
      }      
      while ((c = in->get()) == '.' || isdigit(c)) {
	 s[i++] = c;
	 if (c == '.') {
	    f = 1;
	 }
      }
      char e = tolower(c);
      if (e == 'e' || e == 'f' || e == 'd' || e == 'g') {
	 f = 1;
	 s[i++] = c;
	 if ((c = in->get()) == '-' || c == '+' || isdigit(c)) {
	    s[i++] = c;
	 } else {
	    return (currtoken = UNKNOWNTOKEN);
	 }
	 while (isdigit((c = in->get()))) {
	    s[i++] = c;
	 }
      }
      s[i] = '\0';
      in->putback(c);
      if (f == 0) {       
	 value.SetType(TypedValue::VAR_INT);
	 value.Set(Int(atoi(s)));
      } else {      
	 value.SetType(TypedValue::VAR_REAL);
	 value.Set(Real(atof(s)));
      }
      
      return (currtoken = NUM);
   }
   
   /* name? */
   if (isalpha(c) || c == '_') {
      int l = 0;	
      namebuf[l++] = char(c);	
      while ((c = in->get()), isalnum(c) || c == '_') {
	 namebuf[l++] = char(c);
	 if (l ==  namebuflen) {
	    IncNameBuf();
	 }
      }
      namebuf[l] = '\0';
      in->putback(c);
      return (currtoken = NAME);
   }
   
   switch (c) {
    case '^': 
      return (currtoken = EXP);
    case '*':
      return (currtoken = MULT);
    case '/':  
      if ((c = in->get()) == '*') {
	 for (;; c = in->get()) {
	    if (c == '*' && (c = in->get()) == '/') {
	       goto start_parsing;
	    }	    
	 }	 
      } else {
	 in->putback(c);	            
	 return (currtoken = DIV);
      }      
    case '-':
      return (currtoken = MINUS);
    case '+':
      return (currtoken = PLUS);
    case '>':
      if ((c = in->get()), c == '=') {
	 return (currtoken = GE);
      }
      in->putback(c);
      return (currtoken = GT);
    case '=':
      if ((c = in->get()), c == '=') {	    
	 return (currtoken = EQ);
      }	 
      in->putback(c);
      return (currtoken = ASSIGN);
    case '<':
      if ((c = in->get()), c == '=') {	 	    
	 return (currtoken = LE);
      }	 
      in->putback(c);
      return (currtoken = LT);
    case '!':
      if ((c = in->get()), c == '=') {
	 return (currtoken = NE);
      }
      in->putback(c);
      return (currtoken = NOT);
    case '&':
      if ((c = in->get()), c != '&') {
	 return (currtoken = UNKNOWNTOKEN);
      }
      return (currtoken = AND);
    case '|':
      if ((c = in->get()), c != '|') {	    
	 return (currtoken = UNKNOWNTOKEN);
      }	 
      return (currtoken = OR);
    case '~':
      if ((c = in->get()), c != '|') {	    
	 return (currtoken = UNKNOWNTOKEN);
      }	 
      return (currtoken = XOR); 	    
    case '(':
      return (currtoken = OBR);
    case ')':	 
      return (currtoken = CBR);
    case '[':
      return (currtoken = OPGIN);
    case ']':
      return (currtoken = CPGIN);
    case ';':
      return (currtoken = STMTSEP);
    case ',':	 
      return (currtoken = ARGSEP);
    default:
      return (currtoken = UNKNOWNTOKEN);
   }
}

void 
MathParser::IncNameBuf(void) 
{
   int oldlen = namebuflen;
   namebuflen *= 2;
   char* s = NULL;
   SAFENEWARR(s, char, (namebuflen+1), MPmm);
   strncpy(s, namebuf, oldlen);
   s[oldlen] = '\0';
   SAFEDELETEARR(namebuf, MPmm);
   namebuf = s;
}

TypedValue 
MathParser::logical(void) 
{
   TypedValue d = relational();
   while (1) {      
      switch (currtoken) {
       case AND:
	 GetToken();
	 d = d && relational();
	 break;
       case OR:
	 GetToken();
	 d = d || relational();
	 break;
       case XOR: {	  
	  GetToken();
	  TypedValue e = relational();
	  d = !(d && e) && (d || e);
	  break;
       }
       default:      
	 return d;
      }      
   }      
}

TypedValue 
MathParser::relational(void) 
{
   TypedValue d = binary();
   while (1) {      
      switch (currtoken) {	 
       case GT:
	 GetToken();
	 d = d > binary();
	 break;
       case GE:
	 GetToken();
	 d =  d >= binary();
	 break;
       case EQ:
	 GetToken();
	 d = d == binary();
	 break;
       case LE:
	 GetToken();
	 d = d <= binary();
	 break;
       case LT:
	 GetToken();
	 d = d < binary();
	 break;
       case NE:
	 GetToken();
	 d = d != binary();
	 break;
       default:	
	 return d;
      }
   }
}

TypedValue 
MathParser::binary(void) {
   TypedValue d = mult();
   while (1) {      
      switch (currtoken) {
       case PLUS:
	 GetToken();
	 d += mult();
	 break;
       case MINUS:
	 GetToken();
	 d -= mult();
	 break;
       default:
	 return d;
      }      
   }
}

TypedValue 
MathParser::mult(void)
{
   TypedValue d = power();
   while (1) {      
      switch(currtoken) {
       case MULT: 
	 GetToken();
	 d *= power();
	 break;
       case DIV: {	  
	  GetToken();
	  TypedValue e = power();
	  if (e == 0.) {
	     DEBUGCERR("divide by zero in mult()" << endl);	  
	     THROW(ErrGeneric(this, "divide by zero in mult()"));
	  }
	  d /= e;
	  break;
       }	 
       default: 
	 return d;
      }
   }
}

TypedValue 
MathParser::power(void) {
   TypedValue d = unary();
   while (1) {
      switch(currtoken) {
       case EXP:
	 GetToken();	
	 powerstack.Push(d);
	 d = unary();
	 break;
       default: {
	  TypedValue e(Real(0));
	  while (powerstack.Pop(e)) {	    
	     if (e < 0. && d <= 0.) {
		DEBUGCERR("can't compute " << e << '^' << d << " in power()" << endl);
		THROW(ErrGeneric(this, "invalid operands in power()"));
	     }	    
	     if (d.GetType() == TypedValue::VAR_INT
		 && e.GetType() == TypedValue::VAR_INT) {
		Int i = d.GetInt();
#if 1
		Int j = e.GetInt();
		Int r = j;
		if (i == 0) {
		   r = 1;
		} else if (i < 0) {
		   r = 0;
		} else {
		   for (Int k = i-1; k-- > 0; ) {
		      r *= j;
		   }
		}
		d = TypedValue(r);
#else /* 0 */
		d = TypedValue(Int(pow(e.GetInt(), i)));
#endif /* 0 */
	     } else {
		Real r = d.GetReal();
		d.SetType(TypedValue::VAR_REAL);
		d = TypedValue(Real(pow(e.GetReal(), r)));
	     }	     
	  }	  	  	  
	  
	  return d;
       }
      }
   }
}


TypedValue 
MathParser::unary(void) {
   switch(currtoken) {
    case MINUS:
      GetToken();
      return -expr();
    case PLUS:
      GetToken();
      return expr();
    case NOT:
      GetToken();
      return !expr();
    default:
      return expr();
   }
}

TypedValue 
MathParser::evalfunc(mathfuncs* f) 
{
   if (f->nargs == 0) {
      return TypedValue((*(f->f.f0))());
   }
   Real d[max_nargs];
   for (int i = 0; i < f->nargs-1; i++) {
      d[i] = stmtlist().GetReal();
      if (currtoken != ARGSEP) {	 
	 DEBUGCERR("argument separator expected in evalfunc()" << endl);
	 THROW(ErrGeneric(this, "argument separator expected"));
      }
      GetToken();
   }
   d[f->nargs-1] = stmtlist().GetReal();
   
   if (f->t != NULL) {
      if ((*(f->t))(d)) {
	 DEBUGCERR("error in function " << f->fname << " (msg: "
		   << f->errmsg << ')' << " in evalfunc()" << endl);
	 THROW(ErrGeneric(this, f->fname, ": error ", f->errmsg));
      }
   }
   
   switch (f->nargs) {
    case 1:
      return TypedValue((*(f->f.f1))(d[0]));
    case 2:
      return TypedValue((*(f->f.f2))(d[0], d[1]));
    default:
      DEBUGCERR("unsupported feature in evalfunc()" << endl);
      THROW(ErrGeneric(this, "unsupported number of arguments"));
   }
   return 0;
}

TypedValue 
MathParser::expr(void) 
{
   if (currtoken == NUM) {
      GetToken();
      return value;
   }
   
   if (currtoken == OBR) {
      GetToken();
      TypedValue d = stmtlist();
      if (currtoken != CBR) {
	 DEBUGCERR("closing parenthesis expected in expr()" << endl);
	 THROW(ErrGeneric(this, "closing parenthesis expected"));
      }
      GetToken();
      return d;
   }

   if (currtoken == OPGIN) {
       TypedValue d = readplugin();
       if (currtoken != CPGIN) {
	  DEBUGCERR("closing plugin expected in expr()" << endl);
	  THROW(ErrGeneric(this, "closing plugin expected"));
       }
       GetToken();
       return d;
   }
   
   if (currtoken == NAME) {
      GetToken();	 
      if (currtoken == OBR) {
	 /* in futuro ci potranno essere magari i dati strutturati */
	 if (!IsFunc(namebuf)) {
	    THROW(ErrGeneric(this, "user-defined functions not supported yet!"));
	 }
	 
	 mathfuncs* f = GetFunc(namebuf);
	 if (f == NULL) {
	    THROW(ErrGeneric(this, "function '", namebuf, "' not found"));
	 }
	 GetToken();
	 TypedValue d = evalfunc(f);
	 if (currtoken != CBR) {	   
	    DEBUGCERR("closing parenthesis expected after function '"
		      << f->fname << "' in expr()" << endl);
	    THROW(ErrGeneric(this, "closing parenthesis expected after function '",
			     f->fname, "' in expr()"));
	 }
	 GetToken();
	 return d;	    
      } else {
	 NamedValue* v = table.Get(namebuf);
	 if (v != NULL) {	       
	    return v->GetVal();
	 }

	 /* try to get a plugin */

      }
           
      DEBUGCERR("unknown name '" << namebuf << "' in expr()" << endl);
      THROW(ErrGeneric(this, "unknown name '", namebuf, "'"));
   }
   
   /* invalid expr */
   if (currtoken != ENDOFFILE) {
      DEBUGCERR("unknown token in expr()" << currtoken << endl);
      THROW(ErrGeneric(this, "unknown token"));
   }
   return TypedValue(0.);
}

TypedValue 
MathParser::stmt(void) 
{
   if (currtoken == NAME) {
      /* declaration? */
      if (IsType(namebuf)) {
	 TypedValue::Type type = GetType(namebuf);
	 ASSERT(type != TypedValue::VAR_UNKNOWN);
	 
	 if (GetToken() != NAME) {
	    cerr << "name expected after type in declaration in stmt()" << endl;
	    THROW(ErrGeneric(this, "name expected after type in declaration"));
	 }
	 
	 if (IsKeyWord(namebuf)) {
	    THROW(ErrGeneric(this, "name '", namebuf, "' is a keyword"));
	 }

	 /* with assign? */
	 if (GetToken() == ASSIGN) {
	    /* faccio una copia del nome! */
	    char* varname = NULL;
	    SAFESTRDUP(varname, namebuf, MPmm);
	    
	    GetToken();
	    TypedValue d = logical();

	    /* cerco la variabile */
	    NamedValue* v = table.Get(varname);
	    
	    if (v == NULL) {
	       /* Se non c'e' la inserisco */
	       v = table.Put(varname, TypedValue(type));
	       ((Var *)v)->SetVal(d);
	    } else {
	       /* altrimenti, se la posso ridefinire, mi limito 
		* ad assegnarle il nuovo valore */
	       if (redefine_vars) {
		  if (v->IsVar()) {
		     ((Var *)v)->SetVal(d);
		  } else {
		     cerr << "cannot redefine a non-var named value" << endl;
		     THROW(MathParser::ErrGeneric(this, "cannot redefine non-var named value '", v->GetName(), "'"));
		  }
	       } else {
		  /* altrimenti la reinserisco, cosi' da provocare l'errore
		   * di table */
		  v = table.Put(varname, TypedValue(type));
	       }
	    }
	    /* distruggo il temporaneo */
	    SAFEDELETEARR(varname, MPmm);
	    return v->GetVal();
	    
	 } else if (currtoken == STMTSEP) {
	    NamedValue* v = table.Get(namebuf);
	    if (v == NULL || !redefine_vars) {
	       /* se la var non esiste, la inserisco;
		* se invece esiste e non vale la ridefinizione,
		* tento di inserirla comunque, cosi' table da' errore */
	       v = table.Put(namebuf, TypedValue(type));
	    }

	    return v->GetVal();
	 }	    
      } else {
	 /* assignment? */
	 NamedValue* v = table.Get(namebuf);
	 if (v != NULL) {
	    if (GetToken() == ASSIGN) {
	       GetToken();
	       TypedValue d = logical();
	       if (v->IsVar()) {
	          ((Var *)v)->SetVal(d);
	       } else {
		  cerr << "cannot assign a non-val named value" << endl;
		  THROW(MathParser::ErrGeneric(this, "cannot assign non-val named value '", v->GetName(), "'"));
	       }
	       return v->GetVal();
	    } else {
	       TokenPush(currtoken);
	       currtoken = NAME;		 
	    }
	 } /* else could be a function or a variable */	  	
      }
   } 
   return logical();
}

TypedValue 
MathParser::readplugin(void)
{
   THROW(MathParser::ErrGeneric(this, "plugins not supported yet"));
}
   
TypedValue 
MathParser::stmtlist(void) 
{
   TypedValue d = stmt();
   if (currtoken == STMTSEP) {
      GetToken();
      return stmtlist();
   }
   return d;
}


const int default_namebuflen = 127;

MathParser::MathParser(const InputStream& strm, Table& t, int redefine_vars)
: PlugIns(NULL),
table(t),
redefine_vars(redefine_vars),
in((InputStream*)&strm), 
namebuf(NULL), 
namebuflen(default_namebuflen),
value(Real(0)), 
varlist(NULL), 
tokenlist(NULL), 
powerstack()
{
   DEBUGCOUTFNAME("MathParser::MathParser");
   
   SAFENEWARR(namebuf, char, (namebuflen+1), MPmm);
   
   time_t tm;
   time(&tm);
   srand(tm);
}

MathParser::MathParser(Table& t, int redefine_vars)
: PlugIns(NULL),
table(t), 
redefine_vars(redefine_vars),
in(NULL),
namebuf(NULL), 
namebuflen(default_namebuflen),
value(Real(0)), 
varlist(NULL), 
tokenlist(NULL), 
powerstack()
{
   DEBUGCOUTFNAME("MathParser::MathParser");

   SAFENEWARR(namebuf, char, (namebuflen+1), MPmm);
   
   time_t tm;
   time(&tm);
   srand(tm);
}

void 
MathParser::InsertSym(const char* const s, const Real& v)
{
   NamedValue* var = table.Put(s, TypedValue(v));
   if (var == NULL) {
      THROW(ErrGeneric(this, "error while adding real var '", s, ""));
   }
}

void 
MathParser::InsertSym(const char* const s, const Int& v) 
{
   NamedValue* var = table.Put(s, TypedValue(v));
   if (var == NULL) {  
      THROW(ErrGeneric(this, "error while adding integer var '", s, ""));
   }
}
   
MathParser::~MathParser(void) 
{
   DEBUGCOUTFNAME("MathParser::~MathParser");

   if (namebuf != NULL) {
      SAFEDELETEARR(namebuf, MPmm);
   }
}
     
Real 
MathParser::GetLastStmt(Real d, Token t)
{
   if (GetToken() == t) {     
      return d;
   }      
   while (currtoken != ENDOFFILE && currtoken != t) {
      d = stmtlist().GetReal();
   }   
   return d;
}

Real 
MathParser::GetLastStmt(const InputStream& strm, Real d, Token t) 
{
   const InputStream* p = in;
   in = (InputStream*)&strm;
   d =  GetLastStmt(d, t);
   in = (InputStream*)p;
   return d;
}

Real 
MathParser::Get(Real d)
{
   GetToken();
   d = stmt().GetReal();
   if (currtoken != STMTSEP && currtoken != ENDOFFILE) {     
      DEBUGCERR("stmtsep expected in Get()" << endl);
      THROW(ErrGeneric(this, "statement separator expected"));
   }
   return d;
}

Real 
MathParser::Get(const InputStream& strm, Real d)
{
   const InputStream* p = in;
   in = (InputStream*)&strm;
   GetToken();
   if (currtoken != STMTSEP && currtoken != ARGSEP) {
      d = stmt().GetReal();
   }
   if (currtoken == STMTSEP) {
      in->putback(';');
   } else if (currtoken == ARGSEP) {
      in->putback(',');
   } else {
      DEBUGCERR("separator expected in Get()" << endl);
      THROW(ErrGeneric(this, "separator expected"));
   }
   in = (InputStream*)p;
   return d;
}

void 
MathParser::GetForever(ostream& out, const char* const sep)
{
   do {
      out << Get(0.) << sep;	
   } while (currtoken == STMTSEP);
}

void 
MathParser::GetForever(const InputStream& strm, ostream& out, const char* const sep)
{
   const InputStream* p = in;
   in = (InputStream*)&strm;
   GetForever(out);
   in = (InputStream*)p;
}

int 
MathParser::RegisterPlugIn(const char *name, MathParser::PlugIn * (*constructor)(MathParser&, void *), void *arg)
{
   PlugInRegister *p = NULL;
   SAFENEW(p, PlugInRegister, MPmm);
   SAFESTRDUP(p->name, name, MPmm);
   p->constructor = constructor;
   p->arg = arg;
   p->next = PlugIns;
   PlugIns = p;
   
   return 0;
}

