/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

#include <ac/float.h>
#include <mathp.h>


/* helper per le funzioni built-in */
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
	Real a = d[0] - int(d[0]/M_PI)*M_PI;
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
atanh_t(Real* d)
{
	if (d[0] >= 1. || d[0] <= -1.) {
		return 1;
	}
	return 0;
}

Int
log_t(Real* d)
{
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
	return -1. + 2.*(Real(rand())/Real(RAND_MAX));
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
	silent_cout(s << std::endl);
	return s;
}

Real
mp_stop(Real s, Real v)
{
	if (s != 0.) {
		if (int(v) == 0) {
			silent_cout("mp_stop(SUCCESS)" << std::endl);
			throw NoErr();
		} else {
			silent_cout("mp_stop(FAILURE)" << std::endl);
			throw ErrGeneric();
		}
	}
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
	Real a = d[0] - int(d[0]/M_PI)*M_PI;
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

Real
mp_sign(Real s)
{
	return copysign(1., s);
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
mp_sramp(Real s, Real s_max)
{
	if (s <= 0.) {
		return 0.;
	}
	if (s >= s_max) {
		return s_max;
	}
	return s;
}

Real
mp_par(Real s)
{
	if (s > 0.) {
		return s*s;
	}
	return 0.;
}

Real max_t (Real a, Real b)
{
	return std::max(a, b);
}

Real min_t (Real a, Real b)
{
	return std::min(a, b);
}

/* tabella statica delle funzioni built-in */
static MathFunc_t FuncTable[] = {
	{ "asin",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)asin) },		asin_t,		"" },
	{ "acos",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)acos) },		asin_t,		"" },
	{ "atan",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)atan) },		0,		"" },
	{ "actan",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_actg) },	0,		"" },
	{ "atan2",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)atan2) },	0,		"" },
	{ "actan2",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)mp_actg2) },	0,		"" },

	{ "cos",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)cos) },		0,		"" },
	{ "sin",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)sin) },		0,		"" },
	{ "tan",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)tan) },		tan_t,		"" },
	{ "ctan",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_ctg) },	mp_ctg_t,	"" },

	{ "cosh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)cosh) },		0,		"" },
	{ "sinh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)sinh) },		0,		"" },
	{ "tanh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)tanh) },		0,		"" },
	{ "ctanh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_ctgh) },	mp_ctgh_t,	"" },

#ifdef __USE_XOPEN
	{ "acosh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)acosh) },	acosh_t,	"" },
	{ "asinh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)asinh) },	0,		"" },
	{ "atanh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)atanh) },	atanh_t,	"" },
	{ "actanh",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_actgh) },	0,		"" },
#endif /* __USE_XOPEN */

	{ "exp",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)exp) },		0,		"" },
	{ "log",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)log) },		log_t,		"" },
	{ "log10",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)log10) },	log_t,		"" },

	{ "sqrt",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)sqrt) },		log_t,		"" },

	{ "abs",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)fabs) },		0,		"" },
	{ "sign",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_sign) },	0,		"" },
	{ "copysign",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)copysign) },	0,		"" },
	{ "max",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)max_t) },	0,		"" },
	{ "min",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)min_t) },	0,		"" },
	{ "floor",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)floor) },	0,		"" },
	{ "ceil",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)ceil) },		0,		"" },

#ifdef __USE_XOPEN
	{ "round",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)rint) },		0,		"" },
#endif /* __USE_XOPEN */

	{ "rand",	0,	{ (MathFunc_0args_t)((MathFunc_0args_t)mp_rand) },	0,		"" },
	{ "random",	0,	{ (MathFunc_0args_t)((MathFunc_0args_t)mp_rndm) },	0,		"" },
	{ "seed",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_srnd) },	0,		"" },

	{ "step",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_step) },	0,		"" },
	{ "ramp",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_ramp) },	0,		"" },
	{ "sramp",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)mp_sramp) },	0,		"" },
	{ "par",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_par) },	0,		"" },

	{ "print",	1,	{ (MathFunc_0args_t)((MathFunc_1args_t)mp_prnt) },	0,		"" },
	{ "stop",	2,	{ (MathFunc_0args_t)((MathFunc_2args_t)mp_stop) },	0,		"" },

     /* add more as needed */
	{ 0,		0,	{ (MathFunc_0args_t)0 },				0,		0 }
};


/* tipi delle variabili */
struct TypeName_t {
	const char* name;
	TypedValue::Type type;
};

static TypeName_t TypeNames[] = {
	{ "integer",	TypedValue::VAR_INT },
	{ "real",	TypedValue::VAR_REAL },

	{ NULL,		TypedValue::VAR_UNKNOWN }
};

struct typemodifiernames {
	const char* name;
	TypedValue::TypeModifier type;
};

static typemodifiernames TypeModifierNames[] = {
	{ "const",	TypedValue::MOD_CONST },
	{ NULL,		TypedValue::MOD_UNKNOWN }
};


TypedValue::TypedValue(const Int& i, bool isConst)
: type(TypedValue::VAR_INT), bConst(isConst)
{
	v.i = i;
}

TypedValue::TypedValue(const Real& r, bool isConst)
: type(TypedValue::VAR_REAL), bConst(isConst)
{
	v.r = r;
}

TypedValue::TypedValue(const TypedValue::Type t, bool isConst)
: type(t), bConst(isConst)
{
	switch (type) {
	case TypedValue::VAR_INT:
		v.i = 0;
		break;
	case TypedValue::VAR_REAL:
		v.r = 0.;
		break;
	default:
		throw ErrUnknownType();
	}
}

TypedValue::TypedValue(const TypedValue& var)
: type(var.type), bConst(var.bConst)
{
	switch (type) {
	case TypedValue::VAR_INT:
		v.i = var.v.i;
		break;
	case TypedValue::VAR_REAL:
		v.r = var.v.r;
		break;
	default:
		throw ErrUnknownType();
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
		throw ErrUnknownType();
	}
	bConst = var.Const();
	return *this;
}

TypedValue::Type
TypedValue::GetType(void) const
{
	return type;
}

bool
TypedValue::Const(void) const
{
	return bConst;
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
		throw ErrUnknownType();
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
		throw ErrUnknownType();
	}
	return 0.;
}

void
TypedValue::SetType(TypedValue::Type t, bool isConst)
{
	type = t;
	bConst = isConst;
}

void
TypedValue::SetConst(bool isConst)
{
	bConst = isConst;
}

const TypedValue&
TypedValue::Set(const Int& i)
{
	switch (GetType()) {
	case TypedValue::VAR_INT:
		v.i = i;
		break;
	case TypedValue::VAR_REAL:
		v.r = Real(i);
		break;
	default:
		throw ErrUnknownType();
	}
	return *this;
}

const TypedValue&
TypedValue::Set(const Real& r)
{
	switch (GetType()) {
	case TypedValue::VAR_INT:
		v.i = Int(r);
		break;
	case TypedValue::VAR_REAL:
		v.r = r;
		break;
	default:
		throw ErrUnknownType();
	}
	return *this;
}

bool
TypedValue::operator && (const TypedValue& v) const
{
	return (GetReal() && v.GetReal());
}

bool
TypedValue::operator || (const TypedValue& v) const
{
	return (GetReal() || v.GetReal());
}

bool
TypedValue::operator > (const TypedValue& v) const
{
	return (GetReal() > v.GetReal());
}

bool
TypedValue::operator >= (const TypedValue& v) const
{
	return (GetReal() >= v.GetReal());
}

bool
TypedValue::operator == (const TypedValue& v) const
{
	return (GetReal() == v.GetReal());
}

bool
TypedValue::operator <= (const TypedValue& v) const
{
	return (GetReal() <= v.GetReal());
}

bool
TypedValue::operator < (const TypedValue& v) const
{
	return (GetReal() < v.GetReal());
}

bool
TypedValue::operator != (const TypedValue& v) const
{
	return (GetReal() != v.GetReal());
}

TypedValue
TypedValue::operator + (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return TypedValue(GetInt() + v.GetInt());
	}
	return TypedValue(GetReal() + v.GetReal());
}

TypedValue
TypedValue::operator - (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return TypedValue(GetInt() - v.GetInt());
	}
	return TypedValue(GetReal() - v.GetReal());
}

TypedValue
TypedValue::operator * (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return TypedValue(GetInt()*v.GetInt());
	}
	return TypedValue(GetReal()*v.GetReal());
}

TypedValue
TypedValue::operator / (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return TypedValue(GetInt()/v.GetInt());
	}
	return TypedValue(GetReal()/v.GetReal());
}

const TypedValue&
TypedValue::operator += (const TypedValue& v)
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return Set(GetInt() + v.GetInt());
	}
	Real d = GetReal() + v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

const TypedValue&
TypedValue::operator -= (const TypedValue& v)
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return Set(GetInt() - v.GetInt());
	}
	Real d = GetReal() - v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

const TypedValue&
TypedValue::operator *= (const TypedValue& v)
{
	if ((GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_INT))
	{
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
			&& (v.GetType() == TypedValue::VAR_INT))
	{
		return Set(GetInt()/v.GetInt());
	}
	Real d = GetReal()/v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

bool
operator ! (const TypedValue& v)
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
		throw TypedValue::ErrUnknownType();
	}
	return 0;
}

TypedValue
operator + (const TypedValue& v)
{
	return v;
}

std::ostream&
operator << (std::ostream& out, const TypedValue& v)
{
	switch (v.GetType()) {
	case TypedValue::VAR_INT:
		return out << v.GetInt();
	case TypedValue::VAR_REAL:
		return out << v.GetReal();
	default:
		throw TypedValue::ErrUnknownType();
	}
	return out;
}

/* classe per la memorizzazione delle variabili */

NamedValue::NamedValue(const char *const s)
: name(NULL)
{
	AllocName(s);
}

NamedValue::~NamedValue(void)
{
	ASSERT(name != NULL);
	SAFEDELETEARR(name);
}

void
NamedValue::AllocName(const char* const s)
{
	ASSERT(s != NULL);
	SAFESTRDUP(name, s);
}

bool
NamedValue::IsVar(void) const
{
	return false;
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

bool
Var::IsVar(void) const
{
	return true;
}

TypedValue::Type
Var::GetType(void) const
{
	return value.GetType();
}

bool
Var::Const(void) const
{
	return value.Const();
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
		throw TypedValue::ErrUnknownValue();
	}
}

void
MathParser::trim_arg(char *const s)
{
	int i, l;

	for (i = 0; isspace(s[i]); ++i) {
		NO_OP;
	}

	l = strlen(s+i);
	if (i > 0) {
		memmove(s, s+i, l+1);
	}

	for (i = l-1; isspace(s[i]); --i) {
		NO_OP;
	}
	s[i+1] = '\0';
}

MathParser::PlugInVar::PlugInVar(const char *const s, MathParser::PlugIn *p)
: NamedValue(s), pgin(p)
{
	NO_OP;
}

MathParser::PlugInVar::~PlugInVar(void)
{
	if (pgin) {
		SAFEDELETE(pgin);
	}
}

TypedValue::Type
MathParser::PlugInVar::GetType(void) const
{
	return pgin->GetType();
}

bool
MathParser::PlugInVar::Const(void) const
{
	return true;
}

TypedValue
MathParser::PlugInVar::GetVal(void) const
{
	return pgin->GetVal();
}

MathParser::ErrGeneric::ErrGeneric(void)
{
	NO_OP;
}

MathParser::ErrGeneric::ErrGeneric(MathParser* p, const char* const s)
{
	silent_cerr(s << " at line " << p->GetLineNumber() << std::endl);
}

MathParser::ErrGeneric::ErrGeneric(MathParser* p,
		const char* const s1,
		const char* const s2,
		const char* const s3)
{
	silent_cerr("MathParser - " << s1 << s2 << s3
			<< " at line " << p->GetLineNumber() << std::endl);
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
	SAFESTRDUP(name, s);
}

MathParser::TokenList::TokenList(const TypedValue& v)
: t(NUM), value(v), name(NULL), next(NULL)
{
	NO_OP;
}

MathParser::TokenList::~TokenList(void)
{
	if (t == NAME) {
		SAFEDELETEARR(name);
	}
}

void
MathParser::TokenPush(enum Token t)
{
	TokenList* p = 0;

	if (t == NUM) {
		SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(value));
	} else if (t == NAME) {
		SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(namebuf));
	} else {
		SAFENEWWITHCONSTRUCTOR(p, TokenList, TokenList(t));
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
	SAFEDELETE(p);

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
		throw ErrGeneric(this, "var '", s, "' already defined!" );
	}

	Var* v = 0;
	switch (t) {
	case TypedValue::VAR_INT:
		SAFENEWWITHCONSTRUCTOR(v, Var, Var(s, Int(d)));
		break;
	case TypedValue::VAR_REAL:
		SAFENEWWITHCONSTRUCTOR(v, Var, Var(s, d));
		break;
	default:
		throw TypedValue::ErrUnknownType();
	}
	VarList* p = 0;
	SAFENEW(p, VarList);
	p->var = v;
	p->next = varlist;
	varlist = p;
	return v;
}

MathParser::NameSpace::NameSpace(const char *const n)
: name(0)
{
	ASSERT(n != 0);
	SAFESTRDUP(name, n);
}

MathParser::NameSpace::~NameSpace(void)
{
	if (name) {
		SAFEDELETEARR(name);
	}
}

const char *
MathParser::NameSpace::sGetName(void) const
{
	return name;
}

MathParser::StaticNameSpace::StaticNameSpace(const char *const n, MathFunc_t *f)
: MathParser::NameSpace(n), func(f)
{
	NO_OP;
}

MathParser::StaticNameSpace::~StaticNameSpace(void)
{
	NO_OP;
}

bool
MathParser::StaticNameSpace::IsFunc(const char* const s) const
{
	for (int i = 0; func[i].fname != NULL; i++) {
		if (strcmp(s, func[i].fname) == 0) {
			return true;
		}
	}
	return false;
}

MathFunc_t*
MathParser::StaticNameSpace::GetFunc(const char* const s) const
{
	for (int i = 0; func[i].fname != NULL; i++) {
		if (strcmp(s, func[i].fname) == 0) {
			return &func[i];
		}
	}
	return NULL;
}

TypedValue
MathParser::StaticNameSpace::EvalFunc(MathFunc_t *f, Real *d) const
{
	switch (f->nargs) {
	case 0:
		return TypedValue((*(f->f.f0))());

	case 1:
		return TypedValue((*(f->f.f1))(d[0]));

	case 2:
		return TypedValue((*(f->f.f2))(d[0], d[1]));

	default:
		throw ErrGeneric();
	}
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

TypedValue::TypeModifier
MathParser::GetTypeModifier(const char* const s) const
{
	for (Int i = 0; TypeModifierNames[i].name != NULL; i++) {
		if (strcmp(s, TypeModifierNames[i].name) == 0) {
			return TypeModifierNames[i].type;
		}
	}
	return TypedValue::MOD_UNKNOWN;
}

bool
MathParser::IsType(const char* const s) const
{
	for (int i = 0; TypeNames[i].name != NULL; i++) {
		if (strcmp(s, TypeNames[i].name) == 0) {
			return true;
		}
	}
	return false;
}

bool
MathParser::IsTypeModifier(const char* const s) const
{
	for (int i = 0; TypeModifierNames[i].name != NULL; i++) {
		if (strcmp(s, TypeModifierNames[i].name) == 0) {
			return true;
		}
	}
	return false;
}

bool
MathParser::IsKeyWord(MathParser::NameSpace *ns, const char* const s) const
{
	if (IsTypeModifier(s)) {
		return true;
	}
	if (IsType(s)) {
		return true;
	}
	if (ns->IsFunc(s)) {
		return true;
	}
	return false;
}

enum MathParser::Token
MathParser::GetToken(void)
{
	ASSERT(in != NULL);

	if (TokenPop()) {
		/* se lo trova! */
		return currtoken;
	}

	int c = 0;

start_parsing:;
	/* skip spaces */
	while ((c = in->get()), isspace(c)) {
		NO_OP;
	};

	if (c == EOF || in->eof()) {
		return (currtoken = ENDOFFILE);
	}

	if (c == ONE_LINE_REMARK) {
		for (c = in->get(); c != '\n'; c = in->get()) {
			if (c == '\\') {
				c = in->get();
				if (c == '\r') {
					c = in->get();
				}
			}
		}
		goto start_parsing;
	}

	/* number? */
	if (c == '.' || isdigit(c)) {
		/* lot of space ... */
		static char s[256];
		bool f = false;
		int i = 0;

		/* FIXME: need to check for overflow */

		s[i++] = c;

		if (c == '.') {
			f = true;
		}
		while ((c = in->get()) == '.' || isdigit(c)) {
			s[i++] = c;
			if (c == '.') {
				if (f) {
					return (currtoken = UNKNOWNTOKEN);
				}
				f = true;
			}
		}
		char e = tolower(c);
		if (e == 'e' || e == 'f' || e == 'd' || e == 'g') {
			f = true;
			s[i++] = 'e';
			if ((c = in->get()) == '-' || c == '+') {
				s[i++] = c;
				c = in->get();
			}
			if (isdigit(c)) {
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
		char *endptr = NULL;
		if (!f) {
			value.SetType(TypedValue::VAR_INT);
#ifdef HAVE_STRTOL
			value.Set(Int(strtol(s, &endptr, 10)));
#else /* !HAVE_STRTOL */
			value.Set(Int(atoi(s)));
#endif /* !HAVE_STRTOL */
		} else {
			value.SetType(TypedValue::VAR_REAL);
#ifdef HAVE_STRTOD
			value.Set(Real(strtod(s, &endptr)));
#else /* !HAVE_STRTOD */
			value.Set(Real(atof(s)));
#endif /* !HAVE_STRTOD */
		}

		if (endptr && endptr[0] != '\0') {
			return (currtoken = UNKNOWNTOKEN);
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
			for (c = in->get();; c = in->get()) {
				if (c == '*') {
end_of_comment:;
					c = in->get();
					if (c == '/') {
						goto start_parsing;
					}
				} else if (c == '/') {
					c = in->get();
					if (c == '*') {
						silent_cerr("warning: '/*' "
							"inside a comment at line "
							<< GetLineNumber()
							<< std::endl);
						goto end_of_comment;
					}
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

	case ':':
		if ((c = in->get()), c != ':') {
			return (currtoken = UNKNOWNTOKEN);
		}
		return (currtoken = NAMESPACESEP);

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
	SAFENEWARR(s, char, namebuflen+1);
	strncpy(s, namebuf, oldlen);
	s[oldlen] = '\0';
	SAFEDELETEARR(namebuf);
	namebuf = s;
}

TypedValue
MathParser::logical(void)
{
	TypedValue d = relational();
	while (true) {
		switch (currtoken) {
		case AND:
			GetToken();
			d = (d && relational());
			break;

		case OR:
			GetToken();
			d = (d || relational());
			break;

		case XOR: {
			GetToken();
			TypedValue e = relational();
			d = ((!(d && e)) && (d || e));
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
	while (true) {
		switch (currtoken) {
		case GT:
			GetToken();
			d = (d > binary());
			break;

		case GE:
			GetToken();
			d = (d >= binary());
			break;

		case EQ:
			GetToken();
			d = (d == binary());
			break;

		case LE:
			GetToken();
			d = (d <= binary());
			break;

		case LT:
			GetToken();
			d = (d < binary());
			break;

		case NE:
			GetToken();
			d = (d != binary());
			break;

		default:
			return d;
		}
	}
}

TypedValue
MathParser::binary(void)
{
	TypedValue d = mult();
	while (true) {
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
	while (true) {
		switch (currtoken) {
		case MULT:
			GetToken();
			d *= power();
			break;
		
		case DIV: {
			GetToken();
			TypedValue e = power();
			if (e == 0.) {
				throw ErrGeneric(this, "divide by zero in mult()");
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
MathParser::power(void)
{
	TypedValue d = unary();

	if (currtoken == EXP) {
		GetToken();

		/*
		 * Per l'esponente chiamo di nuovo power cosi' richiama unary;
		 * se per caso dopo unary c'e' di nuovo un esponente,
		 * l'associazione avviene correttamente da destra:
		 *
		 * 	d^e1^e2 == d^(e1^e2)
		 */
		TypedValue e = power();

		if (d < 0. && e <= 0.) {
			DEBUGCERR("can't compute " << d << '^'
					<< e << " in power()" << std::endl);
			throw ErrGeneric(this, "invalid operands in power()");
		}

		/*
		 * Se sono entrambi interi, uso la sequenza di prodotti
		 * (maggiore accuratezza? comunque va in overflow
		 * correttamente)
		 */
		if (e.GetType() == TypedValue::VAR_INT
				&& d.GetType() == TypedValue::VAR_INT)
		{
			Int i = e.GetInt();
			Int j = d.GetInt();
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
		/*
		 * Altrimenti li forzo entrambi a reale e uso pow
		 * (ottimizzata di suo)
		 */
		} else {
			Real r = e.GetReal();
			Real b = d.GetReal();
			d.SetType(TypedValue::VAR_REAL);
			d = TypedValue(Real(pow(b, r)));
		}
	}

	return d;
}


TypedValue
MathParser::unary(void)
{
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
MathParser::evalfunc(MathParser::NameSpace *ns, MathFunc_t* f)
{
	Real d[max_nargs];
	for (int i = 0; i < f->nargs; i++) {
		d[i] = stmtlist().GetReal();
		if (i != f->nargs - 1) {
			if (currtoken != ARGSEP) {
				throw ErrGeneric(this,
						"argument separator expected");
			}
			GetToken();
		}
	}

	if (f->t != NULL) {
		if ((*(f->t))(d)) {
			DEBUGCERR("error in function " << f->fname << " "
					"(msg: " << f->errmsg << ")"
					<< " in evalfunc()" << std::endl);
			throw ErrGeneric(this, f->fname, ": error ", f->errmsg);
		}
	}

	return ns->EvalFunc(f, d);
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
			throw ErrGeneric(this, "closing parenthesis expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == OPGIN) {
		TypedValue d = readplugin();
		if (currtoken != CPGIN) {
			throw ErrGeneric(this, "closing plugin expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == NAME) {
		GetToken();
		MathParser::NameSpace *currNameSpace = &defaultNameSpace;
		if (currtoken == NAMESPACESEP) {
			std::string name(namebuf);
			NameSpaceMap::iterator i = nameSpaceMap.find(name);
			if (i == nameSpaceMap.end()) {
				throw ErrGeneric(this, "unable to find namespace \"", namebuf, "\"");
			}
			currNameSpace = i->second;
			GetToken();
			if (currtoken != NAME) {
				throw ErrGeneric(this, "name expected after namespace");
			}
			GetToken();
		}

		if (currtoken == OBR) {
			/* in futuro ci potranno essere magari i dati strutturati */
			if (!currNameSpace->IsFunc(namebuf)) {
				throw ErrGeneric(this, "user-defined functions not supported yet!");
			}

			MathFunc_t* f = currNameSpace->GetFunc(namebuf);
			if (f == NULL) {
				throw ErrGeneric(this, "function '", namebuf, "' not found");
			}
			GetToken();
			TypedValue d = evalfunc(currNameSpace, f);
			if (currtoken != CBR) {
				throw ErrGeneric(this, "closing parenthesis "
						"expected after function "
						"\"", f->fname, "\" in expr()");
			}
			GetToken();
			return d;
			
		} else {
			NamedValue* v = table.Get(namebuf);
			if (v != NULL) {
				return v->GetVal();
			}
		}

		throw ErrGeneric(this, "unknown name \"", namebuf, "\"");
	}

	/* invalid expr */
	if (currtoken != ENDOFFILE) {
		throw ErrGeneric(this, "unknown token");
	}

	return TypedValue(0.);
}

TypedValue
MathParser::stmt(void)
{
	if (currtoken == NAME) {
		bool isTypeModifier = false;
		bool isConst = false;

		if (IsTypeModifier(namebuf)) {
			isTypeModifier = true;
			TypedValue::TypeModifier mod = GetTypeModifier(namebuf);
			ASSERT(mod != TypedValue::MOD_UNKNOWN);

			if (mod == TypedValue::MOD_CONST) {
				isConst = true;
			}

			if (GetToken() != NAME || !IsType(namebuf)) {
				throw ErrGeneric(this, "type expected "
						"after type modifier "
						"in declaration");
			}
		}

		/* declaration? */
		if (isTypeModifier || IsType(namebuf)) {
			TypedValue::Type type = GetType(namebuf);
			ASSERT(type != TypedValue::VAR_UNKNOWN);
			
			if (GetToken() != NAME) {
				throw ErrGeneric(this, "name expected "
						"after type in declaration");
			}

			/* FIXME: need to specialize symbol table for namespaces */
			if (IsKeyWord(&defaultNameSpace, namebuf)) {
				throw ErrGeneric(this, "name '", namebuf, "' "
						"is a keyword");
			}

			/* with assign? */
			if (GetToken() == ASSIGN) {
				/* faccio una copia del nome! */
				char* varname = NULL;
				SAFESTRDUP(varname, namebuf);

				GetToken();
				TypedValue d = logical();

				/* cerco la variabile */
				NamedValue* v = table.Get(varname);

				if (v == NULL) {
					/* Se non c'e' la inserisco */
					v = table.Put(varname, TypedValue(type, isConst));
					((Var *)v)->SetVal(d);

				} else {
					/* altrimenti, se la posso ridefinire, mi limito
					 * ad assegnarle il nuovo valore */
					if (redefine_vars) {
						if (v->Const()) {
							silent_cerr("cannot redefine a const named value"
									<< std::endl);
							throw MathParser::ErrGeneric(this,
									"cannot redefine "
									"a const named value "
									"\"", v->GetName(), "\"");
						}

						if (v->IsVar()) {
							((Var *)v)->SetVal(d);
						
						} else {
				   			throw MathParser::ErrGeneric(this,
									"cannot redefine "
									"non-var named value "
									"\"", v->GetName(), "\"");
						}
					} else {
						/* altrimenti la reinserisco,
						 * cosi' da provocare l'errore
						 * di table */
						v = table.Put(varname, TypedValue(type));
					}
				}
				
				/* distruggo il temporaneo */
				SAFEDELETEARR(varname);
				return v->GetVal();

			} else if (currtoken == STMTSEP) {
				NamedValue* v = table.Get(namebuf);
				if (v == NULL || !redefine_vars) {
					/* se la var non esiste, la inserisco;
					 * se invece esiste e non vale 
					 * la ridefinizione, tento
					 * di inserirla comunque, cosi'
					 * table da' errore */
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
					if (v->Const()) {
			      			throw MathParser::ErrGeneric(this,
		      						"cannot assign const named value "
								"\"", v->GetName(), "\"");
			 		}

			 		if (v->IsVar()) {
						((Var *)v)->SetVal(d);
					} else {
						throw MathParser::ErrGeneric(this,
								"cannot assign non-var named value "
								"\"", v->GetName(), "\"");
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
	/*
	 * parse plugin:
	 * 	- arg[0]: nome plugin
	 * 	- arg[1]: nome variabile
	 * 	- arg[2]->arg[n]: dati da passare al costrutture
	 */
	char **argv = NULL, **tmp;
	char c, buf[1024];
	int argc = 0;
	unsigned int i = 0, in_quotes = 0;

	/*
	 * inizializzo l'array degli argomenti
	 */
	SAFENEWARR(argv, char *, 1);
	argv[0] = NULL;

	/*
	 * parserizzo la stringa:
	 * <plugin> ::= '[' <type> ',' <var_name> <list_of_args> ']'
	 * <type> ::= <registered_type>
	 * <var_name> ::= <legal_var_name>
	 * <list_of_args> ::= ',' <arg> <list_of_args> | ''
	 * <arg> ::= <legal_string> (no unescaped ']' or ','!)
	 */
	while ((c = in->get()), !in->eof()) {
		switch (c) {
		case '\\':
			c = in->get();
			if (in_quotes) {
				switch (c) {
				case 'n':
					c = '\n';
					break;

				default:
					in->putback(c);
					c = '\\';
					break;
				}
			}
			buf[i++] = c;
			break;

		case '"':
			if (in_quotes == 0) {
				in_quotes = 1;
				break;
			}
			in_quotes = 0;
			while (isspace((c = in->get()))) {
				NO_OP;
			}
			if (c != ',' && c != ']') {
				silent_cerr("need a separator "
					"after closing quotes" << std::endl);
				throw ErrGeneric();
			}
			in->putback(c);
			break;

		case ',':
		case ']':
			if (in_quotes) {
				buf[i++] = c;
				break;
			}
			buf[i] = '\0';
			tmp = argv;
			argv = NULL;
			SAFENEWARR(argv, char *, argc+2);
			memcpy(argv, tmp, sizeof(char *)*argc);
			SAFEDELETEARR(tmp);
			argv[argc] = NULL;
			trim_arg(buf);
			SAFESTRDUP(argv[argc], buf);
			++argc;
			argv[argc] = NULL;
			if (c == ']') {
				goto last_arg;
			}
			i = 0;
			break;

		default:
			buf[i++] = c;
			break;
		}

		/*
		 * FIXME: rendere dinamico il buffer ...
		 */
		if (i >= sizeof(buf)) {
			silent_cerr("buffer overflow" << std::endl);
			throw ErrGeneric();
		}
	}

last_arg:
	if (in->eof()) {
		silent_cerr("eof encountered while parsing plugin"
			<< std::endl);
		throw ErrGeneric();
	}

	/*
	 * put the close plugin token back
	 */
	in->putback(c);
	buf[i] = '\0';

	/*
	 * argomenti comuni a tutti i plugin
	 */
	char *pginname = argv[0];
	char *varname = argv[1];
	trim_arg(pginname);
	trim_arg(varname);

	/*
	 * verifiche di validita' argomenti
	 */
	if (pginname == NULL || *pginname == '\0') {
		silent_cerr("illegal or missing plugin name" << std::endl);
		throw ErrGeneric();
	}

	if (varname == NULL || *varname == '\0') {
		silent_cerr("illegal or missing plugin variable name"
			<< std::endl);
		throw ErrGeneric();
	}

	/*
	 * verifica esistenza nome
	 */
	NamedValue* v = table.Get(varname);
	if (v != NULL) {
		silent_cerr("variable " << varname << " already defined"
			<< std::endl);
		throw ErrGeneric();
	}

	/*
	 * ricerca registrazione plugin
	 */
	for (struct PlugInRegister *p = PlugIns; p != NULL; p = p->next) {
		if (strcasecmp(p->name, pginname) != 0) {
			continue;
		}
#ifdef DEBUG
		for (int i = 0; argv[i] != NULL; i++) {
			silent_cout("argv[" << i << "]=" << argv[i]
					<< std::endl);
		}
#endif /* DEBUG */

		/*
		 * costruisce il plugin e gli fa interpretare gli argomenti
		 */
		MathParser::PlugIn *pgin = (*p->constructor)(*this, p->arg);
		pgin->Read(argc-2, argv+2);

		/*
		 * riporta il parser nello stato corretto
		 */
		GetToken();

		/*
		 * costruisce la variabile, la inserisce nella tabella
		 * e ne ritorna il valore (prima esecuzione)
		 */
		SAFENEWWITHCONSTRUCTOR(v, PlugInVar,
				PlugInVar(varname, pgin));
		table.Put(v);

		/*
		 * pulizia ...
		 */
		for (int i = 0; argv[i] != NULL; i++) {
			SAFEDELETEARR(argv[i]);
		}
		SAFEDELETEARR(argv);

		return v->GetVal();
	}

	/*
	 * si arriva qui solo se il plugin non e' stato registrato
	 */
	silent_cerr("plugin '" << pginname << "' not supported" << std::endl);
	throw ErrGeneric();
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
defaultNameSpace("default", FuncTable),
namebuf(NULL),
namebuflen(default_namebuflen),
value(Real(0)),
varlist(NULL),
tokenlist(NULL),
powerstack()
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	if (RegisterNameSpace(&defaultNameSpace)) {
		throw ErrGeneric(this, "unable to register namespace "
				"\"", defaultNameSpace.sGetName(), "\"");
	}

	time_t tm;
	time(&tm);
	srand(tm);
}

MathParser::MathParser(Table& t, int redefine_vars)
: PlugIns(NULL),
table(t),
redefine_vars(redefine_vars),
in(NULL),
defaultNameSpace("default", FuncTable),
namebuf(NULL),
namebuflen(default_namebuflen),
value(Real(0)),
varlist(NULL),
tokenlist(NULL),
powerstack()
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	if (RegisterNameSpace(&defaultNameSpace)) {
		throw ErrGeneric(this, "unable to register namespace "
				"\"", defaultNameSpace.sGetName(), "\"");
	}

	time_t tm;
	time(&tm);
	srand(tm);
}

NamedValue *
MathParser::InsertSym(const char* const s, const Real& v, int redefine)
{
	/* cerco la variabile */
	NamedValue* var = table.Get(s);

	if (var == NULL) {
		/* Se non c'e' la inserisco */
		var = table.Put(s, TypedValue(v));
	} else {
		/* altrimenti, se la posso ridefinire, mi limito
		 * ad assegnarle il nuovo valore */
		if (redefine) {
			if (var->IsVar()) {
				((Var *)var)->SetVal(TypedValue(v));
			} else {
				throw MathParser::ErrGeneric(this,
						"cannot redefine "
						"non-var " "named value "
						"\"", var->GetName(), "\"");
			}

		} else {
			/* altrimenti la reinserisco, cosi' 
			 * da provocare l'errore di table */
			var = table.Put(s, TypedValue(v));
		}
	}

	if (var == NULL) {
		throw ErrGeneric(this, "error while adding real var '", s, "");
	}

	return var;
}

NamedValue *
MathParser::InsertSym(const char* const s, const Int& v, int redefine)
{
	/* cerco la variabile */
	NamedValue* var = table.Get(s);

	if (var == NULL) {
		/* Se non c'e' la inserisco */
		var = table.Put(s, TypedValue(v));

	} else {
		/* altrimenti, se la posso ridefinire, mi limito
		 * ad assegnarle il nuovo valore */
		if (redefine) {
			if (var->IsVar()) {
				((Var *)var)->SetVal(TypedValue(v));

			} else {
				throw MathParser::ErrGeneric(this,
						"cannot redefine "
						"non-var named value "
						"\"", var->GetName(), "\"");
			}

		} else {
			/* altrimenti la reinserisco, cosi' 
			 * da provocare l'errore di table */
			var = table.Put(s, TypedValue(v));
		}
	}

	if (var == NULL) {
		throw ErrGeneric(this, "error while adding integer var '", s, "");
	}

	return var;
}

MathParser::~MathParser(void)
{
	DEBUGCOUTFNAME("MathParser::~MathParser");

	if (namebuf != NULL) {
		SAFEDELETEARR(namebuf);
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
	d = GetLastStmt(d, t);
	in = (InputStream*)p;
	return d;
}

Real
MathParser::Get(Real d)
{
	GetToken();
	d = stmt().GetReal();
	if (currtoken != STMTSEP && currtoken != ENDOFFILE) {
		throw ErrGeneric(this, "statement separator expected");
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
		throw ErrGeneric(this, "separator expected");
	}
	in = (InputStream*)p;
	return d;
}

void
MathParser::GetForever(std::ostream& out, const char* const sep)
{
	do {
		out << Get(0.) << sep;
	} while (currtoken == STMTSEP);
}

void
MathParser::GetForever(const InputStream& strm, std::ostream& out,
		const char* const sep)
{
 	const InputStream* p = in;
	in = (InputStream*)&strm;
	GetForever(out);
	in = (InputStream*)p;
}

int
MathParser::RegisterPlugIn(const char *name,
		MathParser::PlugIn * (*constructor)(MathParser&, void *),
		void *arg)
{
	PlugInRegister *p = NULL;
	SAFENEW(p, PlugInRegister);
	SAFESTRDUP(p->name, name);
	p->constructor = constructor;
	p->arg = arg;
	p->next = PlugIns;
	PlugIns = p;

	return 0;
}

int
MathParser::RegisterNameSpace(MathParser::NameSpace *ns)
{
	ASSERT(ns != 0);

	silent_cout("registering namespace \"" << ns->sGetName() << "\""
			<< std::endl);

	std::string name(ns->sGetName());

	if (nameSpaceMap.find(name) != nameSpaceMap.end()) {
		return 1;
	}

	nameSpaceMap[name] = ns;

	return 0;
}
