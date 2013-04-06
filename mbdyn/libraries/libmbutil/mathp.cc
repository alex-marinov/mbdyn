/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <cfloat>
#include <cstdlib>
#include <climits>
#include <limits>

#include "mathp.h"
#include "parser.h"


/* helper per le funzioni built-in */
typedef double (*mp_f1_f)(double);
typedef double (*mp_f2_f)(double, double);

template <class Tin, class Tout, mp_f1_f F>
static int
mp_func_1(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	Tout *out = dynamic_cast<Tout *>(args[0]);
	ASSERT(out != 0);

	Tin *arg1 = dynamic_cast<Tin *>(args[1]);
	ASSERT(arg1 != 0);

	*out = F((*arg1)());

	return 0;
}

template <class Tin, class Tout, mp_f2_f F>
static int
mp_func_2(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);

	Tout *out = dynamic_cast<Tout *>(args[0]);
	ASSERT(out != 0);

	Tin *arg1 = dynamic_cast<Tin *>(args[1]);
	ASSERT(arg1 != 0);

	Tin *arg2 = dynamic_cast<Tin *>(args[2]);
	ASSERT(arg2 != 0);

	Real a1 = (*arg1)();
	Real a2 = (*arg2)();
	*out = F(a1, a2);

	return 0;
}

static int
mp_asin_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	const Real a1 = (*arg1)();
	if (a1 > 1. || a1 < -1.) {
		return 1;
	}

	return 0;
}

static int
mp_acos_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	const Real a1 = (*arg1)();
	if (a1 > 1. || a1 < -1.) {
		return 1;
	}

	return 0;
}

static int
mp_tan_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	Real a = a1;
	a -= int(a1/M_PI)*M_PI;
	if (fabs(fabs(a) - M_PI_2) < std::numeric_limits<double>::epsilon()) {
		return 1;
	}

	return 0;
}

static int
mp_acosh_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	if ((*arg1)() <= 1.) {
		return 1;
	}

	return 0;
}

static int
mp_atanh_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	if (a1 >= 1. || a1 <= -1.) {
		return 1;
	}

	return 0;
}

static int
mp_greater_than_0_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	if ((*arg1)() <= 0.) {
		return 1;
	}

	return 0;
}

static int
mp_greater_than_or_equal_to_0_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	if ((*arg1)() < 0.) {
		return 1;
	}

	return 0;
}

static int
mp_rand(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 0);
	ASSERT(args[0]->Type() == MathParser::AT_INT);

	MathParser::MathArgInt_t* out = dynamic_cast<MathParser::MathArgInt_t*>(args[0]);
	ASSERT(out != 0);

	*out = rand();

	return 0;
}

static int
mp_rndm(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 0);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t* out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	*out = -1. + 2.*(Real(rand())/Real(RAND_MAX));

	return 0;
}

static int
mp_srnd(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_VOID);
	ASSERT(args[1]->Type() == MathParser::AT_INT);

	MathParser::MathArgInt_t *arg1 = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(arg1 != 0);

	srand((unsigned int)(*arg1)());

	return 0;
}

std::ostream&
operator << (std::ostream& out, const MathParser::MathArgVoid_t& /* v */ )
{
	return out;
}

std::ostream&
operator << (std::ostream& out, const MathParser::MathArgBool_t& v)
{
	return out << v();
}

std::ostream&
operator << (std::ostream& out, const MathParser::MathArgInt_t& v)
{
	return out << v();
}

std::ostream&
operator << (std::ostream& out, const MathParser::MathArgReal_t& v)
{
	return out << v();
}

std::ostream&
operator << (std::ostream& out, const MathParser::MathArgString_t& v)
{
	return out << v();
}

static int
mp_print(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);

	switch (args[1]->Type()) {
	case MathParser::AT_BOOL:
		silent_cout((*dynamic_cast<MathParser::MathArgBool_t*>(args[1])) << std::endl);
		break;

	case MathParser::AT_INT:
		silent_cout((*dynamic_cast<MathParser::MathArgInt_t*>(args[1])) << std::endl);
		break;

	case MathParser::AT_REAL:
		silent_cout((*dynamic_cast<MathParser::MathArgReal_t*>(args[1])) << std::endl);
		break;

	case MathParser::AT_STRING:
		silent_cout((*dynamic_cast<MathParser::MathArgString_t*>(args[1])) << std::endl);
		break;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return 0;
}

static int
mp_stop(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_VOID);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[2]->Type() == MathParser::AT_INT);

	MathParser::MathArgInt_t *s = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(s != 0);

	MathParser::MathArgInt_t *v = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(v != 0);

	if ((*s)() != 0) {
		if ((*v)() == 0) {
			silent_cout("mp_stop(SUCCESS)" << std::endl);
			throw NoErr(MBDYN_EXCEPT_ARGS);

		} else {
			silent_cout("mp_stop(FAILURE)" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	return 0;
}

static int
mp_ctg(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t *>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	*out = 1./tan((*arg1)());

	return 0;
}

static int
mp_ctg_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	Real a = a1;
	a -= int(a1/M_PI)*M_PI;
	if (fabs(a) < std::numeric_limits<double>::epsilon()) {
		return 1;
	}

	return 0;
}

static int
mp_actg(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	*out = atan2(1., (*arg1)());

	return 0;
}

static int
mp_actg2(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t*>(args[2]);
	ASSERT(arg2 != 0);

	Real a1 = (*arg1)();
	Real a2 = (*arg2)();

	*out = atan2(a2, a1);

	return 0;
}

static int
mp_ctgh(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	*out = 1./tanh((*arg1)());

	return 0;
}

static int
mp_ctgh_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	if (fabs((*arg1)()) < std::numeric_limits<double>::epsilon()) {
		return 1;
	}

	return 0;
}

static int
mp_sign(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgInt_t *out = dynamic_cast<MathParser::MathArgInt_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	if (a1 == 0.) {
		*out = 0;
	} else {
		*out = Int(copysign(1., a1));
	}

	return 0;
}

static int
mp_max(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t*>(args[2]);
	ASSERT(arg2 != 0);

	Real a1 = (*arg1)();
	Real a2 = (*arg2)();
	*out = std::max(a1, a2);

	return 0;
}

static int
mp_min(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t*>(args[2]);
	ASSERT(arg2 != 0);

	Real a1 = (*arg1)();
	Real a2 = (*arg2)();
	*out = std::min(a1, a2);

	return 0;
}

#ifdef __USE_XOPEN
static int
mp_actgh(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	*out = atanh(1./(*arg1)());

	return 0;
}
#endif

static int
mp_step(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	if (a1 > 0.) {
		*out = 1.;

	} else if (a1 < 0.) {
		*out = 0.;

	} else {
		*out = .5;
	}

	return 0;
}

static int
mp_ramp(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	if (a1 > 0.) {
		*out = a1;

	} else {
		*out = 0.;
	}

	return 0;
}

static int
mp_sramp(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t*>(args[2]);
	ASSERT(arg2 != 0);

	Real a1 = (*arg1)();
	Real a2 = (*arg2)();
	if (a1 < 0.) {
		*out = 0.;

	} else if (a1 > a2) {
		*out = a2;

	} else {
		*out = a1;
	}

	return 0;
}

static int
mp_par(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	Real a1 = (*arg1)();
	if (a1 > 0.) {
		*out = a1*a1;

	} else {
		*out = 0.;
	}

	return 0;
}

enum mp_in_e {
	IN_LL,
	IN_LE,
	IN_EL,
	IN_EE
};

template<mp_in_e IN>
static int
mp_in(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 3);
	ASSERT(args[0]->Type() == MathParser::AT_REAL);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);
	ASSERT(args[2]->Type() == MathParser::AT_REAL);
	ASSERT(args[3]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	MathParser::MathArgReal_t *arg2 = dynamic_cast<MathParser::MathArgReal_t*>(args[2]);
	ASSERT(arg2 != 0);

	MathParser::MathArgReal_t *arg3 = dynamic_cast<MathParser::MathArgReal_t*>(args[3]);
	ASSERT(arg3 != 0);

	Real l = (*arg1)();
	Real x = (*arg2)();
	Real u = (*arg3)();

	switch (IN) {
	case IN_LL:
	case IN_LE:
		if (x <= l) {
			*out = 0;
			return 0;
		}
		break;

	case IN_EL:
	case IN_EE:
		if (x < l) {
			*out = 0;
			return 0;
		}
		break;
	}

	switch (IN) {
	case IN_LL:
	case IN_EL:
		if (x >= u) {
			*out = 0;
			return 0;
		}
		break;

	case IN_LE:
	case IN_EE:
		if (x > u) {
			*out = 0;
			return 0;
		}
		break;
	}

	*out = 1;

	return 0;
}

static int
mp_cast(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[0]->Type() == MathParser::AT_ANY);

	MathParser::MathArgAny_t *arg1 = dynamic_cast<MathParser::MathArgAny_t*>(args[1]);
	ASSERT(arg1 != 0);

	if (args[0]->Type() == MathParser::AT_BOOL) {
		MathParser::MathArgBool_t *out = dynamic_cast<MathParser::MathArgBool_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetBool();

	} else if (args[0]->Type() == MathParser::AT_INT) {
		MathParser::MathArgInt_t *out = dynamic_cast<MathParser::MathArgInt_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetInt();

	} else if (args[0]->Type() == MathParser::AT_REAL) {
		MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetReal();

	} else if (args[0]->Type() == MathParser::AT_STRING) {
		MathParser::MathArgString_t *out = dynamic_cast<MathParser::MathArgString_t*>(args[0]);
		ASSERT(out != 0);
		TypedValue val(TypedValue::VAR_STRING);
		val.Cast((*arg1)());
		(*out)() = val.GetString();

	} else {
	}

	return 0;
}

/* tipi delle variabili */
struct TypeName_t {
	const char* name;
	TypedValue::Type type;
};

static TypeName_t TypeNames[] = {
	{ "bool",	TypedValue::VAR_BOOL },
	{ "integer",	TypedValue::VAR_INT },
	{ "real",	TypedValue::VAR_REAL },
	{ "string",	TypedValue::VAR_STRING },

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

struct declarationmodifiernames {
	const char* name;
	MathParser::DeclarationModifier type;
};

static declarationmodifiernames DeclarationModifierNames[] = {
	{ "ifndef",	MathParser::DMOD_IFNDEF },
	{ NULL,		MathParser::DMOD_UNKNOWN }
};

TypedValue::ErrWrongType::ErrWrongType(const char *file, int line, const char *func,
	const TypedValue::Type& to,
	const TypedValue::Type& from)
: MBDynErrBase(file, line, func, std::string("cannot cast \"") + GetTypeName(from) + "\" into \"" + GetTypeName(to) + "\"")
{
	NO_OP;
}

TypedValue::TypedValue(void)
: type(TypedValue::VAR_UNKNOWN), bConst(false)
{
	NO_OP;
}

TypedValue::~TypedValue(void)
{
	NO_OP;
}

TypedValue::TypedValue(const bool& b, bool isConst)
: type(TypedValue::VAR_BOOL), bConst(isConst)
{
	v.i = b;
}

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

TypedValue::TypedValue(const std::string& s, bool isConst)
: type(TypedValue::VAR_STRING), bConst(isConst), s(s)
{
	NO_OP;
}

TypedValue::TypedValue(const TypedValue::Type t, bool isConst)
: type(t), bConst(isConst)
{
	switch (type) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
		v.i = 0;
		break;

	case TypedValue::VAR_REAL:
		v.r = 0.;
		break;

	case TypedValue::VAR_STRING:
		s = "";
		break;

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

// NOTE: TypedValue does not inherit const'ness
TypedValue::TypedValue(const TypedValue& var)
: type(var.type), bConst(false)
{
	switch (type) {
	case TypedValue::VAR_BOOL:
		ASSERT(var.v.i == 0 || var.v.i == 1);
		// fallthru
	case TypedValue::VAR_INT:
		v.i = var.v.i;
		break;

	case TypedValue::VAR_REAL:
		v.r = var.v.r;
		break;

	case TypedValue::VAR_STRING:
		s = var.s;
		break;

	case TypedValue::VAR_UNKNOWN:
		/* allow initialization from unspecified type right now */
		break;

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

const TypedValue&
TypedValue::operator = (const TypedValue& var)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	bConst = false;
	if (type == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (var.type) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)var.GetInt());
			Set(std::string(buf));
			break;

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)var.GetReal());
			Set(std::string(buf));
			break;

		case TypedValue::VAR_STRING:
			this->s = var.GetString();
			break;

		default:
			throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
		}

	} else {
		switch (var.type) {
		case TypedValue::VAR_BOOL:
			type = var.type;
			Set(var.GetBool());
			break;

		case TypedValue::VAR_INT:
			type = var.type;
			Set(var.GetInt());
			break;

		case TypedValue::VAR_REAL:
			type = var.type;
			Set(var.GetReal());
			break;

		case TypedValue::VAR_STRING:
			if (type != TypedValue::VAR_UNKNOWN) {
				throw ErrWrongType(MBDYN_EXCEPT_ARGS, type, var.type);
			}
			type = TypedValue::VAR_STRING;
			s = var.s;
			break;

		case TypedValue::VAR_UNKNOWN:
			if (type != TypedValue::VAR_UNKNOWN) {
				throw ErrWrongType(MBDYN_EXCEPT_ARGS, type, var.type);
			}
			break;

		default:
			throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
		}
	}

	// NOTE: TypedValue does not inherit const'ness
	return *this;
}

const TypedValue&
TypedValue::Cast(const TypedValue& var)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if (type == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (var.type) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)var.GetInt());
			Set(std::string(buf));
			break;

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)var.GetReal());
			Set(std::string(buf));
			break;

		case TypedValue::VAR_STRING:
			this->s = var.GetString();
			break;

		default:
			throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
		}

	} else {
		switch (type) {
		case TypedValue::VAR_BOOL:
			switch (var.type) {
			case TypedValue::VAR_BOOL:
				break;

			case TypedValue::VAR_INT:
			case TypedValue::VAR_REAL:
				silent_cerr("Warning: implicit cast "
					"from " << GetTypeName(var.type)
					<< " to " << GetTypeName(type)
					<< " may alter its value"
					<< " at line " << mbdyn_get_line_data()
					<< std::endl);
				break;

			default:
				throw ErrWrongType(MBDYN_EXCEPT_ARGS, type, var.type);
			}
			Set(var.GetBool());
			break;

		case TypedValue::VAR_INT:
			switch (var.type) {
			case TypedValue::VAR_BOOL:
			case TypedValue::VAR_INT:
				break;

			case TypedValue::VAR_REAL:
				silent_cerr("Warning: implicit cast "
					"from " << GetTypeName(var.type)
					<< " to " << GetTypeName(type)
					<< " may alter its value"
					<< " at line " << mbdyn_get_line_data()
					<< std::endl);
				break;

			default:
				throw ErrWrongType(MBDYN_EXCEPT_ARGS, type, var.type);
			}
			Set(var.GetInt());
			break;

		case TypedValue::VAR_REAL:
			switch (var.type) {
			case TypedValue::VAR_BOOL:
			case TypedValue::VAR_INT:
			case TypedValue::VAR_REAL:
				break;

			default:
				throw ErrWrongType(MBDYN_EXCEPT_ARGS, type, var.type);
			}
			Set(var.GetReal());
			break;

		default:
			throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
		}
	}

	bConst = var.Const();
	return *this;
}

TypedValue::Type
TypedValue::GetType(void) const
{
	return type;
}

const char *const
TypedValue::GetTypeName(TypedValue::Type t)
{
	switch (t) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
	case TypedValue::VAR_REAL:
	case TypedValue::VAR_STRING:
		return TypeNames[t].name;

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

const char *const
TypedValue::GetTypeName(void) const
{
	return GetTypeName(type);
}

bool
TypedValue::Const(void) const
{
	return bConst;
}

bool
TypedValue::GetBool(void) const
{
	switch (type) {
	case TypedValue::VAR_BOOL:
		return v.i;

	case TypedValue::VAR_INT:
		return v.i ? 1 : 0;

	case TypedValue::VAR_REAL:
		return v.r ? 1 : 0;

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, VAR_BOOL, type);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

Int
TypedValue::GetInt(void) const
{
	switch (type) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
		return v.i;

	case TypedValue::VAR_REAL:
		return Int(v.r);

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, VAR_INT, type);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

Real
TypedValue::GetReal(void) const
{
	switch (type) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
		return Real(v.i);

	case TypedValue::VAR_REAL:
 		return v.r;

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, VAR_REAL, type);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

const std::string &
TypedValue::GetString(void) const
{
	switch (type) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
	case TypedValue::VAR_REAL:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, VAR_STRING, type);

	case TypedValue::VAR_STRING:
 		return s;

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}
}

void
TypedValue::SetType(TypedValue::Type t, bool isConst)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	type = t;
	bConst = isConst;
}

void
TypedValue::SetConst(bool isConst, bool bForce)
{
	if (Const() && !isConst && !bForce) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	bConst = isConst;
}

const TypedValue&
TypedValue::Set(const bool& b)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	switch (GetType()) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
		v.i = b ? 1 : 0;
		break;

	case TypedValue::VAR_REAL:
		v.r = b ? 1. : 0.;
		break;

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), VAR_BOOL);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}

	return *this;
}

const TypedValue&
TypedValue::Set(const Int& i)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	switch (GetType()) {
	case TypedValue::VAR_BOOL:
		v.i = i ? 1 : 0;
		break;

	case TypedValue::VAR_INT:
		v.i = i;
		break;

	case TypedValue::VAR_REAL:
		v.r = Real(i);
		break;

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), VAR_INT);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}

	return *this;
}

const TypedValue&
TypedValue::Set(const Real& r)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	switch (GetType()) {
	case TypedValue::VAR_BOOL:
		v.i = r ? 1 : 0;
		break;

	case TypedValue::VAR_INT:
		v.i = Int(r);
		break;

	case TypedValue::VAR_REAL:
		v.r = r;
		break;

	case TypedValue::VAR_STRING:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), VAR_REAL);

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
	}

	return *this;
}

const TypedValue&
TypedValue::Set(const std::string& s)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	switch (GetType()) {
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
	case TypedValue::VAR_REAL:
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), VAR_STRING);

	case TypedValue::VAR_STRING:
		this->s = s;
		break;

	default:
		throw ErrUnknownType(MBDYN_EXCEPT_ARGS);
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
	if (GetType() == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (v.GetType()) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)v.GetInt());
			return TypedValue(GetString() + buf);

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)v.GetReal());
			return TypedValue(GetString() + buf);

		case TypedValue::VAR_STRING:
			return TypedValue(GetString() + v.GetString());

		default:
			throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), v.GetType());
		}
	}

	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return TypedValue(GetInt() + v.GetInt());
	}

	return TypedValue(GetReal() + v.GetReal());
}

TypedValue
TypedValue::operator - (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return TypedValue(GetInt() - v.GetInt());
	}

	return TypedValue(GetReal() - v.GetReal());
}

TypedValue
TypedValue::operator * (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return TypedValue(GetInt()*v.GetInt());
	}

	return TypedValue(GetReal()*v.GetReal());
}

TypedValue
TypedValue::operator / (const TypedValue& v) const
{
	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return TypedValue(GetInt()/v.GetInt());
	}

	return TypedValue(GetReal()/v.GetReal());
}

TypedValue
TypedValue::operator % (const TypedValue& v) const
{
	// bool is implicitly cast to Int
	if ((GetType() != TypedValue::VAR_BOOL && GetType() != TypedValue::VAR_INT)
		|| (v.GetType() != TypedValue::VAR_BOOL && v.GetType() != TypedValue::VAR_INT))
	{
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), v.GetType());
	}

	return TypedValue(GetInt()%v.GetInt());
}

const TypedValue&
TypedValue::operator += (const TypedValue& v)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if (GetType() == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (v.GetType()) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)v.GetInt());
			this->s += buf;
			break;

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)v.GetReal());
			this->s += buf;
			break;

		case TypedValue::VAR_STRING:
			this->s += v.GetString();
			break;

		default:
			throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), v.GetType());
		}

		return *this;
	}

	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return Set(GetInt() + v.GetInt());
	}

	Real d = GetReal() + v.GetReal();
	type = TypedValue::VAR_REAL;

	return Set(d);
}

const TypedValue&
TypedValue::operator -= (const TypedValue& v)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return Set(GetInt() - v.GetInt());
	}
	Real d = GetReal() - v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

const TypedValue&
TypedValue::operator *= (const TypedValue& v)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return Set(GetInt()*v.GetInt());
	}
	Real d = GetReal()*v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

const TypedValue&
TypedValue::operator /= (const TypedValue& v)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if ((GetType() == TypedValue::VAR_BOOL || GetType() == TypedValue::VAR_INT)
			&& (v.GetType() == TypedValue::VAR_BOOL || v.GetType() == TypedValue::VAR_INT))
	{
		// bool is implicitly cast to Int
		return Set(GetInt()/v.GetInt());
	}
	Real d = GetReal()/v.GetReal();
	type = TypedValue::VAR_REAL;
	return Set(d);
}

const TypedValue&
TypedValue::operator %= (const TypedValue& v)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	// bool is implicitly cast to Int
	if ((GetType() != TypedValue::VAR_BOOL && GetType() != TypedValue::VAR_INT)
			|| (v.GetType() != TypedValue::VAR_BOOL && v.GetType() != TypedValue::VAR_INT))
	{
		throw ErrWrongType(MBDYN_EXCEPT_ARGS, GetType(), v.GetType());
	}

	Int i = GetInt()%v.GetInt();
	type = TypedValue::VAR_INT;
	return Set(i);
}

bool
operator ! (const TypedValue& v)
{
	if (v.GetType() == TypedValue::VAR_STRING) {
		return v.GetString().empty();
	}

	return (!v.GetReal());
}

TypedValue
operator - (const TypedValue& v)
{
	switch (v.GetType()) {
	case TypedValue::VAR_BOOL:
		// bool is implicitly cast to Int
	case TypedValue::VAR_INT:
		return TypedValue(-v.GetInt());

	case TypedValue::VAR_REAL:
		return TypedValue(-v.GetReal());

	case TypedValue::VAR_STRING:
		throw TypedValue::ErrWrongType(MBDYN_EXCEPT_ARGS);

	default:
		throw TypedValue::ErrUnknownType(MBDYN_EXCEPT_ARGS);
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
	case TypedValue::VAR_BOOL:
	case TypedValue::VAR_INT:
		return out << v.GetInt();

	case TypedValue::VAR_REAL:
		return out << v.GetReal();

	case TypedValue::VAR_STRING:
		return out << v.GetString();

	default:
		throw TypedValue::ErrUnknownType(MBDYN_EXCEPT_ARGS);
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

const char *const
NamedValue::GetTypeName(void) const
{
	return TypedValue::GetTypeName(GetType());
}

Var::Var(const char* const s, const TypedValue& v)
: NamedValue(s), value(v)
{
	value.SetConst(v.Const(), true);
}

Var::Var(const char* const s, const bool& v)
: NamedValue(s), value(v)
{
	NO_OP;
}

Var::Var(const char* const s, const Int& v)
: NamedValue(s), value(v)
{
	NO_OP;
}

Var::Var(const char* const s, const Real& v)
: NamedValue(s), value(v)
{
	NO_OP;
}

Var::Var(const char* const s, const std::string& v)
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
Var::SetVal(const bool& v)
{
	value.Set(v);
}

void
Var::SetVal(const Int& v)
{
	value.Set(v);
}

void
Var::SetVal(const Real& v)
{
	value.Set(v);
}

void
Var::SetVal(const std::string& v)
{
	value.Set(v);
}

void
Var::SetVal(const TypedValue& v)
{
	value = v;
}

void
Var::Cast(const TypedValue& v)
{
	value.Cast(v);
}

void
MathParser::trim_arg(char *const s)
{
	int i, l;

	for (i = 0; isspace(s[i]); ++i) {
		NO_OP;
	}

	l = strlen(&s[i]);
	if (i > 0) {
		memmove(s, &s[i], l + 1);
	}

	for (i = l - 1; isspace(s[i]); --i) {
		NO_OP;
	}
	s[i + 1] = '\0';
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

MathParser::ErrGeneric::ErrGeneric(MBDYN_EXCEPT_ARGS_DECL_NODEF) : 
	MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};

MathParser::ErrGeneric::ErrGeneric(MathParser* p, MBDYN_EXCEPT_ARGS_DECL_NODEF) :
	MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {
	silent_cerr(what() << " at line " << p->GetLineNumber() << std::endl);
}

MathParser::ErrGeneric::ErrGeneric(MathParser* p,
	MBDYN_EXCEPT_ARGS_DECL_NODEF, 
	const char* const s2,
	const char* const s3)
: MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU + s2 + s3)
{
	silent_cerr("MathParser - " << r << s2 << s3
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
	table = T;
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
MathParser::TokenPush(Token t)
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

MathParser::NameSpace::NameSpace(const std::string& name)
: name(name)
{
	NO_OP;
}

MathParser::NameSpace::~NameSpace(void)
{
	NO_OP;
}

const std::string&
MathParser::NameSpace::sGetName(void) const
{
	return name;
}

MathParser::StaticNameSpace::StaticNameSpace()
: MathParser::NameSpace("default")
{
	MathFunc_t	*f;

	// asin
	f = new MathFunc_t;
	f->fname = std::string("asin");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, asin>;
	f->t = mp_asin_t;
	f->errmsg = std::string("invalid arg to asin()");

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// acos
	f = new MathFunc_t;
	f->fname = std::string("acos");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, acos>;
	f->t = mp_acos_t;
	f->errmsg = std::string("invalid arg to acos()");

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// atan
	f = new MathFunc_t;
	f->fname = std::string("atan");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, atan>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// actan
	f = new MathFunc_t;
	f->fname = std::string("actan");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_actg;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// atan2
	f = new MathFunc_t;
	f->fname = std::string("atan2");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_func_2<MathParser::MathArgReal_t, MathParser::MathArgReal_t, atan2>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// actan2
	f = new MathFunc_t;
	f->fname = std::string("actan2");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_actg2;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// cos
	f = new MathFunc_t;
	f->fname = std::string("cos");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, cos>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// sin
	f = new MathFunc_t;
	f->fname = std::string("sin");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, sin>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// tan
	f = new MathFunc_t;
	f->fname = std::string("tan");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, tan>;
	f->t = mp_tan_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ctan
	f = new MathFunc_t;
	f->fname = std::string("ctan");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_ctg;
	f->t = mp_ctg_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// cosh
	f = new MathFunc_t;
	f->fname = std::string("cosh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, cosh>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// sinh
	f = new MathFunc_t;
	f->fname = std::string("sinh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, sinh>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// tanh
	f = new MathFunc_t;
	f->fname = std::string("tanh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, tanh>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ctanh
	f = new MathFunc_t;
	f->fname = std::string("ctanh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_ctgh;
	f->t = mp_ctgh_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef __USE_XOPEN
	// acosh
	f = new MathFunc_t;
	f->fname = std::string("acosh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, acosh>;
	f->t = mp_acosh_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// asinh
	f = new MathFunc_t;
	f->fname = std::string("asinh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, asinh>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// atanh
	f = new MathFunc_t;
	f->fname = std::string("atanh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, atanh>;
	f->t = mp_atanh_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// actanh
	f = new MathFunc_t;
	f->fname = std::string("actanh");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_actgh;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif /* __USE_XOPEN */

	// exp
	f = new MathFunc_t;
	f->fname = std::string("exp");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, exp>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// log
	f = new MathFunc_t;
	f->fname = std::string("log");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, log>;
	f->t = mp_greater_than_0_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// log10
	f = new MathFunc_t;
	f->fname = std::string("log10");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, log10>;
	f->t = mp_greater_than_0_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// sqrt
	f = new MathFunc_t;
	f->fname = std::string("sqrt");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, sqrt>;
	f->t = mp_greater_than_or_equal_to_0_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// abs
	f = new MathFunc_t;
	f->fname = std::string("abs");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, fabs>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// sign
	f = new MathFunc_t;
	f->fname = std::string("sign");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgInt_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_sign;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// copysign
	f = new MathFunc_t;
	f->fname = std::string("copysign");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_func_2<MathParser::MathArgReal_t, MathParser::MathArgReal_t, copysign>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// max
	f = new MathFunc_t;
	f->fname = std::string("max");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_max;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// min
	f = new MathFunc_t;
	f->fname = std::string("min");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_min;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// floor
	f = new MathFunc_t;
	f->fname = std::string("floor");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, floor>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ceil
	f = new MathFunc_t;
	f->fname = std::string("ceil");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, ceil>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

#ifdef __USE_XOPEN
	// round
	f = new MathFunc_t;
	f->fname = std::string("round");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, rint>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
#endif /* __USE_XOPEN */

	// rand
	f = new MathFunc_t;
	f->fname = std::string("rand");
	f->args.resize(1 + 0);
	f->args[0] = new MathArgInt_t;
	f->f = mp_rand;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// random
	f = new MathFunc_t;
	f->fname = std::string("random");
	f->args.resize(1 + 0);
	f->args[0] = new MathArgReal_t;
	f->f = mp_rndm;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// seed
	f = new MathFunc_t;
	f->fname = std::string("seed");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgVoid_t;
	f->args[1] = new MathArgInt_t;
	f->f = mp_srnd;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// step
	f = new MathFunc_t;
	f->fname = std::string("step");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_step;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// ramp
	f = new MathFunc_t;
	f->fname = std::string("ramp");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_ramp;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// sramp
	f = new MathFunc_t;
	f->fname = std::string("sramp");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->f = mp_sramp;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// par
	f = new MathFunc_t;
	f->fname = std::string("par");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_par;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// in
	f = new MathFunc_t;
	f->fname = std::string("in_ll");
	f->args.resize(1 + 3);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->args[3] = new MathArgReal_t;
	f->f = mp_in<IN_LL>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	f = new MathFunc_t;
	f->fname = std::string("in_le");
	f->args.resize(1 + 3);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->args[3] = new MathArgReal_t;
	f->f = mp_in<IN_LE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	f = new MathFunc_t;
	f->fname = std::string("in_el");
	f->args.resize(1 + 3);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->args[3] = new MathArgReal_t;
	f->f = mp_in<IN_EL>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	f = new MathFunc_t;
	f->fname = std::string("in_ee");
	f->args.resize(1 + 3);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->args[2] = new MathArgReal_t;
	f->args[3] = new MathArgReal_t;
	f->f = mp_in<IN_EE>;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// print
	f = new MathFunc_t;
	f->fname = std::string("print");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgVoid_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_print;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// stop
	f = new MathFunc_t;
	f->fname = std::string("stop");
	f->args.resize(1 + 2);
	f->args[0] = new MathArgVoid_t;
	f->args[1] = new MathArgInt_t;
	f->args[2] = new MathArgInt_t;
	f->f = mp_stop;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	// cast
	{
		struct {
			TypedValue::Type type;
			MathArg_t *arg;
		} data[5] = {
			{ TypedValue::VAR_BOOL },
			{ TypedValue::VAR_INT },
			{ TypedValue::VAR_REAL },
			{ TypedValue::VAR_STRING },
			{ TypedValue::VAR_UNKNOWN }
		};

		data[0].arg = new MathArgBool_t;
		data[1].arg = new MathArgInt_t;
		data[2].arg = new MathArgReal_t;
		data[3].arg = new MathArgString_t;

		for (int i = 0; data[i].type != TypedValue::VAR_UNKNOWN; i++) {
			f = new MathFunc_t;
			f->fname = std::string(TypedValue::GetTypeName(data[i].type));
			f->args.resize(1 + 1);
			f->args[0] = data[i].arg;
			f->args[1] = new MathArgAny_t;
			f->f = mp_cast;
			f->t = 0;

			if (!func.insert(funcType::value_type(f->fname, f)).second) {
				silent_cerr("static namespace: "
					"unable to insert handler "
					"for function " << f->fname << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}
}

MathParser::StaticNameSpace::~StaticNameSpace(void)
{
	for (funcType::iterator f = func.begin(); f != func.end(); ++f) {
		for (MathParser::MathArgs::iterator i = f->second->args.begin();
			i != f->second->args.end();
			++i)
		{
			delete *i;
		}

		delete f->second;
	}
}

bool
MathParser::StaticNameSpace::IsFunc(const std::string& fname) const
{
	if (func.find(fname) != func.end()) {
		return true;
	}

	return false;
}

MathParser::MathFunc_t*
MathParser::StaticNameSpace::GetFunc(const std::string& fname) const
{
	funcType::const_iterator i = func.find(fname);

	if (i != func.end()) {
		return i->second;
	}

	return 0;
}

TypedValue
MathParser::StaticNameSpace::EvalFunc(MathParser::MathFunc_t *f, const MathArgs& args) const
{
	f->f(args);

	switch (args[0]->Type()) {
	case MathParser::AT_VOID:
		return TypedValue(0);

	case MathParser::AT_BOOL:
		return TypedValue((*dynamic_cast<MathArgBool_t*>(args[0]))());

	case MathParser::AT_INT:
		return TypedValue((*dynamic_cast<MathArgInt_t*>(args[0]))());

	case MathParser::AT_REAL:
		return TypedValue((*dynamic_cast<MathArgReal_t*>(args[0]))());

	case MathParser::AT_STRING:
		return TypedValue((*dynamic_cast<MathArgString_t*>(args[0]))());

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

MathParser::DeclarationModifier
MathParser::GetDeclarationModifier(const char* const s) const
{
	for (Int i = 0; DeclarationModifierNames[i].name != NULL; i++) {
		if (strcmp(s, DeclarationModifierNames[i].name) == 0) {
			return DeclarationModifierNames[i].type;
		}
	}

	return DMOD_UNKNOWN;
}

bool
MathParser::IsType(const char* const s) const
{
	return GetType(s) != TypedValue::VAR_UNKNOWN;
}

bool
MathParser::IsTypeModifier(const char* const s) const
{
	return GetTypeModifier(s) != TypedValue::MOD_UNKNOWN;
}

bool
MathParser::IsDeclarationModifier(const char* const s) const
{
	return GetDeclarationModifier(s) != DMOD_UNKNOWN;
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

MathParser::Token
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
			/* a trailing '\' continues the one-line remark
			 * on the following line
			 * FIXME: are we sure we want this? */
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

		s[i++] = char(c);

		if (c == '.') {
			f = true;
		}
		while ((c = in->get()) == '.' || isdigit(c)) {
			s[i++] = char(c);
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
				s[i++] = char(c);
				c = in->get();
			}
			if (isdigit(c)) {
				s[i++] = char(c);
			} else {
			       	return (currtoken = UNKNOWNTOKEN);
			}
			while (isdigit((c = in->get()))) {
				s[i++] = char(c);
			}
		}
		s[i] = '\0';
		in->putback(char(c));
		char *endptr = 0;
		if (!f) {
			value.SetType(TypedValue::VAR_INT);
#ifdef HAVE_STRTOL
			errno = 0;
			long l = strtol(s, &endptr, 10);
			int save_errno = errno;
			if (endptr == s || endptr[0] != '\0') {
				return (currtoken = UNKNOWNTOKEN);
			}

			if (save_errno == ERANGE) {
				// over/under-flow
				if (l == LONG_MIN) {
					throw ErrGeneric(this,
						MBDYN_EXCEPT_ARGS,
						"integer value ",
						std::string(s, endptr - s).c_str(),
						" underflow");
				}

				if (l == LONG_MAX) {
					throw ErrGeneric(this,
						MBDYN_EXCEPT_ARGS,
						"integer value ",
						std::string(s, endptr - s).c_str(),
						" overflow");
				}

				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			value.Set(Int(l));
#else /* !HAVE_STRTOL */
			value.Set(Int(atoi(s)));
#endif /* !HAVE_STRTOL */

		} else {
			value.SetType(TypedValue::VAR_REAL);
#ifdef HAVE_STRTOD
			errno = 0;
			double d = strtod(s, &endptr);
			int save_errno = errno;
			if (endptr == s || endptr[0] != '\0') {
				return (currtoken = UNKNOWNTOKEN);
			}

			if (save_errno == ERANGE) {
				// over/under-flow
				if (std::abs(d) == HUGE_VAL) {
					throw ErrGeneric(this,
						MBDYN_EXCEPT_ARGS,
						"real value ",
						std::string(s, endptr - s).c_str(),
						" overflow");
				}

				if (d == 0.) {
					throw ErrGeneric(this,
						MBDYN_EXCEPT_ARGS,
						"real value ",
						std::string(s, endptr - s).c_str(),
						" underflow");
				}

				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			value.Set(Real(d));
#else /* !HAVE_STRTOD */
			value.Set(Real(atof(s)));
#endif /* !HAVE_STRTOD */
		}

		return (currtoken = NUM);
	}

	/* name? */
	if (isalpha(c) || c == '_') {
		int l = 0;
		namebuf[l++] = char(c);
		while ((c = in->get())) {
			if (!(c == '_'
				|| isalnum(c)
				|| ((currtoken == NAMESPACESEP) && c == ':')))
			{
				break;
			}

			namebuf[l++] = char(c);
			if (l ==  namebuflen) {
				IncNameBuf();
			}
		}
		namebuf[l] = '\0';
		in->putback(char(c));
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
			in->putback(char(c));
			return (currtoken = DIV);
		}

	case '%':
		return (currtoken = MOD);

	case '-':
		return (currtoken = MINUS);

	case '+':
		return (currtoken = PLUS);

	case '>':
		if ((c = in->get()), c == '=') {
			return (currtoken = GE);
		}
		in->putback(char(c));
		return (currtoken = GT);

	case '=':
		if ((c = in->get()), c == '=') {
			return (currtoken = EQ);
		}
		in->putback(char(c));
		return (currtoken = ASSIGN);

	case '<':
		if ((c = in->get()), c == '=') {
			return (currtoken = LE);
		}
		in->putback(char(c));
		return (currtoken = LT);

	case '!':
		if ((c = in->get()), c == '=') {
			return (currtoken = NE);
		}
		in->putback(char(c));
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

	case '"': {
		int l = 0;
		while ((c = in->get()) != '"') {
			if (c == '\\') {
				c = in->get();
				if (c == '\0') {
					return (currtoken = UNKNOWNTOKEN);
				}
				if (c == EOF || in->eof()) {
					return (currtoken = ENDOFFILE);
				}
			}
			namebuf[l++] = char(c);
			if (l == namebuflen) {
				IncNameBuf();
			}
		}
		namebuf[l] = '\0';
		value.SetType(TypedValue::VAR_STRING);
		value.Set(std::string(namebuf));
		return (currtoken = NUM);
		}

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
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "divide by zero in mult()");
			}
			d /= e;
			break;
		}

		case MOD: {
			GetToken();
			TypedValue e = power();
			if (e == 0.) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "divide by zero in mult()");
			}
			d %= e;
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
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "invalid operands in power()");
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
MathParser::evalfunc(MathParser::NameSpace *ns, MathParser::MathFunc_t* f)
{
	MathArgs args(f->args.size());

	for (unsigned i = 0; i < f->args.size(); i++) {
		args[i] = f->args[i]->Copy();
	}

	for (unsigned i = 1; i < args.size(); i++) {
		switch (args[i]->Type()) {
		case MathParser::AT_ANY:
			if (currtoken == CBR) {
				if (!args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "argument expected");
				}
				(*dynamic_cast<MathArgAny_t*>(args[i]))() = (*dynamic_cast<MathArgAny_t*>(f->args[i]))();
				args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgAny_t*>(args[i]))() = stmtlist();
			}
			break;

		case MathParser::AT_BOOL:
			if (currtoken == CBR) {
				if (!args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "bool argument expected");
				}
				(*dynamic_cast<MathArgBool_t*>(args[i]))() = (*dynamic_cast<MathArgBool_t*>(f->args[i]))();
				args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgBool_t*>(args[i]))() = stmtlist().GetBool();
			}
			break;

		case MathParser::AT_INT:
			if (currtoken == CBR) {
				if (!args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "integer argument expected");
				}
				(*dynamic_cast<MathArgInt_t*>(args[i]))() = (*dynamic_cast<MathArgInt_t*>(f->args[i]))();
				args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgInt_t*>(args[i]))() = stmtlist().GetInt();
			}
			break;

		case MathParser::AT_REAL:
			if (currtoken == CBR) {
				if (!args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "real argument expected");
				}
				(*dynamic_cast<MathArgReal_t*>(args[i]))() = (*dynamic_cast<MathArgReal_t*>(f->args[i]))();
				args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgReal_t*>(args[i]))() = stmtlist().GetReal();
			}
			break;

		case MathParser::AT_STRING:
			if (currtoken == CBR) {
				if (!args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "string argument expected");
				}
				(*dynamic_cast<MathArgString_t*>(args[i]))() = (*dynamic_cast<MathArgString_t*>(f->args[i]))();
				args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgString_t*>(args[i]))() = stmtlist().GetString();
			}
			break;

		case MathParser::AT_PRIVATE:
			/* ignore */
			break;

		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i < args.size() - 1) {
			if (args[i + 1]->Type() != AT_PRIVATE) {
				switch (currtoken) {
				case CBR:
					if (!args[i + 1]->IsFlag(MathParser::AF_OPTIONAL)) {
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"mandatory argument expected");
					}
					break;

				case ARGSEP:
					GetToken();
					break;

				default:
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS,
						"argument separator expected");
				}
			}
		}
	}

	if (f->t != 0) {
		if (f->t(args)) {
			DEBUGCERR("error in function "
				<< ns->sGetName() << "::" << f->fname 
				<< " " "(msg: " << f->errmsg << ")"
				<< " in evalfunc()" << std::endl);
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, f->fname.c_str(), ": error ", f->errmsg.c_str());
		}
	}

	TypedValue val = ns->EvalFunc(f, args);

	for (unsigned i = 0; i < args.size(); i++) {
		delete args[i];
	}

	return val;
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
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == OPGIN) {
		TypedValue d = readplugin();
		if (currtoken != CPGIN) {
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing plugin expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == NAME) {
		GetToken();
		MathParser::NameSpace *currNameSpace = defaultNameSpace;
		if (currtoken == NAMESPACESEP) {
			std::string name(namebuf);
			NameSpaceMap::iterator i = nameSpaceMap.find(name);
			if (i == nameSpaceMap.end()) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf, "\"");
			}
			currNameSpace = i->second;
			GetToken();
			if (currtoken != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
			}
			GetToken();
		}

		if (currtoken == OBR) {
			/* in futuro ci potranno essere magari i dati strutturati */
			if (!currNameSpace->IsFunc(namebuf)) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf, "' not found; user-defined functions not supported yet!");
			}

			MathParser::MathFunc_t* f = currNameSpace->GetFunc(namebuf);
			if (f == NULL) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf, "' not found");
			}
			GetToken();
			TypedValue d = evalfunc(currNameSpace, f);
			if (currtoken != CBR) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis "
						"expected after function "
						"\"", f->fname.c_str(), "\" in expr()");
			}
			GetToken();
			return d;
			
		} else {
			NamedValue* v = table.Get(namebuf);
			if (v != NULL) {
				return v->GetVal();
			}
		}

		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unknown name \"", namebuf, "\"");
	}

	/* invalid expr */
	if (currtoken != ENDOFFILE) {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unknown token");
	}

	return TypedValue(0.);
}

TypedValue
MathParser::stmt(void)
{
	if (currtoken == NAME) {
		bool isIfndef = false;
		bool isConst = false;

		DeclarationModifier declarationmodifier = GetDeclarationModifier(namebuf);
		if (declarationmodifier != DMOD_UNKNOWN) {
			switch (declarationmodifier) {
			case DMOD_IFNDEF:
				isIfndef = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled definition modifier ", namebuf, "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type (modifier) expected "
					"after definition modifier in declaration");
			}
		}

		TypedValue::TypeModifier typemodifier = GetTypeModifier(namebuf);
		if (typemodifier != TypedValue::MOD_UNKNOWN) {
			switch (typemodifier) {
			case TypedValue::MOD_CONST:
				isConst = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled type modifier ", namebuf, "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type expected "
					"after type modifier in declaration");
			}
		}

		/* declaration? */
		TypedValue::Type type = GetType(namebuf);
		if (type != TypedValue::VAR_UNKNOWN) {
			switch (GetToken()) {
			case OBR:
				// explicit cast?
				TokenPush(currtoken);
				currtoken = NAME;
				return logical();
				
			case NAME:
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected "
					"after type in declaration");
			}

			/* FIXME: need to specialize symbol table for namespaces */
			if (IsKeyWord(defaultNameSpace, namebuf)) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name '", namebuf, "' "
						"is a keyword");
			}

			/* with assign? */
			if (GetToken() == ASSIGN) {
				/* faccio una copia del nome! */
				char* varname = NULL;
				SAFESTRDUP(varname, namebuf);

				GetToken();
				TypedValue d = logical();
				/* set const'ness of newly gotten value */
				if (isConst) {
					d.SetConst();
				}

				/* cerco la variabile */
				NamedValue* v = table.Get(varname);

				if (v == NULL) {
					/* create new var with assigned type */
					TypedValue newvar(type);
					/* assign new var, so that it internally
					 * takes care of casting, while it inherits
					 * const'ness from d */
					newvar.Cast(d);
					v = table.Put(varname, newvar);

					if (isIfndef) {
						silent_cerr("warning, ifndef variable " << v->GetTypeName() << " \"" << v->GetName()
							<< "\" not yet defined; set to \"" << newvar << "\" at line " << mbdyn_get_line_data() << std::endl);
					}

				} else {
					/* altrimenti, se la posso ridefinire, mi limito
					 * ad assegnarle il nuovo valore */
					if (!bRedefineVars && !isIfndef) {
				   		throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"var \"", v->GetName(), "\"");
					}

					if (v->Const() && !isIfndef) {
						// TODO: check redefinition of const'ness
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"a const named value "
							"\"", v->GetName(), "\"");
					}

					if (!v->IsVar()) {
			   			throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"non-var named value "
							"\"", v->GetName(), "\"");
					}

					if (!isIfndef) {
						dynamic_cast<Var *>(v)->SetVal(d);

					} else {
						if (v->GetType() != type) {
							silent_cerr("warning, skipping redefinition of \"" << v->GetName()
								<< "\" from " << v->GetTypeName() << " to " << TypedValue::GetTypeName(type) << " at line " << mbdyn_get_line_data() << std::endl);

						} else {
							silent_cerr("warning, skipping redefinition of " << v->GetTypeName() << " \"" << v->GetName()
								<< "\" at line " << mbdyn_get_line_data() << std::endl);
						}
					}
				}
				
				/* distruggo il temporaneo */
				SAFEDELETEARR(varname);
				return v->GetVal();

			} else if (currtoken == STMTSEP) {
				NamedValue* v = table.Get(namebuf);
				if (v == NULL || !bRedefineVars) {
					if (isConst) {
						/* cannot insert a const var
						 * with no value */
			      			throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
		      						"cannot create const named value "
								"\"", namebuf, "\" with no value");
					}
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
			if (declarationmodifier != DMOD_UNKNOWN) {
			      	throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
		      			"definition modifier without type");
			}

			if (typemodifier != TypedValue::MOD_UNKNOWN) {
			      	throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
		      			"type modifier without type");
			}

			/* assignment? */
			NamedValue* v = table.Get(namebuf);
			if (v != NULL) {
				if (GetToken() == ASSIGN) {
					GetToken();
					TypedValue d = logical();
					if (v->Const()) {
			      			throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
		      						"cannot assign const named value "
								"\"", v->GetName(), "\"");
			 		}

			 		if (!v->IsVar()) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
								"cannot assign non-var named value "
								"\"", v->GetName(), "\"");
			 		}
					dynamic_cast<Var *>(v)->Cast(d);
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
	std::vector<char *> argv(1);
	char c, buf[1024];
	int argc = 0;
	unsigned int i = 0, in_quotes = 0;

	/*
	 * inizializzo l'array degli argomenti
	 */
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
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
			argv.resize(argc + 2);
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
			silent_cerr("MathParser::readplugin(): buffer overflow" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

last_arg:
	if (in->eof()) {
		silent_cerr("eof encountered while parsing plugin"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (varname == NULL || *varname == '\0') {
		silent_cerr("illegal or missing plugin variable name"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	/*
	 * verifica esistenza nome
	 */
	NamedValue* v = table.Get(varname);
	if (v != NULL) {
		silent_cerr("variable " << varname << " already defined"
			<< std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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
		pgin->Read(argc - 2, &argv[2]);

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

		return v->GetVal();
	}

	/*
	 * si arriva qui solo se il plugin non e' stato registrato
	 */
	silent_cerr("plugin '" << pginname << "' not supported" << std::endl);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

MathParser::MathParser(const InputStream& strm, Table& t, bool bRedefineVars)
: PlugIns(0),
table(t),
bRedefineVars(bRedefineVars),
in(const_cast<InputStream*>(&strm)),
defaultNameSpace(0),
namebuf(0),
namebuflen(default_namebuflen),
value(),
currtoken(UNKNOWNTOKEN),
tokenlist(NULL)
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	defaultNameSpace = new StaticNameSpace();
	if (RegisterNameSpace(defaultNameSpace)) {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to register namespace "
				"\"", defaultNameSpace->sGetName().c_str(), "\"");
	}

	time_t tm;
	time(&tm);
	srand(tm);
}

MathParser::MathParser(Table& t, bool bRedefineVars)
: PlugIns(0),
table(t),
bRedefineVars(bRedefineVars),
in(0),
defaultNameSpace(0),
namebuf(0),
namebuflen(default_namebuflen),
value(),
currtoken(UNKNOWNTOKEN),
tokenlist(0)
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	defaultNameSpace = new StaticNameSpace();
	if (RegisterNameSpace(defaultNameSpace)) {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to register namespace "
				"\"", defaultNameSpace->sGetName().c_str(), "\"");
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
				throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
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
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "error while adding real var '", s, "");
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
				throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
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
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "error while adding integer var '", s, "");
	}

	return var;
}

MathParser::~MathParser(void)
{
	DEBUGCOUTFNAME("MathParser::~MathParser");

	if (namebuf != NULL) {
		SAFEDELETEARR(namebuf);
	}

	while (PlugIns) {
		PlugInRegister *next = PlugIns->next;
		SAFEDELETEARR(PlugIns->name);
		SAFEDELETE(PlugIns);
		PlugIns = next;
	}

	for (NameSpaceMap::iterator i = nameSpaceMap.begin(); i != nameSpaceMap.end(); ++i) {
		delete i->second;
	}
}

Real
MathParser::GetLastStmt(Real d, Token t)
{
	if (GetToken() == t) {
		return d;
	}

	for (;;) {
		d = stmtlist().GetReal();
		if (currtoken == ENDOFFILE || currtoken == t) {
			break;
		}
		if (currtoken != STMTSEP) {
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "statement separator expected");
		}
	}

	return d;
}

Real
MathParser::GetLastStmt(const InputStream& strm, Real d, Token t)
{
	const InputStream* save_in = in;
	Token save_currtoken = currtoken;

	in = const_cast<InputStream *>(&strm);
	currtoken = UNKNOWNTOKEN;

	d = GetLastStmt(d, t);

	in = const_cast<InputStream *>(save_in);
	currtoken = save_currtoken;

	return d;
}

Real
MathParser::Get(Real d)
{
	TypedValue v(d);
	v = Get(v);
	return v.GetReal();
}

Real
MathParser::Get(const InputStream& strm, Real d)
{
	TypedValue v(d);
	v = Get(strm, v);
	return v.GetReal();
}

TypedValue
MathParser::Get(const TypedValue& /* v */ )
{
	GetToken();
	TypedValue vv = stmt();
	if (currtoken != STMTSEP && currtoken != ENDOFFILE) {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "statement separator expected");
	}
	return vv;
}

TypedValue
MathParser::Get(const InputStream& strm, const TypedValue& v)
{
	const InputStream* p = in;
	in = (InputStream*)&strm;
	GetToken();
	TypedValue vv = v;
	if (currtoken != STMTSEP && currtoken != ARGSEP) {
		vv = stmt();
	}
	if (currtoken == STMTSEP) {
		in->putback(';');

	} else if (currtoken == ARGSEP) {
		in->putback(',');

	} else {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "separator expected");
	}
	in = (InputStream*)p;

	return vv;
}

void
MathParser::GetForever(std::ostream& out, const char* const sep)
{
	do {
		TypedValue val(0.);
		out << Get(val) << sep;
	} while (currtoken == STMTSEP);
}

void
MathParser::GetForever(const InputStream& strm, std::ostream& out,
		const char* const /* sep */ )
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
	pedantic_cout("registering plugin \"" << name << "\"" << std::endl);

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

	pedantic_cout("MathParser::RegisterNameSpace: "
		"registering namespace \"" << ns->sGetName() << "\""
		<< std::endl);

	if (nameSpaceMap.find(ns->sGetName()) != nameSpaceMap.end()) {
		return 1;
	}

	nameSpaceMap[ns->sGetName()] = ns;

	return 0;
}

MathParser::NameSpace *
MathParser::GetNameSpace(const std::string& name) const
{
	NameSpaceMap::const_iterator i = nameSpaceMap.find(name);
	if (i == nameSpaceMap.end()) {
		return 0;
	}

	return i->second;
}

