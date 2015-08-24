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

/*
 * With the contribution of Ankit Aggarwal <ankit.ankit.aggarwal@gmail.com>
 * during Google Summer of Code 2015
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cerrno>
#include <cfloat>
#include <cstdlib>
#include <climits>
#include <limits>
#include <sstream>

#include "mathp.h"
#include "parser.h"

#ifdef USE_EE
#include "evaluator_impl.h"
#endif // USE_EE


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

	switch (args[0]->Type()) {
	case MathParser::AT_BOOL: {
		MathParser::MathArgBool_t *out = dynamic_cast<MathParser::MathArgBool_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetBool();
		} break;

	case MathParser::AT_INT: {
		MathParser::MathArgInt_t *out = dynamic_cast<MathParser::MathArgInt_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetInt();
		} break;

	case MathParser::AT_REAL: {
		MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
		ASSERT(out != 0);
		(*out)() = (*arg1)().GetReal();
		} break;

	case MathParser::AT_STRING: {
		MathParser::MathArgString_t *out = dynamic_cast<MathParser::MathArgString_t*>(args[0]);
		ASSERT(out != 0);
		TypedValue val(TypedValue::VAR_STRING);
		val.Cast((*arg1)());
		(*out)() = val.GetString();
		} break;

	default:
		// add support for future types
		return 1;
	}

	return 0;
}

/* tipi delle variabili */
struct TypeName_t {
	const char* name;
	TypedValue::Type type;
};

static const TypeName_t TypeNames[] = {
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

static const typemodifiernames TypeModifierNames[] = {
	{ "const",	TypedValue::MOD_CONST },
	{ NULL,		TypedValue::MOD_UNKNOWN }
};

struct declarationmodifiernames {
	const char* name;
	MathParser::DeclarationModifier type;
};

static const declarationmodifiernames DeclarationModifierNames[] = {
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

TypedValue&
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

TypedValue&
TypedValue::Cast(const TypedValue& var, bool bErr)
{
	if (Const()) {
		throw ErrConstraintViolation(MBDYN_EXCEPT_ARGS);
	}

	if (type == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (var.type) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			// TODO: use std::ostringstream?
			snprintf(buf, sizeof(buf), "%ld", (long)var.GetInt());
			Set(std::string(buf));
			break;

		case TypedValue::VAR_REAL:
			// TODO: use std::ostringstream?
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
				if (bErr) {
					Real r = var.GetReal();
					if (r != 0. && r != 1.) {
						silent_cerr(" Error: implicit cast "
							"of " << var
							<< " from " << GetTypeName(var.type)
							<< " to " << GetTypeName(type)
							<< " alters its value"
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

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
				if (bErr) {
					if (var.GetReal() != Real(var.GetInt())) {
						silent_cerr(" Error: implicit cast "
							"of " << var
							<< " from " << GetTypeName(var.type)
							<< " to " << GetTypeName(type)
							<< " alters its value"
							<< std::endl);
						throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					}
				}

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

	case TypedValue::VAR_REAL: {
		// FIXME: precision?
		// make sure there is a (trailing) '.'
		std::ostringstream os;
		os << v.GetReal();
		if (os.str().find('.') == std::string::npos) {
			os << '.';
		}
		return out << os.str();
		}

	case TypedValue::VAR_STRING:
		return out << '"' << v.GetString() << '"';

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

bool
Var::MayChange(void) const
{
	return !value.Const();
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
Var::Cast(const TypedValue& v, bool bErr)
{
	value.Cast(v, bErr);
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

bool
MathParser::PlugInVar::MayChange(void) const
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
	silent_cerr(" MathParser - " << r << s2 << s3
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

void
MathParser::TokenPush(Token t)
{
	TokenVal tv;
	tv.m_t = t;
	if (t == NUM) {
		tv.m_v = value;
	} else if (t == NAME) {
		tv.m_v = TypedValue(namebuf);
	}
	TokenStack.push(tv);
}

int
MathParser::TokenPop(void)
{
	if (TokenStack.empty()) {
		return 0;
	}

	TokenVal tv = TokenStack.top();
	TokenStack.pop();

	currtoken = tv.m_t;
	if (currtoken == NUM) {
		value = tv.m_v;
	} else if (currtoken == NAME) {
		// note: namebuf is expected to be large enough for value
		namebuf = tv.m_v.GetString();
	}

	return 1;
}

const MathParser::NameSpaceMap&
MathParser::GetNameSpaceMap(void) const
{
	return nameSpaceMap;
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

static bool sns = 0;

MathParser::StaticNameSpace::StaticNameSpace(Table *pTable)
: MathParser::NameSpace("default"), m_pTable(pTable)
{
	// make sure there's only one StaticNameSpace
	if (sns) {
		// error
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	sns = true;

	MathFunc_t	*f;

	// asin
	f = new MathFunc_t;
	f->fname = std::string("asin");
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
	f->ns = this;
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
			f->ns = this;
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
		delete f->second;
	}
}

bool
MathParser::StaticNameSpace::IsFunc(const std::string& fname) const
{
#if 0
	for (funcType::const_iterator f = func.begin(); f != func.end(); ++f) {
		silent_cerr("*** " << f->second->fname << std::endl);
	}
#endif

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
		return new MathParser::MathFunc_t(*i->second);
	}

	return 0;
}

TypedValue
MathParser::StaticNameSpace::EvalFunc(MathParser::MathFunc_t *f) const
{
	f->f(f->args);

	switch (f->args[0]->Type()) {
	case MathParser::AT_VOID:
		return TypedValue(0);

	case MathParser::AT_BOOL:
		return TypedValue((*dynamic_cast<MathArgBool_t*>(f->args[0]))());

	case MathParser::AT_INT:
		return TypedValue((*dynamic_cast<MathArgInt_t*>(f->args[0]))());

	case MathParser::AT_REAL:
		return TypedValue((*dynamic_cast<MathArgReal_t*>(f->args[0]))());

	case MathParser::AT_STRING:
		return TypedValue((*dynamic_cast<MathArgString_t*>(f->args[0]))());

	default:
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

Table*
MathParser::StaticNameSpace::GetTable(void)
{
	return m_pTable;
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
		// lot of space...
		char s[BUFSIZ];
		bool f = false;
		unsigned i = 0;

		// FIXME: need to check for overflow

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
			if (i >= sizeof(s)) {
				// buffer about to overflow
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "value too long");
			}
		}
		if (std::strchr("efdgEFDG", c) != 0) {
			f = true;
			// use 'e' because strtod only understands 'e' or 'E'
			s[i++] = 'e';
			if (i >= sizeof(s)) {
				// buffer about to overflow
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "value too long");
			}
			if ((c = in->get()) == '-' || c == '+') {
				s[i++] = char(c);
				if (i >= sizeof(s)) {
					// buffer about to overflow
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "value too long");
				}
				c = in->get();
			}
			if (isdigit(c)) {
				s[i++] = char(c);
				if (i >= sizeof(s)) {
					// buffer about to overflow
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "value too long");
				}

			} else {
				return (currtoken = UNKNOWNTOKEN);
			}
			while (isdigit((c = in->get()))) {
				s[i++] = char(c);
				if (i >= sizeof(s)) {
					// buffer about to overflow
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "value too long");
				}
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
				silent_cerr(" MathParser - unable to parse \"" << s << "\" as integer"
					<< " at line " << GetLineNumber() << std::endl);
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
				silent_cerr(" MathParser - unable to parse \"" << s << "\" as real"
					<< " at line " << GetLineNumber() << std::endl);
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
		namebuf.clear();
		namebuf.push_back(char(c));
		while ((c = in->get())) {
			if (!(c == '_'
				|| isalnum(c)
				|| ((currtoken == NAMESPACESEP) && c == ':')))
			{
				break;
			}

			namebuf.push_back(char(c));
		}
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
		namebuf.clear();
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
			namebuf.push_back(char(c));
		}
		value.SetType(TypedValue::VAR_STRING);
		value.Set(namebuf);
		return (currtoken = NUM);
		}

	default:
		return (currtoken = UNKNOWNTOKEN);
	}
}

bool
MathParser::bNameValidate(const std::string& s) const
{
	std::string::const_iterator i = s.begin();
	if (*i != '_' && !isalpha(*i)) {
		return false;
	}

	while (++i != s.end()) {
		if (*i != '_' && !isalnum(*i)) {
			return false;
		}
	}

	return true;
}

#ifndef USE_EE
TypedValue
MathParser::logical(void)
{
	return logical_int(relational());
}

TypedValue
MathParser::logical(TypedValue d)
{
	return logical_int(relational(d));
}

TypedValue
MathParser::logical_int(TypedValue d)
{
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
	return relational_int(binary());
}

TypedValue
MathParser::relational(TypedValue d)
{
	return relational_int(binary(d));
}

TypedValue
MathParser::relational_int(TypedValue d)
{
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
	return binary_int(mult());
}

TypedValue
MathParser::binary(TypedValue d)
{
	return binary_int(mult(d));
}

TypedValue
MathParser::binary_int(TypedValue d)
{
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
	return mult_int(power());
}

TypedValue
MathParser::mult(TypedValue d)
{
	return mult_int(power(d));
}

TypedValue
MathParser::mult_int(TypedValue d)
{
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
	return power_int(unary());
}

TypedValue
MathParser::power(TypedValue d)
{
	return power_int(d);
}

TypedValue
MathParser::power_int(TypedValue d)
{
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
			DEBUGCERR("can't compute (" << d << ")^("
					<< e << ") in power()" << std::endl);
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
			d = TypedValue(Real(std::pow(b, r)));
		}
	}

	return d;
}


TypedValue
MathParser::unary(void)
{
	switch (currtoken) {
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
	for (unsigned i = 1; i < f->args.size(); i++) {
		switch (f->args[i]->Type()) {
		case MathParser::AT_ANY:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "argument expected");
				}
				(*dynamic_cast<MathArgAny_t*>(f->args[i]))() = (*dynamic_cast<MathArgAny_t*>(f->args[i]))();
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgAny_t*>(f->args[i]))() = stmtlist();
			}
			break;

		case MathParser::AT_BOOL:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "bool argument expected");
				}
				(*dynamic_cast<MathArgBool_t*>(f->args[i]))() = (*dynamic_cast<MathArgBool_t*>(f->args[i]))();
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgBool_t*>(f->args[i]))() = stmtlist().GetBool();
			}
			break;

		case MathParser::AT_INT:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "integer argument expected");
				}
				(*dynamic_cast<MathArgInt_t*>(f->args[i]))() = (*dynamic_cast<MathArgInt_t*>(f->args[i]))();
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgInt_t*>(f->args[i]))() = stmtlist().GetInt();
			}
			break;

		case MathParser::AT_REAL:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "real argument expected");
				}
				(*dynamic_cast<MathArgReal_t*>(f->args[i]))() = (*dynamic_cast<MathArgReal_t*>(f->args[i]))();
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgReal_t*>(f->args[i]))() = stmtlist().GetReal();
			}
			break;

		case MathParser::AT_STRING:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "string argument expected");
				}
				(*dynamic_cast<MathArgString_t*>(f->args[i]))() = (*dynamic_cast<MathArgString_t*>(f->args[i]))();
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				(*dynamic_cast<MathArgString_t*>(f->args[i]))() = stmtlist().GetString();
			}
			break;

		case MathParser::AT_PRIVATE:
			/* ignore */
			break;

		default:
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i < f->args.size() - 1) {
			if (f->args[i + 1]->Type() != AT_PRIVATE) {
				switch (currtoken) {
				case CBR:
					if (!f->args[i + 1]->IsFlag(MathParser::AF_OPTIONAL)) {
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
		if (f->t(f->args)) {
			DEBUGCERR("error in function "
				<< ns->sGetName() << "::" << f->fname 
				<< " " "(msg: " << f->errmsg << ")"
				<< " in evalfunc()" << std::endl);
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, f->fname.c_str(), ": error ", f->errmsg.c_str());
		}
	}

	TypedValue val = ns->EvalFunc(f);

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
		std::string name(namebuf);
		MathParser::NameSpace *currNameSpace = defaultNameSpace;
		Table *currTable = &table;

		GetToken();
		if (currtoken == NAMESPACESEP) {
			NameSpaceMap::iterator i = nameSpaceMap.find(name);
			if (i == nameSpaceMap.end()) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
			}
			currNameSpace = i->second;
			currTable = currNameSpace->GetTable();
			GetToken();
			if (currtoken != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
			}
			name += "::";
			name += namebuf;
			GetToken();
		}

		if (currtoken == OBR) {
			/* in futuro ci potranno essere magari i dati strutturati */
			MathParser::MathFunc_t* f = currNameSpace->GetFunc(namebuf);
			if (f == NULL) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf.c_str(), "' not found");
			}
			GetToken();
			TypedValue d = evalfunc(currNameSpace, f);
			delete f;
			if (currtoken != CBR) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis "
						"expected after function "
						"\"", f->fname.c_str(), "\" in expr()");
			}
			GetToken();
			return d;
			
		} else {
			NamedValue* v = 0;
			if (currTable) {
				v = currTable->Get(namebuf);
			}

			if (v != NULL) {
				return v->GetVal();
			}
		}

		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unknown name \"", name.c_str(), "\"");
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

		DeclarationModifier declarationmodifier = GetDeclarationModifier(namebuf.c_str());
		if (declarationmodifier != DMOD_UNKNOWN) {
			switch (declarationmodifier) {
			case DMOD_IFNDEF:
				isIfndef = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled definition modifier ", namebuf.c_str(), "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type (modifier) expected "
					"after definition modifier in declaration");
			}
		}

		TypedValue::TypeModifier typemodifier = GetTypeModifier(namebuf.c_str());
		if (typemodifier != TypedValue::MOD_UNKNOWN) {
			switch (typemodifier) {
			case TypedValue::MOD_CONST:
				isConst = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled type modifier ", namebuf.c_str(), "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type expected "
					"after type modifier in declaration");
			}
		}

		/* declaration? */
		TypedValue::Type type = GetType(namebuf.c_str());
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

			if (IsKeyWord(defaultNameSpace, namebuf.c_str())) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name '", namebuf.c_str(), "' "
						"is a keyword");
			}

			MathParser::NameSpace *currNameSpace = defaultNameSpace;
			Table *currTable = &table;
			std::string name(namebuf);

			GetToken();
			if (currtoken == NAMESPACESEP) {
				NameSpaceMap::iterator i = nameSpaceMap.find(name);
				if (i == nameSpaceMap.end()) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
				}
				currNameSpace = i->second;
				currTable = currNameSpace->GetTable();
				GetToken();
				if (currtoken != NAME) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
				}
				name += "::";
				name += namebuf;
				GetToken();
			}

			/* with assign? */
			if (currtoken == ASSIGN) {
				if (currTable == 0) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
				}

				/* faccio una copia del nome! */
				std::string varname(namebuf);

				GetToken();
				TypedValue d = logical();
				/* set const'ness of newly gotten value */
				if (isConst) {
					d.SetConst();
				}
				/* cerco la variabile */
				NamedValue* v = currTable->Get(varname.c_str());

				if (v == NULL) {
					/* create new var with assigned type */
					TypedValue newvar(type);
					/* assign new var, so that it internally
					 * takes care of casting, while it inherits
					 * const'ness from d */
					newvar.Cast(d);
					v = currTable->Put(varname.c_str(), newvar);

					if (isIfndef) {
						silent_cerr("warning, ifndef variable " << v->GetTypeName() << " \"" << name
							<< "\" not yet defined; set to \"" << newvar << "\" at line " << mbdyn_get_line_data() << std::endl);
					}

				} else {
					/* altrimenti, se la posso ridefinire, mi limito
					 * ad assegnarle il nuovo valore */
					if (!bRedefineVars && !isIfndef) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"var \"", name.c_str(), "\"");
					}

					if (v->Const() && !isIfndef) {
						// TODO: check redefinition of const'ness
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"a const named value "
							"\"", name.c_str(), "\"");
					}

					if (!v->IsVar()) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"non-var named value "
							"\"", name.c_str(), "\"");
					}

					if (!isIfndef) {
						dynamic_cast<Var *>(v)->SetVal(d);

					} else {
						if (v->GetType() != type) {
							silent_cerr("warning, skipping redefinition of \"" << name.c_str() << "\""
								<< " from " << v->GetTypeName() << " to " << TypedValue::GetTypeName(type)
								<< " (orig=" << v->GetVal() << ", unchanged; new=" << d << ")"
								<< " at line " << mbdyn_get_line_data() << std::endl);

						} else {
							silent_cerr("warning, skipping redefinition of " << v->GetTypeName() << " \"" << name.c_str() << "\""
								<< " (orig=" << v->GetVal() << ", unchanged; new=" << d << ")"
								<< " at line " << mbdyn_get_line_data() << std::endl);
						}
					}
				}
				return v->GetVal();

			} else if (currtoken == STMTSEP) {
				if (currTable == 0) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
				}
				
				NamedValue* v = currTable->Get(namebuf);
				if (v == NULL || (!bRedefineVars && !isIfndef)) {
					if (isConst) {
						/* cannot insert a const var
						 * with no value */
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot create const named value "
							"\"", namebuf.c_str(), "\" with no value");
					}
					/* se la var non esiste, la inserisco;
					 * se invece esiste e non vale 
					 * la ridefinizione, tento
					 * di inserirla comunque, cosi'
					 * table da' errore */
					v = currTable->Put(namebuf, TypedValue(type));
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

			MathParser::NameSpace *currNameSpace = defaultNameSpace;
			Table *currTable = &table;
			std::string name(namebuf);

			GetToken();
			if (currtoken == NAMESPACESEP) {
				NameSpaceMap::iterator i = nameSpaceMap.find(name);
				if (i == nameSpaceMap.end()) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
				}
				currNameSpace = i->second;
				currTable = currNameSpace->GetTable();
				GetToken();
				if (currtoken != NAME) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
				}
				name += "::";
				name += namebuf;
				GetToken();
			}

			if (currtoken == OBR) {
				/* in futuro ci potranno essere magari i dati strutturati */
				MathParser::MathFunc_t* f = currNameSpace->GetFunc(namebuf);
				if (f == NULL) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf.c_str(), "' not found");
				}
				GetToken();
				TypedValue d = evalfunc(currNameSpace, f);
				delete f;
				if (currtoken != CBR) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis "
							"expected after function "
							"\"", f->fname.c_str(), "\" in expr()");
				}
				GetToken();

				return logical(d);
			}

			/* assignment? */
			if (currTable == 0) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
			}
			NamedValue* v = currTable->Get(namebuf);
			if (v == NULL) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "var \"", namebuf.c_str(), "\" not found");
			}

			if (currtoken == ASSIGN) {
				GetToken();
				TypedValue d = logical();
				if (v->Const()) {
					throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
						"cannot assign const named value "
						"\"", name.c_str(), "\"");
				}

				if (!v->IsVar()) {
					throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot assign non-var named value "
							"\"", name.c_str(), "\"");
				}
				dynamic_cast<Var *>(v)->Cast(d);
				return v->GetVal();

			} else {
				// NOTE: fails if <name> is actually <namespace>::<name>
				// ASSERT(currtoken != NAME);
				// TokenPush(currtoken);
				// currtoken = NAME;

				// NOTE: fails if <name> is not the complete <stmt>
				// return v->GetVal();

				return logical(v->GetVal());
			}
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
	silent_cerr(" plugin '" << pginname << "' not supported" << std::endl);
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

#endif // ! USE_EE


MathParser::MathParser(const InputStream& strm, Table& t, bool bRedefineVars)
: PlugIns(0),
table(t),
bRedefineVars(bRedefineVars),
in(const_cast<InputStream*>(&strm)),
defaultNameSpace(0),
value(),
currtoken(UNKNOWNTOKEN)
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	// namebuf.resize(4);

	defaultNameSpace = new StaticNameSpace(&t);
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
value(),
currtoken(UNKNOWNTOKEN)
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	// namebuf.resize(4);

	defaultNameSpace = new StaticNameSpace(&t);
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

	if (var == 0) {
		/* Se non c'e' la inserisco */
		var = table.Put(s, TypedValue(v));
	} else {
		/* altrimenti, se la posso ridefinire, mi limito
		 * ad assegnarle il nuovo valore */
		if (redefine) {
			if (var->IsVar()) {
				dynamic_cast<Var *>(var)->SetVal(TypedValue(v));
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
				dynamic_cast<Var *>(var)->SetVal(TypedValue(v));

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
#ifdef USE_EE
		ExpressionElement *e = stmtlist();
		d = e->Eval().GetReal();
		if (pedantic_out) {
			std::cout << "GetLastStmt: \"", e->Output(std::cout) << "\" = " << d << std::endl;
		}
		delete e;
#else // ! USE_EE
		d = stmtlist().GetReal();
#endif // ! USE_EE
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

#ifdef USE_EE
ExpressionElement *
MathParser::GetExpr(void)
{
	if (GetToken() == ARGSEP) {
		return new EE_Value(0.);
	}

	ExpressionElement *e = 0;
	for (;;) {
		e = stmtlist();
		if (pedantic_out) {
			std::cout << "GetExpr: \"", e->Output(std::cout) << "\"" << std::endl;
		}

		if (currtoken == ENDOFFILE || currtoken == ARGSEP) {
			break;
		}

		if (currtoken != STMTSEP) {
			delete e;
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "statement separator expected");
		}
	}

	// NOTE: caller must explicitly delete it
	return e;
}

ExpressionElement *
MathParser::GetExpr(const InputStream& strm)
{
	const InputStream *save_in = in;
	Token save_currtoken = currtoken;

	in = const_cast<InputStream *>(&strm);
	currtoken = UNKNOWNTOKEN;

	ExpressionElement *e = GetExpr();

	in = const_cast<InputStream *>(save_in);
	currtoken = save_currtoken;

	// NOTE: caller must explicitly delete it
	return e;
}
#endif // USE_EE

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
#ifdef USE_EE
	ExpressionElement *e = stmt();
	TypedValue vv = e->Eval();
	if (pedantic_out) {
		std::cout << "Get: \"", e->Output(std::cout) << "\" = " << vv << std::endl;
	}
	delete e;
#else // ! USE_EE
	TypedValue vv = stmt();
#endif // ! USE_EE
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
#ifdef USE_EE
		ExpressionElement *e = stmt();
		vv = e->Eval();
		if (pedantic_out) {
			std::cout << "Get: \"", e->Output(std::cout) << "\" = " << vv << std::endl;
		}
		delete e;
#else // ! USE_EE
		vv = stmt();
#endif // ! USE_EE
	}
	if (currtoken == STMTSEP) {
		in->putback(';');

	} else if (currtoken == ARGSEP) {
		in->putback(',');

	} else {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "separator expected");
	}
	in = const_cast<InputStream *>(p);

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
	InputStream *p = in;
	in = const_cast<InputStream *>(&strm);
	GetForever(out);
	in = p;
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

#ifdef USE_EE
/*
 * Evaluator code 
 * Each function return ExpressionElement * for constructing a tree
 */
ExpressionElement *
MathParser::logical(void)
{
	return logical_int(relational());
}

ExpressionElement *
MathParser::logical(ExpressionElement * d)
{
	return logical_int(relational(d));
}

ExpressionElement *
MathParser::logical_int(ExpressionElement * d)
{
	while (true) {
		switch (currtoken) {
		case AND:
			GetToken();
			// d = new EE_AND(d , relational()); 
			d = EECreate<EE_AND>(d , relational()); 
			break;

		case OR:
			GetToken();
			// d = new EE_OR(d, relational()); 
			d = EECreate<EE_OR>(d , relational()); 
			break;

		case XOR:
			GetToken();
			// d = new EE_XOR(d, relational());
			d = EECreate<EE_XOR>(d , relational()); 
			break;

		default:
			return d;
		}
	}
}

ExpressionElement *
MathParser::relational(void)
{
	return relational_int(binary());
}

ExpressionElement *
MathParser::relational(ExpressionElement * d)
{
	return relational_int(binary(d));
}

ExpressionElement *
MathParser::relational_int(ExpressionElement *  d)
{
	while (true) {
		switch (currtoken) {
		case GT:
			GetToken();
			// d = new EE_Greater(d, binary());
			d = EECreate<EE_Greater>(d , binary()); 
			break;

		case GE:
			GetToken();
			// d = new EE_Greater_Equal(d, binary());
			d = EECreate<EE_Greater_Equal>(d , binary()); 
			break;

		case EQ:
			GetToken();
			// d = new EE_Equal_Equal(d, binary());
			d = EECreate<EE_Equal_Equal>(d , binary()); 
			break;

		case LE:
			GetToken();
			// d = new EE_Lesser_Equal(d, binary());
			d = EECreate<EE_Lesser_Equal>(d , binary()); 
			break;

		case LT:
			GetToken();
			// d = new EE_Lesser(d, binary());
			d = EECreate<EE_Lesser>(d , binary()); 
			break;

		case NE:
			GetToken();
			// d = new EE_Not_Equal(d, binary());
			d = EECreate<EE_Not_Equal>(d , binary()); 
			break;

		default:
			return d;
		}
	}
}

ExpressionElement *
MathParser::binary(void)
{
	return binary_int(mult());
}

ExpressionElement *
MathParser::binary(ExpressionElement * d)
{
	return binary_int(mult(d));
}

ExpressionElement *
MathParser::binary_int(ExpressionElement * d)
{
	while (true) {
		switch (currtoken) {
		case PLUS:
			GetToken();
			// d = new EE_Plus(d, mult());
			d = EECreate<EE_Plus>(d , mult()); 
			break;

		case MINUS:
			GetToken();
			// d = new EE_Minus(d, mult());
			d = EECreate<EE_Minus>(d , mult()); 
			break;

		default:
			return d;
		}
	}
}

ExpressionElement *
MathParser::mult(void)
{
	return mult_int(power());
}

ExpressionElement *
MathParser::mult(ExpressionElement * d)
{
	return mult_int(power(d));
}

ExpressionElement *
MathParser::mult_int(ExpressionElement * d)
{
	while (true) {
		switch (currtoken) {
		case MULT:
			GetToken();
			// d = new EE_Multiply(d, power());
			d = EECreate<EE_Multiply>(d , power()); 
			break;
		
		case DIV:
			GetToken();
			// d = new EE_Divide(d, power());
			d = EECreate<EE_Divide>(d , power()); 
			break;

		case MOD:
			GetToken();
			// d = new EE_Modulus(d, power());
			d = EECreate<EE_Modulus>(d , power()); 
			break;

		default:
			return d;
		}
	}
}


ExpressionElement *
MathParser::power(void)
{
	return power_int(unary());
}

ExpressionElement *
MathParser::power(ExpressionElement* d)
{
	return power_int(d);
}

ExpressionElement *
MathParser::power_int(ExpressionElement * d)
{
	if (currtoken == EXP) {
		GetToken();

		/*
		 * Per l'esponente chiamo di nuovo power cosi' richiama unary;
		 * se per caso dopo unary c'e' di nuovo un esponente,
		 * l'associazione avviene correttamente da destra:
		 *
		 * 	d^e1^e2 == d^(e1^e2)
		 */

		// d = new EE_Power(d, power());
		d = EECreate<EE_Power>(d , power()); 
	}

	return d;
}


ExpressionElement *
MathParser::unary(void)
{
	switch (currtoken) {
	case MINUS:
		GetToken();
		// return new EE_Unary_minus(expr());
		return EECreate<EE_Unary_minus>(expr()); 

	case PLUS:
		GetToken();
		return expr();

	case NOT:
		GetToken();
		// return new EE_NOT(expr());
		return EECreate<EE_NOT>(expr()); 

	default:
		return expr();
	}
}

ExpressionElement *
MathParser::parsefunc(MathFunc_t* f)
{
	for (unsigned i = 1; i < f->args.size(); i++) {
		switch (f->args[i]->Type()) {
		case MathParser::AT_ANY:
		case MathParser::AT_BOOL:
		case MathParser::AT_INT:
		case MathParser::AT_REAL:
		case MathParser::AT_STRING:
			if (currtoken == CBR) {
				if (!f->args[i]->IsFlag(MathParser::AF_OPTIONAL)) {
					switch (f->args[i]->Type()) {
					case MathParser::AT_ANY:
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "argument expected");
					case MathParser::AT_BOOL:
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "bool argument expected");
					case MathParser::AT_INT:
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "integer argument expected");
					case MathParser::AT_REAL:
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "real argument expected");
					case MathParser::AT_STRING:
						throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "string argument expected");
					}
				}
				f->args[i]->SetFlag(MathParser::AF_OPTIONAL_NON_PRESENT);

			} else {
				ExpressionElement *ee = stmtlist();
				f->args[i]->SetExpr(ee);
				if (dynamic_cast<EE_Value *>(ee)) {
					f->args[i]->SetFlag(AF_CONST);
				}
			}
			break;

		case MathParser::AT_PRIVATE:
			/* ignore */
			break;

		default:
			/* NOTE: will need to deal with future types */
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if (i < f->args.size() - 1) {
			if (f->args[i + 1]->Type() != MathParser::AT_PRIVATE) {
				switch (currtoken) {
				case CBR:
					if (!f->args[i + 1]->IsFlag(MathParser::AF_OPTIONAL)) {
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

	return new EE_Func(this, f);
}

ExpressionElement *
MathParser::expr(void)
{
	if (currtoken == NUM) {
		GetToken();
		return new EE_Value(value);
	}

	if (currtoken == OBR) {
		GetToken();
		ExpressionElement *d = stmtlist();
		if (currtoken != CBR) {
			delete d;
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == OPGIN) {
		ExpressionElement *d = readplugin();
		if (currtoken != CPGIN) {
			delete d;
			throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing plugin expected");
		}
		GetToken();
		return d;
	}

	if (currtoken == NAME) {
		std::string name(namebuf);
		MathParser::NameSpace *currNameSpace = defaultNameSpace;
		Table *currTable = &table;

		GetToken();
		if (currtoken == NAMESPACESEP) {
			NameSpaceMap::iterator i = nameSpaceMap.find(name);
			if (i == nameSpaceMap.end()) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
			}
			currNameSpace = i->second;
			currTable = currNameSpace->GetTable();
			GetToken();
			if (currtoken != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
			}
			name += "::";
			name += namebuf;
			GetToken();
		}

		if (currtoken == OBR) {
			/* function handling 
			 * in futuro ci potranno essere magari i dati strutturati */
			MathFunc_t* f = currNameSpace->GetFunc(namebuf);
			if (f == 0) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf.c_str(), "' not found");
			}
			GetToken();
			ExpressionElement *e = parsefunc(f);
			if (currtoken != CBR) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis "
						"expected after function "
						"\"", f->fname.c_str(), "\" in expr()");
			}
			GetToken();
			return e;
			
		} else {
			NamedValue* v = 0;
			if (currTable) {
				v = currTable->Get(namebuf);
			}

			if (v != NULL) {
				if (v->MayChange() || !ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
					return new EE_Var(v, currNameSpace);
				}
				return new EE_Value(v->GetVal());
			}
		}

		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unknown name \"", name.c_str(), "\"");
	}

	/* invalid expr */
	if (currtoken != ENDOFFILE) {
		throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unknown token");
	}

	// FIXME: is this ok?
	TypedValue a(0., true); 
	return new EE_Value(a);
}


ExpressionElement *
MathParser::stmt(void)
{
	if (currtoken == NAME) {
		bool bIsIfndef = false;
		bool bIsConst = false;

		DeclarationModifier declarationmodifier = GetDeclarationModifier(namebuf.c_str());
		if (declarationmodifier != DMOD_UNKNOWN) {
			switch (declarationmodifier) {
			case DMOD_IFNDEF:
				bIsIfndef = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled definition modifier ", namebuf.c_str(), "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type (modifier) expected "
					"after definition modifier in declaration");
			}
		}

		TypedValue::TypeModifier typemodifier = GetTypeModifier(namebuf.c_str());
		if (typemodifier != TypedValue::MOD_UNKNOWN) {
			switch (typemodifier) {
			case TypedValue::MOD_CONST:
				bIsConst = true;
				break;

			default:
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unhandled type modifier ", namebuf.c_str(), "");
			}

			if (GetToken() != NAME) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "type expected "
					"after type modifier in declaration");
			}
		}

		/* declaration? */
		TypedValue::Type type = GetType(namebuf.c_str());
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

			if (IsKeyWord(defaultNameSpace, namebuf.c_str())) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name '", namebuf.c_str(), "' "
						"is a keyword");
			}

			MathParser::NameSpace *currNameSpace = defaultNameSpace;
			Table *currTable = &table;
			std::string name(namebuf);

			GetToken();
			if (currtoken == NAMESPACESEP) {
				NameSpaceMap::iterator i = nameSpaceMap.find(name);
				if (i == nameSpaceMap.end()) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
				}
				currNameSpace = i->second;
				currTable = currNameSpace->GetTable();
				GetToken();
				if (currtoken != NAME) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
				}
				name += "::";
				name += namebuf;
				GetToken();
			}

			/* TODO: this needs a *LOT* of work
			 *
			 * - if the definition is without assignment?
			 * - if the definition is with assignment and it is const?
			 * - when can we safely evaluate the value?

				(const integer i = 10)	-> declaration, assignment, referentiation ::= EE_Declare, EE_Assign?
				(integer i = 10)	-> declaration, assignment, referentiation ::= EE_Declare, EE_Assign?
				(integer i)		-> declaration ::= EE_Declare?
				(i = 10)		-> assigment, referentiation ::= EE_Assign?
				(i)			-> referentiation ::= EE_Var

			 * if we pre-declare, and only use EE_Assign, we implicitly
			 * turn stuff like "real A = B + C" into "real A" once,
			 * followed by multiple "A = B + C", which may be good
			 * but is not consistent with earlier behavior
			 *
			 * Alternatively, we need EE_Declare (and possibly EE_Declare_Assign),
			 * although its only purpose is to fail when repeatedly invoked.
			 *
			 * Another point: an expression like

				((integer i = 10) && i)

			 * is now perfectly valid, since "i" is created before it is referenced
			 * after "&&".
			 *
			 * However, in the new implementation as it is now (pre-declaration)
			 * it works fine, but it would fail if we introduce EE_Declare,
			 * because "i" would not exist during parsing, it would only exist
			 * during evaluation.
			 *
			 * Pre-declaration seems good, but we need to find a way
			 * to intercept repeated declarations.  For instance:
			 * - pre-declaration and assignment of const
			 * - pre-declaration, deferred assignment of non-const
			 * - evaluation flag for all? Fail if evaluated more than once.
			 */

			/* with assign? */
			if (currtoken == ASSIGN) {
				if (currTable == 0) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
				}

				/* faccio una copia del nome! */
				std::string varname(namebuf);

				GetToken();
				ExpressionElement *e = logical();

				/* cerco la variabile */
				NamedValue* v = currTable->Get(varname.c_str());

				if (v == 0) {
					/* create new var with assigned type */
					TypedValue newvar(type);
					/* assign new var, so that it internally
					 * takes care of casting, while it inherits
					 * const'ness from d */
					TypedValue d(e->Eval());
					delete e;
					e = 0;

					if (bIsConst) {
						d.SetConst();
					}
					newvar.Cast(d);

					v = currTable->Put(varname.c_str(), newvar);

					if (bIsIfndef) {
						silent_cerr("warning, ifndef variable " << v->GetTypeName() << " \"" << name
							<< "\" not yet defined; set to \"" << newvar << "\" at line " << mbdyn_get_line_data() << std::endl);
					}

				} else {
					/* altrimenti, se la posso ridefinire, mi limito
					 * ad assegnarle il nuovo valore */
					if (!bRedefineVars && !bIsIfndef) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"var \"", name.c_str(), "\"");
					}

					if (v->Const() && !bIsIfndef) {
						// TODO: check redefinition of const'ness
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"a const named value "
							"\"", name.c_str(), "\"");
					}

					if (!v->IsVar()) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot redefine "
							"non-var named value "
							"\"", name.c_str(), "\"");
					}

					if (!bIsIfndef) {
						dynamic_cast<Var *>(v)->SetVal(e->Eval());

					} else {
						if (v->GetType() != type) {
							silent_cerr("warning, skipping redefinition of \"" << name.c_str() << "\""
								<< " from " << v->GetTypeName() << " to " << TypedValue::GetTypeName(type)
								<< " (orig=" << v->GetVal() << ", unchanged; new=" << e->Eval() << ")"
								<< " at line " << mbdyn_get_line_data() << std::endl);

						} else {
							silent_cerr("warning, skipping redefinition of " << v->GetTypeName() << " \"" << name.c_str() << "\""
								<< " (orig=" << v->GetVal() << ", unchanged; new=" << e->Eval() << ")"
								<< " at line " << mbdyn_get_line_data() << std::endl);
						}
					}

					delete e;
					e = 0;
				}

				if (bIsConst && ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
					return new EE_Value(dynamic_cast<Var *>(v)->GetVal());

				} else {
					return new EE_Var(v, currNameSpace);
				}

			} else if (currtoken == STMTSEP) {
				if (currTable == 0) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
				}
				
				NamedValue* v = currTable->Get(namebuf);
				if (v == NULL || (!bRedefineVars && !bIsIfndef)) {
					if (bIsConst) {
						/* cannot insert a const var
						 * with no value */
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot create const named value "
							"\"", namebuf.c_str(), "\" with no value");
					}

					/* se la var non esiste, la inserisco;
					 * se invece esiste e non vale 
					 * la ridefinizione, tento
					 * di inserirla comunque, cosi'
					 * table da' errore */
					v = currTable->Put(namebuf, TypedValue(type));
				}

				return new EE_Var(v, currNameSpace);
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

			MathParser::NameSpace *currNameSpace = defaultNameSpace;
			Table *currTable = &table;
			std::string name(namebuf);

			GetToken();
			if (currtoken == NAMESPACESEP) {
				NameSpaceMap::iterator i = nameSpaceMap.find(name);
				if (i == nameSpaceMap.end()) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "unable to find namespace \"", namebuf.c_str(), "\"");
				}
				currNameSpace = i->second;
				currTable = currNameSpace->GetTable();
				GetToken();
				if (currtoken != NAME) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "name expected after namespace");
				}
				name += "::";
				name += namebuf;
				GetToken();
			}

			if (currtoken == OBR) {
				/* in futuro ci potranno essere magari i dati strutturati */
				MathFunc_t* f = currNameSpace->GetFunc(namebuf);
				if (f == 0) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "function '", namebuf.c_str(), "' not found");
				}
				GetToken();
				ExpressionElement *e = parsefunc(f);
				if (currtoken != CBR) {
					throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "closing parenthesis "
							"expected after function "
							"\"", f->fname.c_str(), "\" in expr()");
				}
				GetToken();

				return logical(e);
			}

			/* assignment? */
			if (currTable == 0) {
				throw ErrGeneric(this, MBDYN_EXCEPT_ARGS, "namespace \"", currNameSpace->sGetName().c_str(), "\" does not support variables");
			}
			NamedValue* v = currTable->Get(namebuf);
			if (v != NULL) {
				if (currtoken == ASSIGN) {
					GetToken();
					ExpressionElement *e = logical();
					if (v->Const()) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
							"cannot assign const named value "
							"\"", name.c_str(), "\"");
					}

					if (!v->IsVar()) {
						throw MathParser::ErrGeneric(this, MBDYN_EXCEPT_ARGS,
								"cannot assign non-var named value "
								"\"", name.c_str(), "\"");
					}

					// FIXME: we need to define a EE_Assign that does the operation below
					dynamic_cast<Var *>(v)->Cast(e->Eval());
					delete e;
					return new EE_Var(v, currNameSpace);

				} else {
					// NOTE: fails if <name> is actually <namespace>::<name>
					// ASSERT(currtoken != NAME);
					// TokenPush(currtoken);
					// currtoken = NAME;

					// NOTE: fails if <name> is not the complete <stmt>
					// return v->GetVal();
					ExpressionElement *e;
					if (v->MayChange() || !ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
						e = new EE_Var(v, currNameSpace);
					} else {
						e = new EE_Value(v->GetVal());
					}
					return logical(e);
				}

			} else {
				ASSERT(0);
			}
		}
	}
	return logical();
}

ExpressionElement *
MathParser::stmtlist(void)
{
	ExpressionElement * d = stmt();
	if (currtoken == STMTSEP) {
		GetToken();
		return d = new EE_StmtList(d, stmtlist());
	}
	//return new EE_StmtList(NULL , d);
	return d;
}

ExpressionElement *
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
	struct PlugInRegister *p = 0;
	for (p = PlugIns; p != 0; p = p->next) {
		if (strcasecmp(p->name, pginname) == 0) {
			break;
		}
	}

	if (p == 0) {
		silent_cerr(" plugin '" << pginname << "' not supported" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

	return new EE_Var(v, defaultNameSpace);
}

#endif // USE_EE
