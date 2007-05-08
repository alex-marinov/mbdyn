/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
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

	*out = F((*arg1)(), (*arg2)());

	return 0;
}

static int
mp_asin_t(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);
	ASSERT(args[1]->Type() == MathParser::AT_REAL);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t *>(args[1]);
	ASSERT(arg1 != 0);

	if ((*arg1)() > 1. || (*arg1)() < -1.) {
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

	if ((*arg1)() < 0.) {
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

	Real a = (*arg1)() - int((*arg1)()/M_PI)*M_PI;
	if (fabs(fabs(a) - M_PI_2) < DBL_EPSILON) {
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

	if ((*arg1)() >= 1. || (*arg1)() <= -1.) {
		return 1;
	}

	return 0;
}

static int
mp_log_t(const MathParser::MathArgs& args)
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

static std::ostream&
operator << (std::ostream& out, const MathParser::MathArgVoid_t& v)
{
	return out;
}

static std::ostream&
operator << (std::ostream& out, const MathParser::MathArgInt_t& v)
{
	return out << v();
}

static std::ostream&
operator << (std::ostream& out, const MathParser::MathArgReal_t& v)
{
	return out << v();
}

static std::ostream&
operator << (std::ostream& out, const MathParser::MathArgString_t& v)
{
	return out << v();
}

static int
mp_print(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 1);

	switch (args[1]->Type()) {
	case MathParser::AT_INT:
		silent_cout((*dynamic_cast<MathParser::MathArgInt_t*>(args[1]))() << std::endl);
		break;

	case MathParser::AT_REAL:
		silent_cout((*dynamic_cast<MathParser::MathArgReal_t*>(args[1]))() << std::endl);
		break;

	case MathParser::AT_STRING:
		silent_cout((*dynamic_cast<MathParser::MathArgString_t*>(args[1]))() << std::endl);
		break;

	default:
		throw ErrGeneric();
	}

	return 0;
}

static int
mp_stop(const MathParser::MathArgs& args)
{
	ASSERT(args.size() == 1 + 2);
	ASSERT(args[0]->Type() == MathParser::AT_VOID);
	ASSERT(args[1]->Type() == MathParser::AT_INT);
	ASSERT(args[1]->Type() == MathParser::AT_INT);

	MathParser::MathArgInt_t *s = dynamic_cast<MathParser::MathArgInt_t *>(args[1]);
	ASSERT(s != 0);

	MathParser::MathArgInt_t *v = dynamic_cast<MathParser::MathArgInt_t *>(args[2]);
	ASSERT(v != 0);

	if ((*s)() != 0) {
		if ((*v)() == 0) {
			silent_cout("mp_stop(SUCCESS)" << std::endl);
			throw NoErr();

		} else {
			silent_cout("mp_stop(FAILURE)" << std::endl);
			throw ErrGeneric();
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

	Real a = (*arg1)() - int((*arg1)()/M_PI)*M_PI;
	if (fabs(a) < DBL_EPSILON) {
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

	*out = atan2((*arg2)(), (*arg1)());

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

	if (fabs((*arg1)()) < DBL_EPSILON) {
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

	MathParser::MathArgReal_t *out = dynamic_cast<MathParser::MathArgReal_t*>(args[0]);
	ASSERT(out != 0);

	MathParser::MathArgReal_t *arg1 = dynamic_cast<MathParser::MathArgReal_t*>(args[1]);
	ASSERT(arg1 != 0);

	*out = copysign(1., (*arg1)());

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

	*out = std::max((*arg1)(), (*arg2)());

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

	*out = std::min((*arg1)(), (*arg2)());

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

	if ((*arg1)() >= 0.) {
		*out = 1.;

	} else if ((*arg1)() == 0.) {
		*out = .5;

	} else {
		*out = 0.;
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

	if ((*arg1)() > 0.) {
		*out = (*arg1)();

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

	if ((*arg1)() < 0.) {
		*out = 0.;

	} else if ((*arg1)() > (*arg2)()) {
		*out = (*arg2)();

	} else {
		*out = (*arg1)();
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

	if ((*arg1)() > 0.) {
		*out = (*arg1)()*(*arg1)();

	} else {
		*out = 0.;
	}

	return 0;
}

/* tipi delle variabili */
struct TypeName_t {
	const char* name;
	TypedValue::Type type;
};

static TypeName_t TypeNames[] = {
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


TypedValue::TypedValue(void)
: type(TypedValue::VAR_UNKNOWN), bConst(false)
{
	NO_OP;
}

TypedValue::~TypedValue(void)
{
	NO_OP;
}

TypedValue::TypedValue(const Int& i, bool isConst)
: type(TypedValue::VAR_INT), bConst(isConst)
{
	v.i = i;
}

TypedValue::TypedValue(const bool& b, bool isConst)
: type(TypedValue::VAR_INT), bConst(isConst)
{
	v.i = b;
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

	case TypedValue::VAR_STRING:
		s = var.s;
		break;

	case TypedValue::VAR_UNKNOWN:
		/* allow initialization from unspecified type right now */
		break;

	default:
		throw ErrUnknownType();
	}
}

const TypedValue&
TypedValue::operator = (const TypedValue& var)
{
	if (Const()) {
		throw ErrConstraintViolation();
	}

	if (type == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (var.type) {
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)var.GetInt());
			Set(buf);
			break;

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)var.GetReal());
			Set(buf);
			break;

		case TypedValue::VAR_STRING:
			this->s = var.GetString();
			break;

		default:
			throw ErrUnknownType();
		}

	} else {
		switch (var.type) {
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
				throw ErrWrongType();
			}
			type = TypedValue::VAR_STRING;
			s = var.s;
			break;

		case TypedValue::VAR_UNKNOWN:
			if (type != TypedValue::VAR_UNKNOWN) {
				throw ErrWrongType();
			}
			break;

		default:
			throw ErrUnknownType();
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

	case TypedValue::VAR_STRING:
		throw ErrWrongType();

	default:
		throw ErrUnknownType();
	}
}

Real
TypedValue::GetReal(void) const
{
	switch (type) {
	case TypedValue::VAR_INT:
		return Real(v.i);

	case TypedValue::VAR_REAL:
 		return v.r;

	case TypedValue::VAR_STRING:
		throw ErrWrongType();

	default:
		throw ErrUnknownType();
	}
}

const std::string &
TypedValue::GetString(void) const
{
	switch (type) {
	case TypedValue::VAR_INT:
	case TypedValue::VAR_REAL:
		throw ErrWrongType();

	case TypedValue::VAR_STRING:
 		return s;

	default:
		throw ErrUnknownType();
	}
}

void
TypedValue::SetType(TypedValue::Type t, bool isConst)
{
	if (Const()) {
		throw ErrConstraintViolation();
	}

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
	if (Const()) {
		throw ErrConstraintViolation();
	}

	switch (GetType()) {
	case TypedValue::VAR_INT:
		v.i = i;
		break;

	case TypedValue::VAR_REAL:
		v.r = Real(i);
		break;

	case TypedValue::VAR_STRING:
		throw ErrWrongType();

	default:
		throw ErrUnknownType();
	}

	return *this;
}

const TypedValue&
TypedValue::Set(const Real& r)
{
	if (Const()) {
		throw ErrConstraintViolation();
	}

	switch (GetType()) {
	case TypedValue::VAR_INT:
		v.i = Int(r);
		break;

	case TypedValue::VAR_REAL:
		v.r = r;
		break;

	case TypedValue::VAR_STRING:
		throw ErrWrongType();

	default:
		throw ErrUnknownType();
	}

	return *this;
}

const TypedValue&
TypedValue::Set(const std::string& s)
{
	if (Const()) {
		throw ErrConstraintViolation();
	}

	switch (GetType()) {
	case TypedValue::VAR_INT:
	case TypedValue::VAR_REAL:
		throw ErrWrongType();

	case TypedValue::VAR_STRING:
		this->s = s;
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
	if (GetType() == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (v.GetType()) {
		case TypedValue::VAR_INT:
			snprintf(buf, sizeof(buf), "%ld", (long)v.GetInt());
			return TypedValue(GetString() + buf);

		case TypedValue::VAR_REAL:
			snprintf(buf, sizeof(buf), "%e", (double)v.GetReal());
			return TypedValue(GetString() + buf);

		case TypedValue::VAR_STRING:
			return TypedValue(GetString() + v.GetString());

		default:
			throw ErrWrongType();
		}
	}

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
	if (Const()) {
		throw ErrConstraintViolation();
	}

	if (GetType() == TypedValue::VAR_STRING) {
		char buf[BUFSIZ];

		switch (v.GetType()) {
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
			throw ErrWrongType();
		}

		return *this;
	}

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
	if (Const()) {
		throw ErrConstraintViolation();
	}

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
	if (Const()) {
		throw ErrConstraintViolation();
	}

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
	if (Const()) {
		throw ErrConstraintViolation();
	}

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
	if (v.GetType() == TypedValue::VAR_STRING) {
		return v.GetString().empty();
	}

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

	case TypedValue::VAR_STRING:
		throw TypedValue::ErrWrongType();

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

	case TypedValue::VAR_STRING:
		return out << v.GetString();

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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
	}

	// log
	f = new MathFunc_t;
	f->fname = std::string("log");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, log>;
	f->t = mp_log_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// log10
	f = new MathFunc_t;
	f->fname = std::string("log10");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, log10>;
	f->t = mp_log_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
	}

	// sqrt
	f = new MathFunc_t;
	f->fname = std::string("sqrt");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_func_1<MathParser::MathArgReal_t, MathParser::MathArgReal_t, sqrt>;
	f->t = mp_log_t;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
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
		throw ErrGeneric();
	}

	// sign
	f = new MathFunc_t;
	f->fname = std::string("sign");
	f->args.resize(1 + 1);
	f->args[0] = new MathArgReal_t;
	f->args[1] = new MathArgReal_t;
	f->f = mp_sign;
	f->t = 0;

	if (!func.insert(funcType::value_type(f->fname, f)).second) {
		silent_cerr("static namespace: "
			"unable to insert handler "
			"for function " << f->fname << std::endl);
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
	}

	// sramp
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
		throw ErrGeneric();
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
		throw ErrGeneric();
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
		throw ErrGeneric();
	}
}

MathParser::StaticNameSpace::~StaticNameSpace(void)
{
	for (funcType::iterator f = func.begin(); f != func.end(); f++) {
		for (MathParser::MathArgs::iterator i = f->second->args.begin();
			i != f->second->args.end();
			i++)
		{
			delete *i;
		}

		delete f->second;
	}
}

bool
MathParser::StaticNameSpace::IsFunc(const char* const s) const
{
	if (func.find(std::string(s)) != func.end()) {
		return true;
	}

	return false;
}

MathParser::MathFunc_t*
MathParser::StaticNameSpace::GetFunc(const char* const s) const
{
	funcType::const_iterator i = func.find(std::string(s));

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

	case MathParser::AT_INT:
		return TypedValue((*dynamic_cast<MathArgInt_t*>(args[0]))());

	case MathParser::AT_REAL:
		return TypedValue((*dynamic_cast<MathArgReal_t*>(args[0]))());

	case MathParser::AT_STRING:
		return TypedValue((*dynamic_cast<MathArgString_t*>(args[0]))());

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
		while ((c = in->get()), isalnum(c)
			|| c == '_'
			|| ((currtoken == NAMESPACESEP) && c == ':'))
		{
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
		value.Set(namebuf);
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
MathParser::evalfunc(MathParser::NameSpace *ns, MathParser::MathFunc_t* f)
{
	MathArgs args = f->args;

	for (unsigned i = 1; i < args.size(); i++) {
		switch (args[i]->Type()) {
		case MathParser::AT_INT:
			(*dynamic_cast<MathArgInt_t*>(args[i]))() = stmtlist().GetInt();
			break;

		case MathParser::AT_REAL:
			(*dynamic_cast<MathArgReal_t*>(args[i]))() = stmtlist().GetReal();
			break;

		case MathParser::AT_PRIVATE:
			/* ignore */
			break;

		default:
			throw ErrGeneric();
		}

		if (i < args.size() - 1) {
			if (args[i + 1]->Type() != AT_PRIVATE) {
				if (currtoken != ARGSEP) {
					throw ErrGeneric(this,
						"argument separator expected");
				}
				GetToken();
			}
		}
	}
			
	if (f->t != 0) {
		if (f->t(args)) {
			DEBUGCERR("error in function "
				<< ns->sGetName() << "::" << f->fname 
				<< " " "(msg: " << f->errmsg << ")"
				<< " in evalfunc()" << std::endl);
			throw ErrGeneric(this, f->fname.c_str(), ": error ", f->errmsg.c_str());
		}
	}

	return ns->EvalFunc(f, args);
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
		MathParser::NameSpace *currNameSpace = defaultNameSpace;
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

			MathParser::MathFunc_t* f = currNameSpace->GetFunc(namebuf);
			if (f == NULL) {
				throw ErrGeneric(this, "function '", namebuf, "' not found");
			}
			GetToken();
			TypedValue d = evalfunc(currNameSpace, f);
			if (currtoken != CBR) {
				throw ErrGeneric(this, "closing parenthesis "
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
			if (IsKeyWord(defaultNameSpace, namebuf)) {
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
					newvar = d;
					v = table.Put(varname, newvar);

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

						if (!v->IsVar()) {
				   			throw MathParser::ErrGeneric(this,
								"cannot redefine "
								"non-var named value "
								"\"", v->GetName(), "\"");
						}
						dynamic_cast<Var *>(v)->SetVal(d);

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
					if (isConst) {
						/* cannot insert a const var
						 * with no value */
			      			throw MathParser::ErrGeneric(this,
		      						"cannot create const named value "
								"\"", v->GetName(), "\" with no value");
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

			 		if (!v->IsVar()) {
						throw MathParser::ErrGeneric(this,
								"cannot assign non-var named value "
								"\"", v->GetName(), "\"");
			 		}
					dynamic_cast<Var *>(v)->SetVal(d);
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
defaultNameSpace(0),
namebuf(NULL),
namebuflen(default_namebuflen),
value(Real(0)),
varlist(NULL),
tokenlist(NULL),
powerstack()
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	defaultNameSpace = new StaticNameSpace();
	if (RegisterNameSpace(defaultNameSpace)) {
		throw ErrGeneric(this, "unable to register namespace "
				"\"", defaultNameSpace->sGetName(), "\"");
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
defaultNameSpace(0),
namebuf(NULL),
namebuflen(default_namebuflen),
value(Real(0)),
varlist(NULL),
tokenlist(NULL),
powerstack()
{
	DEBUGCOUTFNAME("MathParser::MathParser");

	SAFENEWARR(namebuf, char, namebuflen + 1);

	defaultNameSpace = new StaticNameSpace();
	if (RegisterNameSpace(defaultNameSpace)) {
		throw ErrGeneric(this, "unable to register namespace "
				"\"", defaultNameSpace->sGetName(), "\"");
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

	while (PlugIns) {
		PlugInRegister *next = PlugIns->next;
		SAFEDELETEARR(PlugIns->name);
		SAFEDELETE(PlugIns);
		PlugIns = next;
	}

	for (NameSpaceMap::iterator i = nameSpaceMap.begin(); i != nameSpaceMap.end(); i++) {
		delete i->second;
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
		throw ErrGeneric(this, "statement separator expected");
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
		throw ErrGeneric(this, "separator expected");
	}
	in = (InputStream*)p;

	return vv;
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

	pedantic_cout("registering namespace \"" << ns->sGetName() << "\""
			<< std::endl);

	std::string name(ns->sGetName());

	if (nameSpaceMap.find(name) != nameSpaceMap.end()) {
		return 1;
	}

	nameSpaceMap[name] = ns;

	return 0;
}
