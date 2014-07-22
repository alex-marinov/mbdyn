/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef MATHTYP_H
#define MATHTYP_H

#include <string>
#include "except.h"

/* definizione dei tipi */
typedef double Real;
typedef int Int;

/* valori con tipo */
class TypedValue {
public:
	enum Type {
		VAR_UNKNOWN = -1,

		VAR_BOOL,
		VAR_INT,
		VAR_REAL,
		VAR_STRING,

		VAR_LAST
	};

	enum TypeModifier {
		MOD_UNKNOWN = -1,

		MOD_CONST,

		MOD_LAST
	};

	class ErrUnknownType : public MBDynErrBase {
	public:
		ErrUnknownType(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrWrongType : public MBDynErrBase {
	public:
		ErrWrongType(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
		ErrWrongType(const char *file, int line, const char *func,
			const TypedValue::Type& to,
			const TypedValue::Type& from);
	};
	class ErrUnknownValue : public MBDynErrBase {
	public:
		ErrUnknownValue(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrConstraintViolation : public MBDynErrBase {
	public:
		ErrConstraintViolation(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

protected:
	TypedValue::Type type;
	bool bConst;
	union {
		Int i;
		Real r;
	} v;
	std::string s;

public:
	TypedValue(void);
	~TypedValue(void);
	TypedValue(const bool& b, bool isConst = false);
	TypedValue(const Int& i, bool isConst = false);
	TypedValue(const Real& r, bool isConst = false);
	TypedValue(const std::string& s, bool isConst = false);
	TypedValue(const TypedValue::Type t, bool isConst = false);
	TypedValue(const TypedValue& var);

	const TypedValue& operator = (const TypedValue& var);
	const TypedValue& Cast(const TypedValue& var);

	TypedValue::Type GetType(void) const;
	const char *const GetTypeName(void) const;
	static const char *const GetTypeName(TypedValue::Type t);
	bool Const(void) const;
	bool GetBool(void) const;
	Int GetInt(void) const;
	Real GetReal(void) const;
	const std::string& GetString(void) const;

	void SetType(TypedValue::Type t, bool isConst = false);
	void SetConst(bool isConst = true, bool bForce = false);
	const TypedValue& Set(const bool& b);
	const TypedValue& Set(const Int& i);
	const TypedValue& Set(const Real& r);
	const TypedValue& Set(const std::string& s);

	bool operator && (const TypedValue& v) const;
	bool operator || (const TypedValue& v) const;
	bool operator > (const TypedValue& v) const;
	bool operator >= (const TypedValue& v) const;
	bool operator == (const TypedValue& v) const;
	bool operator <= (const TypedValue& v) const;
	bool operator < (const TypedValue& v) const;
	bool operator != (const TypedValue& v) const;

	TypedValue operator + (const TypedValue& v) const;
	TypedValue operator - (const TypedValue& v) const;
	TypedValue operator * (const TypedValue& v) const;
	TypedValue operator / (const TypedValue& v) const;
	TypedValue operator % (const TypedValue& v) const;

	const TypedValue& operator += (const TypedValue& v);
	const TypedValue& operator -= (const TypedValue& v);
	const TypedValue& operator *= (const TypedValue& v);
	const TypedValue& operator /= (const TypedValue& v);
	const TypedValue& operator %= (const TypedValue& v);
};

extern bool operator ! (const TypedValue& v);
extern TypedValue operator - (const TypedValue& v);
extern TypedValue operator + (const TypedValue& v);
extern std::ostream& operator << (std::ostream& out, const TypedValue& v);


/* classe per la memorizzazione delle variabili */
class NamedValue {
private:
	char *name;

	void AllocName(const char *const s);

public:
	NamedValue(const char *const s);
	virtual ~NamedValue(void);

	virtual bool IsVar(void) const;

	const char *GetName(void) const;
	virtual TypedValue::Type GetType(void) const = 0;
	virtual const char *const GetTypeName(void) const;
	virtual bool Const(void) const = 0;
	virtual TypedValue GetVal(void) const = 0;
};

class Var : public NamedValue {
private:
	TypedValue value;

public:
	Var(const char* const s, const TypedValue& v);
	Var(const char* const s, const bool& b);
	Var(const char* const s, const Int& v);
	Var(const char* const s, const Real& v);
	Var(const char* const s, const std::string& v);
	~Var(void);

	bool IsVar(void) const;

	TypedValue::Type GetType(void) const;
	bool Const(void) const;
	TypedValue GetVal(void) const;

	void SetVal(const bool& b);
	void SetVal(const Int& v);
	void SetVal(const Real& v);
	void SetVal(const std::string& v);
	void SetVal(const TypedValue& v);
	void Cast(const TypedValue& v);
};

#endif /* MATHTYP_H */

