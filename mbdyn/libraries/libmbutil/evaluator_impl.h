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

#ifndef EVALUATOR_IMPL_H
#define EVALUATOR_IMPL_H

#include "evaluator.h"
#include <iostream>
#include <sstream>
#include <cmath>
#include <stdio.h>
#include <limits>

// for storing constant values like "1","1.23"
class EE_Value : public ExpressionElement {
private:
	TypedValue m_Val;

public:
	EE_Value(TypedValue Val) : m_Val(Val) { m_Val.SetConst(true); };
	// nothing to destroy
	~EE_Value(void) {};

	TypedValue
	Eval(void) const
	{
		return m_Val;
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << m_Val;
	};
};

// for storing non const variable like "a", "b"
// FIXME: probably needs to know about its namespace
class EE_Var : public ExpressionElement {
private:
	NamedValue *m_Var;
	// NOTE: probably the string would suffice, but just in case...
	MathParser::NameSpace *m_ns;

public:
	EE_Var(NamedValue *var, MathParser::NameSpace *ns = 0) : m_Var(var), m_ns(ns) {};
	// nothing to destroy (Var belongs to the symbol table)
	virtual ~EE_Var(void) {};

	virtual TypedValue
	Eval(void) const
	{
		// NOTE: this should not be const,
		// otherwise it would have been turned into a EE_Value
		return m_Var->GetVal();
	};

	virtual std::ostream&
	Output(std::ostream& out) const
	{
		if (m_ns) {
			out << m_ns->sGetName() << "::";
		}

		return out << m_Var->GetName();
	};
};

class EE_Plus : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Plus(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Plus(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		// FIXME: check integer overflow?
		return m_pEE1->Eval() + m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " + ", m_pEE2->Output(out) << ")";
	};
};

class EE_Minus : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Minus(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Minus(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		// FIXME: check integer overflow?
		return m_pEE1->Eval() - m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " - ", m_pEE2->Output(out) << ")";
	};
};

class EE_Modulus : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

#if 0 // TODO: check correctness whenever possible
	bool
	Check_int(const TypedValue& value)
	{
		switch (value.GetType()) {
		case TypedValue::VAR_BOOL:
		case TypedValue::VAR_INT:
			if (value.GetInt() == 0) {
				return false;
			}
			break;

		default:
			return false;
		}
	};
#endif

public:
	EE_Modulus(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Modulus(void) { delete m_pEE1; delete m_pEE2; };

#if 0 // TODO: check correctness whenever possible
	bool
	Check(void) const
	{
		if (dynamic_cast<EE_Value *>(m_pEE2)) {
			return Check_int(m_pEE2->Eval());
		}
		return true;
	};
#endif

	TypedValue
	Eval(void) const
	{
		// FIXME: the check on types should be done while parsing, if possible...
		TypedValue a = m_pEE2->Eval();
		if (a.GetType() == TypedValue::VAR_INT || a.GetType() == TypedValue::VAR_BOOL) {
			ASSERT(a.GetInt() != 0);
			return m_pEE1->Eval() % a;
		}
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " % ", m_pEE2->Output(out) << ")";
	};
};

class EE_Multiply : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Multiply(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Multiply(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() * m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " * ", m_pEE2->Output(out) << ")";
	};
};

class EE_Divide : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Divide(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Divide(void) { delete m_pEE1; delete m_pEE2; };
	// need to take care of zero denominator

	TypedValue
	Eval(void) const
	{
		TypedValue den = m_pEE2->Eval();
		if (den.GetType() == TypedValue::VAR_BOOL) {
			ASSERT(den.GetBool() != 0);

		} else if (den.GetType() == TypedValue::VAR_INT) {
			ASSERT(den.GetInt() != 0);

		} else if (den.GetType() == TypedValue::VAR_REAL) {
			Real value = den.GetReal();
			if (std::abs(value) < std::numeric_limits<double>::epsilon()) {
				// need to throw exception over here
				std::cout << "denominator cannot be zero" << std::endl;
			}
		}

		return m_pEE1->Eval() / den;
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " / ", m_pEE2->Output(out) << ")";
	};
};

class EE_Unary_minus : public ExpressionElement {
private:
	ExpressionElement *m_pEE1;

public:
	EE_Unary_minus(ExpressionElement *pEE1) : m_pEE1(pEE1) {};
	~EE_Unary_minus(void) { delete m_pEE1; };

	TypedValue
	Eval(void) const
	{
		return -m_pEE1->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(-", m_pEE1->Output(out) << ")";
	};
};

class EE_AND : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_AND(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_AND(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() && m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " && ", m_pEE2->Output(out) << ")";
	};
};

class EE_OR : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_OR(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_OR(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() || m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " || ", m_pEE2->Output(out) << ")";
	};
};

class EE_NOT : public ExpressionElement {
private:
	ExpressionElement *m_pEE1;

public:
	EE_NOT(ExpressionElement *pEE1) : m_pEE1(pEE1) {};
	~EE_NOT(void) { delete m_pEE1; };

	TypedValue
	Eval(void) const
	{
		return  !(m_pEE1->Eval());
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(!", m_pEE1->Output(out) << ")";
	};
};

// need to implement ~= as overloaded operator in mathp.cc
class EE_XOR : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_XOR(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_XOR(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		TypedValue v1 = m_pEE1->Eval();
		TypedValue v2 = m_pEE2->Eval();
		return ((!(v1 && v2)) && (v1 || v2));
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " ~| ", m_pEE2->Output(out) << ")";
	};
};

class EE_Greater : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Greater(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Greater(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() > m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " > ", m_pEE2->Output(out) << ")";
	};
};

class EE_Greater_Equal : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Greater_Equal(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Greater_Equal(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() >= m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " >= ", m_pEE2->Output(out) << ")";
	};
};

class EE_Lesser : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Lesser(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Lesser(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() < m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " < ", m_pEE2->Output(out) << ")";
	};
};

class EE_Lesser_Equal : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Lesser_Equal(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Lesser_Equal(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return m_pEE1->Eval() <= m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " <= ", m_pEE2->Output(out) << ")";
	};
};

class EE_Equal_Equal : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Equal_Equal(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Equal_Equal(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return (m_pEE1->Eval()) == (m_pEE2->Eval());
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " == ", m_pEE2->Output(out) << ")";
	};
};

class EE_Not_Equal : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Not_Equal(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Not_Equal(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		return (m_pEE1->Eval()) != (m_pEE2->Eval());
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " != ", m_pEE2->Output(out) << ")";
	};
};

class EE_Power : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_Power(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_Power(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		TypedValue a = m_pEE1->Eval();
		TypedValue b = m_pEE2->Eval();
		if (a.GetType() == TypedValue::VAR_STRING || b.GetType() == TypedValue::VAR_STRING) {
			silent_cerr(" power undefined between string " << a << " and " << b << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if ((a.GetType() == TypedValue::VAR_INT && b.GetType() == TypedValue::VAR_INT)
			|| (a.GetType() == TypedValue::VAR_BOOL && b.GetType() == TypedValue::VAR_BOOL))
		{
			Real value = std::pow(a.GetInt(), b.GetInt());
			if (Real(Int(value)) != value) {
				silent_cerr(" power(" << a.GetInt() << ", " << b.GetInt() << ") overflows (use explicit/implicit cast?)" << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			return value;

		} else {
			return std::pow(a.GetReal(), b.GetReal());
		}
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " ^ ", m_pEE2->Output(out) << ")";
	};
};

#if 1
// FIXME: needs review
// for storing value that is which is evaluated first then assigned

class EE_Assign : public ExpressionElement {
private:
	mutable Var *m_Var;
	MathParser::NameSpace *m_ns;
	const ExpressionElement *m_pEE;

public:
	EE_Assign(Var *var, MathParser::NameSpace *ns, ExpressionElement *pEE) : m_Var(var), m_ns(ns), m_pEE(pEE) {};
	~EE_Assign(void) { delete m_pEE; };

	TypedValue
	Eval(void) const
	{
		TypedValue v(m_pEE->Eval());
		m_Var->SetVal(v);
		return v;
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		out << "(";
		if (m_ns) {
			out << m_ns->sGetName() << "::";
		}

		return out << m_Var->GetName() << " = ", m_pEE->Output(out) << ")";
	};
};

/* the variable is pre-declared */
class EE_DeclareAssign : public ExpressionElement {
private:
	mutable Var *m_Var;
	MathParser::NameSpace *m_ns;
	const ExpressionElement *m_pEE;
	mutable bool m_bEvaluated;

public:
	EE_DeclareAssign(Var *var, MathParser::NameSpace *ns, ExpressionElement *pEE = 0) : m_Var(var), m_ns(ns), m_pEE(pEE), m_bEvaluated(false) {};
	~EE_DeclareAssign(void) { if (m_pEE) { delete m_pEE; } };

	TypedValue
	Eval(void) const
	{
		if (m_bEvaluated) {
			silent_cerr(" variable " << m_Var->GetName() << " already declared" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		m_bEvaluated = true;

		if (m_pEE) {
			TypedValue v(m_pEE->Eval());
			m_Var->SetVal(v);
			return v;
		}

		return TypedValue(0, true);
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		out << "(";
		if (m_Var->Const()) {
			out << "const ";
		}

		out << m_Var->GetTypeName() << " ";
		
		if (m_ns) {
			out << m_ns->sGetName() << "::";
		}

		return out << m_Var->GetName() << " = ", m_pEE->Output(out) << ")";
	};
};
#endif

// for storing like <stmt>;<stmt>
class EE_StmtList : public ExpressionElement {
private:
	ExpressionElement *m_pEE1, *m_pEE2;

public:
	EE_StmtList(ExpressionElement *pEE1, ExpressionElement *pEE2) : m_pEE1(pEE1), m_pEE2(pEE2) {};
	~EE_StmtList(void) { delete m_pEE1; delete m_pEE2; };

	TypedValue
	Eval(void) const
	{
		// FIXME: the first expression must exist
		ASSERT(m_pEE1 != 0);

		if (m_pEE1 != 0) {
			m_pEE1->Eval();
		}

		return m_pEE2->Eval();
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		return out << "(", m_pEE1->Output(out) << " ; ", m_pEE2->Output(out) << ")";
	};
};

class EE_Func : public ExpressionElement {
private:
	// pointer to MathParser only needed for exception handling
	MathParser *m_p;
	MathParser::MathFunc_t *m_f;

public:
	EE_Func(MathParser *p, MathParser::MathFunc_t *f)
	: m_p(p), m_f(f)
	{
	};

	~EE_Func(void)
	{
		for (unsigned i = 0; i < m_f->args.size(); ++i) {
			const ExpressionElement *ee = m_f->args[i]->GetExpr();
			if (ee) {
				delete ee;
			}
		}
		delete m_f;
	};

	TypedValue
	Eval(void) const
	{
		// sure args[0] does not need to be evaluated: it's the return value
		for (unsigned i = 1; i < m_f->args.size(); i++) {
			m_f->args[i]->Eval();
		}

		if (m_f->t != 0) {
			if (m_f->t(m_f->args)) {
				DEBUGCERR("error in function "
					<< m_f->ns->sGetName() << "::" << m_f->fname
					<< " " "(msg: " << m_f->errmsg << ")"
					<< " in EE_Func()" << std::endl);
				throw MathParser::ErrGeneric(m_p, MBDYN_EXCEPT_ARGS, m_f->fname.c_str(), ": error ", m_f->errmsg.c_str());
			}
		}

		TypedValue val = m_f->ns->EvalFunc(m_f);

		return val;
	};

	std::ostream&
	Output(std::ostream& out) const
	{
		if (m_f->ns) {
			out << m_f->ns->sGetName() << "::";
		}

		out << m_f->fname << "(";
		const ExpressionElement *ee = m_f->args[1]->GetExpr();
		if (ee != 0) {
			ee->Output(out);
			for (int i = 2; i < m_f->args.size(); ++i) {
				ee = m_f->args[i]->GetExpr();
				if (ee == 0) {
					break;
				}
				out << ", ", ee->Output(out);
			}
		}
		out << ")";

		return out;
	};
};

// unary operator construction helper with optimization
template <class T>
ExpressionElement *
EECreate(ExpressionElement *e1)
{
	ExpressionElement *out = new T(e1);
	if (ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
		if (dynamic_cast<EE_Value *>(e1)) {
			ExpressionElement *tmp = new EE_Value(out->Eval());
			delete out;
			return tmp;
		}
	}
	return out;
}

// binary operator construction helper with optimization
template <class T>
ExpressionElement *
EECreate(ExpressionElement *e1, ExpressionElement *e2)
{
	ExpressionElement *out = new T(e1, e2);
	if (ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
		if (dynamic_cast<EE_Value *>(e1) && dynamic_cast<EE_Value *>(e2)) {
			ExpressionElement *tmp = new EE_Value(out->Eval());
			delete out;
			return tmp;
		}
	}
	return out;
}

// ternary operator construction helper with optimization
template <class T>
ExpressionElement *
EECreate(ExpressionElement *e1, ExpressionElement *e2, ExpressionElement *e3)
{
	ExpressionElement *out = new T(e1, e2, e3);
	if (ExpressionElement::IsFlag(ExpressionElement::EE_CONSTIFY)) {
		if (dynamic_cast<EE_Value *>(e1) && dynamic_cast<EE_Value *>(e2) && dynamic_cast<EE_Value *>(e3)) {
			ExpressionElement *tmp = new EE_Value(out->Eval());
			delete out;
			return tmp;
		}
	}
	return out;
}

#endif // EVALUATOR_IMPL_H
