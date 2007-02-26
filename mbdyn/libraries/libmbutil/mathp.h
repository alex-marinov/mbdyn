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

#ifndef MATHP_H
#define MATHP_H

#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <ac/math.h>
#include <time.h>
#include <ac/iostream>
#include <map>
#include <string>

#include <myassert.h>
#include <mynewmem.h>
#include <except.h>

#include <mathtyp.h>
#include <table.h>
#include <input.h>

#include <stack.h>

typedef Real (*MathFunc_0args_t)(void);
typedef Real (*MathFunc_1args_t)(Real);
typedef Real (*MathFunc_2args_t)(Real, Real);
typedef Int (*MathFuncTest_t)(Real *);

/* struttura delle funzioni built-in */
struct MathFunc_t {
	const char* fname;		/* nome */
	Int nargs;			/* numero di argomenti */
	union {
		MathFunc_0args_t	f0;
		MathFunc_1args_t	f1;
		MathFunc_2args_t	f2;
		/* add more as needed */
	} f;				/* puntatore a funzione */
	MathFuncTest_t	t;		/* puntatore a funzione di test */
	const char* errmsg;		/* messaggio di errore */
};


/* massimo numero di argomenti */
const Int max_nargs = 2;		/* keep it updated */

const char ONE_LINE_REMARK = '#';		/* carattere di inizio commento */

class MathParser {
public:
	/* Namespace */
	class NameSpace {
		const char	*name;

	public:
		NameSpace(const char *n);
		virtual ~NameSpace(void);
		virtual const char *sGetName(void) const;
		virtual bool IsFunc(const char* const s) const = 0;
		virtual MathFunc_t* GetFunc(const char* const s) const = 0;
		virtual TypedValue EvalFunc(MathFunc_t *f, Real *d) const = 0;
	};

	/* Static namespace */
	class StaticNameSpace : public MathParser::NameSpace {
	private:
		MathFunc_t	*func;

	public:
		StaticNameSpace(const char *n, MathFunc_t f[]);
		~StaticNameSpace(void);

		bool IsFunc(const char* const s) const;
		MathFunc_t* GetFunc(const char* const s) const;
		virtual TypedValue EvalFunc(MathFunc_t *f, Real *d) const;
	};

	/* Prototipo dei plugin */
	class PlugIn {
	protected:
		MathParser& mp;
 
	/*
	 * Idea base: 
	 * un plugin e' un qualcosa che viene collegato al parser
	 * dopo la sua costruzione.  Il plugin e' dotato di 'tipo'
	 * (sName).  Quando il parser incontra una parentesi quadra
	 * che si apre ('['), si comporta come segue:
	 *   - esegue il parsing dello stream fino alla chiusura
	 *     della parentesi (']'), 
	 *   - separando i token in base alle virgole (','), 
	 *   - quindi legge il primo token che contiene il nome del
	 *     plugin
	 *   - legge il secondo token con il quale definisce una
	 *     metavariabile a cui il risultato della creazione
	 *     del plugin viene assegnato.
	 * Il risultato e' una metavariabile, a cui e' associata una
	 * certa operazione.  Ogni qualvolta la variabile viene usata,
	 * questo provoca l'esecuzione del plugin.
	 */
	public:
		PlugIn(MathParser& mp) : mp(mp) {};
		virtual ~PlugIn() {};
		virtual const char *sName(void) const = 0;
		virtual int Read(int argc, char *argv[]) = 0;
		virtual TypedValue::Type GetType(void) const = 0;
		virtual TypedValue GetVal(void) const = 0;
	};

protected:
	class PlugInVar : public NamedValue {
	private:
		MathParser::PlugIn *pgin;

	public:
		PlugInVar(const char *const s, MathParser::PlugIn *p);
		~PlugInVar(void);

		TypedValue::Type GetType(void) const;
		bool Const(void) const;
		TypedValue GetVal(void) const;
	};

	struct PlugInRegister {
		const char *name;
		MathParser::PlugIn * (*constructor)(MathParser&, void *);
		void *arg;
		struct PlugInRegister *next; 

		PlugInRegister(void)
		: name(NULL), constructor(NULL), arg(NULL) {};
	} *PlugIns;

public:

	/* Gestione degli errori */
	class ErrGeneric {
	public:
		ErrGeneric(void);
		ErrGeneric(MathParser* p, const char* const s);
		ErrGeneric(MathParser* p, const char* const s1, 
				const char* const s2, const char* const s3);
	};
   
	Table&   table;      /* symbol table */
	int redefine_vars;  /* redefine_vars flag */
   
	/* gioca con table e stream di ingresso */
	Table& GetSymbolTable(void) const;
	void PutSymbolTable(Table& T);

	InputStream* in;     /* stream in ingresso */
	int GetLineNumber(void) const;

public:
   
	/* i token che il lex riconosce */
	enum Token {
		ENDOFFILE = -2,
		UNKNOWNTOKEN = -1,
	
		MP_INT = TypedValue::VAR_INT,
		MP_REAL = TypedValue::VAR_REAL,	
	
		NUM,		/* Numero */
		NAME, 		/* Nome */
		EXP,		/* '^'	: Elevamento a potenza */
		MULT,		/* '*'	: Moltiplicazione */
		DIV,		/* '/'	: Divisione */
		MINUS,		/* '-'	: Meno */
		PLUS,		/* '+'	: Piu' */
		GT,		/* '>'	: Maggiore di */
		GE,		/* '>='	: Maggiore o uguale */
		EQ,		/* '=='	: Uguale */
		LE,		/* '<='	: Minore o uguale */
		LT,		/* '<'	: Minore di */
		NE,		/* '!='	: Diverso da */
		NOT,		/* '!'	: Negazione (operatore logico) */
		AND,		/* '&&'	: AND (operatore logico) */
		OR,		/* '||'	: OR (operatore logico) */
		XOR,		/* '~|'	: XOR, o OR esclusivo (op. logico) */
		
		OBR,		/* '('	: Parentesi aperta */
		CBR,		/* ')'	: Parentesi chiusa */
		OPGIN,		/* '['	: Apertura di plugin statement */
		CPGIN,		/* ']'	: Chiusura di plugin statement */
		STMTSEP,	/* ';'	: Separatore di statements */
		ARGSEP,		/* ','	: Separatore di argomenti */
		NAMESPACESEP,	/* '::'	: Separatore di namespace */
		ASSIGN,		/* '='	: Assegnazione */
		
		LASTTOKEN	
	};
   
protected:

	StaticNameSpace* defaultNameSpace;
	typedef std::map<std::string, NameSpace *> NameSpaceMap;
	NameSpaceMap nameSpaceMap;

	/* buffer statico reallocabile per leggere nomi */
	char* namebuf;
	Int namebuflen;

	/* valore numerico dotato di tipo */
	TypedValue value;
   
	/* token corrente */
	enum Token currtoken;
      
	/* inizio della lista */
	VarList* varlist;

	/* lista dei token; implementa una stack rudimentale */
	struct TokenList {
		enum Token t;
		TypedValue value;     
		char* name;

		TokenList* next;
      
		TokenList(Token t);
		TokenList(const char* const s);
		TokenList(const TypedValue& v);
		~TokenList(void);
	};
	
	/* cima della stack */
	TokenList* tokenlist;

	/* operazioni sulla stack */
	void TokenPush(enum Token t);
	int TokenPop(void);

	/* lista di valori per l'elevamento a potenza (associa da destra) */
	Stack<TypedValue> powerstack;

	/* operazioni sulle var, se si usa la lista */
	NamedValue* GetVar(const char* const s);
	Var* NewVar(const char* const s, TypedValue::Type t,
			const Real& d = 0.);

	TypedValue::Type GetType(const char* const s) const;
	TypedValue::TypeModifier GetTypeModifier(const char*) const;

	bool IsType(const char* const s) const;
	bool IsTypeModifier(const char* const s) const;
	bool IsKeyWord(NameSpace *ns, const char* const s) const;
   
	/* lexer */
	enum Token GetToken(void);

	/* aumenta il buffer se necessario */
	void IncNameBuf(void);

	/*
	 * funzioni la cui chiamata ricorsiva esprime la precedenza
	 * degli operatori
	 */
	TypedValue logical(void);
	TypedValue relational(void);
	TypedValue binary(void);
	TypedValue mult(void);
	TypedValue unary(void);
	TypedValue power(void);

	/* helper per expr, che valuta le funzioni built-in */
	TypedValue evalfunc(MathParser::NameSpace *ns, MathFunc_t* f);

	/* valuta le espressioni */
	TypedValue expr(void);

	/* valuta gli statement e le liste di statement */
	TypedValue stmt(void);
	TypedValue stmtlist(void);
	TypedValue readplugin(void);
	void trim_arg(char *const s);

public:   
	MathParser(const InputStream& strm, Table& t, int redefine_vars = 0);
	MathParser(Table& t, int redefine_vars = 0);
	~MathParser(void);

	NamedValue *
	InsertSym(const char* const s, const Real& v, int redefine = 0);
	NamedValue *
	InsertSym(const char* const s, const Int& v, int redefine = 0);   
     
	/*
	 * interpreta una sequenza di stmt e restituisce il valore
	 * dell'ultimo
	 */
	Real GetLastStmt(Real d = 0., Token t = ARGSEP);
	Real GetLastStmt(const InputStream& strm, Real d = 0.,
			Token t = ARGSEP);

	/* interpreta uno stmt e ne restitutisce il valore */
	Real Get(Real d = 0.);
	Real Get(const InputStream& strm, Real d = 0.);

	/* modalita' calcolatrice: elabora e stampa ogni stmt */
	void GetForever(std::ostream& out, const char* const sep = "\n");
	void GetForever(const InputStream& strm, std::ostream& out,
			const char* const sep = "\n");

	int RegisterPlugIn(const char *name,
			MathParser::PlugIn * (*)(MathParser&, void *),
			void *arg);

	int RegisterNameSpace(NameSpace *ns);
};

#endif /* MATHP_H */

