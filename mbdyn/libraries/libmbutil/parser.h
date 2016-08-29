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

/* Parser per l'ingresso dati - parte generale */

/* Si compone di tre diverse strutture di scansione, 
 * piu' le strutture di memorizzazione di stringhe e variabili.
 * 
 * La prima struttura e' data dal LowParser, che riconosce gli elementi 
 * della sintassi. E' data da:
 * <statement_list>::=
 *   <statement> ; <statement_list>
 *   epsilon
 * <statement>::=
 *   <description>
 *   <description> : <arg_list>
 * <arg_list>::=
 *   <arg>
 *   <arg> , <arg_list>
 * <arg>::=
 *   <word>
 *   <number>
 * ecc. ecc.
 * 
 * La seconda struttura e' data dall'HighParser, che riconosce la sintassi 
 * vera e propria. In alternativa al LowParser, qualora sia atteso un valore
 * numerico esprimibile mediante un'espressione regolare 
 * (espressione matematica), e' possibile invocare la terza struttura,
 * il MathParser. Questo analizza espressioni anche complesse e multiple,
 * se correttamente racchiuse tra parentesi.
 * 
 * L'HighParser deve necessariamente riconoscere una parola chiave nel campo
 * <description>, mentre puo' trovare parole qualsiasi nel campo <arg>
 * qualora sia attesa una stringa.
 * 
 * Le parole chiave vengono fornite all'HighParser attraverso la KeyTable, 
 * ovvero una lista di parole impaccate (senza spazi). L'uso consigliato e':
 * 
 *   const char sKeyWords[] = { "keyword0",
 *                              "keyword1",
 *                              "...",
 *                              "keywordN"};
 * 
 *   enum KeyWords { KEYWORD0 = 0,
 *                   KEYWORD1,
 *                   ...,
 *                   KEYWORDN,
 *                   LASTKEYWORD};
 * 
 *   KeyTable K((int)LASTKEYWORD, sKeyWords);
 * 
 * Il MathParser usa una tabella di simboli, ovvero nomi (dotati di tipo) 
 * a cui e' associato un valore. La tabella e' esterna e quindi puo' essere 
 * conservata ed utilizzata in seguito conservando in memoria i nomi
 * definiti in precedenza.
 * 
 * A questo punto si puo' generare la tabella dei simboli:
 * 
 *   int iSymbolTableInitialSize = 10;
 *   Table T(iSymbolTableInitialSize);
 * 
 * Quindi si crea il MathParser:
 * 
 *   MathParser Math(T);
 * 
 * Infine si genera l'HighParser:
 * 
 *   HighParser HP(Math, K, StreamIn);
 * 
 * dove StreamIn e' l'istream da cui avviene la lettura.
 */
            

#ifndef PARSER_H
#define PARSER_H

#include <iostream>
#include <fstream>
#include <limits>
#include <ac/f2c.h>

#include <string.h>
#include <ctype.h>

#include <stdlib.h>
#include <unistd.h>

#include "myassert.h"
#include "input.h"
#include "mathp.h"
#include "matvec3.h"
#include "matvec3n.h"
#include "matvec6.h"
#include "mbsleep.h"
#include "ltstrcase.h"


/* Classi dichiarate */
class LowParser;
class KeyTable;
class HighParser;


const unsigned int iDefaultBufSize =
#ifdef BUFSIZ
	BUFSIZ
#else /* ! BUFSIZ */
	8192
#endif /* ! BUFSIZ */
;


/* LowParser - begin */

class LowParser {
	friend class HighParser;

public:
	enum Token {
		UNKNOWN,

		WORD,
		COMMA = ',',
		COLON = ':',
		SEMICOLON = ';',
		NUMBER,
		ENDOFFILE,

		LASTTOKEN
	};
   
private:
	HighParser& HP;
	enum Token CurrToken;
	char *sCurrWordBuf;
	unsigned iBufSize;
	doublereal dCurrNumber;

	void PackWords(InputStream& In);
   
public:
	LowParser(HighParser& hp);
	~LowParser(void);
	Token GetToken(InputStream& In);
	doublereal dGetReal(void) const;
	integer iGetInt(void) const;
	char* sGetWord(void);
};

/* LowParser - end */


/* KeyTable - begin */

class KeyTable {
private:
	char* const* sKeyWords;
	const KeyTable *oldKey;
	HighParser& HP;
   
public:      
	KeyTable(HighParser& hp, const char* const sTable[]);   
	virtual ~KeyTable(void);
	int Find(const char* sToFind) const;
};

/* KeyTable - end */


/* DescRead - begin */

/* prototype of the functional object: reads a description */
struct DescRead {
public:
	virtual ~DescRead(void);
	virtual bool Read(HighParser& HP) = 0;
};

/* description registration function: call to register one */
extern bool
SetDescData(const std::string &s, DescRead *rf);

/* function that reads a description */
extern bool
ReadDescription(HighParser& HP, const std::string& desc);

/* DescRead - end */


/* HighParser - begin */

class HighParser {
public:   
	class ErrInvalidCallToGetDescription : public MBDynErrBase {
	public:
		ErrInvalidCallToGetDescription(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrKeyWordExpected : public MBDynErrBase {
	public:
		ErrKeyWordExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrSemicolonExpected : public MBDynErrBase {
	public:
		ErrSemicolonExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrColonExpected : public MBDynErrBase {
	public:
		ErrColonExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrMissingSeparator : public MBDynErrBase {
	public:
		ErrMissingSeparator(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrIntegerExpected : public MBDynErrBase {
	public:
		ErrIntegerExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrRealExpected : public MBDynErrBase {
	public:
		ErrRealExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrStringExpected : public MBDynErrBase {
	public:
		ErrStringExpected(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrIllegalDelimiter : public MBDynErrBase {
	public:
		ErrIllegalDelimiter(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrValueOutOfRange : public MBDynErrBase {
	public:
		ErrValueOutOfRange(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

public:      
	enum Token {
		UNKNOWN = -1,
		
		DESCRIPTION,
		FIRSTARG,
		ARG,
		LASTARG,
		NOARGS,
		WORD,	
		NUMBER,
		STRING,
		ENDOFFILE,

		LASTITEM
	};
   
	enum Delims {
		UNKNOWNDELIM = -1,
		
		PLAINBRACKETS,
		SQUAREBRACKETS,
		CURLYBRACKETS,
		SINGLEQUOTE,
		DOUBLEQUOTE,
		DEFAULTDELIM,

		LASTDELIM
	};
   
	const char ESCAPE_CHAR;
   
public:
	struct ErrOut {
		const char* sFileName;
		const char* sPathName;
		unsigned int iLineNumber;
	};

	struct WordSet {
		virtual ~WordSet(void) { NO_OP; };
		virtual bool IsWord(const std::string &s) const = 0;
	};
   
	enum {
		NONE		= 0x00U,
		EATSPACES	= 0x01U,
		ESCAPE		= 0x02U,
		LOWER		= 0x04U,
		UPPER		= 0x08U
	};

	template <class T>
	struct range_base {
		virtual bool check(const T& value) const { return false; };
	};

	template <class T>
	struct range_any : public range_base<T> {
		virtual bool check(const T& value) const { return true; };
	};

	template<class T>
	struct range_ge : public range_base<T> {
		T m_lower;
		range_ge(const T& l) : m_lower(l) {};

		bool check(const T& value) const {
			return (value >= m_lower);
		};
	};

	template<class T>
	struct range_gt : public range_base<T> {
		T m_lower;
		range_gt(const T& l) : m_lower(l) {};

		bool check(const T& value) const {
			return (value > m_lower);
		};
	};

	template<class T>
	struct range_le : public range_base<T> {
		T m_upper;
		range_le(const T& u) : m_upper(u) {};

		bool check(const T& value) const {
			return (value <= m_upper);
		};
	};

	template<class T>
	struct range_lt : public range_base<T> {
		T m_upper;
		range_lt(const T& u) : m_upper(u) {};

		bool check(const T& value) const {
			return (value < m_upper);
		};
	};

	template<class T>
	struct range_2_base : public range_base<T> {
		T m_lower, m_upper;
		range_2_base(const T& l, const T& u) : m_lower(l), m_upper(u) {};
	};

	template<class T>
	struct range_ge_le : public range_2_base<T> {
		range_ge_le(const T& l, const T& u) : range_2_base<T>(l, u) {};

		bool check(const T& value) const {
			return ((value >= range_2_base<T>::m_lower) && (value <= range_2_base<T>::m_upper));
		};
	};

	template<class T>
	struct range_gt_le : public range_2_base<T> {
		range_gt_le(const T& l, const T& u) : range_2_base<T>(l, u) {};

		bool check(const T& value) const {
			return ((value > range_2_base<T>::m_lower) && (value <= range_2_base<T>::m_upper));
		};
	};

	template<class T>
	struct range_ge_lt : public range_2_base<T> {
		range_ge_lt(const T& l, const T& u) : range_2_base<T>(l, u) {};

		bool check(const T& value) const {
			return ((value >= range_2_base<T>::m_lower) && (value < range_2_base<T>::m_upper));
		};
	};

	template<class T>
	struct range_gt_lt : public range_2_base<T> {
		range_gt_lt(const T& l, const T& u) : range_2_base<T>(l, u) {};

		bool check(const T& value) const {
			return ((value > range_2_base<T>::m_lower) && (value < range_2_base<T>::m_upper));
		};
	};

protected:
	/* Parser di basso livello, per semplice lettura dei tipi piu' comuni */
	LowParser LowP;
     
	/* Stream in ingresso */
	InputStream* pIn;
	std::ifstream* pf;

	/* Buffer per le stringhe */
	char sStringBuf[iDefaultBufSize];
	char sStringBufWithSpaces[iDefaultBufSize];

	/* Parser delle espressioni matematiche, 
	 * usato per acquisire valori sicuramente numerici */
	MathParser& MathP;

	/* Tabella dei simboli da riconoscere; puo' essere cambiata */
	const KeyTable* KeyT;

	/* Token di basso ed alto livello */
	LowParser::Token CurrLowToken;
	Token CurrToken;

	virtual HighParser::Token FirstToken(void);
	virtual void NextToken(const char* sFuncName);

	int iGetDescription_int(const char* const s);
	virtual void Eof(void);

	virtual void
	SetDelims(enum Delims Del, char &cLdelim, char &cRdelim) const;

	int ParseWord(unsigned flags = HighParser::NONE);
	void PutbackWord(void);

public:   
	HighParser(MathParser& MP, InputStream& streamIn);
	virtual ~HighParser(void);
	/* Attacca una nuova KeyTable (e ritorna la vecchia) */
	virtual const KeyTable* PutKeyTable(const KeyTable& KT);
	/* Numero di linea corrente */   
	virtual int GetLineNumber(void) const;
	/* Numero di nome file e linea corrente */   
	virtual HighParser::ErrOut GetLineData(void) const;
	/* Restituisce il math parser */
	virtual MathParser& GetMathParser(void);
	/* "Chiude" i flussi */
	virtual void Close(void);
	/* verifica se il token successivo e' una description (ambiguo ...) */
	bool IsDescription(void) const;
	/* ha appena trovato una description */
	Token GotDescription(void);
	/* Legge una parola chiave */
	int GetDescription(void);
	/* si attende una descrizione */
	virtual void ExpectDescription(void);
	/* si attende una lista di argomenti */
	virtual void ExpectArg(void);
	/* 1 se trova la keyword sKeyWord */
	virtual bool IsKeyWord(const char* sKeyWord);
	/* numero della keyword trovata */
	virtual int IsKeyWord(void);
	/* 1 se e' atteso un argomento */
	virtual bool IsArg(void);
	/* 1 se e' atteso un argomento */
	virtual bool IsStringWithDelims(enum Delims Del = DEFAULTDELIM);
	/* se l'argomento successivo e' una parola in un WordSet, la ritorna */
	virtual const char *IsWord(const HighParser::WordSet& ws);
	/* Se ha letto un ";" lo rimette a posto */
	virtual void PutBackSemicolon(void);
	/* legge un booleano con il mathpar */
	virtual bool GetBool(bool bDefval = false);
	/* legge un "yes"/"no" */
	virtual bool GetYesNo(bool& bRet);
	/* legge un "yes"/"no" o booleano con il mathpar */
	virtual bool GetYesNoOrBool(bool bDefval = false);
	/* legge un intero con il mathpar */
	virtual integer GetInt(integer iDefval = 0);
	template <class Range>
	integer GetInt(integer iDefval, Range r);
	/* legge un reale con il mathpar */
	virtual doublereal GetReal(const doublereal& dDefval = 0.0);
	template <class Range>
	doublereal GetReal(const doublereal& dDefval, Range r);
	/* legge una string con il mathpar */
	virtual std::string GetString(const std::string& sDefVal);
	/* legge un valore tipizzato con il mathpar */
	virtual TypedValue GetValue(const TypedValue& v);
	template <class Range>
	TypedValue GetValue(const TypedValue& v, Range r);
	/* legge un timeout */
	virtual mbsleep_t GetTimeout(const mbsleep_t& DefVal);
	/* legge una keyword */
	virtual int GetWord(void);
	/* legge una stringa */
	virtual const char* GetString(unsigned flags = HighParser::NONE);
	/* stringa delimitata */
	virtual const char* GetStringWithDelims(enum Delims Del = DEFAULTDELIM, bool escape = true); 
};

/* Le funzioni:
 *   ExpectDescription()
 *   ExpectArg()
 * informano il parser di cio' che e' atteso; di default il costruttore
 * setta ExpectDescription().
 * 
 * Le funzioni:
 *   GetDescription()
 *   IsKeyWord()
 *   GetWord()
 * restituiscono un intero che corrisponde alla posizione occupata nella
 * KeyTable dalle parole corrispondenti, oppure -1 se la parola non e'
 * trovata. Si noti che IsKeyWord(), in caso di esito negativo, ripristina 
 * l'istream. Tutte preparano poi il parser per la lettura successiva.
 * 
 * La funzione 
 *   IsKeyWord(const char*) 
 * restituisce 0 se non trova la parola e ripristina l'istream, altrimenti 
 * restituisce 1 e prepara il parser alla lettura successiva.
 * 
 * Le funzioni 
 *   GetInt(), 
 *   GetReal(), 
 *   GetString(), 
 *   GetStringWithDelims(enum Delims)
 * restituiscono i valori attesi e preparano il prser alla lettura successiva.
 */

/* HighParser - end */
   
extern std::ostream&
operator << (std::ostream& out, const HighParser::ErrOut& err);

extern HighParser::ErrOut
mbdyn_get_line_data(void);

extern std::ostream&
mbdyn_print_line_data(std::ostream& out);

template <class Range>
integer
HighParser::GetInt(int iDefVal, Range range)
{
	TypedValue v(iDefVal);
	v = GetValue(v);
	integer val = v.GetInt();
	if (!range.check(val)) {
		throw ErrValueOutOfRange(MBDYN_EXCEPT_ARGS);
	}
	return val;
}


template <class Range>
doublereal
HighParser::GetReal(const doublereal& dDefVal, Range range)
{
	TypedValue v(dDefVal);
	v = GetValue(v);
	doublereal val = v.GetReal();
	if (!range.check(val)) {
		throw ErrValueOutOfRange(MBDYN_EXCEPT_ARGS);
	}
	return val;
}

template <class Range>
TypedValue
HighParser::GetValue(const TypedValue& vDefVal, Range range)
{
	const char sFuncName[] = "HighParser::GetValue()";

	if (CurrToken != HighParser::ARG) {
		silent_cerr("Parser error in "
			<< sFuncName << ", arg expected at line "
			<< GetLineData() << std::endl);
		throw HighParser::ErrIntegerExpected(MBDYN_EXCEPT_ARGS);
	}

	TypedValue v(vDefVal);

	try {
		v = MathP.Get(*pIn, v);
	}
	catch (TypedValue::ErrWrongType& e) {
		silent_cerr(sFuncName << ": " << e.what() << " at line "
			<< GetLineData() << std::endl);
		throw e;
	}
	catch (MathParser::ErrGeneric& e) {
		silent_cerr(sFuncName << ": error return from MathParser at line "
			<< GetLineData() << std::endl);
		throw e;
	}
	catch (...) {
		throw;
	}

	NextToken(sFuncName);

	if (!range.check(v)) {
		throw ErrValueOutOfRange(MBDYN_EXCEPT_ARGS);
	}

	return v;
}

#endif /* PARSER_H */

