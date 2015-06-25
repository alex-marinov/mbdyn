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


#ifndef PARSINC_H
#define PARSINC_H

#include <fstream>
#include <stack>

#include "parser.h"

/* IncludeParser - begin */

class IncludeParser : public HighParser {
	friend struct IncludeDR;
protected:

	/*
	 * Struttura e dati per la stack di flussi in ingresso, usata
	 * per consentire l'inclusione multipla di files
	 * ATTENZIONE: non impedisce l'inclusione ricorsiva, e quindi il loop
	 */
	struct MyInput {
		std::ifstream* pfile;
		InputStream* pis;

#ifdef USE_INCLUDE_PARSER
		char* sPath;
		char* sFile;

		MyInput(std::ifstream* pf = NULL, InputStream* pi = NULL,
			char* sp = NULL, char* sfile = NULL)
		: pfile(pf), pis(pi), sPath(sp), sFile(sfile) {
			NO_OP;
		};

#else /* !USE_INCLUDE_PARSER */
		MyInput(std::ifstream* pf = NULL, InputStream* pi = NULL)
		: pfile(pf), pis(pi) {
			NO_OP;
		};
#endif /* !USE_INCLUDE_PARSER */
	};

	std::stack<MyInput *> myinput;
	char* sCurrPath;
	char* sInitialPath;
	char* sCurrFile;

	flag fCheckStack(void);
	bool Include_int(void);
	virtual void Eof(void);

public:
	IncludeParser(MathParser& MP, InputStream& streamIn,
		const char *initial_file = "initial file");
	virtual ~IncludeParser(void);

	virtual void Close(void);                  /* "Chiude" i flussi */

	virtual const char* GetFileName(enum Delims Del = DEFAULTDELIM);

	virtual HighParser::ErrOut GetLineData(void) const;
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
 *   GatFileName(enum Delims)
 * restituiscono i valori attesi e preparano il parser alla lettura successiva.
 */

/* IncludeParser - end */

#endif /* PARSINC_H */

