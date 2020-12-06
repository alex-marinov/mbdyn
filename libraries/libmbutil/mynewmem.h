/* $Header$ */
/*
 * This library comes with MBDyn (C), a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

/******************************************************************************

	Allocazione dinamica di memoria controllata con ASSERT e Memory Manager


	Scritto da
	Pierangelo Masarati
	il 05/04/1997


Se si usa #define DEBUG sono attivate solo le routines di ASSERT;

Se si usa anche #define DEBUG_MEMMANAGER, viene attivato il memory manager.

Se si aggiunge #define DEBUG_MEMMANAGER_DELETEBUTKEEPLOG, quando un blocco
viene cancellato, ne rimane traccia nel log.

Se si aggiunge #define DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE, quando un blocco
viene cancellato in realta' la memoria non viene mai rilasciata, solo se ne
perde il riferimento


E' necessario aggiungere al makefile i seguenti files:

	<myassert.cc>
	<mynewmem.cc>

Inoltre e' necessario mettere questa intestazione in tutti i files che
usano la gestione sicura della memoria:

	#define DEBUG
	#define DEBUG_MEMMANAGER
	[ #define DEBUG_MEMMANAGER_DELETEBUTKEEPLOG ||
		#define DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE ]
	#include "mynewmem.h"



Le macro di gestione dinamica della memoria sono:

	SAFENEW(item* pit, item);
	SAFENEW_(item* pit, item, clMemMan mm);
		<pit>  e' un puntatore al tipo <item> e deve essere <pit> = NULL
		<item> e' un tipo valido
		<mm>   e' un oggetto della classe clMemMan
					 se #ifdef DEBUG_MEMMANAGER = TRUE,
					 altrimenti non viene usato
        ###
        ### ATTENZIONE: QUANDO SI USA LA MACRO "SAFENEW", O <item> E' UN TIPO
        ### SEMPLICE, IL CUI COSTRUTTORE NON CREA EFFETTI COLLATERALI
        ### (ALLORA SI PUO' USARE IL COSTRUTTORE AL POSTO DI <item>), 
        ### OPPURE SI DESIDERA UTILIZZARE IL COSTRUTTORE DI DEFAULT.
        ### SE SI DEVE UTILIZZARE UN COSTRUTTORE COMPLESSO, E' NECESSARIO
        ### USARE LA MACRO "SAFENEWWITHCONSTRUCTOR" 
        ### (ED E' CONSIGLIABILE USARLA SEMPRE PER TIPI COMPLESSI, ALTRIMENTI 
        ### PUO' ESSERE ESEGUITO PIU' VOLTE IL COSTRUTTORE
        ###

	SAFENEWWCONS(item* pit, item, constructor);
	SAFENEWWITHCONSTRUCTOR_(item* pit, item, constructor, clMemMan mm);
		<pit>  e' un puntatore al tipo <item> e deve essere <pit> = NULL
		<item> e' un tipo valido
                <constructor> ' il costruttore di <item>
		<mm>   e' un oggetto della classe clMemMan
					 se #ifdef DEBUG_MEMMANAGER = TRUE,
					 altrimenti non viene usato

	SAFENEWARR(item* pit, item, type size);
	SAFENEWARR_(item* pit, item, type size, clMemMan mm);
		come SAFENEW, solo che viene allocato un array di oggetti
		di tipo <item> dalle dimensioni determinate da <size>, 
		di tipo indefinito.



	SAFEDELETE(item* pit);
	SAFEDELETE_(item* pit, clMemMan mm);
		come sopra; <pit> deve essere definito e deve essere 
		contenuto in <mm> in qualita' di puntatore all'inizio 
		di un blocco se #ifdef DEBUG_MEMMANAGER = TRUE

	SAFEDELETEARR(item* pit);
	SAFEDELETEARR_(item* pit, clMemMan mm);
		come SAFEDELETE, solo che viene cancellata un'array di oggetti.

	SAFEDELETEANDFILLARR(item* pit);
	SAFEDELETEANDFILLARR_(item* pit, clMemMan mm);
		come SAFEDELETEARR, ma se #ifdef DEBUG_MEMMANAGER = TRUE,
		prima di deallocare riempe la memoria con il carattere di free;
		deve essere usato solo quando ogni oggetto dell'array e' stato
                distrutto con l'apposito distruttore, oppure se gli oggetti
		non hanno distruttore (verificare come viene distrutto un 
		array di oggetti con distruttore)


Uso:

		item* pit = NULL;
		clMemMan m;

		// per allocare:
		SAFENEW(pit, item, m);

		if( pit == NULL) // Handle error
		// L'errore puo' essere gestito in modo migliore attraverso
		// il casting dell'operatore <new> in modo che non sia 
		// possibile mancare l'allocazione. 
		// In ogni caso la mancata allocazione viene segnalata 
		// su <cerr>.

		// per deallocare:
		SAFEDELETE(pit, m);

		// In modo analogo si usano
		SAFENEWARR(pit, item, size, m);
		SAFEDELETEARR(pit, m);
		SAFEDELETEANDFILLARR(pit, m);


Nota:   per comportamenti particolari, agire sul codice
	delle funzioni clMemMan::remove, clMemMan::removeAndFill
	e delle macro SAFEDELETE(), SAFEDELETEARR()
	per evitare il delete fisico della memoria o per mantenere un log
	delle allocazioni avvenute.



class clMemMan :

La classe clMemMan e' un importante strumento di debugging.
Consente di mantenere informazioni relative alla memoria allocata
attraverso le macro e quindi di verificare in qualsiasi momento
se la memoria utilizzata e' effettivamente disponibile.

L'inserimento e la cancellazione dei record avviene automaticamente
tramite le macro. L'utente dispone di due strumenti di verifica.
Il primo consiste nel dump della struttura della memoria. In esso
i records sono ordinati in base al valore del puntatore a cui fanno
riferimento. Questo rende possibile una verifica manuale.
Il secondo strumento consiste nella possibilita' di marcare
i puntatori desiderati con un riferimento qualora siano	effettiva-
mente	disponibili. Quindi un dump dei riferimenti marcati consente
di verificare sia la presenza di garbage (memoria non piu' puntata
e quindi non referenziabile) che di dangling referencies (puntatori
a memoria non piu' disponibile).



Le funzioni pubbliche della classe sono:

flag clMemMan::fIsBlock(const void* pvBlock, size_t sizeBlock = 1) const;
	se <pvBlock> e' un puntatore all'inizio di un blocco
	e <sizeBlock> e' la corretta dimensione del blocco ritorna TRUE,
	altrimenti FALSE

flag clMemMan::fIsPointerToBlock(const void* pvBlock) const;
	se <pvBlock> e' un puntatore all'inizio di un blocco
	ritorna TRUE,	altrimenti FALSE

flag clMemMan::fIsValid(const void* pvValid, size_t sizeValid = 1) const;
	se <pvValid> e' una valida locazione all'interno di un blocco e
	<(void*)pvalid+sizeValid-1> non eccede il blocco ritorna TRUE,
	altrimenti FALSE

size_t clMemMan::sizeOfBlock(const void* pvSizeOf) const;
	se <pvSizeOf> e' un puntatore all'inizio di un blocco,
	ne ritorna la	dimensione

flag clMemMan::fIsArray(const void* pvIsArray) const;
	se <pvIsArray> e' un puntatore all'inizio di un blocco,
	ritorna TRUE se e' un'array, altrimenti FALSE

eStatus clMemMan::eBlockStatus(const void* pvBStatus) const;
	se <pvIsArray> e' un puntatore all'inizio di un blocco,
	ne ritorna lo stato con un valore appartenente all'enum <eStatus>

void clMemMan::ClearRefs(void);
	azzera i flag di riferimento a tutti i puntatori

void clMemMan::PutRef(const void* pvRef);
	se <pvRef> e' un puntatore all'inizio di un blocco,
	ne setta il riferimento

flag clMemMan::fIsRefd(const void* pvIsRefd) const;
	se <pvRef> e' un puntatore all'inizio di un blocco,
	ritorna TRUE se il riferimento e' settato

void clMemMan::DumpRef(ostream& rout) const;
	invia la stampa dei blocchi suddivisi in base ai riferimenti
	ad un ostream

ostream& operator << (ostream& rout, const clMemMan& rm);
	invia la stampa dei blocchi con le relative proprieta' ad un ostream


Esempio completo:

			// ...

		#define DEBUG
		#define DEBUG_MEMMANAGER
		#include "mynewmem.h"

			// ...

		// Uso dei riferimenti:

		#ifdef DEBUG_MEMMANAGER
		clMemMan m("Test"); // Memory Manager
		#endif

		void* pb = NULL;
		SAFENEW(pb, char, m); // Allocazione di un byte

		#ifdef DEBUG_MEMMANAGER
		m.ClearRefs();       		// Pulisce i Ref
		m.PutRef(pb);        		// Setta il Ref di pb
		if (m.fIsRefd(pb)) { 		// Ritorna TRUE
			m.DumpRef(cout);  	// scrive i dati
		}
		#endif

		SAFEDELETE(pb, m);    	// Deallocazione
		if (pb) {              	// Ritorna FALSE, perche' SAFEDELETE
			NULL;           // pone comunque pb = NULL
		}

		// Uso delle utilities:

		SAFENEWARR(pb, char, 10, m); // Allocazione di 10 bytes

		#ifdef DEBUG_MEMMANAGER
		if (m.fIsPointerToBlock(pb), 10) {	// Ritorna TRUE
			size_t size = m.sizeOfBlock(pb);// Dim. del blocco
			int iItems = 1;
			if (m.fIsArray(pb)) {	// Se array calcola le dim.
				iItems = size/sizeof(*pb);
				cout << "\nNumero di elementi: " 
					<< iItems << endl;
			}
			void* pbTmp = pb;
			while (size--) {	// Azzera la memoria
						// (operazione da non fare
				*pbTmp++ = '\0';// se devono essere eseguiti
						// i distruttori
			}
		}
		#endif
		SAFEDELETEARR(pb, m);          	// Libera la memoria


Si noti come, togliendo #define DEBUG_MEMMANAGER, si abbia il codice
di debug semplificato, con gli ASSERT ed altri controlli, mentre,
togliendo anche #define DEBUG, si abbia in modo autometico il codice
in versione release, senza alcun overhead

******************************************************************************/

#ifndef MYNEWMEM_H
#define MYNEWMEM_H


#include <myassert.h>
#include <except.h>

#ifdef DEBUG

/* Dichiarazione di funzioni usate nella versione debug */
extern void _Safenew(const char* file, int line, int flag = 0);
extern void _Safenewfill(void* pv, size_t size, char fill);

const char cDebugAlloc = 'A';
const char cDebugFree  = 'F';

#ifdef DEBUG_MEMMANAGER

/* Dichiarazione delle macro di debug con memory manager */
#define SAFENEW_(pnt, item, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		try{(pnt) = new item;} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__); \
			throw(ba); \
		} \
		(memman).add((void*)(pnt), sizeof(item)); \
	} while(0)

/* Dichiarazione delle macro di debug con memory manager */
#define SAFENEWWITHCONSTRUCTOR_(pnt, item, constructor, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		try{(pnt) = new constructor;} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__); \
			throw(ba); \
           	} \
           	(memman).add((void*)(pnt), sizeof(item)); \
	} while (0)

/* Attenzione: questa operazione e' lecita solo se
 * non e' stato eseguito un costruttore
 * _Safenewfill(pnt, sizeof(item), cDebugAlloc); */
#define SAFENEWARR_(pnt, item, sz, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		ASSERT((sz) != 0);	   \
		try{(pnt) = new item[sz];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
		(memman).add((void*)(pnt), sizeof(item)*(sz), 1); \
		_Safenewfill((void*)pnt, sizeof(item)*(sz), cDebugAlloc); \
	} while (0)

/* questa e' sicura anche se e' stato eseguito un costruttore */
#define SAFENEWARRNOFILL_(pnt, item, sz, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		ASSERT(sz); \
		try{(pnt) = new item[sz];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
		(memman).add((void*)(pnt), sizeof(item)*(sz), 1); \
	} while (0)

#define SAFESTRDUP_(pnt, src, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT((src)); \
		unsigned int l = strlen((src))+1; \
		try{(pnt) = new char[l];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
		strcpy((char *)(pnt), (src)); \
		(memman).add((void*)(pnt), l, 1); \
	} while (0)

/* Qui il fill e' lecito perche' per le arrays
 * non sono eseguiti i costruttori */
 
#ifndef DEBUG_MEMMANAGER_DELETEBUTKEEPLOG

#ifndef DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE

#define SAFEDELETE_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).remove((void*)(pnt)); \
		delete (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).remove((void*)(pnt), 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEANDFILLARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).remove((void*)(pnt), 1, 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#endif /* !DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE */

#endif /* !DEBUG_MEMMANAGER_DELETEBUTKEEPLOG */

#ifdef DEBUG_MEMMANAGER_DELETEBUTKEEPLOG

#define SAFEDELETE_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButKeepLog((void*)(pnt)); \
		delete (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButKeepLog((void*)(pnt), 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEANDFILLARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButKeepLog((void*)(pnt), 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#endif /* DEBUG_MEMMANAGER_DELETEBUTKEEPLOG */

#ifdef DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE

#define SAFEDELETE_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButNotRelease((void*)(pnt)); \
		delete (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButNotRelease((void*)(pnt), 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEANDFILLARR_(pnt, memman) \
	do { \
		ASSERT(pnt); \
		(memman).removeButNotRelease((void*)(pnt), 1); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#endif /* DEBUG_MEMMANAGER_DELETEBUTNOTRELEASE */

/* enum degli stati della memoria */
enum eStatus { 
	UNKNOWN, 
	ALLOCATED, 
	FREED, 
	FREEDBUTNOTRELEASED 
};

/* struttura dei dati di un blocco di memoria */
struct stMemBlock {
	void* pv;
	size_t size;
	eStatus eSt;
	flag fArr;
	flag fRef;
	
	stMemBlock(void* _pv = NULL, size_t _size = 0, eStatus _eSt = UNKNOWN, 
		flag _fArr = 0, flag _fRef = 0)
		: pv(_pv), size(_size), eSt(_eSt), fArr(_fArr), fRef(_fRef) {
	    	NO_OP;
	};
};

/* classe del memory manager */
class clMemMan {
	friend std::ostream& operator << (std::ostream& rout, const clMemMan& rm);

public:
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrNotFound : public MBDynErrBase {
	public:
		ErrNotFound(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

private:
	struct stList {
		stMemBlock stMB;
		stList* pstNext;
		stList(stMemBlock stIn) : stMB(stIn), pstNext(NULL) {
			NO_OP;
		};
	};
	
	stList* pstRoot;
	char* sName;
	
	stList* pstFindElem(const void* pvToFind) const;
	stList* pstFindPrev(const void* pvToFindPrev) const;
	
	enum eRemoveMode { 
		RELEASE,
		DELBUTKEEPLOG,
		DELBUTNOTRELEASE
	};
	
	void _remove(const void* pvToRemove, eRemoveMode eMode, flag fArr, 
		flag fFill);
		
public:
	clMemMan(char* sName = NULL);
	~clMemMan(void);
       
	flag fIsBlock(const void* pvBlock, size_t sizeBlock = 1) const;
	flag fIsPointerToBlock(const void* pvBlock) const;
	flag fIsValid(const void* pvValid, size_t sizeValid = 1) const;
	size_t sizeOfBlock(const void* pvSizeOf) const;
	flag fIsArray(const void* pvIsArray) const;
	eStatus eBlockStatus(const void* pvBStatus) const;
       
	void ClearRefs(void);
	void PutRef(const void* pvRef);
	flag fIsRefd(const void* pvIsRefd) const;
	std::ostream& DumpRef(std::ostream& rout) const;
       
	void add(const void* pvIn, size_t sizeIn, flag fArr = 0);
	inline void remove(const void* pvToRemove, flag fArr = 0, flag fFill = 0) {
		_remove(pvToRemove, RELEASE, fArr, fFill);
	};
       
	inline void removeButNotRelease(const void* pvToRemove, flag fArr = 0, flag fFill = 0) {
		_remove(pvToRemove, DELBUTNOTRELEASE, fArr, fFill); 
	};
       
	inline void removeButKeepLog(const void* pvToRemove, flag fArr = 0, flag fFill = 0) {
		_remove(pvToRemove, DELBUTKEEPLOG, fArr, fFill);
	};
};

extern clMemMan defaultMemoryManager;

/* funzione di stream del memory manager */
extern std::ostream& operator << (std::ostream& rout, const clMemMan& rm);

#else /* !DEBUG_MEMMANAGER */

/* dichiarazione delle macro di debug senza memory manager */
#define SAFENEW_(pnt, item, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		try{(pnt) = new item;} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__); \
			throw(ba); \
		} \
	} while (0)

/* Dichiarazione delle macro di debug senza memory manager */
#define SAFENEWWITHCONSTRUCTOR_(pnt, item, constructor, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		try{(pnt) = new constructor;} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__); \
			throw(ba); \
		} \
	} while (0)

/* Attenzione: questa operazione e' lecita solo se
 * non e' stato eseguito un costruttore
 * _Safenewfill(pnt, sizeof(item), cDebugAlloc) */

#define SAFENEWARR_(pnt, item, sz, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		ASSERT((sz) != 0);		   \
		try{(pnt) = new item[sz];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
		_Safenewfill(pnt, sizeof(item)*(sz), cDebugAlloc); \
	} while (0) 

/* questa e' sicura anche se e' stato eseguito un costruttore */
#define SAFENEWARRNOFILL_(pnt, item, sz, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT(sizeof(item)); \
		ASSERT(sz); \
		try{(pnt) = new item[sz];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
	} while (0) 

#define SAFESTRDUP_(pnt, src, memman) \
	do { \
		ASSERT(!(pnt)); \
		ASSERT((src)); \
		unsigned int l = strlen((src))+1; \
		try{(pnt) = new char[l];} \
		catch(std::bad_alloc& ba) { \
			_Safenew(__FILE__, __LINE__, 1); \
			throw(ba); \
		} \
		strcpy((char *)(pnt), (src)); \
	} while (0)

/* Qui il fill e' lecito perche' per le arrays
 * non sono eseguiti i costruttori */

#define SAFEDELETE_(pnt, memman) \
	do { \
		ASSERT((pnt) != 0); \
		delete (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEARR_(pnt, memman) \
	do { \
		ASSERT((pnt) != 0); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEANDFILLARR_(pnt, memman) \
	do { \
		ASSERT((pnt) != 0); \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)      

#endif /* !DEBUG_MEMMANAGER */

#else /* !DEBUG */

/* dichiarazione delle macro nella versione da release */
#define SAFENEW_(pnt, item, memman) \
	(pnt) = new item

/* dichiarazione delle macro nella versione da release */
#define SAFENEWWITHCONSTRUCTOR_(pnt, item, constructor, memman) \
	(pnt) = new constructor

#define SAFENEWARR_(pnt, item, sz, memman) \
	(pnt) = new item[sz]

#define SAFENEWARRNOFILL_(pnt, item, sz, memman) \
	(pnt) = new item[sz]

#define SAFESTRDUP_(pnt, src, memman) \
	do { \
		unsigned int l = strlen((src))+1; \
		(pnt) = new char[l]; \
		strcpy((char *)(pnt), (src)); \
	} while (0)

#define SAFEDELETE_(pnt, memman) \
	do { \
		delete (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEARR_(pnt, memman) \
	do { \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define SAFEDELETEANDFILLARR_(pnt, memman) \
	do { \
		delete[] (pnt); \
		(pnt) = NULL; \
	} while (0)

#define defaultMemoryManager	0

#endif /* !DEBUG */

#define SAFENEW(pnt, item) \
	SAFENEW_(pnt, item, defaultMemoryManager)

#define SAFENEWWITHCONSTRUCTOR(pnt, item, constructor) \
	SAFENEWWITHCONSTRUCTOR_(pnt, item, constructor, defaultMemoryManager)

#define SAFENEWARR(pnt, item, sz) \
	SAFENEWARR_(pnt, item, sz, defaultMemoryManager)

#define SAFENEWARRNOFILL(pnt, item, sz) \
	SAFENEWARRNOFILL_(pnt, item, sz, defaultMemoryManager)

#define SAFESTRDUP(pnt, src) \
	SAFESTRDUP_(pnt, src, defaultMemoryManager)

#define SAFEDELETE(pnt) \
	SAFEDELETE_(pnt, defaultMemoryManager)

#define SAFEDELETEARR(pnt) \
	SAFEDELETEARR_(pnt, defaultMemoryManager)

#define SAFEDELETEANDFILLARR(pnt) \
	SAFEDELETEANDFILLARR_(pnt, defaultMemoryManager)

#endif /* MYNEWMEM_H */

