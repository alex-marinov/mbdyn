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
            

#ifndef MBPAR_H
#define MBPAR_H

#include <parsinc.h>

#include <llist.h>
/* Per reference frame */
#if defined(USE_STRUCT_NODES)
#include <reffrm.h>
/* Riferimento assoluto */
extern const ReferenceFrame AbsRefFrame;
#endif /* defined(USE_STRUCT_NODES) */

/* Per nodi idraulici */
#if defined(USE_HYDRAULIC_NODES)
#include <hfluid.h>
#endif /* defined(USE_HYDRAULIC_NODES) */

#if defined(USE_AERODYNAMIC_ELEMS)
#include <aerodata.h>
#endif /* USE_AERODYNAMIC_ELEMS */

#include "constltp.h"

/* Deals with license and disclaimer output */
extern void mbdyn_license(std::ostream& out);
extern void mbdyn_warranty(std::ostream& out);

/* MBDynParser - begin */

class MBDynParser : public IncludeParser {
public:
	class ErrGeneric {};
	class ErrReferenceAlreadyDefined {};
	class ErrReferenceUndefined {};
   
public:
	enum Frame {
		UNKNOWNFRAME = 0,
		
		GLOBAL,
		NODE,
		LOCAL,
		REFERENCE,
		
		LASTFRAME
	};
   
protected:      
	/* Struttura e dati per la linked list di reference frames */
#if defined(USE_STRUCT_NODES)   
	HardDestructor<ReferenceFrame> RFHD;
	MyLList<ReferenceFrame> RF;
	
	Frame GetRef(ReferenceFrame& rf);
	
	void Reference_int(void);
#endif /* USE_STRUCT_NODES */
 
	/* Struttura e dati per la linked list di hydraulic fluids */
#if defined(USE_HYDRAULIC_NODES)   
	HardDestructor<HydraulicFluid> HFHD;
	MyLList<HydraulicFluid> HF;
	
	void HydraulicFluid_int(void);
#endif /* USE_HYDRAULIC_NODES */
 
	/* Struttura e dati per la linked list di c81 data */
#if defined(USE_AERODYNAMIC_ELEMS)   
	HardDestructor<C81Data> ADHD;
	MyLList<C81Data> AD;
	
	void C81Data_int(void);
#endif /* USE_AERODYNAMIC_ELEMS */
	
	HardDestructor<ConstitutiveLaw1D> C1DHD;
	MyLList<ConstitutiveLaw1D> C1D;
	HardDestructor<ConstitutiveLaw3D> C3DHD;
	MyLList<ConstitutiveLaw3D> C3D;
	HardDestructor<ConstitutiveLaw6D> C6DHD;
	MyLList<ConstitutiveLaw6D> C6D;

	void ConstitutiveLaw_int(void);

	/* Drives */
	HardDestructor<DriveCaller> DCHD;
	MyLList<DriveCaller> DC;

	void DriveCaller_int(void);

	/* Legge una parola chiave */
	bool GetDescription_int(const char *s);

	DataManager *pDM;

public:
	MBDynParser(MathParser& MP, InputStream& streamIn,
			const char *initial_file);
	~MBDynParser(void);

	void SetDataManager(DataManager *pdm);
	
	
	/*
	 * Lettura di posizioni, vettori e matrici di rotazione
	 * relative ed assolute rispetto ad un riferimento
	 */
#if defined(USE_STRUCT_NODES)
	Vec3 GetPosRel(const ReferenceFrame& rf);
	Vec3 GetPosAbs(const ReferenceFrame& rf);
	Vec3 GetVelRel(const ReferenceFrame& rf, const Vec3& x);
	Vec3 GetVelAbs(const ReferenceFrame& rf, const Vec3& x);
	Vec3 GetOmeRel(const ReferenceFrame& rf);
	Vec3 GetOmeAbs(const ReferenceFrame& rf);
	Vec3 GetVecRel(const ReferenceFrame& rf);
	Vec3 GetVecAbs(const ReferenceFrame& rf);
	Mat3x3 GetMatRel(const ReferenceFrame& rf);
	Mat3x3 GetMatAbs(const ReferenceFrame& rf);
	Mat3x3 GetRotRel(const ReferenceFrame& rf);
	Mat3x3 GetRotAbs(const ReferenceFrame& rf);

	void OutputFrames(std::ostream& out) const;
#endif /* USE_STRUCT_NODES */
   
#if defined(USE_HYDRAULIC_NODES)
	HydraulicFluid* GetHydraulicFluid(void);
#endif /* USE_HYDRAULIC_NODES */

#if defined(USE_AERODYNAMIC_ELEMS)
	const c81_data* GetC81Data(integer profile);
#endif /* USE_AERODYNAMIC_ELEMS */

	ConstitutiveLaw1D* GetConstLaw1D(ConstLawType::Type& clt);
	ConstitutiveLaw3D* GetConstLaw3D(ConstLawType::Type& clt);
	ConstitutiveLaw6D* GetConstLaw6D(ConstLawType::Type& clt);
	DriveCaller *MBDynParser::GetDriveCaller(void);
};

/*
 * Le funzioni:
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
 *   GetFileName(enum Delims)
 * restituiscono i valori attesi e preparano il parser alla lettura successiva.
 */

/* MBDynParser - end */
   
#endif /* MBPAR_H */

