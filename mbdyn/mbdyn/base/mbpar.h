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


#ifndef MBPAR_H
#define MBPAR_H

#include <stack>

#include "parsinc.h"

/* Per reference frame */
#include "reffrm.h"
/* Riferimento assoluto */
extern const ReferenceFrame AbsRefFrame;

/* Per nodi idraulici */
#include "hfluid.h"

/* Aerodynamic data */
#include "aerodata.h"

#include "constltp.h"
#include "ScalarFunctions.h"

/* Deals with license and disclaimer output */
extern void mbdyn_license(void);
extern void mbdyn_warranty(void);

/* MBDynParser - begin */

class MBDynParser : public IncludeParser {
public:
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrReferenceAlreadyDefined : public MBDynErrBase {
	public:
		ErrReferenceAlreadyDefined(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrReferenceUndefined : public MBDynErrBase {
	public:
		ErrReferenceUndefined(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

public:
	enum Frame {
		UNKNOWNFRAME = 0,

		GLOBAL,
		NODE,
		LOCAL,
		REFERENCE,

		OTHER_POSITION,
		OTHER_ORIENTATION,
		OTHER_NODE,

		LASTFRAME
	};

	/* Struttura e dati per la linked list di reference frames */
public:
	typedef std::map<unsigned, const ReferenceFrame *> RFType;
	const RFType& GetReferenceFrameContainer(void) const;

protected:
	RFType RF;

	Frame GetRef(ReferenceFrame& rf);

	bool Reference_int(void);
friend struct RefFrameDR;

	/* Struttura e dati per la linked list di hydraulic fluids */
public:
	typedef std::map<unsigned, const HydraulicFluid *> HFType;
	const HFType& GetHydraulicFluidContainer(void) const;

protected:
	HFType HF;

	bool HydraulicFluid_int(void);
friend struct HFluidDR;

	/* Struttura e dati per la linked list di c81 data */
public:
	typedef std::map<unsigned, C81Data *> ADType;
	const ADType& GetC81DataContainer(void) const;

protected:
	ADType AD;

	bool C81Data_int(void);
friend struct C81DataDR;

	typedef std::map<unsigned, const ConstitutiveLaw1D *> CL1DType;
	typedef std::map<unsigned, const ConstitutiveLaw3D *> CL3DType;
	typedef std::map<unsigned, const ConstitutiveLaw6D *> CL6DType;
	CL1DType CL1D;
	CL3DType CL3D;
	CL6DType CL6D;

	bool ConstitutiveLaw_int(void);
friend struct ConstLawDR;

	/* Drives */
public:
	typedef std::map<unsigned, const DriveCaller *> DCType;
	const DCType& GetDriveCallerContainer(void) const;

protected:
	DCType DC;

	bool DriveCaller_int(void);
friend struct DriveCallerDR;

public:
	/* Template drives */
	typedef std::map<unsigned, const TplDriveCaller<doublereal> *> DC1DType;
	typedef std::map<unsigned, const TplDriveCaller<Vec3> *> DC3DType;
	typedef std::map<unsigned, const TplDriveCaller<Vec6> *> DC6DType;
	typedef std::map<unsigned, const TplDriveCaller<Mat3x3> *> DC3x3DType;
	typedef std::map<unsigned, const TplDriveCaller<Mat6x6> *> DC6x6DType;

	const DC1DType& GetDriveCaller1DContainer(void) const;
	const DC3DType& GetDriveCaller3DContainer(void) const;
	const DC6DType& GetDriveCaller6DContainer(void) const;
	const DC3x3DType& GetDriveCaller3x3DContainer(void) const;
	const DC6x6DType& GetDriveCaller6x6DContainer(void) const;

protected:
	DC1DType DC1D;
	DC3DType DC3D;
	DC6DType DC6D;
	DC3x3DType DC3x3D;
	DC6x6DType DC6x6D;

	bool TplDriveCaller_int(void);
friend struct TplDriveCallerDR;

	/* Scalar functions */
	typedef std::map<std::string, const BasicScalarFunction *> SFType;
	SFType SF;

	bool ScalarFunction_int(void);
friend struct SFuncDR;

	/* Dynamic modules */
	bool moduleInitialized;
	bool ModuleLoad_int(void);
friend struct ModuleLoadDR;

	DataManager *pDM;

public:
	class Manip {
	public:
		Manip(void) {};
		virtual ~Manip(void) {};
	};

private:
	std::stack<const Manip *> manip;

public:
	void PopManip(void);
	void PushManip(const Manip *);
	const Manip *GetManip(void) const;
	bool bEmptyManip(void) const;

public:
	template <class T>
	class TplManip : public Manip {
	public:
		virtual ~TplManip(void);
		virtual T Get(const T& t) const = 0;
	};

private:
	/* not allowed */
	MBDynParser(const MBDynParser&);

public:
	MBDynParser(MathParser& MP, InputStream& streamIn,
			const char *initial_file);
	~MBDynParser(void);

	void SetDataManager(DataManager *pdm);
	DataManager *GetDataManager(void) const;

	// Vec*, Mat*x* moved here from parser.h
	/* vettore Vec3 */
	virtual Vec3 GetVec3(void);
	/* vettore Vec3 */
	virtual Vec3 GetVec3(const Vec3& vDef);
	/* matrice R mediante i due vettori */
	virtual Mat3x3 GetMatR2vec(void);
	/* matrice 3x3 simmetrica */
	virtual Mat3x3 GetMat3x3Sym(void);
	/* matrice 3x3 arbitraria */
	virtual Mat3x3 GetMat3x3(void);
	virtual Mat3x3 GetMat3x3(const Mat3x3& mDef);

	virtual Vec6 GetVec6(void);
	virtual Vec6 GetVec6(const Vec6& vDef);
	virtual Mat6x6 GetMat6x6(void);
	virtual Mat6x6 GetMat6x6(const Mat6x6& mDef);

public:
	virtual void GetMat6xN(Mat3xN& m1, Mat3xN& m2, integer iNumCols);

	/*
	 * Lettura di posizioni, vettori e matrici di rotazione
	 * relative ed assolute rispetto ad un riferimento
	 */
	Vec3 GetPosRel(const ReferenceFrame& rf);
	Vec3 GetPosRel(const ReferenceFrame& rf, const ReferenceFrame& other_rf, const Vec3& other_X);
	Vec3 GetPosAbs(const ReferenceFrame& rf);
	Vec3 GetVelRel(const ReferenceFrame& rf, const Vec3& x);
	Vec3 GetVelAbs(const ReferenceFrame& rf, const Vec3& x);
	Vec3 GetOmeRel(const ReferenceFrame& rf);
	Vec3 GetOmeAbs(const ReferenceFrame& rf);
	Vec3 GetVecRel(const ReferenceFrame& rf);
	Vec3 GetVecAbs(const ReferenceFrame& rf);
	Vec3 GetUnitVecRel(const ReferenceFrame& rf);
	Vec3 GetUnitVecAbs(const ReferenceFrame& rf);
	Mat3x3 GetMatRel(const ReferenceFrame& rf);
	Mat3x3 GetMatAbs(const ReferenceFrame& rf);
	Mat3x3 GetRotRel(const ReferenceFrame& rf);
	Mat3x3 GetRotRel(const ReferenceFrame& rf, const ReferenceFrame& other_rf, const Mat3x3& other_R);
	Mat3x3 GetRotAbs(const ReferenceFrame& rf);

	enum VecMatOpType {
		VM_NONE		= 0,

		VM_POSABS,
		VM_VELABS,
		VM_OMEABS,
		VM_VECABS,
		VM_UNITVECABS,

		VM_MATABS,
		VM_ROTABS,

		VM_MOD_REL	= 0x1000U,
		VM_MOD_OTHER	= 0x2000U,

		VM_POSREL	= VM_POSABS | VM_MOD_REL,
		VM_VELREL	= VM_VELABS | VM_MOD_REL,
		VM_OMEREL	= VM_OMEABS | VM_MOD_REL,
		VM_VECREL	= VM_VECABS | VM_MOD_REL,
		VM_UNITVECREL	= VM_UNITVECABS | VM_MOD_REL,

		VM_MATREL	= VM_MATABS | VM_MOD_REL,
		VM_ROTREL	= VM_ROTABS | VM_MOD_REL,

		VM_LAST
	};

	template <class T>
	class RefTplManip : public TplManip<T> {
	protected:
		MBDynParser& HP;
		const ReferenceFrame& rf;
		VecMatOpType type;

	public:
		RefTplManip(MBDynParser& HP, const ReferenceFrame& rf, VecMatOpType type);
		virtual ~RefTplManip(void);
	};

	class RefVec3Manip : public RefTplManip<Vec3> {
	public:
		RefVec3Manip(MBDynParser& HP, const ReferenceFrame& rf, VecMatOpType type);
		virtual ~RefVec3Manip(void);
		virtual Vec3 Get(const Vec3& v) const;
	};

	class VecRelManip : public RefVec3Manip {
	public:
		VecRelManip(MBDynParser& HP, const ReferenceFrame& rf);
		virtual ~VecRelManip(void);
	};

	class VecAbsManip : public RefVec3Manip {
	public:
		VecAbsManip(MBDynParser& HP, const ReferenceFrame& rf);
		virtual ~VecAbsManip(void);
	};

	virtual doublereal Get(const doublereal& d);
	virtual Vec3 Get(const Vec3& v);
	virtual Mat3x3 Get(const Mat3x3& m);
	virtual Vec6 Get(const Vec6& v);
	virtual Mat6x6 Get(const Mat6x6& m);

	void OutputFrames(std::ostream& out) const;

	HydraulicFluid* GetHydraulicFluid(void);

	const c81_data* GetC81Data(unsigned profile) const;

	ConstitutiveLaw1D* GetConstLaw1D(ConstLawType::Type& clt);
	ConstitutiveLaw3D* GetConstLaw3D(ConstLawType::Type& clt);
	ConstitutiveLaw6D* GetConstLaw6D(ConstLawType::Type& clt);
	const DriveCaller* GetDrive(unsigned uLabel) const;

	DriveCaller* GetDriveCaller(bool bDeferred = false);
	template <class T> const TplDriveCaller<T> *GetTplDrive(unsigned uLabel) const;
	template <class T> TplDriveCaller<T> *GetTplDriveCaller(void);

	const BasicScalarFunction* GetScalarFunction(void);
	const BasicScalarFunction* GetScalarFunction(const std::string &s);
	bool SetScalarFunction(const std::string &s, const BasicScalarFunction *sf);
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

