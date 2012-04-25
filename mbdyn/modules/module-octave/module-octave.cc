/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 AUTHOR: Reinhard Resch <reinhard.resch@accomp.it>
        Copyright (C) 2011(-2012) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <cmath>
#include <cfloat>
#include <string>
#include <cassert>
#include <limits>
#include <set>

// FIXME: there is a conflict between the MBDyn and octave real type
#define real mbdyn_real_type
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#include <matvec3.h>
#include <matvec6.h>
#include <dataman.h>
#include <drive.h>
#include <tpldrive.h>
#include <userelem.h>
#undef real

#include <octave/oct.h>
#include <octave/parse.h>
#include <octave/toplev.h>
#include <octave/octave.h>

#include "module-octave.h"
#include "octave_object.h"

//#define DEBUG

#ifdef DEBUG
#define TRACE(msg) ((void)(std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << ":" << msg << std::endl))
#else
#define TRACE(msg) ((void)0)
#endif

class OctaveInterface {
public:
	enum OctaveCallFlags_t {
		DEFAULT_CALL_FLAGS = 0x0,
		PASS_DATA_MANAGER = 0x1,
		UPDATE_GLOBAL_VARIABLES = 0x2
	};

private:
	explicit OctaveInterface(const DataManager* pDM);
	virtual ~OctaveInterface(void);

public:
	static OctaveInterface* CreateInterface(const DataManager* pDM);
	void Destroy();
	static OctaveInterface* GetInterface(void) { return pOctaveInterface; };
	void UpdateVariables(void);
	const DataManager* GetDataManager(void) const { return pDM; };
	octave_value GetDataManagerInterface(void) const { return pDataManager; };
	bool AddOctaveSearchPath(const std::string& path);
	static bool ConvertMBDynToOctave(const TypedValue& mbValue, octave_value& octValue);
	static bool ConvertMBDynToOctave(const Vec3& mbValue, octave_value& octValue);
	static bool ConvertMBDynToOctave(const Mat3x3& mbValue, octave_value& octValue);
	static bool ConvertOctaveToMBDyn(const Matrix& octValue, doublereal& mbValue);
	static bool ConvertOctaveToMBDyn(const Matrix& octValue, Vec3& mbValue);
	static bool ConvertOctaveToMBDyn(const Matrix& octValue, Vec6& mbValue);
	static bool ConvertOctaveToMBDyn(const Matrix& octValue, Mat3x3& mbValue);
	static bool ConvertOctaveToMBDyn(const Matrix& octValue, Mat6x6& mbValue);
	inline octave_value_list EvalFunction(const std::string& func, const octave_value_list& args, int nargout = 1, int flags = DEFAULT_CALL_FLAGS);
	inline octave_value_list EvalFunctionDerivative(const std::string& func, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
	inline octave_value_list EvalFunctionDerivative(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
	inline octave_value_list EvalFunction(const std::string& func, doublereal dVar, int nargout = 1, int flags = DEFAULT_CALL_FLAGS);
	inline doublereal EvalScalarFunction(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
	inline doublereal EvalScalarFunctionDerivative(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
	template<typename T>
	inline void EvalMatrixFunction(const std::string& func, T& res, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
	template<typename T>
	inline void EvalMatrixFunctionDerivative(const std::string& func, T& res, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
	static bool HaveADPackage(void) { return bHaveADPackage; };
	void AddEmbedFileName(const std::string& strFile);

private:
	template<typename T>
	inline static bool ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows);
	template<typename T>
	inline static bool ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows, int cols);
	bool LoadADPackage(void);

private:
	typedef std::set<std::string> EmbedFileNameSet_t;
	typedef EmbedFileNameSet_t::iterator EmbedFileNameIter_t;
	const DataManager* pDM;
	bool bFirstCall;
	octave_value pDataManager;
	EmbedFileNameSet_t strEmbedFileNames;
	static OctaveInterface* pOctaveInterface;
	static int iRefCount;
	static const std::string strADFunc;
	static bool bHaveADPackage;
};

class MBDynInterface : public octave_object {
public:
	// An instance of this class might be created directly from within octave
	// In this case pInterface must be defined
	explicit MBDynInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface());
	virtual ~MBDynInterface(void);
	virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;

protected:
	BEGIN_METHOD_TABLE_DECLARE()
		METHOD_DECLARE(GetVersion)
	END_METHOD_TABLE_DECLARE()

private:
	DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
	DECLARE_OCTAVE_ALLOCATOR

protected:
	OctaveInterface* GetInterface(void) const { return pInterface; };

private:
	OctaveInterface* const pInterface;
};

class StructNodeInterface : public MBDynInterface {
public:
	explicit StructNodeInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const StructNode* pNode = 0);
	virtual ~StructNodeInterface();
	virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;

protected:
	BEGIN_METHOD_TABLE_DECLARE()
		METHOD_DECLARE(iGetFirstIndex)
		METHOD_DECLARE(iGetFirstPositionIndex)
		METHOD_DECLARE(iGetFirstMomentumIndex)
		METHOD_DECLARE(GetLabel)
		METHOD_DECLARE(GetXCurr)
		METHOD_DECLARE(GetXPrev)
		METHOD_DECLARE(GetgCurr)
		METHOD_DECLARE(GetgRef)
		METHOD_DECLARE(GetgPCurr)
		METHOD_DECLARE(GetgPRef)
		METHOD_DECLARE(GetRCurr)
		METHOD_DECLARE(GetRPrev)
		METHOD_DECLARE(GetRRef)
		METHOD_DECLARE(GetVCurr)
		METHOD_DECLARE(GetVPrev)
		METHOD_DECLARE(GetWCurr)
		METHOD_DECLARE(GetWPrev)
		METHOD_DECLARE(GetXPPCurr)
		METHOD_DECLARE(GetXPPPrev)
		METHOD_DECLARE(GetWPCurr)
		METHOD_DECLARE(GetWPPrev)
	END_METHOD_TABLE_DECLARE()

private:
	DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
	DECLARE_OCTAVE_ALLOCATOR

private:
	octave_value GetVec3(const Vec3& v, const octave_value_list& args) const;
	octave_value GetMat3x3(const Mat3x3& m, const octave_value_list& args) const;

private:
	const StructNode *const pNode;
};

class DataManagerInterface : public MBDynInterface {
public:
	// An instance of this class might be created directly from within octave
	// In this case pInterface must be defined
	explicit DataManagerInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface());
	virtual ~DataManagerInterface();
	virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;

protected:
	BEGIN_METHOD_TABLE_DECLARE()
		METHOD_DECLARE(GetVariable)
		METHOD_DECLARE(GetStructNodePos)
		METHOD_DECLARE(GetStructNode)
	END_METHOD_TABLE_DECLARE()

private:
	DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
	DECLARE_OCTAVE_ALLOCATOR

private:
	const DataManager* GetDataManager(void);
	const Table& GetSymbolTable(void);
};

class OctaveDriveCaller : public DriveCaller {
public:
	explicit OctaveDriveCaller(const std::string& strFunc, OctaveInterface* pInterface, int iFlags);
	virtual ~OctaveDriveCaller(void);
 
	/* Copia */
	virtual DriveCaller* pCopy(void) const;
 
	/* Scrive il contributo del DriveCaller al file di restart */   
	virtual std::ostream& Restart(std::ostream& out) const;
 
	inline doublereal dGet(const doublereal& dVar) const;
	inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;
private:
	int iFlags;
	const std::string strFunc;
	OctaveInterface* pInterface;
};

template <class T>
class OctaveTplDriveCaller : public TplDriveCaller<T> {
public:
	OctaveTplDriveCaller(const std::string& strFunction, OctaveInterface* pInterface, int iFlags);
	~OctaveTplDriveCaller(void);
	virtual TplDriveCaller<T>* pCopy(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	virtual std::ostream& Restart_int(std::ostream& out) const;
	virtual inline T Get(const doublereal& dVar) const;
	virtual inline T Get(void) const;
	virtual inline bool bIsDifferentiable(void) const;
	virtual inline T GetP(void) const;
	virtual inline int getNDrives(void) const;

private:
	const std::string strFunction;
	OctaveInterface* const pInterface;
	const int iFlags;
};

class DerivativeDriveCaller : public DriveCaller {
public:
	explicit DerivativeDriveCaller(DriveCaller* pDriveCaller);
	virtual ~DerivativeDriveCaller(void);

	/* Copia */
	virtual DriveCaller* pCopy(void) const;

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const;

	inline doublereal dGet(const doublereal& dVar) const;
	inline doublereal dGet(void) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
	virtual doublereal dGetP(const doublereal& dVar) const;
	virtual inline doublereal dGetP(void) const;

private:
	DriveOwner pDriveCaller;
};

class OctaveScalarFunction : public DifferentiableScalarFunction {
public:
	OctaveScalarFunction(const std::string& strFunc, OctaveInterface* pInterface, int iFlags);
	virtual ~OctaveScalarFunction(void);
	virtual doublereal operator()(const doublereal x) const;
	virtual doublereal ComputeDiff(const doublereal t, const integer order = 1) const;

private:
	const std::string strFunc;
	OctaveInterface* const pInterface;
	const int iFlags;
};

class OctaveBaseDCR {
public:
	OctaveBaseDCR(void);
	virtual ~OctaveBaseDCR(void);
	void Read(const DataManager* pDM, MBDynParser& HP, bool bDefered = false) const;
	OctaveInterface* GetInterface(void) const { return pInterface; };
	const std::string& GetFunction(void) const { return strFunction; };
	int GetFlags(void) const { return iFlags; };

private:
	// FIXME: mutable needed for ScalarFunctionRead
	mutable OctaveInterface* pInterface;
	mutable std::string strFunction;
	mutable int iFlags;
};

class OctaveDCR : public DriveCallerRead, public OctaveBaseDCR {
public:
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

template <class T>
class OctaveTDCR : public TplDriveCallerRead<T>, public OctaveBaseDCR {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP);
};

class DerivativeDCR : public DriveCallerRead {
public:
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

class OctaveSFR : public ScalarFunctionRead, public OctaveBaseDCR {
public:
	const BasicScalarFunction *
	Read(DataManager* pDM, MBDynParser& HP) const;
};

class OctaveElement: virtual public Elem, public UserDefinedElem {
public:
	OctaveElement(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~OctaveElement(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
   	OctaveInterface* GetInterface(void) const { return dcr.GetInterface(); };
   	int GetFlags(void) const { return dcr.GetFlags(); };
   	const std::string& GetClass(void) const { return dcr.GetFunction(); };

private:
   	octave_value octObject;
   	octave_value mbdObject;
   	OctaveBaseDCR dcr;
   	enum {
   		JACOBIAN_NO = 0,
   		JACOBIAN_FULL = 1,
   		JACOBIAN_SPARSE = 2
   	} haveJacobian;
   	static const std::string strWorkSpaceDim;
   	static const std::string striGetNumDof;
   	static const std::string strAssRes;
   	static const std::string strAssJac;
   	static const std::string strIsMethod;
};

class OctaveElementInterface: public MBDynInterface {
public:
	explicit OctaveElementInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), OctaveElement* pElem = 0);
	virtual ~OctaveElementInterface(void);
	virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;

protected:
	BEGIN_METHOD_TABLE_DECLARE()
		METHOD_DECLARE(iGetFirstIndex)
	END_METHOD_TABLE_DECLARE()

private:
	DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
	DECLARE_OCTAVE_ALLOCATOR

private:
	OctaveElement* const pElem;
};

OctaveInterface* OctaveInterface::pOctaveInterface = 0;
int OctaveInterface::iRefCount = 0;
const std::string OctaveInterface::strADFunc("mbdyn_derivative");
bool OctaveInterface::bHaveADPackage = false;

OctaveInterface::OctaveInterface(const DataManager* pDM)
: pDM(pDM),
bFirstCall(true),
pDataManager(new DataManagerInterface(this))
{
	TRACE("constructor");

	ASSERT(pOctaveInterface == 0);

	pOctaveInterface = this;

	int argc = 1;
	if (silent_output) {
		argc++;
	}
	if (pedantic_output) {
		argc++;
	}
	string_vector argv(argc);

	argc = 0;
	argv(argc) = "octave";
	if (silent_output) {
		argv(++argc) = "-q";
	}
	if (pedantic_output) {
		argv(++argc) = "-V";
	}

	octave_main(argv.nelem(), argv.c_str_vec(), 1);

	LoadADPackage();
}

OctaveInterface::~OctaveInterface(void)
{
	TRACE("destructor");

	ASSERT(this == pOctaveInterface);

	do_octave_atexit();
	pOctaveInterface = 0;
}

bool
OctaveInterface::LoadADPackage(void)
{
	octave_value_list args;
	args.append(octave_value("load"));
	args.append(octave_value("ad"));

	feval("pkg", args, 0);

	bHaveADPackage = (error_state == 0);

	if (!bHaveADPackage) {
		silent_cerr("warning: octave package for automatic forward differentiation is not available" << std::endl);
		error_state = 0; // ignore error
	}

	return bHaveADPackage;
}

OctaveInterface *
OctaveInterface::CreateInterface(const DataManager* pDM)
{
	++iRefCount;

	if (pOctaveInterface) {
		return pOctaveInterface;
	}

	return new OctaveInterface(pDM);
}

void
OctaveInterface::Destroy(void)
{
	ASSERT(iRefCount >= 1);

	if (0 == --iRefCount) {
		delete this;
	}
}

void
OctaveInterface::UpdateVariables(void)
{
	typedef Table::VM::const_iterator iterator;

	const Table& symbolTable = pDM->GetMathParser().GetSymbolTable();

	for (iterator it = symbolTable.begin(); it != symbolTable.end(); ++it)
	{
		const std::string& mbName(it->first);
		const NamedValue*const namedValue = it->second;
		const TypedValue mbValue(namedValue->GetVal());
		octave_value octValue;

		if (!ConvertMBDynToOctave(mbValue, octValue)) {
			silent_cerr("octave error: data type \"" << mbValue.GetType() << "\" of variable \"" << mbName << "\": not handled in switch statement " << std::endl);
			ASSERT(0);
    			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		set_global_value(mbName, octValue);
	}
}

bool
OctaveInterface::AddOctaveSearchPath(const std::string& path)
{
	feval("addpath", octave_value(path));

	return error_state == 0;
}

bool
OctaveInterface::ConvertMBDynToOctave(const TypedValue& mbValue, octave_value& octValue)
{
	switch (mbValue.GetType()) {
	case TypedValue::VAR_BOOL:
		octValue = mbValue.GetBool();
		break;
	case TypedValue::VAR_INT:
		octValue = mbValue.GetInt();
		break;
	case TypedValue::VAR_REAL:
		octValue = mbValue.GetReal();
		break;
	case TypedValue::VAR_STRING:
		octValue = mbValue.GetString();
		break;
	default:
		return false;
	}

	return true;
}

bool
OctaveInterface::ConvertMBDynToOctave(const Vec3& mbValue, octave_value& octValue)
{
	ColumnVector V(3);

	for (int i = 0; i < 3; ++i) {
		V(i) = mbValue(i + 1);
	}

	octValue = V;

	return true;
}

bool
OctaveInterface::ConvertMBDynToOctave(const Mat3x3& mbValue, octave_value& octValue)
{
	Matrix M(3, 3);

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			M(i, j) = mbValue(i + 1, j + 1);
		}
	}

	octValue = M;

	return true;
}

template<typename T>
bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows)
{
	if (octValue.rows() != rows || octValue.columns() != 1) {
		silent_cerr("octave error: invalid matrix size " << octValue.rows() << "x" << octValue.columns() << " expected " << rows << "x1" << std::endl);
		return false;
	}

	for (int i = 0; i < rows; ++i) {
		mbValue(i + 1) = octValue(i, 0);
	}

	return true;
}

template<typename T>
bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows, int cols)
{
	if (octValue.rows() != rows || octValue.columns() != cols) {
		silent_cerr("octave error: invalid matrix size " << octValue.rows() << "x" << octValue.columns() << " expected " << rows << "x" << cols << std::endl);
		return false;
	}

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			mbValue(i + 1, j + 1) = octValue(i, j);
		}
	}

	return true;
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, doublereal& mbValue)
{
	if (octValue.rows() != 1 || octValue.columns() != 1) {
		silent_cerr("octave error: invalid matrix size " << octValue.rows() << "x" << octValue.columns() << " expected 1x1" << std::endl);
		return false;
	}

	mbValue = octValue(0, 0);

	return true;
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, Vec3& mbValue)
{
	return ConvertOctaveToMBDyn(octValue, mbValue, 3);
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, Vec6& mbValue)
{
	return ConvertOctaveToMBDyn(octValue, mbValue, 6);
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, Mat3x3& mbValue)
{
	return ConvertOctaveToMBDyn(octValue, mbValue, 3, 3);
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const Matrix& octValue, Mat6x6& mbValue)
{
	return ConvertOctaveToMBDyn(octValue, mbValue, 6, 6);
}

octave_value_list
OctaveInterface::EvalFunction(const std::string& func, doublereal dVar, int nargout, int flags)
{
	octave_value_list args;

	args.append(octave_value(dVar));

	if (flags & PASS_DATA_MANAGER) {
		args.append(pDataManager);
	}

	return EvalFunction(func, args, nargout, flags);
}

octave_value_list
OctaveInterface::EvalFunction(const std::string& func, const octave_value_list& args, int nargout, int flags)
{
	if ((flags & UPDATE_GLOBAL_VARIABLES) || bFirstCall) {
		UpdateVariables();
	}

	if (bFirstCall) {
		for (EmbedFileNameIter_t pFile = strEmbedFileNames.begin();
			pFile != strEmbedFileNames.end(); ++pFile )
		{
			feval("source", octave_value(*pFile));

			if (error_state) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	bFirstCall = false;

	octave_value_list ans = feval(func, args, nargout);

	if (error_state) {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (ans.length() < nargout) {
		silent_cerr("octave error: function \"" << func << "\" returned less than " << nargout << " values" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	for (int i = 0; i < nargout; ++i) {
		if (!ans(i).is_defined()) {
			silent_cerr("octave error: result " << i + 1 << "of function \"" << func << "\" is undefined" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	return ans;
}

octave_value_list
OctaveInterface::EvalFunctionDerivative(const std::string& func, const octave_value_list& args, int flags)
{
	TRACE("func=" << func);
	TRACE("flags=" << flags);

	octave_value_list derivArgs;
	derivArgs.append(func);
	derivArgs.append(args);

	return EvalFunction(strADFunc, derivArgs, 2, flags);
}

octave_value_list
OctaveInterface::EvalFunctionDerivative(const std::string& func, doublereal dVar, int flags)
{
	TRACE("func=" << func);
	TRACE("flags=" << flags);

	octave_value_list args = octave_value(dVar);

	if (flags & PASS_DATA_MANAGER) {
		args.append(pDataManager);
	}

	return EvalFunctionDerivative(func, args, flags);
}

inline doublereal
OctaveInterface::EvalScalarFunction(const std::string& func, doublereal dVar, int flags)
{
	octave_value_list ans = EvalFunction(func, dVar, 1, flags);

	if (!ans(0).is_real_scalar()) {
		silent_cerr("octave error: result of function \"" << func << "\" is not a scalar value" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return ans(0).scalar_value();
}

doublereal
OctaveInterface::EvalScalarFunctionDerivative(const std::string& func, doublereal dVar, int flags)
{
	TRACE("func=" << func);
	TRACE("flags=" << flags);

	octave_value_list ans = EvalFunctionDerivative(func, dVar, flags);

	ASSERT(ans.length() == 2);

	if (!ans(1).is_real_scalar()) {
		silent_cerr("octave error: derivative of function \"" << func << "\" is not a scalar value" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return ans(1).scalar_value();
}

template <typename T>
inline void
OctaveInterface::EvalMatrixFunction(const std::string& func, T& res, doublereal dVar, int flags)
{
	octave_value_list ans = EvalFunction(func, dVar, 1, flags);

	if (!ans(0).is_real_matrix()) {
		silent_cerr("octave error: result of function \"" << func << "\" is not a matrix value" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Matrix A(ans(0).matrix_value());

	if (!ConvertOctaveToMBDyn(A, res)) {
		silent_cerr("octave error: conversion of octave matrix failed!" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

template <typename T>
inline void
OctaveInterface::EvalMatrixFunctionDerivative(const std::string& func, T& res, doublereal dVar, int flags)
{
	octave_value_list ans = EvalFunctionDerivative(func, dVar, flags);

	ASSERT(ans.length() == 2);

	if (!ans(0).is_real_matrix()) {
		silent_cerr("octave error: result of function \"" << func << "\" is not a matrix value" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (!ans(1).is_real_matrix()) {
		silent_cerr("octave error: result of derivative of function \"" << func << "\" is not a matrix value" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	Matrix A(ans(1).matrix_value());

	if (!ConvertOctaveToMBDyn(A, res)) {
		silent_cerr("octave error: conversion of octave matrix failed!" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

void
OctaveInterface::AddEmbedFileName(const std::string& strFile)
{
	strEmbedFileNames.insert(strFile);
}

MBDynInterface::MBDynInterface(OctaveInterface* pInterface)
: pInterface(pInterface)
{
	TRACE("constructor");
}

MBDynInterface::~MBDynInterface(void)
{
	TRACE("destructor");
}

void
MBDynInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
	os << "MBDyn " << VERSION << " interface" << std::endl;
}

METHOD_DEFINE(MBDynInterface, GetVersion, args, nargout)
{
	octave_value octValue;

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	octValue = VERSION;

	return octValue;
}

BEGIN_METHOD_TABLE(MBDynInterface, octave_object)
	METHOD_DISPATCH(MBDynInterface, GetVersion)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(MBDynInterface);
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(MBDynInterface, "MBDyn", "MBDyn");

StructNodeInterface::StructNodeInterface(OctaveInterface* pInterface, const StructNode* pNode)
: MBDynInterface(pInterface), pNode(pNode)
{
	TRACE("constructor");
}

StructNodeInterface::~StructNodeInterface(void)
{
	TRACE("destructor");
}

void
StructNodeInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
	os << "StructNode " << pNode->GetLabel() << ": " << pNode->GetName() << std::endl;
}

octave_value
StructNodeInterface::GetVec3(const Vec3& V, const octave_value_list& args) const
{
	octave_value octValue;

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	if (!GetInterface()->ConvertMBDynToOctave(V, octValue)) {
		error("could not convert data!\n");
	}

	return octValue;
}

octave_value
StructNodeInterface::GetMat3x3(const Mat3x3& M, const octave_value_list& args) const
{
	octave_value octValue;

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	if (!GetInterface()->ConvertMBDynToOctave(M, octValue)) {
		error("could not convert data!\n");
	}

	return octValue;
}

METHOD_DEFINE(StructNodeInterface, iGetFirstIndex, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octave_value();
	}

	return octave_value(pNode->iGetFirstIndex());
}

METHOD_DEFINE(StructNodeInterface, iGetFirstPositionIndex, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octave_value();
	}

	return octave_value(pNode->iGetFirstPositionIndex());
}

METHOD_DEFINE(StructNodeInterface, iGetFirstMomentumIndex, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octave_value();
	}

	return octave_value(pNode->iGetFirstMomentumIndex());
}

METHOD_DEFINE(StructNodeInterface, GetLabel, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octave_value();
	}

	return octave_value(pNode->GetLabel());
}

METHOD_DEFINE(StructNodeInterface, GetXCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetXCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetXPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetXPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetgCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgRef, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetgRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgPCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetgPCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgPRef, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetgPRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetMat3x3(pNode->GetRCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetMat3x3(pNode->GetRPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRRef, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetMat3x3(pNode->GetRRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetVCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetVCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetVPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetVPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetWCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	return GetVec3(pNode->GetWPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetXPPCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (!pNode->bComputeAccelerations()) {
		error("accelerations are not available for node %u!\n", pNode->GetLabel());
		return octave_value();
	}

	return GetVec3(pNode->GetXPPCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetXPPPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (!pNode->bComputeAccelerations()) {
		error("accelerations are not available for node %u!\n", pNode->GetLabel());
		return octave_value();
	}

	return GetVec3(pNode->GetXPPPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPCurr, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (!pNode->bComputeAccelerations()) {
		error("accelerations are not available for node %u!\n", pNode->GetLabel());
		return octave_value();
	}

	return GetVec3(pNode->GetWPCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPPrev, args, nargout)
{
	if (pNode == 0) {
		error("node data is not available!\n");
		return octave_value();
	}

	if (!pNode->bComputeAccelerations()) {
		error("accelerations are not available for node %u!\n", pNode->GetLabel());
		return octave_value();
	}

	return GetVec3(pNode->GetWPPrev(), args);
}

BEGIN_METHOD_TABLE(StructNodeInterface, MBDynInterface)
	METHOD_DISPATCH(StructNodeInterface, iGetFirstIndex)
	METHOD_DISPATCH(StructNodeInterface, iGetFirstPositionIndex)
	METHOD_DISPATCH(StructNodeInterface, iGetFirstMomentumIndex)
	METHOD_DISPATCH(StructNodeInterface, GetLabel)
	METHOD_DISPATCH(StructNodeInterface, GetXCurr)
	METHOD_DISPATCH(StructNodeInterface, GetXPrev)
	METHOD_DISPATCH(StructNodeInterface, GetgCurr)
	METHOD_DISPATCH(StructNodeInterface, GetgRef)
	METHOD_DISPATCH(StructNodeInterface, GetgPCurr)
	METHOD_DISPATCH(StructNodeInterface, GetgPRef)
	METHOD_DISPATCH(StructNodeInterface, GetRCurr)
	METHOD_DISPATCH(StructNodeInterface, GetRPrev)
	METHOD_DISPATCH(StructNodeInterface, GetRRef)
	METHOD_DISPATCH(StructNodeInterface, GetVCurr)
	METHOD_DISPATCH(StructNodeInterface, GetVPrev)
	METHOD_DISPATCH(StructNodeInterface, GetWCurr)
	METHOD_DISPATCH(StructNodeInterface, GetWPrev)
	METHOD_DISPATCH(StructNodeInterface, GetXPPCurr)
	METHOD_DISPATCH(StructNodeInterface, GetXPPPrev)
	METHOD_DISPATCH(StructNodeInterface, GetWPCurr)
	METHOD_DISPATCH(StructNodeInterface, GetWPPrev)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(StructNodeInterface);
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(StructNodeInterface, "StructNode", "StructNode");

DataManagerInterface::DataManagerInterface(OctaveInterface* pInterface)
	:MBDynInterface(pInterface)
{
	TRACE("constructor");

	ASSERT(pInterface != 0);
}

DataManagerInterface::~DataManagerInterface()
{
	TRACE("destructor");
}

void DataManagerInterface::print(std::ostream& os, bool pr_as_read_syntax)const
{
	os << "a MBDyn DataManager interface" << std::endl;
}

METHOD_DEFINE(DataManagerInterface, GetVariable, args, nargout)
{
	octave_value octValue;

	if (args.length() < 1) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	const std::string varName = args(0).string_value();

	if (error_state) {
		error("invalid argument type for argument number one!\n");
		return octValue;
	}

	const NamedValue* const pVar = GetSymbolTable().Get(varName.c_str());

	if (pVar == 0) {
		error("Variable \"%s\" not found in symbol table!\n", varName.c_str());
		return octValue;
	}

	const TypedValue mbValue = pVar->GetVal();

	if (!GetInterface()->ConvertMBDynToOctave(mbValue, octValue)) {
		error("Failed to convert MBDyn variable \"%s\" value to octave value!\n", varName.c_str());
	}

	return octValue;
}

METHOD_DEFINE(DataManagerInterface, GetStructNodePos, args, nargout)
{
	octave_value octValue;

	if (args.length() < 1) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	integer iNodeLabel = args(0).int32_scalar_value();

	if (error_state) {
		error("invalid argument type for argument number one!\n");
		return octValue;
	}

	Node* pNode = GetDataManager()->pFindNode(Node::STRUCTURAL, iNodeLabel);

	if (pNode == 0) {
		error("could not find node %ld!\n", (long)iNodeLabel);
		return octValue;
	}

	StructNode* pStructNode = dynamic_cast<StructNode*>(pNode);

	if (pStructNode == 0) {
		error("invalid type of node %ld\n", (long)iNodeLabel);
		return octValue;
	}

	const Vec3& X = pStructNode->GetXCurr();

	ColumnVector octX(3);

	for (int i = 0; i < 3; ++i) {
		octX(i) = X(i + 1);
	}

	octValue = octX;

	return octValue;
}

METHOD_DEFINE(DataManagerInterface, GetStructNode, args, nargout)
{
	octave_value octValue;

	if (args.length() < 1) {
		error("invalid number of arguments!\n");
		return octValue;
	}

	integer iNodeLabel = args(0).int32_scalar_value();

	if (error_state) {
		error("invalid argument type for argument number one (node id)!\n");
		return octValue;
	}

	Node* pNode = GetDataManager()->pFindNode(Node::STRUCTURAL, iNodeLabel);

	if (pNode == 0) {
		error("could not find node %ld!\n", (long)iNodeLabel);
		return octValue;
	}

	StructNode* pStructNode = dynamic_cast<StructNode*>(pNode);

	if (pStructNode == 0) {
		error("invalid type of node %ld\n", (long)iNodeLabel);
		return octValue;
	}

	octValue = new StructNodeInterface(GetInterface(), pStructNode);

	return octValue;
}

const DataManager* DataManagerInterface::GetDataManager()
{
	ASSERT(GetInterface() != 0);

	return GetInterface()->GetDataManager();
}

const Table& DataManagerInterface::GetSymbolTable()
{
	return GetDataManager()->GetMathParser().GetSymbolTable();
}

BEGIN_METHOD_TABLE(DataManagerInterface, MBDynInterface)
	METHOD_DISPATCH(DataManagerInterface, GetVariable)
	METHOD_DISPATCH(DataManagerInterface, GetStructNodePos)
	METHOD_DISPATCH(DataManagerInterface, GetStructNode)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(DataManagerInterface);
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(DataManagerInterface, "DataManager", "DataManager");

OctaveDriveCaller::OctaveDriveCaller(const std::string& strFunc, OctaveInterface* pInterface, int iFlags)
: DriveCaller(0),
  iFlags(iFlags),
  strFunc(strFunc),
  pInterface(pInterface)
{
	TRACE("constructor");
	TRACE("strFunc=" << strFunc);
	TRACE("iFlags=" <<  iFlags);
}

OctaveDriveCaller::~OctaveDriveCaller(void)
{
	TRACE("destructor");
}

DriveCaller *
OctaveDriveCaller::pCopy(void) const
{
	return new OctaveDriveCaller(strFunc, pInterface, iFlags);
}

std::ostream&
OctaveDriveCaller::Restart(std::ostream& out) const
{
	return out << "octave, \"" << strFunc << "\"";
}

inline doublereal 
OctaveDriveCaller::dGet(const doublereal& dVar) const 
{
    return pInterface->EvalScalarFunction(strFunc, dVar, iFlags);
}

inline doublereal
OctaveDriveCaller::dGet(void) const
{
	return pInterface->EvalScalarFunction(strFunc, pInterface->GetDataManager()->dGetTime(), iFlags);
}

inline bool
OctaveDriveCaller::bIsDifferentiable(void) const
{
	return pInterface->HaveADPackage();
}

inline doublereal 
OctaveDriveCaller::dGetP(const doublereal& dVar) const
{
	if (!bIsDifferentiable()) {
		return 0.;
	}

	return pInterface->EvalScalarFunctionDerivative(strFunc, dVar, iFlags);
}

inline doublereal 
OctaveDriveCaller::dGetP(void) const
{
	if (!bIsDifferentiable()) {
		return 0.;
	}

	const doublereal t = pInterface->GetDataManager()->dGetTime();

	return pInterface->EvalScalarFunctionDerivative(strFunc, t, iFlags);
}

template <class T>
OctaveTplDriveCaller<T>::OctaveTplDriveCaller(const std::string& strFunction, OctaveInterface* pInterface, int iFlags)
	:strFunction(strFunction), pInterface(pInterface), iFlags(iFlags)
{
	TRACE("constructor");
};

template <class T>
OctaveTplDriveCaller<T>::~OctaveTplDriveCaller(void)
{
	TRACE("destructor");
};


template <class T>
TplDriveCaller<T>* OctaveTplDriveCaller<T>::pCopy(void) const
{
	return new OctaveTplDriveCaller(strFunction, pInterface, iFlags);
};

template <class T>
std::ostream& OctaveTplDriveCaller<T>::Restart(std::ostream& out) const
{
	return out;
};

template <class T>
std::ostream& OctaveTplDriveCaller<T>::Restart_int(std::ostream& out) const
{
	return out;
};

template <class T>
T OctaveTplDriveCaller<T>::Get(const doublereal& dVar) const
{
	T X;
	pInterface->EvalMatrixFunction(strFunction, X, dVar, iFlags);
	return X;
};

template <class T>
T OctaveTplDriveCaller<T>::Get(void) const
{
	doublereal t = pInterface->GetDataManager()->dGetTime();
	T X;
	pInterface->EvalMatrixFunction(strFunction, X, t, iFlags);
	return X;
};

template <class T>
bool OctaveTplDriveCaller<T>::bIsDifferentiable(void) const
{
	return pInterface->HaveADPackage();
};

template <class T>
T OctaveTplDriveCaller<T>::GetP(void) const
{
	const doublereal t = pInterface->GetDataManager()->dGetTime();
	T XP;
	pInterface->EvalMatrixFunctionDerivative(strFunction, XP, t, iFlags);
	return XP;
};

template <class T>
int OctaveTplDriveCaller<T>::getNDrives(void) const
{
	return 0;
};

template <class T>
TplDriveCaller<T> *
OctaveTDCR<T>::Read(const DataManager* pDM, MBDynParser& HP)
{
	OctaveBaseDCR::Read(pDM, HP);

	return new OctaveTplDriveCaller<T>(GetFunction(), GetInterface(), GetFlags());
};

DerivativeDriveCaller::DerivativeDriveCaller(DriveCaller* pDriveCaller)
	:DriveCaller(0), pDriveCaller(pDriveCaller)
{
	TRACE("constructor");
	ASSERT(pDriveCaller->bIsDifferentiable());
}

DerivativeDriveCaller::~DerivativeDriveCaller(void)
{
	TRACE("destructor");
}

DriveCaller* DerivativeDriveCaller::pCopy(void) const
{
	return new DerivativeDriveCaller(pDriveCaller.pGetDriveCaller()->pCopy());
}

std::ostream& DerivativeDriveCaller::Restart(std::ostream& out) const
{
	return out;
}

doublereal DerivativeDriveCaller::dGet(const doublereal& dVar) const
{
	return pDriveCaller.dGetP(dVar);
}

doublereal DerivativeDriveCaller::dGet(void) const
{
	return pDriveCaller.dGetP();
}

bool DerivativeDriveCaller::bIsDifferentiable(void) const
{
	return false;
}

doublereal DerivativeDriveCaller::dGetP(const doublereal& dVar) const
{
	return 0.;
}

doublereal DerivativeDriveCaller::dGetP(void) const
{
	return 0.;
}

OctaveScalarFunction::OctaveScalarFunction(const std::string& strFunc, OctaveInterface* pInterface, int iFlags)
	:strFunc(strFunc), pInterface(pInterface), iFlags(iFlags)
{

}

OctaveScalarFunction::~OctaveScalarFunction(void)
{

}

doublereal OctaveScalarFunction::operator()(const doublereal x) const
{
	return pInterface->EvalScalarFunction(strFunc, x, iFlags);
}

doublereal OctaveScalarFunction::ComputeDiff(const doublereal t, const integer order) const
{
	if (order != 1) {
		silent_cerr("octave scalar function \"" << strFunc << "\" derivative of order " << order << " not supported" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pInterface->EvalScalarFunctionDerivative(strFunc, t, iFlags);
}

OctaveBaseDCR::OctaveBaseDCR()
	:pInterface(0), iFlags(OctaveInterface::DEFAULT_CALL_FLAGS)
{
	TRACE("constructor");
	TRACE("iFlags=" <<  iFlags);
}

OctaveBaseDCR::~OctaveBaseDCR()
{
	TRACE("destructor");

	if (pInterface) {
		pInterface->Destroy();
	}
}

void OctaveBaseDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDefered)const
{
	if (pInterface == 0) {
	   pInterface = OctaveInterface::CreateInterface(pDM);
	}

	strFunction = HP.GetStringWithDelims(HighParser::DOUBLEQUOTE);

	iFlags = OctaveInterface::DEFAULT_CALL_FLAGS;

	if (HP.IsKeyWord("update" "global" "variables") && HP.GetYesNoOrBool()) {
		iFlags |= OctaveInterface::UPDATE_GLOBAL_VARIABLES;
	}

	if (HP.IsKeyWord("pass" "data" "manager") && HP.GetYesNoOrBool()) {
		iFlags |= OctaveInterface::PASS_DATA_MANAGER;
	}

	if (HP.IsKeyWord("embed" "octave") && HP.GetYesNoOrBool()) {
		IncludeParser::ErrOut sFileInfo = HP.GetLineData();
		pInterface->AddEmbedFileName(sFileInfo.sFileName);
	}

	if (HP.IsKeyWord("octave" "search" "path")) {
		while (HP.IsStringWithDelims(HighParser::DOUBLEQUOTE)) {
			std::string path = HP.GetFileName(HighParser::DOUBLEQUOTE);
			pedantic_cout("octave: adding \"" << path << "\" to octave search path at line " << HP.GetLineData() << std::endl);
			if (!pInterface->AddOctaveSearchPath(path)) {
				silent_cerr("octave error: addpath(\"" << path << "\") failed at line " << HP.GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
		}
	}

	TRACE("strFunction=" << strFunction);
	TRACE("iFlags=" << iFlags);
}

DriveCaller *
OctaveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	OctaveBaseDCR::Read(pDM, HP, bDeferred);

	return new OctaveDriveCaller(GetFunction(), GetInterface(), GetFlags());
};

DriveCaller *
DerivativeDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	DriveCaller* pDriveCaller = HP.GetDriveCaller(bDeferred);

	if (!pDriveCaller->bIsDifferentiable()) {
		silent_cerr("DriveCaller " << pDriveCaller->GetLabel() << " " << pDriveCaller->GetName() << " is not differentiable at line " << HP.GetLineData() << std::endl);

		delete pDriveCaller;

		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return new DerivativeDriveCaller(pDriveCaller);
};


const BasicScalarFunction *
OctaveSFR::Read(DataManager* pDM, MBDynParser& HP) const
{
	OctaveBaseDCR::Read(pDM, HP);
	return new OctaveScalarFunction(GetFunction(), GetInterface(), GetFlags());
};


const std::string OctaveElement::strWorkSpaceDim("WorkSpaceDim");
const std::string OctaveElement::striGetNumDof("iGetNumDof");
const std::string OctaveElement::strAssRes("AssRes");
const std::string OctaveElement::strAssJac("AssJac");
const std::string OctaveElement::strIsMethod("ismethod");

OctaveElement::OctaveElement(
	unsigned uLabel, const DofOwner *pDO,
	DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
  UserDefinedElem(uLabel, pDO),
  haveJacobian(JACOBIAN_NO)
{
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout("Module: 	octave\n"
			<< std::endl);

		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}

	dcr.Read(pDM, HP);

	mbdObject = new OctaveElementInterface(GetInterface(), this);

	octave_value_list args;

	while (HP.IsArg()) {
		const TypedValue vDef;
		const TypedValue varValue(HP.GetValue(vDef));
		octave_value octValue;
		GetInterface()->ConvertMBDynToOctave(varValue, octValue);
		args.append(octValue);
	}

	args.append(GetInterface()->GetDataManagerInterface());
	args.append(mbdObject);

	const int flags = OctaveInterface::UPDATE_GLOBAL_VARIABLES | OctaveInterface::PASS_DATA_MANAGER;

	octave_value_list ans = GetInterface()->EvalFunction(GetClass(), args, 1, flags);

	ASSERT(ans.length() == 1);

	octObject = ans(0);

	if (!octObject.is_object()) {
		silent_cerr("octave(" << GetLabel() << "): result of constructor is not an object at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	args.resize(0);
	args.append(octObject);
	args.append(octave_value(strAssJac));

	ans = GetInterface()->EvalFunction(strIsMethod, args, 1);

	ASSERT(ans.length() == 1);

	haveJacobian = ans(0).bool_value(true) ? JACOBIAN_FULL : JACOBIAN_NO;

	if (error_state) {
		silent_cerr("octave(" << GetLabel() << "): unexpected error in function " << strIsMethod << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

OctaveElement::~OctaveElement(void)
{
	GetInterface()->Destroy();
}

void
OctaveElement::Output(OutputHandler& OH) const
{
	// should do something useful
	NO_OP;
}

void
OctaveElement::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	octave_value_list args(octObject);
	args.append(GetInterface()->GetDataManagerInterface());
	args.append(mbdObject);

	octave_value_list ans = GetInterface()->EvalFunction(strWorkSpaceDim, args, 2, GetFlags());

	ASSERT(ans.length() == 2);

	*piNumRows = ans(0).int32_scalar_value();
	*piNumCols = ans(1).int32_scalar_value();

	if (error_state) {
		silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strWorkSpaceDim << " returned a incorrect data type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

VariableSubMatrixHandler&
OctaveElement::AssJac(VariableSubMatrixHandler& WorkMatVar,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	if (haveJacobian == JACOBIAN_NO) {
		WorkMatVar.SetNullMatrix();
		return WorkMatVar;
	}

	octave_value_list args(octObject);

	args.append(octave_value(dCoef));
	args.append(GetInterface()->GetDataManagerInterface());
	args.append(mbdObject);

	octave_value_list ans = GetInterface()->EvalFunction(strAssJac, args, 3, GetFlags());

	ASSERT(ans.length() == 3);

	Matrix Jac = ans(0).matrix_value();
	int32NDArray ridx = ans(1).int32_array_value();
	int32NDArray cidx = ans(2).int32_array_value();

	if (error_state) {
		silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssJac << " returned an incorrect data type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNumRows, iNumCols;

	WorkSpaceDim(&iNumRows, &iNumCols);

	FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

	WorkMat.ResizeReset(iNumRows, iNumCols);

	if (Jac.rows() != iNumRows
		|| Jac.columns() != iNumCols
		|| ridx.length() != iNumRows
		|| cidx.length() != iNumCols)
	{
		silent_cerr("octave(" << GetLabel() << "):"
			<< " rows(Jac)=" << Jac.rows()
			<< " columns(Jac)= " << Jac.columns()
			<< " length(ridx)= " << ridx.length()
			<< " length(cidx)=" << cidx.length()
			<< " are not consistent with iNumRows=" << iNumRows
			<< " and iNumCols=" << iNumCols << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNumDof = GetInterface()->GetDataManager()->iGetNumDofs();

	for (int i = 0; i < iNumRows; ++i) {
		if (int(ridx(i)) <= 0 || int(ridx(i)) > iNumDof) {
			silent_cerr("octave(" << GetLabel() << "):"
				<< " function " << GetClass() << "." << strAssJac
				<< ": row index " << ridx(i)
				<< " out of range [1:" << iNumDof << "]" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		WorkMat.PutRowIndex(i + 1, ridx(i));
	}

	for (int j = 0; j < iNumCols; ++j) {
		if (int(cidx(j)) <= 0 || int(cidx(j)) > iNumDof) {
			silent_cerr("octave(" << GetLabel() << "):"
				<< " function " << GetClass() << "." << strAssJac
				<< ": column index " << cidx(j)
				<< " out of range [1:" << iNumDof << "]" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		WorkMat.PutColIndex(j + 1, cidx(j));
	}

	for (int i = 0; i < iNumRows; ++i) {
		for (int j = 0; j < iNumCols; ++j) {
			WorkMat.PutCoef(i + 1, j + 1, Jac(i, j));
		}
	}

	return WorkMatVar;
}

SubVectorHandler&
OctaveElement::AssRes(SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{
	int iNumRows, iNumCols;

	WorkSpaceDim(&iNumRows, &iNumCols);

	WorkVec.ResizeReset(iNumRows);

	octave_value_list args(octObject);
	args.append(octave_value(dCoef));
	args.append(GetInterface()->GetDataManagerInterface());
	args.append(mbdObject);

	octave_value_list ans = GetInterface()->EvalFunction(strAssRes, args, 2, GetFlags());

	ASSERT(ans.length() == 2);

	ColumnVector f = ans(0).column_vector_value();
	int32NDArray ridx = ans(1).int32_array_value();

	if (error_state) {
		silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << " returned an incorrect data type" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (ridx.length() != iNumRows || f.length() != iNumRows) {
		silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": length(f)=" << f.length() << " length(ridx)=" << ridx.length() << " is not equal to iNumRows=" << iNumRows << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	integer iNumDof = GetInterface()->GetDataManager()->iGetNumDofs();

	for (octave_idx_type i = 0; i < f.length(); ++i) {
		if (int(ridx(i)) <= 0 || int(ridx(i)) > iNumDof) {
			silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": row index " << ridx(i) << " out of range [1:" << iNumDof << "]" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		WorkVec.PutItem(i + 1, ridx(i), f(i));
	}

	return WorkVec;
}

unsigned int
OctaveElement::iGetNumPrivData(void) const
{
	return 0;
}

int
OctaveElement::iGetNumConnectedNodes(void) const
{
	return 0;
}

void
OctaveElement::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	connectedNodes.resize(0);
}

void
OctaveElement::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
OctaveElement::Restart(std::ostream& out) const
{
	return out << "# OctaveElement: not implemented" << std::endl;
}

unsigned int
OctaveElement::iGetInitialNumDof(void) const
{
	return 0;
}

void
OctaveElement::InitialWorkSpaceDim(
	integer* piNumRows,
	integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
OctaveElement::InitialAssJac(
	VariableSubMatrixHandler& WorkMat,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkMat.SetNullMatrix();

	return WorkMat;
}

SubVectorHandler&
OctaveElement::InitialAssRes(
	SubVectorHandler& WorkVec,
	const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);

	WorkVec.ResizeReset(0);

	return WorkVec;
}

OctaveElementInterface::OctaveElementInterface(OctaveInterface* pInterface, OctaveElement* pElem)
	:MBDynInterface(pInterface),
	 pElem(pElem)
{

}

OctaveElementInterface::~OctaveElementInterface()
{

}

void
OctaveElementInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
	os << "MBDynElement(" << (pElem ? pElem->GetLabel() : (unsigned)(-1)) << ")" << std::endl;
}

METHOD_DEFINE(OctaveElementInterface, iGetFirstIndex, args, nargout)
{
	if (args.length() != 0) {
		error("invalid number of arguments!\n");
		return octave_value();
	}

	return octave_value(pElem->iGetFirstIndex());
}

BEGIN_METHOD_TABLE(OctaveElementInterface, MBDynInterface)
	METHOD_DISPATCH(OctaveElementInterface, iGetFirstIndex)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(OctaveElementInterface);
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(OctaveElementInterface, "MBDynElement", "MBDynElement");


bool
mbdyn_octave_set(void)
{
	DriveCallerRead	*rf = new OctaveDCR;
	if (!SetDriveData("octave", rf)) {
		delete rf;
		return false;
	}

	OctaveTDCR<doublereal>* const rf1D = new OctaveTDCR<doublereal>;
	if (!SetDC1D("octave", rf1D)) {
		delete rf1D;
		return false;
	}

	OctaveTDCR<Vec3>* const rf3D = new OctaveTDCR<Vec3>;
	if (!SetDC3D("octave", rf3D)) {
		delete rf3D;
		return false;
	}

	OctaveTDCR<Vec6>* const rf6D = new OctaveTDCR<Vec6>;
	if (!SetDC6D("octave", rf6D)) {
		delete rf6D;
		return false;
	}

	OctaveTDCR<Mat3x3>* const rf3x3D = new OctaveTDCR<Mat3x3>;
	if (!SetDC3x3D("octave", rf3x3D)) {
		delete rf3x3D;
		return false;
	}

	OctaveTDCR<Mat6x6>* const rf6x6D = new OctaveTDCR<Mat6x6>;
	if (!SetDC6x6D("octave", rf6x6D)) {
		delete rf6x6D;
		return false;
	}

	rf = new DerivativeDCR;
	if (!SetDriveData("derivative", rf)) {
		delete rf;
		return false;
	}

	ScalarFunctionRead* rsf = new OctaveSFR;
	if (!SetSF("octave", rsf)) {
		delete rsf;
		return false;
	}

	UserDefinedElemRead *urf = new UDERead<OctaveElement>;
	if (!SetUDE("octave", urf)) {
		delete urf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!mbdyn_octave_set()) {
		silent_cerr("octave: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // ! STATIC_MODULES
