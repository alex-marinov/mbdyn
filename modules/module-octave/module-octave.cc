/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2011(-2021) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#include <iostream>

#ifdef USE_OCTAVE

#if OCTAVE_MAJOR_VERSION < 5
// FIXME: there is a conflict between the MBDyn and octave real type
#define real mbdyn_real_type
#endif

#include <mbdefs.h>
#include <matvec3.h>
#include <matvec6.h>
#include <dataman.h>
#include <drive.h>
#include <tpldrive.h>
#include <userelem.h>
#include <constltp.h>

#if OCTAVE_MAJOR_VERSION < 5
#undef real
#endif

#include <cmath>
#include <cfloat>
#include <string>
#include <cassert>
#include <limits>
#include <map>
#include <sstream>
#include <stdint.h>

#if OCTAVE_MAJOR_VERSION < 5
// note: undef MPI_Comm and MPI_Info because of a conflict
// between some mpi implementations and hdf5 (used by NetCDF)
#undef MPI_Comm
#undef MPI_Info
#endif

#include <octave/oct.h>
#include <octave/parse.h>
#if (OCTAVE_MAJOR_VERSION >= 4 && OCTAVE_MINOR_VERSION >= 2 || OCTAVE_MAJOR_VERSION > 4)
#include <octave/interpreter.h>
#else
#include <octave/toplev.h>
#endif
#include <octave/octave.h>
#if !(OCTAVE_MAJOR_VERSION >= 4 && OCTAVE_MINOR_VERSION >= 4 || OCTAVE_MAJOR_VERSION > 4)
#include <octave/oct-alloc.h>
#endif
#include "module-octave.h"
#include "octave_object.h"

#if OCTAVE_MAJOR_VERSION < 5
#define isinteger() is_integer_type()
#define isreal() is_real_type()
#define isobject() is_object()
#define numel() length()
#define iscell() is_cell()
#endif

#ifdef DEBUG
#define TRACE(msg) ((void)(std::cerr << __FILE__ << ":" << __LINE__ << ":" << __PRETTY_FUNCTION__ << ":" << msg << std::endl))
#define ASSERT(expr) assert(expr)
#else
#define TRACE(msg) ((void)0)
#endif

namespace oct {

class OctaveInterface {
public:
        enum OctaveCallFlags_t {
                DEFAULT_CALL_FLAGS              = 0x0,
                PASS_DATA_MANAGER               = 0x1,
                UPDATE_OCTAVE_VARIABLES = 0x2,
                UPDATE_MBDYN_VARIABLES  = 0x4,
                UPDATE_GLOBAL_VARIABLES = UPDATE_OCTAVE_VARIABLES | UPDATE_MBDYN_VARIABLES,
                OPTIONAL_OUTPUT_ARGS    = 0x8
        };

private:
        explicit OctaveInterface(const DataManager* pDM, MBDynParser* pHP);
        virtual ~OctaveInterface(void);

public:
        static OctaveInterface* CreateInterface(const DataManager* pDM, MBDynParser* pHP);
        void AddRef() {
                ++iRefCount;
                TRACE("octave info: increase reference count to:" << iRefCount << "\n");                                
        }
        void Destroy();
        static OctaveInterface* GetInterface(void);
        void UpdateOctaveVariables(void);
        void UpdateMBDynVariables(void);
        const DataManager* GetDataManager(void) const{ return pDM; }
        MBDynParser* GetMBDynParser(void) const { return pHP; }
        const octave_value& GetDataManagerInterface(void) const { return octDM; }
        const octave_value& GetMBDynParserInterface(void) const { return octHP; }
        bool AddOctaveSearchPath(const std::string& path);
        static bool ConvertMBDynToOctave(const TypedValue& mbValue, octave_value& octValue);
        static bool ConvertMBDynToOctave(doublereal mbValue, octave_value& octValue); // useful mainly for templates
        static bool ConvertMBDynToOctave(const Vec3& mbValue, octave_value& octValue);
        static bool ConvertMBDynToOctave(const Vec6& mbValue, octave_value& octValue);
        static bool ConvertMBDynToOctave(const Mat3x3& mbValue, octave_value& octValue);
        static bool ConvertMBDynToOctave(const Mat6x6& mbValue, octave_value& octValue);
        static bool ConvertOctaveToMBDyn(const ColumnVector& octValue, doublereal& mbValue);
        static bool ConvertOctaveToMBDyn(const ColumnVector& octValue, Vec3& mbValue);
        static bool ConvertOctaveToMBDyn(const ColumnVector& octValue, Vec6& mbValue);
        static bool ConvertOctaveToMBDyn(const Matrix& octValue, doublereal& mbValue);
        static bool ConvertOctaveToMBDyn(const Matrix& octValue, Vec3& mbValue);
        static bool ConvertOctaveToMBDyn(const Matrix& octValue, Vec6& mbValue);
        static bool ConvertOctaveToMBDyn(const Matrix& octValue, Mat3x3& mbValue);
        static bool ConvertOctaveToMBDyn(const Matrix& octValue, Mat6x6& mbValue);
        static bool ConvertOctaveToMBDyn(const octave_value& octValue, TypedValue& mbValue);
        static bool ConvertOctaveToMBDyn(const octave_value& octValue, doublereal& mbValue);
        inline octave_value_list EvalFunction(const std::string& func, const octave_value_list& args, int nargout = 1, int flags = DEFAULT_CALL_FLAGS);
        inline octave_value_list EvalFunctionDerivative(const std::string& func, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
        inline octave_value_list EvalFunctionDerivative(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
        inline octave_value_list EvalFunction(const std::string& func, doublereal dVar, int nargout = 1, int flags = DEFAULT_CALL_FLAGS);
        inline doublereal EvalScalarFunction(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
        inline doublereal EvalScalarFunction(const std::string& func, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
        inline doublereal EvalScalarFunctionDerivative(const std::string& func, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
        inline doublereal EvalScalarFunctionDerivative(const std::string& func, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
        template<typename T>
        inline void EvalMatrixFunction(const std::string& func, T& res, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
        template <typename T>
        inline void EvalMatrixFunction(const std::string& func, T& res, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
        template<typename T>
        inline void EvalMatrixFunctionDerivative(const std::string& func, T& res, doublereal dVar, int flags = DEFAULT_CALL_FLAGS);
        template<typename T>
        inline void EvalMatrixFunctionDerivative(const std::string& func, T& res, const octave_value_list& args, int flags = DEFAULT_CALL_FLAGS);
        static bool HaveADPackage(void) { return bHaveADPackage; };
        void AddEmbedFileName(const std::string& strFile);
        bool bHaveMethod(const octave_value& octObject, const std::string& strClass, const std::string& strName);
        inline octave_value_list MakeArgList(doublereal dVar, const octave_value_list& args, int iFlags) const;
private:
        template<typename T>
        inline static bool ConvertOctaveToMBDyn(const ColumnVector& octValue, T& mbValue, int rows);
        template<typename T>
        inline static bool ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows);
        template<typename T>
        inline static bool ConvertOctaveToMBDyn(const Matrix& octValue, T& mbValue, int rows, int cols);
        bool LoadADPackage(void);
        template<typename T>
        inline static bool ConvertMBDynToOctave(const T& mbValue, octave_value& octValue, int rows);
        template<typename T>
        inline static bool ConvertMBDynToOctave(const T& mbValue, octave_value& octValue, int rows, int cols);
        static void exit(int);

private:
#if OCTAVE_MAJOR_VERSION >= 5
        octave::interpreter interpreter;
#endif
        typedef std::map<std::string, bool> EmbedFileNameMap_t;
        typedef EmbedFileNameMap_t::iterator EmbedFileNameIter_t;
        bool bFirstCall;
        bool bEmbedFileDirty;
        const octave_value octDM;
        const octave_value octHP;
        const DataManager* const pDM;
        MBDynParser* const pHP;
        EmbedFileNameMap_t strEmbedFileNames;
        static OctaveInterface* pOctaveInterface;
        static int iRefCount;
        static const std::string strADFunc;
        static const std::string strIsMethod;
        static bool bHaveADPackage;
};

class MBDynInterface : public octave_object {
public:
        // An instance of this class might be created directly from within octave
        // In this case pInterface must be defined
        explicit MBDynInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), bool bAddRef = true);
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
        const bool bAddRef;
};

class ConstVectorHandlerInterface : public MBDynInterface {
public:
        explicit ConstVectorHandlerInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const VectorHandler* pX=0);
        virtual ~ConstVectorHandlerInterface();
        void Set(const VectorHandler* pX){ this->pX = const_cast<VectorHandler*>(pX); }
        virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;
        virtual octave_value operator()(const octave_value_list& idx) const;
        virtual dim_vector dims (void) const;
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(dGetCoef)
                METHOD_DECLARE(GetVec)
                METHOD_DECLARE(iGetSize)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

protected:
        VectorHandler* pX;
};

class VectorHandlerInterface : public ConstVectorHandlerInterface {
public:
        explicit VectorHandlerInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), VectorHandler* X=0);
        virtual ~VectorHandlerInterface();
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(PutCoef)
                METHOD_DECLARE(PutVec)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

};

class OStreamInterface : public MBDynInterface {
public:
        explicit OStreamInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), std::ostream* pOS = 0);
        virtual ~OStreamInterface();
        void Set(std::ostream* pOS){ this->pOS = pOS; }
        std::ostream* Get()const{ return pOS; }
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(printf)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        std::ostream* pOS;
        static const std::string strsprintf;
};

class SimulationEntityInterface: public MBDynInterface {
public:
        explicit SimulationEntityInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface());
        virtual ~SimulationEntityInterface();
        virtual const SimulationEntity* Get()const=0;
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(iGetNumPrivData)
                METHOD_DECLARE(iGetPrivDataIdx)
                METHOD_DECLARE(dGetPrivData)
        END_METHOD_TABLE_DECLARE()
};

class NodeInterface: public SimulationEntityInterface {
public:
        explicit NodeInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface());
        virtual ~NodeInterface();
        virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;
        virtual const Node* Get()const=0;
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(GetLabel)
                METHOD_DECLARE(GetName)
                METHOD_DECLARE(iGetFirstIndex)
                METHOD_DECLARE(iGetFirstRowIndex)
                METHOD_DECLARE(iGetFirstColIndex)
                METHOD_DECLARE(dGetDofValue)
                METHOD_DECLARE(dGetDofValuePrev)
        END_METHOD_TABLE_DECLARE()
};

class ScalarNodeInterface: public NodeInterface {
public:
        explicit ScalarNodeInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), ScalarNode* pNode = 0);
        virtual ~ScalarNodeInterface();
        virtual const ScalarNode* Get()const;

protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(SetX)
                METHOD_DECLARE(dGetX)
                METHOD_DECLARE(SetXPrime)
                METHOD_DECLARE(dGetXPrime)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        ScalarNode* const pNode;
};

class StructDispNodeBaseInterface : public NodeInterface {
public:
        explicit StructDispNodeBaseInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface());
        virtual ~StructDispNodeBaseInterface();
        virtual const StructDispNode* Get()const=0;

protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(iGetFirstPositionIndex)
                METHOD_DECLARE(iGetFirstMomentumIndex)
                METHOD_DECLARE(GetLabel)
                METHOD_DECLARE(GetXCurr)
                METHOD_DECLARE(GetXPrev)
                METHOD_DECLARE(GetVCurr)
                METHOD_DECLARE(GetVPrev)
                METHOD_DECLARE(GetXPPCurr)
                METHOD_DECLARE(GetXPPPrev)
        END_METHOD_TABLE_DECLARE()

protected:
        octave_value GetVec3(const Vec3& v, const octave_value_list& args) const;
        octave_value GetMat3x3(const Mat3x3& m, const octave_value_list& args) const;
};

class StructDispNodeInterface: public StructDispNodeBaseInterface {
public:
        explicit StructDispNodeInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const StructDispNode* pNode = 0);
        virtual ~StructDispNodeInterface();
        virtual const StructDispNode* Get()const;

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        const StructDispNode* const pNode;
};

class StructNodeInterface : public StructDispNodeBaseInterface {
public:
        explicit StructNodeInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const StructNode* pNode = 0);
        virtual ~StructNodeInterface();
        virtual const StructNode* Get()const;
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(GetgCurr)
                METHOD_DECLARE(GetgRef)
                METHOD_DECLARE(GetgPCurr)
                METHOD_DECLARE(GetgPRef)
                METHOD_DECLARE(GetRCurr)
                METHOD_DECLARE(GetRPrev)
                METHOD_DECLARE(GetRRef)
                METHOD_DECLARE(GetWCurr)
                METHOD_DECLARE(GetWPrev)
                METHOD_DECLARE(GetWRef)
                METHOD_DECLARE(GetWPCurr)
                METHOD_DECLARE(GetWPPrev)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        const StructNode* const pNode;
};

class DataManagerInterface : public MBDynInterface {
public:
        // An instance of this class might be created directly from within octave
        // In this case pInterface must be defined
        explicit DataManagerInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const DataManager* pDM = 0, bool bAddRef = true);
        virtual ~DataManagerInterface();
        virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;
        inline const DataManager* GetDataManager(void)const;
        inline const Table& GetSymbolTable(void)const;

protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(GetVariable)
                METHOD_DECLARE(GetStructNodePos)
                METHOD_DECLARE(GetStructNode)
                METHOD_DECLARE(pFindNode)
                METHOD_DECLARE(ReadNode)
                METHOD_DECLARE(dGetTime)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        Node* pFindNode_(const octave_value_list& args, Node::Type& Type) const;
        Node::Type GetNodeType(const std::string& strType) const;
        NodeInterface* CreateNodeInterface(Node* pNode, Node::Type Type) const;

private:
        const DataManager* const pDM;
};

class MBDynParserInterface : public MBDynInterface {
public:
        // An instance of this class might be created directly from within octave
        // In this case pInterface must be defined
        explicit MBDynParserInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), MBDynParser* pHP=0, bool bAddRef=true);
        virtual ~MBDynParserInterface();
        virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;
        MBDynParser* Get()const{ return pHP; }

protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(IsKeyWord)
                METHOD_DECLARE(IsArg)
                METHOD_DECLARE(IsStringWithDelims)
                METHOD_DECLARE(GetReal)
                METHOD_DECLARE(GetInt)
                METHOD_DECLARE(GetBool)
                METHOD_DECLARE(GetYesNoOrBool)
                METHOD_DECLARE(GetString)
                METHOD_DECLARE(GetStringWithDelims)
                METHOD_DECLARE(GetValue)
                METHOD_DECLARE(GetPosRel)
                METHOD_DECLARE(GetRotRel)
                METHOD_DECLARE(GetLineData)
                METHOD_DECLARE(GetDriveCaller)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

private:
        const StructDispNode* GetStructNode(const octave_value& arg)const;
        static bool GetDelimsFromString(const std::string& strDelims, HighParser::Delims& delims);
private:
        MBDynParser* pHP;
        static const struct MBDynStringDelims {
                HighParser::Delims value;
                char name[15];
        } mbStringDelims[6];
};

class DriveCallerInterface : public MBDynInterface {
public:
        explicit DriveCallerInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), const DriveCaller* pDC=0);
        virtual ~DriveCallerInterface();
        void Set(const DriveCaller* pDC){ DC.Set(pDC); }
protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(dGet)
                METHOD_DECLARE(dGetP)
        END_METHOD_TABLE_DECLARE()

private:
        DECLARE_OV_TYPEID_FUNCTIONS_AND_DATA
        DECLARE_OCTAVE_ALLOCATOR

protected:
        DriveOwner DC;
};

class OctaveDriveCaller : public DriveCaller {
public:
        explicit
        OctaveDriveCaller(const std::string& strFunc,
                        OctaveInterface* pInterface,
                        int iFlags,
                        const octave_value_list& args);
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
        inline octave_value_list MakeArgList(doublereal dVar) const;

private:
        const int iFlags;
        const std::string strFunc;
        OctaveInterface* const pInterface;
        const octave_value_list args;
};

template <class T>
class OctaveTplDriveCaller : public TplDriveCaller<T> {
public:
        OctaveTplDriveCaller(const std::string& strFunction, OctaveInterface* pInterface, int iFlags, const octave_value_list& args);
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
        inline octave_value_list MakeArgList(doublereal dVar) const;

private:
        const std::string strFunction;
        OctaveInterface* const pInterface;
        const int iFlags;
        const octave_value_list args;
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
        OctaveScalarFunction(const std::string& strFunc, OctaveInterface* pInterface, int iFlags, const octave_value_list& args);
        virtual ~OctaveScalarFunction(void);
        virtual doublereal operator()(const doublereal x) const;
        virtual doublereal ComputeDiff(const doublereal t, const integer order = 1) const;

private:
        inline octave_value_list MakeArgList(doublereal dVar) const;

private:
        const std::string strFunc;
        OctaveInterface* const pInterface;
        const int iFlags;
        const octave_value_list args;
};

class OctaveConstitutiveLawBase
{
public:
        explicit inline OctaveConstitutiveLawBase(const std::string& strClass, OctaveInterface* pInterface, int iFlags);
        virtual ~OctaveConstitutiveLawBase();
        inline bool bHaveMethod(const std::string& strName) const;
        OctaveInterface* GetInterface()const{ return pInterface; }
        const std::string& GetClass()const{ return strClass; }
        int GetFlags()const{ return iFlags; }

protected:
        octave_value octObject;

private:
        const std::string strClass;
        OctaveInterface* const pInterface;
        const int iFlags;
        static const std::string strGetConstLawType;
        static const std::string strUpdate;
};

template <class T, class Tder>
class OctaveConstitutiveLaw
: public ConstitutiveLaw<T, Tder>, private OctaveConstitutiveLawBase {
        typedef ConstitutiveLaw<T, Tder> Base_t;
public:
        OctaveConstitutiveLaw(const std::string& strClass, OctaveInterface* pInterface, int iFlags);
        virtual ~OctaveConstitutiveLaw(void);
        ConstLawType::Type GetConstLawType(void) const;
        virtual ConstitutiveLaw<T, Tder>* pCopy(void) const;
        virtual std::ostream& Restart(std::ostream& out) const;
        virtual void Update(const T& mbEps, const T& mbEpsPrime);

private:
        mutable ConstLawType::Type clType;
};

class OctaveBaseDCR {
public:
        OctaveBaseDCR(void);
        virtual ~OctaveBaseDCR(void);
        void Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred = false) const;
        OctaveInterface* GetInterface(void) const { return pInterface; };
        const std::string& GetFunction(void) const { return strFunction; };
        int GetFlags(void) const { return iFlags; };

private:
        // FIXME: mutable needed for ScalarFunctionRead
        mutable OctaveInterface* pInterface;
        mutable std::string strFunction;
        mutable int iFlags;
};

class OctaveFunctionDCR: public OctaveBaseDCR {
public:
        void Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred = false) const;
        const octave_value_list& GetArgs()const{ return args; }

private:
        mutable octave_value_list args;
};

class OctaveDCR : public DriveCallerRead, public OctaveFunctionDCR {
public:
        virtual DriveCaller *
        Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

template <class T>
class OctaveTDCR : public TplDriveCallerRead<T>, public OctaveFunctionDCR {
public:
        virtual TplDriveCaller<T> *
        Read(const DataManager* pDM, MBDynParser& HP);
};

class DerivativeDCR : public DriveCallerRead {
public:
        virtual DriveCaller *
        Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

class OctaveSFR : public ScalarFunctionRead, public OctaveFunctionDCR {
public:
        const BasicScalarFunction *
        Read(DataManager* pDM, MBDynParser& HP) const;
};

template <class T, class Tder>
struct OctaveCLR : public ConstitutiveLawRead<T, Tder>, public OctaveBaseDCR {
        virtual ConstitutiveLaw<T, Tder> *
        Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
                ConstitutiveLaw<T, Tder>* pCL = 0;

                OctaveBaseDCR::Read(pDM, HP);

                typedef OctaveConstitutiveLaw<T, Tder> L;
                SAFENEWWITHCONSTRUCTOR(pCL, L, L(GetFunction(), GetInterface(), GetFlags()));

                CLType = pCL->GetConstLawType();

                return pCL;
        };
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
        virtual unsigned int iGetPrivDataIdx(const char *s) const;
        virtual doublereal dGetPrivData(unsigned int i) const;
        int iGetNumConnectedNodes(void) const;
        void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
        void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
                SimulationEntity::Hints *ph);
        virtual unsigned int iGetNumDof(void) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
        virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
        virtual void AfterConvergence(const VectorHandler& X,
                        const VectorHandler& XP);
        std::ostream& Restart(std::ostream& out) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual void
        InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        VariableSubMatrixHandler&
        InitialAssJac(VariableSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr);
        SubVectorHandler&
        InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
        virtual void SetInitialValue(VectorHandler& X);
        OctaveInterface* GetInterface(void) const { return dcr.GetInterface(); };
        int GetFlags(void) const { return dcr.GetFlags(); };
        const std::string& GetClass(void) const { return dcr.GetFunction(); };

private:
        inline VariableSubMatrixHandler&
        AssMatrix(VariableSubMatrixHandler& WorkMatVar, const octave_value_list& ans, bool bInitial);
        inline VariableSubMatrixHandler&
        AssMatrix(VariableSubMatrixHandler& WorkMatVar, const Matrix& Jac, const int32NDArray& ridx, const int32NDArray& cidx, bool bSparse, bool bInitial);
        bool bHaveMethod(const std::string& strName)const;
private:
        octave_value octObject;
        octave_value mbdObject;
        OctaveBaseDCR dcr;
        enum OctaveMethods_t {
                HAVE_DEFAULT               = 0x0000,
                HAVE_JACOBIAN              = 0x0001,
                HAVE_UPDATE                        = 0x0002,
                HAVE_SET_VALUE             = 0x0004,
                HAVE_PRIVATE_DOF           = 0x0008,
                HAVE_INITIAL_ASSEMBLY  = 0x0010,
                HAVE_SET_INITIAL_VALUE = 0x0020,
                HAVE_AFTER_CONVERGENCE = 0x0040,
                HAVE_PRIVATE_DATA          = 0x0080,
                HAVE_CONNECTED_NODES   = 0x0100,
                HAVE_OUTPUT                        = 0x0200,
                HAVE_DESCRIBE_DOF          = 0x0400,
                HAVE_DESCRIBE_EQ           = 0x0800,
                HAVE_AFTER_PREDICT         = 0x1000,
                HAVE_RESTART               = 0x2000
        };
        int haveMethod;
        octave_object_ptr<ConstVectorHandlerInterface> X;
        octave_object_ptr<ConstVectorHandlerInterface> XP;
        octave_object_ptr<OStreamInterface> OS;
        static const std::string strWorkSpaceDim;
        static const std::string striGetNumDof;
        static const std::string strAssRes;
        static const std::string strAssJac;
        static const std::string strUpdate;
        static const std::string strSetValue;
        static const std::string striGetInitialNumDof;
        static const std::string strSetInitialValue;
        static const std::string strInitialAssRes;
        static const std::string strInitialAssJac;
        static const std::string strInitialWorkSpaceDim;
        static const std::string strGetDofType;
        static const std::string strGetEqType;
        static const std::string strAfterConvergence;
        static const std::string striGetNumPrivData;
        static const std::string striGetPrivDataIdx;
        static const std::string strdGetPrivData;
        static const std::string striGetNumConnectedNodes;
        static const std::string strGetConnectedNodes;
        static const std::string strOutput;
        static const std::string strDescribeDof;
        static const std::string strDescribeEq;
        static const std::string strAfterPredict;
        static const std::string strRestart;
};

class OctaveElementInterface: public MBDynInterface {
public:
        explicit OctaveElementInterface(OctaveInterface* pInterface = OctaveInterface::GetInterface(), OctaveElement* pElem = 0);
        virtual ~OctaveElementInterface(void);
        virtual void print(std::ostream& os, bool pr_as_read_syntax = false) const;

protected:
        BEGIN_METHOD_TABLE_DECLARE()
                METHOD_DECLARE(GetLabel)
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
const std::string OctaveInterface::strIsMethod("ismethod");
bool OctaveInterface::bHaveADPackage = false;

OctaveInterface::OctaveInterface(const DataManager* pDM, MBDynParser* pHP)
: bFirstCall(true),
bEmbedFileDirty(false),
octDM(new DataManagerInterface(this, pDM, false)),
octHP(new MBDynParserInterface(this, pHP, false)),
pDM(pDM),
pHP(pHP)
{
        TRACE("constructor");

        ASSERT(iRefCount == 0);
        
        ASSERT(pDM != 0);
        ASSERT(pHP != 0);
        ASSERT(pOctaveInterface == 0);

        pOctaveInterface = this;

#if OCTAVE_MAJOR_VERSION >= 5
        int status = interpreter.execute ();

        if (status != 0) {
                silent_cerr("module-octave: creating embedded Octave interpreter failed!" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

#else
        const int nmax_args = 4;
        const int nmax_char = 18;
        char args[nmax_char] = "octave-cli"; //FIXME: this might be unsafe
        char *argv[nmax_args] = { &args[0] } ;
        int argc = 1;

        char* p = strchr(args, '\0') + 1;

        if (silent_output) {
                strncat(p, "-q", &args[nmax_char] - p - 2);
                if (argc < nmax_args - 1)
                        argv[argc++] = p;
                p = strchr(p, '\0') + 1;
        }

        if (pedantic_output) {
                strncat(p, "-V", &args[nmax_char] - p - 2);
                if (argc < nmax_args - 1)
                        argv[argc++] = p;
                p = strchr(p, '\0') + 1;
        }

#if defined(DEBUG)
        int i;
        char** pp;
        for (i = 0, pp = argv; *pp; ++pp, ++i) {
                silent_cerr("octave>argv[" << i << "]=\"" << *pp << "\"" << std::endl);
        }
#endif
        octave_main(argc, argv, 1);
#endif

#if OCTAVE_MAJOR_VERSION >= 5
        octave::feval("page_screen_output", octave_value(false), 0); // turn off pager
#else
        feval("page_screen_output", octave_value(false), 0); // turn off pager
#endif

#if OCTAVE_MAJOR_VERSION < 6
        if ((error_state != 0)) {
                silent_cerr("warning: page_screen_output failed" << std::endl);
                error_state = 0; // ignore error
        }
#endif

        LoadADPackage();
}

void OctaveInterface::exit(int)
{
  silent_cerr("octave interface has been uninitialized" << std::endl);
}

OctaveInterface::~OctaveInterface(void)
{
        TRACE("destructor");

        ASSERT(this == pOctaveInterface);

#if !(OCTAVE_MAJOR_VERSION >= 4 && OCTAVE_MINOR_VERSION >= 4 || OCTAVE_MAJOR_VERSION > 4)
        octave_exit = &OctaveInterface::exit;
#elif OCTAVE_MAJOR_VERSION < 5
#if defined(HAVE_DO_OCTAVE_ATEXIT)
        do_octave_atexit();
#elif defined(HAVE_CLEAN_UP_AND_EXIT)
        clean_up_and_exit(0, true);
#endif        
#elif OCTAVE_MAJOR_VERSION >= 6
        interpreter.shutdown();
#endif

        pOctaveInterface = 0;
}

bool
OctaveInterface::LoadADPackage(void)
{
#if OCTAVE_MAJOR_VERSION >= 4 && OCTAVE_MINOR_VERSION >= 2 || OCTAVE_MAJOR_VERSION > 4
        // AD package is broken on this octave version
        // Use finite differences in var/mbdyn_derivative.m if AD is not available
        bHaveADPackage = true;
#else
        octave_value_list args;
        args.append(octave_value("load"));
        args.append(octave_value("ad"));

        feval("pkg", args, 0);

        bHaveADPackage = true; // Use finite differences in var/mbdyn_derivative.m if AD is not available

#if OCTAVE_MAJOR_VERSION < 6
        if ((error_state != 0)) {
            silent_cerr("warning: octave package for automatic forward differentiation is not available" << std::endl);
            error_state = 0; // ignore error
            bHaveADPackage = false;
        }
#endif
#endif
        return bHaveADPackage;
}

OctaveInterface *
OctaveInterface::CreateInterface(const DataManager* pDM, MBDynParser* pHP)
{
        if (!pOctaveInterface) {
                pOctaveInterface = new OctaveInterface(pDM, pHP);
        }
        
        pOctaveInterface->AddRef();
        
        return pOctaveInterface;
}

void
OctaveInterface::Destroy(void)
{
        ASSERT(iRefCount >= 1);

        if (0 == --iRefCount) {
                TRACE("octave info: Octave interface is released\n"); 
                delete this;
        }
        
        TRACE("octave info: Octave remaining count:" << iRefCount << "\n");        
}

OctaveInterface* OctaveInterface::GetInterface(void) {
        return pOctaveInterface;
}

void
OctaveInterface::UpdateOctaveVariables(void)
{
        typedef Table::VM::const_iterator iterator;

        const Table& symbolTable = GetDataManager()->GetMathParser().GetSymbolTable();

        for (iterator it = symbolTable.begin(); it != symbolTable.end(); ++it)
        {
                const std::string& mbName(it->first);
                const NamedValue*const namedValue = it->second;
                const TypedValue mbValue(namedValue->GetVal());
                octave_value octValue;

                if (!ConvertMBDynToOctave(mbValue, octValue)) {
                        silent_cerr("octave error: data type \"" << mbValue.GetTypeName() << "\" of variable \"" << mbName << "\": not handled in switch statement " << std::endl);
                        ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

#if OCTAVE_MAJOR_VERSION >= 6
                interpreter.global_assign(mbName, octValue);
#elif OCTAVE_MAJOR_VERSION >= 5
                interpreter.get_symbol_table().global_assign(mbName, octValue);
#else
                set_global_value(mbName, octValue);
#endif
        }
}

void
OctaveInterface::UpdateMBDynVariables(void)
{
        typedef Table::VM::const_iterator iterator;

        Table& symbolTable = GetDataManager()->GetMathParser().GetSymbolTable();

        for (iterator it = symbolTable.begin(); it != symbolTable.end(); ++it)
        {
                Var* const varValue = dynamic_cast<Var*>(it->second);

                if (!varValue || varValue->Const()) {
                        continue;
                }

                const std::string& mbName(it->first);

#if OCTAVE_MAJOR_VERSION >= 6
                const octave_value octValue(interpreter.global_varval(mbName));                
#elif OCTAVE_MAJOR_VERSION >= 5
                const octave_value octValue(interpreter.get_symbol_table().global_varval(mbName));
#else
                const octave_value octValue(get_global_value(mbName));
#endif

                if (!octValue.is_defined()) {
                        pedantic_cerr("octave warning: global variable " << mbName << " is not defined in octave" << std::endl);
                        continue;
                }

                TypedValue mbValue;

                if (!ConvertOctaveToMBDyn(octValue, mbValue)) {
                        silent_cerr("octave error: data type \"" << octValue.type_name() << "\" of variable \"" << mbName << "\": can not be converted into MBDyn format " << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                varValue->SetVal(mbValue);
        }
}

bool
OctaveInterface::AddOctaveSearchPath(const std::string& path)
{
#if OCTAVE_MAJOR_VERSION >= 5
        octave::feval("addpath", octave_value(path));
#else
        feval("addpath", octave_value(path));
#endif
#if OCTAVE_MAJOR_VERSION < 6
        return error_state == 0;
#else
        return true;
#endif
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
OctaveInterface::ConvertMBDynToOctave(doublereal mbValue, octave_value& octValue)
{
        // This function is especially useful for template based constitutive laws
        octValue = mbValue;

        return true;
}

bool
OctaveInterface::ConvertMBDynToOctave(const Vec3& mbValue, octave_value& octValue)
{
        return ConvertMBDynToOctave(mbValue, octValue, 3);
}

bool
OctaveInterface::ConvertMBDynToOctave(const Vec6& mbValue, octave_value& octValue)
{
        return ConvertMBDynToOctave(mbValue, octValue, 6);
}

bool
OctaveInterface::ConvertMBDynToOctave(const Mat3x3& mbValue, octave_value& octValue)
{
        return ConvertMBDynToOctave(mbValue, octValue, 3, 3);
}

bool
OctaveInterface::ConvertMBDynToOctave(const Mat6x6& mbValue, octave_value& octValue)
{
        return ConvertMBDynToOctave(mbValue, octValue, 6, 6);
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const ColumnVector& octValue, doublereal& mbValue)
{
        if (octValue.numel() != 1) {
                silent_cerr("octave error: invalid vector size " << octValue.numel() << " expected 1" << std::endl);
                return false;
        }

        mbValue = octValue(0);

        return true;
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const ColumnVector& octValue, Vec3& mbValue)
{
        return ConvertOctaveToMBDyn(octValue, mbValue, 3);
}

bool
OctaveInterface::ConvertOctaveToMBDyn(const ColumnVector& octValue, Vec6& mbValue)
{
        return ConvertOctaveToMBDyn(octValue, mbValue, 6);
}

template<typename T>
bool
OctaveInterface::ConvertOctaveToMBDyn(const ColumnVector& octValue, T& mbValue, int rows)
{
        if (octValue.numel() != rows) {
                silent_cerr("octave error: invalid vector size " << octValue.numel() << " expected " << rows << std::endl);
                return false;
        }

        for (int i = 0; i < rows; ++i) {
                mbValue(i + 1) = octValue(i);
        }

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

template<typename T>
bool
OctaveInterface::ConvertMBDynToOctave(const T& mbValue, octave_value& octValue, int rows)
{
        ColumnVector V(rows);

        for (int i = 0; i < rows; ++i) {
                V(i) = mbValue(i + 1);
        }

        octValue = V;

        return true;
}

template<typename T>
bool
OctaveInterface::ConvertMBDynToOctave(const T& mbValue, octave_value& octValue, int rows, int cols)
{
        Matrix M(rows, cols);

        for (int i = 0; i < rows; ++i) {
                for (int j = 0; j < cols; ++j) {
                        M(i, j) = mbValue(i + 1, j + 1);
                }
        }

        octValue = M;

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

bool
OctaveInterface::ConvertOctaveToMBDyn(const octave_value& octValue, TypedValue& mbValue)
{
        if (!octValue.is_scalar_type()) {
                return false;
        }

        if (octValue.is_bool_scalar()) {
                mbValue = TypedValue(octValue.bool_value());
        } else if (octValue.is_int32_type()) {
                mbValue = TypedValue(static_cast<int32_t>(octValue.int32_scalar_value()));
        } else if (octValue.is_real_scalar()) {
                mbValue = TypedValue(octValue.scalar_value());
        } else if (octValue.is_string()) {
                mbValue = TypedValue(octValue.string_value());
        } else {
                return false;
        }

        return true;
}

bool OctaveInterface::ConvertOctaveToMBDyn(const octave_value& octValue, doublereal& mbValue)
{
        if (!octValue.is_real_scalar()) {
                return false;
        }

        mbValue = octValue.scalar_value();

        return true;
}

octave_value_list
OctaveInterface::EvalFunction(const std::string& func, doublereal dVar, int nargout, int flags)
{
        octave_value_list args = octave_value(dVar);

        if (flags & PASS_DATA_MANAGER) {
                args.append(octDM);
        }

        return EvalFunction(func, args, nargout, flags);
}

octave_value_list
OctaveInterface::EvalFunction(const std::string& func, const octave_value_list& args, int nargout, int flags)
{
        if ((flags & UPDATE_OCTAVE_VARIABLES) || bFirstCall) {
                UpdateOctaveVariables();
        }

        if (bFirstCall) {
            octave_value_list args;
            args.append(octave_value("load"));
            args.append(octave_value("mbdyn_util_oct"));

#if OCTAVE_MAJOR_VERSION >= 5
            octave::feval("pkg", args, 0);
#else
            feval("pkg", args, 0);
#endif

#if OCTAVE_MAJOR_VERSION < 6
            if ((error_state != 0)) {
                silent_cerr("warning: octave package mbdyn_util_oct has not been installed" << std::endl);

                error_state = 0;
            }
#endif
            // octave .m files by default are installed here
            if (!AddOctaveSearchPath(OCTAVEPATH)) {
                    silent_cerr("OctaveInterface error: addpath(\"" << OCTAVEPATH << "\") failed" << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        bFirstCall = false;

        if ( bEmbedFileDirty ) {
                for (EmbedFileNameIter_t pFile = strEmbedFileNames.begin();
                        pFile != strEmbedFileNames.end(); ++pFile )
                {
                        if (!pFile->second) {
#if OCTAVE_MAJOR_VERSION >= 5
                                try {
                                        octave::feval("source", octave_value(pFile->first));
                                } catch (const octave::exit_exception& err) {
                                        silent_cerr("module-octave: Octave interpreter exited with status "
                                                    << err.exit_status() << std::endl);
                                        throw NoErr(MBDYN_EXCEPT_ARGS);
                                } catch (const octave::execution_exception& err) {
                                        silent_cerr("module-octave: Error encountered in file \""
                                                    << pFile->first << "\"" << std::endl);
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }
#else
                                feval("source", octave_value(pFile->first));
#endif
#if OCTAVE_MAJOR_VERSION < 6
                                if (error_state) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }
#endif
                                pFile->second = true;
                        }
                }
                bEmbedFileDirty = false;
        }

#if OCTAVE_MAJOR_VERSION >= 5
        octave_value_list ans;

        try {
                ans = octave::feval(func, args, nargout);
        } catch (const octave::exit_exception& err) {
                silent_cerr("module-octave: Octave interpreter exited with status "
                            << err.exit_status() << std::endl);
                throw NoErr(MBDYN_EXCEPT_ARGS);
        } catch (const std::exception& err) {
                silent_cerr("Exception in Octave function \"" << func << "\": " << err.what() << "\n");
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
#else
        octave_value_list ans = feval(func, args, nargout);
#endif
#if OCTAVE_MAJOR_VERSION < 6
        if (error_state) {
                // An error message has been displayed by octave
                // There is no need to output further error information
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
#endif
        if (flags & OPTIONAL_OUTPUT_ARGS) {
                if (ans.length() > nargout) {
                                silent_cerr("octave error: function \"" << func << "\" returned " << ans.length() << " output parameters\n"
                                                "expected maximum " << nargout  << " output parameters" << std::endl);
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }
        else if (ans.length() != nargout) {
                silent_cerr("octave error: function \"" << func << "\" returned " << ans.length() << " output parameters\n"
                                "expected " << nargout  << " output parameters" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        for (int i = 0; i < ans.length(); ++i) {
                if (!ans(i).is_defined()) {
                        silent_cerr("octave error: result " << i + 1 << " of function \"" << func << "\" is undefined" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }

        if (flags & UPDATE_MBDYN_VARIABLES) {
                UpdateMBDynVariables();
        }

        return ans;
}

octave_value_list
OctaveInterface::EvalFunctionDerivative(const std::string& func, const octave_value_list& args, int flags)
{
        TRACE("func=" << func);
        TRACE("flags=" << flags);

        octave_value_list derivArgs = octave_value(func);
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
                args.append(octDM);
        }

        return EvalFunctionDerivative(func, args, flags);
}

inline doublereal
OctaveInterface::EvalScalarFunction(const std::string& func, doublereal dVar, int flags)
{
        return EvalScalarFunction(func, octave_value(dVar), flags);
}

inline doublereal
OctaveInterface::EvalScalarFunction(const std::string& func, const octave_value_list& args, int flags)
{
        const octave_value_list ans = EvalFunction(func, args, 1, flags);

        if (!ans(0).is_real_scalar()) {
                silent_cerr("octave error: result of function \"" << func << "\" is not a scalar value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return ans(0).scalar_value();
}

inline doublereal
OctaveInterface::EvalScalarFunctionDerivative(const std::string& func, doublereal dVar, int flags)
{
        return EvalScalarFunctionDerivative(func, octave_value(dVar), flags);
}

inline doublereal
OctaveInterface::EvalScalarFunctionDerivative(const std::string& func, const octave_value_list& args, int flags)
{
        TRACE("func=" << func);
        TRACE("flags=" << flags);

        octave_value_list ans = EvalFunctionDerivative(func, args, flags);

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
        EvalMatrixFunction(func, res, octave_value(dVar), flags);
}

template <typename T>
inline void
OctaveInterface::EvalMatrixFunction(const std::string& func, T& res, const octave_value_list& args, int flags)
{
        const octave_value_list ans = EvalFunction(func, args, 1, flags);

        ASSERT(ans.length() == 1);
        // accept also scalar values for scalar template drive caller
        if (!(ans(0).is_real_matrix() || ans(0).is_real_scalar())) {
                silent_cerr("octave error: result of function \"" << func << "\" is not a matrix value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Matrix A(ans(0).matrix_value());

        if (!ConvertOctaveToMBDyn(A, res)) {
                silent_cerr("octave error: conversion of octave matrix failed!" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

template <typename T>
inline void
OctaveInterface::EvalMatrixFunctionDerivative(const std::string& func, T& res, doublereal dVar, int flags)
{
        EvalMatrixFunctionDerivative(func, res, octave_value(dVar), flags);
}

template<typename T>
inline void
OctaveInterface::EvalMatrixFunctionDerivative(const std::string& func, T& res, const octave_value_list& args, int flags)
{
        const octave_value_list ans = EvalFunctionDerivative(func, args, flags);

        ASSERT(ans.length() == 2);

        if (!(ans(0).is_real_matrix() || ans(0).is_real_scalar())) {
                silent_cerr("octave error: result of function \"" << func << "\" is not a matrix value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!(ans(1).is_real_matrix() || ans(1).is_real_scalar())) {
                silent_cerr("octave error: result of derivative of function \"" << func << "\" is not a matrix value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Matrix A(ans(1).matrix_value());

        if (!ConvertOctaveToMBDyn(A, res)) {
                silent_cerr("octave error: conversion of octave matrix failed!" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

void
OctaveInterface::AddEmbedFileName(const std::string& strFile)
{
        strEmbedFileNames.insert(std::pair<std::string, bool>(strFile,false));
        bEmbedFileDirty = true;
}

bool OctaveInterface::bHaveMethod(const octave_value& octObject, const std::string& strClass, const std::string& strName)
{
        octave_value_list args(octObject);
        args.append(octave_value(strName));

        octave_value_list ans = EvalFunction(strIsMethod, args, 1);

        ASSERT(ans.length() == 1);

        if (!ans(0).is_bool_scalar()) {
                silent_cerr("octave error: unexpected error in function " << strIsMethod << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return ans(0).bool_value(true);
}

octave_value_list OctaveInterface::MakeArgList(doublereal dVar, const octave_value_list& args, int iFlags) const
{
        octave_value_list fargs = octave_value(dVar);
        fargs.append(args);

        if (iFlags & PASS_DATA_MANAGER) {
                fargs.append(GetDataManagerInterface());
        }

        return fargs;
}

MBDynInterface::MBDynInterface(OctaveInterface* pInterface, bool bAddRef)
: pInterface(pInterface), bAddRef(bAddRef)
{
        TRACE("constructor");
        ASSERT(pInterface != 0);
        
        if (bAddRef) {
                pInterface->AddRef(); // Avoid invalid memory access in octave-6.2.1
        }
}

MBDynInterface::~MBDynInterface(void)
{
        TRACE("destructor");
        
        if (bAddRef) {
                pInterface->Destroy(); // Avoid invalid memory access in octave-6.2.1
        }
}

void
MBDynInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
        os << type_name() << ": MBDyn version" << VERSION << std::endl;
}

METHOD_DEFINE(MBDynInterface, GetVersion, args, nargout)
{
        octave_value octValue;

        if (args.length() != 0) {
                error("%s: invalid number of arguments: %ld",
                                type_name().c_str(),
                                long(args.length()));
                return octValue;
        }

        octValue = VERSION;

        return octValue;
}

BEGIN_METHOD_TABLE(MBDynInterface, octave_object)
        METHOD_DISPATCH(MBDynInterface, GetVersion)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(MBDynInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(MBDynInterface, "MBDyn", "MBDyn");

ConstVectorHandlerInterface::ConstVectorHandlerInterface(OctaveInterface* pInterface, const VectorHandler* X)
        :MBDynInterface(pInterface),
         pX(const_cast<VectorHandler*>(X))
{
        TRACE("constructor");
}

ConstVectorHandlerInterface::~ConstVectorHandlerInterface()
{
        TRACE("destructor");
}

void
ConstVectorHandlerInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return;
        }

        for ( int i = 1; i <= pX->iGetSize(); ++i ) {
                os << '\t' << pX->dGetCoef(i) << std::endl;
        }
}

octave_value ConstVectorHandlerInterface::operator()(const octave_value_list& idx) const
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (idx.length() != 1) {
                error("%s: invalid number of indices %ld", type_name().c_str(), long(idx.length()));
                return octave_value();
        }

        const integer iSize = pX->iGetSize();

        ColumnVector X;

        if (idx(0).is_magic_colon()) {
                X.resize(iSize);

                for (int i = 0; i < X.numel(); ++i) {
                        X(i) = (*pX)(i + 1);
                }
        } else {
                if (!(idx(0).is_range()
                                || (idx(0).isinteger()
                                        && (idx(0).rows() == 1 || idx(0).columns() == 1)))) {
                        error("%s: invalid index type %ldx%ld (%s)\n"
                                        "expected integer vector",
                                        type_name().c_str(),
                              static_cast<long>(idx(0).rows()),
                              static_cast<long>(idx(0).columns()),
                                        idx(0).type_name().c_str());
                        return octave_value();
                }

                const int32NDArray iRow(idx(0).int32_array_value());

#if OCTAVE_MAJOR_VERSION < 6
                if (error_state) {
                        return octave_value();
                }
#endif

                X.resize(iRow.numel());

                for (int i = 0; i < iRow.numel(); ++i) {
                        if (int32_t(iRow(i)) < 1 || int32_t(iRow(i)) > iSize) {
                                error("%s: index %d out of range [%d-%d]",
                                                type_name().c_str(),
                                                int32_t(iRow(i)),
                                                1,
                                                int32_t(iSize));

                                return octave_value();
                        }

                        X(i) = (*pX)(int32_t(iRow(i)));
                }
        }

        return octave_value(X);
}

dim_vector ConstVectorHandlerInterface::dims() const
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return dim_vector(-1, -1);
        }

        return dim_vector(pX->iGetSize(), 1);
}

METHOD_DEFINE(ConstVectorHandlerInterface, dGetCoef, args, nargout)
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if ( args.length() != 1 ) {
                error("%s: invalid number of arguments %ld\n"
                          "one argument expected",
                          type_name().c_str(),
                          long(args.length()));
                return octave_value();
        }

        if ( !(args(0).isinteger() && args(0).is_scalar_type()) ) {
                error("%s: invalid argument type \"%s\"\n"
                          "iRow must be an integer",
                          type_name().c_str(),
                          args(0).type_name().c_str());
                return octave_value();
        }

        const integer iRow = args(0).int32_scalar_value();

        if ( iRow < 1 || iRow > pX->iGetSize() ) {
                error("VectorHandler: index %ld out of range [%ld-%ld]", long(iRow), 1L, long(pX->iGetSize()));
                return octave_value();
        }

        return octave_value(pX->dGetCoef(iRow));
}

METHOD_DEFINE(ConstVectorHandlerInterface, GetVec, args, nargout)
{
        if (pX == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments %ld\n"
                          "one argument expected",
                          type_name().c_str(),
                          long(args.length()));
                return octave_value();
        }

        if (!(args(0).isinteger() && args(0).columns() == 1)) {
                error("%s: invalid argument type \"%s\"\n"
                          "iRow must be an integer column vector",
                          type_name().c_str(),
                          args(0).type_name().c_str());
                return octave_value();
        }

        const int32NDArray iRow(args(0).int32_array_value());

        ColumnVector X(iRow.numel());

        const int32_t iSize = pX->iGetSize();

        for (int i = 0; i < iRow.numel(); ++i) {
                if (int32_t(iRow(i)) < 1 || int32_t(iRow(i)) > iSize) {
                        error("%s: index %d out of range [%d-%d]",
                                        type_name().c_str(),
                                        int32_t(iRow(i)),
                                        1,
                                        iSize);

                        return octave_value();
                }

                X(i) = pX->dGetCoef(iRow(i));
        }

        return octave_value(X);
}

METHOD_DEFINE(ConstVectorHandlerInterface, iGetSize, args, nargout)
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if ( args.length() != 0 ) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<integer>(pX->iGetSize()));
}

BEGIN_METHOD_TABLE(ConstVectorHandlerInterface, MBDynInterface)
        METHOD_DISPATCH(ConstVectorHandlerInterface, dGetCoef)
        METHOD_DISPATCH(ConstVectorHandlerInterface, GetVec)
        METHOD_DISPATCH(ConstVectorHandlerInterface, iGetSize)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(ConstVectorHandlerInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(ConstVectorHandlerInterface, "ConstVectorHandler", "ConstVectorHandler");

VectorHandlerInterface::VectorHandlerInterface(OctaveInterface* pInterface, VectorHandler* X)
        :ConstVectorHandlerInterface(pInterface, X)
{
        TRACE("constructor");
}

VectorHandlerInterface::~VectorHandlerInterface()
{
        TRACE("destructor");
}

METHOD_DEFINE(VectorHandlerInterface, PutCoef, args, nargout)
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if ( args.length() != 2 ) {
                error("%s: invalid number of arguments %ld\n"
                                "two arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        if (!(args(0).is_scalar_type() && args(0).isinteger())) {
                error("%s: invalid argument type (%s)\n"
                                "iRow must be an integer",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        const integer iRow = args(0).int32_scalar_value();

        if ( iRow < 1 || iRow > pX->iGetSize() ) {
                error("%s: index %ld out of range [%ld-%ld]",
                                type_name().c_str(),
                                long(iRow),
                                1L,
                                long(pX->iGetSize()));
                return octave_value();
        }

        if ( !args(1).is_scalar_type() ) {
                error("%s: invalid argument type \"%s\"\n"
                                "dCoef must be a scalar value",
                                type_name().c_str(),
                                args(1).type_name().c_str());
                return octave_value();
        }

        const doublereal dCoef = args(1).scalar_value();

        pX->PutCoef(iRow, dCoef);

        return octave_value();
}

METHOD_DEFINE(VectorHandlerInterface, PutVec, args, nargout)
{
        if ( pX == 0 ) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if ( args.length() != 2 ) {
                error("%s: invalid number of arguments %ld\n"
                                "two arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        if (!(args(0).isinteger() && args(0).columns() == 1)) {
                error("%s: invalid argument type %ldx%ld (%s)\n"
                                "iRow must be an integer column vector",
                                type_name().c_str(),
                                long(args(0).rows()),
                                long(args(0).columns()),
                                args(0).type_name().c_str());

                return octave_value();
        }

        const int32NDArray iRow(args(0).int32_array_value());

        if (!(args(1).is_real_matrix()
                        && args(1).columns() == 1
                        && args(1).rows() == iRow.numel())) {
                error("%s: invalid argument type %ldx%ld (%s)\n"
                                "X must be a real column vector "
                                "with the same number of rows like iRow",
                                type_name().c_str(),
                                long(args(1).rows()),
                                long(args(1).columns()),
                                args(1).type_name().c_str());

                return octave_value();
        }

        const ColumnVector X(args(1).column_vector_value());

        ASSERT(X.numel() == iRow.numel());

        const int32_t iSize = pX->iGetSize();

        for (int i = 0; i < iRow.numel(); ++i) {
                if (int32_t(iRow(i)) < 1 || int32_t(iRow(i)) > iSize) {
                        error("%s: index %d out of range [%d-%d]",
                                        type_name().c_str(),
                                        int32_t(iRow(i)),
                                        1,
                                        iSize);
                        return octave_value();
                }

                pX->PutCoef(iRow(i), X(i));
        }

        return octave_value();
}

BEGIN_METHOD_TABLE(VectorHandlerInterface, ConstVectorHandlerInterface)
        METHOD_DISPATCH(VectorHandlerInterface, PutCoef)
        METHOD_DISPATCH(VectorHandlerInterface, PutVec)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(VectorHandlerInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(VectorHandlerInterface, "VectorHandler", "VectorHandler");

const std::string OStreamInterface::strsprintf("sprintf");

OStreamInterface::OStreamInterface(OctaveInterface* pInterface, std::ostream* pOS)
: MBDynInterface(pInterface),
  pOS(pOS)
{

}

OStreamInterface::~OStreamInterface()
{

}

METHOD_DEFINE(OStreamInterface, printf, args, nargout)
{
        if (!pOS) {
                error("ostream: not connected");
                return octave_value();
        }

#if OCTAVE_MAJOR_VERSION >= 5
        const octave_value_list ans = octave::feval(strsprintf, args, 1);
#else
        const octave_value_list ans = feval(strsprintf, args, 1);
#endif

#if OCTAVE_MAJOR_VERSION < 6
        if (error_state) {
                return octave_value();
        }
#endif
        
        if (!(ans.length() >= 1 && ans(0).is_string())) { // sprintf returns more than one output argument
                error("ostream: %s failed", strsprintf.c_str());
                return octave_value();
        }

        const std::string str(ans(0).string_value());

        for (std::string::const_iterator p = str.begin(); p != str.end(); ++p) {
                *pOS << *p;

                if (*p == '\n') // conform to default behavior in C and octave
                        pOS->flush();
        }

        return octave_value(octave_int<size_t>(str.length()));
}

BEGIN_METHOD_TABLE(OStreamInterface, MBDynInterface)
        METHOD_DISPATCH(OStreamInterface, printf)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(OStreamInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(OStreamInterface, "ostream", "ostream");


SimulationEntityInterface::SimulationEntityInterface(OctaveInterface* pInterface)
: MBDynInterface(pInterface)
{

}

SimulationEntityInterface::~SimulationEntityInterface()
{

}

METHOD_DEFINE(SimulationEntityInterface, iGetNumPrivData, args, nargout)
{
        const SimulationEntity* const pSimulationEntity = Get();

        if (!pSimulationEntity) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments\n"
                                "expected no arguments", type_name().c_str());

                return octave_value();
        }

        return octave_value(octave_int<unsigned int>(pSimulationEntity->iGetNumPrivData()));
}

METHOD_DEFINE(SimulationEntityInterface, iGetPrivDataIdx, args, nargout)
{
        const SimulationEntity* pSimulationEntity = Get();

        if (!pSimulationEntity) {
                error("%s: not connected",
                                type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments\n"
                                "expected one argument",
                                type_name().c_str());

                return octave_value();
        }

        if ( !args(0).is_string() ) {
                error("%s: invalid argument type (%s)\n"
                                "expected string",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        const std::string name = args(0).string_value();

        return octave_value(octave_int<unsigned int>(pSimulationEntity->iGetPrivDataIdx(name.c_str())));
}

METHOD_DEFINE(SimulationEntityInterface, dGetPrivData, args, nargout)
{
        const SimulationEntity* pSimulationEntity = Get();

        if (!pSimulationEntity) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments\n"
                                "expected one arguments", type_name().c_str());

                return octave_value();
        }

        if (!(args(0).isinteger() && args(0).is_scalar_type())) {
                error("%s: invalid argument type (%s)\n"
                                "expected integer scalar",
                                type_name().c_str(),
                                args(0).type_name().c_str());
                return octave_value();
        }

        const int32_t idx = args(0).int32_scalar_value();

        return octave_value(pSimulationEntity->dGetPrivData(idx));
}

BEGIN_METHOD_TABLE(SimulationEntityInterface, MBDynInterface)
        METHOD_DISPATCH(SimulationEntityInterface, iGetNumPrivData)
        METHOD_DISPATCH(SimulationEntityInterface, iGetPrivDataIdx)
        METHOD_DISPATCH(SimulationEntityInterface, dGetPrivData)
END_METHOD_TABLE()

NodeInterface::NodeInterface(OctaveInterface* pInterface)
: SimulationEntityInterface(pInterface)
{

}

NodeInterface::~NodeInterface()
{

}

void NodeInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
        const Node* const pNode = Get();

        if (!pNode) {
                error("%s: not connected", type_name().c_str());
                return;
        }

        os << type_name() << "(" << pNode->GetLabel() << "):" << pNode->GetName() << std::endl;
}

METHOD_DEFINE(NodeInterface, GetLabel, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(octave_int<unsigned int>(pNode->GetLabel()));
}

METHOD_DEFINE(NodeInterface, GetName, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(pNode->GetName());
}

METHOD_DEFINE(NodeInterface, iGetFirstIndex, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->iGetFirstIndex()));
}

METHOD_DEFINE(NodeInterface, iGetFirstRowIndex, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->iGetFirstRowIndex()));
}

METHOD_DEFINE(NodeInterface, iGetFirstColIndex, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->iGetFirstColIndex()));
}

METHOD_DEFINE(NodeInterface, dGetDofValue, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 2) {
                error("%s: invalid number of arguments %ld\n"
                                "%ld arguments expected",
                                type_name().c_str(),
                                long(args.length()),
                                2L);

                return octave_value();
        }

        if (!(args(0).isinteger() && args(0).is_scalar_type())) {
                error("%s: invalid argument type (%s) for argument 1\n"
                                "integer scalar expected",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        const int32_t iDof = args(0).int32_scalar_value();

        if (!(args(1).isinteger() && args(1).is_scalar_type())) {
                error("%s: invalid argument type (%s) for argument 2\n"
                                "integer scalar expected",
                                type_name().c_str(),
                                args(1).type_name().c_str());

                return octave_value();
        }

        const int32_t iOrder = args(1).int32_scalar_value();

        return octave_value(pNode->dGetDofValue(iDof, iOrder));
}

METHOD_DEFINE(NodeInterface, dGetDofValuePrev, args, nargout)
{
        const Node* const pNode = Get();

        if (pNode == 0) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 2) {
                error("%s: invalid number of arguments %ld\n"
                                "%ld arguments expected",
                                type_name().c_str(),
                                long(args.length()),
                                2L);

                return octave_value();
        }

        if (!(args(0).isinteger() && args(0).is_scalar_type())) {
                error("%s: invalid argument type (%s) for argument 1\n"
                                "integer scalar expected",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        const int32_t iDof = args(0).int32_scalar_value();

        if (!(args(1).isinteger() && args(1).is_scalar_type())) {
                error("%s: invalid argument type (%s) for argument 2\n"
                                "integer scalar expected",
                                type_name().c_str(),
                                args(1).type_name().c_str());

                return octave_value();
        }

        const int32_t iOrder = args(1).int32_scalar_value();

        return octave_value(pNode->dGetDofValuePrev(iDof, iOrder));
}

BEGIN_METHOD_TABLE(NodeInterface, SimulationEntityInterface)
        METHOD_DISPATCH(NodeInterface, GetLabel)
        METHOD_DISPATCH(NodeInterface, GetName)
        METHOD_DISPATCH(NodeInterface, iGetFirstIndex)
        METHOD_DISPATCH(NodeInterface, iGetFirstRowIndex)
        METHOD_DISPATCH(NodeInterface, iGetFirstColIndex)
        METHOD_DISPATCH(NodeInterface, dGetDofValue)
        METHOD_DISPATCH(NodeInterface, dGetDofValuePrev)
END_METHOD_TABLE()

ScalarNodeInterface::ScalarNodeInterface(OctaveInterface* pInterface, ScalarNode* pNode)
: NodeInterface(pInterface),
pNode(pNode)
{

}

ScalarNodeInterface::~ScalarNodeInterface()
{

}

const ScalarNode* ScalarNodeInterface::Get()const
{
        return pNode;
}

METHOD_DEFINE(ScalarNodeInterface, SetX, args, nargout)
{
        if (!pNode) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected one argument",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        if (!args(0).is_real_scalar()) {
                error("%s: invalid argument type (%s)\n"
                                "expected real scalar",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        pNode->SetX(args(0).scalar_value());

        return octave_value();
}

METHOD_DEFINE(ScalarNodeInterface, dGetX, args, nargout)
{
        if (!pNode) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "expected no argument",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(pNode->dGetX());
}

METHOD_DEFINE(ScalarNodeInterface, SetXPrime, args, nargout)
{
        if (!pNode) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected one argument",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        if (!args(0).is_real_scalar()) {
                error("%s: invalid argument type (%s)\n"
                                "expected real scalar",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        pNode->SetXPrime(args(0).scalar_value());

        return octave_value();
}

METHOD_DEFINE(ScalarNodeInterface, dGetXPrime, args, nargout)
{
        if (!pNode) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "expected no argument",
                                type_name().c_str(),
                                long(args.length()));

                return octave_value();
        }

        return octave_value(pNode->dGetXPrime());
}

BEGIN_METHOD_TABLE(ScalarNodeInterface, NodeInterface)
        METHOD_DISPATCH(ScalarNodeInterface, SetX)
        METHOD_DISPATCH(ScalarNodeInterface, dGetX)
        METHOD_DISPATCH(ScalarNodeInterface, SetXPrime)
        METHOD_DISPATCH(ScalarNodeInterface, dGetXPrime)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(ScalarNodeInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(ScalarNodeInterface, "ScalarNode", "ScalarNode");

StructDispNodeBaseInterface::StructDispNodeBaseInterface(OctaveInterface* pInterface)
: NodeInterface(pInterface)
{
        TRACE("constructor");
}

StructDispNodeBaseInterface::~StructDispNodeBaseInterface()
{
        TRACE("destructor");
}

octave_value
StructDispNodeBaseInterface::GetVec3(const Vec3& V, const octave_value_list& args) const
{
        octave_value octValue;

        if (args.length() != 0) {
                error("StructNode: invalid number of arguments %ld"
                                "no arguments expected", long(args.length()));
                return octValue;
        }

        if (!GetInterface()->ConvertMBDynToOctave(V, octValue)) {
                error("StructNode: could not convert data");
        }

        return octValue;
}

octave_value
StructDispNodeBaseInterface::GetMat3x3(const Mat3x3& M, const octave_value_list& args) const
{
        octave_value octValue;

        if (args.length() != 0) {
                error("StructNode: invalid number of arguments %ld"
                                "no arguments expected", long(args.length()));
                return octValue;
        }

        if (!GetInterface()->ConvertMBDynToOctave(M, octValue)) {
                error("StructNode: could not convert data");
        }

        return octValue;
}

METHOD_DEFINE(StructDispNodeBaseInterface, iGetFirstPositionIndex, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (args.length() != 0) {
                error("StructNode: invalid number of arguments %ld"
                                "no arguments expected", long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->iGetFirstPositionIndex()));
}

METHOD_DEFINE(StructDispNodeBaseInterface, iGetFirstMomentumIndex, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (args.length() != 0) {
                error("StructNode: invalid number of arguments %ld"
                                "no arguments expected", long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->iGetFirstMomentumIndex()));
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetLabel, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (args.length() != 0) {
                error("StructNode: invalid number of arguments %ld\n"
                                "no arguments expected", long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<integer>(pNode->GetLabel()));
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetXCurr, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetXCurr(), args);
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetXPrev, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetXPrev(), args);
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetVCurr, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetVCurr(), args);
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetVPrev, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetVPrev(), args);
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetXPPCurr, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (!pNode->bComputeAccelerations()) {
                error("StructNode: accelerations are not available for node %u", pNode->GetLabel());
                return octave_value();
        }

        return GetVec3(pNode->GetXPPCurr(), args);
}

METHOD_DEFINE(StructDispNodeBaseInterface, GetXPPPrev, args, nargout)
{
        const StructDispNode* const pNode = Get();

        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (!pNode->bComputeAccelerations()) {
                error("StructNode: accelerations are not available for node %u", pNode->GetLabel());
                return octave_value();
        }

        return GetVec3(pNode->GetXPPPrev(), args);
}

BEGIN_METHOD_TABLE(StructDispNodeBaseInterface, NodeInterface)
        METHOD_DISPATCH(StructDispNodeBaseInterface, iGetFirstIndex)
        METHOD_DISPATCH(StructDispNodeBaseInterface, iGetFirstPositionIndex)
        METHOD_DISPATCH(StructDispNodeBaseInterface, iGetFirstMomentumIndex)
        METHOD_DISPATCH(StructDispNodeBaseInterface, GetLabel)
        METHOD_DISPATCH(StructDispNodeBaseInterface, GetXCurr)
        METHOD_DISPATCH(StructDispNodeBaseInterface, GetXPrev)
        METHOD_DISPATCH(StructNodeInterface, GetVCurr)
        METHOD_DISPATCH(StructNodeInterface, GetVPrev)
        METHOD_DISPATCH(StructNodeInterface, GetXPPCurr)
        METHOD_DISPATCH(StructNodeInterface, GetXPPPrev)
END_METHOD_TABLE()

StructDispNodeInterface::StructDispNodeInterface(OctaveInterface* pInterface, const StructDispNode* pNode)
: StructDispNodeBaseInterface(pInterface),
  pNode(pNode)
{

}

StructDispNodeInterface::~StructDispNodeInterface()
{

}

const StructDispNode* StructDispNodeInterface::Get()const
{
        return pNode;
}

DEFINE_OCTAVE_ALLOCATOR(StructDispNodeInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(StructDispNodeInterface, "StructDispNode", "StructDispNode");

StructNodeInterface::StructNodeInterface(OctaveInterface* pInterface, const StructNode* pNode)
: StructDispNodeBaseInterface(pInterface),
pNode(pNode)
{
        TRACE("constructor");
}

StructNodeInterface::~StructNodeInterface(void)
{
        TRACE("destructor");
}

const StructNode*
StructNodeInterface::Get()const
{
        return pNode;
}

METHOD_DEFINE(StructNodeInterface, GetgCurr, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetgCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgRef, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetgRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgPCurr, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetgPCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetgPRef, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetgPRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRCurr, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetMat3x3(pNode->GetRCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRPrev, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetMat3x3(pNode->GetRPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetRRef, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetMat3x3(pNode->GetRRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWCurr, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetWCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPrev, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetWPrev(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWRef, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        return GetVec3(pNode->GetWRef(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPCurr, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (!pNode->bComputeAccelerations()) {
                error("StructNode: accelerations are not available for node %u", pNode->GetLabel());
                return octave_value();
        }

        return GetVec3(pNode->GetWPCurr(), args);
}

METHOD_DEFINE(StructNodeInterface, GetWPPrev, args, nargout)
{
        if (pNode == 0) {
                error("StructNode: not connected");
                return octave_value();
        }

        if (!pNode->bComputeAccelerations()) {
                error("StructNode: accelerations are not available for node %u", pNode->GetLabel());
                return octave_value();
        }

        return GetVec3(pNode->GetWPPrev(), args);
}

BEGIN_METHOD_TABLE(StructNodeInterface, StructDispNodeBaseInterface)
        METHOD_DISPATCH(StructNodeInterface, GetgCurr)
        METHOD_DISPATCH(StructNodeInterface, GetgRef)
        METHOD_DISPATCH(StructNodeInterface, GetgPCurr)
        METHOD_DISPATCH(StructNodeInterface, GetgPRef)
        METHOD_DISPATCH(StructNodeInterface, GetRCurr)
        METHOD_DISPATCH(StructNodeInterface, GetRPrev)
        METHOD_DISPATCH(StructNodeInterface, GetRRef)
        METHOD_DISPATCH(StructNodeInterface, GetWCurr)
        METHOD_DISPATCH(StructNodeInterface, GetWPrev)
        METHOD_DISPATCH(StructNodeInterface, GetWRef)
        METHOD_DISPATCH(StructNodeInterface, GetWPCurr)
        METHOD_DISPATCH(StructNodeInterface, GetWPPrev)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(StructNodeInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(StructNodeInterface, "StructNode", "StructNode");

DataManagerInterface::DataManagerInterface(OctaveInterface* pInterface, const DataManager* pDM, bool bAddRef)
        :MBDynInterface(pInterface, bAddRef),
         pDM(pDM)
{
        TRACE("constructor");

        ASSERT(pInterface != 0);
        ASSERT(pDM != 0);
}

DataManagerInterface::~DataManagerInterface()
{
        TRACE("destructor");
}

void DataManagerInterface::print(std::ostream& os, bool pr_as_read_syntax)const
{
        os << "DataManager" << std::endl
                        << "iTotDofs=" << GetDataManager()->iGetNumDofs() << std::endl;
}

METHOD_DEFINE(DataManagerInterface, GetVariable, args, nargout)
{
        octave_value octValue;

        if (args.length() != 1) {
                error("DataManager: invalid number of arguments %ld\n"
                                "one argument expected", long(args.length()));
                return octValue;
        }

        if (!args(0).is_string()) {
                error("DataManager: invalid argument type \"%s\"\n"
                                "varName must be a string", args(0).type_name().c_str());
                return octValue;
        }

        const std::string varName = args(0).string_value();

        const NamedValue* const pVar = GetSymbolTable().Get(varName.c_str());

        if (pVar == 0) {
                error("DataManager: variable \"%s\" not found in symbol table", varName.c_str());
                return octValue;
        }

        const TypedValue mbValue = pVar->GetVal();

        if (!GetInterface()->ConvertMBDynToOctave(mbValue, octValue)) {
                error("%s: failed to convert MBDyn variable %s"
                                " of type %s value into octave value",
                                type_name().c_str(),
                                varName.c_str(),
                                mbValue.GetTypeName());
        }

        return octValue;
}

METHOD_DEFINE(DataManagerInterface, GetStructNodePos, args, nargout)
{
        octave_value octValue;

        if (args.length() != 1) {
                error("DataManager: invalid number of arguments\n"
                                "one argument expected");
                return octValue;
        }

        if (!(args(0).is_scalar_type() && args(0).isinteger())) {
                error("DataManager: invalid argument type \"%s\"\n"
                                "iNodeLabel must be an integer", args(0).type_name().c_str());
                return octValue;
        }

        integer iNodeLabel = args(0).int32_scalar_value();

        Node* pNode = GetDataManager()->pFindNode(Node::STRUCTURAL, iNodeLabel);

        if (pNode == 0) {
                error("DataManager: invalid node label\n"
                                "could not find node %ld", long(iNodeLabel));
                return octValue;
        }

        StructDispNode* pStructNode = dynamic_cast<StructDispNode*>(pNode);

        if (pStructNode == 0) {
                error("DataManager: invalid node type\n"
                                "node %ld is not a structural displacement node", long(iNodeLabel));
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

        if (args.length() != 1) {
                error("DataManager: invalid number of arguments\n"
                                "one argument expected");
                return octValue;
        }

        if (!(args(0).is_scalar_type() && args(0).isinteger())) {
                error("DataManager: invalid argument type \"%s\"\n"
                                "iNodeLabel must be an integer", args(0).type_name().c_str());
                return octValue;
        }

        integer iNodeLabel = static_cast<int32_t>(args(0).int32_scalar_value());

        Node* pNode = GetDataManager()->pFindNode(Node::STRUCTURAL, iNodeLabel);

        if (pNode == 0) {
                error("DataManager: could not find node %ld", long(iNodeLabel));
                return octValue;
        }

        StructNode* pStructNode = dynamic_cast<StructNode*>(pNode);

        if (pStructNode) {
                octValue = new StructNodeInterface(GetInterface(), pStructNode);
                return octValue;
        }

        StructDispNode* pStructDispNode = dynamic_cast<StructDispNode*>(pNode);

        if (!pStructDispNode) {
                error("%s: node %ld is not a structural node", type_name().c_str(), long(iNodeLabel));
                return octValue;
        }

        octValue = new StructDispNodeInterface(GetInterface(), pStructDispNode);

        return octValue;
}

METHOD_DEFINE(DataManagerInterface, pFindNode, args, nargout)
{
        Node::Type Type;

        Node* const pNode = pFindNode_(args, Type);

        if (!pNode) {
                return octave_value(); // error has been called
        }

        NodeInterface* const pNodeInterface = CreateNodeInterface(pNode, Type);

        if (pNodeInterface == 0) {
                error("%s: unknown node type %ld", type_name().c_str(), long(Type));
                return octave_value();
        }

        return octave_value(pNodeInterface);
}

METHOD_DEFINE(DataManagerInterface, ReadNode, args, nargout)
{
        if (args.length() != 2) {
                error("%s: invalid number of arguments (%ld)\n"
                                "arguments expected: %ld",
                                type_name().c_str(),
                                long(args.length()),
                                2L);

                return octave_value();
        }

        const MBDynParserInterface* const pHPI =
                        dynamic_cast<const MBDynParserInterface*>(&args(0).get_rep());

        if (pHPI == 0) {
                error("%s: invalid argument type (%s)\n"
                                "expected MBDynParser",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return octave_value();
        }

        MBDynParser* const pHP = pHPI->Get();

        if (pHP == 0) {
                error("%s: not connected", pHPI->type_name().c_str());
                return octave_value();
        }

        if (!args(1).is_string()) {
                error("%s: invalid argument type (%s)\n"
                                "expected string",
                                type_name().c_str(),
                                args(1).type_name().c_str());

                return octave_value();
        }

        const std::string strType(args(1).string_value());

        const Node::Type Type = GetNodeType(strType);

        if (Type == Node::UNKNOWN) {
                error("%s: invalid node type (%s)", type_name().c_str(), strType.c_str());
                return octave_value();
        }

        Node* pNode = 0;

        try {
                // FIXME: This function should not be called from a drive caller
                // FIXME: because DataManager is assumed to be constant in this case?
                pNode = const_cast<DataManager*>(pDM)->ReadNode(*pHP, Type);
        } catch(...) { // do not throw across call back boundaries
                error("%s: ReadNode failed", type_name().c_str());
                return octave_value();
        }

        NodeInterface* const pNodeInterface = CreateNodeInterface(pNode, Type);

        if (pNodeInterface == 0) {
                error("%s: unknown node type %ld", type_name().c_str(), long(Type));
                return octave_value();
        }

        return octave_value(pNodeInterface);
}

METHOD_DEFINE(DataManagerInterface, dGetTime, args, nargout)
{
        if (args.length() != 0) {
                error("DataManager: invalid number of arguments\n"
                                "no argument expected");
                return octave_value();
        }

        return octave_value(GetDataManager()->dGetTime());
}

const DataManager* DataManagerInterface::GetDataManager()const
{
        ASSERT(GetInterface() != 0);
        ASSERT(pDM != 0);

        return pDM;
}

const Table& DataManagerInterface::GetSymbolTable()const
{
        return GetDataManager()->GetMathParser().GetSymbolTable();
}

Node* DataManagerInterface::pFindNode_(const octave_value_list& args, Node::Type& Type) const
{
        if (args.length() != 2) {
                error("%s: invalid number of arguments\n"
                                "arguments expected 2", type_name().c_str());

                return 0;
        }

        if (!(args(0).is_string())) {
                error("%s: invalid argument type (%s)\n"
                                "Type must be a string",
                                type_name().c_str(),
                                args(0).type_name().c_str());

                return 0;
        }

        const std::string strType(args(0).string_value());

        Type = GetNodeType(strType);

        if (Type == Node::UNKNOWN) {
                error("%s: unknown node type (%s)", type_name().c_str(), strType.c_str());
                return 0;
        }

        if (!(args(1).is_scalar_type() && args(1).isinteger())) {
                error("%s: invalid argument type (%s)\n"
                                "iNodeLabel must be an integer",
                                type_name().c_str(),
                                args(1).type_name().c_str());

                return 0;
        }

        const int32_t iNodeLabel = args(1).int32_scalar_value();

        Node* pNode = GetDataManager()->pFindNode(Type, iNodeLabel);

        if (pNode == 0) {
                error("%s: could not find node %ld", type_name().c_str(), long(iNodeLabel));
                return 0;
        }

        return pNode;
}

Node::Type DataManagerInterface::GetNodeType(const std::string& strType) const
{
        static const struct {
                char name[11];
                Node::Type type;
        }
        nodeTypes[] = {
                        { "ABSTRACT",   Node::ABSTRACT },
                        { "STRUCTURAL", Node::STRUCTURAL },
                        { "ELECTRIC",   Node::ELECTRIC },
                        { "THERMAL",    Node::THERMAL },
                        { "PARAMETER",	Node::PARAMETER },
                        { "HYDRAULIC",	Node::HYDRAULIC }
        };

        static const int count = sizeof(nodeTypes) / sizeof(nodeTypes[0]);

        for (int i = 0; i < count; ++i) {
                if (strType == nodeTypes[i].name) {
                        return nodeTypes[i].type;
                }
        }

        return Node::UNKNOWN;
}

NodeInterface* DataManagerInterface::CreateNodeInterface(Node* pNode, Node::Type Type) const
{
        switch (Type) {
                case Node::STRUCTURAL:
                {
                        StructNode* const pStructNode = dynamic_cast<StructNode*>(pNode);

                        if (pStructNode == 0) {
                                StructDispNode* pStructDispNode = dynamic_cast<StructDispNode*>(pNode);
                                if (pStructDispNode == 0) {
                                        ASSERT(0);
                                        return 0;
                                }

                                return new StructDispNodeInterface(GetInterface(), pStructDispNode);
                        }
                        return new StructNodeInterface(GetInterface(), pStructNode);
                }
                case Node::ABSTRACT:
                case Node::ELECTRIC:
                case Node::THERMAL:
                case Node::PARAMETER:
                case Node::HYDRAULIC:
                {
                        ScalarNode* const pScalarNode = dynamic_cast<ScalarNode*>(pNode);
                        ASSERT(pScalarNode != 0);
                        return new ScalarNodeInterface(GetInterface(), pScalarNode);
                }
                default:
                        ASSERT(0); // FIXME: add StructuralDisplacementNode
                        return 0;
        }
}

BEGIN_METHOD_TABLE(DataManagerInterface, MBDynInterface)
        METHOD_DISPATCH(DataManagerInterface, GetVariable)
        METHOD_DISPATCH(DataManagerInterface, GetStructNodePos)
        METHOD_DISPATCH(DataManagerInterface, GetStructNode)
        METHOD_DISPATCH(DataManagerInterface, pFindNode)
        METHOD_DISPATCH(DataManagerInterface, ReadNode)
        METHOD_DISPATCH(DataManagerInterface, dGetTime)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(DataManagerInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(DataManagerInterface, "DataManager", "DataManager");

const MBDynParserInterface::MBDynStringDelims
MBDynParserInterface::mbStringDelims[6] = {
                { HighParser::PLAINBRACKETS,	"PLAINBRACKETS" },
                { HighParser::SQUAREBRACKETS,	"SQUAREBRACKETS" },
                { HighParser::CURLYBRACKETS,	"CURLYBRACKETS" },
                { HighParser::SINGLEQUOTE,		"SINGLEQUOTE" },
                { HighParser::DOUBLEQUOTE,		"DOUBLEQUOTE" },
                { HighParser::DEFAULTDELIM,		"DEFAULTDELIM" }
};

MBDynParserInterface::MBDynParserInterface(OctaveInterface* pInterface, MBDynParser* pHP, bool bAddRef)
        :MBDynInterface(pInterface, bAddRef),
         pHP(pHP)
{
        if (!pHP) {
                if (pInterface) {
                        ASSERT(0); // should not happen
                        this->pHP = &pInterface->GetDataManager()->GetMBDynParser();
                }
        }
}

MBDynParserInterface::~MBDynParserInterface()
{

}

void MBDynParserInterface::print(std::ostream& os, bool pr_as_read_syntax) const
{
        os << "MBDynParser";

        if (pHP) {
                const IncludeParser::ErrOut lineData(pHP->GetLineData());
                os << " at " << lineData.sFileName << ":" << lineData.iLineNumber;
        }

        os << std::endl;
}

METHOD_DEFINE(MBDynParserInterface, IsKeyWord, args, nargout)
{
        if (!pHP) {
                error("MBDynParser: not connected");
                return octave_value();
        }

        if (args.length() != 1) {
                error("MBDynParser: invalid number of arguments %ld\n"
                                "expected %ld arguments", long(args.length()), 1L);
                return octave_value();
        }

        if (!args(0).is_string()) {
                error("MBDynParser: invalid argument type (%s)\n"
                                "expected string", args(0).type_name().c_str());
                return octave_value();
        }

        const std::string strKeyWord = args(0).string_value();

        bool bVal;

        try {
                bVal = pHP->IsKeyWord(strKeyWord.c_str());
        } catch(...) {
                error("%s: IsKeyWord failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(bVal);
}

METHOD_DEFINE(MBDynParserInterface, IsArg, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments %ld\n"
                                "no arguments expected",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        }

        return octave_value(pHP->IsArg());
}

METHOD_DEFINE(MBDynParserInterface, IsStringWithDelims, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() > 1) {
                error("%s: invalid number of arguments %ld\n"
                                "one argument expected",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        }

        HighParser::Delims delims = HighParser::DEFAULTDELIM;

        if (args.length() >= 1) {
                if (!args(0).is_string()) {
                        error("%s: invalid argument type %s\n"
                                        "expected string",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }

                const std::string strDelims(args(0).string_value());

                if (!GetDelimsFromString(strDelims, delims)) {
                        error("%s: invalid argument %s",
                                        type_name().c_str(),
                                        strDelims.c_str());
                        return octave_value();
                }
        }

        return octave_value(pHP->IsStringWithDelims(delims));
}

METHOD_DEFINE(MBDynParserInterface, GetReal, args, nargout)
{
        if (!pHP) {
                error("MBDynParser: not connected");
                return octave_value();
        }

        doublereal dDefVal = 0.0;

        if (args.length() > 1) {
                error("MBDynParser: invalid number of arguments %ld\n"
                                "expected 0-1 arguments", long(args.length()));
                return octave_value();
        } else if (args.length() >= 1) {
                if (!args(0).is_real_scalar()) {
                        error("MBDynParser: invalid argument type (%s)\n"
                                        "expected real scalar", args(0).type_name().c_str());
                        return octave_value();
                }
                dDefVal = args(0).scalar_value();
        }

        doublereal dVal;

        try {
                dVal = pHP->GetReal(dDefVal);
        } catch(...) {
                error("%s: GetReal failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(dVal);
}

METHOD_DEFINE(MBDynParserInterface, GetInt, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        integer iDefVal = 0;

        if (args.length() > 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected 0-1 arguments",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        } else if (args.length() >= 1) {
                if (!args(0).isinteger() && args(0).is_scalar_type()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected integer scalar",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }
                iDefVal = args(0).int_value();
        }

        integer iVal;

        try {
                iVal = pHP->GetInt(iDefVal);
        } catch(...) {
                error("%s: GetInt failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(iVal);
}

METHOD_DEFINE(MBDynParserInterface, GetBool, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        bool bDefVal = 0;

        if (args.length() > 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected 0-1 arguments",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        } else if (args.length() >= 1) {
                if (!args(0).is_bool_scalar()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected bool scalar",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }
                bDefVal = args(0).bool_value();
        }

        bool bVal;

        try {
                bVal = pHP->GetBool(bDefVal);
        } catch(...) {
                error("%s: GetBool failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(bVal);
}

METHOD_DEFINE(MBDynParserInterface, GetYesNoOrBool, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        bool bDefVal = 0;

        if (args.length() > 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected 0-1 arguments",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        } else if (args.length() >= 1) {
                if (!args(0).is_bool_scalar()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected bool scalar",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }
                bDefVal = args(0).bool_value();
        }

        bool bVal;

        try {
                bVal = pHP->GetYesNoOrBool(bDefVal);
        } catch(...) {
                error("%s: GetYesNoOrBool failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(bVal);
}

METHOD_DEFINE(MBDynParserInterface, GetString, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        std::string strDefVal;

        if (args.length() > 1) {
                error("%s: invalid number of arguments %ld\n"
                                "expected 0-1 arguments",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        } else if (args.length() >= 1) {
                if (!args(0).is_string()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected string",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }
                strDefVal = args(0).string_value();
        }

        std::string strVal;

        try {
                strVal = pHP->GetString(strDefVal);
        } catch(...) {
                error("%s: GetString failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(strVal);
}

METHOD_DEFINE(MBDynParserInterface, GetStringWithDelims, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        std::string strDelims;

        if (args.length() > 2) {
                error("%s: invalid number of arguments %ld\n"
                                "expected 0-2 arguments",
                                type_name().c_str(),
                                long(args.length()));
                return octave_value();
        }

        HighParser::Delims delims = HighParser::DEFAULTDELIM;

        if (args.length() >= 1) {
                if (!args(0).is_string()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected string",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }

                std::string strDelims(args(0).string_value());

                if (!GetDelimsFromString(strDelims, delims)) {
                        error("%s: invalid argument %s",
                                        type_name().c_str(),
                                        strDelims.c_str());
                        return octave_value();
                }
        }

        bool bEscape = true;

        if (args.length() >= 2) {
                if (!args(1).is_bool_scalar()) {
                        error("%s: invalid argument type (%s)\n"
                                        "expected bool",
                                        type_name().c_str(),
                                        args(1).type_name().c_str());
                        return octave_value();
                }

                bEscape = args(1).bool_value();
        }

        std::string strVal;

        try {
                strVal = pHP->GetStringWithDelims(delims, bEscape);
        } catch(...) {
                error("%s: GetStringWithDelims failed", type_name().c_str());
                return octave_value();
        }

        return octave_value(strVal);
}

METHOD_DEFINE(MBDynParserInterface, GetPosRel, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments\n"
                                "expected one argument", type_name().c_str());

                return octave_value();
        }

        const StructDispNode* const pNode = GetStructNode(args(0));

        if (pNode == 0) {
                return octave_value();
        }

        Vec3 mbX;

        try {
                mbX = pHP->GetPosRel(ReferenceFrame(pNode));
        } catch(...) { // do not throw across call back boundaries
                error("%s: GetPosRel failed", type_name().c_str());
                return octave_value();
        }

        octave_value octX;

        if (!GetInterface()->ConvertMBDynToOctave(mbX, octX)) {
                error("%s: could not convert data", type_name().c_str());
                return octave_value();
        }

        return octX;
}

METHOD_DEFINE(MBDynParserInterface, GetRotRel, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 1) {
                error("%s: invalid number of arguments\n"
                                "expected one argument", type_name().c_str());

                return octave_value();
        }

        const StructDispNode* const pNode = GetStructNode(args(0));

        if (pNode == 0) {
                return octave_value();
        }

        Mat3x3 mbR;

        try {
                mbR = pHP->GetRotRel(ReferenceFrame(pNode));
        } catch(...) { // do not throw across call back boundaries
                error("%s: GetRotRel failed", type_name().c_str());
                return octave_value();
        }

        octave_value octR;

        if (!GetInterface()->ConvertMBDynToOctave(mbR, octR)) {
                error("%s: could not convert data", type_name().c_str());
                return octave_value();
        }

        return octR;
}

METHOD_DEFINE(MBDynParserInterface, GetValue, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() > 1) {
                error("%s: invalid number of arguments\n"
                                "expected 0-1 arguments", type_name().c_str());
                return octave_value();
        }

        TypedValue mbDefValue;

        if (args.length() >= 1) {
                if (!GetInterface()->ConvertOctaveToMBDyn(args(0), mbDefValue)) {
                        error("%s: could not convert octave data type %s to MBDyn",
                                        type_name().c_str(),
                                        args(0).type_name().c_str());
                        return octave_value();
                }
        }

        TypedValue mbValue;

        try {
                mbValue = pHP->GetValue(mbDefValue);
        } catch (...) {
                error("%s: GetValue failed", type_name().c_str());
                return octave_value();
        }

        octave_value octValue;

        if (!GetInterface()->ConvertMBDynToOctave(mbValue, octValue)) {
                error("%s: could not convert data type %s into octave value",
                                type_name().c_str(),
                                mbValue.GetTypeName());
        }

        return octValue;
}

METHOD_DEFINE(MBDynParserInterface, GetLineData, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments\n"
                                "expected no argument", type_name().c_str());

                return octave_value();
        }

        std::ostringstream os;

        os << pHP->GetLineData();

        return octave_value(os.str());
}

const StructDispNode* MBDynParserInterface::GetStructNode(const octave_value& arg) const
{
        const StructDispNodeBaseInterface* const pNodeInterface =
                        dynamic_cast<const StructDispNodeBaseInterface*>(&arg.get_rep());

        if (pNodeInterface == 0) {
                error("%s: invalid argument type (%s)\n"
                                "expected StructNode or StructDispNode",
                                type_name().c_str(),
                                arg.type_name().c_str());

                return 0;
        }

        const StructDispNode* const pNode = pNodeInterface->Get();

        if (pNode == 0) {
                error("%s: not connected", pNodeInterface->type_name().c_str());
                return 0;
        }

        return pNode;
}

bool MBDynParserInterface::GetDelimsFromString(const std::string& strDelims, HighParser::Delims& delims)
{
        const int count = sizeof(mbStringDelims) / sizeof(mbStringDelims[0]);

        for (int i = 0; i < count; ++i) {
                if (strDelims == mbStringDelims[i].name) {
                        delims = mbStringDelims[i].value;
                        return true;
                }
        }

        delims = HighParser::DEFAULTDELIM;

        return false;
}

METHOD_DEFINE(MBDynParserInterface, GetDriveCaller, args, nargout)
{
        if (!pHP) {
                error("%s: not connected", type_name().c_str());
                return octave_value();
        }

        if (args.length() != 0) {
                error("%s: invalid number of arguments\n"
                                "expected no argument", type_name().c_str());

                return octave_value();
        }

        DriveCaller* pDC = 0;

        try {
                pDC = pHP->GetDriveCaller();
        } catch (...) {
                error("%s: GetDriveCaller failed", type_name().c_str());
                return octave_value();
        }

        DriveCallerInterface* pInterface = new DriveCallerInterface(GetInterface(), pDC);

        return octave_value(pInterface);
}

BEGIN_METHOD_TABLE(MBDynParserInterface, MBDynInterface)
        METHOD_DISPATCH(MBDynParserInterface, IsKeyWord)
        METHOD_DISPATCH(MBDynParserInterface, IsArg)
        METHOD_DISPATCH(MBDynParserInterface, IsStringWithDelims)
        METHOD_DISPATCH(MBDynParserInterface, GetReal)
        METHOD_DISPATCH(MBDynParserInterface, GetInt)
        METHOD_DISPATCH(MBDynParserInterface, GetBool)
        METHOD_DISPATCH(MBDynParserInterface, GetYesNoOrBool)
        METHOD_DISPATCH(MBDynParserInterface, GetString)
        METHOD_DISPATCH(MBDynParserInterface, GetStringWithDelims)
        METHOD_DISPATCH(MBDynParserInterface, GetValue)
        METHOD_DISPATCH(MBDynParserInterface, GetPosRel)
        METHOD_DISPATCH(MBDynParserInterface, GetRotRel)
        METHOD_DISPATCH(MBDynParserInterface, GetLineData)
        METHOD_DISPATCH(MBDynParserInterface, GetDriveCaller)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(MBDynParserInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(MBDynParserInterface, "MBDynParser", "MBDynParser");

DriveCallerInterface::DriveCallerInterface(OctaveInterface* pInterface, const DriveCaller* pDC)
:MBDynInterface(pInterface), DC(pDC)
{

}

DriveCallerInterface::~DriveCallerInterface()
{

}


METHOD_DEFINE(DriveCallerInterface, dGet, args, nargout)
{
        octave_value y;

        if (!DC.pGetDriveCaller()) {
                error("%s: not connected", type_name().c_str());
                return y;
        }

        switch (args.length())
        {
        case 0:
                y = DC.dGet();
                break;

        case 1:
        {
                doublereal x;

                if (!GetInterface()->ConvertOctaveToMBDyn(args(0), x)) {
                        error("%s: invalid argument type %s", type_name().c_str(), args(0).type_name().c_str());
                        break;
                }

                y = DC.dGet(x);
                break;
        }
        default:
                error("%s: invalid number of arguments\n"
                                "expected zero or one argument", type_name().c_str());
        };

        return y;
}

METHOD_DEFINE(DriveCallerInterface, dGetP, args, nargout)
{
        octave_value y;

        if (!DC.pGetDriveCaller()) {
                error("%s: not connected", type_name().c_str());
                return y;
        }

        if (!DC.pGetDriveCaller()->bIsDifferentiable()) {
                error("%s: drive caller(%d %s) is not differentiable", type_name().c_str(), DC.pGetDriveCaller()->GetLabel(), DC.pGetDriveCaller()->GetName().c_str());
                return y;
        }

        switch (args.length())
        {
        case 0:
                y = DC.dGetP();
                break;

        case 1:
        {
                doublereal x;

                if (!GetInterface()->ConvertOctaveToMBDyn(args(0), x)) {
                        error("%s: invalid argument type %s", type_name().c_str(), args(0).type_name().c_str());
                        break;
                }

                y = DC.dGetP(x);
                break;
        }
        default:
                error("%s: invalid number of arguments\n"
                                "expected zero or one argument", type_name().c_str());
        };

        return y;
}

BEGIN_METHOD_TABLE(DriveCallerInterface, MBDynInterface)
        METHOD_DISPATCH(DriveCallerInterface, dGet)
        METHOD_DISPATCH(DriveCallerInterface, dGetP)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(DriveCallerInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(DriveCallerInterface, "DriveCaller", "DriveCaller");

OctaveDriveCaller::OctaveDriveCaller(const std::string& strFunc, OctaveInterface* pInterface, int iFlags, const octave_value_list& args)
: DriveCaller(0),
  iFlags(iFlags),
  strFunc(strFunc),
  pInterface(pInterface),
  args(args)
{
        TRACE("constructor");
        TRACE("strFunc=" << strFunc);
        TRACE("iFlags=" <<  iFlags);
        
        if (pInterface) {
                pInterface->AddRef();
        }
}

OctaveDriveCaller::~OctaveDriveCaller(void)
{
        TRACE("destructor");

        if (pInterface) {
                pInterface->Destroy();
        }
}

DriveCaller *
OctaveDriveCaller::pCopy(void) const
{
        return new OctaveDriveCaller(strFunc, pInterface, iFlags, args);
}

std::ostream&
OctaveDriveCaller::Restart(std::ostream& out) const
{
        return out << "octave, \"" << strFunc << "\"";
}

inline doublereal
OctaveDriveCaller::dGet(const doublereal& dVar) const
{
        return pInterface->EvalScalarFunction(strFunc, MakeArgList(dVar), iFlags);
}

inline doublereal
OctaveDriveCaller::dGet(void) const
{
        return dGet(pInterface->GetDataManager()->dGetTime());
}

inline bool
OctaveDriveCaller::bIsDifferentiable(void) const
{
        return pInterface->HaveADPackage();
}

inline doublereal
OctaveDriveCaller::dGetP(const doublereal& dVar) const
{
        return pInterface->EvalScalarFunctionDerivative(strFunc, MakeArgList(dVar), iFlags);
}

inline doublereal
OctaveDriveCaller::dGetP(void) const
{
        if (!bIsDifferentiable()) {
                return 0.;
        }

        return dGetP(pInterface->GetDataManager()->dGetTime());
}

octave_value_list OctaveDriveCaller::MakeArgList(doublereal dVar) const
{
        return pInterface->MakeArgList(dVar, args, iFlags);
}

template <class T>
OctaveTplDriveCaller<T>::OctaveTplDriveCaller(const std::string& strFunction, OctaveInterface* pInterface, int iFlags, const octave_value_list& args)
        :strFunction(strFunction),
         pInterface(pInterface),
         iFlags(iFlags),
         args(args)
{
        TRACE("constructor");

        if (pInterface) {
                pInterface->AddRef();
        }
};

template <class T>
OctaveTplDriveCaller<T>::~OctaveTplDriveCaller(void)
{
        TRACE("destructor");
        
        if (pInterface) {
                pInterface->Destroy();
        }
};


template <class T>
TplDriveCaller<T>* OctaveTplDriveCaller<T>::pCopy(void) const
{
        return new OctaveTplDriveCaller(strFunction, pInterface, iFlags, args);
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
        pInterface->EvalMatrixFunction(strFunction, X, MakeArgList(dVar), iFlags);
        return X;
};

template <class T>
T OctaveTplDriveCaller<T>::Get(void) const
{
        return Get(pInterface->GetDataManager()->dGetTime());
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
        pInterface->EvalMatrixFunctionDerivative(strFunction, XP, MakeArgList(t), iFlags);
        return XP;
};

template <class T>
int OctaveTplDriveCaller<T>::getNDrives(void) const
{
        return 0;
};

template <class T>
octave_value_list OctaveTplDriveCaller<T>::MakeArgList(doublereal dVar) const
{
        return pInterface->MakeArgList(dVar, args, iFlags);
}

template <class T>
TplDriveCaller<T> *
OctaveTDCR<T>::Read(const DataManager* pDM, MBDynParser& HP)
{
        OctaveFunctionDCR::Read(pDM, HP);

        return new OctaveTplDriveCaller<T>(GetFunction(), GetInterface(), GetFlags(), GetArgs());
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

OctaveScalarFunction::OctaveScalarFunction(const std::string& strFunc, OctaveInterface* pInterface, int iFlags, const octave_value_list& args)
        :strFunc(strFunc),
         pInterface(pInterface),
         iFlags(iFlags),
         args(args)
{

}

OctaveScalarFunction::~OctaveScalarFunction(void)
{

}

doublereal OctaveScalarFunction::operator()(const doublereal x) const
{
        return pInterface->EvalScalarFunction(strFunc, MakeArgList(x), iFlags);
}

doublereal OctaveScalarFunction::ComputeDiff(const doublereal t, const integer order) const
{
        if (order != 1) {
                silent_cerr("octave scalar function \"" << strFunc << "\" derivative of order " << order << " not supported" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return pInterface->EvalScalarFunctionDerivative(strFunc, MakeArgList(t), iFlags);
}

octave_value_list OctaveScalarFunction::MakeArgList(doublereal dVar) const
{
        return pInterface->MakeArgList(dVar, args, iFlags);
}

OctaveConstitutiveLawBase::OctaveConstitutiveLawBase(const std::string& strClass, OctaveInterface* pInterface, int iFlags)
                : strClass(strClass),
                  pInterface(pInterface),
                  iFlags(iFlags)
{
        pInterface->AddRef();
}

OctaveConstitutiveLawBase::~OctaveConstitutiveLawBase()
{
        pInterface->Destroy();
}

bool OctaveConstitutiveLawBase::bHaveMethod(const std::string& strName) const
{
        return GetInterface()->bHaveMethod(octObject, strClass, strName);
}

const std::string OctaveConstitutiveLawBase::strGetConstLawType("GetConstLawType");
const std::string OctaveConstitutiveLawBase::strUpdate("Update");

template <class T, class Tder>
OctaveConstitutiveLaw<T, Tder>::OctaveConstitutiveLaw(const std::string& strClass, OctaveInterface* pInterface, int iFlags)
: OctaveConstitutiveLawBase(strClass, pInterface, iFlags),
  clType(ConstLawType::UNKNOWN)
{
        octave_value_list args(GetInterface()->GetDataManagerInterface());
        args.append(GetInterface()->GetMBDynParserInterface());

        const int flags = GetFlags() | OctaveInterface::UPDATE_OCTAVE_VARIABLES;

        octave_value_list ans = GetInterface()->EvalFunction(GetClass(), args, 1, flags);

        ASSERT(ans.length() == 1);

        octObject = ans(0);

        if (!octObject.isobject()) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel() << "): result of constructor ("
                                << ans(0).type_name() << ") is not an object at line "
                                << GetInterface()->GetMBDynParser()->GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!bHaveMethod(strGetConstLawType)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): method " << GetClass() << "." << strGetConstLawType
                                << " is undefined" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!bHaveMethod(strUpdate)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): method " << GetClass() << "." << strUpdate
                                << " is undefined" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
};

template <class T, class Tder>
OctaveConstitutiveLaw<T, Tder>::~OctaveConstitutiveLaw(void) {
        NO_OP;
};

template <class T, class Tder>
ConstLawType::Type OctaveConstitutiveLaw<T, Tder>::GetConstLawType(void) const {

        switch (clType) {
        case ConstLawType::UNKNOWN:
                break;
        default:
                return clType;
        }

        static const struct {
                ConstLawType::Type value;
                char name[13];
        } types[] = {
                        { ConstLawType::ELASTIC,		"ELASTIC" },
                        { ConstLawType::VISCOUS,		"VISCOUS" },
                        { ConstLawType::VISCOELASTIC,	"VISCOELASTIC" }
        };

        static const int count = sizeof(types)/sizeof(types[0]);

        octave_value_list args(octObject);
        const octave_value_list ans = GetInterface()->EvalFunction(strGetConstLawType, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!ans(0).is_string()) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): " << GetClass() << "." << strGetConstLawType
                                << " returned a invalid data type (" << ans(0).type_name() << std::endl
                                << "expected string" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const std::string strCLType(ans(0).string_value());

        for (int i = 0; i < count; ++i) {
                if (strCLType == types[i].name) {
                        clType = types[i].value;
                        return clType;
                }
        }

        silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                        << "): " << GetClass() << "." << strGetConstLawType
                        << " returned a invalid value (" << strCLType << ")" << std::endl);

        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
};

template <class T, class Tder>
ConstitutiveLaw<T, Tder>* OctaveConstitutiveLaw<T, Tder>::pCopy(void) const {
        Base_t* pCL = NULL;

        typedef OctaveConstitutiveLaw<T, Tder> cl;
        SAFENEWWITHCONSTRUCTOR(pCL, cl, cl(GetClass(), GetInterface(), GetFlags()));
        return pCL;
};

template <class T, class Tder>
std::ostream& OctaveConstitutiveLaw<T, Tder>::Restart(std::ostream& out) const {
        return out << "# octave constitutive law not implemented" << std::endl;
};

template <class T, class Tder>
void OctaveConstitutiveLaw<T, Tder>::Update(const T& mbEps, const T& mbEpsPrime) {

        octave_value octEps, octEpsPrime;

        if (!GetInterface()->ConvertMBDynToOctave(mbEps, octEps)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): failed to convert MBDyn data type to octave_value" << std::endl);

                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!GetInterface()->ConvertMBDynToOctave(mbEpsPrime, octEpsPrime)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): failed to convert MBDyn data type to octave_value" << std::endl);

                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octave_value_list args(octObject);
        args.append(octEps);
        args.append(octEpsPrime);

        const octave_value_list ans = GetInterface()->EvalFunction(strUpdate, args, 4, GetFlags());

        ASSERT(ans.length() == 4);

        if (!(ans(0).isobject() && ans(0).type_id() == octObject.type_id())) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                        << "): output argument 1 is not an instance of " << octObject.type_name() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octObject = ans(0);

        if (!((ans(1).is_real_matrix() || ans(1).is_real_scalar()) && ans(1).columns() == 1)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument F is invalid (" << ans(1).type_name() << ")" << std::endl
                                << "expected real vector" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const ColumnVector octF(ans(1).column_vector_value());

        if (!(ans(2).is_real_matrix() || ans(2).is_real_scalar())) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument FDE is invalid (" << ans(2).type_name() << ")" << std::endl
                                << "expected real matrix" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Matrix octFDE(ans(2).matrix_value());

        if (!(ans(3).is_real_matrix() || ans(3).is_real_scalar())) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument FDEPrime is invalid (" << ans(3).type_name() << ")" << std::endl
                                << "expected real matrix" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Matrix octFDEPrime(ans(3).matrix_value());

        Base_t::Epsilon = mbEps;
        Base_t::EpsilonPrime = mbEpsPrime;

        if (!GetInterface()->ConvertOctaveToMBDyn(octF, Base_t::F)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument F is invalid "
                                << octF.rows() << "x" << octF.columns() << " (" << ans(1).type_name() << ")" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!GetInterface()->ConvertOctaveToMBDyn(octFDE, Base_t::FDE)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument FDE is invalid "
                                << octFDE.rows() << "x" << octFDE.columns() << " (" << ans(2).type_name() << ")" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!GetInterface()->ConvertOctaveToMBDyn(octFDEPrime, Base_t::FDEPrime)) {
                silent_cerr("octave constitutive law(" << Base_t::GetLabel()
                                << "): data type of output argument FDEPrime is invalid "
                                << octFDEPrime.rows() << "x" << octFDEPrime.columns() << " (" << ans(3).type_name() << ")" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
};

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

void OctaveBaseDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)const
{
        if (pInterface) {
                pInterface->Destroy();
        }
        
        pInterface = OctaveInterface::CreateInterface(pDM, &HP);

        strFunction = HP.GetStringWithDelims(HighParser::DOUBLEQUOTE);

        iFlags = OctaveInterface::DEFAULT_CALL_FLAGS;

        if (HP.IsKeyWord("update" "octave" "variables") && HP.GetYesNoOrBool()) {
                iFlags |= OctaveInterface::UPDATE_OCTAVE_VARIABLES;
        }

        if (HP.IsKeyWord("update" "mbdyn" "variables") && HP.GetYesNoOrBool()) {
                iFlags |= OctaveInterface::UPDATE_MBDYN_VARIABLES;
        }

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

void OctaveFunctionDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) const
{
        args.resize(0);

        OctaveBaseDCR::Read(pDM, HP, bDeferred);

        if (HP.IsKeyWord("arguments")) {
                integer iNumArgs = HP.GetInt();
                while (iNumArgs-- > 0) {
                        TypedValue mbVal(HP.GetValue(TypedValue()));
                        octave_value octVal;

                        if (!GetInterface()->ConvertMBDynToOctave(mbVal, octVal)) {
                                silent_cerr("octave error: could not convert MBDyn data type " << mbVal.GetTypeName() << " to octave" << std::endl);
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                        }

                        args.append(octVal);
                }
        }
}

DriveCaller *
OctaveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
        OctaveFunctionDCR::Read(pDM, HP, bDeferred);

        return new OctaveDriveCaller(GetFunction(), GetInterface(), GetFlags(), GetArgs());
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
        OctaveFunctionDCR::Read(pDM, HP);
        return new OctaveScalarFunction(GetFunction(), GetInterface(), GetFlags(), GetArgs());
};


const std::string OctaveElement::strWorkSpaceDim("WorkSpaceDim");
const std::string OctaveElement::striGetNumDof("iGetNumDof");
const std::string OctaveElement::strAssRes("AssRes");
const std::string OctaveElement::strAssJac("AssJac");
const std::string OctaveElement::strUpdate("Update");
const std::string OctaveElement::strSetValue("SetValue");
const std::string OctaveElement::striGetInitialNumDof("iGetInitialNumDof");
const std::string OctaveElement::strSetInitialValue("SetInitialValue");
const std::string OctaveElement::strInitialAssRes("InitialAssRes");
const std::string OctaveElement::strInitialAssJac("InitialAssJac");
const std::string OctaveElement::strInitialWorkSpaceDim("InitialWorkSpaceDim");
const std::string OctaveElement::strGetDofType("GetDofType");
const std::string OctaveElement::strGetEqType("GetEqType");
const std::string OctaveElement::strAfterConvergence("AfterConvergence");
const std::string OctaveElement::striGetNumPrivData("iGetNumPrivData");
const std::string OctaveElement::striGetPrivDataIdx("iGetPrivDataIdx");
const std::string OctaveElement::strdGetPrivData("dGetPrivData");
const std::string OctaveElement::striGetNumConnectedNodes("iGetNumConnectedNodes");
const std::string OctaveElement::strGetConnectedNodes("GetConnectedNodes");
const std::string OctaveElement::strOutput("Output");
const std::string OctaveElement::strDescribeDof("DescribeDof");
const std::string OctaveElement::strDescribeEq("DescribeEq");
const std::string OctaveElement::strAfterPredict("AfterPredict");
const std::string OctaveElement::strRestart("Restart");

OctaveElement::OctaveElement(
        unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)),
  UserDefinedElem(uLabel, pDO),
  haveMethod(HAVE_DEFAULT)
{
        // help
        if (HP.IsKeyWord("help")) {
                silent_cout("Module:    octave\n"
                        << std::endl);

                if (!HP.IsArg()) {
                        /*
                         * Exit quietly if nothing else is provided
                         */
                        throw NoErr(MBDYN_EXCEPT_ARGS);
                }
        }

        dcr.Read(pDM, HP);

        X = new ConstVectorHandlerInterface(GetInterface());
        XP = new ConstVectorHandlerInterface(GetInterface());
        mbdObject = new OctaveElementInterface(GetInterface(), this);
        OS = new OStreamInterface(GetInterface());

        octave_value_list args(mbdObject);

        args.append(GetInterface()->GetDataManagerInterface());
        args.append(GetInterface()->GetMBDynParserInterface());
        args.append(OS);
        const int flags = OctaveInterface::UPDATE_OCTAVE_VARIABLES;

        OS->Set(&pDM->GetLogFile());

        octave_value_list ans = GetInterface()->EvalFunction(GetClass(), args, 1, flags);

        OS->Set(0);

        ASSERT(ans.length() == 1);

        octObject = ans(0);

        if (!octObject.isobject()) {
                silent_cerr("octave(" << GetLabel() << "): result of constructor ("
                                << ans(0).type_name() << ") is not an object at line "
                                << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

        if (bHaveMethod(strAssJac)) {
                haveMethod |= HAVE_JACOBIAN;
        } else {
                pedantic_cerr("octave waring: method " << strAssJac << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strUpdate)) {
                haveMethod |= HAVE_UPDATE;
        } else {
                pedantic_cerr("octave warning: method " << strUpdate << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strAfterConvergence)) {
                haveMethod |= HAVE_AFTER_CONVERGENCE;
        } else {
                pedantic_cerr("octave warning: method " << strAfterConvergence << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strSetValue)) {
                haveMethod |= HAVE_SET_VALUE;
        } else {
                pedantic_cerr("octave warning: method " << strSetValue << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(striGetNumDof)) {
                if (!bHaveMethod(strGetDofType)) {
                        silent_cerr("octave error: method " << strGetDofType << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!bHaveMethod(strGetEqType)) {
                        silent_cerr("octave error: method " << strGetEqType << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                haveMethod |= HAVE_PRIVATE_DOF;

                if (bHaveMethod(strDescribeDof)) {
                        haveMethod |= HAVE_DESCRIBE_DOF;
                }

                if (bHaveMethod(strDescribeEq)) {
                        haveMethod |= HAVE_DESCRIBE_EQ;
                }
        } else {
                pedantic_cerr("octave warning: method " << striGetNumDof << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(striGetNumPrivData)) {
                if (!bHaveMethod(striGetPrivDataIdx)) {
                        silent_cerr("octave error: method " << striGetPrivDataIdx << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!bHaveMethod(strdGetPrivData)) {
                        silent_cerr("octave error: method " << strdGetPrivData << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                haveMethod |= HAVE_PRIVATE_DATA;
        } else {
                pedantic_cerr("octave warning: method " << striGetNumPrivData << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(striGetNumConnectedNodes)) {
                if (!bHaveMethod(strGetConnectedNodes)) {
                        silent_cerr("octave error: method " << strGetConnectedNodes << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                haveMethod |= HAVE_CONNECTED_NODES;
        } else {
                pedantic_cerr("octave warning: method " << striGetNumConnectedNodes << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strInitialAssRes)) {
                if (!bHaveMethod(strInitialAssJac)) {
                        silent_cerr("octave error: method " << strInitialAssJac << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!bHaveMethod(strInitialWorkSpaceDim)) {
                        silent_cerr("octave error: method " << strInitialWorkSpaceDim << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (!bHaveMethod(striGetInitialNumDof)) {
                        silent_cerr("octave error: method " << striGetInitialNumDof << " not found in class " << GetClass() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                haveMethod |= HAVE_INITIAL_ASSEMBLY;

                if (bHaveMethod(strSetInitialValue)) {
                        haveMethod |= HAVE_SET_INITIAL_VALUE;
                } else {
                        pedantic_cerr("octave warning: method " << strSetInitialValue << " not found in class " << GetClass() << std::endl);
                }
        } else {
                pedantic_cerr("octave warning: method " << strInitialAssRes
                                          << " not found in class " << GetClass() << std::endl
                                          << " initial assembly has been disabled" << std::endl);
        }

        if (bHaveMethod(strAfterPredict)) {
                haveMethod |= HAVE_AFTER_PREDICT;
        } else {
                pedantic_cerr("octave warning: method " << strAfterPredict << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strRestart)) {
                haveMethod |= HAVE_RESTART;
        } else {
                pedantic_cerr("octave warning: method " << strRestart << " not found in class " << GetClass() << std::endl);
        }

        if (bHaveMethod(strOutput)) {
                haveMethod |= HAVE_OUTPUT;
        } else {
                pedantic_cerr("octave warning: method " << strOutput << " not found in class " << GetClass() << std::endl);
        }
}

OctaveElement::~OctaveElement(void)
{

}

void
OctaveElement::Output(OutputHandler& OH) const
{
        if (!(bToBeOutput() && OH.UseText(OutputHandler::LOADABLE))) {
                return;
        }

        if (!(haveMethod & HAVE_OUTPUT)) {
                return;
        }

        OS->Set(&OH.Loadable());

        octave_value_list args(octObject);
        args.append(OS);

        GetInterface()->EvalFunction(strOutput, args, 0, GetFlags());

        OS->Set(0);
}

void
OctaveElement::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(strWorkSpaceDim, args, 2, GetFlags());

        ASSERT(ans.length() == 2);

        if (!(ans(0).is_scalar_type()
                        && ans(0).isinteger()
                        && ans(1).is_scalar_type()
                        && ans(1).isinteger())) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strWorkSpaceDim
                                << " must return two integer values\nreturned(" << ans(0).type_name() << ", " << ans(1).type_name() << ")" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        *piNumRows = static_cast<int32_t>(ans(0).int32_scalar_value());
        *piNumCols = static_cast<int32_t>(ans(1).int32_scalar_value());
}

VariableSubMatrixHandler&
OctaveElement::AssJac(VariableSubMatrixHandler& WorkMatVar,
        doublereal dCoef,
        const VectorHandler& XCurr,
        const VectorHandler& XPrimeCurr)
{
        if ( !(haveMethod & HAVE_JACOBIAN) ) {
                WorkMatVar.SetNullMatrix();
                return WorkMatVar;
        }

        integer iNumRows, iNumCols;

        WorkSpaceDim(&iNumRows, &iNumCols);

        octave_value_list args(octObject);

        X->Set(&XCurr);
        XP->Set(&XPrimeCurr);

        args.append(octave_value(dCoef));
        args.append(X);
        args.append(XP);

        const octave_value_list ans = GetInterface()->EvalFunction(strAssJac, args, 5, GetFlags() | OctaveInterface::OPTIONAL_OUTPUT_ARGS);

        ASSERT(ans.length() <= 5);

        X->Set(0);
        XP->Set(0);

        if (ans.length() >= 5) {
                if (!(ans(4).isobject() && ans(4).type_id() == octObject.type_id())) {
                        silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": output argument 5 is not an instance of " << octObject.type_name() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
                octObject = ans(4);
        }

        return AssMatrix(WorkMatVar, ans, false);
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

        X->Set(&XCurr);
        XP->Set(&XPrimeCurr);

        octave_value_list args(octObject);
        args.append(octave_value(dCoef));
        args.append(X);
        args.append(XP);

        octave_value_list ans = GetInterface()->EvalFunction(strAssRes, args, 3, GetFlags() | OctaveInterface::OPTIONAL_OUTPUT_ARGS);

        X->Set(0);
        XP->Set(0);

        ASSERT(ans.length() <= 3);

        if (ans.length() < 2) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": expected 2-3 output arguments" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        ASSERT(ans.length() >= 2);
        if (!(ans(0).is_matrix_type()
                        && ans(0).isreal()
                        && ans(0).columns() == 1)) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strAssRes
                                << ": output argument f (" << ans(0).type_name()
                                << ") must be a real column vector" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const ColumnVector f = ans(0).column_vector_value();

        if (!(ans(1).isinteger()
                        && ans(1).is_matrix_type()
                        && ans(1).columns() == 1)) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strAssRes
                                << ": output argument ridx (" << ans(1).type_name()
                                << ") must be an integer column vector" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32NDArray ridx = ans(1).int32_array_value();

        if (ridx.numel() != iNumRows || f.numel() != iNumRows) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes
                                << "\nlength(f)=" << f.numel()
                                << " length(ridx)=" << ridx.numel()
                                << " is not equal to iNumRows=" << iNumRows << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const integer iNumDof = GetInterface()->GetDataManager()->iGetNumDofs();

        for (octave_idx_type i = 0; i < f.numel(); ++i) {
                if (int(ridx(i)) <= 0 || int(ridx(i)) > iNumDof) {
                        silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes
                                        << ": row index " << ridx(i)
                                        << " out of range [1:" << iNumDof << "]" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                WorkVec.PutItem(i + 1, ridx(i), f(i));
        }

        if (ans.length() >= 3) {
                if (!(ans(2).isobject() && ans(2).type_id() == octObject.type_id())) {
                        silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": output argument 3 is not an instance of " << octObject.type_name() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
                octObject = ans(2);
        }

        return WorkVec;
}

unsigned int
OctaveElement::iGetNumPrivData(void) const
{
        if ( !(haveMethod & HAVE_PRIVATE_DATA) ) {
                return 0u;
        }

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(striGetNumPrivData, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if ( !(ans(0).isinteger() && ans(0).is_scalar_type()) ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetNumPrivData
                                << " returned (" << ans(0).type_name() << ")\n"
                                << "expected integer scalar" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32_t iNumPrivData = ans(0).int32_scalar_value();

        if ( iNumPrivData < 0 ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetNumPrivData
                                << " returned a negative value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return iNumPrivData;
}

unsigned int
OctaveElement::iGetPrivDataIdx(const char *s) const
{
        if ( !(haveMethod & HAVE_PRIVATE_DATA) ) {
                ASSERT(0); // We should not get here!
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octave_value_list args(octObject);
        args.append(octave_value(s));

        octave_value_list ans = GetInterface()->EvalFunction(striGetPrivDataIdx, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if ( !(ans(0).isinteger() && ans(0).is_scalar_type()) ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetPrivDataIdx
                                << " returned (" << ans(0).type_name() << ")\n"
                                << "expected integer scalar" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32_t iPrivDataIdx = ans(0).int32_scalar_value();

        if ( iPrivDataIdx <= 0 ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetPrivDataIdx
                                << " returned a negative or zero value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return iPrivDataIdx;
}

doublereal
OctaveElement::dGetPrivData(unsigned int i) const
{
        if ( !(haveMethod & HAVE_PRIVATE_DATA) ) {
                ASSERT(0); // We should not get here!
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octave_value_list args(octObject);
        args.append(octave_int32(i));

        octave_value_list ans = GetInterface()->EvalFunction(strdGetPrivData, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if ( !(ans(0).is_real_scalar()) ) {
                silent_cerr("octave error: method " << GetClass() << "." << strdGetPrivData
                                << " returned (" << ans(0).type_name() << ")\n"
                                << "expected real scalar" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return ans(0).scalar_value();
}

int
OctaveElement::iGetNumConnectedNodes(void) const
{
        if (!(haveMethod & HAVE_CONNECTED_NODES)) {
                return 0u;
        }

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(striGetNumConnectedNodes, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if ( !(ans(0).isinteger() && ans(0).is_scalar_type()) ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetNumConnectedNodes
                                << " returned (" << ans(0).type_name() << ")\n"
                                << "expected integer scalar" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32_t iNumConnectedNodes = ans(0).int32_scalar_value();

        if ( iNumConnectedNodes < 0 ) {
                silent_cerr("octave error: method " << GetClass() << "." << striGetNumConnectedNodes
                                << " returned a negative value" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return iNumConnectedNodes;
}

void
OctaveElement::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
        if (!(haveMethod & HAVE_CONNECTED_NODES)) {
                connectedNodes.resize(0);
                return;
        }

        connectedNodes.resize(iGetNumConnectedNodes());

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(strGetConnectedNodes, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!ans(0).iscell()) {
                silent_cerr("octave error: method " << GetClass() << "." << strGetConnectedNodes
                                << " returned (" << ans(0).type_name() << ")\n"
                                << "expected cell" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Cell nodes(ans(0).cell_value());

        if (size_t(nodes.numel()) != connectedNodes.size()) {
                silent_cerr("octave error: method " << GetClass() << "." << strGetConnectedNodes
                                << " returned an object of size " << nodes.numel()
                                << " whereas " << GetClass() << "." << striGetNumConnectedNodes
                                << " returned "<< connectedNodes.size() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        for (size_t i = 0; i < connectedNodes.size(); ++i) {
                const NodeInterface* pNode = dynamic_cast<const NodeInterface*>(&nodes(i).get_rep());

                if (pNode == 0) {
                        silent_cerr("octave error: cell(" << i + 1 << ") (data type " << nodes(i).type_name() << ") is not a node" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
                connectedNodes[i] = pNode->Get();
        }
}

void
OctaveElement::SetValue(DataManager *pDM,
        VectorHandler& X, VectorHandler& XP,
        SimulationEntity::Hints *ph)
{
        if (!(haveMethod & HAVE_SET_VALUE)) {
                return;
        }

        octave_value_list args(octObject);
        args.append(octave_value(new VectorHandlerInterface(GetInterface(), &X)));
        args.append(octave_value(new VectorHandlerInterface(GetInterface(), &XP)));

        GetInterface()->EvalFunction(strSetValue, args, 0, GetFlags());
}

unsigned int OctaveElement::iGetNumDof(void) const
{
        if (!(haveMethod & HAVE_PRIVATE_DOF)) {
                return 0u;
        }

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(striGetNumDof, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!(ans(0).isinteger() && ans(0).is_scalar_type())){
                silent_cerr("octave error: method " << GetClass() << "." << striGetNumDof
                                << " returned (" << ans(0).type_name()
                                << ")\nmethod must return an integer" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return static_cast<int32_t>(ans(0).int32_scalar_value());
}

DofOrder::Order OctaveElement::GetDofType(unsigned int i) const
{
        if (!(haveMethod & HAVE_PRIVATE_DOF)) {
                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octave_value_list args(octObject);
        args.append(octave_int32(i + 1)); // FIXME: Why is this index zero based?

        octave_value_list ans = GetInterface()->EvalFunction(strGetDofType, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (ans(0).isinteger() && ans(0).is_scalar_type()) {
                const int32_t order = static_cast<int32_t>(ans(0).int32_scalar_value());

                switch ( order ) {
                case DofOrder::ALGEBRAIC:
                case DofOrder::DIFFERENTIAL:
                        return static_cast<DofOrder::Order>(order);
                default:
                        silent_cerr("octave error: method " << GetClass() << "." << strGetDofType
                                        << " returned an illegal value" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }

        if ( !ans(0).is_string() ) {
                silent_cerr("octave error: method " << GetClass() << "." << strGetDofType
                                << " returned (" << ans(0).type_name()
                                << ")\nmethod must return integer or string values" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const std::string order = ans(0).string_value();

        if ( order == "ALGEBRAIC" ) {
                return DofOrder::ALGEBRAIC;
        } else if ( order == "DIFFERENTIAL" ) {
                return DofOrder::DIFFERENTIAL;
        } else {
                silent_cerr("octave error: method " << GetClass() << "." << strGetDofType
                                << " returned an illegal value (" << order << ")" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

DofOrder::Order OctaveElement::GetEqType(unsigned int i) const
{
        if (!(haveMethod & HAVE_PRIVATE_DOF)) {
                ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octave_value_list args(octObject);
        args.append(octave_int32(i + 1)); // FIXME: Why is this index zero based?

        octave_value_list ans = GetInterface()->EvalFunction(strGetEqType, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (ans(0).isinteger() && ans(0).is_scalar_type()) {
                const int32_t order = static_cast<int32_t>(ans(0).int32_scalar_value());

                switch ( order ) {
                case DofOrder::ALGEBRAIC:
                case DofOrder::DIFFERENTIAL:
                        return static_cast<DofOrder::Order>(order);
                default:
                        silent_cerr("octave error: method " << GetClass() << "." << strGetEqType
                                        << " returned an illegal value" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
        }

        if ( !ans(0).is_string() ) {
                silent_cerr("octave error: method " << GetClass() << "." << strGetEqType
                                << " returned (" << ans(0).type_name()
                                << ")\nmethod must return integer or string values" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const std::string order = ans(0).string_value();

        if ( order == "ALGEBRAIC" ) {
                return DofOrder::ALGEBRAIC;
        } else if ( order == "DIFFERENTIAL" ) {
                return DofOrder::DIFFERENTIAL;
        } else {
                silent_cerr("octave error: method " << GetClass() << "." << strGetEqType
                                << " returned an illegal value (" << order << ")" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
}

std::ostream& OctaveElement::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
        if (!(haveMethod & HAVE_DESCRIBE_DOF)) {
                return out;
        }

        OS->Set(&out);

        octave_value_list args(octObject);
        args.append(OS);
        args.append(octave_value(prefix));
        args.append(octave_value(bInitial));

        GetInterface()->EvalFunction(strDescribeDof, args, 0, GetFlags());

        OS->Set(0);

        return out;
}

std::ostream& OctaveElement::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
        if (!(haveMethod & HAVE_DESCRIBE_EQ)) {
                return out;
        }

        OS->Set(&out);

        octave_value_list args(octObject);
        args.append(OS);
        args.append(octave_value(prefix));
        args.append(octave_value(bInitial));

        GetInterface()->EvalFunction(strDescribeEq, args, 0, GetFlags());

        OS->Set(0);

        return out;
}

void OctaveElement::Update(const VectorHandler& XCurr,const VectorHandler& XPrimeCurr)
{
        if (!(haveMethod & HAVE_UPDATE)) {
                return;
        }

        X->Set(&XCurr);
        XP->Set(&XPrimeCurr);

        octave_value_list args(octObject);
        args.append(X);
        args.append(XP);

        octave_value_list ans = GetInterface()->EvalFunction(strUpdate, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!(ans(0).isobject() && ans(0).type_id() == octObject.type_id())) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strUpdate << ": output argument 1 is not an instance of " << octObject.type_name() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octObject = ans(0);

        X->Set(0);
        XP->Set(0);
}

void OctaveElement::AfterPredict(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
{
        if (!(haveMethod & HAVE_AFTER_PREDICT)) {
                return;
        }

        X->Set(&XCurr);
        XP->Set(&XPrimeCurr);

        octave_value_list args(octObject);
        args.append(X);
        args.append(XP);

        octave_value_list ans = GetInterface()->EvalFunction(strAfterPredict, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!(ans(0).isobject() && ans(0).type_id() == octObject.type_id())) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAfterPredict << ": output argument 1 is not an instance of " << octObject.type_name() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octObject = ans(0);

        X->Set(0);
        XP->Set(0);
}

void OctaveElement::AfterConvergence(const VectorHandler& X,
                        const VectorHandler& XP)
{
        if ( !(haveMethod & HAVE_AFTER_CONVERGENCE) ) {
                return;
        }

        this->X->Set(&X);
        this->XP->Set(&XP);

        octave_value_list args(octObject);
        args.append(this->X);
        args.append(this->XP);

        octave_value_list ans = GetInterface()->EvalFunction(strAfterConvergence, args, 1, GetFlags());

        this->X->Set(0);
        this->XP->Set(0);

        ASSERT(ans.length() == 1);

        if ( !ans(0).isobject() ) {
                silent_cerr("octave error: return value of " << strAfterConvergence << " ("
                                << ans(0).type_name() << ") is not an object" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        octObject = ans(0);
}

std::ostream&
OctaveElement::Restart(std::ostream& out) const
{
        if (!(haveMethod & HAVE_RESTART)) {
                out << "# OctaveElement(" << GetClass() << "): not implemented" << std::endl;
                return out;
        }

        OS->Set(&out);

        octave_value_list args(octObject);
        args.append(OS);

        GetInterface()->EvalFunction(strRestart, args, 0, GetFlags());

        OS->Set(0);

        return out;
}

unsigned int
OctaveElement::iGetInitialNumDof(void) const
{
        if (!(haveMethod & HAVE_INITIAL_ASSEMBLY)) {
                return 0u;
        }

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(striGetInitialNumDof, args, 1, GetFlags());

        ASSERT(ans.length() == 1);

        if (!(ans(0).isinteger() && ans(0).is_scalar_type())){
                silent_cerr("octave error: method " << GetClass() << "." << striGetInitialNumDof
                                << " returned (" << ans(0).type_name()
                                << ")\nmethod must return an integer" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return static_cast<int32_t>(ans(0).int32_scalar_value());
}

void
OctaveElement::InitialWorkSpaceDim(
        integer* piNumRows,
        integer* piNumCols) const
{
        if (!(haveMethod & HAVE_INITIAL_ASSEMBLY)) {
                *piNumRows = 0;
                *piNumCols = 0;
                return;
        }

        octave_value_list args(octObject);

        octave_value_list ans = GetInterface()->EvalFunction(strInitialWorkSpaceDim, args, 2, GetFlags());

        ASSERT(ans.length() == 2);

        if (!(ans(0).isinteger()
                        && ans(0).is_scalar_type()
                        && ans(1).isinteger()
                        && ans(1).is_scalar_type())) {
                silent_cerr("octave error: method " << GetClass() << "." << strInitialWorkSpaceDim
                                << " must return two integer values\n"
                                << "returned (" << ans(0).type_name() << ", " << ans(1).type_name() << ")" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        *piNumRows = static_cast<int32_t>(ans(0).int32_scalar_value());
        *piNumCols = static_cast<int32_t>(ans(1).int32_scalar_value());
}

VariableSubMatrixHandler&
OctaveElement::InitialAssJac(
        VariableSubMatrixHandler& WorkMatVar,
        const VectorHandler& XCurr)
{
        if (!(haveMethod & HAVE_INITIAL_ASSEMBLY)) {
                WorkMatVar.SetNullMatrix();
                return WorkMatVar;
        }

        octave_value_list args(octObject);

        X->Set(&XCurr);

        args.append(X);

        const octave_value_list ans = GetInterface()->EvalFunction(strInitialAssJac, args, 4, GetFlags() | OctaveInterface::OPTIONAL_OUTPUT_ARGS);

        ASSERT(ans.length() <= 4);

        X->Set(0);

        return AssMatrix(WorkMatVar, ans, true);
}

SubVectorHandler&
OctaveElement::InitialAssRes(
        SubVectorHandler& WorkVec,
        const VectorHandler& XCurr)
{
        if ( !(haveMethod & HAVE_INITIAL_ASSEMBLY) ) {
                WorkVec.ResizeReset(0);
                return WorkVec;
        }

        int iNumRows, iNumCols;

        InitialWorkSpaceDim(&iNumRows, &iNumCols);

        WorkVec.ResizeReset(iNumRows);

        X->Set(&XCurr);

        octave_value_list args(octObject);
        args.append(X);

        octave_value_list ans = GetInterface()->EvalFunction(strInitialAssRes, args, 2, GetFlags());

        X->Set(0);

        ASSERT(ans.length() == 2);

        if (!(ans(0).is_matrix_type()
                        && ans(0).isreal()
                        && ans(0).columns() == 1)) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strInitialAssRes
                                << " output argument f (" << ans(0).type_name() << ") must be a real column vector" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const ColumnVector f = ans(0).column_vector_value();

        if (!(ans(1).isinteger()
                        && ans(1).is_matrix_type()
                        && ans(1).columns() == 1)) {
                silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strInitialAssRes
                                << ": output argument ridx (" << ans(1).type_name() << ") must be an integer column vector" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32NDArray ridx = ans(1).int32_array_value();

        if (ridx.numel() != iNumRows || f.numel() != iNumRows) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strInitialAssRes
                                << ": length(f)=" << f.numel()
                                << " length(ridx)=" << ridx.numel()
                                << " is not equal to iNumRows=" << iNumRows << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const integer iNumDof = GetInterface()->GetDataManager()->iGetNumDofs();

        for (octave_idx_type i = 0; i < f.numel(); ++i) {
                if (int(ridx(i)) <= 0 || int(ridx(i)) > iNumDof) {
                        silent_cerr("octave(" << GetLabel() << "): function " << GetClass() << "." << strAssRes << ": row index " << ridx(i) << " out of range [1:" << iNumDof << "]" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                WorkVec.PutItem(i + 1, ridx(i), f(i));
        }

        return WorkVec;
}

void OctaveElement::SetInitialValue(VectorHandler& X)
{
        if (!(haveMethod & HAVE_SET_INITIAL_VALUE)) {
                return;
        }

        octave_value_list args(octObject);
        args.append(octave_value(new VectorHandlerInterface(GetInterface(), &X)));

        GetInterface()->EvalFunction(strSetInitialValue, args, 0, GetFlags());
}

bool OctaveElement::bHaveMethod(const std::string& strName)const
{
        return GetInterface()->bHaveMethod(octObject, GetClass(), strName);
}

VariableSubMatrixHandler&
OctaveElement::AssMatrix(VariableSubMatrixHandler& WorkMatVar, const octave_value_list& ans, bool bInitial)
{
        const std::string& strFunction = bInitial ? strInitialAssJac : strAssJac;

        ASSERT(ans.length() <= 5);

        if (ans.length() < 3) {
                silent_cerr("octave(" << GetLabel()
                                        << "): function " << GetClass() << "." << strFunction
                                        << " returned " << ans.length() << " output arguments\n"
                                        << "expected 3-4 output arguments" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        ASSERT(ans.length() >= 3);

        if (!ans(0).is_real_matrix()) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strFunction
                                << " output argument Jac (" << ans(0).type_name()
                                << ") must be a real matrix" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const Matrix Jac = ans(0).matrix_value();

        if (!(ans(1).isinteger() && ans(1).is_matrix_type())) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strFunction
                                << " output argument ridx (" << ans(1).type_name()
                                << ") must be a vector of integer values" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32NDArray ridx = ans(1).int32_array_value();

        if (!(ans(2).isinteger() && ans(2).is_matrix_type())) {
                silent_cerr("octave(" << GetLabel()
                                << "): function " << GetClass() << "." << strFunction
                                << " output argument cidx (" << ans(2).type_name()
                                << ") must be a vector of integer values" << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const int32NDArray cidx = ans(2).int32_array_value();

        bool bSparse = false;

        if (ans.length() >= 4) {
                if (!ans(3).is_bool_scalar()) {
                        silent_cerr("octave(" << GetLabel()
                                        << "): invalid data type (" << ans(3).type_name() << ")\n"
                                        << "expected bool scalar" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                bSparse = ans(3).bool_value();
        }

        return AssMatrix(WorkMatVar, Jac, ridx, cidx, bSparse, bInitial);
}

VariableSubMatrixHandler&
OctaveElement::AssMatrix(VariableSubMatrixHandler& WorkMatVar, const Matrix& Jac, const int32NDArray& ridx, const int32NDArray& cidx, bool bSparse, bool bInitial)
{
        const std::string& strFunction = bInitial ? strInitialAssJac : strAssJac;

        const integer iNumRows = Jac.rows();
        const integer iNumCols = Jac.cols();

        const integer iNumDof = GetInterface()->GetDataManager()->iGetNumDofs();

        if (bSparse) {
                if (iNumCols != 1
                                || ridx.numel() != iNumRows
                                || cidx.numel() != iNumRows) {
                        silent_cerr("octave(" << GetLabel() << "):"
                                        << " rows(Jac)=" << Jac.rows()
                                        << " columns(Jac)= " << Jac.columns()
                                        << " length(ridx)= " << ridx.numel()
                                        << " length(cidx)=" << cidx.numel()
                                        << " are not consistent for a sparse submatrix" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                SparseSubMatrixHandler& WorkMat = WorkMatVar.SetSparse();
                WorkMat.ResizeReset(iNumRows, 1);

                for (int i = 0; i < iNumRows; ++i) {
                        if (int32_t(ridx(i)) <= 0 || int32_t(ridx(i)) > iNumDof) {
                                silent_cerr("octave(" << GetLabel() << "):"
                                                << " function " << GetClass() << "." << strFunction
                                                << ": row index " << ridx(i)
                                                << " out of range [1:" << iNumDof << "]" << std::endl);
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                        }

                        if (int32_t(cidx(i)) <= 0 || int32_t(cidx(i)) > iNumDof) {
                                silent_cerr("octave(" << GetLabel() << "):"
                                                << " function " << GetClass() << "." << strFunction
                                                << ": column index " << cidx(i)
                                                << " out of range [1:" << iNumDof << "]" << std::endl);
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                        }

                        WorkMat.PutItem(i + 1, ridx(i), cidx(i), Jac(i, 0));
                }
        } else {
                FullSubMatrixHandler& WorkMat = WorkMatVar.SetFull();

                WorkMat.ResizeReset(iNumRows, iNumCols);

                if (ridx.numel() != iNumRows
                                || cidx.numel() != iNumCols) {
                        silent_cerr("octave(" << GetLabel() << "):"
                                        << " rows(Jac)=" << Jac.rows()
                                        << " columns(Jac)= " << Jac.columns()
                                        << " length(ridx)= " << ridx.numel()
                                        << " length(cidx)=" << cidx.numel()
                                        << " are not consistent for a full submatrix" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                for (int i = 0; i < iNumRows; ++i) {
                        if (int32_t(ridx(i)) <= 0 || int32_t(ridx(i)) > iNumDof) {
                                silent_cerr("octave(" << GetLabel() << "):"
                                                << " function " << GetClass() << "." << strFunction
                                                << ": row index " << ridx(i)
                                                << " out of range [1:" << iNumDof << "]" << std::endl);
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                        }

                        WorkMat.PutRowIndex(i + 1, ridx(i));
                }

                for (int j = 0; j < iNumCols; ++j) {
                        if (int32_t(cidx(j)) <= 0 || int32_t(cidx(j)) > iNumDof) {
                                silent_cerr("octave(" << GetLabel() << "):"
                                                << " function " << GetClass() << "." << strFunction
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
        }

        return WorkMatVar;
}

OctaveElementInterface::OctaveElementInterface(OctaveInterface* pInterface, OctaveElement* pElem)
: MBDynInterface(pInterface),
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

METHOD_DEFINE(OctaveElementInterface, GetLabel, args, nargout)
{
        if (args.length() != 0) {
                error("OctaveElement: invalid number of arguments %ld\n"
                                "expected no arguments", long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<unsigned int>(pElem->GetLabel()));
}

METHOD_DEFINE(OctaveElementInterface, iGetFirstIndex, args, nargout)
{
        if (args.length() != 0) {
                error("OctaveElement: invalid number of arguments %ld\n"
                                "expected no arguments", long(args.length()));
                return octave_value();
        }

        return octave_value(octave_int<integer>(pElem->iGetFirstIndex()));
}

BEGIN_METHOD_TABLE(OctaveElementInterface, MBDynInterface)
        METHOD_DISPATCH(OctaveElementInterface, GetLabel)
        METHOD_DISPATCH(OctaveElementInterface, iGetFirstIndex)
END_METHOD_TABLE()

DEFINE_OCTAVE_ALLOCATOR(OctaveElementInterface)
DEFINE_OV_TYPEID_FUNCTIONS_AND_DATA(OctaveElementInterface, "MBDynElement", "MBDynElement");

}

#endif // USE_OCTAVE

bool
mbdyn_octave_set(void)
{
#ifdef USE_OCTAVE
        using namespace oct;

        DriveCallerRead	*rf = new OctaveDCR;
        if (!SetDriveCallerData("octave", rf)) {
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
        if (!SetDriveCallerData("derivative", rf)) {
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

        ConstitutiveLawRead<doublereal, doublereal> *rfcl1D = new OctaveCLR<doublereal, doublereal>;
        if (!SetCL1D("octave", rfcl1D)) {
                delete rfcl1D;
                return false;
        }

        ConstitutiveLawRead<Vec3, Mat3x3> *rfcl3D = new OctaveCLR<Vec3, Mat3x3>;
        if (!SetCL3D("octave", rfcl3D)) {
                delete rfcl3D;
                return false;
        }

        ConstitutiveLawRead<Vec6, Mat6x6> *rfcl6D = new OctaveCLR<Vec6, Mat6x6>;
        if (!SetCL6D("octave", rfcl6D)) {
                delete rfcl6D;
                return false;
        }
#else
        pedantic_cerr("warning: MBDyn has been configured without octave support\n"
                                  "warning: module-octave is not available in this version of MBDyn" << std::endl);
#endif

        // Return true also if USE_OCTAVE is not enabled
        // This prevents one assertion to fail in userelem.cc if built as a static module
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
