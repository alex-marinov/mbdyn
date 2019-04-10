/* $Header$ */
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
 * ACKNOWLEDGEMENTS:
 * Support for output with NetCDF is based on a contribution
 * by Patrick Rix <patrick.rix@online.de>
 */

/* gestore dell'output */

#ifndef OUTPUT_H
#define OUTPUT_H

/* se #define DEBUG_COUT l'output avviene su cout anziche' nei files */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <typeinfo>

#if defined(USE_NETCDF)
#if defined(USE_NETCDFC)
#include <netcdfcpp.h>
typedef const NcDim * MBDynNcDim;
typedef NcVar * MBDynNcVar;
typedef NcFile MBDynNcFile;
typedef NcType MBDynNcType;
#define MBDynNcInt ncLong
#define MBDynNcDouble ncDouble
#define MBDynNcChar ncChar
#elif defined(USE_NETCDF4)
#include <netcdf>
typedef netCDF::NcDim MBDynNcDim; // not const because cannot be if not a pointer (in this case)
typedef netCDF::NcVar MBDynNcVar;
typedef netCDF::NcFile MBDynNcFile;
typedef netCDF::NcType MBDynNcType;
#define MBDynNcInt netCDF::NcType::nc_INT /**< replaces long in netcdf4 */
#define MBDynNcDouble netCDF::NcType::nc_DOUBLE
#define MBDynNcChar netCDF::NcType::nc_CHAR
#endif
#define MbNcInt MBDynNcType(MBDynNcInt) /**< creates a NcType object for a int, makes the notation simpler */
#define MbNcDouble MBDynNcType(MBDynNcDouble) /**< creates a NcType object for a double, makes the notation simpler */
#define MbNcChar MBDynNcType(MBDynNcChar) /**< makes the notation simpler */
#endif

#include "myassert.h"
#include "except.h"
#include "solman.h"
#include "filename.h"

/* OutputHandler - begin */

class OutputHandler : public FileName {
public:
	enum OutFiles {
		UNKNOWN			= -1,
		FIRSTFILE		= 0,
		OUTPUT			= 0,	//  0
		STRNODES,
		ELECTRIC,
		ABSTRACT,
		INERTIA,
		JOINTS,				//  5
		FORCES,
		BEAMS,
		ROTORS,
		RESTART,
		RESTARTXSOL,			// 10
		AERODYNAMIC,
		HYDRAULIC,
		PRESNODES,
		LOADABLE,
		GENELS,				// 15
		PARTITION,
		AEROMODALS,
		REFERENCEFRAMES,
		LOG,
		AIRPROPS,			// 20
		PARAMETERS,
		EXTERNALS,
		MODAL,
		NETCDF,
		THERMALNODES,			// 25
		THERMALELEMENTS,
		PLATES,
		GRAVITY,
		DOFSTATS,
		DRIVECALLERS,			// 30
		TRACES,
		EIGENANALYSIS,			// NOTE: ALWAYS LAST!
		LASTFILE			// 33
	};

private:
	long currentStep;

public:
	inline void SetCurrentStep(long Step) {
		currentStep = Step;
	};
	inline void IncCurrentStep(void) {
		currentStep++;
#if defined(USE_NETCDFC)
		ncStart1[0] = this->GetCurrentStep();
#elif defined(USE_NETCDF4)
		ncStart1[0] = ncStart1x3[0] = ncStart1x3x3[0] = this->GetCurrentStep();
#endif  /* USE_NETCDF4 */
	};
       	inline long GetCurrentStep(void) const {
		return currentStep;
	};

private:

	// flag values
	enum {
		OUTPUT_NONE			= 0x00U,
		OUTPUT_USE_DEFAULT_PRECISION	= 0x01U,
		OUTPUT_USE_SCIENTIFIC		= 0x02U,

		OUTPUT_MAY_USE_TEXT		= 0x10U,
		OUTPUT_USE_TEXT			= 0x20U,
		OUTPUT_MAY_USE_NETCDF		= 0x40U,
		OUTPUT_USE_NETCDF		= 0x80U,

// REMEMBER TO MODIFY OUTPUT_PRIVATE AND OUTPUT_MASK WHEN THE ABOVE MASKS
// BECOME GRATER THAN OUTPUT_PRIVATE !

		LAST
	};

	/* Aggiungere qui i files che si desidera avere a disposizione */
	struct {
		std::ofstream*	pof;
		unsigned	flags;
	} OutData[LASTFILE];

	// NetCDF dimensions and global attributes related to the binary file
#ifdef USE_NETCDF
	MBDynNcDim m_DimTime;
	MBDynNcDim m_DimV1;
	MBDynNcDim m_DimV3;
	MBDynNcFile *m_pBinFile;   /* ! one ! binary NetCDF data file */
#endif /* USE_NETCDF */

	/* handlers to streams */
	std::ofstream ofOutput;      		/*  0 */
	std::ofstream ofStrNodes;
	std::ofstream ofElectric;
	std::ofstream ofAbstract;
	std::ofstream ofInertia;
	std::ofstream ofJoints;      		/*  5 */
	std::ofstream ofForces;
	std::ofstream ofBeams;
	std::ofstream ofRotors;
	std::ofstream ofRestart;
	std::ofstream ofRestartXSol; 		/* 10 */
	std::ofstream ofAerodynamic;
	std::ofstream ofHydraulic;
	std::ofstream ofPresNodes;
	std::ofstream ofLoadable;
	std::ofstream ofGenels;			/* 15 */
	std::ofstream ofPartition;
	std::ofstream ofAeroModals;
	std::ofstream ofReferenceFrames;	/* 20 */
	std::ofstream ofLog;
	std::ofstream ofAirProps;
	std::ofstream ofParameters;
	std::ofstream ofExternals;
	std::ofstream ofModal;			/* 25 */
	std::ofstream ofThermalNodes;
	std::ofstream ofThermalElements;
	std::ofstream ofPlates;
	std::ofstream ofGravity;
	std::ofstream ofDofStats;		/* 30 */
	std::ofstream ofDriveCallers;
	std::ofstream ofTraces;
	std::ofstream ofEigenanalysis;

	int iCurrWidth;
	int iCurrPrecision;
	int nCurrRestartFile;

	// private because we know we're using valid out index
	bool IsOpen(int out) const;
	bool UseDefaultPrecision(int out) const;
	bool UseScientific(int out) const;

	bool UseText(int out) const;
	bool UseNetCDF(int out) const;

	// Pseudo-constructor
	void OutputHandler_int(void);

public:
	OutputHandler(void);

	OutputHandler(const char* sFName, int iExtNum = -1);

	void Init(const char* sFName, int iExtNum = -1);

	virtual ~OutputHandler(void);

	/* Aggiungere qui le funzioni che aprono i singoli stream */
	bool Open(const OutputHandler::OutFiles out);

	/* Overload for eigenanalysis text output */
	bool Open(const int out, const std::string& postfix);
	bool IsOpen(const OutputHandler::OutFiles out) const;
	bool UseDefaultPrecision(const OutputHandler::OutFiles out) const;
	bool UseScientific(const OutputHandler::OutFiles out) const;

	void SetText(const OutputHandler::OutFiles out);
	void ClearText(void);
	void ClearText(const OutputHandler::OutFiles out);
	bool UseText(const OutputHandler::OutFiles out) const;

	void SetNetCDF(const OutputHandler::OutFiles out);
	void ClearNetCDF(void);
	void ClearNetCDF(const OutputHandler::OutFiles out);
	bool UseNetCDF(const OutputHandler::OutFiles out) const;

	bool Close(const OutputHandler::OutFiles out);

	bool OutputOpen(void);
	bool RestartOpen(bool openResXSol = false);

	bool PartitionOpen(void);
	bool LogOpen(void);

	/* Aggiungere qui le funzioni che ritornano gli stream desiderati */
	inline std::ostream& Get(const OutputHandler::OutFiles f);

	inline std::ostream& Output(void) const;
	inline std::ostream& StrNodes(void) const;
	inline std::ostream& Electric(void) const;
	inline std::ostream& Abstract(void) const;
	inline std::ostream& Inertia(void) const;
	inline std::ostream& Joints(void) const;
	inline std::ostream& Forces(void) const;
	inline std::ostream& Beams(void) const;
	inline std::ostream& Rotors(void) const;
	inline std::ostream& Restart(void) const;
	inline std::ostream& RestartXSol(void) const;
	inline std::ostream& Aerodynamic(void) const;
	inline std::ostream& Hydraulic(void) const;
	inline std::ostream& PresNodes(void) const;
	inline std::ostream& Loadable(void) const;
	inline std::ostream& Genels(void) const;
	inline std::ostream& Partition(void) const;
	inline std::ostream& AeroModals(void) const;
	inline std::ostream& ReferenceFrames(void) const;
	inline std::ostream& Log(void) const;
	inline std::ostream& AirProps(void) const;
	inline std::ostream& Parameters(void) const;
	inline std::ostream& Externals(void) const;
	inline std::ostream& Modal(void) const;
	inline std::ostream& ThermalNodes(void) const;
	inline std::ostream& ThermalElements(void) const;
	inline std::ostream& Plates(void) const;
	inline std::ostream& Gravity(void) const;
	inline std::ostream& DofStats(void) const;
	inline std::ostream& DriveCallers(void) const;
	inline std::ostream& Traces(void) const;
	inline std::ostream& Eigenanalysis(void) const;

	inline int iW(void) const;
	inline int iP(void) const;

	void SetWidth(int iNewWidth);
	void SetPrecision(int iNewPrecision);
	void SetExceptions(std::ios::iostate flags);

#ifdef USE_NETCDF
	inline MBDynNcFile * pGetBinFile(void) const;

	struct AttrVal {
		std::string attr;
		std::string val;
		AttrVal(void) { NO_OP; };
		AttrVal(const std::string& attr, const std::string& val) : attr(attr), val(val) { NO_OP; };
	};

	typedef std::vector<OutputHandler::AttrVal> AttrValVec;
	typedef std::vector<MBDynNcDim> NcDimVec;

	MBDynNcDim 
	CreateDim(const std::string& name, integer size = -1);

	MBDynNcDim 
	GetDim(const std::string& name) const;

	inline MBDynNcDim DimTime(void) const;
	inline MBDynNcDim DimV1(void) const;
	inline MBDynNcDim DimV3(void) const;

	std::vector<size_t> ncStart1;
#if defined(USE_NETCDF4)
	std::vector<size_t> ncCount1;	
	std::vector<size_t> ncStart1x3;
	std::vector<size_t> ncCount1x3;
	std::vector<size_t> ncStart1x3x3;
	std::vector<size_t> ncCount1x3x3;
#endif  /* USE_NETCDF4 */

	MBDynNcVar
	CreateVar(const std::string& name, const MBDynNcType& type,
		const AttrValVec& attrs, const NcDimVec& dims);
	
	void
	WriteNcVar(const MBDynNcVar&, const Vec3&);

	void
	WriteNcVar(const MBDynNcVar&, const Mat3x3&);

	template <class Tvar>
	void
	WriteNcVar(const MBDynNcVar&, const Tvar&);
	
	MBDynNcVar
	CreateVar(const std::string& name, const std::string& type);

	template <class T>
	MBDynNcVar
	CreateVar(const std::string& name,
		const std::string& units, const std::string& description);

	MBDynNcVar
	CreateRotationVar(const std::string& name_prefix,
		const std::string& name_postfix,
		OrientationDescription od,
		const std::string& description);
#endif /* USE_NETCDF */
}; /* End class OutputHandler */

#ifdef USE_NETCDF
template <class T>
MBDynNcVar
OutputHandler::CreateVar(const std::string& name,
	const std::string& units, const std::string& description)
{
	AttrValVec attrs(3);
	NcDimVec dims(1);

	attrs[0] = AttrVal("units", units);
	attrs[2] = AttrVal("description", description);
	dims[0] = DimTime();
	
	if (typeid(T) == typeid(integer)) {
		attrs[1] = AttrVal("type", "integer");
		return CreateVar(name, MbNcInt, attrs, dims); // put this one inside if because couldn't figure out how to assign type inside if while declaring it outside of if.. alternative would be to create a pointer and change CreateVar function to copy type parameter instead of passing by reference (because type cannot be delete in CreateVar since addVar passes it by reference to netcdf)...
	} else if (typeid(T) == typeid(doublereal)) {
		attrs[1] = AttrVal("type", "doublereal");
	return CreateVar(name, MbNcDouble, attrs, dims); // see comment above
	} else if (typeid(T) == typeid(Vec3)) {
		attrs[1] = AttrVal("type", "Vec3");
		dims.resize(2);
		dims[1] = DimV3();
		return CreateVar(name, MbNcDouble, attrs, dims); // see comment above
	} else {
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
}

inline MBDynNcFile *
OutputHandler::pGetBinFile(void) const
{
	return m_pBinFile;
}

inline MBDynNcDim 
OutputHandler::DimTime(void) const
{
	return m_DimTime;
}

inline MBDynNcDim 
OutputHandler::DimV1(void) const
{
	return m_DimV1;
}

inline MBDynNcDim 
OutputHandler::DimV3(void) const
{
	return m_DimV3;
}
#endif /* USE_NETCDF */

inline std::ostream&
OutputHandler::Get(const OutputHandler::OutFiles f)
{
	ASSERT(f > -1 && f < LASTFILE);
	ASSERT(IsOpen(f));
	return *(OutData[f].pof);
}

inline std::ostream&
OutputHandler::Output(void) const
{
#ifdef DEBUG_COUT
	return const_cast<std::ostream &>(cout);
#else
	ASSERT(IsOpen(OUTPUT));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofOutput));
#endif
}

inline std::ostream&
OutputHandler::StrNodes(void) const
{
	ASSERT(IsOpen(STRNODES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofStrNodes));
}

inline std::ostream&
OutputHandler::Electric(void) const
{
	ASSERT(IsOpen(ELECTRIC));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofElectric));
}

inline std::ostream&
OutputHandler::ThermalNodes(void) const
{
	ASSERT(IsOpen(THERMALNODES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofThermalNodes));
}

inline std::ostream&
OutputHandler::ThermalElements(void) const
{
	ASSERT(IsOpen(THERMALELEMENTS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofThermalElements));
}

inline std::ostream&
OutputHandler::Abstract(void) const
{
	ASSERT(IsOpen(ABSTRACT));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofAbstract));
}

inline std::ostream&
OutputHandler::Inertia(void) const
{
	ASSERT(IsOpen(INERTIA));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofInertia));
}

inline std::ostream&
OutputHandler::Joints(void) const
{
	ASSERT(IsOpen(JOINTS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofJoints));
}

inline std::ostream&
OutputHandler::Forces(void) const
{
	ASSERT(IsOpen(FORCES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofForces));
}

inline std::ostream&
OutputHandler::Beams(void) const
{
	ASSERT(IsOpen(BEAMS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofBeams));
}

inline std::ostream&
OutputHandler::Rotors(void) const
{
	ASSERT(IsOpen(ROTORS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofRotors));
}

inline std::ostream&
OutputHandler::Restart(void) const
{
	ASSERT(IsOpen(RESTART));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofRestart));
}

inline std::ostream&
OutputHandler::RestartXSol(void) const
{
	ASSERT(IsOpen(RESTART));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofRestartXSol));
}

inline std::ostream&
OutputHandler::Aerodynamic(void) const
{
	ASSERT(IsOpen(AERODYNAMIC));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofAerodynamic));
}

inline std::ostream&
OutputHandler::Hydraulic(void) const
{
	ASSERT(IsOpen(HYDRAULIC));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofHydraulic));
}

inline std::ostream&
OutputHandler::PresNodes(void) const
{
	ASSERT(IsOpen(PRESNODES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofPresNodes));
}

inline std::ostream&
OutputHandler::Loadable(void) const
{
	ASSERT(IsOpen(LOADABLE));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofLoadable));
}

inline std::ostream&
OutputHandler::Genels(void) const
{
	ASSERT(IsOpen(GENELS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofGenels));
}

inline std::ostream&
OutputHandler::Partition(void) const
{
	ASSERT(IsOpen(PARTITION));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofPartition));
}

inline std::ostream&
OutputHandler::AeroModals(void) const
{
	ASSERT(IsOpen(AEROMODALS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofAeroModals));
}

inline std::ostream&
OutputHandler::ReferenceFrames(void) const
{
	ASSERT(IsOpen(REFERENCEFRAMES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofReferenceFrames));
}

inline std::ostream&
OutputHandler::Log(void) const
{
#ifdef DEBUG_COUT
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(cout));
#else
	ASSERT(IsOpen(LOG));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofLog));
#endif
}

inline std::ostream&
OutputHandler::AirProps(void) const
{
	ASSERT(IsOpen(AIRPROPS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofAirProps));
}

inline std::ostream&
OutputHandler::Parameters(void) const
{
	ASSERT(IsOpen(PARAMETERS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofParameters));
}

inline std::ostream&
OutputHandler::Externals(void) const
{
	ASSERT(IsOpen(EXTERNALS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofExternals));
}

inline std::ostream&
OutputHandler::Modal(void) const
{
	ASSERT(IsOpen(MODAL));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofModal));
}

inline std::ostream&
OutputHandler::Plates(void) const
{
	ASSERT(IsOpen(PLATES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofPlates));
}

inline std::ostream&
OutputHandler::Gravity(void) const
{
	ASSERT(IsOpen(GRAVITY));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofGravity));
}

inline std::ostream&
OutputHandler::DofStats(void) const
{
	ASSERT(IsOpen(DOFSTATS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofDofStats));
}

inline std::ostream&
OutputHandler::DriveCallers(void) const
{
	ASSERT(IsOpen(DRIVECALLERS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofDriveCallers));
}

inline std::ostream&
OutputHandler::Traces(void) const
{
	ASSERT(IsOpen(TRACES));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofTraces));
}

inline std::ostream&
OutputHandler::Eigenanalysis(void) const
{
	ASSERT(IsOpen(EIGENANALYSIS));
	return const_cast<std::ostream &>(dynamic_cast<const std::ostream &>(ofEigenanalysis));
}


inline int
OutputHandler::iW(void) const
{
	return iCurrWidth;
}

inline int
OutputHandler::iP(void) const
{
	return iCurrPrecision;
}

/* OutputHandler - end */


/* ToBeOutput - begin */

const flag fDefaultOut = 1;

class ToBeOutput {
public:
	enum {
		OUTPUT = 0x1U,

		// use OUTPUT_MASK to isolate public output flags
		OUTPUT_MASK = 0xFFU,

		// reserve values up to OUTPUT_PRIVATE for public use;
		// reserved output flags can start from OUTPUT_PRIVATE up
		OUTPUT_PRIVATE = 0x100U,

		// use OUTPUT_PRIVATE_MASK to isolate private output flags
		OUTPUT_PRIVATE_MASK = ~OUTPUT_MASK
	};

protected:
	flag fOutput;

public:
	ToBeOutput(flag fOut = fDefaultOut);
	virtual ~ToBeOutput(void);

	virtual void OutputPrepare(OutputHandler &OH);

	/* Regular output */
	virtual void Output(OutputHandler& OH) const;

	/* Output of perturbed solution (modes ...) */
	virtual void Output(OutputHandler& OH,
	const VectorHandler& X, const VectorHandler& XP) const;

	virtual flag fToBeOutput(void) const;
	virtual bool bToBeOutput(void) const;
	virtual void SetOutputFlag(flag f = flag(1));
};

/* ToBeOutput - end */

class Traceable {
public:
	enum {
		TRACE 				= 0x01U,
		TRACE_PUBLIC_MASK	= 0x0FU,
		TRACE_PRIVATE		= 0x10U,
		TRACE_PRIVATE_MASK	= ~TRACE_PUBLIC_MASK
	};

	Traceable(flag fTrace = 0);
	virtual ~Traceable(void);

	virtual void Trace(OutputHandler& OH) const=0;
	virtual flag fToBeTraced(void) const;
	virtual void SetTraceFlag(flag f = TRACE);

private:
	flag fTrace;
};

#endif /* OUTPUT_H */
