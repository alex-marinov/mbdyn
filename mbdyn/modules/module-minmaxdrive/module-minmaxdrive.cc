#ifdef HAVE_CONFIG_H
#include "mbconfig.h"
#endif

#include <vector>
#include <drive.h>

#include <myassert.h>
#include <except.h>

#include <strnode.h>
#include <elem.h>
#include <mynewmem.h>
#include <dataman.h>

#include "module-minmaxdrive.h"

class MinMaxDriveCaller : public DriveCaller
{
protected:
	MinMaxDriveCaller(const std::vector<DriveOwner>& drives);
public:
	virtual ~MinMaxDriveCaller();
	bool bIsDifferentiable(void) const;
	virtual std::ostream& Restart(std::ostream& out) const;
	doublereal dGet(void) const;
	doublereal dGetP(void) const;
	doublereal dGet(const doublereal& dVar) const;
	virtual doublereal dGetP(const doublereal& dVar) const;

protected:
	virtual bool bCompare(doublereal lhs, doublereal rhs) const = 0;
	virtual const char *c_str(void) const = 0;
	typedef std::vector<DriveOwner>::const_iterator iterator;
	const std::vector<DriveOwner> drives;
};

class MinDriveCaller : public MinMaxDriveCaller
{
public:
	MinDriveCaller(const std::vector<DriveOwner>& drives);
	virtual ~MinDriveCaller(void);
	virtual DriveCaller* pCopy(void) const;

protected:
	virtual bool bCompare(doublereal lhs, doublereal rhs) const;
	virtual const char *c_str(void) const;
};

class MaxDriveCaller : public MinMaxDriveCaller
{
public:
	MaxDriveCaller(const std::vector<DriveOwner>& drives);
	virtual ~MaxDriveCaller(void);
	virtual DriveCaller* pCopy(void) const;

protected:
	virtual bool bCompare(doublereal lhs, doublereal rhs) const;
	virtual const char *c_str(void) const;
};

struct MinMaxDriveDCR : public DriveCallerRead {
	const enum Type {
		MMD_MIN,
		MMD_MAX
	} eType;

	const char *c_str(void) const;

	explicit MinMaxDriveDCR(enum Type type);
	DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
};

MinMaxDriveCaller::MinMaxDriveCaller(const std::vector<DriveOwner>& drives)
: DriveCaller(0),
 drives(drives)
{
	NO_OP;
}

MinMaxDriveCaller::~MinMaxDriveCaller()
{
	NO_OP;
}

doublereal
MinMaxDriveCaller::dGet(void) const
{
	iterator i = drives.begin();
	doublereal dVal = i++->dGet();

	for ( ; i != drives.end(); ++i) {
		const doublereal dTmp = i->dGet();

		if (bCompare(dTmp, dVal)) {
			dVal = dTmp;
		}
	}

	return dVal;
}

doublereal MinMaxDriveCaller::dGetP(void) const
{
	iterator i = drives.begin();
	iterator iMin = i;
	doublereal dVal = i++->dGet();

	for ( ; i != drives.end(); ++i) {
		const doublereal dTmp = i->dGet();

		if (bCompare(dTmp, dVal)) {
			dVal = dTmp;
			iMin = i;
		}
	}

	return iMin->dGetP();
}

doublereal
MinMaxDriveCaller::dGet(const doublereal& dVar) const
{
	iterator i = drives.begin();
	doublereal dVal = i++->dGet(dVar);

	for ( ; i != drives.end(); ++i) {
		const doublereal dTmp = i->dGet(dVar);

		if (bCompare(dTmp, dVal)) {
			dVal = dTmp;
		}
	}

	return dVal;
}

doublereal MinMaxDriveCaller::dGetP(const doublereal& dVar) const
{
	iterator i = drives.begin();
	iterator iMin = i;
	doublereal dVal = i++->dGet(dVar);

	for ( ; i != drives.end(); ++i) {
		const doublereal dTmp = i->dGet(dVar);

		if (bCompare(dTmp, dVal)) {
			dVal = dTmp;
			iMin = i;
		}
	}

	return iMin->dGetP(dVar);
}

bool MinMaxDriveCaller::bIsDifferentiable(void) const
{
	for (iterator i = drives.begin(); i != drives.end(); ++i) {
		if (!i->bIsDifferentiable()) {
			return false;
		}
	}

	return true;
}

/* Restart */
std::ostream&
MinMaxDriveCaller::Restart(std::ostream& out) const
{
	out << c_str() << ", "
		<< drives.size() << ", ";

	for (iterator i = drives.begin(); i != drives.end(); ++i) {
		i->pGetDriveCaller()->Restart(out);

		if (drives.end() - i > 1) {
			out << ", ";
		}
	}

	return out;
}

MinDriveCaller::MinDriveCaller(const std::vector<DriveOwner>& drives)
: MinMaxDriveCaller(drives)
{
	NO_OP;
};

MinDriveCaller::~MinDriveCaller(void)
{
	NO_OP;
}

DriveCaller*
MinDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
			MinDriveCaller,
			MinDriveCaller(drives));

	return pDC;
}

bool MinDriveCaller::bCompare(doublereal lhs, doublereal rhs) const
{
	return lhs < rhs;
}

const char *MinDriveCaller::c_str(void) const
{
	return "min";
}

MaxDriveCaller::MaxDriveCaller(const std::vector<DriveOwner>& drives)
: MinMaxDriveCaller(drives)
{
	NO_OP;
};

MaxDriveCaller::~MaxDriveCaller(void)
{
	NO_OP;
}

DriveCaller*
MaxDriveCaller::pCopy(void) const
{
	DriveCaller* pDC = 0;

	SAFENEWWITHCONSTRUCTOR(pDC,
			MaxDriveCaller,
			MaxDriveCaller(drives));

	return pDC;
}

bool MaxDriveCaller::bCompare(doublereal lhs, doublereal rhs) const
{
	return lhs > rhs;
}

const char *MaxDriveCaller::c_str(void) const
{
	return "max";
}

MinMaxDriveDCR::MinMaxDriveDCR(enum Type type)
: eType(type)
{
	NO_OP;
}

const char *
MinMaxDriveDCR::c_str(void) const
{
	switch (eType) {
	case MMD_MIN:
		return "min";
	case MMD_MAX:
		return "max";
	}

	ASSERT(0);
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

DriveCaller *
MinMaxDriveDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
{
	DriveCaller *pDC = 0;

	const integer iNumDrives = HP.GetInt();

	if (iNumDrives < 1) {
		silent_cerr("at least one drive caller expected at line " << HP.GetLineData() << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	std::vector<DriveOwner> drives;

	drives.reserve(iNumDrives);

	for (int i = 0; i < iNumDrives; ++i) {
		drives.push_back(HP.GetDriveCaller());
	}

	switch (eType) {
	case MMD_MIN:
		SAFENEWWITHCONSTRUCTOR(pDC,
			MinDriveCaller,
			MinDriveCaller(drives));
		break;

	case MMD_MAX:
		SAFENEWWITHCONSTRUCTOR(pDC,
			MaxDriveCaller,
			MaxDriveCaller(drives));
		break;

	default:
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pDC;
}

bool
minmaxdrive_set()
{
	DriveCallerRead	*rf = new MinMaxDriveDCR(MinMaxDriveDCR::MMD_MIN);

	if (!SetDriveCallerData("min", rf)) {
		delete rf;
		return false;
	}

	rf = new MinMaxDriveDCR(MinMaxDriveDCR::MMD_MAX);

	if (!SetDriveCallerData("max", rf)) {
		delete rf;
		return false;
	}

	return true;
}

#ifndef STATIC_MODULES

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	if (!minmaxdrive_set()) {
		silent_cerr("minmaxdrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);
		return -1;
	}

	return 0;
}

#endif // ! STATIC_MODULES
