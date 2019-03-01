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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>

#include "dataman.h"
#include "drive.h"

#include <boost/random.hpp>

static boost::mt19937 gen;
static std::string sInSeedFileName;
static std::string sOutSeedFileName;

class BoostRandomDriveCaller : public DriveCaller {
private:
	boost::normal_distribution<double> m_dist;
	mutable boost::variate_generator<boost::mt19937&, 
		boost::normal_distribution<double> > m_var_nor;

	mutable double m_dCurrTime;
	mutable double m_dCurrValue;

public:
	BoostRandomDriveCaller(const DriveHandler* pDH,
		doublereal dMean, doublereal dVariance);
	virtual ~BoostRandomDriveCaller(void);
 
	/* Copia */
	virtual DriveCaller* pCopy(void) const;
 
	/* Scrive il contributo del DriveCaller al file di restart */   
	virtual std::ostream& Restart(std::ostream& out) const;
 
	inline doublereal dGet(const doublereal& dVar) const;

	/* this is about drives that are differentiable */
	virtual bool bIsDifferentiable(void) const;
};

BoostRandomDriveCaller::BoostRandomDriveCaller(const DriveHandler* pDH,
	doublereal dMean, doublereal dVariance)
: DriveCaller(pDH),
m_dist(dMean, dVariance),
m_var_nor(::gen, m_dist),
m_dCurrTime(pDH->dGetTime()),
m_dCurrValue(m_var_nor())
{
	NO_OP;
}

BoostRandomDriveCaller::~BoostRandomDriveCaller(void)
{
	NO_OP;
}

DriveCaller *
BoostRandomDriveCaller::pCopy(void) const
{
	return new BoostRandomDriveCaller(pDrvHdl,
		m_dist.mean(), m_dist.sigma());
}

std::ostream&
BoostRandomDriveCaller::Restart(std::ostream& out) const
{
	return out << "boost random, "
		<< m_dist.mean() << ", "
		<< m_dist.sigma();
}

inline doublereal 
BoostRandomDriveCaller::dGet(const doublereal& dVar) const 
{
	if (dVar > m_dCurrTime) {
		m_dCurrTime = dVar;
		m_dCurrValue = m_var_nor();
	}

	return m_dCurrValue;
}

inline bool
BoostRandomDriveCaller::bIsDifferentiable(void) const
{
	return false;
}

/* prototype of the functional object: reads a drive caller */
struct BoostRandomDCR : public DriveCallerRead {
	virtual DriveCaller *
	Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred) {
		doublereal dMean = HP.GetReal();
		doublereal dVariance = HP.GetReal();
		return new BoostRandomDriveCaller(pDM->pGetDrvHdl(),
			dMean, dVariance);
	};

	virtual ~BoostRandomDCR(void) {
		if (!::sOutSeedFileName.empty()) {
			std::ofstream out(::sOutSeedFileName.c_str());
			if (out) {
				out << ::gen;
			}
		}
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
#endif
	MBDynParser	*pHP = (MBDynParser *)php;

	bool bGotSeed(false);
	bool bGotSeedInputFileName(false);
	int iSeed;

	while (pHP->IsArg()) {
		if (pHP->IsKeyWord("seed")) {
			if (bGotSeedInputFileName) {
				silent_cerr("BoostRandom::module_init(): previously parsed \"seed input file name\" incompatible with \"seed\" at line " << pHP->GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			iSeed = pHP->GetInt();
			bGotSeed = true;

		} else if (pHP->IsKeyWord("seed" "input" "file" "name")) {
			if (bGotSeed) {
				silent_cerr("BoostRandom::module_init(): previously parsed \"seed\" incompatible with \"seed input file name\" at line " << pHP->GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

			const char *sTmp = pHP->GetFileName();
			if (sTmp == 0) {
				silent_cerr("BoostRandom::module_init(): unable to parse \"seed input file name\" at line " << pHP->GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			::sInSeedFileName = sTmp;
			bGotSeedInputFileName = true;

		} else if (pHP->IsKeyWord("seed" "output" "file" "name")) {
			const char *sTmp = pHP->GetFileName();
			if (sTmp == 0) {
				silent_cerr("BoostRandom::module_init(): unable to parse \"seed output file name\" at line " << pHP->GetLineData() << std::endl);
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			::sOutSeedFileName = sTmp;

		} else {
			silent_cerr("BoostRandom::module_init(): unknown arg at line " << pHP->GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	}

	DriveCallerRead	*rf = new BoostRandomDCR;

	if (!SetDriveCallerData("boost" "random", rf)) {
		delete rf;

		silent_cerr("BoostRandomDrive: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	if (!::sInSeedFileName.empty()) {
		std::ifstream in(::sInSeedFileName.c_str());
		if (in) {
			in >> ::gen;
		}

	} else if (bGotSeed) {
		::gen.seed(iSeed);
	}

	return 0;
}

