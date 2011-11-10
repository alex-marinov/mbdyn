/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2011
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

#include "myassert.h"
#include "mynewmem.h"

#include "dataman.h"
#include "tpldrive_impl.h"
#include "reffrm.h"

/* ZeroTplDriveCaller<T> defined in "tpldrive_impl.h",
 * which includes "tpldrive.h" */

template <class T>
class NullTDCR : public TplDriveCallerRead<T> {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		TplDriveCaller<T>* pTplDC = 0;

		SAFENEW(pTplDC, ZeroTplDriveCaller<T>);

		return pTplDC;
	};
};

template <class T>
class ZeroTDCR : public NullTDCR<T> {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		silent_cerr("\"zero\" template drive caller "
			"at line " << HP.GetLineData() << " is deprecated; "
			"use \"null\" instead" << std::endl);

		return NullTDCR<T>::Read(pDM, HP);
	};
};

/* ZeroTplDriveCaller - end */


/* SingleTplDriveCaller - begin */

#if 0
template <class T>
class SingleTplDriveCaller : public TplDriveCaller<T>, public DriveOwner {
protected:
	T t;

public:
	SingleTplDriveCaller(const DriveCaller* pDC, const T& x)
	: DriveOwner(pDC), t(const_cast<T&>(x)) {
		NO_OP;
	};

	~SingleTplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		typedef SingleTplDriveCaller<T> dc;
		TplDriveCaller<T>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy(), t));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "single, ",
		Write(out, t, ", ") << ", ";
		return pGetDriveCaller()->Restart(out);
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		Write(out, t, ", ") << ", ";
		return pGetDriveCaller()->Restart(out);
	};

	inline T Get(void) const {
		return t*dGet();
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		return DriveOwner::bIsDifferentiable();
	};

	inline T GetP(void) const {
		return t*dGetP();
	};

	inline int getNDrives(void) const {
		return 1;
	};
};

/* Nota: in caso scalare, viene semplificata la classe in modo da
 *       usare solo il drive senza pesatura che viene assunta unitaria */

template<>
class SingleTplDriveCaller<doublereal>
: public TplDriveCaller<doublereal>, public DriveOwner {
public:
	SingleTplDriveCaller(const DriveCaller* pDC, const doublereal& = 0.)
	: DriveOwner(pDC) {
		NO_OP;
	};

	~SingleTplDriveCaller(void) {
		NO_OP;
	};

	/* copia */
	virtual TplDriveCaller<doublereal>* pCopy(void) const {
		TplDriveCaller<doublereal>* pDC = 0;

		typedef SingleTplDriveCaller<doublereal> dc;
		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(pGetDriveCaller()->pCopy()));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "single, ";
		return pGetDriveCaller()->Restart(out);
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		return pGetDriveCaller()->Restart(out);
	};

	inline doublereal Get(void) const {
		return dGet();
	};

	inline bool bIsDifferentiable(void) const {
		return DriveOwner::bIsDifferentiable();
	};

	inline doublereal GetP(void) const {
		return dGetP();
	};

	inline int getNDrives(void) const {
		return 1;
	};
};
#endif

template <class T>
class SingleTDCR : public TplDriveCallerRead<T> {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		T t(mb_zero<T>());

		t = GetT(HP, t);

		DriveCaller* pDC = HP.GetDriveCaller();

		TplDriveCaller<T>* pTplDC = 0;

		SAFENEWWITHCONSTRUCTOR(pTplDC,
			SingleTplDriveCaller<T>,
			SingleTplDriveCaller<T>(pDC, t));

		return pTplDC;
	};
};

/* SingleTplDriveCaller - end */


/* CompTplDriveCaller - begin */

template <class T>
class CompTplDriveCaller : public TplDriveCaller<T>, public DriveOwner {
protected:
	std::vector<DriveCaller *> m_dc;

public:
	CompTplDriveCaller(std::vector<DriveCaller *>& dc)
	: m_dc(dc)
	{
		if (typeid(T) == typeid(Vec3)) {
			if (dc.size() != 3) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else if (typeid(T) == typeid(Vec6)) {
			if (dc.size() != 6) {
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}

		} else {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	~CompTplDriveCaller(void) {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			delete m_dc[i];
		}
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		typedef CompTplDriveCaller<T> dc;
		TplDriveCaller<T>* pDC = 0;

		std::vector<DriveCaller *> tmpdc(m_dc.size());

		for (unsigned i = 0; i < m_dc.size(); i++) {
			tmpdc[i] = m_dc[i]->pCopy();
		}

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(tmpdc));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "component";

		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	inline T Get(const doublereal& dVar) const {
		T t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i + 1) = m_dc[i]->dGet(dVar);
		}

		return t;
	};

	inline T Get(void) const {
		T t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i + 1) = m_dc[i]->dGet();
		}

		return t;
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			if (!m_dc[i]->bIsDifferentiable()) {
				return false;
			}
		}

		return true;
	};

	inline T GetP(void) const {
		T t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i + 1) = m_dc[i]->dGetP();
		}

		return t;
	};

	inline int getNDrives(void) const {
		return m_dc.size();
	};
};

template <>
class CompTplDriveCaller<Mat3x3> : public TplDriveCaller<Mat3x3>, public DriveOwner {
protected:
	std::vector<DriveCaller *> m_dc;

public:
	CompTplDriveCaller(std::vector<DriveCaller *>& dc)
	: m_dc(dc)
	{
		if (dc.size() != 9) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	~CompTplDriveCaller(void) {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			delete m_dc[i];
		}
	};

	/* copia */
	virtual TplDriveCaller<Mat3x3>* pCopy(void) const {
		typedef CompTplDriveCaller<Mat3x3> dc;
		TplDriveCaller<Mat3x3>* pDC = 0;

		std::vector<DriveCaller *> tmpdc(m_dc.size());

		for (unsigned i = 0; i < m_dc.size(); i++) {
			tmpdc[i] = m_dc[i]->pCopy();
		}

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(tmpdc));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "component";

		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	inline Mat3x3 Get(const doublereal& dVar) const {
		Mat3x3 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/3 + 1, i%3 + 1) = m_dc[i]->dGet(dVar);
		}

		return t;
	};

	inline Mat3x3 Get(void) const {
		Mat3x3 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/3 + 1, i%3 + 1) = m_dc[i]->dGet();
		}

		return t;
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			if (!m_dc[i]->bIsDifferentiable()) {
				return false;
			}
		}

		return true;
	};

	inline Mat3x3 GetP(void) const {
		Mat3x3 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/3 + 1, i%3 + 1) = m_dc[i]->dGetP();
		}

		return t;
	};

	inline int getNDrives(void) const {
		return m_dc.size();
	};
};

template <>
class CompTplDriveCaller<Mat6x6> : public TplDriveCaller<Mat6x6>, public DriveOwner {
protected:
	std::vector<DriveCaller *> m_dc;

public:
	CompTplDriveCaller(std::vector<DriveCaller *>& dc)
	: m_dc(dc)
	{
		if (dc.size() != 36) {
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	};

	~CompTplDriveCaller(void) {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			delete m_dc[i];
		}
	};

	/* copia */
	virtual TplDriveCaller<Mat6x6>* pCopy(void) const {
		typedef CompTplDriveCaller<Mat6x6> dc;
		TplDriveCaller<Mat6x6>* pDC = 0;

		std::vector<DriveCaller *> tmpdc(m_dc.size());

		for (unsigned i = 0; i < m_dc.size(); i++) {
			tmpdc[i] = m_dc[i]->pCopy();
		}

		SAFENEWWITHCONSTRUCTOR(pDC, dc, dc(tmpdc));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "component";

		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i]->Restart(out);
		}

		return out;
	};

	inline Mat6x6 Get(const doublereal& dVar) const {
		Mat6x6 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/6 + 1, i%6 + 1) = m_dc[i]->dGet(dVar);
		}

		return t;
	};

	inline Mat6x6 Get(void) const {
		Mat6x6 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/6 + 1, i%6 + 1) = m_dc[i]->dGet();
		}

		return t;
	};

	/* this is about drives that are differentiable */
	inline bool bIsDifferentiable(void) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			if (!m_dc[i]->bIsDifferentiable()) {
				return false;
			}
		}

		return true;
	};

	inline Mat6x6 GetP(void) const {
		Mat6x6 t;

		for (unsigned i = 0; i < m_dc.size(); i++) {
			t(i/6 + 1, i%6 + 1) = m_dc[i]->dGetP();
		}

		return t;
	};

	inline int getNDrives(void) const {
		return m_dc.size();
	};
};

template <class T>
class CompTDCR : public TplDriveCallerRead<T> {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		std::vector<DriveCaller *> dc;
		unsigned nr = 0, nc = 1;
		if (typeid(T) == typeid(Vec3)) {
			nr = 3;

		} else if (typeid(T) == typeid(Vec6)) {
			nr = 6;

		} else if (typeid(T) == typeid(Mat3x3)) {
			nr = 3;
			nc = 3;

		} else if (typeid(T) == typeid(Mat6x6)) {
			nr = 6;
			nc = 6;

		} else {
			silent_cerr("component template drive used with unknown type" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		dc.resize(nr*nc);

		if (nc > 1 && HP.IsKeyWord("sym")) {
			for (unsigned ir = 0; ir < nr; ir++) {
				for (unsigned ic = ir; ic < nc; ic++) {
					if (HP.IsKeyWord("inactive")) {
						dc[nc*ir + ic] = 0;
						SAFENEW(dc[nc*ir + ic], NullDriveCaller);

					} else {
						dc[nc*ir + ic] = HP.GetDriveCaller();
					}

					if (ic > ir) {
						dc[nc*ic + ir] = dc[nc*ir + ic]->pCopy();
					}
				}
			}

		} else if (nc > 1 && HP.IsKeyWord("diag")) {
			for (unsigned ir = 0; ir < nr; ir++) {
				if (HP.IsKeyWord("inactive")) {
					dc[nc*ir + ir] = 0;
					SAFENEW(dc[nc*ir + ir], NullDriveCaller);

				} else {
					dc[nc*ir + ir] = HP.GetDriveCaller();
				}

				for (unsigned ic = ir + 1; ic < nc; ic++) {
					dc[nc*ir + ic] = 0;
					SAFENEW(dc[nc*ir + ic], NullDriveCaller);
					dc[nc*ic + ir] = dc[nc*ir + ic]->pCopy();
				}
			}

		} else {
			for (unsigned i = 0; i < dc.size(); i++) {
				if (HP.IsKeyWord("inactive")) {
					dc[i] = 0;
					SAFENEW(dc[i], NullDriveCaller);

				} else {
					dc[i] = HP.GetDriveCaller();
				}
			}
		}

		TplDriveCaller<T>* pTplDC = 0;

		SAFENEWWITHCONSTRUCTOR(pTplDC,
			CompTplDriveCaller<T>,
			CompTplDriveCaller<T>(dc));

		return pTplDC;
	};
};

/* CompTplDriveCaller - end */

/* ArrayTplDriveCaller - begin */

template <class T>
struct DrivesArray {
	DriveCaller* pDriveCaller;
	T t;
	DrivesArray(void) : pDriveCaller(0), t() {};
};

template <class T>
class ArrayTplDriveCaller : public TplDriveCaller<T> {
protected:
	std::vector<DrivesArray<T> > m_dc;

public:
	ArrayTplDriveCaller(std::vector<DrivesArray<T> >& dc)
	: m_dc(dc) {
		ASSERT(!dc.empty());
	};

	~ArrayTplDriveCaller(void) {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			SAFEDELETE(m_dc[i].pDriveCaller);
		}
	};

	/* copia */
	virtual TplDriveCaller<T>* pCopy(void) const {
		std::vector<DrivesArray<T> > dc(m_dc.size());
		for (unsigned i = 0; i < m_dc.size(); i++) {
			dc[i].pDriveCaller = m_dc[i].pDriveCaller->pCopy();
			dc[i].t = m_dc[i].t;
		}

		typedef ArrayTplDriveCaller<T> dc_t;
		TplDriveCaller<T>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc_t, dc_t(dc));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "array, " << m_dc.size();
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ",
				Write(out, m_dc[i].t, ", ") << ", ",
				m_dc[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ",
				Write(out, m_dc[i].t, ", ") << ", ",
				m_dc[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	inline T Get(const doublereal& dVar) const {
		T v = mb_zero<T>();
		for (unsigned i = 0; i < m_dc.size(); i++) {
			v += (m_dc[i].t)*(m_dc[i].pDriveCaller->dGet(dVar));
		}
		return v;
	};

	inline T Get(void) const {
		T v = mb_zero<T>();
		for (unsigned i = 0; i < m_dc.size(); i++) {
			v += (m_dc[i].t)*(m_dc[i].pDriveCaller->dGet());
		}
		return v;
	};

	inline bool bIsDifferentiable(void) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			if (!m_dc[i].pDriveCaller->bIsDifferentiable()) {
				return false;
			}
		}
		return true;
	};

	inline T GetP(void) const {
		T v = mb_zero<T>();
		for (unsigned i = 0; i < m_dc.size(); i++) {
			v += (m_dc[i].t)*(m_dc[i].pDriveCaller->dGetP());
		}
		return v;
	};

	inline int getNDrives(void) const {
		return m_dc.size();
	};
};

template<>
class ArrayTplDriveCaller<doublereal> : public TplDriveCaller<doublereal> {
protected:
	std::vector<DrivesArray<doublereal> > m_dc;

public:
	ArrayTplDriveCaller(std::vector<DrivesArray<doublereal> > dc)
	: m_dc(dc) {
		ASSERT(!m_dc.empty());
	};

	virtual ~ArrayTplDriveCaller(void) {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			SAFEDELETE(m_dc[i].pDriveCaller);
		}
	};

	/* copia */
	virtual TplDriveCaller<doublereal>* pCopy(void) const {
		std::vector<DrivesArray<doublereal> > dc(m_dc.size());
		for (unsigned i = 0; i < m_dc.size(); i++) {
			dc[i].pDriveCaller = m_dc[i].pDriveCaller->pCopy();
			dc[i].t = m_dc[i].t;
		}

		typedef ArrayTplDriveCaller<doublereal> dc_t;
		TplDriveCaller<doublereal>* pDC = 0;

		SAFENEWWITHCONSTRUCTOR(pDC, dc_t, dc_t(dc));

		return pDC;
	};

	/* Scrive il contributo del DriveCaller al file di restart */
	virtual std::ostream& Restart(std::ostream& out) const {
		out << "array, " << m_dc.size();
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	virtual std::ostream& Restart_int(std::ostream& out) const {
		for (unsigned i = 0; i < m_dc.size(); i++) {
			out << ", ", m_dc[i].pDriveCaller->Restart(out);
		}
		return out;
	};

	inline doublereal Get(const doublereal& dVar) const {
		doublereal v = 0.;
		for (unsigned i = 0; i < m_dc.size(); i++) {
			v += m_dc[i].pDriveCaller->dGet(dVar);
		}
		return v;
	};

	inline doublereal Get(void) const {
		doublereal v = 0.;
		for (unsigned i = 0; i < m_dc.size(); i++) {
			v += m_dc[i].pDriveCaller->dGet();
		}
		return v;
	};

	inline int getNDrives(void) const {
		return m_dc.size();
	};
};

/* Nota: di questa classe non viene scritta esplicitamente la versione
 *       per reali in quanto il coefficiente moltiplicativo
 *       puo' essere usato per un'ulteriore pesatura del drive */

template <class T>
class ArrayTDCR : public SingleTDCR<T> {
public:
	virtual TplDriveCaller<T> *
	Read(const DataManager* pDM, MBDynParser& HP) {
		unsigned short int iNumDr = HP.GetInt();
		if (iNumDr == 0) {
			silent_cerr("At least one drive is required "
				"in array template drive" << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);

		} else if (iNumDr == 1) {
			return SingleTDCR<T>::Read(pDM, HP);
		} /* else */

		std::vector<DrivesArray<T> > dc(iNumDr);

		for (unsigned short int i = 0; i < iNumDr; i++) {
			T t(mb_zero<T>());
			dc[i].t = GetT(HP, t);
			dc[i].pDriveCaller = HP.GetDriveCaller();
		}

		TplDriveCaller<T>* pTplDC = 0;

		SAFENEWWITHCONSTRUCTOR(pTplDC,
			ArrayTplDriveCaller<T>,
			ArrayTplDriveCaller<T>(dc));

		return pTplDC;
	};
};

/* ArrayTplDriveCaller - end */


extern doublereal GetT(MBDynParser& HP, const doublereal& t);

template <class T> T GetT(MBDynParser& HP, const T& t)
{
	return HP.Get(t);
}

/* template drive caller containers */
typedef std::map<std::string, TplDriveCallerRead<doublereal> *, ltstrcase> DC1DFuncMapType;
typedef std::map<std::string, TplDriveCallerRead<Vec3> *, ltstrcase> DC3DFuncMapType;
typedef std::map<std::string, TplDriveCallerRead<Vec6> *, ltstrcase> DC6DFuncMapType;

typedef std::map<std::string, TplDriveCallerRead<Mat3x3> *, ltstrcase> DC3x3DFuncMapType;
typedef std::map<std::string, TplDriveCallerRead<Mat6x6> *, ltstrcase> DC6x6DFuncMapType;

static DC1DFuncMapType DC1DFuncMap;
static DC3DFuncMapType DC3DFuncMap;
static DC6DFuncMapType DC6DFuncMap;

static DC3x3DFuncMapType DC3x3DFuncMap;
static DC6x6DFuncMapType DC6x6DFuncMap;

struct DC1DWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::DC1DFuncMap.find(std::string(s)) != ::DC1DFuncMap.end();
	};
};

struct DC3DWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::DC3DFuncMap.find(std::string(s)) != ::DC3DFuncMap.end();
	};
};

struct DC6DWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::DC6DFuncMap.find(std::string(s)) != ::DC6DFuncMap.end();
	};
};

struct DC3x3DWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::DC3x3DFuncMap.find(std::string(s)) != ::DC3x3DFuncMap.end();
	};
};

struct DC6x6DWordSetType : public HighParser::WordSet {
	bool IsWord(const std::string& s) const {
		return ::DC6x6DFuncMap.find(std::string(s)) != ::DC6x6DFuncMap.end();
	};
};

static DC1DWordSetType DC1DWordSet;
static DC3DWordSetType DC3DWordSet;
static DC6DWordSetType DC6DWordSet;

static DC3x3DWordSetType DC3x3DWordSet;
static DC6x6DWordSetType DC6x6DWordSet;

/* template drive caller registration functions: call to register one */
bool
SetDC1D(const char *name, TplDriveCallerRead<doublereal> *rf)
{
	pedantic_cout("registering template drive caller 1D \"" << name << "\""
		<< std::endl );
	return DC1DFuncMap.insert(DC1DFuncMapType::value_type(name, rf)).second;
}

bool
SetDC3D(const char *name, TplDriveCallerRead<Vec3> *rf)
{
	pedantic_cout("registering template drive caller 3D \"" << name << "\""
		<< std::endl );
	return DC3DFuncMap.insert(DC3DFuncMapType::value_type(name, rf)).second;
}

bool
SetDC6D(const char *name, TplDriveCallerRead<Vec6> *rf)
{
	pedantic_cout("registering template drive caller 6D \"" << name << "\""
		<< std::endl );
	return DC6DFuncMap.insert(DC6DFuncMapType::value_type(name, rf)).second;
}

bool
SetDC3x3D(const char *name, TplDriveCallerRead<Mat3x3> *rf)
{
	pedantic_cout("registering template drive caller 3x3D \"" << name << "\""
		<< std::endl );
	return DC3x3DFuncMap.insert(DC3x3DFuncMapType::value_type(name, rf)).second;
}

bool
SetDC6x6D(const char *name, TplDriveCallerRead<Mat6x6> *rf)
{
	pedantic_cout("registering template drive caller 6x6D \"" << name << "\""
		<< std::endl );
	return DC6x6DFuncMap.insert(DC6x6DFuncMapType::value_type(name, rf)).second;
}

/* functions that read a template drive caller */
TplDriveCaller<doublereal> *
ReadDC1D(const DataManager* pDM, MBDynParser& HP)
{
	const char *s = HP.IsWord(DC1DWordSet);
	if (s == 0) {
		s = "single";
	}

	DC1DFuncMapType::iterator func = DC1DFuncMap.find(std::string(s));
	if (func == DC1DFuncMap.end()) {
		silent_cerr("unknown template drive caller 1D type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

TplDriveCaller<Vec3> *
ReadDC3D(const DataManager* pDM, MBDynParser& HP)
{
	const char *s = HP.IsWord(DC3DWordSet);
	if (s == 0) {
		s = "single";
	}

	DC3DFuncMapType::iterator func = DC3DFuncMap.find(std::string(s));
	if (func == DC3DFuncMap.end()) {
		silent_cerr("unknown template drive caller 3D type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

TplDriveCaller<Vec6> *
ReadDC6D(const DataManager* pDM, MBDynParser& HP)
{
	const char *s = HP.IsWord(DC6DWordSet);
	if (s == 0) {
		s = "single";
	}

	DC6DFuncMapType::iterator func = DC6DFuncMap.find(std::string(s));
	if (func == DC6DFuncMap.end()) {
		silent_cerr("unknown template drive caller 6D type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

TplDriveCaller<Mat3x3> *
ReadDC3x3D(const DataManager* pDM, MBDynParser& HP)
{
	const char *s = HP.IsWord(DC3x3DWordSet);
	if (s == 0) {
		s = "single";
	}

	DC3x3DFuncMapType::iterator func = DC3x3DFuncMap.find(std::string(s));
	if (func == DC3x3DFuncMap.end()) {
		silent_cerr("unknown template drive caller 3x3D type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

TplDriveCaller<Mat6x6> *
ReadDC6x6D(const DataManager* pDM, MBDynParser& HP)
{
	const char *s = HP.IsWord(DC6x6DWordSet);
	if (s == 0) {
		s = "single";
	}

	DC6x6DFuncMapType::iterator func = DC6x6DFuncMap.find(std::string(s));
	if (func == DC6x6DFuncMap.end()) {
		silent_cerr("unknown template drive caller 6x6D type \"" << s << "\" "
			"at line " << HP.GetLineData() << std::endl);
		throw DataManager::ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return func->second->Read(pDM, HP);
}

TplDriveCaller<Vec3> *
ReadDCVecRel(const DataManager* pDM, MBDynParser& HP, const ReferenceFrame& rf)
{
	MBDynParser::VecRelManip manip(HP, rf);

	HP.PushManip(&manip);
	TplDriveCaller<Vec3>* pDC = ReadDC3D(pDM, HP);
	if (HP.GetManip() == &manip) {
		HP.PopManip();
	}

	return pDC;
}

TplDriveCaller<Vec3> *
ReadDCVecAbs(const DataManager* pDM, MBDynParser& HP, const ReferenceFrame& rf)
{
	MBDynParser::VecAbsManip manip(HP, rf);

	HP.PushManip(&manip);
	TplDriveCaller<Vec3>* pDC = ReadDC3D(pDM, HP);
	if (HP.GetManip() == &manip) {
		HP.PopManip();
	}

	return pDC;
}

static unsigned done;

void
InitTplDC(void)
{
	if (::done++ > 0) {
		return;
	}

	/* null */
	SetDC1D("null", new NullTDCR<doublereal>);
	SetDC3D("null", new NullTDCR<Vec3>);
	SetDC6D("null", new NullTDCR<Vec6>);

	SetDC3x3D("null", new NullTDCR<Mat3x3>);
	SetDC6x6D("null", new NullTDCR<Mat6x6>);

	/* zero (deprecated) */
	SetDC1D("zero", new ZeroTDCR<doublereal>);
	SetDC3D("zero", new ZeroTDCR<Vec3>);
	SetDC6D("zero", new ZeroTDCR<Vec6>);

	SetDC3x3D("zero", new ZeroTDCR<Mat3x3>);
	SetDC6x6D("zero", new ZeroTDCR<Mat6x6>);

	/* single */
	SetDC1D("single", new SingleTDCR<doublereal>);
	SetDC3D("single", new SingleTDCR<Vec3>);
	SetDC6D("single", new SingleTDCR<Vec6>);

	SetDC3x3D("single", new SingleTDCR<Mat3x3>);
	SetDC6x6D("single", new SingleTDCR<Mat6x6>);

	/* component */
	/* in the scalar case, "single" and "component" are identical */
	SetDC1D("component", new SingleTDCR<doublereal>);
	SetDC3D("component", new CompTDCR<Vec3>);
	SetDC6D("component", new CompTDCR<Vec6>);

	SetDC3x3D("component", new CompTDCR<Mat3x3>);
	SetDC6x6D("component", new CompTDCR<Mat6x6>);

	/* array */
	SetDC1D("array", new ArrayTDCR<doublereal>);
	SetDC3D("array", new ArrayTDCR<Vec3>);
	SetDC6D("array", new ArrayTDCR<Vec6>);

	SetDC3x3D("array", new ArrayTDCR<Mat3x3>);
	SetDC6x6D("array", new ArrayTDCR<Mat6x6>);
}

void
DestroyTplDC(void)
{
	if (::done == 0) {
		silent_cerr("DestroyTplDC() called once too many" << std::endl);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	if (--::done > 0) {
		return;
	}

	/* free stuff */
	for (DC1DFuncMapType::iterator i = DC1DFuncMap.begin(); i != DC1DFuncMap.end(); ++i) {
		delete i->second;
	}
	DC1DFuncMap.clear();

	for (DC3DFuncMapType::iterator i = DC3DFuncMap.begin(); i != DC3DFuncMap.end(); ++i) {
		delete i->second;
	}
	DC3DFuncMap.clear();

	for (DC6DFuncMapType::iterator i = DC6DFuncMap.begin(); i != DC6DFuncMap.end(); ++i) {
		delete i->second;
	}
	DC6DFuncMap.clear();

	for (DC3x3DFuncMapType::iterator i = DC3x3DFuncMap.begin(); i != DC3x3DFuncMap.end(); ++i) {
		delete i->second;
	}
	DC3x3DFuncMap.clear();

	for (DC6x6DFuncMapType::iterator i = DC6x6DFuncMap.begin(); i != DC6x6DFuncMap.end(); ++i) {
		delete i->second;
	}
	DC6x6DFuncMap.clear();
}

