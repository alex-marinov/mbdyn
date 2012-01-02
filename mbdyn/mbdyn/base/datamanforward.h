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

#ifndef DATAMANFORWARD_H
#define DATAMANFORWARD_H

/* used by maps to compare strings case-insensitive */
struct ltstrcase {
	/* case-insensitive string comparison */
	bool operator()(const std::string& s1, const std::string& s2) const {
		return strcasecmp(s1.c_str(), s2.c_str()) < 0;
	};
};

class DataManagerErrors {
public:
	class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrAssemblyDiverged: public MBDynErrBase {
	public:
		ErrAssemblyDiverged(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrAssemblyMaxIters: public MBDynErrBase {
	public:
		ErrAssemblyMaxIters(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrElemNotAllowedInAssembly: MBDynErrBase {
	public:
		ErrElemNotAllowedInAssembly(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrUnknownElem: public MBDynErrBase {
	public:
		ErrUnknownElem(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrUnknownFunction: public MBDynErrBase {
	public:
		ErrUnknownFunction(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrUnknownNode: public MBDynErrBase {
	public:
		ErrUnknownNode(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrMissingNodes: public MBDynErrBase {
	public:
		ErrMissingNodes(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
	class ErrNeedDataManager: public MBDynErrBase {
	public:
		ErrNeedDataManager(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
};

class DataManager;


#endif /* DATAMANFORWARD_H */
