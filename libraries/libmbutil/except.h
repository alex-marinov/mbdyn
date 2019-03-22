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

#ifndef EXCEPT_H
#define EXCEPT_H

#include <exception>
#include <stdio.h>

#include <iostream>

// Use this macro instead of the required set of args when declaring the constructor
// of classes derived from Err_base
#define MBDYN_EXCEPT_ARGS_DECL_NOOPT \
	const char *file, int line, const char *func
#define MBDYN_EXCEPT_ARGS_DECL \
	MBDYN_EXCEPT_ARGS_DECL_NOOPT , const std::string r = std::string()
// Use this macro instead of the required set of args when coding, not inline, the constructor
// of classes derived from Err_base
#define MBDYN_EXCEPT_ARGS_DECL_NOOPT_NODEF \
	const char *file, int line, const char *func
#define MBDYN_EXCEPT_ARGS_DECL_NODEF \
	MBDYN_EXCEPT_ARGS_DECL_NOOPT_NODEF , const std::string r
// Use this macro to pass the required set of args thru to Err_base from the constructor
// of derived classes
#define MBDYN_EXCEPT_ARGS_NOOPT_PASSTHRU \
	file, line, func
#define MBDYN_EXCEPT_ARGS_PASSTHRU \
	MBDYN_EXCEPT_ARGS_NOOPT_PASSTHRU , r
// Use this macro to pass the required set of args to error classes derived from Err_base
#if __GNUC__ >= 2
#define MBDYN_EXCEPT_ARGS \
	__FILE__ , __LINE__ , __PRETTY_FUNCTION__
#else // ! __GNUC__
// FIXME: need to detect whether __func__ (C99) is available
#define MBDYN_EXCEPT_ARGS \
	__FILE__ , __LINE__ , "(unknown)"
#endif // ! __GNUC__

class MBDynErrBase : public std::exception {
private:
	std::string s;

public:
	MBDynErrBase(MBDYN_EXCEPT_ARGS_DECL);
	virtual ~MBDynErrBase(void) noexcept {};
	void Set(const std::string& s);
	const char * what(void) const noexcept;
};


class NoErr : public MBDynErrBase {
public:
	NoErr(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrGeneric : public MBDynErrBase {
public:
	ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrInterrupted : public MBDynErrBase {
public:
	ErrInterrupted(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
  
class ErrOutOfRange : public MBDynErrBase {
public:
	ErrOutOfRange(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrIndexOutOfRange : public ErrOutOfRange {
protected:
	ErrIndexOutOfRange(const char *idx_type, int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL);
	void WriteMsg(const char *idx_type, int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL);
public:
	ErrIndexOutOfRange(int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL);
};
class ErrRowIndexOutOfRange : public ErrIndexOutOfRange {
public:
	ErrRowIndexOutOfRange(int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL)
		: ErrIndexOutOfRange("row ", idx, imin, imax, MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrColIndexOutOfRange : public ErrIndexOutOfRange {
public:
	ErrColIndexOutOfRange(int idx, int imin, int imax, MBDYN_EXCEPT_ARGS_DECL)
		: ErrIndexOutOfRange("col ", idx, imin, imax, MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrDivideByZero : public MBDynErrBase {
public:
	ErrDivideByZero(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrNullNorm : public MBDynErrBase {
public:
	ErrNullNorm(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrMemory : public MBDynErrBase {
public: 
	ErrMemory(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
  
class EndOfFile : public MBDynErrBase {
public:
	EndOfFile(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrFile : public MBDynErrBase {
public:
	ErrFile(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrFileSystem : public MBDynErrBase {
public:
	ErrFileSystem(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
  
class ErrNotAvailableYet : public MBDynErrBase {
public:
	ErrNotAvailableYet(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
class ErrNotImplementedYet : public MBDynErrBase {
public:
	ErrNotImplementedYet(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
};
  
#endif /* EXCEPT_H */

