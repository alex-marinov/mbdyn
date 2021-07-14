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

/* socket driver */

#ifndef BUFMOD_H
#define BUFMOD_H

typedef std::map<std::string, bool> TypeMap_t;


/* BufCast - begin */

class BufCast {
protected:
	size_t m_offset;

public:
	BufCast(size_t offset);
	virtual ~BufCast(void);

	virtual size_t size(void) const = 0;
	virtual size_t offset(void) const = 0;
	virtual doublereal cast(const void *p) const = 0;
	virtual void uncast(void *pTo, doublereal d) const = 0;
	virtual BufCast *copy(size_t offset) const = 0;
};

template <class T>
class TBufCast : public BufCast {
public:
	TBufCast(size_t offset) : BufCast(offset) { NO_OP; };

	size_t size(void) const {
		return sizeof(T);
	};

	size_t offset(void) const {
		return m_offset;
	};

	// "cast" in the sense that whatever type comes from the stream,
	// it is cast into a doublereal
	doublereal cast(const void *pFrom) const {
		const char *p = &((const char *)pFrom)[m_offset];
		return doublereal(*((T *)p));
	};

	// "uncast" in the sense that a doublereal is transformed into
	// the correct type for the stream
	void uncast(void *pTo, doublereal d) const {
		char *p = &((char *)pTo)[m_offset];
		((T *)p)[0] = T(d);
	};

	BufCast *copy(size_t offset) const {
		return new TBufCast<T>(offset);
	};
};

template <typename T>
static T
mbswap(const T in)
{
	const char *pin = (const char *)&in;
	T out = T();
	char *pout = (char *)&out;

	for (unsigned int i = 0; i < sizeof(T)/2; i++) {
		pout[i] = pin[sizeof(T) - 1 - i];
		pout[sizeof(T) - 1 - i] = pin[i];
	}

	return out;
}

template <int8_t>
static int8_t
mbswap(const int8_t in)
{
	return in;
}

template <uint8_t>
static uint8_t
mbswap(const uint8_t in)
{
	return in;
}

template <class T>
class TBufCastHToN : public TBufCast<T> {
public:
	TBufCastHToN(size_t offset) : TBufCast<T>(offset) {};

	// "cast" in the sense that whatever type comes from the stream,
	// it is cast into a doublereal
	doublereal cast(const void *pFrom) const {
		const char *p = &((const char *)pFrom)[TBufCast<T>::m_offset];
		return mbswap<T>(*((T *)p));
	};

	// "uncast" in the sense that a doublereal is transformed into
	// the correct type for the stream
	void uncast(void *pTo, doublereal d) const {
		char *p = &((char *)pTo)[TBufCast<T>::m_offset];
		((T *)p)[0] = mbswap<T>(T(d));
	};

	BufCast *copy(size_t offset) const {
		return new TBufCastHToN<T>(offset);
	};
};

extern void
ReadBufCast(HighParser& HP, std::vector<BufCast *>& data);

extern bool
bIsLittleEndian(void);

extern void
SwapMapInit(TypeMap_t& swapmap);

/* BufCast - end */

#endif // BUFMOD_H
