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

//#include "unistd.h"
//#include "ac/f2c.h"

#include "parser.h"
#include "bufmod.h"

// FIXME: conditional build?
#include <arpa/inet.h>
// Some systems require the inclusion of <netinet/in.h> instead of <arpa/inet.h>.
// #include <netinet/in.h>

BufCast::BufCast(size_t offset)
: m_offset(offset)
{
	NO_OP;
}

BufCast::~BufCast(void)
{
	NO_OP;
}

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
	T out;
	char *pout = (char *)&out;

	for (int i = 0; i < sizeof(T)/2; i++) {
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

typedef std::map<std::string, bool> TypeMap_t;

static BufCast *
ReadOneBufCast(HighParser& HP, size_t& offset, TypeMap_t& swapmap, bool bNoSkip = false)
{
	BufCast *pBC(0);

	if (HP.IsKeyWord("int8_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(int8_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<int8_t>(offset);
			
		} else {
			pBC = new TBufCast<int8_t>(offset);
		}
		offset += sizeof(int8_t);

	} else if (HP.IsKeyWord("uint8_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(uint8_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<uint8_t>(offset);
			
		} else {
			pBC = new TBufCast<uint8_t>(offset);
		}
		offset += sizeof(uint8_t);

	} else if (HP.IsKeyWord("int16_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(int16_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<int16_t>(offset);
			
		} else {
			pBC = new TBufCast<int16_t>(offset);
		}
		offset += sizeof(int16_t);

	} else if (HP.IsKeyWord("uint16_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(uint16_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<uint16_t>(offset);
			
		} else {
			pBC = new TBufCast<uint16_t>(offset);
		}
		offset += sizeof(uint16_t);

	} else if (HP.IsKeyWord("int32_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(int32_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<int32_t>(offset);
			
		} else {
			pBC = new TBufCast<int32_t>(offset);
		}
		offset += sizeof(int32_t);

	} else if (HP.IsKeyWord("uint32_t")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(uint32_t).name());
		if (i->second) {
			pBC = new TBufCastHToN<uint32_t>(offset);
			
		} else {
			pBC = new TBufCast<uint32_t>(offset);
		}
		offset += sizeof(uint32_t);

	} else if (HP.IsKeyWord("float")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(float).name());
		if (i->second) {
			pBC = new TBufCastHToN<float>(offset);
			
		} else {
			pBC = new TBufCast<float>(offset);
		}
		offset += sizeof(float);

	} else if (HP.IsKeyWord("double")) {
		TypeMap_t::const_iterator i = swapmap.find(typeid(double).name());
		if (i->second) {
			pBC = new TBufCastHToN<double>(offset);
			
		} else {
			pBC = new TBufCast<double>(offset);
		}
		offset += sizeof(double);

	} else if (HP.IsKeyWord("skip") && !bNoSkip) {
		integer skip = HP.GetInt();
		if (skip < 0) {
			silent_cerr("ReadOneBufCast: invalid number of bytes " << skip
				<< " to be skipped at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		offset += skip;

#if 0
	} else if (HP.IsKeyWord("offset")) {
		integer ioffset = HP.GetInt();
		if (ioffset < 0) {
			silent_cerr("ReadOneBufCast: invalid number of bytes " << ioffset
				<< " as offset at line " << HP.GetLineData() << std::endl);
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		offset = ioffset;
#endif

	} else {
		ASSERT(0);
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}

	return pBC;
}

static bool
bIsLittleEndian(void)
{
	const int i = 1;
	return (*(char *)&i) == 0;
}

static void
SwapMapInit(TypeMap_t& swapmap)
{
	swapmap.insert(TypeMap_t::value_type(typeid(int8_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(uint8_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(int16_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(uint16_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(int32_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(uint32_t).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(float).name(), false));
	swapmap.insert(TypeMap_t::value_type(typeid(double).name(), false));
}

struct TypeWordSet_t : public HighParser::WordSet {
	const TypeMap_t& m_swapmap;
	TypeWordSet_t(const TypeMap_t& swapmap) : m_swapmap(swapmap) { NO_OP; };
	bool IsWord(const std::string& s) const {
		return m_swapmap.find(s) != m_swapmap.end();
	};
};

void
ReadBufCast(HighParser& HP, std::vector<BufCast *>& data)
{
	TypeMap_t swapmap;
	SwapMapInit(swapmap);

	if (HP.IsKeyWord("swap")) {
		TypeWordSet_t typewordset(swapmap);
		const char *s = HP.IsWord(typewordset);
		if (s != 0) {
			do {
				TypeMap_t::iterator i = swapmap.find(std::string(s));
				i->second = true;
			} while ((s = HP.IsWord(typewordset)) != 0);

		} else {
			bool bSwap(false);
			if (HP.IsKeyWord("detect")) {
				bSwap = bIsLittleEndian();

			} else {
				bSwap = HP.GetYesNoOrBool();
			}

			if (bSwap) {
				for (TypeMap_t::iterator i = swapmap.begin(); i != swapmap.end(); ++i) {
					i->second = true;
				}
			}
		}
	}

	if (HP.IsKeyWord("all")) {
		size_t size(0), offset(0);
		BufCast *pBC = 0;
		pBC = ReadOneBufCast(HP, size, swapmap, true);

		data[0] = pBC;
		for (size_t i = 1; i < data.size(); i++) {
			offset += size;
			data[i] = data[i - 1]->copy(offset);
		}

	} else {
		size_t offset(0);
		for (size_t i = 0; i < data.size(); i++) {
retry:;
			BufCast *pBC = 0;
			pBC = ReadOneBufCast(HP, offset, swapmap);

			if (pBC == 0) {
				// got skip
				goto retry;
			}
			data[i] = pBC;
		}
	}
}

