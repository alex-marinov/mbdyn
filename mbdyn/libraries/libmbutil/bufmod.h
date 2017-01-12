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

extern void
ReadBufCast(HighParser& HP, std::vector<BufCast *>& data);

/* BufCast - end */

#endif // BUFMOD_H

