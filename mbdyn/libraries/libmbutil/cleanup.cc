/* $Header$*/
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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
#include "cleanup.h"

#include <stack>

struct cleanup {
	mbdyn_cleanup_f handler;
	void *data;

public:
	cleanup(mbdyn_cleanup_f handler);
	~cleanup(void);
};

cleanup::cleanup(mbdyn_cleanup_f handler)
: handler(handler), data(0)
{
	NO_OP;
}

cleanup::~cleanup(void)
{
	NO_OP;
}

static std::stack<cleanup *> c;

extern "C" int
mbdyn_cleanup_register(mbdyn_cleanup_f handler, void ***datapp)
{

	pedantic_cout("mbdyn_cleanup_register: " << (void *)handler << ":" << (void *)datapp << std::endl);

	cleanup *p = new cleanup(handler);
	if (datapp != 0) {
		*datapp = &p->data;
	}
	c.push(p);

	return 0;
}

extern "C" int
mbdyn_cleanup(void)
{
	while (!c.empty()) {
		cleanup *p = c.top();

		pedantic_cout("mbdyn_cleanup: " << (void *)p->handler << ":" << (void *)p->data << std::endl);

		p->handler(p->data);
		c.pop();
		delete p;
	}

	return 0;
}

extern "C" void
mbdyn_cleanup_destroy(void)
{
	while (!c.empty()) {
		cleanup *p = c.top();
		c.pop();
		delete p;
	}
}

