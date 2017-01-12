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

#ifndef PRIVPGIN_H
#define PRIVPGIN_H

#include "mathp.h"
#include "dataman.h"

class PrivPlugIn : public MathParser::PlugIn {
	friend std::ostream& operator << (std::ostream& out, const PrivPlugIn& p);

protected:
	SimulationEntity *pSE;
	unsigned int iIndex;
	std::string sIndexName;
	DataManager *pDM;

public:
	PrivPlugIn(MathParser& mp, DataManager *pDM);
	virtual ~PrivPlugIn(void);
	virtual const char *sName(void) const = 0;
	int Read(int argc, char *argv[]);
	TypedValue::Type GetType(void) const;
	TypedValue GetVal(void) const;

protected:
	unsigned int ReadLabel(const char* s);
	virtual void ReadSE(unsigned int uLabel, const char *s) = 0;
	void ReadIndex(unsigned int iMaxIndex, const char *s);
	virtual std::ostream& Err(std::ostream& out) const = 0;
};

class NodePrivPlugIn : public PrivPlugIn {
public:
	NodePrivPlugIn(MathParser& mp, DataManager *pDM);
	virtual ~NodePrivPlugIn(void);
	const char *sName(void) const;

protected:
	virtual void ReadSE(unsigned int uLabel, const char *s);
	virtual std::ostream& Err(std::ostream& out) const;
};

class ElemPrivPlugIn : public PrivPlugIn {
public:
	ElemPrivPlugIn(MathParser& mp, DataManager *pDM);
	virtual ~ElemPrivPlugIn(void);
	const char *sName(void) const;

protected:
	virtual void ReadSE(unsigned int uLabel, const char *s);
	virtual std::ostream& Err(std::ostream& out) const;
};

extern std::ostream& 
operator << (std::ostream& out, const PrivPlugIn& p);

extern MathParser::PlugIn *
node_priv_plugin(MathParser& mp, void *arg);

extern MathParser::PlugIn *
elem_priv_plugin(MathParser& mp, void *arg);

#endif /* PRIVPGIN_H */

