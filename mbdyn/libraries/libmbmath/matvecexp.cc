/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
 * 
 * This code is a partial merge of HmFe and MBDyn.
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 * Teodoro Merlini  <merlini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */


#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "matvecexp.h"

const VecExp ZeroExp(0.);
const MatExp EyeExp(0.);

ostream&
VecExp::Write(ostream& out, const char* sFill) const
{
	return out 
		<< vec.dGet(1) << sFill
		<< vec.dGet(2) << sFill
		<< vec.dGet(3) << sFill
		<< mom.dGet(1) << sFill
		<< mom.dGet(2) << sFill
		<< mom.dGet(3);
}

VecExp
operator - (const VecExp& v)
{
	return VecExp(-v.Vec(), -v.Mom());
}

ostream&
operator << (ostream& out, const VecExp& v)
{
	const Vec3& vec = v.Vec();
	const Vec3& mom = v.Mom();
	return out 
		<< vec.dGet(1) << " "
		<< vec.dGet(2) << " "
		<< vec.dGet(3) << " " 
		<< mom.dGet(1) << " "
		<< mom.dGet(2) << " "
		<< mom.dGet(3);    
}

ostream&
Write(ostream& out, const VecExp& v, const char* sFill)
{
	return v.Write(out, sFill);
}

ostream&
MatExp::Write(ostream& out, const char* sFill, const char* sFill2) const
{
	char* sF2 = (char*)sFill2; 
	if (sFill2 == NULL) {      
		sF2 = (char*)sFill;
	}   
	
	return out 
		<< vec.dGet(1,1) << sFill
		<< vec.dGet(1,2) << sFill
		<< vec.dGet(1,3) << sFill
		<< mom.dGet(1,1) << sFill
		<< mom.dGet(1,2) << sFill
		<< mom.dGet(1,3) << sF2
		<< vec.dGet(2,1) << sFill
		<< vec.dGet(2,2) << sFill
		<< vec.dGet(2,3) << sFill
		<< mom.dGet(2,1) << sFill
		<< mom.dGet(2,2) << sFill
		<< mom.dGet(2,3) << sF2
		<< vec.dGet(3,1) << sFill
		<< vec.dGet(3,2) << sFill
		<< vec.dGet(3,3) << sFill
		<< mom.dGet(3,1) << sFill
		<< mom.dGet(3,2) << sFill
		<< mom.dGet(3,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(1,1) << sFill
		<< vec.dGet(1,2) << sFill
		<< vec.dGet(1,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(2,1) << sFill
		<< vec.dGet(2,2) << sFill
		<< vec.dGet(2,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(3,1) << sFill
		<< vec.dGet(3,2) << sFill
		<< vec.dGet(3,3);
}

ostream&
operator << (ostream& out, const MatExp& m)
{
	const Mat3x3& vec = m.Vec();
	const Mat3x3& mom = m.Mom();
	return out 
		<< vec.dGet(1, 1) << " "
		<< vec.dGet(1, 2) << " "
		<< vec.dGet(1, 3) << " " 
		<< mom.dGet(1, 1) << " "
		<< mom.dGet(1, 2) << " "
		<< mom.dGet(1, 3) << endl
		<< vec.dGet(2, 1) << " "
		<< vec.dGet(2, 2) << " "
		<< vec.dGet(2, 3) << " " 
		<< mom.dGet(2, 1) << " "
		<< mom.dGet(2, 2) << " "
		<< mom.dGet(2, 3) << endl
		<< vec.dGet(3, 1) << " "
		<< vec.dGet(3, 2) << " "
		<< vec.dGet(3, 3) << " " 
		<< mom.dGet(3, 1) << " "
		<< mom.dGet(3, 2) << " "
		<< mom.dGet(3, 3) << endl
		<< 0. << " "
		<< 0. << " "
		<< 0. << " " 
		<< vec.dGet(1, 1) << " "
		<< vec.dGet(1, 2) << " "
		<< vec.dGet(1, 3) << endl
		<< 0. << " "
		<< 0. << " " 
		<< 0. << " " 
		<< vec.dGet(2, 1) << " "
		<< vec.dGet(2, 2) << " "
		<< vec.dGet(2, 3) << endl
		<< 0. << " "
		<< 0. << " "
		<< 0. << " " 
		<< vec.dGet(3, 1) << " "
		<< vec.dGet(3, 2) << " "
		<< vec.dGet(3, 3);     
}

ostream&
Write(ostream& out, const MatExp& m, const char* sFill, const char* sFill2)
{
	return m.Write(out, sFill, sFill2);
}

#if 0
VecExp
MultRV(const VecExp& v, const Mat3x3& R)
{
	return VecExp(R*v.Vec(),R*v.Mom());
}

MatExp
MultRM(const MatExp& m, const Mat3x3& R)
{
	return MatExp(R*m.GetA(), R*m.GetXA());
}


MatExp
MultMRt(const MatExp& m, const Mat3x3& R)
{
	Mat3x3 Rt = R.Transpose();
	return MatExp(m.GetA()*Rt, m.GetXA()*Rt);
}

MatExp
MultRMRt(const MatExp& m, const Mat3x3& R)
{   
	Mat3x3 Rt = R.Transpose();
	return MatExp(R*m.GetA()*Rt, R*m.GetXA()*Rt);
}
#endif /* 0 */

