/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <matvecexp.h>

const VecExp ZeroExp(0.);
const MatExp EyeExp(0.);

ostream& VecExp::Write(ostream& out, const char* sFill) const
{
   return out 
     << x.dGet(1) << sFill
     << x.dGet(2) << sFill
     << x.dGet(3) << sFill
     << g.dGet(1) << sFill
     << g.dGet(2) << sFill
     << g.dGet(3);
}


VecExp operator + (const VecExp& v)
{
   return v;
}


VecExp operator - (const VecExp& v)
{
   return VecExp(-v.GetX(), -v.GetG());
}


ostream& operator << (ostream& out, const VecExp& v)
{
   const Vec3& x = v.GetX();
   const Vec3& g = v.GetG();
   return out 
     << x.dGet(1) << " " << x.dGet(2) << " " << x.dGet(3) << " " 
     << g.dGet(1) << " " << g.dGet(2) << " " << g.dGet(3);    
}


ostream& Write(ostream& out, const VecExp& v, const char* sFill)
{
   return v.Write(out, sFill);
}


ostream& MatExp::Write(ostream& out,
		       const char* sFill, 
		       const char* sFill2) const
{
   char* sF2 = (char*)sFill2; 
   if (sFill2 == NULL) {      
      sF2 = (char*)sFill;
   }   
   
   return out 
     << a.dGet(1,1) << sFill
     << a.dGet(1,2) << sFill
     << a.dGet(1,3) << sFill
     << xa.dGet(1,1) << sFill
     << xa.dGet(1,2) << sFill
     << xa.dGet(1,3) << sF2
     << a.dGet(2,1) << sFill
     << a.dGet(2,2) << sFill
     << a.dGet(2,3) << sFill
     << xa.dGet(2,1) << sFill
     << xa.dGet(2,2) << sFill
     << xa.dGet(2,3) << sF2
     << a.dGet(3,1) << sFill
     << a.dGet(3,2) << sFill
     << a.dGet(3,3) << sFill
     << xa.dGet(3,1) << sFill
     << xa.dGet(3,2) << sFill
     << xa.dGet(3,3) << sF2
     << 0. << sFill
     << 0. << sFill
     << 0. << sFill
     << a.dGet(1,1) << sFill
     << a.dGet(1,2) << sFill
     << a.dGet(1,3) << sF2
     << 0. << sFill
     << 0. << sFill
     << 0. << sFill
     << a.dGet(2,1) << sFill
     << a.dGet(2,2) << sFill
     << a.dGet(2,3) << sF2
     << 0. << sFill
     << 0. << sFill
     << 0. << sFill
     << a.dGet(3,1) << sFill
     << a.dGet(3,2) << sFill
     << a.dGet(3,3);
}


ostream& operator << (ostream& out, const MatExp& m)
{
   const Mat3x3& a = m.GetA();
   const Mat3x3& xa = m.GetXA();
   return out 
     << a.dGet(1, 1)  << " " << a.dGet(1, 2)  << " " << a.dGet(1,3)  << " " 
     << xa.dGet(1, 1) << " " << xa.dGet(1, 2) << " " << xa.dGet(1,3) << endl
     << a.dGet(2, 1)  << " " << a.dGet(2, 2)  << " " << a.dGet(2,3)  << " " 
     << xa.dGet(2, 1) << " " << xa.dGet(2, 2) << " " << xa.dGet(2,3) << endl
     << a.dGet(3, 1)  << " " << a.dGet(3, 2)  << " " << a.dGet(3,3)  << " " 
     << xa.dGet(3, 1) << " " << xa.dGet(3, 2) << " " << xa.dGet(3,3) << endl
     << 0.            << " " << 0.            << " " << 0.           << " " 
     << a.dGet(1, 1)  << " " << a.dGet(1, 2)  << " " << a.dGet(1,3)  << endl
     << 0.            << " " << 0.            << " " << 0.           << " " 
     << a.dGet(2, 1)  << " " << a.dGet(2, 2)  << " " << a.dGet(2,3)  << endl
     << 0.            << " " << 0.            << " " << 0.           << " " 
     << a.dGet(3, 1)  << " " << a.dGet(3, 2)  << " " << a.dGet(3,3);     
}


ostream& Write(ostream& out,
	       const MatExp& m,
	       const char* sFill,
	       const char* sFill2)
{
   return m.Write(out, sFill, sFill2);
}


VecExp MultRV(const VecExp& v, const Mat3x3& R)
{
   return VecExp(R*v.GetX(),R*v.GetG());
}


MatExp MultRM(const MatExp& m, const Mat3x3& R)
{
   return MatExp(R*m.GetA(), R*m.GetXA());
}


MatExp MultMRt(const MatExp& m, const Mat3x3& R)
{
   Mat3x3 Rt = R.Transpose();
   return MatExp(m.GetA()*Rt, m.GetXA()*Rt);
}


MatExp MultRMRt(const MatExp& m, const Mat3x3& R)
{   
   Mat3x3 Rt = R.Transpose();
   return MatExp(R*m.GetA()*Rt, R*m.GetXA()*Rt);
}

