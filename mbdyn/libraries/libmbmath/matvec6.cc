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

#include <matvec6.h>

const Vec6 Zero6(0.);
const Mat6x6 Zero6x6(0.);
const Mat6x6 Eye6(1.);


ostream& Vec6::Write(ostream& out, const char* sFill) const
{
   return out 
     << v[0].dGet(1) << sFill
     << v[0].dGet(2) << sFill
     << v[0].dGet(3) << sFill
     << v[1].dGet(1) << sFill
     << v[1].dGet(2) << sFill
     << v[1].dGet(3);
}


ostream& operator << (ostream& out, const Vec6& v)
{
   const Vec3& v1 = v.GetVec1();
   const Vec3& v2 = v.GetVec2();
   return out 
     << v1.dGet(1) << " " << v1.dGet(2) << " " << v1.dGet(3) << " " 
     << v2.dGet(1) << " " << v2.dGet(2) << " " << v2.dGet(3);    
}


ostream& Write(ostream& out, const Vec6& v, const char* sFill)
{
   return v.Write(out, sFill);
}


ostream& Mat6x6::Write(ostream& out,
		       const char* sFill, 
		       const char* sFill2) const
{
   char* sF2 = (char*)sFill2; 
   if (sFill2 == NULL) {      
      sF2 = (char*)sFill;
   }   
   
   return out 
     << m[0][0].dGet(1,1) << sFill
     << m[0][0].dGet(1,2) << sFill
     << m[0][0].dGet(1,3) << sFill
     << m[0][1].dGet(1,1) << sFill
     << m[0][1].dGet(1,2) << sFill
     << m[0][1].dGet(1,3) << sF2
     << m[0][0].dGet(2,1) << sFill
     << m[0][0].dGet(2,2) << sFill
     << m[0][0].dGet(2,3) << sFill
     << m[0][1].dGet(2,1) << sFill
     << m[0][1].dGet(2,2) << sFill
     << m[0][1].dGet(2,3) << sF2
     << m[0][0].dGet(3,1) << sFill
     << m[0][0].dGet(3,2) << sFill
     << m[0][0].dGet(3,3) << sFill
     << m[0][1].dGet(3,1) << sFill
     << m[0][1].dGet(3,2) << sFill
     << m[0][1].dGet(3,3) << sF2
     << m[1][0].dGet(1,1) << sFill
     << m[1][0].dGet(1,2) << sFill
     << m[1][0].dGet(1,3) << sFill
     << m[1][1].dGet(1,1) << sFill
     << m[1][1].dGet(1,2) << sFill
     << m[1][1].dGet(1,3) << sF2
     << m[1][0].dGet(2,1) << sFill
     << m[1][0].dGet(2,2) << sFill
     << m[1][0].dGet(2,3) << sFill
     << m[1][1].dGet(2,1) << sFill
     << m[1][1].dGet(2,2) << sFill
     << m[1][1].dGet(2,3) << sF2
     << m[1][0].dGet(3,1) << sFill
     << m[1][0].dGet(3,2) << sFill
     << m[1][0].dGet(3,3) << sFill
     << m[1][1].dGet(3,1) << sFill
     << m[1][1].dGet(3,2) << sFill
     << m[1][1].dGet(3,3);
}


ostream& operator << (ostream& out, const Mat6x6& m)
{
   const Mat3x3& m11 = m.GetMat11();
   const Mat3x3& m12 = m.GetMat12();
   const Mat3x3& m21 = m.GetMat21();
   const Mat3x3& m22 = m.GetMat22();
   return out 
     << m11.dGet(1, 1) << " " << m11.dGet(1, 2) << " " << m11.dGet(1,3) << " " 
     << m12.dGet(1, 1) << " " << m12.dGet(1, 2) << " " << m12.dGet(1,3) << endl
     << m11.dGet(2, 1) << " " << m11.dGet(2, 2) << " " << m11.dGet(2,3) << " " 
     << m12.dGet(2, 1) << " " << m12.dGet(2, 2) << " " << m12.dGet(2,3) << endl
     << m11.dGet(3, 1) << " " << m11.dGet(3, 2) << " " << m11.dGet(3,3) << " " 
     << m12.dGet(3, 1) << " " << m12.dGet(3, 2) << " " << m12.dGet(3,3) << endl
     << m21.dGet(1, 1) << " " << m21.dGet(1, 2) << " " << m21.dGet(1,3) << " " 
     << m22.dGet(1, 1) << " " << m22.dGet(1, 2) << " " << m22.dGet(1,3) << endl
     << m21.dGet(2, 1) << " " << m21.dGet(2, 2) << " " << m21.dGet(2,3) << " " 
     << m22.dGet(2, 1) << " " << m22.dGet(2, 2) << " " << m22.dGet(2,3) << endl
     << m21.dGet(3, 1) << " " << m21.dGet(3, 2) << " " << m21.dGet(3,3) << " " 
     << m22.dGet(3, 1) << " " << m22.dGet(3, 2) << " " << m22.dGet(3,3);     
}


ostream& Write(ostream& out,
	       const Mat6x6& m,
	       const char* sFill,
	       const char* sFill2)
{
   return m.Write(out, sFill, sFill2);
}


Vec6 MultRV(const Vec6& v, const Mat3x3& R)
{
   return Vec6(R*v.GetVec1(),R*v.GetVec2());
}


Mat6x6 MultRM(const Mat6x6& m, const Mat3x3& R)
{
   return Mat6x6(R*m.GetMat11(), R*m.GetMat21(),
		 R*m.GetMat12(), R*m.GetMat22());
}


Mat6x6 MultMRt(const Mat6x6& m, const Mat3x3& R)
{
   Mat3x3 Rt = R.Transpose();
   return Mat6x6(m.GetMat11()*Rt, m.GetMat21()*Rt, 
		 m.GetMat12()*Rt, m.GetMat22()*Rt);   
}


Mat6x6 MultRMRt(const Mat6x6& m, const Mat3x3& R)
{   
   Mat3x3 Rt = R.Transpose();
   return Mat6x6(R*m.GetMat11()*Rt, R*m.GetMat21()*Rt, 
		 R*m.GetMat12()*Rt, R*m.GetMat22()*Rt);   
}


/* esegue l'operazione:
 * [I   0] [     ]
 * [     ] [  m  ]
 * [vx  I] [     ] */
Mat6x6 MultVM(const Mat6x6& m, const Vec3& v)
{
   return Mat6x6(m.GetMat11(), 
		 v.Cross(m.GetMat11())+m.GetMat21(),
		 m.GetMat12(),
		 v.Cross(m.GetMat12())+m.GetMat22());
}


/* esegue l'operazione:
 * [I  vx] [     ]
 * [     ] [  m  ]
 * [0   I] [     ] */
Mat6x6 MultVtM(const Mat6x6& m, const Vec3& v)
{
   return Mat6x6(m.GetMat11()+v.Cross(m.GetMat21()),
		 m.GetMat21(),
		 m.GetMat12()+v.Cross(m.GetMat22()),
		 m.GetMat22());
}


/* esegue l'operazione:
 * [     ] [I  vx] 
 * [  m  ] [     ] 
 * [     ] [0   I] */
Mat6x6 MultMV(const Mat6x6& m, const Vec3& v)
{
   return Mat6x6(m.GetMat11(),
		 m.GetMat21(),
		 m.GetMat11()*Mat3x3(v)+m.GetMat12(),
		 m.GetMat21()*Mat3x3(v)+m.GetMat22());
}


/* esegue l'operazione:
 * [     ] [I   0] 
 * [  m  ] [     ] 
 * [     ] [vx  I] */
Mat6x6 MultMVt(const Mat6x6& m, const Vec3& v)
{
   return Mat6x6(m.GetMat11()+m.GetMat12()*Mat3x3(v),
		 m.GetMat21()+m.GetMat22()*Mat3x3(v),
		 m.GetMat12(),
		 m.GetMat22());
}
