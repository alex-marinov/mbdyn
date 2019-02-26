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

#include <matvec6.h>

// NOTE: do not use Zero3, Zero3x3 or mb_zero<>()
// because they might be not initialized yet
const Vec6 Zero6(0., 0., 0., 0., 0., 0.);
const Mat6x6 Zero6x6(0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0.,
		0., 0., 0., 0., 0., 0.);
const Mat6x6 Eye6(1., 0., 0., 0., 0., 0.,
		0., 1., 0., 0., 0., 0.,
		0., 0., 1., 0., 0., 0.,
		0., 0., 0., 1., 0., 0.,
		0., 0., 0., 0., 1., 0.,
		0., 0., 0., 0., 0., 1.);


std::ostream&
Vec6::Write(std::ostream& out, const char* sFill) const
{
   return out 
     << v[0].dGet(1) << sFill
     << v[0].dGet(2) << sFill
     << v[0].dGet(3) << sFill
     << v[1].dGet(1) << sFill
     << v[1].dGet(2) << sFill
     << v[1].dGet(3);
}

void
Vec6::Reset(void)
{
	v[0].Reset();
	v[1].Reset();
}

void
Mat6x6::Reset(void)
{
	m[0][0].Reset();
	m[0][1].Reset();
	m[1][0].Reset();
	m[1][1].Reset();
}

Vec6 operator + (const Vec6& v)
{
   return v;
}


Vec6 operator - (const Vec6& v)
{
   return Vec6(-v.GetVec1(), -v.GetVec2());
}


std::ostream& 
operator << (std::ostream& out, const Vec6& v)
{
   const Vec3& v1 = v.GetVec1();
   const Vec3& v2 = v.GetVec2();
   return out 
     << v1.dGet(1) << " " << v1.dGet(2) << " " << v1.dGet(3) << " " 
     << v2.dGet(1) << " " << v2.dGet(2) << " " << v2.dGet(3);    
}


std::ostream&
Write(std::ostream& out, const Vec6& v, const char* sFill)
{
   return v.Write(out, sFill);
}


std::ostream&
Mat6x6::Write(std::ostream& out, const char* sFill, const char* sFill2) const
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


std::ostream&
operator << (std::ostream& out, const Mat6x6& m)
{
   const Mat3x3& m11 = m.GetMat11();
   const Mat3x3& m12 = m.GetMat12();
   const Mat3x3& m21 = m.GetMat21();
   const Mat3x3& m22 = m.GetMat22();
   return out 
     << m11.dGet(1, 1) << " " << m11.dGet(1, 2) << " " << m11.dGet(1,3) << " " 
     << m12.dGet(1, 1) << " " << m12.dGet(1, 2) << " " << m12.dGet(1,3) << std::endl
     << m11.dGet(2, 1) << " " << m11.dGet(2, 2) << " " << m11.dGet(2,3) << " " 
     << m12.dGet(2, 1) << " " << m12.dGet(2, 2) << " " << m12.dGet(2,3) << std::endl
     << m11.dGet(3, 1) << " " << m11.dGet(3, 2) << " " << m11.dGet(3,3) << " " 
     << m12.dGet(3, 1) << " " << m12.dGet(3, 2) << " " << m12.dGet(3,3) << std::endl
     << m21.dGet(1, 1) << " " << m21.dGet(1, 2) << " " << m21.dGet(1,3) << " " 
     << m22.dGet(1, 1) << " " << m22.dGet(1, 2) << " " << m22.dGet(1,3) << std::endl
     << m21.dGet(2, 1) << " " << m21.dGet(2, 2) << " " << m21.dGet(2,3) << " " 
     << m22.dGet(2, 1) << " " << m22.dGet(2, 2) << " " << m22.dGet(2,3) << std::endl
     << m21.dGet(3, 1) << " " << m21.dGet(3, 2) << " " << m21.dGet(3,3) << " " 
     << m22.dGet(3, 1) << " " << m22.dGet(3, 2) << " " << m22.dGet(3,3);     
}


std::ostream&
Write(std::ostream& out, const Mat6x6& m, const char* sFill, 
		const char* sFill2)
{
   return m.Write(out, sFill, sFill2);
}


Vec6 MultRV(const Vec6& v, const Mat3x3& R)
{
   return Vec6(R*v.GetVec1(), R*v.GetVec2());
}


Mat6x6 MultRM(const Mat6x6& m, const Mat3x3& R)
{
   return Mat6x6(R*m.GetMat11(), R*m.GetMat21(),
		 R*m.GetMat12(), R*m.GetMat22());
}


Mat6x6 MultMRt(const Mat6x6& m, const Mat3x3& R)
{
   return Mat6x6(m.GetMat11().MulMT(R), m.GetMat21().MulMT(R), 
		 m.GetMat12().MulMT(R), m.GetMat22().MulMT(R));   
}


Mat6x6 MultRMRt(const Mat6x6& m, const Mat3x3& R)
{   
   return Mat6x6(R*m.GetMat11().MulMT(R), R*m.GetMat21().MulMT(R), 
		 R*m.GetMat12().MulMT(R), R*m.GetMat22().MulMT(R));   
}

Mat6x6 MultRMRt(const Mat6x6& m, const Mat3x3& R, const doublereal& c)
{
   Mat3x3 Rc(R*c);
   return Mat6x6(R*m.GetMat11().MulMT(Rc), R*m.GetMat21().MulMT(Rc), 
		 R*m.GetMat12().MulMT(Rc), R*m.GetMat22().MulMT(Rc));   
}


/* esegue l'operazione:
 * [I   0] [     ]
 * [     ] [  m  ]
 * [vx  I] [     ] */
Mat6x6 MultVCrossMat(const Mat6x6& m, const Vec3& v)
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
Mat6x6 MultVCrossTMat(const Mat6x6& m, const Vec3& v)
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
Mat6x6 MultMatVCross(const Mat6x6& m, const Vec3& v)
{
	Mat3x3 vCross(MatCross, v);

	return Mat6x6(m.GetMat11(),
		m.GetMat21(),
		m.GetMat11()*vCross + m.GetMat12(),
		m.GetMat21()*vCross + m.GetMat22());
}


/* esegue l'operazione:
 * [     ] [I   0] 
 * [  m  ] [     ] 
 * [     ] [vx  I] */
Mat6x6 MultMatVCrossT(const Mat6x6& m, const Vec3& v)
{
	Mat3x3 vCross(MatCross, v);

	return Mat6x6(m.GetMat11() + m.GetMat12()*vCross,
		m.GetMat21() + m.GetMat22()*vCross,
		m.GetMat12(),
		m.GetMat22());
}
