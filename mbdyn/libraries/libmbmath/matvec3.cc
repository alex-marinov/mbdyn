/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* vettori 3 e matrici 3x3 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/math.h>
#include <ac/float.h>

#include <matvec3.h>

/* noteworthy constant */
const Mat3x3 Eye3(1., 0., 0., 0., 1., 0., 0., 0., 1.);
const Mat3x3 Zero3x3(0., 0., 0., 0., 0., 0., 0., 0., 0.);
const Vec3 Zero3(0., 0., 0.);


/* Vec3 - begin */

Mat3x3 Vec3::Tens(const Vec3& v) const
{
   return Mat3x3(pdVec[V1]*v.pdVec[V1],
		 pdVec[V2]*v.pdVec[V1],
		 pdVec[V3]*v.pdVec[V1],
		 pdVec[V1]*v.pdVec[V2],
		 pdVec[V2]*v.pdVec[V2],
		 pdVec[V3]*v.pdVec[V2],
		 pdVec[V1]*v.pdVec[V3],
		 pdVec[V2]*v.pdVec[V3],
		 pdVec[V3]*v.pdVec[V3]);
}

/* Prodotto "tensore".  Restituisce se stesso per se stesso */
Mat3x3 Vec3::Tens(void) const
{
	return Tens(*this);
}     

/* Prodotto vettore per matrice */
Mat3x3 Vec3::Cross(const Mat3x3& m) const {
   return Mat3x3(pdVec[V2]*m.pdMat[M31]-pdVec[V3]*m.pdMat[M21],
		 pdVec[V3]*m.pdMat[M11]-pdVec[V1]*m.pdMat[M31],
		 pdVec[V1]*m.pdMat[M21]-pdVec[V2]*m.pdMat[M11],
		 pdVec[V2]*m.pdMat[M32]-pdVec[V3]*m.pdMat[M22],
		 pdVec[V3]*m.pdMat[M12]-pdVec[V1]*m.pdMat[M32],
		 pdVec[V1]*m.pdMat[M22]-pdVec[V2]*m.pdMat[M12],
		 pdVec[V2]*m.pdMat[M33]-pdVec[V3]*m.pdMat[M23],
		 pdVec[V3]*m.pdMat[M13]-pdVec[V1]*m.pdMat[M33],
		 pdVec[V1]*m.pdMat[M23]-pdVec[V2]*m.pdMat[M13]);
}

/* Vec3 - end */

/* Mat3x3 - begin */

/* Prodotto matrice per matrice */
Mat3x3 Mat3x3::operator * (const Mat3x3& m) const
{
   return Mat3x3(pdMat[M11]*m.pdMat[M11]+pdMat[M12]*m.pdMat[M21]+pdMat[M13]*m.pdMat[M31],
		 pdMat[M21]*m.pdMat[M11]+pdMat[M22]*m.pdMat[M21]+pdMat[M23]*m.pdMat[M31],
		 pdMat[M31]*m.pdMat[M11]+pdMat[M32]*m.pdMat[M21]+pdMat[M33]*m.pdMat[M31],
		 pdMat[M11]*m.pdMat[M12]+pdMat[M12]*m.pdMat[M22]+pdMat[M13]*m.pdMat[M32],
		 pdMat[M21]*m.pdMat[M12]+pdMat[M22]*m.pdMat[M22]+pdMat[M23]*m.pdMat[M32],
		 pdMat[M31]*m.pdMat[M12]+pdMat[M32]*m.pdMat[M22]+pdMat[M33]*m.pdMat[M32],
		 pdMat[M11]*m.pdMat[M13]+pdMat[M12]*m.pdMat[M23]+pdMat[M13]*m.pdMat[M33],
		 pdMat[M21]*m.pdMat[M13]+pdMat[M22]*m.pdMat[M23]+pdMat[M23]*m.pdMat[M33],
		 pdMat[M31]*m.pdMat[M13]+pdMat[M32]*m.pdMat[M23]+pdMat[M33]*m.pdMat[M33]);
}


/* inversione */
doublereal Mat3x3::dDet(void) const
{
   doublereal* p = (doublereal*)pdMat;

   return p[M11]*(p[M22]*p[M33]-p[M23]*p[M32])
     +p[M12]*(p[M23]*p[M31]-p[M21]*p[M33])
     +p[M13]*(p[M21]*p[M32]-p[M22]*p[M31]);
}
   
/* inversione */
Mat3x3 Mat3x3::Inv(const doublereal &d) const
{
   ASSERT(fabs(d) > DBL_EPSILON);

   doublereal* p = (doublereal*)pdMat;

   return Mat3x3((p[M22]*p[M33]-p[M23]*p[M32])/d,
		 (p[M23]*p[M31]-p[M21]*p[M33])/d,
		 (p[M21]*p[M32]-p[M22]*p[M31])/d,
		 (p[M13]*p[M32]-p[M12]*p[M33])/d,
		 (p[M11]*p[M33]-p[M13]*p[M31])/d,
		 (p[M12]*p[M31]-p[M11]*p[M32])/d,
		 (p[M12]*p[M23]-p[M13]*p[M22])/d,
		 (p[M13]*p[M21]-p[M11]*p[M23])/d,
		 (p[M11]*p[M22]-p[M12]*p[M21])/d);
}

/* inversione */
Mat3x3 Mat3x3::Inv(void) const
{
   doublereal d = dDet();
   if (fabs(d) < DBL_EPSILON) {
      silent_cerr("matrix is singular" << std::endl);
      throw MatrixHandler::ErrMatrixIsSingular();
   }
   
   return Inv(d);
}


/* soluzione */
Vec3 Mat3x3::Inv(const doublereal& d, const Vec3& v) const
{
   doublereal* p = (doublereal*)pdMat;
   doublereal* pv = v.pGetVec();

   ASSERT(fabs(d) > DBL_EPSILON);
   
   return Vec3((pv[V1]*(p[M22]*p[M33]-p[M23]*p[M32])
		+pv[V2]*(p[M13]*p[M32]-p[M12]*p[M33])
		+pv[V3]*(p[M12]*p[M23]-p[M13]*p[M22]))/d,
	       (pv[V1]*(p[M23]*p[M31]-p[M21]*p[M33])
		+pv[V2]*(p[M11]*p[M33]-p[M13]*p[M31])
		+pv[V3]*(p[M13]*p[M21]-p[M11]*p[M23]))/d,
	       (pv[V1]*(p[M21]*p[M32]-p[M22]*p[M31])
		+pv[V2]*(p[M12]*p[M31]-p[M11]*p[M32])
		+pv[V3]*(p[M11]*p[M22]-p[M12]*p[M21]))/d);
}

/* soluzione */
Vec3 Mat3x3::Inv(const Vec3& v) const
{
   doublereal d = dDet();
   
   if (fabs(d) < DBL_EPSILON) {
      silent_cerr("matrix is singular" << std::endl);
      throw ErrGeneric();
   }

   return Inv(d, v);
}

/* Mat3x3 - end */

/* Manipolatori */
_MatR_Manip MatR;
_MatG_Manip MatG;
_MatGm1_Manip MatGm1;


Vec3 operator - (const Vec3& v)
{
   return Vec3(-v.pdVec[V1],
	       -v.pdVec[V2],
	       -v.pdVec[V3]);
}


Mat3x3 operator - (const Mat3x3& m)
{
   doublereal* pdMat = m.pGetMat();
   return Mat3x3(-pdMat[M11],
		 -pdMat[M21],
		 -pdMat[M31],
		 -pdMat[M12],
		 -pdMat[M22],
		 -pdMat[M32],
		 -pdMat[M13],
		 -pdMat[M23],
		 -pdMat[M33]);
}


const char sForm[] = "%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e%15.6e";
const char sDefFill[] = " ";

/* output di matrici */

std::ostream& 
operator << (std::ostream& out, const Mat3x3& m)
{
   doublereal* pd = m.pGetMat();
#ifdef HAVE_FORM_IN_OSTREAM
   out.form(sForm,
		   pd[M11], pd[M12], pd[M13],
		   pd[M21], pd[M22], pd[M23],
		   pd[M31], pd[M32], pd[M33]);
#else /* !HAVE_FORM_IN_OSTREAM */
   out << pd[M11] << sDefFill << pd[M12] << sDefFill << pd[M13] << sDefFill
     << pd[M21] << sDefFill << pd[M22] << sDefFill << pd[M23] << sDefFill
     << pd[M31] << sDefFill << pd[M32] << sDefFill << pd[M33];
#endif /* !HAVE_FORM_IN_OSTREAM */
   
   return out;
}


std::ostream& 
Write(std::ostream& out, const Mat3x3& m, const char* s, const char* s2)
{
   return m.Write(out, s, s2);
}


/* output di vettori */

std::ostream&
operator << (std::ostream& out, const Vec3& v)
{
   doublereal* pd = v.pGetVec();
   
// out.form(sForm, *pd, *(pd+1), *(pd+2));

   out << pd[0] << sDefFill << pd[1] << sDefFill << pd[2];
   
   return out;
}


std::ostream& 
Write(std::ostream& out, const Vec3& v, const char* s)
{
   return v.Write(out, s);
}


/* Output di matrici */
std::ostream&
Mat3x3::Write(std::ostream& out, const char* sFill, const char* sFill2) const
{
   char* sF2 = (char*)sFill2; 
   if (sFill2 == NULL) {
      sF2 = (char*)sFill;
   }
   out 
     << pdMat[M11] << sFill << pdMat[M12] << sFill << pdMat[M13] << sF2
     << pdMat[M21] << sFill << pdMat[M22] << sFill << pdMat[M23] << sF2
     << pdMat[M31] << sFill << pdMat[M32] << sFill << pdMat[M33];
   
   return out;
}


std::ostream&
Vec3::Write(std::ostream& out, const char* sFill) const
{
   out << pdVec[V1] << sFill << pdVec[V2] << sFill << pdVec[V3];
   
   return out;
}


std::ostream&
Write(std::ostream& out, const doublereal& d, const char*)
{
   return out << d;
}


/* calcolo dei parametri di rotazione a partire dalla matrice R */

Vec3 gparam(const Mat3x3& m)
{
   /* test di singolarita' */
   doublereal d = 1.+m.Trace();
   
   if (d == 0.) {
      silent_cerr("gparam(): divide by zero,"
        " probably due to singularity in rotation parameters" << std::endl);
      throw ErrGeneric();
   }
   
   return m.Ax()*(4./d);
}


/* Calcolo della matrice R a partire da due vettori sghembi */

Mat3x3 MatR2vec(unsigned short int ia, const Vec3& va, 
		unsigned short int ib, const Vec3& vb)
{
   const char sFuncName[] = "MatR2vec()";
   
   ASSERT(ia >= 1 && ia <= 3); 
   ASSERT(va.Dot() > DBL_EPSILON);
   ASSERT(ib >= 1 && ib <= 3); 
   ASSERT(vb.Dot() > DBL_EPSILON);
   
   Vec3 r[3];

   DEBUGCOUT(sFuncName << ": ia = " << ia << " (" << va << "),"
	     << " ib = " << ib << " (" << vb << ")" << std::endl);
   
   if (ia < 1 || ia > 3) {
      silent_cerr(sFuncName << ": first index is illegal" 
	      << std::endl);
      throw ErrGeneric();
   }
   
   int i1 = ia-1;
   int i2 = ia%3;
   int i3 = (ia+1)%3;
   
   if (ib == (ia%3)+1) {
      r[i1] = va/va.Norm();
      r[i3] = r[i1].Cross(vb);
      doublereal d = r[i3].Dot();
      if (d <= DBL_EPSILON) {
	 silent_cerr(sFuncName << ": vectors must be distinct" 
		 << std::endl);
	 throw ErrGeneric();
      }	
      d = sqrt(d);
      r[i3] /= d;
      r[i2] = r[i3].Cross(r[i1]);
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << std::endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else if (ib == ((ia+1)%3+1)) {
      r[i1] = va/va.Norm();	
      r[i2] = vb.Cross(r[i1]);
      doublereal d = r[i2].Dot();
      if (d <= DBL_EPSILON) {
	 silent_cerr(sFuncName << ": vectors must be distinct" 
		 << std::endl);
	 throw ErrGeneric();
      }	
      d = sqrt(d);
      r[i2] /= d;
      r[i3] = r[i1].Cross(r[i2]);  
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << std::endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else {
      silent_cerr(sFuncName << ": second index is illegal" << std::endl);
      throw ErrGeneric();
   }
   
   return Zero3x3; // phony call, not reachable
}


/* Subroutine che calcola gli angoli di Eulero a partire dalla matrice R
      SUBROUTINE RPYPAR(ALFA, PHI)
      
      IMPLICIT NONE
      REAL*8 PHI(3), ALFA(3,3)

      REAL*8 RADEGR
      PARAMETER (RADEGR = 180.D0/3.141592653589793D0)

      REAL*8 CA, SA

      PHI(1) = ATAN2(- ALFA(2,3), ALFA(3,3))
      CA = COS(PHI(1))
      SA = SIN(PHI(1))
      PHI(1) = PHI(1)*RADEGR
      PHI(2) = ATAN2(ALFA(1,3), - SA*ALFA(2,3) + CA*ALFA(3,3))*RADEGR
      PHI(3) = ATAN2(CA*ALFA(2,1) + SA*ALFA(3,1),
     &               CA*ALFA(2,2) + SA*ALFA(3,2))*RADEGR
      RETURN
      END

 */

const doublereal dRaDegr = 180./M_PI;

Vec3 MatR2EulerAngles(const Mat3x3& R)
{  
   doublereal dAlpha = atan2(-R.dGet(2,3), R.dGet(3,3));
   doublereal dCosAlpha = cos(dAlpha);
   doublereal dSinAlpha = sin(dAlpha);
      
   return Vec3(dAlpha,
	       atan2(R.dGet(1,3), 
		     dCosAlpha*R.dGet(3,3)-dSinAlpha*R.dGet(2,3)),
	       atan2(dCosAlpha*R.dGet(2,1)+dSinAlpha*R.dGet(3,1),
		     dCosAlpha*R.dGet(2,2)+dSinAlpha*R.dGet(3,2)));   
}

void MatR2EulerParams(const Mat3x3& R, doublereal& e0, Vec3& e)
{
   doublereal t = R.Tr();
   doublereal T[4];
   T[0] = 1. + t;
   T[1] = 1. - t + 2.*R.dGet(1, 1);
   T[2] = 1. - t + 2.*R.dGet(2, 2);
   T[3] = 1. - t + 2.*R.dGet(3, 3);
   
   int k = 0;
   for (int i = 1; i <= 3; i++) {
      if (fabs(T[i]) >= fabs(T[k])) {
	 k = i;       
      }
   }
   
   switch (k) {
    case 0:   
      e0 = .5*sqrt(T[0]);
      e = Vec3(R.dGet(3, 2)-R.dGet(2, 3),
	       R.dGet(1, 3)-R.dGet(3, 1),
	       R.dGet(2, 1)-R.dGet(1, 2))/(4.*e0);
      break;
      
    case 1: {
       doublereal e1 = .5*sqrt(T[1]);       
       e0 = (R.dGet(3, 2) - R.dGet(2, 3))/(4.*e1);
       e = Vec3(T[1],
		R.dGet(2, 1)+R.dGet(1, 2),
		R.dGet(3, 1)+R.dGet(1, 3))/(4.*e1);
       break;
    }
      
    case 2: {
       doublereal e2 = .5*sqrt(T[2]);
       e0 = (R.dGet(1, 3) - R.dGet(3, 1))/(4.*e2);
       e = Vec3(R.dGet(1, 2)+R.dGet(2, 1),
		T[2],
		R.dGet(3, 2)+R.dGet(2, 3))/(4.*e2);
       break;
    }
      
    case 3: {
       doublereal e3 = .5*sqrt(T[3]);
       e0 = (R.dGet(2, 1) - R.dGet(1, 2))/(4.*e3);
       e = Vec3(R.dGet(1, 3)+R.dGet(3, 1),
		R.dGet(2, 3)+R.dGet(3, 2),
		T[3])/(4.*e3);
       break;
    }
   }
}

Mat3x3 EulerAngles2MatR(const Vec3& v)
{
   doublereal d = v.dGet(1);
   doublereal dCosAlpha(cos(d));
   doublereal dSinAlpha(sin(d));
   d = v.dGet(2);
   doublereal dCosBeta(cos(d));
   doublereal dSinBeta(sin(d));
   d = v.dGet(3);
   doublereal dCosGamma(cos(d));
   doublereal dSinGamma(sin(d));
   
   return Mat3x3(dCosBeta*dCosGamma,
		 dCosAlpha*dSinGamma+dSinAlpha*dSinBeta*dCosGamma,
		 dSinAlpha*dSinGamma-dCosAlpha*dSinBeta*dCosGamma,
		 -dCosBeta*dSinGamma,
		 dCosAlpha*dCosGamma-dSinAlpha*dSinBeta*dSinGamma,
		 dSinAlpha*dCosGamma+dCosAlpha*dSinBeta*dSinGamma,
		 dSinBeta,
		 -dSinAlpha*dCosBeta,
		 dCosAlpha*dCosBeta);
};
