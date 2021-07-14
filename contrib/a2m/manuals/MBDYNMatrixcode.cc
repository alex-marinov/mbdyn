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
	     << " ib = " << ib << " (" << vb << ")" << endl);
   
   if (ia < 1 || ia > 3) {
      cerr << endl << sFuncName << ": first index is illegal" << endl;
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
   }
   
   int i1 = ia-1;
   int i2 = ia%3;
   int i3 = (ia+1)%3;
   
   if (ib == (ia%3)+1) {
      r[i1] = va/va.Norm();
      r[i3] = r[i1].Cross(vb);
      double d = r[i3].Dot();
      if (d <= DBL_EPSILON) {
	 cerr << endl << sFuncName << ": vectors must be distinct" << endl;
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }	
      d = sqrt(d);
      r[i3] /= d;
      r[i2] = r[i3].Cross(r[i1]);
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else if (ib == ((ia+1)%3+1)) {
      r[i1] = va/va.Norm();	
      r[i2] = vb.Cross(r[i1]);
      double d = r[i2].Dot();
      if (d <= DBL_EPSILON) {
	 cerr << endl << sFuncName << ": vectors must be distinct" << endl;
	 throw ErrGeneric(MBDYN_EXCEPT_ARGS);
      }	
      d = sqrt(d);
      r[i2] /= d;
      r[i3] = r[i1].Cross(r[i2]);  
      
      DEBUGCOUT("R = " << Mat3x3(r[0], r[1], r[2]) << endl);
      
      return Mat3x3(r[0], r[1], r[2]);
   } else {
      cerr << endl << sFuncName << ": second index is illegal" << endl;
      throw ErrGeneric(MBDYN_EXCEPT_ARGS);
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

Vec3 EulerAngles(const Mat3x3& R)
{  
   doublereal dAlpha = atan2(-R.dGet(2,3), R.dGet(3,3));
   doublereal dCosAlpha = cos(dAlpha);
   doublereal dSinAlpha = sin(dAlpha);
      
   return Vec3(dAlpha*dRaDegr,
	       atan2(R.dGet(1,3), 
		     dCosAlpha*R.dGet(3,3)-dSinAlpha*R.dGet(2,3))*dRaDegr,
	       atan2(dCosAlpha*R.dGet(2,1)+dSinAlpha*R.dGet(3,1),
		     dCosAlpha*R.dGet(2,2)+dSinAlpha*R.dGet(3,2))*dRaDegr);   
}

void EulerParams(const Mat3x3& R, doublereal& e0, Vec3& e)
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

Mat3x3 RFromEulerAngles(const Vec3& v)
{
   doublereal d;
   doublereal dCosAlpha(cos((d = v.dGet(1))));
   doublereal dSinAlpha(sin(d));
   doublereal dCosBeta(cos((d = v.dGet(2))));
   doublereal dSinBeta(sin(d));
   doublereal dCosGamma(cos((d = v.dGet(3))));
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
