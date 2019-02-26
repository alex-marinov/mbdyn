/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/

//                                 MATHEM.CC

// Contiene le funzioni relative alle matrici e alle operazioni
// nel piano

#include <mathem.h>

Mat6x6 KMatrix  (double L, double Ixx, double Iyy, double Izz,
		 double A, double E, double G, double Asy, double Asz,
		 double* Cmatrix, double Cratio)
{
   double Py,Pz;
   Mat6x6* K=new Mat6x6;
   Py=12*E*Izz*Asy/(G*A*L*L);
   Pz=12*E*Iyy*Asz/(G*A*L*L);
#if 0
   (*K)[0][0]=E*A/L;
   (*K)[1][1]=12*E*Izz/(L*L*L*(1+Py));
   (*K)[1][5]=-6*E*Izz/(L*L*(1+Py));
   (*K)[2][2]=12*E*Iyy/(L*L*L*(1+Pz));
   (*K)[2][4]=6*E*Iyy/(L*L*(1+Pz));
   (*K)[3][3]=G*Ixx/L;
   (*K)[4][4]=(4+Pz)*E*Iyy/(L*(1+Pz));
   (*K)[5][5]=(4+Py)*E*Izz/(L*(1+Py));
#endif
   (*K)[0][0]=E*A;
   (*K)[1][1]=G*A;
   (*K)[2][2]=G*A;
   (*K)[3][3]=G*Ixx;
   (*K)[4][4]=E*Iyy;
   (*K)[5][5]=E*Izz;
   return (*K);
}

void Gauss (Vector& X, Matrix& AA, Vector& R)
{
   if (Det(AA)==0) {
      cout << "Il sistema non ammette soluzione univoca." << endl;
      return;
   }
   double sum=0;
   unsigned int ord=AA.rows();
   Vector Ra(ord);
   for (int k=0;k<ord;k++) Ra[k]=R[k];
   Matrix A(ord,ord+1);
   Vector XORD(ord);
   Vector temp(ord);
   Vector temp2(ord);
   Vector BUF(ord+1);
   for (int i=0; i<ord; i++)
     for (int j=0; j<ord; j++) A[i][j]=AA[i][j];
   for (int i=0; i<ord; i++) A[i][ord]=Ra[i];
   Matrix cprimo(ord,ord+1);
   cprimo=A;
   Vector rik(ord);
   for (int u=1;u<=ord;u++) rik[u-1]=(double) u;
   for (int I=1;I<=ord;I++) {
      if (A[I-1][I-1]==0) {
	 unsigned int lp=I;
	 unsigned int flag=0;
	 while ((flag==0)&(lp<(ord))) {
	    lp=lp+1;
	    if (A[I-1][lp-1]!=0) flag=1;
	 }
	 if (flag==1) {
	    for (int i=1;i<=ord;i++) {
	       temp[i-1]=A[i-1][I-1];
	       temp2[i-1]=cprimo[i-1][I-1];
	       A[i-1][I-1]=A[i-1][lp-1];
	       A[i-1][lp-1]=temp[i-1];
	       cprimo[i-1][I-1]=cprimo[i-1][lp-1];
	       cprimo[i-1][lp-1]=temp2[i-1];
	    }
	    double t=rik[I-1];
	    rik[I-1]=rik[lp-1];
	    rik[lp-1]=t;
	 }
      }
      for (int K=I;K<=(ord+1);K++) cprimo[I-1][K-1]=(A[I-1][K-1])/(A[I-1][I-1]);
      for (int J=I+1;J<=ord;J++)
	for (int K=1;K<=(ord+1);K++)
	   cprimo[J-1][K-1]=A[J-1][K-1]-cprimo[I-1][K-1]*A[J-1][I-1];
      A=cprimo;
   }
   X[ord-1]=A[ord-1][ord];
   for (int I=(ord-1);I>0;I--) {
      sum=0;
      for (int K=I+1;K<=ord;K++) sum=sum+A[I-1][K-1]*X[K-1];
      X[I-1]=A[I-1][ord]-sum;
   }
   for (int i=1;i<=ord;i++) XORD[rik[i-1]-1]=X[i-1];
   for (int i=0;i<ord;i++) X[i]=XORD[i];
}

// Matrice inversa

Matrix Inv(Matrix& R) {
   unsigned int r=R.rows();
   unsigned int c=R.columns();
   Vector BUF(r);// Vettore che contiene la colonna i-esima
   Vector XBUF(r); // Vettore che contiene la soluzione per la colonna i-esima
   Matrix TM(R); // Matrice che conterra' l'inversa
   Matrix RI(r);
   /* Rende RI identità */
   for (int k; k<r; k++)
     for (int j; j<r; j++)
       if (k==j) RI[k][j]=1; else RI[k][j]=0;
   if (r != c) cout << "Inversione non possibile!" << endl;
   else if (Det(R)==0) cout << "Matrice a determinante nullo!" << endl;
   else {
      // Procedura di inversione
      for (int I=0;I<c;I++) {
	 BUF.Clear(); // Svuota il Vettore buffer
	 for (int i=0;i<r;i++) BUF[i]=RI[i][I]; // Copia la colonna I-esima
	 Gauss(XBUF,R,BUF);
	 for (int i=0;i<r;i++) TM[i][I]=XBUF[i];
      }
   }
   return TM;
}

Matrix operator / (Matrix& A, Matrix& B)
{
   unsigned int r=A.rows();
   unsigned int c=A.columns();
   Matrix P (r,c);
   if (Det(B)==0) cout << "Divisione impossibile" << endl; exit(0);
   P=A*Inv(B);
   return P;
}

double square (double p)
{
   return (p*p);
}


Mat3x3 MomentToInertia(Vec3 Diag, Vec3 Sym)
{
   /* La funzione a partire dai valori di Ixx,Iyy,Izz,.. restituisce
    * la matrice contenente le proprietà inerziali */
   Mat3x3 I;
   I[0][0]=Diag[0];
   I[1][1]=Diag[1];
   I[2][2]=Diag[2];
   I[0][1]=Sym[0];
   I[0][2]=Sym[1];
   I[1][0]=Sym[0];
   I[1][2]=Sym[2];
   I[2][0]=Sym[1];
   I[2][1]=Sym[2];
   return I;
}

double a360tan2(double y, double x)
{
   double t;
   t=atan2(y,x);
   if (t<0) t=((2*pi)-fabs(t));
   return t;
}

Vec3 gparam (const Mat3x3& R)
{ 
   doublereal d = 1+R[0][0]+R[1][1]+R[2][2];
   if (d == 0) {
      cout << "gparam () : divide by 0, singularity in rotation parameters"
	<< endl;
      exit (0);
   }
   d = 2./d;
   Vec3 Answer (d*(R[1][2]-R[2][1]),d*(R[2][0]-R[0][2]),
		d*(R[0][1]-R[1][0]));
   return Answer;
}


Mat3x3 MatR2vec (unsigned short int ia, const Vec3& va,
		 unsigned short int ib, const Vec3& vb)
{
   Vec3 r[3];

   const char* sFuncName = "MatR2vec";
   
   int i1 = ia-1;
   int i2 = ia%3;
   int i3 = (ia+1)%3;
   
   if (ib == (ia%3)+1) {
      r[i1] = va/va.Norm();
      r[i3] = r[i1].Cross(vb);
      double d = r[i3].Dot();

      /* CONTROLLA CHE LA NORMA DEL VETTORE PROIEZIONE SIA != 0 */
      if (d <= DBL_EPSILON) {
	 cerr << endl << sFuncName << ": vectors must be distinct" << endl;
	 exit(-3);
      }	

      d = sqrt(d);
      r[i3] /= d;
      r[i2] = r[i3].Cross(r[i1]);
      return Mat3x3(r[0], r[1], r[2]);
   } 
   else if (ib == ((ia+1)%3+1)) {
      r[i1] = va/va.Norm();	
      r[i2] = vb.Cross(r[i1]);
      double d = r[i2].Dot();

      /* CONTROLLA CHE LA NORMA DEL VETTORE PROIEZIONE SIA != 0 */
      if (d <= DBL_EPSILON) {
	 cerr << endl << sFuncName << ": vectors must be distinct" << endl;
	 exit(-3);
      }	

      d = sqrt(d);
      r[i2] /= d;
      r[i3] = r[i1].Cross(r[i2]);  
      return Mat3x3(r[0], r[1], r[2]);
   } else {
      cerr << endl << sFuncName << ": second index is illegal" << endl;
      exit(-3);
   }
   
   Mat3x3 A;
   return A; // phony call, not reachable
}

Vec3 EulerAngles (const Mat3x3& R)
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

Mat3x3 RFromEulerAngles (const Vec3& v)
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
}

