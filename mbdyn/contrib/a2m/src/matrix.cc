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

//                                MATRIX.CC

#include <matrix.h>
#include <debug.h>

extern Boolean EXTENDED_MATRIX_DISPLAY;


/* VECTOR */

Vector::Vector (unsigned int n=8) : sz(n), dato (new double[n])
{
   for (int i=0; i<n; i++) dato[i]=0;
}

Vector::Vector (const Vector& v) : sz(v.sz), dato(new double[v.sz])
{
   copy(v);
}

void Vector::copy (const Vector& v)
{
   unsigned min_size= (sz < v.sz? sz : v.sz);
   for (int i=0; i< min_size;i++)
     dato[i]=v.dato[i];
}

Vector::~Vector () { delete dato; }

double& Vector::operator [] (unsigned int i) const
{
   if (i>(size()-1)) {
      cout << "Fail to address:" << i << ", " << (size()-1) << endl;
      DEBUGCOUT ("Errore di indirizzamento in vettore");
   }
   return dato[(i)];
}

Vector& Vector::operator = (const Vector& v)
{
   sz=v.sz;
   dato=new double[sz];
   for (int j=0; j<sz; j++) dato[j]=v.dato[j];
   return *this;
}

unsigned int Vector::size() const { return sz; }

ostream& operator << (ostream& out, const Vector& r)
{
   for (int j=0;j<r.size();j++) out << r.dato[j] << ", ";
   out << endl;
   return out;
}

ostream& Vector::Write (ostream& out, const char* fill) const
{
   /* Check if vector is null */
   Boolean IsNull=Y;
   for (int i=0; i<sz; i++) if (dato[i] != 0) IsNull=N;
   if ((IsNull) && (EXTENDED_MATRIX_DISPLAY==N)) out << "null";
   else {
   for (int i=0; i<(sz-1); i++) out << dato[i] << fill;
      out << dato[sz-1];
   }
   return out;
}

istream& operator >> (istream& in, Vector& r)
{
   for (int j=0; j<r.size();j++) in >> r.dato[j];
   return in;
}

Vector operator + (const Vector& a, const Vector& b)
{
   if (a.size() != b.size()) return a;
   Vector sum(a.size());
   for (int i=0;i<a.size();i++)
     sum.dato[i]=(a.dato[i]+b.dato[i]);
   return sum;
}

Vector operator - (const Vector& a, const Vector& b)
{
   if (a.size() != b.size()) return a;
   Vector diff(a.size());
   for (int i=0; i<a.size();i++)
     diff.dato[i]=(a.dato[i]-b.dato[i]);
   return diff;
}

double operator * (const Vector& a, const Vector& b)
{
   /*La funzione restituisce il prodotto scalare tra i vettori*/
   if (a.size()!=b.size()) return 0;
   double R;
   R=(Module(a)*Module(b))*Cos_b2Vector(a,b);
   return R;
}

Vector operator * (const Vector& a, const double q)
{
   Vector R(a.size());
   for (int i=0; i<a.size();i++)
     R.dato[i]=a.dato[i]*q;
   return R;
}

Vector operator * (const double q, const Vector& a)
{
   return (a*q);
}

Vector operator / (const Vector& a, const double q)
{
   return a*(1/q);
}

Vector& operator /= (Vector& a, const double q)
{
   for (int i=0;i<a.size();i++)
     a[i] /= q;
   return (a);
}


Vector operator ^ (const Vector& a, const Vector& b)
{
   /*La funzione restituisce il prodotto vettore, ora solo in 3D*/
   double c1,c2,c3;
   c1= (a[1]*b[2])-(a[2]*b[1]);
   c2= -(a[0]*b[2])-(a[2]*b[0]);
   c3= (a[0]*b[1])-(b[0]*a[1]);
   Vector R(3);
   R[0]=c1; R[1]=c2; R[2]=c3;
   return a;
}

void Vector::Clear()
{
   for (int j=0; j<size();j++) dato[j]=0;
}


// Sottovettori di 2 e 3 elementi

void Vec3::Set(double i, double j, double k)
{
   dato[0]=i; dato[1]=j; dato[2]=k;
}

void Vec2::Set(double i, double j)
{
   dato[0]=i; dato[1]=j;
}



/* MATRIX */

Matrix::Matrix (unsigned int n, unsigned int m) : rws(n),clms(m)
{
   data=new Vector* [n];
   for (int i=0; i<rws; i++) data[i]=new Vector(m);
}

Matrix::Matrix (const Matrix& R) : rws(R.rws),clms(R.clms)
{
   data=new Vector*[rws];
   for (int i=0;i<rws;i++) data[i]=new Vector(clms);
   for (int i=0;i<rws;i++) (*data[i])=(*R.data[i]);
}

Matrix::~Matrix ()
{
   // for (int i=0;i<rws;i++) delete data[i];
}

Vector& Matrix::operator [] (unsigned int j) const 
{
   if (j>(rows()-1)) DEBUGCOUT ("Errore di indirizzamento in matrice");
   return *data[j];
}

Matrix Matrix::tr() 
{
   Matrix MIRROR (columns(),rows());
   for (int i=0; i<columns(); i++)
     for (int j=0; j<rows();j++)
       MIRROR[i][j]=(*data[j])[i];
   for (int i=0; i<rows(); i++)
     for (int j=0; j<columns();j++)
       (*data[i])[j]=MIRROR[j][i];
   return (MIRROR);
}

ostream& operator << (ostream& out, const Matrix& r)
{
   r.Write(out,", ",", ");
   return out;
}

istream& operator >> (istream& in, Matrix& r)
{
   cout << "Introdurre matrice (" << r.rows() << "," << r.columns() << ")"
        << " elementi" << endl;
   for (int i=0;i<r.rows();i++) {
      in >> (*r.data[i]);
   }
   return in;
}

Matrix& Matrix::operator = (const Matrix& R) {
   rws=R.rws;
   clms=R.clms;
   Vector BUF(clms);
   data=new Vector*[rws];
   for (int i=0;i<rws;i++) data[i]=new Vector(clms);
   for (int i=0;i<rws;i++) (*data[i])=R[i];
   return *this;
}

unsigned int Matrix::rows() const{ return rws; }
unsigned int Matrix::columns() const { return clms; }

void Matrix::Display() const
{
   for (int i=0; i<rows(); i++) {
     for (int j=0;j<columns();j++) printf ("%-9.5lg ",(*data[i])[j]);
   printf ("\n");
   }
}

void Matrix::Clear()
{
   for (int i=0;i<rows();i++) {
      for (int j=0;j<columns();j++)
	(*data[i])[j]=0;
   }
}

ostream& Matrix::Write(ostream& out, const char* fill,
		       const char* fill2,
		       const char* indent) const
{
   enum Matrix_type {
      GENERIC_MATRIX,
	NULL_MATRIX,
	DIAG_MATRIX,
	EYE_MATRIX
   };
   
   Boolean TARGET=Y;
   Matrix_type TYPE=GENERIC_MATRIX;
   
   if (EXTENDED_MATRIX_DISPLAY==N)
     {
	/* DIAG MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ((i!=j) && ((*data[i])[j]!=0)) TARGET=N;
	if (TARGET==Y) TYPE=DIAG_MATRIX;
	/* NULL MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ( (*data[i])[j]!=0 ) TARGET=N;
	if (TARGET==Y) TYPE=NULL_MATRIX;
	/* EYE MATRIX */
	if (TYPE==DIAG_MATRIX) {
	   TARGET=Y;
	   for (int i=0;i<rows();i++)
	      if ((*data[i])[i]!=1) TARGET=N;
	   if (TARGET==Y) TYPE=EYE_MATRIX;
	}
	/* VISUALIZZAZIONE SELETTIVA */
	switch (TYPE) {
	 case NULL_MATRIX: out << indent << "null"; break;
	 case EYE_MATRIX: out << indent << "eye"; break;
	 case DIAG_MATRIX: out << indent << "diag";
	   for (int i=0;i<rows();i++) 
	     out << fill << (*(data[i]))[i];
	   break;
	}
        if (TYPE!=GENERIC_MATRIX) return out;
     }
   for (int i=0;i<rows();i++) {
      out << indent;
      out << "",(*data[i]).Write(out,fill);
      if (i<rows()-1) out << "" << fill2;
   }
   return out;
}

Vector Matrix::GetVec (unsigned int n) const
{
   // Ritorna la riga n-esima : indice a base 1
   return (*data[(n-1)]);
}

Matrix operator * (const Matrix& A, const Matrix& B)
{
   unsigned int w=B.columns();
   unsigned int h=A.rows();
   double psum=0;
   if (A.columns()!=B.rows()) DEBUGCOUT ("Errore in moltiplicazione matrice");
   Matrix P(h,w);
   for (int i=0;i<h;i++)
     for (int j=0;j<w;j++) {
	psum=0;
	for (int k=0; k<A.columns();k++) psum=psum+A[i][k]*B[k][j];
	P[i][j]=psum;
     }
   return P;
}

Matrix operator * (const double value, const Matrix& B)
{
   unsigned int r=B.rows();
   unsigned int c=B.columns();
   Matrix P(r,c);
   for (int i=0;i<B.rows();i++)
     for (int j=0;j<B.columns();j++)
       P[i][j]=(B[i][j])*value;
   return P;
}

Matrix operator * (const Matrix& B, const double value)
{
   Matrix P(B.rows(),B.columns());
   P=value*B;
   return P;
}

Vector operator * (const Matrix& M, const Vector& V)
{
   if (V.size()!=M.columns()) DEBUGCOUT ("Matrice x vettore = dimensioni diverse");
   Vector result (V.size());
   Matrix A = M;
   Matrix B (M.columns(),1);
   for (int i=0;i<V.size();i++) B[i][0]=V[i];
   Matrix C (M.columns(),1);
   C = A*B;
   for (int i=0;i<V.size();i++) result[i]=C[i][0];
   return result;
}

Vector operator * (const Vector& V, const Matrix& M)
{
   if (V.size()!=M.rows()) DEBUGCOUT ("Matrice x vettore = dimensioni diverse!");
   Vector result (V.size());
   Matrix A = M;
   Matrix B (1,M.columns());
   for (int i=0;i<V.size();i++) B[0][i]=V[i];
   Matrix C (1,M.columns());
   C=B*A;
   for (int i=0;i<V.size();i++) result[i]=C[0][i];
   return result;
}

Matrix operator + (const Matrix& A, const Matrix& B)
{
   unsigned int r=A.rows();
   unsigned int c=A.columns();
   Matrix P(r,c);
   for (int i=0;i<r;i++)
     for (int j=0;j<c;j++)
       P[i][j]=A[i][j]+B[i][j];
   return P;
}

Matrix operator - (const Matrix& A, const Matrix& B)
{   
   unsigned int r=A.rows();
   unsigned int c=A.columns();
   Matrix P(r,c);
   P=A+(-1)*B;
   return P;
}

// Estrazione di sottomatrice (matrice, riga iniz, colonna iniz, dimensioni)

Matrix SubMat (const Matrix& R, unsigned int sr=1, unsigned int sc=1,
	       unsigned int rl=0, unsigned int cl=0)
{
   unsigned int r=R.rows();
   unsigned int c=R.columns();
   sr--; sc--;
   if (rl==0) rl=r-sr;
   if (cl==0) cl=c-sc;
   Matrix SUB(rl,cl);
   for (int i=0;i<rl;i++)
     for (int j=0;j<cl;j++)
       SUB[i][j]=R[sr+i][sc+j];
   return SUB;
}

// Interpolazione di matrici

Matrix Interp (const Matrix& A, const Matrix& B, double P=0.5)
{
   if ((A.rows() != B.rows()) | (A.columns() != B.columns())) {
      DEBUGCOUT ("Interpolazione non possibile !");
      Mat3x3 R = Zero3x3;
      return R;
   }
   Matrix SUB(A.rows(),A.columns());
   /* Interpolazione lineare a una posizione compresa tra 0 e 1 */
   SUB=P*B+(1-P)*A;
   return SUB;
}

Vector Interp (const Vector& A, const Vector& B, double P=0.5)
{
   if (A.size()!=B.size()) {
      DEBUGCOUT ("Interpolazione non possibile!");
      Vector R(A.size());
      return R;
   }
   Vector R(A.size());
   R=P*B+(1-P)*A;
   return R;
}

Vec3 Interp (const Vec3& A, const Vec3& B, double P=0.5)
{
  Vector C(3),D(3);
  Vec3 E;
  E=Interp ( (Vector) A, (Vector) B,P);
  return E;
}

Mat3x3 Interp (const Mat3x3& A, const Mat3x3& B, double P=0.5)
{
   Mat3x3 SUB;
   /* Solo interpolazione lineare per ora */
   SUB=P*B+(1-P)*A;
   return SUB;
}

// Determinante di una matrice quadrata

double Det (Matrix& R)
{
   unsigned int r=R.rows();
   unsigned int c=R.columns();
   int I,idx,idy;
   double ValRet=0;
   if (r!=c) {
      return 0;
   }
   else {
      if (r==2) {
	 return ((R[0][0]*R[1][1]-R[0][1]*R[1][0]));
      }
      for (I=0;I<c;I++) {
	 Matrix SUBMAT (r-1,c-1);
	 for (idy=0;idy<I;idy++) {
	    for (idx=0;idx<(r-1);idx++) {
	       SUBMAT[idx][idy]=R[idx+1][idy];
	    }
	 }
	 for (idy=I+1;idy<c;idy++) {
	    for (idx=0;idx<(r-1);idx++) {
	       SUBMAT[idx][idy-1]=R[idx+1][idy];
	    }
	 }
	 double coef=R[0][I];
	 if ((I%2)!=0) coef=-1*coef;
	 ValRet=ValRet+coef*Det(SUBMAT);
      }
   }
   return ValRet;
}

Vec3& Vec3::operator = (const Vector& R) {
   if (R.size()!=3) cout << "ILLEGAL OPERATION" << endl;
   for (int i=0;i<3;i++) dato[i]=R.dato[i];
   return *this;
}
   
Mat3x3& Mat3x3::operator = (const Matrix& R) {
   rws=R.rws;
   clms=R.clms;
   Vector BUF(clms);
   data=new Vector*[rws];
   for (int i=0;i<rws;i++) data[i]=new Vector(clms);
   for (int i=0;i<rws;i++) (*data[i])=R[i];
   return *this;
}

Mat6x6& Mat6x6::operator = (const Matrix& R) {
   rws=R.rws;
   clms=R.clms;
   Vector BUF(clms);
   data=new Vector*[rws];
   for (int i=0;i<rws;i++) data[i]=new Vector(clms);
   for (int i=0;i<rws;i++) (*data[i])=R[i];
   return *this;
}

/* Rwrite writes the rotation matrix in the form
 * 1, Vec3, 2, Vec3 */

ostream& Mat3x3::RWrite(ostream& out, const char* fill,
		       const char* fill2,
		       const char* indent) const
{
   enum Matrix_type {
      GENERIC_MATRIX,
	NULL_MATRIX,
	DIAG_MATRIX,
	EYE_MATRIX
   };
   
   Boolean TARGET=Y;
   Matrix_type TYPE=GENERIC_MATRIX;
   
   if (EXTENDED_MATRIX_DISPLAY==N)
     {
	/* DIAG MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ((i!=j) && ((*data[i])[j]!=0)) TARGET=N;
	if (TARGET==Y) TYPE=DIAG_MATRIX;
	/* NULL MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ( (*data[i])[j]!=0 ) TARGET=N;
	if (TARGET==Y) TYPE=NULL_MATRIX;
	/* EYE MATRIX */
	if (TYPE==DIAG_MATRIX) {
	   TARGET=Y;
	   for (int i=0;i<rows();i++)
	      if ((*data[i])[i]!=1) TARGET=N;
	   if (TARGET==Y) TYPE=EYE_MATRIX;
	}
	/* VISUALIZZAZIONE SELETTIVA */
	switch (TYPE) {
	 case NULL_MATRIX: out << indent << "null"; break;
	 case EYE_MATRIX: out << indent << "eye"; break;
	 case DIAG_MATRIX: out << indent << "diag";
	   for (int i=0;i<rows();i++) 
	     out << fill << (*(data[i]))[i];
	   break;
	}
        if (TYPE!=GENERIC_MATRIX) return out;
     }
     out << indent << "1, ",(*data[0]).Write(out,fill) << ", "
     << fill2 << indent << "2, ",(*data[1]).Write(out,fill);
   return out;
}


/* FUNZIONI MISTE */

double Module (Vector P)
{
   /* La funzione restituisce il modulo di un vettore */
   double m=0;
   for (int i=0;i<P.size();i++)
     m=m+(P.dato[i]*P.dato[i]);
   m=sqrt(m);
   return m;
}

Vector Direction (Vector P)
{
   /* La funzione restituisce i coseni direttori di un vettore */
   double rho=Module(P);
   Vector R (P.size());
   if (rho==0) return R;
   /* Se il modulo è nullo è impossibile ottenere la direzione! */
   for (int i=0;i<P.size();i++)
     R[i]=(P.dato[i]/rho);
   return R;
}

double Cos_b2Vector (Vector V1,Vector V2)
{
   double A=0;
   Vector D1,D2;
   D1=Direction(V1); D2=Direction(V2);
   for (int i=0;i<V1.size();i++)
     A=A+D1.dato[i]*D2.dato[i];
   return A;
}

/* Prodotto vettore tra il vettore stesso e il vettore v */
Vec3 Vec3::Cross(const Vec3& v) const {     
      return Vec3(dato[1]*v.dato[2]-dato[2]*v.dato[1],
		  dato[2]*v.dato[0]-dato[0]*v.dato[2],
		  dato[0]*v.dato[1]-dato[1]*v.dato[0]);
};
   
/* Prodotto vettore per matrice. */ 

   
/* Prodotto scalare per un vettore e sè stesso */
doublereal Vec3::Dot(const Vec3& v) const {
      return 
	dato[0]*v.dato[0]+
	dato[1]*v.dato[1]+
	dato[2]*v.dato[2];
};
   
doublereal Vec3::Dot(void) const { 
   return 
	dato[0]*dato[0]+
	dato[1]*dato[1]+
	dato[2]*dato[2];
};
   
/* Norma del vettore */
doublereal Vec3::Norm(void) const { 
      return 
	sqrt(dato[0]*dato[0]+
	     dato[1]*dato[1]+
	     dato[2]*dato[2]);
};

double Vector::dGet(int i) const
{
   return (dato[i-1]);
}

double Matrix::dGet(int i, int j) const
{
   return (*(data[i-1]))[j-1];
}

MBDyn_entity::MBDyn_entity ()
{
   Mode=NUL;
   RS=0;
}

MBDyn_entity::MBDyn_entity (Mode_Type ENT)
{
   Mode=ENT;
   RS=0;
}

MBDyn_entity::MBDyn_entity (Mode_Type ENT, Id RF)
{
   Mode=ENT;
   RS=RF;
}

MBDyn_entity::~MBDyn_entity () {}

ostream& operator << (ostream& out, const MBDyn_entity& P)
{
   switch (P.Mode) {
    case (MBDyn_entity::NUL): out << ""; break;
    case (MBDyn_entity::GLOBAL): out << "reference, global, "; break;
    case (MBDyn_entity::LOCAL): out << "reference, local, "; break;
    case (MBDyn_entity::NODE): out << "reference, node, "; break;
    case (MBDyn_entity::REFERENCED): out << "reference, " << P.RS 
	<< ", "; break;
    default: out << ""; break;
   }
   return out;
}

ostream& MBDyn_entity::Restart (ostream& out) const
{
   out << (*this);
   return out;
}

RVec3::RVec3  (double i, double j, double k) : Vec3(i,j,k) {}
RVec3::RVec3  (Vec3 V) : Vec3 (V), REF (MBDyn_entity::NUL) {}
RVec3::RVec3  (Vec3 V, MBDyn_entity E) : Vec3(V),REF(E) {}
RVec3::RVec3  (Vec3 V, 
	       MBDyn_entity::Mode_Type MTE,
	       Id RF) : Vec3(V), REF(MTE,RF) {}

RVec3::~RVec3 () {}

ostream& RVec3::Restart (ostream& out) const
{
   REF.Restart(out);
   out << "",Write(out,", ");
   return out;
}

RMat3x3::RMat3x3 (Mat3x3 M) : Mat3x3 (M) {}
RMat3x3::RMat3x3 (Mat3x3 M, MBDyn_entity E) : Mat3x3(M),REF(E) {}
RMat3x3::RMat3x3 (Mat3x3 M,
		  MBDyn_entity::Mode_Type MTE,
		  Id RF) : Mat3x3 (M), REF(MTE,RF) {}

RMat3x3::~RMat3x3 () {}

RMat3x3& RMat3x3::operator = (const Matrix& R) {
   rws=R.rws;
   clms=R.clms;
   Vector BUF(clms);
   data=new Vector*[rws];
   for (int i=0;i<rws;i++) data[i]=new Vector(clms);
   for (int i=0;i<rws;i++) (*data[i])=R[i];
   return *this;
}


ostream& RMat3x3::Restart(ostream& out) const
{
   REF.Restart(out);
   out << "",Write(out,", ","\n");
   return out;
}

ostream& operator << (ostream& out, const RVec3& P)
{
   P.Restart(out);
   return out;
}

ostream& operator << (ostream& out, const RMat3x3 P)
{
   P.Restart(out);
   return out;
}

ostream& RVec3::Write (ostream& out, const char* fill) const
{
   /* Check if vector is null */
   REF.Restart(out);
   Boolean IsNull=Y;
   for (int i=0; i<sz; i++) if (dato[i] != 0) IsNull=N;
   if ((IsNull) && (EXTENDED_MATRIX_DISPLAY==N)) out << "null";
   else {
   for (int i=0; i<(sz-1); i++) out << dato[i] << fill;
      out << dato[sz-1];
   }
   return out;
} 

ostream& RMat3x3::Write (ostream& out, const char* fill,
			 const char* fill2,
			 const char* indent) const
{
      enum Matrix_type {
      GENERIC_MATRIX,
	NULL_MATRIX,
	DIAG_MATRIX,
	EYE_MATRIX
   };
   
   Boolean TARGET=Y;
   Matrix_type TYPE=GENERIC_MATRIX;
   
   if (EXTENDED_MATRIX_DISPLAY==N)
     {
	/* DIAG MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ((i!=j) && ((*data[i])[j]!=0)) TARGET=N;
	if (TARGET==Y) TYPE=DIAG_MATRIX;
	/* NULL MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ( (*data[i])[j]!=0 ) TARGET=N;
	if (TARGET==Y) TYPE=NULL_MATRIX;
	/* EYE MATRIX */
	if (TYPE==DIAG_MATRIX) {
	   TARGET=Y;
	   for (int i=0;i<rows();i++)
	      if ((*data[i])[i]!=1) TARGET=N;
	   if (TARGET==Y) TYPE=EYE_MATRIX;
	}
	/* VISUALIZZAZIONE SELETTIVA */
	switch (TYPE) {
	 case NULL_MATRIX: out << indent << REF << "null"; break;
	 case EYE_MATRIX: out << indent << REF << "eye"; break;
	 case DIAG_MATRIX: out << indent << REF << "diag";
	   for (int i=0;i<rows();i++) 
	     out << fill << (*(data[i]))[i];
	   break;
	}
        if (TYPE!=GENERIC_MATRIX) return out;
     }
   REF.Restart(out);
   for (int i=0;i<rows();i++) {
      out << indent;
      out << "",(*data[i]).Write(out,fill);
      if (i<rows()-1) out << "" << fill2;
   }
   return out;
}

ostream& RMat3x3::RWrite(ostream& out, const char* fill,
		       const char* fill2,
		       const char* indent) const
{
   enum Matrix_type {
      GENERIC_MATRIX,
	NULL_MATRIX,
	DIAG_MATRIX,
	EYE_MATRIX
   };
   
   Boolean TARGET=Y;
   Matrix_type TYPE=GENERIC_MATRIX;
   
   if (EXTENDED_MATRIX_DISPLAY==N)
     {
	/* DIAG MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ((i!=j) && ((*data[i])[j]!=0)) TARGET=N;
	if (TARGET==Y) TYPE=DIAG_MATRIX;
	/* NULL MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ( (*data[i])[j]!=0 ) TARGET=N;
	if (TARGET==Y) TYPE=NULL_MATRIX;
	/* EYE MATRIX */
	if (TYPE==DIAG_MATRIX) {
	   TARGET=Y;
	   for (int i=0;i<rows();i++)
	      if ((*data[i])[i]!=1) TARGET=N;
	   if (TARGET==Y) TYPE=EYE_MATRIX;
	}
	/* VISUALIZZAZIONE SELETTIVA */
	switch (TYPE) {
	 case NULL_MATRIX: out << indent << REF << "null"; break;
	 case EYE_MATRIX: out << indent << REF << "eye"; break;
	 case DIAG_MATRIX: out << indent << REF << "diag";
	   for (int i=0;i<rows();i++) 
	     out << fill << (*(data[i]))[i];
	   break;
	}
        if (TYPE!=GENERIC_MATRIX) return out;
     }
     out << indent << REF << "1, ",(*data[0]).Write(out,fill) << ", "
     << fill2 << indent << "2, ",(*data[1]).Write(out,fill);
   return out;
}



// FINE DELLA PARTE DI CODICE SPERIMENTALE

Mat3x3 Unref (const RMat3x3 A)
{
   Mat3x3 R;
   for (int i=0;i<3;i++)
    for (int j=0;j<3;j++)
       R[i][j]=A[i][j];
   return R;
   
}

Vec3 Unref (const RVec3& A)
{
   Vec3 R;
   for (int i=0;i<3;i++)
     R[i]=A[i];
   return R;
}

// MATRICE INTERPOLANTE

RMat3x3 Interp (const RMat3x3& A, const RMat3x3& B, double P=0.5)
{
   if (A.REF!=B.REF) {
      cout << "DEBUG: [RMAT3X3 INTERP] = Incongruent reference systems" << endl;
      exit (-1);
   }
   RMat3x3 SUB;
   /* Solo interpolazione lineare per ora */
   SUB=P*B+(1-P)*A;
   SUB.REF=A.REF;
   return SUB;
}

RVec3 Interp (const RVec3& A, const RVec3& B, double P=0.5)
{
   if (A.REF!=B.REF) {
      cout << "DEBUG: [RMAT3X3 INTERP] = Incongruent reference systems" << endl;
      exit(-1);
   }
   Vec3 MR;
   MR=P*B+(1-P)*A;
   RVec3 R(MR,A.REF);
   return R;
}

MBDyn_entity& MBDyn_entity::operator = (const MBDyn_entity& P)
{
   Mode=P.Mode;
   RS=P.RS;
   return (*this);
}

Boolean MBDyn_entity::operator != (const MBDyn_entity& P) const
{
   if ((Mode!=P.Mode) | (RS!=P.RS)) return Y; else return N;
}

Boolean Vec3::operator != (const Vec3& P) const
{
   unsigned int count=0;
   for (int i=0;i<3;i++)
     if (dato[i]!=P.dato[i]) return Y;
   /* .PHONY CALL, not reachable if the argument is different */
   return N;
}

Boolean Mat3x3::operator != (const Mat3x3& P) const
{
   unsigned int count=0;
   for (int i=0; i<3; i++)
     for (int j=0; j<3; j++)
       if (((*data[i])[j])!=P[i][j]) return Y;
   /* .PHONY CALL, not reachable if the argument is different */
   return N;
}

Boolean Matrix::operator != (const Matrix& P) const
{
   unsigned int count=0;
   for (int i=0;i<rows();i++)
     for (int j=0;j<columns();j++)
       if (((*data[i])[j])!=P[i][j]) return Y;
   return N;
}

RVec3 RMat3x3::GetVec (unsigned int n) const
{
   // Ritorna la riga n-esima : indice a base 1
   Vec3 parteuno;
   MBDyn_entity WS;
   parteuno = (*data[(n-1)]);
   WS=REF;
   RVec3 toreturn(parteuno,WS);
   return (toreturn);
}

ostream& Mat6x6::Write(ostream& out, const char* fill,
		       const char* fill2,
		       const char* indent) const
{
   enum Matrix_type {
      GENERIC_MATRIX,
	NULL_MATRIX,
	DIAG_MATRIX,
	EYE_MATRIX
   };
   
   Boolean TARGET=Y;
   Matrix_type TYPE=GENERIC_MATRIX;
   
   if (EXTENDED_MATRIX_DISPLAY==N)
     {
	/* DIAG MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ((i!=j) && ((*data[i])[j]!=0)) TARGET=N;
	if (TARGET==Y) TYPE=DIAG_MATRIX;
	/* NULL MATRIX */
	for (int i=0;i<rows();i++)
	  for (int j=0;j<columns();j++)
	    if ( (*data[i])[j]!=0 ) TARGET=N;
	if (TARGET==Y) TYPE=NULL_MATRIX;
	/* EYE MATRIX */
	if (TYPE==DIAG_MATRIX) {
	   TARGET=Y;
	   for (int i=0;i<rows();i++)
	      if ((*data[i])[i]!=1) TARGET=N;
	   if (TARGET==Y) TYPE=EYE_MATRIX;
	}
	/* VISUALIZZAZIONE SELETTIVA */
	switch (TYPE) {
	 case NULL_MATRIX: out << indent << "null"; break;
	 case EYE_MATRIX: out << indent << "eye"; break;
	 case DIAG_MATRIX: out << indent << "diag";
	   for (int i=0;i<rows();i++) 
	     out << fill << (*(data[i]))[i];
	   break;
	}
        if (TYPE!=GENERIC_MATRIX) return out;
     }
   for (int i=0;i<rows();i++) {
      out << indent;
      out << "",(*data[i]).Write(out,fill);
      if (i<rows()-1) out << "" << fill2;
   }
   return out;
}
