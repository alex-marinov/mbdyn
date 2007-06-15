/* MATRIX - Libreria contenente funzioni e definizioni di matrice */

#ifndef MATRIX_H
#define MATRIX_H

// Valore 0: Indirizzamento codice, primo elemento = 0;

#include <stdio.h>
#include <iostream>
#include <defs.h>

// Modello di classe vettore

struct Vector {
   Vector (unsigned int n=8);
   Vector (const Vector& );
   ~Vector ();
   void copy(const Vector&);
   double& operator [] (unsigned int) const;
   Vector& operator=(const Vector&);
   ostream& Write (ostream&, const char*) const;
   friend ostream& operator << (ostream&, const Vector&);
   friend istream& operator >> (istream&, Vector&);
   friend Vector operator + (const Vector&, const Vector&);
   friend Vector operator - (const Vector&, const Vector&);
   friend double operator * (const Vector&, const Vector&);
   friend Vector operator * (const Vector&, const double);
   friend Vector operator * (const double, const Vector&);
   friend Vector operator ^ (const Vector&, const Vector&);
   friend Vector operator / (const Vector&, const double);
   friend Vector& operator /= (Vector&,const double);
   void Clear();
   unsigned int size() const;
   double Module() const;   
   // Dati propri della classe
   unsigned int sz;
   double* dato;
   double dGet(int) const;
};

struct Matrix {
   Matrix (unsigned int n=6, unsigned int m=6);
   Matrix (const Matrix&);
   ~Matrix ();
   Vector& operator [] (unsigned int) const;
   Matrix tr();
   friend ostream& operator << (ostream&, const Matrix&);
   friend istream& operator >> (istream&, Matrix&);
   friend Matrix operator * (const Matrix&, const Matrix&);
   friend Matrix operator + (const Matrix&, const Matrix&);
   friend Matrix operator * (const double, const Matrix&);
   friend Matrix operator * (const Matrix&, const double);
   friend Matrix operator - (const Matrix&, const Matrix&);
   friend Matrix operator / (const Matrix&, const Matrix&);
   Boolean operator != (const Matrix&) const;
   // Quest'ultima sopra è dichiarata in math.h, riferimento all'inversione
   Matrix& operator = (const Matrix&);
   unsigned int rows() const;
   unsigned int columns() const;
   void Display () const;
   void Clear();
   /*ostream& Write(ostream&, const char*,const char*) const;*/
   ostream& Write(ostream&, const char*, const char*, const char* c="") const;
   Vector GetVec(unsigned int) const;
   // Dati propri della classe
   Vector** data;
   unsigned int rws;
   unsigned int clms;
   double dGet(int,int) const;
};

struct Eye : public Matrix {
   Eye(unsigned int n=6) : Matrix(n,n)
     {
	for (int i=0;i<n;i++) (*data[i])[i]=1;
     }	  
   ~Eye() {} ;
   Eye& operator = (const Matrix&);
};

struct Eye3x3 : public Eye {
   Eye3x3() : Eye(3) {}
   ~Eye3x3() {}
};

struct Zeros : public Matrix {
   Zeros(unsigned int m=6, unsigned int n=6) : Matrix(m,n) {}
   ~Zeros() {}
   Zeros& operator = (const Matrix&);
};

/* CLASSI DI MATRICI E VETTORI PIU' UTILIZZATI */

struct Vec3 : public Vector {
   Vec3(double i=0, double j=0, double k=0) : Vector(3) {
      dato[0]=i; dato[1]=j; dato[2]=k; }
   ~Vec3() {}
   void Set (double,double,double);
   /* Importate dalla libreria matvec3 di MBDyn */
   Vec3 Cross(const Vec3&) const;
   doublereal Dot(void) const;
   doublereal Dot(const Vec3&) const;
   doublereal Norm(void) const;
   Vec3& operator = (const Vector&);
   Boolean operator != (const Vec3&) const;
};

struct Mat3x3: public Matrix {
   Mat3x3() : Matrix(3,3) {}
   Mat3x3 (Vec3 a1, Vec3 b1, Vec3 c1) : Matrix(3,3) {
      (*data[0])=a1; (*data[1])=b1; (*data[2])=c1; 
   }
   Mat3x3(double c11, double c12, double c13,
	  double c21, double c22, double c23,
	  double c31, double c32, double c33) : Matrix(3,3) {
	  // Costruttore esteso 
	  (*data[0])[0]=c11;
	  (*data[0])[1]=c12;
	  (*data[0])[2]=c13;
	  (*data[1])[0]=c21;
	  (*data[1])[1]=c22;
	  (*data[1])[2]=c23;
	  (*data[2])[0]=c31;
	  (*data[2])[1]=c32;
	  (*data[2])[2]=c33;
   }
   ~Mat3x3() {}
   Mat3x3& operator = (const Matrix&);
   Boolean operator != (const Mat3x3&) const;
   Boolean operator != (const Eye3x3&) const;
   ostream& RWrite(ostream&,const char*, const char*, const char* c="") const;
};

struct Mat6x6: public Matrix {
   Mat6x6() : Matrix (6,6) {}
   ~Mat6x6() {}
   Mat6x6& operator = (const Matrix&);
};

struct Vec2 : public Vector {
   Vec2(double i=0, double j=0) : Vector(2) {
       dato[0]=i; dato[1]=j; }
   ~Vec2() {}
   void Set (double,double);
};

struct Mat2x2: public Matrix {
   Mat2x2() : Matrix(2,2) {}
   ~Mat2x2() {}
};

// Funzioni per la gestione delle sottomatrici

Matrix SubMat(const Matrix&, unsigned int, unsigned int,
	      unsigned int, unsigned int);

// Funzione determinante

double Det(Matrix&);

// Matrice interpolante e vettore interpolante

Matrix Interp(const Matrix&, const Matrix&,double p=0.5);
Vector Interp(const Vector&, const Vector&, double p=0.5);
Vec3 Interp (const Vec3&, const Vec3&, double p=0.5);
Mat3x3 Interp (const Mat3x3&, const Mat3x3&, double p=0.5);

// Funzione miste di utilizzo sui vettori

double Module (Vector);
Vector Direction (Vector);
double Cos_b2Vector(Vector,Vector);

// Entità di MBDYN - struttura contenente la parte di riferimento

struct MBDyn_entity {
   enum Mode_Type {
      NUL,
	GLOBAL,
	LOCAL,
	NODE,
	REFERENCED
   };
   Mode_Type Mode;
   Id RS;
   ostream& Restart (ostream& out) const;
   MBDyn_entity ();
   MBDyn_entity (Mode_Type);
   MBDyn_entity (Mode_Type,Id);
   ~MBDyn_entity ();
   friend ostream& operator << (ostream&, const MBDyn_entity&);
   MBDyn_entity& operator = (const MBDyn_entity&);
   Boolean operator != (const MBDyn_entity&) const ;
};

// PARTE DI CODICE SPERIMENTALE [ ENTITA' REFERENZIATE ]

struct RVec3 : public Vec3 {
   /* Referenced Vec3 : contiene un Vec3+un'entita che determina il rif. */
   RVec3 (): Vec3(0,0,0) {}
   RVec3 (Vec3,MBDyn_entity);
   RVec3 (Vec3,MBDyn_entity::Mode_Type,Id);
   RVec3 (Vec3);
   RVec3 (double i, double j, double k);
   ~RVec3 ();
   ostream& Restart (ostream& out) const;
   friend ostream& operator << (ostream&, const RVec3&);
   ostream& Write (ostream&, const char*) const;
   
   MBDyn_entity REF;
};

/* L'utilità di Rmat è ancora da verificare .. */
struct RMat3x3 : public Mat3x3 {
   /* Referenced Mat3x3 : contiene un Mat3x3 + il riferimento */
   RMat3x3 () : Mat3x3 () {}
   RMat3x3 (Mat3x3,MBDyn_entity);
   RMat3x3 (Mat3x3,MBDyn_entity::Mode_Type,Id);
   RMat3x3 (Mat3x3);
   ~RMat3x3 ();
   ostream& Restart (ostream& out) const;
   friend ostream& operator << (ostream&, const RMat3x3&);
   RVec3 GetVec(unsigned int) const;
   ostream& Write(ostream&, const char*, const char*, const char* c="") const;
   ostream& RWrite(ostream&,const char*, const char*, const char* c="") const;
   RMat3x3& operator = (const Matrix&);
   
   MBDyn_entity REF;
};

// FINE DELLA PARTE DI CODICE SPERIMENTALE

RMat3x3 Interp (const RMat3x3&, const RMat3x3&, double P=0.5);
RVec3 Interp (const RVec3&, const RVec3&, double P=0.5);
Mat3x3 Unref (const RMat3x3);
Vec3 Unref (const RVec3&);

#endif
