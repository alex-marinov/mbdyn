#ifndef MATHEM_H
#define MATHEM_H

#include <defs.h>
#include <matrix.h>
#include <math.h>
#include <float.h>

/* Usage: Gauss (Vector& solution, Matrix A orig, vector Yor) */

Mat3x3 MomentToInertia(Vec3, Vec3);

/* Parametri di rotazione di Eulero a partire dalla matrice R */
Vec3 gparam (const Mat3x3&);

double square(double);
double a360tan2 (double,double);
void Gauss (Vector&, Matrix&, Vector&);
Matrix Inv (Matrix&);
Matrix operator / (Matrix&, Matrix&);
Matrix& KMatrix (double, double, double, double,
		 double, double, double, double, double,
		 double*, double);
Vec3 EulerAngles (const Mat3x3&); 
Mat3x3 RFromEulerAngles (const Vec3&);
Mat3x3 MatR2vec (unsigned short int, const Vec3&,
		 unsigned short int, const Vec3&);

#endif
