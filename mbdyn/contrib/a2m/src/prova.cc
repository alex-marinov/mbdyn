


#include <matrix.h>
#include <iostream.h>
#include <math.h>

Vec3 a(1,0,0),b(50,50,0);
double k1,k2;
Eye PR(3);
Matrix B (3,1);

main ()
{
   const char indent[]="       ";
   cout << "Hello!" << endl;
   cout << "Vectors are:" << a << ", " << b << endl;
   k1=a*b;
   cout << "Prodotto scalare:" << k1 << endl;
   cout << "Coseno angolo:" << Cos_b2Vector(a,b) << endl;
   cout << "Angolo:" << (57.3*acos(Cos_b2Vector(a,b))) << endl;
   Vector Q(3);
   Q=a^b;
   cout << "Prodotto vettore:" << Q << endl;
   cout << "Risultato dell'interpolazione" << endl;
   cout << a << endl << b << endl;
   cout << Interp(a,b) << endl;
   cout << (Interp(a,b)/25.5) << endl;
   cout << "",PR.Write (cout,",","\n",indent);
   cout << endl;
}
