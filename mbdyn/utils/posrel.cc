#include "matvec3.h"
#include "string.h"


int main(int argn, const char* const argv[])
{
   static doublereal d[3];
   while (1) {
      cin >> d[0];
      if (cin) {
	 cin >> d[1] >> d[2];
	 Vec3 x1(d);
	 cin >> d[0] >> d[1] >> d[2];
	 Mat3x3 R1(RFromEulerAngles(Vec3(d)/180.*M_PI));
	 cin >> d[0] >> d[1] >> d[2];
	 Vec3 x2(d);

	 cout << R1.Transpose()*(x2-x1) << endl;
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}
