#include "matvec3.h"
#include "string.h"


int main(int argn, const char* const argv[])
{ 
   if (argn > 1) {
      if (!strcasecmp(argv[1], "-?")
	  || !strcasecmp(argv[1], "-h") 
	  || !strcasecmp(argv[1], "--help")) {
	 cerr << endl << "usage: " << argv[0] << endl << endl
	   << "    reads the Euler angles (in degs) from stdin;" << endl
	   << "    writes the rotation matrix (column-oriented) on standard output" << endl << endl
	   << "part of MBDyn package (Copyright (C) Pierangelo Masarati, 1996)" << endl << endl;
	 exit(EXIT_SUCCESS);
      }
   }   

   static doublereal d[3];
   while (1) {
      cin >> d[0];
      if (cin) {
	 cin >> d[1] >> d[2];
	 cout << RFromEulerAngles(Vec3(d)/180.*M_PI) << endl;
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}
