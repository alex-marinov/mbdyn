#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <matvec3.h>
#include <string.h>


int main(int argn, const char* const argv[])
{
   if (argn > 1) {
      if (!strcasecmp(argv[1], "-?")
	  || !strcasecmp(argv[1], "-h") 
	  || !strcasecmp(argv[1], "--help")) {
	 cerr << endl << "usage: " << argv[0] << endl << endl
	   << "    reads a rotation matrix (row-oriented) from stdin;" << endl
	   << "    writes the Euler angles (in degs) on standard output" << endl << endl
	   << "part of MBDyn package (Copyright (C) Pierangelo Masarati, 1996)" << endl << endl;
	 exit(EXIT_SUCCESS);
      }
   }   

   static doublereal d[9];
   while (1) {
      cin >> d[0];
      if (cin) {
	 Vec3 e;
	 doublereal e0;
	 cin >> d[3] >> d[6] >> d[1] >> d[4] >> d[7] >> d[2] >> d[5] >> d[8];
	 EulerParams(Mat3x3(d, 3), e0, e);
	 cout << e0 << " " << e << endl;
      } else {
	 break;
      }      
   }   
   
   return (EXIT_SUCCESS);
}
