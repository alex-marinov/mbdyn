#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <matvec3.h>
#include <string.h>


int 
main(int argn, const char* const argv[])
{ 
   	if (argn > 1) {
      		if (!strcasecmp(argv[1], "-?")
	  	    || !strcasecmp(argv[1], "-h") 
		    || !strcasecmp(argv[1], "--help")) {
	 		cerr << endl 
				<< "usage: " << argv[0] << endl 
				<< endl
	   			<< "    reads the Euler angles (in degs)"
				" from stdin;" << endl
	   			<< "    writes the rotation matrix"
				" on standard output" << endl
				<< "    (m11, m12, m13,"
				" m21, m22, m23,"
				" m31, m32, m33)" << endl
				<< endl
	   			<< "part of MBDyn package (Copyright (C)"
				" Pierangelo Masarati, 1996-2000)" << endl
				<< endl;
	 		exit(EXIT_SUCCESS);
      		}
   	}   

   	static doublereal d[3];
   	while (1) {
      		cin >> d[0];
      		if (cin) {
	 		cin >> d[1] >> d[2];
	 		cout << EulerAngles2MatR(Vec3(d)/180.*M_PI) << endl;
      		} else {
	 		break;
      		}
   	}
   
   	return EXIT_SUCCESS;
}

