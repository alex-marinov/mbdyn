#ifdef HAVE_CONFIG_H
#include <mbconfig.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>

#include <myassert.h>
#include <mynewmem.h>

#if defined(USE_AERODYNAMIC_ELEMS)
#include <aerodata.h>
#include <aerodc81.h>
#include <c81data.h>
#endif /* USE_AERODYNAMIC_ELEMS */

/*
 * uso: testc81 file alpha mach
 */
int
main(int argc, char *argv[])
{
#if defined(USE_AERODYNAMIC_ELEMS)
	if (argc != 4) {
		cerr << "usage: testc81 file alpha mach" << endl;
		exit(EXIT_SUCCESS);
	}

	ifstream in(argv[1]);
	if (!in) {
		cerr << "unable to open file '" << argv[1] << "'" << endl;
		exit(EXIT_FAILURE);
	}

	c81_data* data = new C81Data(1);

	if (read_c81_data(in, data)) {
		cerr << "unable to read c81 data from file '" 
			<< argv[1] << "'" << endl;
		exit(EXIT_FAILURE);
	}

	double alpha(atof(argv[2])), mach(atof(argv[3]));

	cout 
		<< "alpha: " << alpha << endl
		<< "mach: " << mach << endl
		<< "cl: " << get_c81_coef(data->NML, data->ml, data->NAL, 
				data->al, alpha, mach) << endl
		<< "cd: " << get_c81_coef(data->NMD, data->md, data->NAD, 
				data->ad, alpha, mach) << endl
		<< "cm: " << get_c81_coef(data->NMM, data->mm, data->NAM, 
				data->am, alpha, mach) << endl;

#else /* !USE_AERODYNAMIC_ELEMS */
	cerr << "compile with --with-aero to enable aerodynamic stuff" << endl;
#endif /* !USE_AERODYNAMIC_ELEMS */

	return 0;
}
