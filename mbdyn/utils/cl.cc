#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <input.h>
#include <mathp.h>
#if defined(HAVE_SSTREAM)
#include <sstream>
#elif defined(HAVE_STRSTREAM_H)
#include <strstream.h>
#else /* !HAVE_SSTREAM && !HAVE_STRSTREAM_H */
#error "no sstream and no strstream.h"
#endif /* !HAVE_SSTREAM && !HAVE_STRSTREAM_H */


int 
main(int argc, const char* const argv[])
{
	int verbose = 1;

	if (argc > 1) {
		if (strcmp(argv[1], "-?") == 0 
		    || strcmp(argv[1], "-h") == 0
		    || strcmp(argv[1], "--help") == 0) {
		    	cerr 
				<< "usage: " << argv[0] 
				<< "     reads from stdin" << endl
		   		<< "       " << argv[0] << " {-?|-h|--help}"
		   		" prints this message" << endl
	   			<< "       " << argv[0] 
				<< " <arg1> [<arg2> ...]"
				" evaluates the expressions" << endl;
	 		exit(EXIT_SUCCESS);
      		} else if (strcmp(argv[1], "-s") == 0) {
			verbose = 0;
			argv++;
			argc--;
		}
	}
	
#ifdef USE_TABLE
	Table t(31, 1);
#endif /* USE_TABLE */

	if (argc > 1) {
	 	for (int i = 1; i < argc; i++) {
#if defined(HAVE_SSTREAM)
			std::istringstream in(argv[i]);
#else /* HAVE_STRSTREAM_H */ 
	    		istrstream in(argv[i]);
#endif /* HAVE_STRSTREAM_H */
	    		InputStream In(in);
#ifdef USE_TABLE
	    		MathParser mp(In, t);
#else /* !USE_TABLE */
	    		MathParser mp(In);
#endif /* !USE_TABLE */
			if (verbose) {
	    			cout << "argv[" << i << "] = ";
			}
	    		mp.GetForever(cout, "; ");
	    		cout << endl;
	 	}
	 	exit(EXIT_SUCCESS);
      	}

	InputStream In(cin);
#ifdef USE_TABLE
      	MathParser mp(In, t);
#else /* !USE_TABLE */
      	MathParser mp(In);
#endif /* !USE_TABLE */
      	mp.GetForever(cout, "\n");
      	cout << endl;
      	exit(EXIT_SUCCESS);
}

