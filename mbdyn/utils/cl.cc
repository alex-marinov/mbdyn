#include <mbconfig.h>

#include <input.h>
#include <mathp.h>
#include <strstream.h>

int 
main(int argn, const char* const argv[])
{
	if (argn > 1) {
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
      		}
#ifdef USE_TABLE
	 	Table t(31, 1);		
#endif
		int verbose = 1;
		int arg0 = 1;
		if (strcmp(argv[1], "-s") == 0) {
			verbose = 0;
			arg0++;
		}
	 	for (int i = arg0; i < argn; i++) {
	    		istrstream in(argv[i]);
	    		InputStream In(in);
#ifdef USE_TABLE
	    		MathParser mp(In, t);
#else
	    		MathParser mp(In);
#endif
			if (verbose) {
	    			cout << "argv[" << i << "] = ";
			}
	    		mp.GetForever(cout, "; ");
	    		cout << endl;
	 	}
	 	exit(EXIT_SUCCESS);
      	}
	
#ifdef USE_TABLE
      	Table t(31, 1);
      	InputStream In(cin);
      	MathParser mp(In, t);
#else
      	MathParser mp(In);
#endif
      	mp.GetForever(cout, "\n");
      	cout << endl;
      	exit(EXIT_SUCCESS);
}

