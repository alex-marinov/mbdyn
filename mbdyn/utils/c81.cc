#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <getopt.h>

extern "C" {
#include <aerodc81.h>
}

#include <c81data.h>

/*
 * header NML,NAL,NMD,NAD,NMM,NAM          A30,6I2 
 * ML(1),....,ML(NML)               7x,9F7.0   eventualmente su piu' righe 
 * AL(1)  CL(1,1),....,CL(1,NML)           10F7.0/(7x,9F7.0) 
 * :         :     :       : 
 * AL(NAL)CL(NAL,1),....,CL(NAL,NML)       10F7.0/(7x,9F7.0) 
 * AD(1)  CD(1,1),....,CD(1,NMD)           10F7.0/(7x,9F7.0) 
 * :         :     :       : 
 * AD(NAD)CD(NAD,1),....,CD(NAD,NMD)       10F7.0/(7x,9F7.0) 
 * AM(1)  CM(1,1),....,CL(1,NMM)           10F7.0/(7x,9F7.0) 
 * :         :     :       : 
 * AM(NAM)CM(NAM,1),....,CL(NAM,NMM)       10F7.0/(7x,9F7.0) 
 */

int
main(int argc, char* argv[]) 
{
   int opt;
   int echo = 0;
   
   while (1) {
      opt = getopt(argc, argv, "e");
      
      if (opt == EOF) {
	 break;
      }
      
      switch (opt) {
       case 'e':
	 echo = 1;
	 break;
	 
       default:
	 break;
      }
   }
   
   argc -= optind;
   argv += optind;
   
   if (argc < 1) {
      cerr << "missing filename(s); usage: c81 file [file [...]]" << endl;
      exit(EXIT_FAILURE);
   }
   
   while (argc > 0) {
      c81_data data;
      char buf[1024];      
      
      ifstream in(argv[0]);
      if (!in) {
	 cerr << "unable to open file '" << argv[0] << "'" << endl;
	 exit(EXIT_FAILURE);
      }

      read_c81_data(in, &data);
                  
      in.close();
      
      // echo
      if (echo) {
	 write_c81_data(cout, &data);
      }
      
      argc--;
      argv++;
   }
   
   return 0;
}
