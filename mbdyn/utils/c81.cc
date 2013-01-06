/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "ac/getopt.h"

extern "C" {
#include "aerodc81.h"
}

#include "c81data.h"

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
   
   while (true) {
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

      c81_data_read(in, &data, 1e-6, 0);
                  
      in.close();
      
      // echo
      if (echo) {
	 c81_data_write(cout, &data);
      }
      
      argc--;
      argv++;
   }
   
   return 0;
}
