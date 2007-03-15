#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ActivationFunction.h>
#include <matrix.h>
#include <ann.h>
#include <getopt.h>
#include <string.h>
int TRAINING_MODE = 1;

double mean_value( matrix *MAT, int column ){

	int i;
	double mean = 0.;

	for( i=0; i<MAT->Nrow; i++ ){
		mean += MAT->mat[i][column];
	}
	return( mean/MAT->Nrow );
}

double variance( matrix *MAT, int column ){
	
	int i;
	double mean, var = 0.;

	mean = mean_value( MAT, column);
	for( i=0; i<MAT->Nrow; i++ ){
		var += (MAT->mat[i][column]-mean)*(MAT->mat[i][column]-mean);
	}

	return( var/MAT->Nrow );
}

double maximum( matrix *MAT, int column ){

	int i;
	double MAX;
	
	MAX = MAT->mat[0][column];
	for( i=0; i<MAT->Nrow; i++ ){
		if( MAT->mat[i][column] > MAX ){
			MAX = MAT->mat[i][column];
		}
	}

	return MAX;
}
	
double minimum( matrix *MAT, int column ){

	int i;
	double MIN;
	
	MIN = MAT->mat[0][column];
	for( i=0; i<MAT->Nrow; i++ ){
		if( MAT->mat[i][column] < MIN ){
			MIN = MAT->mat[i][column];
		}
	}

	return MIN;
}
struct option options[] = {
        { "usage",      0, 0, 'u'  },
        { "mode",       1, 0, 'm'  },
        { "mean",    	1, 0, 'M'  },
        { "variance",  	1, 0, 'V'  },
        { "min",    	1, 0, 'N'  },
        { "max",  	1, 0, 'X'  },
        { "file",    	1, 0, 'f'  } 
};
void print_usage( void ){

        fprintf( stdout, "\nUSAGE OPTIONS:\n"
                "  -u, --usage\n"
                "       print usage\n"
                "  -f, --file\n"
                "       input file name\n"
                "  -m, --mode\n"
                "       scaling mode: 1 mean-variance  	(default)\n"
                "                     2 min-max\n" 
                "  -M, --mean\n"
                "       mean value (defaul 0.)\n"
                "  -V, --variance"
                "       variance value (defaul 1.)\n"
                "  -N, --min\n"
                "       minimum value (defaul -1.)\n"
                "  -X, --max\n"
                "       maximum value (defaul 1.)\n"
         );
        exit(0);

}


int MODE = 1;
float MEAN = 0.;
float VARIANCE = 1.;
float MIN = -1.;
float MAX = 1.;
static char *filename = "data.dat";

int main( int argc, char **argv ){

	int opt;
        extern char *optarg;
	matrix MAT, SF;
	int N_sample, i;
	double mean, var, min, max;

        /* 0. Training options */
        do {
                opt = getopt_long( argc, argv, "uf:m:M:V:N:X:", options, NULL  );
                switch( opt ) {
                case 'u':       print_usage();
                                break;   
		case 'f':       filename = strdup( optarg );
                                break;
                case 'm':       MODE = atoi( optarg );
                                break;
                case 'M':       MEAN = atof( optarg );
                                break;
                case 'V':       VARIANCE = atof( optarg );
                                break;
                case 'N':       MIN = atof( optarg );
                                break;
                case 'X':       MAX = atof( optarg );
                                break;
                default:        break;
                }
        }while( opt >= 0  );
	
        N_sample = 0;
        if( ANN_DataRead( &MAT, &N_sample, filename ) ){
                fprintf( stderr, "Error in Input data acquisition\n");
                return 1; 
        }	
	matrix_init( &SF, MAT.Ncolumn, 2 );
	if( MODE == 1 ){
		for( i=0; i<MAT.Ncolumn; i++ ){
			mean = mean_value( &MAT, i );
			var = variance( &MAT, i );
			SF.mat[i][0] = sqrt( VARIANCE/var );
			SF.mat[i][1] = MEAN - SF.mat[i][0]*mean;
		}
	}
	if( MODE == 2 ){
		for( i=0; i<MAT.Ncolumn; i++ ){
			min = minimum( &MAT, i );
			max = maximum( &MAT, i );
			printf("AAAA   %lf   %lf\n",min,max);
			SF.mat[i][0] = (MAX-MIN)/(max-min);
			SF.mat[i][1] = (max*MIN-MAX*min)/(max-min);
		}
	}
	matrix_write( &SF, stdout, W_M_TEXT );
	
	return 0;

}
