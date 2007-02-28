/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2007
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
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
/*
 * Copyright (C) 2007
 *
 * Mattia Mattaboni	<mattaboni@aero.polimi.it>
 */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define FILE_NAME_MAX_LENGTH 100

/* input files' names*/
char fname[FILE_NAME_MAX_LENGTH];
char PathName[FILE_NAME_MAX_LENGTH];
char fnamePAR[FILE_NAME_MAX_LENGTH];
char fnameANN[FILE_NAME_MAX_LENGTH];    /*ANN initilal definition file*/
char fnameTR[FILE_NAME_MAX_LENGTH];     /*ANN trained definition file*/
char fnameINPUT[FILE_NAME_MAX_LENGTH];  /*ANN training input file*/
char fnameOUTPUT[FILE_NAME_MAX_LENGTH]; /*ANN training output file*/

int MAXITER;
double TOLL;
int MODE;
int SAVESTEP;
int PRINTSTEP;
int GR_CHECK;
double dw;


/* Artificial Multi-Layer Feedforward Neural Network */
typedef struct ANN{
	int N_input;     /* network input number*/
	int N_output;    /* network (visible)output number*/
	int N_layer;     /* network layer number*/
	int *N_neuron;   /* neuron number*/
	double alpha;    /* activation function coefficient*/
	double beta;     /* activation function coefficient*/
	double eta;      /* learning rate*/
	double rho;      /* momentum term*/
	int r;           /* number of time delay*/
	double ***W;     /* network synaptic weights*/
	double ***dW;    /* error gradient*/
	double **v;      /* neurons' internal activity*/
	double *****dy;  /* output network gradient*/
	double *yD;
	}ANN;

/*FUNCTIONS' PROTOTYPES:    */

void ANN_Strcat(char *, char *, char *);
int ANN_FileName();
int ANN_Parameters();
int ANN_Initialize( ANN * );
int ANN_TrainingData( double *** , double *** , double *** , int *);
void ANN_Free( ANN * , double *** , double *** , double *** , int );
int ANN_sim( ANN * , double * , double * );
void ANN_print( ANN * );
int ANN_GradientCheck(ANN * , double ** , double ** ,int);
int ANN_GradientCheck2( ANN * , double ** , double ** , int );
int ANN_dW( ANN * , double * );
int ANN_dXdW( ANN * , int , int , int , double *);
int ANN_TrainingEpoch( ANN * , double ** , double ** , double ** , int ,int ,double ***);
void ANN_WeightUpdate( ANN * , double *** );
double ANN_InternalFunction( double , ANN * );
double ANN_InternalFunctionDer( double , ANN * );
int ANN_save( ANN * );


/* MAIN FUNCTION:    */
int main( int argc , char **argv ){
	ANN net;
	int i,j,k,p,q;
	double **INPUT, **NN_OUTPUT, **DES_OUTPUT;
	int N_sample;
	double err1,err2;
	int Niter;
	int CNT = 0;
	double ***DWnew,***DWold,***W1,***W2;

		
	/* input files' names acquisition*/
	if( argc == 2 ){
		strcpy(fname,argv[1]);
	}
	if( argc == 3 ){
		ANN_Strcat(argv[2],argv[1],fname);
	}
	if( !ANN_FileName() ){
		printf("\nError: ANN_FileName@main\n");
		return(0);
	}
	if( !ANN_Parameters() ){
		printf("\nError: ANN_Parameters@main\n");
		return(0);
	}
	/* Artificial Neural Network inizialization*/
	printf("\nLOADING DATA...");
	if( !ANN_Initialize( &net ) ){
		printf("\nError: ANN_initilize@main\n");
		return(0);
	}
	/* Training data acquisition*/
	N_sample = 0;
	if( !ANN_TrainingData( &INPUT,&DES_OUTPUT,&NN_OUTPUT,&N_sample ) ){
		printf("\nError: ANN_TrainingData@main\n");
		return(0);
	}

	/*ANN_GradientCheck(&net,INPUT,DES_OUTPUT,100);*/
	/*ANN_GradientCheck2(&net,INPUT,DES_OUTPUT,N_sample);*/
	ANN_print(&net);
	if( !( DWnew = (double ***)malloc( (net.N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(DWnew)\n");
		return(0);
	}
	if( !( DWold = (double ***)malloc( (net.N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(DWold)\n");
		return(0);
	}
	if( !( W1 = (double ***)malloc( (net.N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(DWold)\n");
		return(0);
	}
	if( !( W2 = (double ***)malloc( (net.N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(DWold)\n");
		return(0);
	}
	for( i=0;i<(net.N_layer+1);i++ ){
		if( !( DWnew[i] = (double **)malloc( net.N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(DWnew)\n");
			return(0);
		}
		if( !( DWold[i] = (double **)malloc( net.N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(DWold)\n");
			return(0);
		}
		if( !( W1[i] = (double **)malloc( net.N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(DWold)\n");
			return(0);
		}
		if( !( W2[i] = (double **)malloc( net.N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(DWold)\n");
			return(0);
		}
		for( j=0;j<net.N_neuron[i];j++ ){
			if( !( DWnew[i][j] = (double *)malloc( net.N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(DWnew)\n");
				return(0);
			}
			if( !( DWold[i][j] = (double *)malloc( net.N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(DWold)\n");
				return(0);
			}
			if( !( W1[i][j] = (double *)malloc( net.N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(DWnew)\n");
				return(0);
			}
			if( !( W2[i][j] = (double *)malloc( net.N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(DWnew)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net.N_layer+1;i++ ){
		for( j=0;j<net.N_neuron[i];j++ ){
			for( k=0;k<net.N_neuron[i+1];k++ ){
				DWold[i][j][k] = 0.;
				DWnew[i][j][k] = 0.;
				W1[i][j][k] = 0.;
				W2[i][j][k] = 0.;
			}
		}
	}

	Niter = 0.;
	err2 = 1000000000;
	do{
		err1 = err2;
		if( GR_CHECK == 1 ){
			if( !ANN_GradientCheck2( &net,INPUT,DES_OUTPUT,N_sample ) ){
				printf("\nError: ANN_sim@ANN_GradientCheck\n");
				return(0);
			}
			printf("\n");
		}
		for( i=0;i<net.N_layer+1;i++ ){
			for( j=0;j<net.N_neuron[i];j++ ){
				for( k=0;k<net.N_neuron[i+1];k++ ){
					net.dW[i][j][k] = 0.;
				}
			}
		}
		for( p=0;p<net.r;p++ ){
			for( q=0;q<net.N_neuron[net.N_layer+1];q++ ){
				for ( i=0; i<net.N_layer+1;i++ ){
					for( j=0;j<net.N_neuron[i];j++ ){
						for( k=0;k<net.N_neuron[i+1];k++ ){
							net.dy[p][q][i][j][k] = 0.;
						}
					}
				}
			}
		}
		for(  i=0;i<(net.r*net.N_neuron[net.N_layer+1]);i++){
			net.yD[i] = 0.;
		}
		Niter++;
		for( i=0;i<net.N_layer+1;i++ ){
			for( j=0;j<net.N_neuron[i];j++ ){
				for( k=0;k<net.N_neuron[i+1];k++ ){
					W2[i][j][k] = W1[i][j][k];
					W1[i][j][k] = net.W[i][j][k];
				}
			}
		}
		if( !ANN_TrainingEpoch( &net , INPUT , DES_OUTPUT , NN_OUTPUT, N_sample , MODE , DWnew ) ){
			printf("\nError: ANN_TrainingEpoch@main\n");
			return(0);
		}
		err2 = 0.;
		for( j=0;j<N_sample;j++ ){
			for( i=0;i<net.N_output;i++ ){
				err2 += 0.5*(DES_OUTPUT[j][i]-NN_OUTPUT[j][i])*(DES_OUTPUT[j][i]-NN_OUTPUT[j][i]);
			}
		}
		if (MODE == 1){
			CNT++;
			while ( err2 >= err1 ){
				printf("\n......    %lf\n",err2-err1);
				CNT = 0;
				net.eta = 0.5*net.eta;
				printf("Network's learning rate decreasing (eta = %lf)\n",net.eta);
				for( i=0;i<net.N_layer+1;i++ ){
					for( j=0;j<net.N_neuron[i];j++ ){
						for( k=0;k<net.N_neuron[i+1];k++ ){
							//net.W[i][j][k] = net.W[i][j][k]-DWold[i][j][k]-DWnew[i][j][k];
							net.W[i][j][k] = W2[i][j][k];
							DWold[i][j][k] = 0.5*DWold[i][j][k];
						}
					}
				}
				for( i=0;i<net.N_layer+1;i++ ){
					for( j=0;j<net.N_neuron[i];j++ ){
						for( k=0;k<net.N_neuron[i+1];k++ ){
							net.W[i][j][k] += DWold[i][j][k];						
							W1[i][j][k] = net.W[i][j][k];
						}
					}
				}
				for( i=0;i<net.N_layer+1;i++ ){
					for( j=0;j<net.N_neuron[i];j++ ){
						for( k=0;k<net.N_neuron[i+1];k++ ){
							net.dW[i][j][k] = 0.;
						}
					}
				}
				for( p=0;p<net.r;p++ ){
					for( q=0;q<net.N_neuron[net.N_layer+1];q++ ){
						for ( i=0; i<net.N_layer+1;i++ ){
							for( j=0;j<net.N_neuron[i];j++ ){
								for( k=0;k<net.N_neuron[i+1];k++ ){
									net.dy[p][q][i][j][k] = 0.;
								}
							}
						}
					}
				}
				for(  i=0;i<(net.r*net.N_neuron[net.N_layer+1]);i++){
					net.yD[i] = 0.;
				}
				if( !ANN_TrainingEpoch( &net , INPUT , DES_OUTPUT , NN_OUTPUT, N_sample , MODE , DWnew) ){
					printf("\nError: ANN_TrainingEpoch@main\n");
					return(0);
				}
				err2 = 0.;
				for( j=0;j<N_sample;j++ ){
					for( i=0;i<net.N_output;i++ ){
						err2 += 0.5*(DES_OUTPUT[j][i]-NN_OUTPUT[j][i])*(DES_OUTPUT[j][i]-NN_OUTPUT[j][i]);
					}
				}
			}	

			for( i=0;i<net.N_layer+1;i++ ){
				for( j=0;j<net.N_neuron[i];j++ ){
					for( k=0;k<net.N_neuron[i+1];k++ ){
						DWold[i][j][k] = DWnew[i][j][k];
					}
				}
			}
			if (CNT == 10){
				net.eta=1.5*net.eta;
				printf("Network's learning rate increasing (eta = %lf)\n",net.eta);
				CNT = 0;
			}
		}
			
		if( !(Niter%PRINTSTEP) ){
			printf("\nTRAINING:    iter:%d       ",Niter);
			printf("err:%1.12lf",err2);
		}

		if( !(Niter%SAVESTEP) ){
			printf("\nSAVING DATA...");
			if( !ANN_save( &net ) ){
				printf("\nError: ANN_save@main\n");
				return(0);
			}
		}

	}while( (err2>TOLL) && (Niter<MAXITER) );

	printf("\nSAVING DATA...");
	if( !ANN_save( &net ) ){
		printf("\nError: ANN_save@main\n");
		return(0);
	}

	/* dynamic memory free*/
	ANN_Free( &net,&INPUT,&DES_OUTPUT,&NN_OUTPUT,N_sample );
	return(1);
}



/* FUNCTION:  */

void ANN_Strcat(char *s1,char *s2,char *s3 ){

	int i;
	int len1,len2;
	char s[FILE_NAME_MAX_LENGTH];

	
	len1 = strlen(s1)-1;
	len2 = strlen(s2)-1;
	strcpy(s,s1);
	
	for( i=0;i<len2;i++ ){
		s[i+len1] = s2[i];
	}
	s[len1+len2] = (char)NULL;
	
	strcpy(s3,s);
}

int ANN_FileName(){

	FILE * fp;
	int len;

	if ( !(fp = fopen( fname , "r" ) ) ){
		printf("\nError: fopen@ANN_FileName\n");
		return(0);
	}


	fgets(PathName,FILE_NAME_MAX_LENGTH,fp);
	len = strlen(PathName);
	PathName[len-1] = (char)NULL;
	ANN_Strcat(PathName,fgets(fnameANN,FILE_NAME_MAX_LENGTH,fp),fnameANN);
	len = strlen(fnameANN);
	fnameANN[len-1] = (char)NULL;
	ANN_Strcat(PathName,fgets(fnameINPUT,FILE_NAME_MAX_LENGTH,fp),fnameINPUT);
	len = strlen(fnameINPUT);
	fnameINPUT[len-1] = (char)NULL;
	ANN_Strcat(PathName,fgets(fnameOUTPUT,FILE_NAME_MAX_LENGTH,fp),fnameOUTPUT);
	len = strlen(fnameOUTPUT);
	fnameOUTPUT[len-1] = (char)NULL;
	ANN_Strcat(PathName,fgets(fnameTR,FILE_NAME_MAX_LENGTH,fp),fnameTR);
	len = strlen(fnameTR);
	fnameTR[len-1] = (char)NULL;
	ANN_Strcat(PathName,fgets(fnamePAR,FILE_NAME_MAX_LENGTH,fp),fnamePAR);
	len = strlen(fnamePAR);
	fnamePAR[len-1] = (char)NULL;

	fclose(fp);
	return(1);
}

int ANN_Parameters(){

	FILE *fp;
	printf("%s",fnamePAR);
	if ( !(fp = fopen( fnamePAR , "r" ) ) ){
		printf("\nError: fopen@ANN_Parameters\n");
		return(0);
	}
	fscanf(fp,"%d",&MAXITER);
	fscanf(fp,"%lf",&TOLL);
	fscanf(fp,"%d",&MODE);
	fscanf(fp,"%d",&SAVESTEP);
	fscanf(fp,"%d",&PRINTSTEP);
	fscanf(fp,"%lf",&dw);
	fscanf(fp,"%d",&GR_CHECK);

	fclose(fp);
	return(1);

}


int ANN_Initialize( ANN *net ){

	int i,j,k,p,q;
	FILE *fp;

	if ( !(fp = fopen( fnameANN , "r" ) ) ){
		printf("\nError: fopen@ANN_initialize\n");
		return(0);
	}


	fscanf(fp,"%d",&(net->N_input));
	fscanf(fp,"%d",&(net->N_output));
	fscanf(fp,"%d",&(net->N_layer));
	fscanf(fp,"%d",&(net->r));


	if( !(net->N_neuron = (int *)malloc( (net->N_layer+2) * sizeof(int) ) ) ){
			printf("\nError: malloc@ANN_initialize(N_neuron)\n");
			return(0);
	}
	net->N_neuron[0] = 0;
	for( i=0;i<net->N_layer+1;i++ ){
		fscanf(fp,"%d",&(net->N_neuron[i+1]));
	}
	net->N_neuron[0] = net->N_input + ( net->r * net->N_neuron[net->N_layer+1] );


	fscanf(fp,"%lf",&(net->alpha));
	fscanf(fp,"%lf",&(net->beta));
	fscanf(fp,"%lf",&(net->eta));
	fscanf(fp,"%lf",&(net->rho));


	if( !( net->W = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(W)\n");
		return(0);
	}
	for( i=0;i<(net->N_layer+1);i++ ){
		if( !( net->W[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(W)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( net->W[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(W)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				fscanf(fp,"%lf",&(net->W[i][j][k]));
			}
		}
	}


	if( !( net->dW = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_initialize(dW)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+1;i++ ){
		if( !( net->dW[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_initialize(dW)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( net->dW[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_initialize(dW)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				net->dW[i][j][k] = 0.;
			}
		}
	}


	if( !( net->v = (double **)malloc( (net->N_layer + 2) * sizeof(double *) ) ) ){
		printf("\nError: malloc@ANN_initialize(v)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+2;i++ ){
		if( !( net->v[i] = (double *)malloc( net->N_neuron[i] * sizeof(double ) ) ) ){
			printf("\nError: malloc@ANN_initialize(v)\n");
			return(0);
		}
	}
	for( i=0;i<net->N_layer+2;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			net->v[i][j] = 0.;
		}
	}


	if( net->r != 0 ){
		if( !( net->dy = (double *****)malloc( net->r * sizeof(double ****) ) ) ){
			printf("\nError: malloc@ANN_initialize(dy)\n");
			return(0);
		}
		for( p=0;p<net->r;p++ ){
			if( !( net->dy[p] = (double ****)malloc( net->N_neuron[net->N_layer+1] * sizeof(double ***) ) ) ){
				printf("\nError: malloc@ANN_initialize(dy)\n");
				return(0);

			}
			for( i=0;i<net->N_neuron[net->N_layer+1];i++ ){
				if( !( net->dy[p][i] = (double ***)malloc( (net->N_layer+1) * sizeof(double **) ) ) ){
					printf("\nError: malloc@ANN_initialize(dy)\n");
					return(0);
				}
				for( j=0;j<net->N_layer+1;j++ ){
					if( !( net->dy[p][i][j] = (double **)malloc( net->N_neuron[j] * sizeof(double *) ) ) ){
						printf("\nError: malloc@ANN_initialize(dy)\n");
						return(0);
					}
					for( k=0;k<net->N_neuron[j];k++ ){
						if( !(net->dy[p][i][j][k] = (double *)malloc( net->N_neuron[j+1] * sizeof(double ) ) ) ){
							printf("\nError: malloc@ANN_initialize(dy)\n");
							return(0);
						}



					}
				}
			}
		}
		for( p=0;p<net->r;p++ ){
			for( q=0;q<net->N_neuron[net->N_layer+1];q++ ){
				for ( i=0; i<net->N_layer+1;i++ ){
					for( j=0;j<net->N_neuron[i];j++ ){
						for( k=0;k<net->N_neuron[i+1];k++ ){
							net->dy[p][q][i][j][k] = 0.;
						}
					}
				}
			}
		}

		if( !( net->yD = (double *)malloc( (net->r * net->N_neuron[net->N_layer+1]) * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_initialize(yD)\n");
			return(0);
		}
		for(  i=0;i<(net->r*net->N_neuron[net->N_layer+1]);i++){
			net->yD[i] = 0.;
		}
	}

	fclose(fp);

	return(1);
}


int ANN_TrainingData( double ***INPUT , double ***DES_OUTPUT , double ***NN_OUTPUT , int * N_s){

	int N_input,N_output,N_sample;
	int i,j;
	FILE *fp;

	if ( !(fp = fopen( fnameINPUT , "r" ) ) ){
		printf("\nError: fopen@ANN_TrainingData(INPUT)\n");
		return(0);
	}

	fscanf(fp,"%d",&N_sample);
	fscanf(fp,"%d",&N_input);

	if( !( *INPUT = (double **)malloc( N_sample * sizeof(double *) ) ) ){
		printf("\nError: malloc@ANN_TrainingData(INPUT)\n");
		return(0);
	}
	for( i=0;i<N_sample;i++ ){
		if( !( (*INPUT)[i] = (double *)malloc( N_input * sizeof(double ) ) ) ){
			printf("\nError: malloc@ANN_TrainingData(INPUT)\n");
			return(0);
		}
	}
	for( i=0;i<N_sample;i++ ){
		for( j=0;j<N_input;j++ ){
			fscanf(fp,"%lf",&(*INPUT)[i][j]);
		}
	}

	fclose(fp);


	if ( !(fp = fopen( fnameOUTPUT , "r" )) ){
		printf("\nError: fopen@ANN_TrainingData(OUTPUT)\n");
		return(0);
	}

	fscanf(fp,"%d",&N_sample);
	fscanf(fp,"%d",&N_output);

	if( !( *DES_OUTPUT = (double **)malloc( N_sample * sizeof(double *) ) ) ){
		printf("\nError: malloc@ANN_TrainingData(DES_OUTPUT)\n");
		return(0);
	}
	for( i=0;i<N_sample;i++ ){
		if( !( (*DES_OUTPUT)[i] = (double *)malloc( N_output * sizeof(double ) ) ) ){
			printf("\nError: malloc@ANN_TrainingData(DES_OUTPUT)\n");
			return(0);
		}
	}
	for( i=0;i<N_sample;i++ ){
		for( j=0;j<N_output;j++ ){
			fscanf(fp,"%lf",&(*DES_OUTPUT)[i][j]);
		}
	}

	fclose(fp);

	if( !( *NN_OUTPUT = (double **)malloc( N_sample * sizeof(double *) ) ) ){
		printf("\nError: malloc@ANN_TrainingData(NN_OUTPUT)\n");
		return(0);
	}
	for( i=0;i<N_sample;i++ ){
		if( !( (*NN_OUTPUT)[i] = (double *)malloc( N_output * sizeof(double ) ) ) ){
			printf("\nError: malloc@ANN_TrainingData(NN_OUTPUT)\n");
			return(0);
		}
	}
	for( i=0;i<N_sample;i++ ){
		for( j=0;j<N_output;j++ ){
			(*NN_OUTPUT)[i][j] = 0.;
		}
	}
	*N_s = N_sample;

	return(1);
}


void ANN_Free( ANN * net , double *** INPUT , double *** DES_OUTPUT , double *** NN_OUTPUT , int N_sample ){

	int i,j,k,p;

	for( i=0;i<N_sample;i++ ){
		free((*INPUT)[i]);
		free((*DES_OUTPUT)[i]);
		free((*NN_OUTPUT)[i]);
	}
	free(*INPUT);
	free(*DES_OUTPUT);
	free(*NN_OUTPUT);

	if( net->r != 0 ){
		for( i=0;i<net->r;i++ ){
			for( j=0;j<net->N_neuron[net->N_layer+1];j++ ){
				for( k=0;k<net->N_layer+1;k++ ){
					for( p=0;p<net->N_neuron[k];p++ ){
						free(net->dy[i][j][k][p]);
					}
					free(net->dy[i][j][k]);
				}
				free(net->dy[i][j]);
			}
			free(net->dy[i]);
		}
		free(net->dy);
		free(net->yD);
	}

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(net->dW[i][j]);
		}
		free(net->dW[i]);
	}
	free(net->dW);

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(net->W[i][j]);
		}
		free(net->W[i]);
	}
	free(net->W);

	for( i=0;i<net->N_layer+2;i++ ){
		free(net->v[i]);
	}
	free(net->v);

	free(net->N_neuron);

}


void ANN_print( ANN *net ){

	int i,j,k;


	printf("\nInput number: %d ",net->N_input);
	printf("\nOutput number: %d ",net->N_output);
	printf("\nHidden layers number: %d ",net->N_layer);

	for( i=0;i<net->N_layer;i++ ){
		printf("\nNeurons number (layer number %d) : %d",i+1,net->N_neuron[i+1]);
	}

	printf("\nLearning rate: %lf ",net->eta);
	printf("\nMomentum term: %lf ",net->rho);

	
	printf("\nActivation function, alpha: %lf ",net->alpha);
	printf("\nActivation fuction, beta: %lf ",net->beta);


	for( i=0;i<net->N_layer+1;i++ ){
		if( i!=net->N_layer )
			printf("\nSynaptic Weight Matrix ( layer number %d ):\n",i+1);
		else
			printf("\nSynaptic Weight Matrix ( visible layer ):\n");
		for( j=0;j<net->N_neuron[i];j++ ){
			printf("\n");
			for( k=0;k<net->N_neuron[i+1];k++ ){
				printf("%lf ",net->W[i][j][k]);
			}
		}
	}
}


int ANN_sim( ANN *net , double *input , double *output ){

	int i,j,k;
	double *IA,*y,*y_p;

	if( !(y_p = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(y_p)\n");
		return(0);
	}
	if( !(IA = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(IA)\n");
		return(0);
	}
	if( !(y = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_sim(y)\n");
		return(0);
	}

	for( i=0;i<net->N_neuron[0];i++ ){
		y_p[i] = input[i];
	}

	for( i=0;i<net->N_neuron[0];i++ )
		net->v[0][i] = input[i];

	for( i=0;i<net->N_layer+1;i++ ){
		if( !( IA = (double *)realloc( IA , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(IA)\n");
			return(0);
		}
		if( !( y = (double *)realloc( y , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(y)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i+1];j++ ){
			IA[j] = 0.;
			for( k=0;k<net->N_neuron[i];k++ ){
				IA[j] += net->W[i][k][j]*y_p[k];
			}

			y[j] = ANN_InternalFunction(IA[j],net);
		}

		for( j=0;j<net->N_neuron[i+1];j++ )
			net->v[i+1][j] = IA[j];

		if( !( y_p = (double *)realloc( y_p , net->N_neuron[i+1]*sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_sim(y_p)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i+1];j++ ){
			y_p[j] = y[j];
		}
	}

	for( i=0;i<net->N_output;i++ ){
		output[i] = y_p[i];
	}

	for( i=0;i<(net->r-1)*(net->N_neuron[net->N_layer+1]);i++ ){
		net->yD[((net->r)*(net->N_neuron[net->N_layer+1]))-1-i] = net->yD[((net->r-1)*(net->N_neuron[net->N_layer+1]))-1-i];
	}
	if( net->r != 0 ){
		for( i=0;i<net->N_neuron[net->N_layer+1];i++ ){
			net->yD[i] = y_p[i];
		}
	}

	free(IA);
	free(y);
	free(y_p);

	return(1);

}



int ANN_GradientCheck(ANN * net , double ** INPUT , double ** DES_OUTPUT , int t ){

	int i,j,k,p,a,I,J,K;
	double E1,E2,MAXerror,error;
	double *input,*output,*e;
	double ***dEdW_FD;

	if( !( input = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck(input)\n");
		return(0);
	}
	if( !( output = (double *)malloc( net->N_output * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck(output)\n");
		return(0);
	}
	if( !( e = (double *)malloc( net->N_output * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck(e)\n");
		return(0);
	}

	if( !( dEdW_FD = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck(dEdw_FD)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+1;i++ ){
		if( !( dEdW_FD[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_GradientCheck(dEdw_FD)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( dEdW_FD[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_GradientCheck(dEdw_FD)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				dEdW_FD[i][j][k] = 0.;
			}
		}
	}

	a = 0;
	for( i=0;i<net->N_input;i++ ){
		input[a] = INPUT[t][i];
		a++;
	}
	for( i=0;i<(net->r * net->N_neuron[net->N_layer+1]);i++ ){
		input[a] = net->yD[i];
		a++;
	}

	for( i=0;i<net->N_output;i++ ){
		output[i] = 0.;
		e[i] = 0.;
	}

	if( !ANN_sim( net,input,output ) ){
		printf("\nError: ANN_sim@ANN_GradientCheck\n");
		return(0);
	}


	for( i=0;i<net->N_output;i++ ){
		e[i] = DES_OUTPUT[t][i]-output[i];
	}
	if( !ANN_dW( net , e ) ){
		printf("\nError: ANN_dW@ANN_GradientCheck\n");
		return(0);
	}

	E1 = 0.;
	for( i=0;i<net->N_output;i++ )
		E1 += 0.5*( DES_OUTPUT[t][i]-output[i] )*( DES_OUTPUT[t][i]-output[i] );

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){

				for( p=0;p<net->r*net->N_neuron[net->N_layer+1];p++ ){
					net->yD[p] = 0.;
				}

				net->W[i][j][k] += dw;
				for( p=0;p<net->N_output;p++ )
					output[p] = 0.;
				if( !ANN_sim( net,input,output ) ){
					printf("\nError: ANN_sim@ANN_GradientCheck\n");
					return(0);
				}
				net->W[i][j][k] -= dw;
				E2 = 0.;
				for( p=0;p<net->N_output;p++ )
					E2 += 0.5*( DES_OUTPUT[t][p]-output[p] )*( DES_OUTPUT[t][p]-output[p] );


				dEdW_FD[i][j][k] = (E2-E1)/dw;
			}
		}
	}

	/*printf("Finite Difference:       Analytic:      Perc error:   ");
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				if( fabs(dEdW_FD[i][j][k]) > 0.0000001){
					printf("\n%lf   %lf   %lf",dEdW_FD[i][j][k],-net->dW[i][j][k]/net->eta,(100*(dEdW_FD[i][j][k]+(net->dW[i][j][k])/net->eta))/dEdW_FD[i][j][k]);
				}
				else{
					printf("\n%lf   %lf",dEdW_FD[i][j][k],-net->dW[i][j][k]);
				}

			}
		}
	}*/

	MAXerror = 0.;
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				if( fabs(dEdW_FD[i][j][k]) > 0.0000001){
					error = (dEdW_FD[i][j][k]+net->dW[i][j][k]/net->eta)/dEdW_FD[i][j][k];
					if( fabs(error) > MAXerror ){
						MAXerror = fabs(error);
						I = i;
						J = j;
						K = k;
					}
				}
			}
		}
	}
	printf("\n\t\tMAX_ERROR: %1.2lf ( FD: %1.2lf,AN: %1.2lf ) ",(100*(dEdW_FD[I][J][K]+net->dW[I][J][K]/net->eta))/dEdW_FD[I][J][K],dEdW_FD[I][J][K],-(net->dW[I][J][K]/net->eta));

	free(e);
	free(input);
	free(output);
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(dEdW_FD[i][j]);
		}
		free(dEdW_FD[i]);
	}
	free(dEdW_FD);
	//getchar();


	return(1);
}


int ANN_GradientCheck2( ANN * net , double ** INPUT , double ** DES_OUTPUT , int N_sample ){

	int i,j,k,p,t,a,I,J,K;
	double EE,Et,E0,MAXerror,error;
	double *input,*output,*e;
	double ***dEdW_FD,***DW;

	if( !( input = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck2(input)\n");
		return(0);
	}
	if( !( output = (double *)malloc( net->N_output * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_GradientCheck2(ioutput)\n");
			return(0);
	}
	if( !( e = (double *)malloc( net->N_output * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_GradientCheck2(e)\n");
			return(0);
	}
	if( !( dEdW_FD = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck2(dEdW_FD)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+1;i++ ){
		if( !( dEdW_FD[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_GradientCheck2(dEdW_FD)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( dEdW_FD[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_GradientCheck2(dEdW_FD)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				dEdW_FD[i][j][k] = 0.;
			}
		}
	}
	for( p=0;p<net->r*net->N_neuron[net->N_layer+1];p++ ){
		net->yD[p] = 0.;
	}
	E0 = 0.;
	for( t=0;t<N_sample;t++ ){
		a = 0;
		for( p=0;p<net->N_input;p++ ){
			input[a] = INPUT[t][a];
			a++;
		}
		for( p=0;p<(net->r)*(net->N_neuron[net->N_layer+1]);p++ ){
			input[a] = net->yD[p];
			a++;
		}

		for( p=0;p<net->N_output;p++ )
			output[p] = 0.;

		if( !ANN_sim( net,input,output ) ){
			printf("\nError: ANN_sim@ANN_GradientCheck2\n");
			return(0);
		}
		Et = 0.;
		for( p=0;p<net->N_output;p++ ){
			Et += (DES_OUTPUT[t][p]-output[p])*(DES_OUTPUT[t][p]-output[p]);
		}
		E0 += 0.5*Et;
	}

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				for( p=0;p<net->r*net->N_neuron[net->N_layer+1];p++ ){
					net->yD[p] = 0.;
				}
				net->W[i][j][k] += dw;
				EE = 0.;
				for( t=0;t<N_sample;t++ ){
					a = 0;
					for( p=0;p<net->N_input;p++ ){
						input[a] = INPUT[t][a];
						a++;
					}
					for( p=0;p<(net->r)*(net->N_neuron[net->N_layer+1]);p++ ){
						input[a] = net->yD[p];
						a++;
					}
					for( p=0;p<net->N_output;p++ )
						output[p] = 0.;
					if( !ANN_sim( net,input,output ) ){
						printf("\nError: ANN_sim@ANN_GradientCheck2\n");
						return(0);
					}
					Et = 0.;
					for( p=0;p<net->N_output;p++ ){
						Et += (DES_OUTPUT[t][p]-output[p])*(DES_OUTPUT[t][p]-output[p]);
					}
					EE += 0.5*Et;
				}
				net->W[i][j][k] -= dw;
				//dEdW_FD[i][j][k] = (EE-E0)/dw;
				dEdW_FD[i][j][k] = EE;
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				for( p=0;p<net->r*net->N_neuron[net->N_layer+1];p++ ){
					net->yD[p] = 0.;
				}
				net->W[i][j][k] -= dw;
				EE = 0.;
				for( t=0;t<N_sample;t++ ){
					a = 0;
					for( p=0;p<net->N_input;p++ ){
						input[a] = INPUT[t][a];
						a++;
					}
					for( p=0;p<(net->r)*(net->N_neuron[net->N_layer+1]);p++ ){
						input[a] = net->yD[p];
						a++;
					}
					for( p=0;p<net->N_output;p++ )
						output[p] = 0.;
					if( !ANN_sim( net,input,output ) ){
						printf("\nError: ANN_sim@ANN_GradientCheck2\n");
						return(0);
					}
					Et = 0.;
					for( p=0;p<net->N_output;p++ ){
						Et += (DES_OUTPUT[t][p]-output[p])*(DES_OUTPUT[t][p]-output[p]);
					}
					EE += 0.5*Et;
				}
				net->W[i][j][k] += dw;
				//dEdW_FD[i][j][k] = (EE-E0)/dw;
				dEdW_FD[i][j][k] = (dEdW_FD[i][j][k] - EE)/(2*dw);
			}
		}
	}

	if( !( DW = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_GradientCheck2(DW)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+1;i++ ){
		if( !( DW[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_GradientCheck2(DW)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( DW[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_GradientCheck2(DW)\n");
				return(0);
			}
		}
	}

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				DW[i][j][k] = 0.;
			}
		}
	}
	for( p=0;p<net->r*net->N_neuron[net->N_layer+1];p++ ){
		net->yD[p] = 0.;
	}
	for( t=0;t<N_sample;t++ ){
		a = 0;
		for( i=0;i<net->N_input;i++ ){
			input[a] = INPUT[t][i];
			a++;
		}
		for( i=0;i<(net->r*net->N_neuron[net->N_layer+1]);i++ ){
			input[a] = net->yD[i];
			a++;
		}

		for( i=0;i<net->N_output;i++ )
			output[i]=0.;
		if( !ANN_sim( net,input,output ) ){
			printf("\nError: ANN_sim@ANN_GradientCheck2\n");
			return(0);
		}
		for( i=0;i<net->N_output;i++ ){
			e[i] = DES_OUTPUT[t][i] - output[i];
		}
		if( !ANN_dW( net , e ) ){
			printf("\nError: ANN_dW@ANN_GradientCheck2\n");
			return(0);
		}
		for( i=0;i<net->N_layer+1;i++ ){
			for( j=0;j<net->N_neuron[i];j++ ){
				for( k=0;k<net->N_neuron[i+1];k++ ){
					DW[i][j][k] += net->dW[i][j][k];
				}
			}
		}


	}


	printf("Finite Difference:       Analytic:      Perc error:   ");
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				if( dEdW_FD[i][j][k] != 0 ){
					printf("\n%lf   %lf   %lf",dEdW_FD[i][j][k],-DW[i][j][k]/net->eta,(100*(dEdW_FD[i][j][k]+DW[i][j][k]/net->eta))/dEdW_FD[i][j][k]);
				}
				else{
					printf("\n%lf   %lf",dEdW_FD[i][j][k],-DW[i][j][k]/net->eta);
				}
			}
		}
	}

	MAXerror = 0.;
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				if( dEdW_FD[i][j][k] != 0 ){
					error = (dEdW_FD[i][j][k]+DW[i][j][k]/net->eta)/dEdW_FD[i][j][k];
					if( fabs(error) > MAXerror ){
						MAXerror = fabs(error);
						I = i;
						J = j;
						K= k;
					}
				}
			}
		}
	}
	printf("\n\t\tMAX_ERROR: %1.2lf ( FD: %1.2lf,AN: %1.2lf ) ",(100*(dEdW_FD[I][J][K]+DW[I][J][K]/net->eta))/dEdW_FD[I][J][K],dEdW_FD[I][J][K],-DW[I][J][K]/net->eta);


	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(DW[i][j]);
		}
		free(DW[i]);
	}
	free(DW);


	free(e);
	free(input);
	free(output);
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(dEdW_FD[i][j]);
		}
		free(dEdW_FD[i]);
	}
	free(dEdW_FD);
	return(1);
}


int ANN_dW( ANN * net , double *e ){

	int i,j,k,p,l,q;

	double *dEdv,*dEdvtmp,*dEdx;
	double *dydv,*dydvtmp,*dydx;
	double ***dEdW;
	double ****dydW;
	double *dXdW,temp;

	if( !( dEdW = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
		printf("\nError: malloc@ANN_dW(dEdW)\n");
		return(0);
	}
	for( i=0;i<net->N_layer+1;i++ ){
		if( !( dEdW[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
			printf("\nError: malloc@ANN_dW(dEdW)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i];j++ ){
			if( !( dEdW[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
				printf("\nError: malloc@ANN_dW(dEdW)\n");
				return(0);
			}
		}
	}
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				dEdW[i][j][k] = 0.;
			}
		}
	}

	/* Gradient error ( visible layer )*/

	if( !( dXdW = (double *)malloc( (net->N_neuron[net->N_layer])*sizeof(double)))){
		printf("\nError: malloc@ANN_dW(dXdW)\n");
		return(0);
	}

	for( k=0;k<net->N_neuron[net->N_layer];k++ ){
		for( l=0;l<net->N_neuron[net->N_layer+1];l++ ){
			for( j=0;j<net->N_output;j++ ){
				if( !ANN_dXdW( net , k , l , net->N_layer , dXdW ) ){
					printf("\nError: ANN_dXdW@ANN_dW\n");
					return(0);
				}
				if( net->N_layer == 0 ){
					temp = net->v[net->N_layer][k]*(l==j);
				}
				else{
					temp = ANN_InternalFunction(net->v[net->N_layer][k],net)*(l==j);
				}
				for( i=0;i<net->N_neuron[net->N_layer];i++ ){
					temp += dXdW[i]*net->W[net->N_layer][i][j];
				}
				temp = temp*ANN_InternalFunctionDer(net->v[net->N_layer+1][j],net);
				dEdW[net->N_layer][k][l] += -e[j]*temp ;
			}
		}
	}
	free(dXdW);

	if( !( dEdv = (double *)malloc( net->N_neuron[net->N_layer+1] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dW(dEdv)\n");
		return(0);
	}
	for( i=0;i<(net->N_neuron[net->N_layer+1] );i++ ){
		dEdv[i] = 0.;
	}
	for( i=0;i<(net->N_output);i++ ){
		dEdv[i] = -e[i]*ANN_InternalFunctionDer(net->v[net->N_layer+1][i],net);
	}

	/* HIDDEN LAYER*/
	if( !( dEdx = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dW(dEdx)\n");
		return(0);
	}
	if( !( dEdvtmp = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dW(dEdvtmp)\n");
		return(0);
	}
	if( !( dXdW = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dW(dXdW)BBBB\n");
		return(0);
	}

	for( k=0;k<net->N_layer;k++ ){
		if( !( dEdx = (double *)realloc( dEdx , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dEdx)\n");
			return(0);
		}
		if( !( dEdvtmp = (double *)realloc( dEdvtmp , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dEdvtmp)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
			dEdx[j] = 0.;
			dEdvtmp[j] = 0.;
		}

		for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
			for( p=0;p<net->N_neuron[net->N_layer+1-k];p++ ){
				dEdx[j] += dEdv[p]*net->W[net->N_layer-k][j][p];
			}
			dEdvtmp[j] = dEdx[j]*ANN_InternalFunctionDer(net->v[net->N_layer-k][j],net);
		}
		if( !( dXdW = (double *)realloc( dXdW , net->N_neuron[net->N_layer-k-1] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dXdW)CCCC\n");
			return(0);
		}
		for( i=0;i<net->N_neuron[net->N_layer-k-1];i++ ){
			for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
				if( !ANN_dXdW( net , i , j , (net->N_layer-k-1) , dXdW ) ){
					printf("\nError: ANN_dXdW@ANN_dW\n");
					return(0);
				}

				if( net->N_layer-k-1 != 0 ){
					temp = ANN_InternalFunction(net->v[net->N_layer-k-1][i],net);
				}
				else{
					temp = (net->v[0][i]);
				}

				for( p=0;p<(net->N_neuron[net->N_layer-k-1]);p++ ){
					temp += net->W[net->N_layer-k-1][p][j]*dXdW[p];
				}
				dEdW[net->N_layer-k-1][i][j] = dEdvtmp[j]*temp;
			}
		}
		if( !( dEdv = (double *)realloc( dEdv , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dEdv)\n");
			return(0);
		}
		for ( p=0;p<(net->N_neuron[net->N_layer-k]);p++ )
			dEdv[p] = dEdvtmp[p];
	}

	free(dEdx);
	free(dEdvtmp);
	free(dXdW);

	free(dEdv);
	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				net->dW[i][j][k] = -net->eta*dEdW[i][j][k] + net->rho*net->dW[i][j][k];
			}
		}
	}

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(dEdW[i][j]);
		}
		free(dEdW[i]);
	}
	free(dEdW);


	/* recurrent network*/
	if( net->r != 0 ){
		if( !( dydW = (double ****)malloc( net->N_neuron[net->N_layer+1] * sizeof(double ***) ) ) ){
			printf("\nError: malloc@ANN_dW(dydW)\n");
			return(0);
		}
		for( p=0;p<(net->N_neuron[net->N_layer+1]);p++ ){
			if( !( dydW[p] = (double ***)malloc( (net->N_layer+1) * sizeof(double **) ) ) ){
				printf("\nError: malloc@ANN_dW(dydW)\n");
				return(0);
			}
			for( k=0;k<net->N_layer+1;k++ ){
				if( !( dydW[p][k] = (double **)malloc( net->N_neuron[k] * sizeof(double *) ) ) ){
					printf("\nError: malloc@ANN_dW(dydW)\n");
					return(0);
				}
				for( j=0;j<net->N_neuron[k];j++ ){
					if( !( dydW[p][k][j] = (double *)malloc( net->N_neuron[k+1] * sizeof(double) ) ) ){
						printf("\nError: malloc@ANN_dW(dydW)\n");
						return(0);
					}
				}
			}
		}
		for( p=0;p<net->N_neuron[net->N_layer+1];p++ ){
			for( k=0;k<net->N_layer+1;k++ ){
				for( i=0;i<net->N_neuron[k];i++){
					for( j=0;j<net->N_neuron[k+1];j++){
						dydW[p][k][i][j] = 0.;
					}
				}
			}
		}
		/*Visible Layer*/
		if( !( dXdW = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dXdW)DDDD\n");
			return(0);
		}
		for( k=0;k<net->N_neuron[net->N_layer];k++ ){
			for( l=0;l<net->N_neuron[net->N_layer+1];l++ ){
				for( j=0;j<net->N_neuron[net->N_layer+1];j++ ){
					if( !ANN_dXdW( net , k , l , net->N_layer , dXdW ) ){
						printf("\nError: ANN_dXdW@ANN_dW\n");
						return(0);
					}
					if( net->N_layer == 0 ){
						temp = net->v[net->N_layer][k]*(l==j);
					}
					else{
						temp = ANN_InternalFunction(net->v[net->N_layer][k],net)*(l==j);
					}
					for( i=0;i<net->N_neuron[net->N_layer];i++ ){
						temp += dXdW[i]*net->W[net->N_layer][i][j];
					}
					temp = temp*ANN_InternalFunctionDer(net->v[net->N_layer+1][j],net);
					dydW[j][net->N_layer][k][l] = temp ;
				}
			}
		}
		free(dXdW);
		/*Hidden Layers*/
		if( !( dXdW = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dXdW)EEEE\n");
			return(0);
		}
		if( !( dydv = (double *)malloc( net->N_neuron[net->N_layer+1] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dydv)\n");
			return(0);
		}
		if( !( dydx = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dydx)\n");
			return(0);
		}
		if( !( dydvtmp = (double *)malloc( net->N_neuron[net->N_layer] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dydvtmp)\n");
			return(0);
		}

		for( q=0;q<net->N_neuron[net->N_layer+1];q++ ){
			for( i=0;i<(net->N_neuron[net->N_layer+1] );i++ ){
				dydv[i] = 0.;
			}
			dydv[q] = ANN_InternalFunctionDer(net->v[net->N_layer+1][q],net);

			for( k=0;k<net->N_layer;k++ ){
				if( !( dydx = (double *)realloc( dydx , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
					printf("\nError: malloc@ANN_dW(dydx)\n");
					return(0);
				}
				if( !( dydvtmp = (double *)realloc( dydvtmp , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
					printf("\nError: malloc@ANN_dW(dydvtmp)\n");
					return(0);
				}
				for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
					dydx[j] = 0.;
					dydvtmp[j] = 0.;
				}

				for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
					for( p=0;p<net->N_neuron[net->N_layer+1-k];p++ ){
						dydx[j] += dydv[p]*net->W[net->N_layer-k][j][p];
					}
					dydvtmp[j] = dydx[j]*ANN_InternalFunctionDer(net->v[net->N_layer-k][j],net);
				}
				if( !( dXdW = (double *)realloc( dXdW , net->N_neuron[net->N_layer-k-1] * sizeof(double) ) ) ){
					printf("\nError: malloc@ANN_dW(dXdW)\n");
					return(0);
				}
				for( i=0;i<net->N_neuron[net->N_layer-k-1];i++ ){
					for( j=0;j<net->N_neuron[net->N_layer-k];j++ ){
						if( !ANN_dXdW( net , i , j , (net->N_layer-k-1) , dXdW ) ){
							printf("\nError: ANN_dXdW@ANN_dW\n");
							return(0);
						}
						if( net->N_layer-k-1 != 0 ){
							temp = ANN_InternalFunction(net->v[net->N_layer-k-1][i],net);
						}
						else{
							temp = (net->v[0][i]);
						}

						for( p=0;p<(net->N_neuron[net->N_layer-k-1]);p++ ){
							temp += net->W[net->N_layer-k-1][p][j]*dXdW[p];
						}
						dydW[q][net->N_layer-1-k][i][j] = dydvtmp[j]*temp;
					}
				}

				if( !( dydv = (double *)realloc( dydv , net->N_neuron[net->N_layer-k] * sizeof(double) ) ) ){
					printf("\nError: malloc@ANN_dW(dydv)\n");
					return(0);
				}
				for ( p=0;p<(net->N_neuron[net->N_layer-k]);p++ )
					dydv[p] = dydvtmp[p];
			}
		}
		free(dydv);
		free(dydx);
		free(dydvtmp);
		free(dXdW);

		for( p=0;p<(net->r-1);p++ ){
			for( q=0;q<net->N_neuron[net->N_layer+1];q++ ){
				for( k=0;k<net->N_layer+1;k++ ){
					for( i=0;i<net->N_neuron[k];i++ ){
						for( j=0;j<net->N_neuron[k+1];j++ ){
							net->dy[net->r-1-p][q][k][i][j] = net->dy[net->r-2-p][q][k][i][j];
						}
					}
				}
			}
		}
		for( q=0;q<net->N_neuron[net->N_layer+1];q++ ){
			for( k=0;k<net->N_layer+1;k++ ){
				for( i=0;i<net->N_neuron[k];i++ ){
					for( j=0;j<net->N_neuron[k+1];j++ ){
						net->dy[0][q][k][i][j] = dydW[q][k][i][j];
					}
				}
			}
		}

		for( q=0;q<net->N_neuron[net->N_layer+1];q++ ){
			for( k=0;k<net->N_layer+1;k++ ){
				for( i=0;i<net->N_neuron[k];i++ ){
					free(dydW[q][k][i]);
				}
				free(dydW[q][k]);
			}
			free(dydW[q]);
		}
		free(dydW);
	}
	return(1);
}


int ANN_dXdW( ANN * net , int I , int J , int N , double *dXdW ){

	int i,j,k;
	double *dxdW,*dvdW;

	if( !( dxdW = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dXdW(dxdW)\n");
		return(0);
	}

	i = 0;
	for( j=0;j<net->N_input;j++ ){
		dxdW[i] = 0.;
		i++;
	}
	for( k=0;k<net->r;k++ ){
		for( j=0;j<net->N_neuron[net->N_layer+1];j++ ){
			dxdW[i] = net->dy[k][j][N][I][J];
			i++;
		}
	}
	if( !( dvdW = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_dW(dvdW)\n");
		return(0);
	}
	for( i=0;i<N;i++ ){
		if( !( dvdW = (double *)realloc( dvdW , net->N_neuron[i+1] * sizeof(double) ) ) ){
			printf("\nError: malloc@ANN_dW(dvdW)\n");
			return(0);
		}

		for( j=0;j<net->N_neuron[i+1];j++ ){
			dvdW[j] = 0.;
			for( k=0;k<net->N_neuron[i];k++ ){
				dvdW[j] += net->W[i][k][j]*dxdW[k];
			}
		}
		if( !( dxdW = (double *)realloc( dxdW,(net->N_neuron[i+1])*sizeof(double)))){
			printf("\nError: malloc@ANN_dW(dxdW)\n");
			return(0);
		}
		for( j=0;j<net->N_neuron[i+1];j++ ){
			dxdW[j] = ANN_InternalFunctionDer( net->v[i+1][j] , net )*dvdW[j];
		}
	}
	free(dvdW);
	for( i=0;i<net->N_neuron[N];i++ ){
		dXdW[i] = dxdW[i];
	}

	free(dxdW);

	return(1);
}


int ANN_TrainingEpoch( ANN * net , double ** INPUT , double ** DES_OUTPUT , double ** NN_OUTPUT, int N_sample ,int TRAINING_MODE , double ***DWtmp){

	int i,j,k,t,a;
	double *input, *output, *e;
	double ***DW;
	double E1,E2;


	if ( TRAINING_MODE == 1 ){
		if( !( DW = (double ***)malloc( (net->N_layer + 1) * sizeof(double **) ) ) ){
			printf("\nError: malloc@ANN_TrainingEpoch(DW)\n");
			return(0);
		}
		for( i=0;i<net->N_layer+1;i++ ){
			if( !( DW[i] = (double **)malloc( net->N_neuron[i] * sizeof(double *) ) ) ){
				printf("\nError: malloc@ANN_TrainingEpoch(DW)\n");
				return(0);
			}
			for( j=0;j<net->N_neuron[i];j++ ){
				if( !( DW[i][j] = (double *)malloc( net->N_neuron[i+1] * sizeof(double) ) ) ){
					printf("\nError: malloc@ANN_TrainingEpoch(DW)\n");
					return(0);
				}
			}
		}
		for( i=0;i<net->N_layer+1;i++ ){
			for( j=0;j<net->N_neuron[i];j++ ){
				for( k=0;k<net->N_neuron[i+1];k++ ){
					DW[i][j][k] = 0.;
				}
			}
		}
	}

	if( !( input = (double *)malloc( net->N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_TrainingEpoch(input)\n");
		return(0);
	}
	if( !( output = (double *)malloc( net->N_output * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_TrainingEpoch(output)\n");
		return(0);
	}
	if( !( e = (double *)malloc( net->N_output * sizeof(double) ) ) ){
		printf("\nError: malloc@ANN_TrainingEpoch(e)\n");
		return(0);
	}

	for( i=0;i<net->N_neuron[0];i++ ){
		input[i] = 0.;
	}
	for( i=0;i<net->N_output;i++ ){
		output[i] = 0.;
		e[i] = 0.;
	}

	E1 = 100000.;
	E2 = 100000.;

	for( t=0;t<N_sample;t++ ){
		E1 = E2;
		a = 0;
		for( i=0;i<net->N_input;i++ ){
			input[a] = INPUT[t][i];
			a++;
		}
		for( i=0;i<(net->r*net->N_neuron[net->N_layer+1]);i++ ){
			input[a] = net->yD[i];
			a++;
		}
		if( !ANN_sim( net,input,output ) ){
			printf("\nError: ANN_sim@ANN_TrainingEpoch\n");
			return(0);
		}
		for( i=0;i<net->N_output;i++ ){
			NN_OUTPUT[t][i] = output[i];
			e[i] = DES_OUTPUT[t][i] - NN_OUTPUT[t][i];
		}

		/*if ( TRAINING_MODE == 2 ){
			E2 = 0.;
			for( i=0;i<net->N_output;i++ ){
				E2 += .5*e[i]*e[i];
			}
			printf(".......E1 = %lf............E2 = %lf\n",E1,E2);
			printf("W(K)=%1.12lf\n",net->W[0][0][0]);
			getchar();
			while ( E2 > E1 ){
				net->eta = 0.5*net->eta;
				printf(".......E1 = %lf............E2 = %lf\n",E1,E2);
				printf("..eta = %lf\n",net->eta);
				getchar();
				for( i=0;i<net->N_layer+1;i++ ){
					for( j=0;j<net->N_neuron[i];j++ ){
						for( k=0;k<net->N_neuron[i+1];k++ ){
							net->dW[i][j][k] = -net->dW[i][j][k];
						}
					}
				}
				ANN_WeightUpdate( net,net->dW );
				printf("W(K-1)=%1.12lf\n",net->W[0][0][0]);
				for( i=0;i<net->N_layer+1;i++ ){
					for( j=0;j<net->N_neuron[i];j++ ){
						for( k=0;k<net->N_neuron[i+1];k++ ){
							net->dW[i][j][k] = -0.5*net->dW[i][j][k];
						}
					}
				}
				ANN_WeightUpdate( net,net->dW );
				printf("W(K)=%1.12lf\n",net->W[0][0][0]);
				if( !ANN_sim( net,input,output ) ){
					printf("\nError: ANN_sim@ANN_TrainingEpoch\n");
					return(0);
				}
				for( i=0;i<net->N_output;i++ ){
					NN_OUTPUT[t][i] = output[i];
					e[i] = DES_OUTPUT[t][i] - NN_OUTPUT[t][i];
				}
				E2 = 0.;
				for( i=0;i<net->N_output;i++ ){
					E2 += .5*e[i]*e[i];
				}
			}
		}*/
				
		if( !ANN_dW( net , e ) ){
			printf("\nError: ANN_dW@ANN_TrainingEpoch\n");
			return(0);
		}

		/*BATCH MODE*/
		if ( TRAINING_MODE == 1 ){
			for( i=0;i<net->N_layer+1;i++ ){
				for( j=0;j<net->N_neuron[i];j++ ){
					for( k=0;k<net->N_neuron[i+1];k++ ){
						DW[i][j][k] += net->dW[i][j][k];
					}
				}
			}
		}

		/*SEQUENTIAL MODE*/
		if ( TRAINING_MODE == 2 ){
			if( GR_CHECK == 1 ){
				printf("time: %d\n",t);
				if( !ANN_GradientCheck( net,INPUT,DES_OUTPUT,t ) ){
					printf("\nError: ANN_sim@ANN_GradientCheck\n");
					return(0);
				}
			}


			ANN_WeightUpdate( net,net->dW );
		}
	}

	if ( TRAINING_MODE == 1 ){
		for( i=0;i<net->N_layer+1;i++ ){
			for( j=0;j<net->N_neuron[i];j++ ){
				for( k=0;k<net->N_neuron[i+1];k++ ){
					DWtmp[i][j][k] = DW[i][j][k];
				}
			}
		}
		ANN_WeightUpdate( net,DW );
		for( i=0;i<net->N_layer+1;i++ ){
			for( j=0;j<net->N_neuron[i];j++ ){
				free(DW[i][j]);
			}
			free(DW[i]);
		}
		free(DW);
	}

	free(input);
	free(output);
	free(e);
	return(1);
}

void ANN_WeightUpdate( ANN * net , double *** DW ){

	int i,j,k;

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				net->W[i][j][k] += DW[i][j][k];
			}
		}
	}

}


int ANN_save( ANN * net){

	int i,j,k;

	FILE *fp;

	if ( !(fp = fopen( fnameTR , "w" )) ){
		printf("\nError: fopen@ANN_save\n");
		return(0);
	}

	fprintf(fp,"%d",net->N_input);
	fprintf(fp,"\n%d",net->N_output);
	fprintf(fp,"\n%d",net->N_layer);
	fprintf(fp,"\n%d\n",net->r);
	for( i=0;i<net->N_layer+1;i++ ){
		fprintf(fp,"%d ",net->N_neuron[i+1]);
	}
	fprintf(fp,"\n\n%1.20lf",net->alpha);
	fprintf(fp,"\n%lf",net->beta);
	fprintf(fp,"\n\n%1.20lf",net->eta);
	fprintf(fp,"\n%lf\n\n",net->rho);

	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			for( k=0;k<net->N_neuron[i+1];k++ ){
				fprintf(fp,"%1.20lf ",net->W[i][j][k]);
			}
			fprintf(fp,"\n");
		}
		fprintf(fp,"\n");
	}

	fclose(fp);

	return(1);
}

double ANN_InternalFunction( double v , ANN * net ){

	double y;

	y = net->alpha*tanh(net->beta*v);

	return(y);
}

double ANN_InternalFunctionDer( double v , ANN * net ){

	double y;

	y = net->alpha*net->beta*( 1 - tanh(net->beta*v)*tanh(net->beta*v));

	return(y);
}
