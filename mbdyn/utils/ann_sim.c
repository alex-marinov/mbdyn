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
char fnameANN[FILE_NAME_MAX_LENGTH];    /*ANN initilal definition file*/
char fnameINPUT[FILE_NAME_MAX_LENGTH];  /*ANN training input file*/
char fnameOUTPUT[FILE_NAME_MAX_LENGTH]; /*ANN training output file*/


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
	int r;        /* number of time delay*/
	double ***W;     /* network synaptic weights*/
	double *yD;
	}ANN;


/*FUNCTIONS' PROTOTYPES:    */

void ANN_Strcat(char *, char *, char *);
int ANN_FileName();
int ANN_Initialize( ANN * );
int ANN_InputData( double *** , double *** , int * , int);
void ANN_Free( ANN * , double *** , double *** , int );
int ANN_sim( ANN * , double * , double * );
void ANN_print( ANN * );
int ANN_OutputSave( double *** , int , int);
double ANN_InternalFunction( double , ANN * );



/* MAIN FUNCTION:    */
int main( int argc , char **argv ){
	ANN net;
	int i,j,p;
	double **INPUT, **NN_OUTPUT;
	double *input, *output;
	int N_sample;

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

	/* Artificial Neural Network inizialization*/
	printf("\nLOADING DATA...");
	if( !ANN_Initialize( &net ) ){
		printf("\nError: ANN_initilize@main\n");
		return(0);
	}
	/* Input data acquisition*/
	N_sample = 0;
	if( !ANN_InputData( &INPUT,&NN_OUTPUT,&N_sample,net.N_output ) ){
		printf("\nError: ANN_InputData@main\n");
		return(0);
	}

	if( !( input = (double *)malloc( net.N_neuron[0] * sizeof(double) ) ) ){
		printf("\nError: malloc@main\n");
		return(0);
	}
	if( !( output = (double *)malloc( net.N_output * sizeof(double) ) ) ){
		printf("\nError: malloc@main\n");
		return(0);
	}

	for( i=0;i<net.N_output;i++ ){
		output[i] = 0.;
	}
	for( i=0;i<net.N_neuron[0];i++ ){
		input[i] = 0.;
	}

	for( i=0;i<N_sample;i++ ){
		p = 0;
		for( j=0;j<net.N_input;j++ ){
			input[p] = INPUT[i][j];
			p++;
		}
		for( j=0;j<(net.r)*(net.N_neuron[net.N_layer+1]);j++ ){
			input[p] = net.yD[j];
			p++;
		}
		if( !ANN_sim(&net,input,output) ){
			printf("\nError: ANN_sim@main\n");
			return(0);
		}

		for( j=0;j<net.N_output;j++ ){
			NN_OUTPUT[i][j] = output[j];
		}

	}

	ANN_OutputSave( &NN_OUTPUT , N_sample , net.N_output);
	ANN_print(&net);

	/* dynamic memory free*/
	ANN_Free( &net,&INPUT,&NN_OUTPUT,N_sample );
	free(input);
	free(output);
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

	fclose(fp);
	return(1);
}

int ANN_Initialize( ANN *net ){

	int i,j,k;
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



	if( net->r != 0 ){
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


int ANN_InputData( double ***INPUT , double ***NN_OUTPUT , int * N_s , int N_output){

	int N_input,N_sample;
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


void ANN_Free( ANN * net , double *** INPUT  , double *** NN_OUTPUT , int N_sample ){

	int i,j;

	for( i=0;i<N_sample;i++ ){
		free((*INPUT)[i]);
		free((*NN_OUTPUT)[i]);
	}
	free(*INPUT);
	free(*NN_OUTPUT);

	if( net->r != 0 ){
		free(net->yD);
	}


	for( i=0;i<net->N_layer+1;i++ ){
		for( j=0;j<net->N_neuron[i];j++ ){
			free(net->W[i][j]);
		}
		free(net->W[i]);
	}
	free(net->W);

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
	printf("\nActivation function, beta: %lf ",net->beta);

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


int ANN_OutputSave( double ***NN_OUTPUT , int N_sample , int N_output){

	int i,j;

	FILE *fp;

	if ( !(fp = fopen( fnameOUTPUT , "w" )) ){
		printf("\nError: fopen@ANN_OutputSave(INPUT)\n");
		return(0);
	}


	for( i=0;i<N_sample;i++ ){
		for( j=0;j<N_output;j++ ){
			fprintf(fp,"%lf ",(*NN_OUTPUT)[i][j]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
	return(1);
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

	if( net->r != 0 ){
		for( i=0;i<(net->r-1)*(net->N_neuron[net->N_layer+1]);i++ ){
			net->yD[((net->r)*(net->N_neuron[net->N_layer+1]))-1-i] = net->yD[((net->r-1)*(net->N_neuron[net->N_layer+1]))-1-i];
		}
		for( i=0;i<net->N_neuron[net->N_layer+1];i++ ){
			net->yD[i] = y_p[i];
		}
	}

	free(IA);
	free(y);
	free(y_p);

	return(1);

}

double ANN_InternalFunction( double v , ANN * net ){

	double y;

	y = net->alpha*tanh(net->beta*v);
	//y = v;

	return(y);
}

