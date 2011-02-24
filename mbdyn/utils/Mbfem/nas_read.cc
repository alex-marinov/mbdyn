#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "nas_read.h"

using namespace std;

int remi;

int iopen(FILE* instrm, char* label){
	
	int key;
	int buf[15]; 
	char chbuf[60];

	// ogni  record e' preceduto e seguito da una parola che inidca di 
	// quanti i byte e composto il messaggio 
	rewind(instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 3) {
		printf("iopen bad key #1: %d\n",buf[1]);
		return -1; 
	}
	fread(buf, sizeof(int), 5, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 7) {
		printf("iopen bad key #2: %d\n",buf[1]);
		return -1; 
	}
	fread(&key, sizeof(int), 1, instrm);
	fread(chbuf, sizeof(char), key, instrm);	
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 2) {
		printf("iopen bad key #3: %d\n",buf[1]);
		return -1; 
	}
	fread(&key, sizeof(int), 1, instrm);
	fread(label, sizeof(char), key, instrm);	
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != -1) {
		printf("iopen bad key #4: %d\n",buf[1]);
		printf("OP2 is probably New!!\n");
		return -1; 
	}
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 0) {
		printf("iopen bad key #5: %d\n",buf[1]);
		return -1; 
	}
	return 0;
}

int iheadr(FILE* instrm, char* name, int* trailer){
	
	int key;
	int buf[10];
 
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 2) {
		if (fread(buf, sizeof(int), 1, instrm) == 0) {
			// eof
			return 1;
		} else {
			printf("iheadr bad key #1: %d\n",buf[1]);
			return -1; 
		}
	}
	// lettura label blocco dati
	fread(&key, sizeof(int), 1, instrm);
	fread(name, sizeof(char), key, instrm);
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != -1) {
		printf("iheadr bad key #2: %d\n",buf[1]);
		return -1; 
	}
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 7) {
		printf("iheadr bad key #3: %d\n",buf[1]);
		return -1; 
	}
	// lettura trailer 
	fread(&key, sizeof(int), 1, instrm);
	fread(trailer, sizeof(int), 7, instrm);
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != -2) {
		printf("iheadr bad key #4: %d\n",buf[1]);
		return -1; 
	}
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != 1) {
		printf("iheadr bad key #5: %d\n",buf[1]);
		return -1; 
	}
	fread(buf, sizeof(int), 3, instrm);
	switch (buf[1]) {
		case 0:
			break;
		case 1:
			printf("iheadr matrix column string\n");
			break;
		case 2:
			printf("iheadr factor matrix string\n");
			break;
		case 3:
			printf("iheadr factor matrix\n");
			break;
		default:
			printf("unknown case\n");
			return -1;	
	}
	fread(buf, sizeof(int), 3, instrm);
	//if (buf[1] != 7) {
	//	printf("iheadr bad key #7: %d\n",buf[1]);
	//	return -1; 
	//}
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(char), key, instrm);
	fread(&key, sizeof(int), 1, instrm);
	fread(buf, sizeof(int), 3, instrm);
	if (buf[1] != -3) {
		printf("iheadr bad key #8: %d\n",buf[1]);
		return -1; 
	}
	return 0;
}

int inasread(FILE* instrm, int* buffer, int bufsize, int iflag, int* nwred) 
{

//     INSTRM--INPUT--INPUT FILE
//     BUFFER--OUTPUT--ARRAY CONTAINING RETURNED WORDS
//     BUFSIZE--INPUT--INTEGER--THE NUMBER OF WORDS TO RETURN
//     IFLAG--INPUT--IF IFLAG IS 0 A NEW RECORD WILL BE READ
//		     IF IFLAG IS 1 A CONTINUATION READ WILL OCCUR
//                   IF IFLAG IS 2 THE REMAINDER OF THE RECORD WILL BE
//                                 SKIPPED
//                   IF IFLAG IS 3 ALL THE RECORD WILL BE
//                                 SKIPPED
//     NWRED--OUTPUT--IF THE RECORD DBLOCK NOT CONTAIN BUFSIZE REMAINING
//                  WORDS THEN NWRED WILL INDICATE THE ACTUAL NUMBER
//                  RETURNED
//     RETURN VALUE--OUTPUT--0 NORMAL RETURN
//                   1 END OF LOGICAL RECORD (NWR WILL BE SET)
//                   2 END OF DATA BLOCK
	int key;
	int buf[10];
	int i;

	switch (iflag) { 
		case 0: 
			// lettura di un nuovo record logico
			// se deve skippare dei dati non fa operazioni differenti
			// li legge in equal modo 
			if (remi >= 0) {
				// � un record logico nuovo
				if (bufsize < 0) { bufsize = -1*bufsize; }

				fread(buf, sizeof(int), 3, instrm);
				if (buf[1] != 1) {
					printf("iread bad key #1: %d\n",buf[1]);
					printf("This point is not the start of a new logical record\n");
					return -1; 
				}
				fread(buf, sizeof(int), 3, instrm);
				if (buf[1] != 0) {
					printf("iread bad key #2: %d\n",buf[1]);
					printf("This record is not a table\n");
					return -1; 
				}
				fread(buf, sizeof(int), 3, instrm);
			}
			if (buf[1] > 0) {
				// buf[1] e' il numeri di dati presenti
				if (buf[1] <= bufsize) {			
					// raggiunge il termine del record prima di aver i
					// riempito il buffer
					fread(&key, sizeof(int), 1, instrm);
					fread(buffer, sizeof(int), buf[1], instrm);
					fread(&key, sizeof(int), 1, instrm);
					// legge la chiusura del record
					*nwred = buf[1];
					fread(&key, sizeof(int), 1, instrm);
					fread(buf, sizeof(int), 1, instrm);
					fread(&key, sizeof(int), 1, instrm);
					if (buf[0] == 0) {
						remi = 0;
						return 2;
					} else if (buf[0] < 0) {
						remi = 0;
						return 1;
					} else {
						remi  = buf[0];
						fread(&key, sizeof(int), 1, instrm);
						return 0;
					}
				} else {
					// ci sono ancora dati che rimangono nel record
					fread(&key, sizeof(int), 1, instrm);
					fread(buffer, sizeof(int), bufsize, instrm);
					*nwred = bufsize;
					remi = buf[1] - bufsize;
					return 0;
				}
			} else if (buf[1] < 0) {
				// fine del record
				*nwred = 0;
				remi = 0;
				return 1;
			} else {
				// buf[1] = 0 fine del file
				*nwred = 0;
				remi = 0;
				return 2;
			}
			break;
		case 1:
			// continua a leggere un record
			if (remi <= bufsize) {
				// raggiunge il termine del record prima
				fread(buffer, sizeof(int), remi, instrm);
				fread(&key, sizeof(int), 1, instrm);
				// legge la chiiusura del record
				fread(&key, sizeof(int), 1, instrm);
				fread(buf, sizeof(int), 1, instrm);
				fread(&key, sizeof(int), 1, instrm);
				*nwred = remi;
				if (buf[0] == 0) {
					remi = 0;
					return 2;
				} else if (buf[0] < 0) {
					remi = 0;
					return 1;
				} else {
					remi  = buf[0];
					fread(&key, sizeof(int), 1, instrm);
					return 0;
				}
			} else {
				// ci sono ancora dati nel record
				fread(buffer, sizeof(int), bufsize, instrm);
				*nwred = bufsize;
				remi  -= bufsize;
				return 0;
			}
			break;
		case 2:
			// salta tante parole a meno che non sia alla fine del record
			if (remi == 0) {
				// legge la chiusura del record
				fread(&key, sizeof(int), 1, instrm);
				fread(buf, sizeof(int), 1, instrm);
				fread(&key, sizeof(int), 1, instrm);
				*nwred = 0;
				if (buf[0] == 0) {
					remi  = 0;
					return 2;
				} else if (buf[0] < 0) {
					remi  = 0;
					return 1;
				} else {
					remi  = buf[0];
					fread(&key, sizeof(int), 1, instrm);
					return 0;
				}	
			} else {
				// raggiunge il termine del record prima
				for (i=0; i < remi; i++) {
					fread(&buf[2], sizeof(int), 1, instrm);
				}
				fread(&key, sizeof(int), 1, instrm);
				// legge la chiusura del record
				fread(&key, sizeof(int), 1, instrm);
				fread(buf, sizeof(int), 1, instrm);
				fread(&key, sizeof(int), 1, instrm);
				*nwred = remi;
				if (buf[0] == 0) {
					remi  = 0;
					return 2;
				} else if (buf[0] < 0) {
					remi  = 0;
					return 1;
				} else {
					remi  = buf[0];
					fread(&key, sizeof(int), 1, instrm);
					return 0;
				}	
			} 
			break;
		case 3:	
			// salta bufsize parole
 			if (remi >= 0) {
 				// il record è nuovo e va saltato tutto
 				if (remi == 0) {
					fread(buf, sizeof(int), 3, instrm);
					if (buf[1] != 1) {
						printf("iread bad key #1: %d\n",buf[1]);
						printf("This point is not the start of a new logical record\n");
						return -1; 
					}
					fread(buf, sizeof(int), 3, instrm);
					if (buf[1] != 0) {
						printf("iread bad key #2: %d\n",buf[1]);
						printf("This record is not a table\n");
						return -1; 
					}
					fread(buf, sizeof(int), 3, instrm);
				} else {
					buf[1] = remi;
				}
				if (buf[1] > 0) {
					if (remi == 0) fread(&key, sizeof(int), 1, instrm);
					for (i=0; i < buf[1]; i++) {
						fread(&buf[2], sizeof(int), 1, instrm);
					}
					fread(&key, sizeof(int), 1, instrm);
					// legge la chiusura del record
					*nwred = buf[1];
					fread(&key, sizeof(int), 1, instrm);
					fread(buf, sizeof(int), 1, instrm);
					fread(&key, sizeof(int), 1, instrm);
					if (buf[0] == 0) {
						remi  = 0;
						return 2;
					} else if (buf[0] < 0) {
						remi  = 0;
						return 1;
					} else {
						remi  = buf[0];
						fread(&key, sizeof(int), 1, instrm);
						return 0;
					}
				} else if (buf[1] < 0) {
					// fine del record
					*nwred = 0;
					remi = 0;
					return 1;
				} else {
					// buf[1] = 0 fine del file
					*nwred = 0;
					remi = 0;
					return 2;
				}
			} else {
				printf("iread: the previous record is not finished\n");
				printf("This point is not the start of a new logical record\n");
				return -1; 
			}	
			break;
		default:
			printf("Unknown iread case: dunno wath to do!!\n");
			return -1;
	}
	// non dovrebbe mai arrivare qui!!	
	return 0;
}

int namecomp(char* label, char** namelist, int listsize) {

	int i = 0;
	while (i < listsize) {
		if (strncmp(label, namelist[i],5) == 0) {
			return i;
		}
		i++;
	}
	return -1;
}

int readBgpdt(FILE* in, int* buffer, int bufsize, map<int,int> iminFem,map<int,int> imaxFem, 
			map<int,int> iminFem_l,map<int,int> imaxFem_l ,
			map<int,Block>* BList,vector<int>* femidnode) {

	int readType;
	int readWord;
	float*  flBuffer = (float*)buffer;
	double* dbBuffer = (double*)buffer;
	int* intbuffer;
	int newbufsize;
	int count;
	int err;
	div_t result;
	int lastn,i;

	readType = 0;
	newbufsize = (int)(floor((double)bufsize/12.0)*12);
	
	err = inasread(in, buffer, newbufsize, readType, &readWord);
	lastn = 0;
	if (err == 2) return err;
	while (err != 1)  {
		// faccio l'output dei dati
		// ogni nodo contiene 12 dati interi
		count = 0;
		result = div(readWord+lastn,12);
		lastn = 0;
		if (result.rem != 0) {
			lastn = result.rem;
		}
		while (count < readWord-lastn) {
			intbuffer = &buffer[count+6];
			dbBuffer = (double*)intbuffer;
			//printf("GRID    %8d %8d %11.6g %11.6g %11.6g\n",buffer[count+2],buffer[count],
			//	dbBuffer[0],dbBuffer[1],dbBuffer[2]);
			int label = buffer[count+2];
			Vec3 pos(dbBuffer[0],dbBuffer[1],dbBuffer[2]);					
			// Assegnazione nel blocco
			for (map<int,int>::iterator it_id = iminFem.begin(); it_id != iminFem.end(); it_id++){
				if ((label>=(it_id->second))&&(label<=imaxFem[it_id->first])){
					// Inserimento nel blocco
					(*BList)[it_id->first].SetNode(label,pos,false,false);
					femidnode->push_back(it_id->first);
				}
				if ((label>=iminFem_l[it_id->first])&&(label<=imaxFem_l[it_id->first])){
					// Inserimento nel blocco - nodo per carico concentrato
					(*BList)[it_id->first].SetNode(label,pos,false,true);
				}
			}
			count += 12;
		}
		readType = 1;
		if (lastn > 0) {
			for(i=0; i < lastn; i++) {
				buffer[i] = buffer[count+i];
			}
		}
		err = inasread(in, &(buffer[lastn]), newbufsize, readType, &readWord);
		if (err == 2) return err;
	}
	// legge i dati rimanenti il cui computo è in readWord
	count = 0;
	while(count < readWord) {
		intbuffer = &buffer[count+6];
		dbBuffer = (double*)intbuffer;
		//printf("GRID    %8d %8d %11.6g %11.6g %11.6g\n",buffer[count+2],buffer[count],
		//	dbBuffer[0],dbBuffer[1],dbBuffer[2]);
		int label = buffer[count+2];
		Vec3 pos(dbBuffer[0],dbBuffer[1],dbBuffer[2]);					
		// Assegnazione nel blocco
		for (map<int,int>::iterator it_id = iminFem.begin(); it_id != iminFem.end(); it_id++){
			if ((label>=(it_id->second))&&(label<=imaxFem[it_id->first])){
				// Inserimento nel blocco
				(*BList)[it_id->first].SetNode(label,pos,false,false);
				femidnode->push_back(it_id->first);
			}
			if ((label>=iminFem_l[it_id->first])&&(label<=imaxFem_l[it_id->first])){
				// Inserimento nel blocco - nodo per carico concentrato
				(*BList)[it_id->first].SetNode(label,pos,false,true);
			}
		}
		count += 12;
	 }
	readType = 3;
	err = inasread(in, buffer, bufsize, readType, &readWord);
	while (err != 2) {
		readType = 3;
		err = inasread(in, buffer, bufsize, readType, &readWord);
	}
	return err;
}

int readop2(char* InFilename,map<int,int> iminFem,map<int,int> imaxFem,
		map<int,int> iminFem_l ,map<int,int> imaxFem_l ,
		map<int,Block>* BList,vector<int>* femidblock) {

	FILE *Infid;
	int err;
	char label[8];
	char dabname[8];
	int  trailer[7];
	int bufsize = 2000;
	int buffer[bufsize];
	//char*   cbuffer = buffer;
	//float*  fbuffer = buffer;
	//double* dbuffer = buffer;
	int readType,readWord;
	int blocknmsz = 1;
	char* blocknames[blocknmsz];
	
	blocknames[0] = (char*)malloc(sizeof(char)*6);
	//malloc(blocknames[1], sizeof(char), 6);
	//malloc(blocknames[2], sizeof(char), 6);
	blocknames[0] = "BGPDT";			
 			
	Infid  = fopen(InFilename, "r");
	
	if (Infid == NULL) {
		err = ferror(Infid);	
		printf("Error #%d in opening the input file!!\n",err);
		return -1;
	}
	 	
	// no head info for post -2 files
	//err = iopen(Infid, label);
	//if (err != 0) {
	//	printf("Error in iopen !!\n");
	//	return -1;
	//}

	while (!feof(Infid)) {
		err = iheadr(Infid, dabname, trailer);
		if (err == -1) {
			printf("Error in iheadr !!\n");
			return -1;
		}
		if (err == 1) {
			fclose(Infid);
			// EOF
			return 0;
		}
		switch (namecomp(dabname,blocknames,blocknmsz)) {
	
			case 0:
				err = readBgpdt(Infid,buffer,bufsize,iminFem,imaxFem,iminFem_l,imaxFem_l,BList,femidblock);
				break;
			default:
				// blocco non interessante
				// salta fino alla fine del blocco 
				readType = 3;
				err = inasread(Infid, buffer, bufsize, readType, &readWord);
				while (err != 2) {
					readType = 3;
					err = inasread(Infid, buffer, bufsize, readType, &readWord);
					if (err == -1) {
						printf("Unknown Error\n");
						fclose(Infid);		
						return -1;
					}
				}
				break;
		}	
	}
	fclose(Infid);
	return 0;
}

FILE* open_nasout(char* InFileName,int n_step,vector<int> step_list){
	FILE* fnas;

	char filenas[100];

	strcpy(filenas,InFileName);
	strcat(filenas,".nas");

	fnas = fopen(filenas,"w");

	fprintf(fnas,"$ NASTRAN cards\n");
	fprintf(fnas,"$ CASE CONTROL DECK\n");
	fprintf(fnas,"SOL 101\nTIME 10000\nCEND\n");
	for (int i_lc = 1; i_lc <= n_step; i_lc++){
		fprintf(fnas,"$STEP %d\n",step_list[i_lc-1]);
		fprintf(fnas,"SUBCASE %d\n",i_lc);
		fprintf(fnas,"ECHO = NONE\nDISPLACEMENT(PLOT) = ALL\nOLOAD(PLOT) = ALL\n");
		fprintf(fnas,"SPCFORCE(PLOT) = ALL\nFORCE(PLOT,CORNER) = ALL\nSTRESS(PLOT,CORNER) = ALL\n");
		fprintf(fnas,"$Edit the constraint set according to your model\n");
		fprintf(fnas,"SPC = 1\n");
		fprintf(fnas,"LOAD = %d\n",i_lc);
	}
	return fnas;
}

void nas_cards(FILE* fnas,map<int,Block> BlList,int n_step, int CID, int start_lc){
	fprintf(fnas,"$ BULK DATA\n");
	double A = 1.;
	char buf[10];
	// Calcolo del numero di schede presenti, per non fare casino coi numeri
	int icards = int(pow(10.,(int(n_step/10)+1)));
	for (int i_lc = 1; i_lc <= n_step; i_lc++){
		int load_type = 0;
		fprintf(fnas,"$ STEP %d - ACCELERATIONS \n",i_lc);
		bool ac_print = false, cl_print = false, ic_print = false;
		for (map<int,Block>::iterator it_bl = BlList.begin(); it_bl != BlList.end(); it_bl++){
			fprintf(fnas,"$ BLOCK %d\n",it_bl->first);
			// Riceve i nodi FEM
			vector<Node> nodes = it_bl->second.GetFemNodes();
			vector<Node>::iterator it_node = nodes.begin();
			// Riceve le accelerazioni sui nodi FEM
			const FullMatrixHandler* Acc = it_bl->second.GetFemAcc(i_lc);
			// Controllino
			if (Acc->iGetNumRows() != nodes.size())
				cout << "Qualcosa è andato storto, pirla!" << endl;
			// Scrittura delle accelerazioni
			for (int i=1; i<=Acc->iGetNumRows(); i++){
				if (sqrt((*Acc)(i,1)*(*Acc)(i,1)+(*Acc)(i,2)*(*Acc)(i,2)+(*Acc)(i,3)*(*Acc)(i,3))>=1.0e-5) {
					if (!ac_print) {
						ac_print = true;
						load_type++;
					}
					fprintf(fnas,"ACCEL1  ");
					fprintf(fnas,"%8d%8d",icards*(i_lc+start_lc-1)+load_type,CID);
					sprintf(buf,"% 8f",A);
					fprintf(fnas,"%8.8s",buf);
					for (int j=0;j<3;j++){
						if ((*Acc)(i,j+1)<0){
							// 7 bit available
							if ((*Acc)(i,j+1)<(-1e6+1)){
								sprintf(buf,"%4.2e",-(*Acc)(i,j+1));
							} else {
								sprintf(buf,"% 8f",-(*Acc)(i,j+1));
							}
						} else {
							// 6 bit available
							if ((*Acc)(i,j+1)>(1e5-1)){
								sprintf(buf,"%4.1e",-(*Acc)(i,j+1));
							} else {
								sprintf(buf,"% 8f",-(*Acc)(i,j+1));
							}
						}
						fprintf(fnas,"%8.8s",buf);
					}
					fprintf(fnas,"                +\n");
					fprintf(fnas,"+       %8d\n",it_node->GetLabel());
				}
				if (it_node+1 != nodes.end()){
					it_node++;
				}
			}
		}
		fprintf(fnas,"$ STEP %d - CONCENTRATED LOADS \n",i_lc);
		for (map<int,Block>::iterator it_bl = BlList.begin(); it_bl != BlList.end(); it_bl++){
			fprintf(fnas,"$ BLOCK %d\n",it_bl->first);
			vector<CLoad> forces = it_bl->second.GetCLoads(i_lc);
			for (vector<CLoad>::iterator it_f = forces.begin(); it_f != forces.end(); it_f++){
				int labtemp = it_f->GetLabFem();
				Vec3 fortemp = it_f->GetLoad();
				bool iscouple = it_f->IsCouple();
				if (sqrt(fortemp(1)*fortemp(1)+fortemp(2)*fortemp(2)+fortemp(3)*fortemp(3))>=1e-5){
					if (!cl_print) {
						cl_print = true;
						load_type++;
					}
					if (iscouple){
						fprintf(fnas,"MOMENT  %8d%8d%8d",icards*(i_lc+start_lc-1)+load_type,labtemp,CID);
						sprintf(buf,"% 8f",A);
						fprintf(fnas,"%8.8s",buf);
						for (int i=0;i<3;i++){
							sprintf(buf,"% 8f",fortemp(i+1));
							fprintf(fnas,"%8.8s",buf);
						}
						fprintf(fnas,"\n");
					} else {						
						fprintf(fnas,"FORCE   %8d%8d%8d",icards*(i_lc+start_lc-1)+load_type,labtemp,CID);
						sprintf(buf,"% 8f",A);
						fprintf(fnas,"%8.8s",buf);
						for (int i=0;i<3;i++){
							sprintf(buf,"% 8f",fortemp(i+1));
							fprintf(fnas,"%8.8s",buf);
						}
						fprintf(fnas,"\n");
					}
				}
			}
		}
		fprintf(fnas,"$ INERTIA COUPLE\n");
		for (map<int,Block>::iterator it_bl = BlList.begin(); it_bl != BlList.end(); it_bl++){
			if (it_bl->second.BeOnBe() && it_bl->second.isMM()){
				fprintf(fnas,"$BLOCK %d\n",it_bl->first);
				vector<Node> nodes = it_bl->second.GetFemNodes();
				vector<Node>::iterator it_node = nodes.begin();
				// Riceve le coppie sui nodi FEM
				const FullMatrixHandler* Cou = it_bl->second.GetFemInCou(i_lc);
				// Scrittura delle accelerazioni
				for (int i=1; i<=Cou->iGetNumRows(); i++){
					if (sqrt((*Cou)(i,1)*(*Cou)(i,1)+(*Cou)(i,2)*(*Cou)(i,2)+(*Cou)(i,3)*(*Cou)(i,3))>=1e-5) {
						if (!ic_print) {
							ic_print = true;
							load_type++;
						}
                                                fprintf(fnas,"MOMENT  %8d%8d%8d",icards*(i_lc+start_lc-1)+load_type,it_node->GetLabel(),CID);
                                                sprintf(buf,"% 8f",A);
                                                fprintf(fnas,"%8.8s",buf);
                                                for (int j=0;j<3;j++){                                                         
						 	sprintf(buf,"% 8f",-(*Cou)(i,j+1));
                                                        fprintf(fnas,"%8.8s",buf);
                                                }
                                                fprintf(fnas,"\n");
					}
					if (it_node+1 != nodes.end()){
						it_node++;
					}
				}
			}
		}
		if (load_type != 0){
			fprintf(fnas,"$ LOAD CARDS\n");
			fprintf(fnas,"LOAD    %8d",i_lc+start_lc-1);
			sprintf(buf,"% 8f",A);
			fprintf(fnas,"%8.8s",buf);
			for (int i = 1; i<=load_type; i++){
				if ((i)%4==0){
					fprintf(fnas,"+\n");
					fprintf(fnas,"        ");
				}
				sprintf(buf,"% 8f",A);
				fprintf(fnas,"%8.8s",buf);
				fprintf(fnas,"%8d",icards*(i_lc+start_lc-1)+i);
			}
			fprintf(fnas,"\n");
		}
	}
	fclose(fnas);
}

void text_dump(char* InFileName,int N_step,int N_block,map<int,Block> BlockList){
	// Scrive un file di testo da leggere con Matlab o simili
	FILE* fout;

	char fileout[100];

	strcpy(fileout,InFileName);
	strcat(fileout,".out");

	fout = fopen(fileout,"w");	
	
	// Prima riga - contiene il numero di blocchi e il numero di step
	fprintf(fout,"%d %d\n",N_block,N_step);
	for (map<int,Block>::iterator it_bl = BlockList.begin(); it_bl != BlockList.end(); it_bl++){
		int ngdl = 3;
		if (it_bl->second.BeOnBe()){
			ngdl = 6;
		}
		fprintf(fout,"%d %d %d %d\n",it_bl->second.GetId(),it_bl->second.GetNmb(),it_bl->second.GetNfem(),ngdl);
		// Nodi Fem
		vector<Node> mb = it_bl->second.GetMbNodes();
		vector<Node> fem = it_bl->second.GetFemNodes();
		// Plottaggio dei nodi mb
		int count_mb = 1;
		for (vector<Node>::iterator it_m = mb.begin(); it_m != mb.end(); it_m++) {
			Vec3 tmp = it_m->GetPos();
			fprintf(fout,"%d %le %le %le ",it_m->GetLabel(),tmp(1),tmp(2),tmp(3));
			// Plottaggio delle accelerazioni
			for (int i_lc = 1; i_lc <= N_step; i_lc++){
				Vec3 tmpacc = it_bl->second.GetMbAcc(i_lc,count_mb);
				fprintf(fout,"%le %le %le ",tmpacc(1),tmpacc(2),tmpacc(3));
				if (ngdl==6){
					Vec3 tmpracc = it_bl->second.GetMbRAcc(i_lc,count_mb);
					fprintf(fout,"%le %le %le ",tmpracc(1),tmpracc(2),tmpracc(3));
				}
			}
			count_mb++;
			fprintf(fout,"\n");
		}
		// Plottaggio dei nodi fem
		int count_fem = 1;
		for (vector<Node>::iterator it_f = fem.begin(); it_f != fem.end(); it_f++) {
			Vec3 tmp = it_f->GetPos();
			fprintf(fout,"%d %le %le %le ",it_f->GetLabel(),tmp(1),tmp(2),tmp(3));
			// Plottaggio delle accelerazioni
			for (int i_lc = 1; i_lc <= N_step; i_lc++){
				Vec3 tmpacc = it_bl->second.GetFemAcc(i_lc,count_fem);
				fprintf(fout,"%le %le %le ",tmpacc(1),tmpacc(2),tmpacc(3));
				if (ngdl==6){
					Vec3 tmpracc = it_bl->second.GetFemRAcc(i_lc,count_fem);
					fprintf(fout,"%le %le %le ",tmpracc(1),tmpracc(2),tmpracc(3));
				}
			}
			count_fem++;
			fprintf(fout,"\n");
		}
		// Plottaggio della matrice di interfaccia
		SpMapMatrixHandler* pH;
		pH = it_bl->second.GetH();
		fprintf(fout,"%d %d\n",pH->iGetNumRows(),pH->iGetNumCols());
		for (int i=0; i < pH->iGetNumRows(); i++){
			for (int j=0; j < pH->iGetNumCols(); j++){
				fprintf(fout,"%le ",pH->operator()(i+1,j+1));
			}
			fprintf(fout,"\n");
		}
	}
	fclose(fout);
}

void readmat(char* InFileName, vector<int> femidblock, map<int,Block>* BlList){
	 
	FILE *Infid;
	int err;
	int buf[5];
	char cname[8];
	 
	Infid = fopen(InFileName,"r");
	if (Infid == NULL) {
		err = ferror(Infid);
	 	printf("Error #%d in opening the input file!!\n",err);
	}

	fread(buf,sizeof(int),5,Infid);
	fread(cname,sizeof(char),8,Infid);

	int Nrow = buf[1], Ncol = buf[2];
	int Ntype = buf[4];
	int nfem = femidblock.size();

	// Controllo: matrice di massa è un vettore colonna
	if (Ncol/nfem != 6){
		std::cout << Ncol << " " << nfem << std::endl;
		std::cerr << "Size not matching - aborting" << std::endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	// Contatori per l'inserimento nei blocchi
	std::map<int,int> nblock;
	for (std::map<int,Block>::iterator it_bl = BlList->begin(); it_bl != BlList->end(); it_bl++){
		nblock[it_bl->first] = 0;
	}
	int row[5];
	fread(row,sizeof(int),5,Infid);
	int N_el = row[4]/Ntype;
	for (int i=0; i<nfem-1; i++){
			double col[6];
			fread(col,sizeof(double),6,Infid);
			// Consider only angular inertia
			if ((*BlList)[femidblock[i]].BeOnBe()){
				Vec3 tmp(col[3],col[4],col[5]);
				(*BlList)[femidblock[i]].SetFemMM(tmp,nblock[femidblock[i]]+1);
				nblock[femidblock[i]]++;
			}
	}
	// Trattiamo l'ultimo record
	// Calcolo dei nodi residui
	int rem = N_el-(nfem-1)*6; 
	double col[rem];
	fread(col,sizeof(double),rem,Infid);
	if ((*BlList)[femidblock[nfem-1]].BeOnBe()){
		switch (rem){
			case 4:{
				Vec3 tmp(col[3], 0., 0.);
				(*BlList)[femidblock[nfem-1]].SetFemMM(tmp,nblock[femidblock[nfem-1]]+1);
				nblock[femidblock[nfem-1]]++;
				break;
				}
			case 5:{
				Vec3 tmp(col[3],col[4], 0.);
				(*BlList)[femidblock[nfem-1]].SetFemMM(tmp,nblock[femidblock[nfem-1]]+1);
				nblock[femidblock[nfem-1]]++;
				break;
				}
			case 6:{
				Vec3 tmp(col[3],col[4],col[5]);
	                        (*BlList)[femidblock[nfem-1]].SetFemMM(tmp,nblock[femidblock[nfem-1]]+1);
				nblock[femidblock.size()-1]++;
				break;
				}
			default:{
				Vec3 tmp(Zero3);
				(*BlList)[femidblock[nfem-1]].SetFemMM(tmp,nblock[femidblock[nfem-1]]+1);
				nblock[femidblock[nfem-1]]++;
				break;
				}
		}
	}
	fclose(Infid);
}

void rdir_dump(char* InFileName,int N_step,int N_block,map<int,Block> BlockList){
	// Scrive un file di testo da leggere con Matlab o simili
	FILE* fout;

	char fileout[100];

	strcpy(fileout,InFileName);
	strcat(fileout,".rdi");

	fout = fopen(fileout,"w");

	// Prima riga - contiene il numero di blocchi e il numero di step
	fprintf(fout,"%d %d\n",N_block,N_step);
	for (map<int,Block>::iterator it_bl = BlockList.begin(); it_bl != BlockList.end(); it_bl++){
		int ngdl = 3;
		if (it_bl->second.BeOnBe()){
			ngdl = 6;
		}
		fprintf(fout,"%d %d %d %d\n",it_bl->second.GetId(),it_bl->second.GetNmb(),it_bl->second.GetNfem(),ngdl);
		//Nodi Fem
		vector<Node> mb = it_bl->second.GetMbNodes();
		vector<Node> fem = it_bl->second.GetFemNodes();
		// Plottaggio dei nodi mb
		int count_mb = 1;
		for (vector<Node>::iterator it_m = mb.begin(); it_m != mb.end(); it_m++) {
			Vec3 tmp = it_m->GetPos();
			fprintf(fout,"%d %le %le %le ",it_m->GetLabel(),tmp(1),tmp(2),tmp(3));
			// Plottaggio degli spostamenti
			for (int i_lc = 1; i_lc <= N_step; i_lc++){
				Vec3 tmppos = it_bl->second.GetMbDisp(i_lc,count_mb);
		 		fprintf(fout,"%le %le %le ",tmppos(1),tmppos(2),tmppos(3));
				if (ngdl==6){
					Vec3 tmprot = it_bl->second.GetMbRot(i_lc,count_mb);
					fprintf(fout,"%le %le %le ",tmprot(1),tmprot(2),tmprot(3));
				}
			}
			count_mb++;
			fprintf(fout,"\n");
		}
		// Plottaggio dei nodi fem
		int count_fem = 1;
		for (vector<Node>::iterator it_f = fem.begin(); it_f != fem.end(); it_f++) {
			Vec3 tmp = it_f->GetPos();
			fprintf(fout,"%d %le %le %le ",it_f->GetLabel(),tmp(1),tmp(2),tmp(3));
			// Plottaggio delle accelerazioni
			for (int i_lc = 1; i_lc <= N_step; i_lc++){
				Vec3 tmppos = it_bl->second.GetFemDisp(i_lc,count_fem);
				fprintf(fout,"%le %le %le ",tmppos(1),tmppos(2),tmppos(3));
				if (ngdl==6){
					Vec3 tmprot = it_bl->second.GetFemRot(i_lc,count_fem);
					fprintf(fout,"%le %le %le ",tmprot(1),tmprot(2),tmprot(3));
				}
			}
			count_fem++;
			fprintf(fout,"\n");
		}
	}
	fclose(fout);
}
	
