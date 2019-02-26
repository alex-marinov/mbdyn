#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "mbdyn_out.h"

using namespace std;

map<int,int> read_mov(char* InFileName,map<int,int> iminMb,
	map<int,int> imaxMb,vector<int> step_list,
	map<int,Block>* BList,bool disp, double scale_f){
	
	FILE* fmov;
	map <int,int> nodeblock;
	
	// Nome del file da aprire
	char filemov[100],buf[1024];
	
	strcpy(filemov,InFileName);
	strcat(filemov,".mov");
	
	fmov = fopen(filemov,"r");
	if (fmov == NULL) {
		cerr << "Invalid path or filename - Aborting" << endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	// Inizio della lettura delle posizioni iniziali dei nodi
	int n1 = -1, N_mb = 0;
	int label;
	bool cont = true;
	while (cont){
		double x,y,z;
		fgets(buf,1024,fmov);
		istringstream iss(buf,istringstream::in);
		iss >> label >> x >> y >> z;
		if (N_mb == 0) n1 = label;
		if ((N_mb != 0)&(label==n1)) cont = false; 	
		else {
			Vec3 pos(x,y,z);					
			// Assegnazione nel blocco
			for (map<int,int>::iterator it_id = iminMb.begin(); it_id != iminMb.end(); it_id++){
				if ((label>=(it_id->second))&&(label<=imaxMb[it_id->first])){
					// Inserimento nel blocco
					(*BList)[it_id->first].SetNode(label,pos*scale_f,true,false);
					nodeblock[label]=it_id->first;
				}
			}
			N_mb++;
		}
	}
	// Ricavate e assegnate ai blocchi, si può procedere alla lettura delle accelerazioni
	rewind(fmov);
	int ini_step = 0;
	// Raggiunge il primo step e inizia la lettura
	int i_step = 1;
	for (vector<int>::iterator it_step = step_list.begin(); it_step != step_list.end(); it_step++){
		int delta_step = (*it_step - ini_step) - 1;
		for (int i=0;i<delta_step*N_mb;i++) fgets(buf,1024,fmov);
		// Inizio la lettura
		for (int i=0;i<N_mb;i++){
			fgets(buf,1024,fmov);
			istringstream iss(buf,istringstream::in);
			int label;
			double dummy;
			Vec3 Pos(Zero3),Rot(Zero3),Acc(Zero3),AngVel(Zero3),AngAcc(Zero3);
			if (iss >> label >> Pos(1) >> Pos(2) >> Pos(3) /*poszioni*/ >> Rot(1) >> Rot(2) >> Rot(3) >> /*rotazioni*/ dummy >> dummy >> dummy >> /*velocità*/ AngVel(1) >> AngVel(2) >> AngVel(3) /*velangolari*/ >> Acc(1) >> Acc(2) >> Acc(3) >> /*AccAngolari*/ AngAcc(1) >> AngAcc(2) >> AngAcc(3)){
			// Inseriamo l'accelerazione nel blocco
				if (nodeblock.find(label)!=nodeblock.end()){
					(*BList)[nodeblock[label]].SetMbAcc(i_step,label,Acc*scale_f,AngVel,AngAcc);
					if (disp){
						(*BList)[nodeblock[label]].SetMbDisp(i_step,label,Pos*scale_f,Rot);
					}
				}
			} else {
				if (nodeblock.find(label)!=nodeblock.end()){
					(*BList)[nodeblock[label]].SetMbAcc(i_step,label,Acc*scale_f,AngVel,AngAcc);
					if (disp){
						(*BList)[nodeblock[label]].SetMbDisp(i_step,label,Pos*scale_f,Rot);
					}
				}
			}
		}
		ini_step = *it_step;
		i_step++;
	}
	fclose(fmov);
	return nodeblock;	
}

void read_frc(char* InFileName,map<int,int> idblock,vector<int> step_list,
		map<int,int> loads,map<int,Block>* BList, double scale_f){
	
	FILE* ffrc;
	
	// Nome del file da aprire
	char filefrc[100],buf[1024];
	
	strcpy(filefrc,InFileName);
	strcat(filefrc,".frc");
	
	ffrc = fopen(filefrc,"r");
	if (ffrc == NULL) {
		cerr << "Invalid path or filename - Aborting" << endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	
	// Acquisizione delle forze presenti nel modello
	// Legge i primi passi e individua il numero di forze
	int n1 = -1, N_load = 0;
	int label;
	bool cont = true;
	while (cont) {
		fgets(buf,1024,ffrc);
		istringstream iss(buf,istringstream::in);
		iss >> label;
		if (N_load==0) n1 = label;
		if ((N_load!=0)&(label==n1)) cont = false; 	
		else N_load++;
	}
	// A questo punto iniziamo ad acquisire i vari load step
	rewind(ffrc);
	int ini_step = 0;
	// Raggiunge il primo step e inizia la lettura
	int i_step = 1;
	for (vector<int>::iterator it_step = step_list.begin(); it_step != step_list.end(); it_step++){
		int delta_step = (*it_step - ini_step) - 1;
		for (int i=0;i<delta_step*N_load;i++) fgets(buf,1024,ffrc);
		// Inizio la lettura
		for (int i=0;i<N_load;i++){
			fgets(buf,1024,ffrc);
			istringstream iss(buf,istringstream::in);
			int node_label,lab_force;
			double dummy;
			Vec3 Load(Zero3);
			if (iss >> lab_force >> node_label >> dummy >> Load(1) >> Load(2) >> Load(3) >> dummy >> dummy >> dummy){
				// E' una forza
				if (loads.find(node_label) != loads.end()){
					// Forza NON è di volume
					(*BList)[idblock[node_label]].SetMbLoad(i_step,node_label,Load,false,loads[node_label]);
				} else (*BList)[idblock[node_label]].SetMbLoad(i_step,node_label,Load,false,0);		
			} else {
				// E' una coppia
				if (loads.find(node_label) != loads.end()){
					// Forza NON è di volume
					(*BList)[idblock[node_label]].SetMbLoad(i_step,node_label,Load*scale_f,true,loads[node_label]);
				} else (*BList)[idblock[node_label]].SetMbLoad(i_step,node_label,Load,true,0);	
			}
		ini_step = *it_step;
		i_step++;
		}
	}			
}
	
	

	


