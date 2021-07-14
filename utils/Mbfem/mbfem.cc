/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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
 * Author: Michele Frumusa <mfrumusa@gmail.com>
 */

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cstdio>
#include <iostream>
#include <sstream>
#include <map>
#include <list>

#include "myassert.h"

#include "block.h"
#include "mbdyn_out.h"
#include "nas_read.h"

using namespace std;

// This is the main file. The main duty of this source is dealing with the input file provided 
// by the user and to call the creation of the map, the interpolation and the output writing.

int
main(int argc, char* ingresso[])
{
	
	// File names of the various input and output	
	char fileFEM[80],fileMB[80],fileMat[100],fileNas[100],buffer[200],fileMass[100];
	FILE* fin;
	FILE* fmat;
	FILE* fnas;

	bool mass_file, mass_treat = false, loads;
	// Step
	int n_step;
	vector<int> step_list;
	// Block
	int N_block;
	double scale_f;
	map<int,Block> BlockList;
	map<int,int> iminMb,imaxMb,iminFem,imaxFem,mbidblock,iminFem_l,imaxFem_l;
	vector<int> femidblock;
	map<int,int> csloads;
	// Output
	bool outnas = false, outdisp = false, outtext = false, outrdir = false;
	int CID, start_lc;

	// Opening input file
	fin = fopen(ingresso[1],"r");
	if (fin!=NULL){
	// Model files
		if ( !fscanf(fin,"%s ",fileMB) )
		{
			cerr << "Unable to read input file - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		if ( fgets(buffer,200,fin) == NULL) {
			cerr << "Unable to read line - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		istringstream issfilefem(buffer,istringstream::in);
		if (issfilefem >> fileFEM >> fileMass) {
			mass_file = true;
		} else {
			mass_file = false;
		}
	// Scale factor - apply to mb 
		if ( !fscanf(fin,"%le ",&scale_f) )
		{
			cerr << "Unable to read scale factor - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
	// Step to analize
		if ( fgets(buffer,100,fin) == NULL) {
			cerr << "Unable to read line - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}

		istringstream isstep(buffer,istringstream::in);
		isstep >> n_step;
		char lst_type;		
		isstep >> lst_type;
		switch (lst_type) {
		case 's': {
			// Step on regular interval
			int st_step,delta_step;			
			isstep >> st_step >> delta_step;
			for (int i_step = 0; i_step < n_step; i_step++) step_list.push_back(st_step+delta_step*i_step);
			break;
			}
		case 'e': {
			// List of step
			for (int i_step = 0; i_step < n_step; i_step++){
				int temp;
				isstep >> temp;
				step_list.push_back(temp);
			}
			break;
			}
		default:
			cerr << "Invalid Step Keyword - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			break;
		}
	// Trattamento delle forze
		if ( !fscanf(fin,"%s ",buffer) )
		{
			cerr << "Unable to read loads flag - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		if (strcmp(buffer,"loads")==0){
			loads = true;
		} else if (strcmp(buffer,"noloads")==0) {
			loads = false;
		}
		if (loads) {
		for (int i_load = 0; i_load < 2; i_load++){
			if ( fgets(buffer,100,fin) == NULL) {
				cerr << "Unable to read line - Aborting" << endl;
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			istringstream issload(buffer,istringstream::in);
			char load_type;
			issload >> load_type;
			switch (load_type) {
				case 'c':{	
					int n_l,temp;
					issload >> n_l; 
					for (int i_l = 0;i_l < n_l;i_l++) {
						issload >> temp;
						csloads[temp] = 1;
					}
					break;
					}
				case 's':{
					int n_l,temp;
					issload >> n_l;
					for (int i_l = 0;i_l < n_l;i_l++) {
						csloads[temp] = 2;
					}
					break;
					}
				default: 
					cerr << "Invalid Load Keyword - Aborting" << endl;
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					break;
			}
		}
		}
	// Blocks
		if ( !fscanf(fin,"%d",&N_block) ) {
			cerr << "Unable to read blocks number - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		// Inizio lettura dei blocchi
		for (int i_bl = 0; i_bl < N_block; i_bl++){
			int id_bl,rbf_ord,near_nodes,n_adj;
			double ref_length;
			char linquad,beam;
			bool lq,bonb;
			// Lettura del tipo di mappa
			if ( !fscanf(fin,"%d ",&id_bl) || !fscanf(fin,"%c %d %d %d %le %c ",&linquad,&rbf_ord,&near_nodes,&n_adj,&ref_length,&beam) ) {
				cerr << "Unable to read from file " << ingresso[1] 
					<< "- Aborting" << endl;
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			switch (linquad) {
				case 'q':
					lq = true;
					break;
				case 'l':
					lq = false;
					break;
				default:
					cerr << "Invalid Linear/Quadratic Keyword - Aborting" << endl;
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					break;
			}
			switch (beam) {
				case 'y':
					// Trattamento diretto delle accelerazioni di inerzia
					bonb = true;
					if (!mass_treat){
						mass_treat = true;
					}
					break;
				case 'n':
					// No accelerazioni di inerzia
					bonb = false;
					break;
				default:
					cerr << "Invalid Rotational Inertia Keyword - Aborting" << endl;
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					break;
			}		
			// Reading block limit
			int fem_l,fem_u,mb_l,mb_u;
			if ( !fscanf(fin,"%d %d ",&iminMb[id_bl],&imaxMb[id_bl]) ||
			     !fscanf(fin,"%d %d ",&iminFem[id_bl],&imaxFem[id_bl]) ||
			     !fscanf(fin,"%d %d ",&iminFem_l[id_bl],&imaxFem_l[id_bl])
			   ) {
				cerr << "Unable to read from file " << ingresso[1] << endl;
				throw ErrGeneric(MBDYN_EXCEPT_ARGS);
			}
			// Creating block
			BlockList[id_bl] = Block(id_bl,n_adj,ref_length*scale_f,n_step,bonb);
			BlockList[id_bl].SetInterpMethod(near_nodes,rbf_ord,lq);
			cout << "Block " << BlockList[id_bl].GetId() << " created" << endl;					  
		}
	// Output Nastran
	char nas_out;
	if ( !fscanf(fin,"%c %d %d ",&nas_out,&CID,&start_lc) ) {
		cerr << "Unable to read from file " << ingresso[1] << endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	if (nas_out=='n'){
		outnas = true;
	} else {
		cerr << "Invalid Output Keyword - Aborting" << endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
	}
	// Output
		if ( fgets(buffer,100,fin) == NULL ) {
			cerr << "Unable to read line - Aborting" << endl;
			throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		}
		int n_out = strlen(buffer)-1;
		for (int i_out = 0; i_out < n_out; i_out++){
			switch (buffer[i_out]) {
				case 'd':
					// Display Result on video
					outdisp = true;
					break;
				case 't':
					// Write a text file for Matlab
					outtext = true;
					break;
				case 'v':
					// Displacement test
					outrdir = true;
					break;
				default:
					cerr << "Invalid Output Keyword - Aborting" << endl;
					throw ErrGeneric(MBDYN_EXCEPT_ARGS);
					break;
			}
		}	  
	} else {
		cerr << "Invalid path or filename - Aborting" << endl;
		throw ErrGeneric(MBDYN_EXCEPT_ARGS);
		return(0);
	}
	// Letto il file, passiamo alla lettura dei file
	cout << "Reading MBDyn file" << endl;
	mbidblock = read_mov(fileMB,iminMb,imaxMb,step_list,&BlockList,outrdir,scale_f);
	if (loads){
		read_frc(fileMB,mbidblock,step_list,csloads,&BlockList,scale_f);
	}
	cout << "Reading op2 file" << endl;
	readop2(fileFEM,iminFem,imaxFem,iminFem_l,imaxFem_l,&BlockList,&femidblock);
	// Lettura eventuale del file contenente la massa lumpata
	if (mass_file&&mass_treat){
		readmat(fileMass,femidblock,&BlockList);
	}
	// Calcolo della Mappa e dell'interpolazione
	for (map<int,Block>::iterator it_bl = BlockList.begin(); it_bl != BlockList.end(); it_bl++){
		it_bl->second.CreateMap();
		it_bl->second.InterpAcc();
		if (loads){
			it_bl->second.InterpCLoads();
		}
		if (outrdir){
			it_bl->second.InterpDisp();
		}
	}
	if (outdisp) {
		for (map<int,Block>::iterator it_bl = BlockList.begin(); it_bl != BlockList.end(); it_bl++){
			it_bl->second.Summary();
			//it_bl->second.ShowAcc();
		}
	}
	if (outnas) {
		// Nastran output
		// Apertura del file
		fnas = open_nasout(ingresso[1],n_step,step_list);
		nas_cards(fnas,BlockList,n_step,CID,start_lc);
	}
	if (outtext) {
		text_dump(ingresso[1],n_step,N_block,BlockList);
	}
	if (outrdir) {
		rdir_dump(ingresso[1],n_step,N_block,BlockList);
	}
}
			
		
	
	

 
