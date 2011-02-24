#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include "block.h"

using namespace std;

// Costruttore vuoto
Block::Block(): BlockId(0),Nmb(0),MbNodes() {}

// Costruttore
Block::Block(int id,int N,double dl,int nstep,bool bonb): BlockId(id),Nmb(0),Nfem(0),Nfem_l(0),Nadj(N),Nstep(nstep),dL(dl),beam(bonb),MbNodes(),FemNodes(),FemMM(0),idnode(),MbAcc(),pIntMap(),H() {}

// Distruttore vuoto
Block::~Block() {}

int Block::GetId() {
	return BlockId;
}
	
int Block::GetNmb() {
	return Nmb;
}

int Block::GetNfem() {
	return Nfem;
}

bool Block::BeOnBe() {
	return beam;
}

bool Block::isMM() {
	if (FemMM != 0) {
		return true;
	} else {
		return false;
	}
}
	
void Block::SetInterpMethod(int Nnodes, int order, bool linquad){
	// Check sul numero di nodi vicini
	int newNnodes = Nnodes;
	if (Nadj != 0){
		int rem = Nnodes%(Nadj+1);
		if (rem != 0){
			// i nodi vicini non sono multipli esatti
			int div = int(Nnodes/(Nadj+1));
			// Aggiornamento del numero di nodi
			newNnodes = (Nadj+1)*(div+1);
		}
	}
	pIntMap = new MLSP(newNnodes,order,linquad,Nadj);
	//(*pIntMap).ShowMe();
}
	
void Block::SetNode(int label, const Vec3& P,bool ismb,bool isCnode){
	if (ismb) {
		MbNodes.push_back(Node(label,P,Nadj,dL));
		// Assegnazione alla mappa nodo-id
		int dim = idnode.size();
		idnode[label]=dim;
		Nmb++;
	} else {
		if (!isCnode) {
			FemNodes.push_back(Node(label,P));
			Nfem++;
		} else {
			FemNodes_l.push_back(Node(label,P));
			Nfem_l++;
		}
	}
}

void Block::ShowNodes(){
	for (vector<Node>::iterator it_n = MbNodes.begin(); it_n != MbNodes.end(); it_n++){
		int label = it_n->GetLabel();
		Vec3 pos = it_n->GetPos();
		cout << label << " : " << pos(1) << "," << pos(2) << "," << pos(3) << endl;
	} 
}

void Block::SetMbAcc(int lc,int label, const Vec3& Acc, const Vec3& AngV,const Vec3& AngA){
	// Ridimensiona la matrice del carico
	std::map<int, FullMatrixHandler*>::iterator m = MbAcc.find(lc);
	if (m == MbAcc.end()) {
		std::pair<std::map<int, FullMatrixHandler*>::iterator, bool> p;
		if (beam){
			p = MbAcc.insert(std::map<int, FullMatrixHandler*>::value_type(lc, new FullMatrixHandler(Nmb*(Nadj+1),6)));
		} else {
			p = MbAcc.insert(std::map<int, FullMatrixHandler*>::value_type(lc, new FullMatrixHandler(Nmb*(Nadj+1),3)));
		}
		if (!p.second) {
			
		}
		m = p.first;	
	} else if (m->second->iGetNumRows() == 0) {
		if (beam){
			m->second->ResizeReset(Nmb*(Nadj+1),6);
		} else {
			m->second->ResizeReset(Nmb*(Nadj+1),3);
		}
	}
	// Assegna i valori del vettore
	int nr = (idnode[label])*(Nadj+1)+1;
	for (int d = 1; d <= 3; d++) m->second->PutCoef(nr,d,Acc(d));
	if (beam) {
		for (int d = 1; d <= 3; d++) m->second->PutCoef(nr,d+3,AngA(d));
	}
	// e poi calcola gli adjoint
	vector<Vec3> delta_acc = MbNodes[idnode[label]].AdjAcc(AngV,AngA);
	for (int i_adj = 0; i_adj < Nadj; i_adj++){
		Vec3 temp = Acc+delta_acc[i_adj];
		for (int d = 1; d <= 3; d++) m->second->PutCoef(nr+i_adj+1,d,temp(d));
		if (beam) {
			for (int d = 1; d <= 3; d++) m->second->PutCoef(nr+i_adj+1,d+3,AngA(d));
		}
	}	
}

void Block::SetMbDisp(int lc,int label, const Vec3& Pos, const Vec3& Rot){
	// Ridimensiona la matrice del carico
	std::map<int, FullMatrixHandler*>::iterator m = MbDisp.find(lc);
	if (m == MbDisp.end()) {
		std::pair<std::map<int, FullMatrixHandler*>::iterator, bool> p;
		if (beam) {
			p = MbDisp.insert(std::map<int, FullMatrixHandler*>::value_type(lc, new FullMatrixHandler(Nmb*(Nadj+1),6)));
		} else {
			p = MbDisp.insert(std::map<int, FullMatrixHandler*>::value_type(lc, new FullMatrixHandler(Nmb*(Nadj+1),3)));
		}
		if (!p.second) {
		
		}
		m = p.first;
	} else if (m->second->iGetNumRows() == 0) {
		if (beam) {	
	       		m->second->ResizeReset(Nmb*(Nadj+1),6);
		} else {
			m->second->ResizeReset(Nmb*(Nadj+1),3);
		}
	}
	// Assegna i valori del vettore
	int nr = (idnode[label])*(Nadj+1)+1;
	for (int d = 1; d <= 3; d++) m->second->PutCoef(nr,d,Pos(d));
	if (beam) {
		for (int d = 1; d <= 3; d++) m->second->PutCoef(nr,d+3,Rot(d));
	}
	// e poi calcola gli adjoint
	vector<Vec3> delta_pos = MbNodes[idnode[label]].AdjPos(Rot);
	for (int i_adj = 0; i_adj < Nadj; i_adj++){
		Vec3 temp = Pos+delta_pos[i_adj];
		for (int d = 1; d <= 3; d++) m->second->PutCoef(nr+i_adj+1,d,temp(d));
		if (beam) {
			for (int d = 1; d <= 3; d++) m->second->PutCoef(nr+i_adj+1,d+3,Rot(d));
		}
	}
}

void Block::ShowAcc(){
	for (map<int,FullMatrixHandler*>::iterator it_lc = MbAcc.begin(); it_lc != MbAcc.end(); it_lc++){
		cout << "Load Case " << it_lc->first << endl;
		cout << "MB: " << endl;
		for (int i_n = 1; i_n <= Nmb*(Nadj+1); i_n++) {
			cout << (*it_lc->second)(i_n,1) << "," << (*it_lc->second)(i_n,2) << "," << (*it_lc->second)(i_n,3) << endl;
		}
		cout << "FEM: " << endl;
		for (int i_nf = 1; i_nf <= Nfem; i_nf++) {
			cout << (*FemAcc[it_lc->first])(i_nf,1) << "," << (*FemAcc[it_lc->first])(i_nf,2) << "," << (*FemAcc[it_lc->first])(i_nf,3) << endl;
		}
	}
}

void Block::SetMbLoad(int lc,int nlabel, const Vec3& Load, bool couple, int type){
	// Allocazione della forza
	switch (type){
		case (0):{
			// Forza di volume
			if (MbVLo[lc]->iGetNumRows()==0) MbVLo[lc]->ResizeReset(Nmb*(Nadj+1),6);
			int nr = (idnode[nlabel])*(Nadj+1)+1;
			if (couple) for (int d = 1; d <= 3; d++) MbVLo[lc]->PutCoef(nr,d+3,Load(d));
			else for (int d = 1; d <= 3; d++) MbVLo[lc]->PutCoef(nr,d,Load(d));
			break;
			}
		case (1):{
			// Carico concentrato
			// recupera la posizione associata al nodo MB
			Vec3 pos = MbNodes[idnode[nlabel]].GetPos();
			MbCLo[lc].push_back(CLoad(nlabel,Load,couple,pos));
			break;
			}
		case (2):{
			// Forza di superficie
			if (MbSLo[lc]->iGetNumRows()==0) MbSLo[lc]->ResizeReset(Nmb*(Nadj+1),6);
			int nr = (idnode[nlabel])*(Nadj+1)+1;
			if (couple) for (int d = 1; d <= 3; d++) MbSLo[lc]->PutCoef(nr,d+3,Load(d));
			else for (int d = 1; d <= 3; d++) MbSLo[lc]->PutCoef(nr,d,Load(d));
			cout << "Not implemented yet - loads will be ignored" << endl;
			break;
			}
	}		
}

void Block::SetFemMM(const Vec3& row,int id){
	if (FemMM  == 0){
		FemMM = new FullMatrixHandler(Nfem,3);
	}
	for (int i=1;i<=3;i++){
		FemMM->PutCoef(id,i,row(i));
	}
}
	

void Block::CreateMap(){
	/* resize della matrice */
	//unsigned int NonZero = (N_MBAdj <= pIM->MaxNds()) ? N_MBAdj :  pIM->MaxNds();
	if (Nfem != 0) { 
		H = new SpMapMatrixHandler(Nfem, Nmb*(Nadj+1));
		pIntMap->Interpolate(FemNodes, MbNodes, H);
	} else {
		H = new SpMapMatrixHandler(1,Nfem);
	} 
	cerr << "Block " << BlockId << ": H matrix computed" << endl;
	/* for (int i=0; i<H.iGetNumRows(); i++){
		for (int j=0; j<H.iGetNumCols(); j++){
			cout << H(i+1,j+1) << "  ";
		}
		cout << endl;
	}	*/ 
}

void Block::InterpAcc(){
	for (map<int,FullMatrixHandler*>::iterator it_st = MbAcc.begin(); it_st != MbAcc.end(); it_st++) {
		// Allocazione della memoria per la matrice
		if (FemAcc[it_st->first] == 0) {
			FemAcc[it_st->first] = new FullMatrixHandler(Nfem, it_st->second->iGetNumCols());
		} else {
			FemAcc[it_st->first]->ResizeReset(Nfem, it_st->second->iGetNumCols());
		}
		// Calcolo delle accelerazioni
		H->MatMatMul(*FemAcc[it_st->first], *it_st->second);
		if (beam && FemMM != 0) {
			// Calcolo delle coppie di inerzia
			// Allocazione dello spazio per la matrice
			if (FemIneCou[it_st->first] == 0) {
				FemIneCou[it_st->first] = new FullMatrixHandler(Nfem, 3);
			} else {
				FemIneCou[it_st->first]->ResizeReset(Nfem, 3);
			}
			// Calcolo delle coppie
			for (int i=1; i<= Nfem; i++) {
				for (int j=1; j<=3; j++) {
					FemIneCou[it_st->first]->PutCoef(i,j,(FemAcc[it_st->first]->operator()(i,j+3))*(FemMM->operator()(i,j)));
				}
			}
		}
	}	
}

void Block::InterpCLoads(){
			
	ANNkd_tree* pTree;
	int labF[Nfem_l];

	unsigned int dim = 3;                   /* dimensioni dello spazio a cui appatengono i punti */
	ANNpointArray  data = annAllocPts(Nfem_l, dim);   /* struttura che contiene le coord dei nodi */
	ANNpoint query = annAllocPt(dim);               /* query point */
	ANNidxArray idx = new ANNidx[1];
	ANNdistArray dist = new ANNdist[1];

	/* Creazione della struttura di ricerca */

	vector<Node>::iterator itFem;
	vector<Node>::const_iterator itFem_end = FemNodes_l.end();

	int n_count = 0;
	for(itFem = FemNodes_l.begin(); itFem != itFem_end; itFem++) {
		Vec3 tmp1 = itFem->GetPos();
		for (unsigned int d=0; d<dim; d++) data[n_count][d] = tmp1(d+1);
		labF[n_count] = itFem->GetLabel();
		n_count++;
	}

	/* crea la stuttura per le ricerche */
	pTree = new ANNkd_tree(data, Nfem_l , dim);   /*struttura dati per la ricerca */

	// Fa il calcolo solo sul primo load case, visto che le forze in gioco saranno sempre le stesse
	for (unsigned int it_f = 0; it_f < (MbCLo[1]).size(); it_f++){
		Vec3 tmp = MbCLo[1][it_f].GetPosMb();
		for(unsigned int d=0; d<dim; d++) query[d] = tmp(d+1);
		/* ricerca del punto più vicino */
		pTree->annkSearch(query, 1, idx, dist);
		for (map< int, vector<CLoad> >::iterator it_lc = MbCLo.begin(); it_lc != MbCLo.end(); it_lc++){
			it_lc->second[it_f].PutLabFem(labF[idx[0]]);
		}
	}
}

void Block::InterpDisp(){
	for (map<int,FullMatrixHandler*>::iterator it_st = MbDisp.begin(); it_st != MbDisp.end(); it_st++) {
		// Allocazione della memoria per la matrice
		if (FemDisp[it_st->first] == 0) {
			FemDisp[it_st->first] = new FullMatrixHandler(Nfem, it_st->second->iGetNumCols());
		} else {
			FemDisp[it_st->first]->ResizeReset(Nfem, it_st->second->iGetNumCols());
		}
		// Calcolo degli spostamenti
		H->MatMatMul(*FemDisp[it_st->first], *it_st->second);
	}
}

FullMatrixHandler* Block::GetFemAcc(int i_lc) {
	return FemAcc[i_lc];
}

Vec3 Block::GetFemAcc(int i_lc,int id){
	Vec3 temp(FemAcc[i_lc]->operator()(id,1),FemAcc[i_lc]->operator()(id,2),FemAcc[i_lc]->operator()(id,3));
	return temp;
}

Vec3 Block::GetMbAcc(int i_lc,int id_0){
	// Diciamo che è una funzione di servizio, e quindi restituisce direttamente la riga utile ignorando gli adj
	int id = (id_0-1)*(Nadj+1)+1;
	Vec3 temp(MbAcc[i_lc]->operator()(id,1),MbAcc[i_lc]->operator()(id,2),MbAcc[i_lc]->operator()(id,3));
	return temp;
}

Vec3 Block::GetFemRAcc(int i_lc,int id){
	if (beam){
		Vec3 temp(FemAcc[i_lc]->operator()(id,4),FemAcc[i_lc]->operator()(id,5),FemAcc[i_lc]->operator()(id,6));
		return temp;
	} else {
		// Return null vector
		Vec3 temp(Zero3);
		return temp;
	}
}

Vec3 Block::GetMbRAcc(int i_lc,int id_0){
        if (beam){
		int id = (id_0-1)*(Nadj+1)+1;
		Vec3 temp(MbAcc[i_lc]->operator()(id,4),MbAcc[i_lc]->operator()(id,5),MbAcc[i_lc]->operator()(id,6));
		return temp;
	} else {
		// Return null vector
		Vec3 temp(Zero3);
		return temp;
	}
}

Vec3 Block::GetMbDisp(int i_lc,int id_0){
	// Diciamo che è una funzione di servizio, e quindi restituisce direttamente la riga utile ignorando gli adj
	int id = (id_0-1)*(Nadj+1)+1;
	Vec3 temp(MbDisp[i_lc]->operator()(id,1),MbDisp[i_lc]->operator()(id,2),MbDisp[i_lc]->operator()(id,3));
	return temp;
}

Vec3 Block::GetFemDisp(int i_lc,int id){
	Vec3 temp(FemDisp[i_lc]->operator()(id,1),FemDisp[i_lc]->operator()(id,2),FemDisp[i_lc]->operator()(id,3));
	return temp;
}

Vec3 Block::GetFemRot(int i_lc,int id){
	if (beam){
		Vec3 temp(FemDisp[i_lc]->operator()(id,4),FemAcc[i_lc]->operator()(id,5),FemAcc[i_lc]->operator()(id,6));
		return temp;
	} else {
		// Return null vector
		Vec3 temp(Zero3);
		return temp;
	}
}

Vec3 Block::GetMbRot(int i_lc,int id_0){
        if (beam){
		int id = (id_0-1)*(Nadj+1)+1;
		Vec3 temp(MbDisp[i_lc]->operator()(id,4),MbAcc[i_lc]->operator()(id,5),MbAcc[i_lc]->operator()(id,6));
		return temp;
	} else {
		// Return null vector
		Vec3 temp(Zero3);
		return temp;
	}
}

FullMatrixHandler* Block::GetFemInCou(int i_lc) {
	return FemIneCou[i_lc];
}

SpMapMatrixHandler* Block::GetH() {
	return H;
}

vector<Node> Block::GetMbNodes() {
	return MbNodes;
}

vector<Node> Block::GetFemNodes() {
	return FemNodes;
}

vector<CLoad> Block::GetCLoads(int i_lc) {
	return MbCLo[i_lc];
}


void Block::Summary(){
	cout << "Block " << BlockId << endl;
	cout << "Mb Nodes: " << Nmb << " (with " << Nadj << " adjoint nodes)" << endl;
	cout << "Fem Nodes: " << Nfem << endl;
	cout << "Volume Loads: ";
	if (MbVLo[1] == 0 || MbVLo[1]->iGetNumRows()==0) {
		cout << "No" << endl;
	} else {
		cout << "Yes" << endl;
	}
	cout << "Surface Loads: ";
	if (MbSLo[1] == 0 || MbSLo[1]->iGetNumRows()==0) {
		cout << "No" << endl;
	} else {
		cout << "Yes" << endl;
	}
	cout << "Concentrated Loads: " << MbCLo[1].size() << endl;
	if (FemMM != 0){
		cout << "Rotational inertia: yes (" << FemMM->iGetNumRows() << "x" << FemMM->iGetNumCols() << ")" << endl;
	}
	cout << "Fem Nodes for concentrated loads: " << Nfem_l << endl;
}
	
				
	

		
			
	
	
