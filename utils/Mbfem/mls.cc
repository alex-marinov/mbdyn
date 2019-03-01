#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <cmath>

#include <myassert.h>

#include "mls.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int
__FC_DECL__(dgelsd)(int* m, int* n, int* nrhs, doublereal* a, int* lda, doublereal* b, int* ldb, double* s, double* rcond, int* rank, double* work, int* lwork, int* iwork, int* info);

#ifdef __cplusplus
}
#endif /* __cplusplus */

class Weight
{
    public:
	Weight(void) {};
	virtual void SetWeight(MyVectorHandler& r_in, SparseMatrixHandler& w_out) = 0;
};


class RBF0: public Weight
{
    public:
	
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
		for (int it = 1; it <= r.iGetSize(); it++) {
			w_out(it,it) = std::pow(1-(r(it)), 2);
		}		
  	};

};		

class RBF2: public Weight
{
    public:
	
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
		for (int it = 1; it <= r.iGetSize(); it++) {
			w_out(it,it) = std::pow(1-r(it), 4)*(4*r(it) + 1);
		}		
	};		


};		

class RBF4: public Weight
{
    public:
	
	inline void SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
		for (int it = 1; it <= r.iGetSize(); it++) {
			w_out(it,it) = std::pow(1-r(it), 6) * (35*std::pow(r(it),2) + 18*r(it) + 3);
		}		
	};

};		


class POrder
{
    public:        
	virtual void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P) = 0;
};
        
class Linear: public POrder
{
    public:
	Linear(void) {};
	    
        	/* 1 x y z */
 	void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P){
    		p(1) = 1.;
		for (unsigned int i=0; i<3; i++)
			p(i+2) = PosFem[i];
		for (unsigned int i=0; i< k; i++){ 
			P(i+1,1) = 1.;
			for (unsigned int j=0; j<3; j++) P(i+1,j+2) = PosMb[id[i]][j];
		}	
	};
};

class Quadratic: public POrder
{
    public:
        Quadratic(void) {};
    
		    /* 1 x y z x^2 xy xz y^2 yz z^2 */
	void SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P){
    		p(1) = 1.;
		for (unsigned int i=0; i<3; i++)
			p(i+2) = PosFem[i];
		for (unsigned int i=0; i<3; i++)
			p(i+5) = PosFem[i] * PosFem[0];
		for (unsigned int i=0; i<2; i++)
			p(i+8)= PosFem[i+1] * PosFem[1];
		p(10) = PosFem[2]*PosFem[2];
		for (unsigned int i=0; i<k; i++){ 
			P(i+1,1) = 1;
			for (unsigned int j=0; j<3; j++) P(i+1,j+2) = PosMb[id[i]][j];
			for (unsigned int j=0; j<3; j++) P(i+1,j+5) = PosMb[id[i]][j] * PosMb[id[i]][0];
			for (unsigned int j=0; j<2; j++) P(i+1,j+8) = PosMb[id[i]][j+1] * PosMb[id[i]][1];
			P(i+1,10) = PosMb[id[i]][2] * PosMb[id[i]][2];
		}	
	};
};

MLSP::MLSP() {}

MLSP::MLSP(unsigned int NodesN, int WgType, bool LinQuad, int Nadj)
:N(NodesN),pp(),pP(),pW(),nadj(Nadj) 
{
 	
	switch (WgType) {
		
	case 0:
		pWg = new(RBF0);
		//std::cout << "C0 RBF" << std::endl;
		break;
	case 1:
		pWg = new(RBF2);
		//std::cout << "C2 RBF" << std::endl;
		break;
	case 2:
		pWg = new(RBF4);
		//std::cout << "C4 RBF" << std::endl;
		break;
	default: 
		std::cerr << "MLSP: Unknown weigths type, switching to C0 RBF\n";
		pWg = new(RBF0);
		break;
	}
		
	if (LinQuad) {
		/* quadratic */
		if (N < 10) {
			// Innalziamo la soglia al primo multiplo esatto
			if (nadj != 0){
				if (10%(nadj+1) != 0){
					int d = int(10/(nadj+1));
					N = (nadj+1)*(d+1);
				} else {
					N = 10;
				}
			} else {
				N = 10;
			}
			std::cerr << "MLSP: for quadratic polyomial basis at least 10 interpolation points are necessary.\n"
				  << "Interpolation points number raised to " << N << std::endl;
		}
		//std::cout << "Quadratic" << std::endl;
		//std::cout << "Near Nodes: " << N << std::endl;
		pOr = new Quadratic();
		pp.ResizeReset(10);
		pP.ResizeReset(N,10);
		pW.ResizeReset(N,N);
	} else {
		/* linear */
		if (N < 4) {
			// Innalziamo la soglia al primo multiplo esatto
			if (nadj != 0){
				if (4%(nadj+1) != 0){
					int d = int(4/(nadj+1));
					N = (nadj+1)*(d+1);
			        } else {
					N = 4;
				}
			} else {
				N = 4;
			}
			std::cerr << "MLSP: for quadratic polyomial basis at least 10 interpolation points are necessary.\n"
				  << "Interpolation points number raised to " << N << std::endl;
		}
		//std::cout << "Linear" << std::endl;
		//std::cout << "Near Nodes: " << N << std::endl; 
		pOr = new Linear();
		pp.ResizeReset(4);
		pP.ResizeReset(N,4);
		pW.ResizeReset(N,N);		/* calcola le matrici A e B */
	}
	
}

void MLSP::ShowMe(){
	std::cout << N << std::endl;
}

void MLSP::Interpolate(std::vector<Node>& FemNodes, std::vector<Node>& MbNodes, SpMapMatrixHandler* H) {

	unsigned int k,MbN;
	
	unsigned int dim = 3; 			/* dimensioni dello spazio a cui appatengono i punti */
	
	unsigned int MbNadj = H->iGetNumCols();           /* il numero di nodi strutturali e' pari alle colonne di H */
	
	ANNpointArray data_complete;	/* struttura che contiene le coord dei nodi */
	ANNidxArray idx_complete;
	ANNdistArray dist_complete;
	
	if (nadj != 0){
		k = N/(nadj+1);
		MbN = MbNadj/(nadj+1);
		data_complete = annAllocPts(MbNadj, dim);   /* struttura che contiene le coord dei nodi */
		idx_complete = new ANNidx[N];
		dist_complete = new ANNdist[N];
	} else {
		k = N;
		MbN = MbNadj;
	}

	ANNpointArray data = annAllocPts(MbN, dim);   /* struttura che contiene le coord dei nodi */
	ANNpoint query = annAllocPt(dim);               /* query point */

	ANNidxArray idx = new ANNidx[k];
	ANNdistArray dist = new ANNdist[k];

	/* inserimento coordinate di NodeList in data */
	 
	std::vector<Node>::iterator itFem, itMb;
	std::vector<Node>::const_iterator itMb_end = MbNodes.end();
		
	unsigned int n_count = 0, n_count_c = 0;
	for(itMb = MbNodes.begin(); itMb != itMb_end; itMb++) {
		Vec3 tmp1 = itMb->GetPos();
		for (unsigned int d=0; d<dim; d++){
			data[n_count][d] = tmp1(d+1);
			if (nadj != 0){
				data_complete[n_count_c][d] = tmp1(d+1);
			}
		}
		n_count++;
		n_count_c++;
		if (nadj != 0){
			std::vector<Vec3> tmp2 = itMb->GetAdjPos();
			for (std::vector<Vec3>::iterator it_adj = tmp2.begin(); it_adj != tmp2.end(); it_adj++){
				for(unsigned int d=0; d<dim; d++) data_complete[n_count_c][d] = it_adj->dGet(d+1);
				n_count_c++;
			}
		}
	}

	/*if (MbN != n_count) {
		std::cerr << "Matrix H dimension is not compatible with "
			<< "node list size !!\n";
		//throw InterpolErr(); 	
	}*/
	
	pTree = new ANNkd_tree(data, MbN, dim);   /*struttura dati per la ricerca */
	
	unsigned int iFemCount = 1; 
	std::vector<Node>::iterator itFem_begin = FemNodes.begin();	
	std::vector<Node>::const_iterator itFem_end = FemNodes.end();	
	/* ciclo principale per ogni nodo aerodinamico */
	for(itFem = itFem_begin; itFem != itFem_end; itFem++) {	  
		Vec3 tmp = itFem->GetPos();
		for(unsigned int d=0; d<dim; d++) query[d] = tmp(d+1);
		/* ricerca dei k punti piu' vicini */
		pTree->annkSearch(query, k, idx, dist);
		// Ciclo di ricerca del numero di nodi presenti nel raggio
		// Associamo ai nodi individuati i nodi adjoint associati
		if (nadj != 0){
			for (unsigned i_node = 0; i_node < k ; i_node++){
				idx_complete[i_node*(nadj+1)] = idx[i_node]*(nadj+1);
				dist_complete[i_node*(nadj+1)] = dist[i_node];
				// Recupero il nodo associato all'id
				std::vector<Vec3> tmp2 = MbNodes[i_node].GetAdjPos();
				int count_a = 1;
				for (std::vector<Vec3>::iterator it_adj = tmp2.begin(); it_adj != tmp2.end(); it_adj++){
					ANNpoint adj = annAllocPt(3);
					for(unsigned int d=0; d<dim; d++) adj[d] = it_adj->dGet(d+1);
					// Calcolo delle distanze
					ANNdist dist_adj = annDist(dim,query,adj);
					// Assegnazione nel vettore completo
					idx_complete[i_node*(nadj+1)+count_a] = idx[i_node]*(nadj+1)+count_a;
					dist_complete[i_node*(nadj+1)+count_a] = dist_adj;
					count_a++;
				}
				//getchar();
			}
			pOr->SetP(query,data_complete,idx_complete, N, pp, pP);
		} else {
			/* inserimento coefficienti nel vettore p e nella matrice P */
			pOr->SetP(query, data, idx, N, pp, pP);
			/* calcolo del vettore delle distanze radiali */
		}
		MyVectorHandler r(N);
		if (nadj != 0){
			double dm = sqrt(dist_complete[N-1]);
			for(unsigned j=0; j<N; j++) r(j+1) = sqrt(dist_complete[j])/(1.05*dm);
		} else {
			double dm = sqrt(dist[N-1]);
			for(unsigned j=0; j<N; j++) r(j+1) = sqrt(dist[j])/(1.05*dm);
		}
		/* calcolo dei pesi con le funzioni RBF */			
		pWg->SetWeight(r,pW);
		/* calcola le matrici A e B */
		int d2 = pP.iGetNumCols();
		FullMatrixHandler A(d2,d2),B(d2,N);;
		pP.MatTMatMul(B,pW);
		B.MatMatMul(A,pP);
		// Calcolo della pseudo inversa
		int ord = N;
		int info;
		double* S = new double[d2];
		double rcond = -1.;
		int rank;
		int nvl = int(2*log(d2/26))+1;
		int lwork =12*d2 + 2*d2*25 + 8*d2*nvl + d2*ord + 676;
		double* work = new double[lwork]; 		
		int* iwork = new int[(3 * d2 * nvl + 11 * d2)];
		// Decomposizione SVD
		__FC_DECL__(dgelsd)(&d2, &d2, &ord, A.pdGetMat(), 
		&d2, B.pdGetMat(), &d2, S, &rcond, &rank, work, &lwork, iwork, &info);
		delete [] S;
		delete [] work;
		delete [] iwork;
		MyVectorHandler h(N);
		B.MatTVecMul(h,pp);
		if (nadj != 0){
			for (unsigned i=0;i<N;i++) H->PutCoef(iFemCount, idx_complete[i]+1,h(i+1));
		} else {
			for (unsigned i=0;i<N;i++) H->PutCoef(iFemCount, idx[i]+1,h(i+1));
		}
		iFemCount++; 
	}
}
