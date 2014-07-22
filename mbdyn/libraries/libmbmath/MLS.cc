/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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
 * Implemented by Michele Frumusa, probably also based on work of others
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_ANN

#include "ac/lapack.h"
#include "MLS.h"

void 
RBF0::SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
                for (int it = 1; it <= r.iGetSize(); it++) {
                	w_out(it,it) = std::pow(1-(r(it)), 2);
		}
}

void
RBF2::SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
		for (int it = 1; it <= r.iGetSize(); it++) {
		 	w_out(it,it) = std::pow(1-r(it), 4)*(4*r(it) + 1);
		}			                 
}

void 
RBF4::SetWeight(MyVectorHandler& r, SparseMatrixHandler& w_out) {
                for (int it = 1; it <= r.iGetSize(); it++) {
			w_out(it,it) = std::pow(1-r(it), 6) * (35*std::pow(r(it),2) + 18*r(it) + 3);
		}
}

void 
Linear::SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P){
		// Filling p vector
		p(1) = 1.;
		for (unsigned int i=0; i<3; i++){
			p(i+2) = PosFem[i];
		}
		// Filling P matrix
		for (unsigned int i=0; i< k; i++){
			P(i+1,1) = 1.;
			for (unsigned int j=0; j<3; j++){
				P(i+1,j+2) = PosMb[id[i]][j];
			}	
		}
}

void
Quadratic::SetP(double* PosFem, double** PosMb, int* id, unsigned int k, MyVectorHandler& p, FullMatrixHandler& P){
                // Filling p vector
		p(1) = 1.;
		for (unsigned int i=0; i<3; i++){
			p(i+2) = PosFem[i];
		}
		for (unsigned int i=0; i<3; i++){
			p(i+5) = PosFem[i] * PosFem[0];
			for (unsigned int i=0; i<2; i++){
				p(i+8)= PosFem[i+1] * PosFem[1];
			}
		}
		p(10) = PosFem[2]*PosFem[2];
		// Filling P matrix
		for (unsigned int i=0; i<k; i++){
			P(i+1,1) = 1;
			for (unsigned int j=0; j<3; j++){
				P(i+1,j+2) = PosMb[id[i]][j];
			}
			for (unsigned int j=0; j<3; j++){
				P(i+1,j+5) = PosMb[id[i]][j] * PosMb[id[i]][0];
			}
			for (unsigned int j=0; j<2; j++){
				P(i+1,j+8) = PosMb[id[i]][j+1] * PosMb[id[i]][1];
			}
			P(i+1,10) = PosMb[id[i]][2] * PosMb[id[i]][2];
		}
		
}

MLSP::MLSP(unsigned int NodesN, int WgType, bool LinQuad)
:N(NodesN), pP(), pW(), pp()
{
	// Weight Matrix
	switch (WgType) {

		        case 0:
				pWg = new(RBF0);
				std::cout << "C0 RBF" << std::endl;
				break;
			case 2:
				pWg = new(RBF2);
				std::cout << "C2 RBF" << std::endl;
				break;
			case 4:
				pWg = new(RBF4);
				std::cout << "C4 RBF" << std::endl;
				break;
	} 
	// p and P matrix
	if (LinQuad) {
		SAFENEWWITHCONSTRUCTOR(pOr , Quadratic , Quadratic() );
		/*SAFENEWWITHCONSTRUCTOR(pp , MyVectorHandler , MyVectorHandler(10) );
		SAFENEWWITHCONSTRUCTOR(pP , FullMatrixHandler , FullMatrixHandler(N,10) );
		SAFENEWWITHCONSTRUCTOR(pW , SpMapMatrixHandler , SpMapMatrixHandler(N,N) );*/
		pp.ResizeReset(10);
		pP.ResizeReset(N,10);
		pW.ResizeReset(N,N);
	} else {
		SAFENEWWITHCONSTRUCTOR(pOr , Linear , Linear() );
		/*SAFENEWWITHCONSTRUCTOR(pp , MyVectorHandler , MyVectorHandler(4) );
		SAFENEWWITHCONSTRUCTOR(pP , FullMatrixHandler , FullMatrixHandler(N,4) );
		SAFENEWWITHCONSTRUCTOR(pW , SpMapMatrixHandler , SpMapMatrixHandler(N,N) );*/
		pp.ResizeReset(4);
		pP.ResizeReset(N,4);
		pW.ResizeReset(N,N);
	}
}

void 
MLSP::Interpolate(const GeometryData& fem_data, const GeometryData& mb_data, SpMapMatrixHandler* pH) 
{
	// Dimension of space
	unsigned dim = 3;

	// Data size
	unsigned mb_size = pH->iGetNumCols();
	unsigned fem_size = pH->iGetNumRows();

	// Data for ANN search
	ANNpointArray data = annAllocPts(mb_size,dim);
	ANNpoint query = annAllocPt(dim);
	ANNidxArray idx = new ANNidx[N];
	ANNdistArray dist = new ANNdist[N];

	// Filling ANN structures with data from model
	for (unsigned i = 0 ; i < mb_size ; i++){
		for (unsigned j = 0 ; j < dim ; j++){
			data[i][j] = mb_data.data[i].X(j+1);
		}
	}

	SAFENEWWITHCONSTRUCTOR(pTree , ANNkd_tree ,  ANNkd_tree(data, mb_size , dim));

	// Starting search
	for (unsigned i = 0; i < fem_size ; i++){
		// Query
		for (unsigned j = 0 ; j < dim ; j++){
			query[j] = fem_data.data[i].X(j+1);
		}
		// Search
		pTree->annkSearch(query, N, idx, dist);
		// Setting p and P
		pOr->SetP(query, data, idx, N, pp, pP);
		// Computing distances
		MyVectorHandler r(N);
		double dm = sqrt(dist[N-1]);
		for (unsigned j = 0; j < N ; j++){
			r(j+1) = sqrt(dist[j])/(1.05*dm);
		}
		// Setting Weight
		pWg->SetWeight(r,pW);
		// Computing Matrix
		int d2 = pP.iGetNumCols();
		FullMatrixHandler A(d2,d2),B(d2,N);;
		pP.MatTMatMul(B,pW);
		B.MatMatMul(A,pP);
		int ord = N;
		int info;
		double* S = new double[d2];
		double rcond = -1.;
		int rank;
		int nvl = int(2*log(d2/26))+1;
		int lwork =12*d2 + 2*d2*25 + 8*d2*nvl + d2*ord + 676;
		double* work = new double[lwork];
		int* iwork = new int[(3 * d2 * nvl + 11 * d2)];
		 __FC_DECL__(dgelsd)(&d2, &d2, &ord, A.pdGetMat(),
		                 &d2, B.pdGetMat(), &d2, S, &rcond, &rank, work, &lwork, iwork, &info);
		delete [] S;
		delete [] work;
		delete [] iwork;
                MyVectorHandler h(N);
		B.MatTVecMul(h,pp);
		for (unsigned int j = 0 ; j < N ; j++){
			pH->operator()(i+1, idx[j]+1) = h(j+1);
		}


	}
}

void 
MLSP::Interpolate_Adj(const GeometryData& fem_data, const GeometryData& mb_data, SpMapMatrixHandler* pH,
		 const std::vector<Vec3>& Adj) 
{
	// Dimension of space
	unsigned dim = 3;

	// Data size
	unsigned mb_size = pH->iGetNumCols();
	unsigned fem_size = pH->iGetNumRows();
	unsigned Nadj = Adj.size();
	unsigned mb_base_size = mb_size/(Nadj+1);

	unsigned k = N/(Nadj+1);

	// Data for ANN search
	ANNpointArray data = annAllocPts(mb_size,dim);
	ANNpointArray data_temp = annAllocPts(mb_base_size,dim);
	ANNpoint query = annAllocPt(dim);
	ANNidxArray idx = new ANNidx[N];
	ANNidxArray idx_temp = new ANNidx[k];
	ANNdistArray dist = new ANNdist[N];
	ANNdistArray dist_temp = new ANNdist[k];

	// Filling ANN structures with data from model
	for (unsigned i = 0 ; i < mb_base_size ; i++){
		for (unsigned j = 0 ; j < dim ; j++){
			data_temp[i][j] = mb_data.data[i].X(j+1);
			data[i*(Nadj+1)][j] = mb_data.data[i].X(j+1);
			for (unsigned i_a = 0; i_a < Nadj ; i_a++){
				for (unsigned d = 0; d < dim ; d++){
					data[i*(Nadj+1)+i_a+1][d] = mb_data.data[i].X(d+1)+Adj[i_a](d+1);
				}
			}	
		}
	}

	SAFENEWWITHCONSTRUCTOR(pTree , ANNkd_tree ,  ANNkd_tree(data_temp, mb_base_size , dim));

	// Starting search
	for (unsigned i = 0; i < fem_size ; i++){
		// Query
		for (unsigned j = 0 ; j < dim ; j++){
			query[j] = fem_data.data[i].X(j+1);
		}
		// Search
		pTree->annkSearch(query, k, idx_temp, dist_temp);
		// Inflating with adjoint nodes
		for (unsigned i_node = 0; i_node < k ; i_node++) {
			idx[i_node*(Nadj+1)] = idx_temp[i_node]*(Nadj+1);
			dist[i_node*(Nadj+1)] = dist_temp[i_node];
			unsigned count_a = 1;
			for (std::vector<Vec3>::const_iterator it_adj = Adj.begin(); 
				it_adj != Adj.end(); ++it_adj)
			{
			     	ANNpoint adj = annAllocPt(dim);
				for (unsigned int d = 0 ; d < dim ; d++) {
					adj[d] = data_temp[idx_temp[i_node]][d]+it_adj->dGet(d+1);
				}
				ANNdist dist_adj = annDist(dim,query,adj);
				idx[i_node*(Nadj+1)+count_a] = idx_temp[i_node]*(Nadj+1)+count_a;
				dist[i_node*(Nadj+1)+count_a] = dist_adj;
			}
		}
		pOr->SetP(query, data, idx, N, pp, pP);
		// Computing distances
		MyVectorHandler r(N);
		double dm = sqrt(dist[N-1]);
		for (unsigned j = 0; j < N ; j++){
			r(j+1) = sqrt(dist[j])/(1.05*dm);
		}
		// Setting Weight
		pWg->SetWeight(r,pW);
		// Computing Matrix
		int d2 = pP.iGetNumCols();
		FullMatrixHandler A(d2,d2),B(d2,N);;
		pP.MatTMatMul(B,pW);
		B.MatMatMul(A,pP);
		int ord = N;
		int info;
		double* S = new double[d2];
		double rcond = -1.;
		int rank;
		int nvl = int(2*log(d2/26))+1;
		int lwork =12*d2 + 2*d2*25 + 8*d2*nvl + d2*ord + 676;
		double* work = new double[lwork];
		int* iwork = new int[(3 * d2 * nvl + 11 * d2)];
		 __FC_DECL__(dgelsd)(&d2, &d2, &ord, A.pdGetMat(),
		                 &d2, B.pdGetMat(), &d2, S, &rcond, &rank, work, &lwork, iwork, &info);
		delete [] S;
		delete [] work;
		delete [] iwork;
                MyVectorHandler h(N);
		B.MatTVecMul(h,pp);
		for (unsigned int j = 0 ; j < N ; j++){
			pH->operator()(i+1, idx[j]+1) = h(j+1);
		}
	}
}

#endif // USE_ANN
