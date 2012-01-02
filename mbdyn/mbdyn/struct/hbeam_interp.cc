/* $Header$ */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2012
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 * Teodoro Merlini  <merlini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

#include "matvecexp.h"
#include "Rot.hh"
#include "hbeam_interp.h"

// Input:
// node_pos	array delle posizioni attuali dei due nodi piu' offset
// node_or	array delle orientazioni attuali dei due nodi
// node_f	array degli offset attuali (R*f_tilde) dei due nodi
// xi		posizione del punto, da 0 a 1
// dexi_des	derivata della posizione del punto rispetto all'ascissa
//		curvilinea
// 
// Output:
// pos		posizione del punto interpolato
// orient		orientazione del punto interpolato
// or_delta_w_or	or_delta in funzione di or_delta nodali
// delta_pos_w_or	delta_pos in funzione di or_delta nodali
// delta_pos_w_pos	delta_pos in funzione di delta_pos nodali
// F		derivata della posizione rispetto all'ascissa curvilinea
// om		assiale di d(orient)/d(s)*or.Transpose()
// delta_om_ws_or	delta_om in funzione di or_delta nodali
// delta_F_ws_or	delta_F in funzione di or_delta nodali
// delta_F_ws_pos	delta_F in funzione di delta_pos nodali
// 
// Nota:
// Dopo l'uscita bisogna convertire gli or_delta nodali 
// in delta_parametri_or

/*
 * Enum puramente locale
 */
enum {
	NODE1 = 0,
	NODE2 = 1,
	NUMNODES = 2	/* serve eventualmente per dimensionare arrays */
};


inline void Compute(const Vec3 *const node_pos,
		const Mat3x3 *const node_or,
		const Vec3 *const node_f,
		const doublereal xi,
		const doublereal dexi_des,
		Vec3 &pos,
		Mat3x3 &orient,
		Vec3 &F,
		Vec3 &om,
		MatExp& A1,
		MatExp& A2,
		MatExp& A12,
		MatExp& A_tilde,
		MatExp& Theta_r_I,
		MatExp& Theta_tilde,
		MatExp& A_xi,
		VecExp& eta_r,
		VecExp& kappa,
		VecExp& eta_tilde
		) {
	A1 = MatExp(node_or[NODE1], node_pos[NODE1].Cross(node_or[NODE1]));
	A2 = MatExp(node_or[NODE2], node_pos[NODE2].Cross(node_or[NODE2]));
	A12 = A2*(A1.Transpose());
	//VecExp e1(0.);
	eta_r = RoTrManip::Helix(A12);
	Theta_r_I = RoTrManip::DRoTr_I(eta_r);
	eta_tilde = eta_r*xi;
	RoTrManip::RoTrAndDRoTr(eta_tilde,A_tilde,Theta_tilde);
	A_xi = A_tilde*A1;
		//extract interpolated orientation and position
		orient = A_xi.GetVec();
		pos = (A_xi.GetMom()*(orient.Transpose())).Ax();

	kappa = Theta_tilde*eta_r*dexi_des;
		//extract interpolated curvature and def. gradient
		om = kappa.GetVec();
		F = kappa.GetMom();
		F -= pos.Cross(om); /* :) */
};

void ComputeInterpolation(const Vec3 *const node_pos,
			const Mat3x3 *const node_or,
			const Vec3 *const node_f,
			const doublereal xi,
			const doublereal dexi_des,
			Vec3 &pos,
			Mat3x3 &orient,
			Vec3 &F,
			Vec3 &om) {
	MatExp A1, A2, A12, A_tilde, Theta_r_I, Theta_tilde, A_xi;
	VecExp eta_r, kappa, eta_tilde;
	Compute(node_pos,node_or, node_f,xi,dexi_des,pos,orient,F,om,
		A1,A2,A12,A_tilde,Theta_r_I,Theta_tilde,A_xi,eta_r, 
		kappa,eta_tilde);
};

void ComputeFullInterpolation(const Vec3 *const node_pos,
			const Mat3x3 *const node_or,
			const Vec3 *const node_f,
			const doublereal xi,
			const doublereal dexi_des,
			Vec3 &pos,
			Mat3x3 &orient,
			Mat3x3 *const or_delta_w_or,
			Mat3x3 *const delta_pos_w_or,
			Mat3x3 *const delta_pos_w_pos,
			Vec3 &F,
			Vec3 &om,
			Mat3x3 *const delta_om_ws_or,
			Mat3x3 *const delta_F_ws_or,
			Mat3x3 *const delta_F_ws_pos) {
	MatExp A1, A2, A12, A_tilde, Theta_r_I, Theta_tilde, A_xi;
	VecExp eta_r, kappa, eta_tilde;
	Compute(node_pos,node_or, node_f,xi,dexi_des,pos,orient,F,om,
		A1,A2,A12,A_tilde,Theta_r_I,Theta_tilde,A_xi,eta_r, 
		kappa,eta_tilde);

	MatExp eta_xi_delta_2(Theta_tilde*Theta_r_I*xi);
	MatExp eta_xi_delta_1(A_tilde);
		eta_xi_delta_1-=Theta_tilde*Theta_r_I*A2*xi;
	//compute traslation wrenches at xi and at the nodes:
	MatExp Trasl_n[NUMNODES];
	Trasl_n[NODE1] = MatExp(Eye3, Mat3x3(MatCross, node_pos[NODE1]));
	Trasl_n[NODE2] = MatExp(Eye3, Mat3x3(MatCross, node_pos[NODE2]));
	MatExp Trasl_xi_bak(Eye3, Mat3x3(MatCross, -pos));
	//compute delta_kappa_i
	MatExp Lambda(RoTrManip::Elle(eta_tilde,eta_r*dexi_des));
	MatExp delta_kappa_2(Lambda*xi);
		delta_kappa_2+=Theta_tilde*dexi_des;
		delta_kappa_2 = delta_kappa_2*Theta_r_I;
	MatExp delta_kappa_1(delta_kappa_2*A2*-1);
	
	//transform eta_xi_delta_i to deltapos_ordelta_xi_i
	MatExp dpos_ord_xi_1(Trasl_xi_bak*eta_xi_delta_1);
	MatExp dpos_ord_xi_2(Trasl_xi_bak*eta_xi_delta_2);
	//project to nodes delta_pos, or_delta:
	MatExp dpos_ord_dpos_ord_1(dpos_ord_xi_1*Trasl_n[NODE1]);
	MatExp dpos_ord_dpos_ord_2(dpos_ord_xi_2*Trasl_n[NODE2]);
	//extract or_delta_w_or, delta_pos_w_or and delta_pos_w_pos
	or_delta_w_or[NODE1] = dpos_ord_dpos_ord_1.GetVec();
	or_delta_w_or[NODE2] = dpos_ord_dpos_ord_2.GetVec();
	delta_pos_w_or[NODE1] = dpos_ord_dpos_ord_1.GetMom();
	delta_pos_w_or[NODE2] = dpos_ord_dpos_ord_2.GetMom();
	delta_pos_w_pos[NODE1] = dpos_ord_dpos_ord_1.GetVec();
	delta_pos_w_pos[NODE2] = dpos_ord_dpos_ord_2.GetVec();
	//project delta_kappa_1 to nodes delta_pos, or_delta:
	MatExp dkappa_dpos_ord_1(delta_kappa_1*Trasl_n[NODE1]);
	MatExp dkappa_dpos_ord_2(delta_kappa_2*Trasl_n[NODE2]);
	//extract delta_om_ws_or, delta_F_ws_or and delta_F_ws_pos
	delta_om_ws_or[NODE1] = dkappa_dpos_ord_1.GetVec();
	delta_om_ws_or[NODE2] = dkappa_dpos_ord_2.GetVec();
	//delta kappa.Mom()
	delta_F_ws_or[NODE1] = dkappa_dpos_ord_1.GetMom();
	delta_F_ws_pos[NODE1] = dkappa_dpos_ord_1.GetVec();
	delta_F_ws_or[NODE2] = dkappa_dpos_ord_2.GetMom();
	delta_F_ws_pos[NODE2] = dkappa_dpos_ord_2.GetVec();
	//x cross delta kappa.Vec()
	delta_F_ws_or[NODE1] -= pos.Cross(dkappa_dpos_ord_1.GetVec());
	delta_F_ws_or[NODE2] -= pos.Cross(dkappa_dpos_ord_2.GetVec());
	//kappa.Vec() x delta_pos
	delta_F_ws_or[NODE1] += kappa.GetVec().Cross(delta_pos_w_or[NODE1]);
	delta_F_ws_pos[NODE1] += kappa.GetVec().Cross(delta_pos_w_pos[NODE1]);
	delta_F_ws_or[NODE2] += kappa.GetVec().Cross(delta_pos_w_or[NODE2]);
	delta_F_ws_pos[NODE2] += kappa.GetVec().Cross(delta_pos_w_pos[NODE2]);

	//tarocco per l'aggiunta dell'offset f
	//le curvature e le posizioni invece sono giuste perche' 
	//viene passato node_pos=x_node+node_or*f
	Mat3x3 fodo;
	
	fodo = node_f[NODE1].Cross(or_delta_w_or[NODE1]);
	delta_pos_w_or[NODE1] -= delta_pos_w_pos[NODE1]*fodo;
	delta_F_ws_or[NODE1] -=	delta_F_ws_pos[NODE1]*fodo;
	
	fodo = node_f[NODE2].Cross(or_delta_w_or[NODE2]);
	delta_pos_w_or[NODE2] -= delta_pos_w_pos[NODE2]*fodo;
	delta_F_ws_or[NODE2] -=	delta_F_ws_pos[NODE2]*fodo;
}

