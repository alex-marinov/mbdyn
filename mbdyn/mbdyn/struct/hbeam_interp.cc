/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
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
 * Copyright (C) 1996-2000
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
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
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

#include "matvecexp.h"
#include "Rot.hh"

// Input:
// node_pos	array delle posizioni attuali dei due nodi
// node_or		array delle orientazioni attuali dei due nodi
// w		array dei pesi associati ai due nodi nel punto di calcolo
// wder		array dei pesi derivati rispetto all'ascissa curvilinea
// 
// Output:
// pos		posizione del punto interpolato
// or		orinetazione del punto interpolato
// or_delta_w_or	or_delta in funzione di or_delta nodali
// delta_pos_w_or	delta_pos in funzione di or_delta nodali
// delta_pos_w_pos	delta_pos in funzione di delta_pos nodali
// F		derivata della posizione rispetto all'ascissa curvilinea
// om		assiale di d(or)/d(s)*or.Transpose()
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

void ComputeInterpolation(const Vec3 *const node_pos,
			const Mat3x3 *const node_or,
			const doublereal *const w,
			const doublereal *const wder,
			Vec3 &pos,
			Mat3x3 &or,
			Mat3x3 *const or_delta_w_or,
			Mat3x3 *const delta_pos_w_or,
			Mat3x3 *const delta_pos_w_pos,
			Vec3 &F,
			Vec3 &om,
			Mat3x3 *const delta_om_ws_or,
			Mat3x3 *const delta_F_ws_or,
			Mat3x3 *const delta_F_ws_pos) {
	MatExp n1(node_or[NODE1],Mat3x3(node_pos[NODE1])*node_or[NODE1]);
	MatExp n2(node_or[NODE2],Mat3x3(node_pos[NODE2])*node_or[NODE2]);
	MatExp n12(n2*(n1.Transpose()));
	//VecExp e1(0.);
	VecExp e2(RoTrManip::Helix(n12));
	VecExp we2(e2*w[NODE2]);
	MatExp E12, Theta12;
	RoTrManip::RoTrAndDRoTr(we2,E12,Theta12);
	MatExp E(E12*n1);
	or = E.GetVec();
	pos = (E.GetMom()*(or.Transpose())).Ax();
	
	
	MatExp Theta2I(RoTrManip::DRoTr_I(e2));
	MatExp d1, d2;
	MatExp Wd1, Wd2;
	d1 = Theta12;
	d2 = Theta12*Theta2I;
	Wd1 = d1*w[NODE1];
	Wd2 = d2*w[NODE2];
	or_delta_w_or[NODE1] = Wd1.GetVec();
	or_delta_w_or[NODE2] = Wd2.GetVec();
	
	delta_pos_w_or[NODE1] = Wd1.GetMom();
	delta_pos_w_or[NODE1] += Wd1.GetVec()*Mat3x3(node_pos[NODE1]);
	delta_pos_w_or[NODE1] -= pos.Cross(Wd1.GetVec());
	delta_pos_w_or[NODE2] = Wd2.GetMom();
	delta_pos_w_or[NODE2] += Wd2.GetVec()*Mat3x3(node_pos[NODE2]);
	delta_pos_w_or[NODE2] -= pos.Cross(Wd2.GetVec());
	
	delta_pos_w_pos[NODE1] = Wd1.GetVec();
	delta_pos_w_pos[NODE2] = Wd2.GetVec();
		
	VecExp kappa(Theta12*(e2*wder[NODE2]));
	om = kappa.GetVec();
	F = kappa.GetMom();
	F -= pos.Cross(om); /* :) */
	
	Wd1 = d1*wder[NODE1];
	Wd2 = d2*wder[NODE2];
	
	delta_om_ws_or[NODE1] = Wd1.GetVec();
	delta_om_ws_or[NODE2] = Wd2.GetVec();
	
	delta_F_ws_or[NODE1] = Wd1.GetVec()*Mat3x3(node_pos[NODE1]);
	delta_F_ws_or[NODE2] = Wd2.GetVec()*Mat3x3(node_pos[NODE2]);
	delta_F_ws_or[NODE1] += Wd1.GetMom();
	delta_F_ws_or[NODE2] += Wd2.GetMom();
	delta_F_ws_or[NODE1] += om.Cross(delta_pos_w_or[NODE1]);
	delta_F_ws_or[NODE2] += om.Cross(delta_pos_w_or[NODE2]);
	delta_F_ws_or[NODE1] -= pos.Cross(delta_om_ws_or[NODE1]);
	delta_F_ws_or[NODE2] -= pos.Cross(delta_om_ws_or[NODE2]);
	
	delta_F_ws_pos[NODE1] = Wd1.GetVec();
	delta_F_ws_pos[NODE2] = Wd2.GetVec();
	delta_F_ws_or[NODE1] += om.Cross(delta_pos_w_pos[NODE1]);
	delta_F_ws_or[NODE2] += om.Cross(delta_pos_w_pos[NODE2]);
	
	
	//hic sunt leones....
// 	Mat3x3 L(RotManip::Elle(we2.GetVec(),(e2.GetVec())*wder[NODE2]));
// 	delta_om_ws_or[NODE1] += L*Theta12.GetVec().Inv()*or_delta_w_or[NODE1];
// 	delta_om_ws_or[NODE2] += L*Theta12.GetVec().Inv()*or_delta_w_or[NODE2];
// 	
// 	L = RotManip::Elle(we2.GetVec(),(e2.GetMom())*wder[NODE2]);
// 	delta_F_ws_or[NODE1] += L*Theta12.GetVec().Inv()*or_delta_w_or[NODE1];
// 	delta_F_ws_or[NODE2] += L*Theta12.GetVec().Inv()*or_delta_w_or[NODE2];


	Mat3x3 L(RotManip::Elle(we2.GetVec(),e2.GetVec()*wder[NODE2]));
	Mat3x3 Theta12vI(RotManip::DRot_I(we2.GetVec()));
	delta_om_ws_or[NODE1] += L*Theta12vI*or_delta_w_or[NODE1];
	delta_om_ws_or[NODE2] += L*Theta12vI*or_delta_w_or[NODE2];
	
// 	delta_F_ws_or[NODE1] += L.GetMom()*Theta12vI*or_delta_w_or[NODE1];
// 	delta_F_ws_or[NODE2] += L.GetMom()*Theta12vI*or_delta_w_or[NODE2];

	Mat3x3 LT(L*Theta12.GetVec());

	Mat3x3 tmp = pos.Cross(or_delta_w_or[NODE1]);
	tmp += delta_pos_w_or[NODE1];
	tmp += delta_pos_w_pos[NODE1];
	tmp += or_delta_w_or[NODE1];
	delta_F_ws_or[NODE1] += LT*tmp;

	tmp = pos.Cross(or_delta_w_or[NODE2]);
	tmp += delta_pos_w_or[NODE2];
	tmp += delta_pos_w_pos[NODE2];
	tmp += or_delta_w_or[NODE2];
	delta_F_ws_or[NODE2] += LT*tmp;
};

