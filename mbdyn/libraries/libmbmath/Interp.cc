
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
	MatExp n1(node_or[0],Mat3x3(node_pos[0])*node_or[0]);
	MatExp n2(node_or[1],Mat3x3(node_pos[1])*node_or[1]);
	MatExp n12(n2*(n1.Transpose()));
	//VecExp e1(0.);
	VecExp e2(RoTrManip::Helix(n12));
	VecExp we2(e2*w[1]);
	MatExp E12, Theta12;
	RoTrManip::RoTrAndDRoTr(we2,E12,Theta12);
	MatExp E(E12*n1);
	or = E.GetVec();
	pos = (E.GetMom()*(or.Transpose())).Ax();
	
	
	MatExp Theta2I(RoTrManip::DRoTr_It(e2).Transpose());
	MatExp d1, d2;
	MatExp Wd1, Wd2;
	d1 = Theta12;
	d2 = Theta12*Theta2I;
	Wd1 = d1*w[0];
	Wd2 = d2*w[1];
	or_delta_w_or[0] = Wd1.GetVec();
	or_delta_w_or[1] = Wd2.GetVec();
	
	Mat3x3 posx(pos);
	
	delta_pos_w_or[0] = Wd1.GetMom();
	delta_pos_w_or[0] += Wd1.GetVec()*Mat3x3(node_pos[0]);
	delta_pos_w_or[0] -= posx*Wd1.GetVec();
	delta_pos_w_or[1] = Wd2.GetMom();
	delta_pos_w_or[1] += Wd2.GetVec()*Mat3x3(node_pos[1]);
	delta_pos_w_or[1] -= posx*Wd2.GetVec();
	
	delta_pos_w_pos[0] = Wd1.GetVec();
	delta_pos_w_pos[1] = Wd2.GetVec();
		
	VecExp kappa(Theta12*(e2*wder[1]));
	om = kappa.GetVec();
	F = kappa.GetMom();
	F -= pos.Cross(om);
	
	Wd1 = d1*wder[0];
	Wd2 = d2*wder[1];
	
	delta_om_ws_or[0] = Wd1.GetVec();
	delta_om_ws_or[1] = Wd2.GetVec();
	
	delta_F_ws_or[0] = Wd1.GetVec()*Mat3x3(node_pos[0]);
	delta_F_ws_or[1] = Wd2.GetVec()*Mat3x3(node_pos[1]);
	delta_F_ws_or[0] += Wd1.GetMom();
	delta_F_ws_or[1] += Wd2.GetMom();
	delta_F_ws_or[0] += om.Cross(delta_pos_w_or[0]);
	delta_F_ws_or[1] += om.Cross(delta_pos_w_or[1]);
	delta_F_ws_or[0] -= pos.Cross(delta_om_ws_or[0]);
	delta_F_ws_or[1] -= pos.Cross(delta_om_ws_or[1]);
	
	delta_F_ws_pos[0] = Wd1.GetVec();
	delta_F_ws_pos[1] = Wd2.GetVec();
	delta_F_ws_or[0] += om.Cross(delta_pos_w_pos[0]);
	delta_F_ws_or[1] += om.Cross(delta_pos_w_pos[1]);
	
	
};
