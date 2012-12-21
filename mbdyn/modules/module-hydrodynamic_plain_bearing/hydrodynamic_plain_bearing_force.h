#ifndef ___HYDRODYNAMIC_PLAIN_BEARING_H__INCLUDED___
#define ___HYDRODYNAMIC_PLAIN_BEARING_H__INCLUDED___

#include <ac/f2c.h>

#ifdef __cplusplus
extern "C" 
{
#endif
	
static const int NBDIRSMAX = 6;
// computes the hydrodynamic bearing force at the bearing
	
//~ SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE(b,d,Psi,eta,omega,e,e_dot,k,eps,eps_dot,delta,SoD,SoV,my,beta)
//~ DOUBLE PRECISION,INTENT(IN) :: b, d, Psi, eta, omega(2), e(2), e_dot(2)
//~ DOUBLE PRECISION,INTENT(OUT) :: k(3), eps, eps_dot, delta, SoD, SoV, my, beta
void __FC_DECL__(hydrodynamic_plain_bearing_force)(
					const doublereal& b,
					const doublereal& d, 
					const doublereal& Psi, 
					const doublereal& eta, 
					const doublereal omega[2], 
					const doublereal e[2], 
					const doublereal e_dot[2], 
					doublereal k[3],	
					doublereal &eps,
					doublereal &eps_dot,
					doublereal& delta,
					doublereal& SoD,
					doublereal& SoV,
					doublereal& my,
					doublereal& beta
					);
					
// computes the hydrodynamic stiffness matrix and damping matrix at the bearing side

//!  Differentiation of hydrodynamic_plain_bearing_force in forward (tangent) mode: (multi-directional mode)
//!   variations   of useful results: k
//!   with respect to varying inputs: e omega e_dot
//SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV(b, d, psi, eta, omega, omegad, e, ed, e_dot, e_dotd, k, kd, eps, eps_dot, delta, sod, sov, my, beta, nbdirs)

//  DOUBLE PRECISION, INTENT(IN) :: b, d, psi, eta, omega(2), e(2), e_dot(2)
//  DOUBLE PRECISION, INTENT(IN) :: omegad(nbdirsmax, 2), ed(nbdirsmax, 2), e_dotd(nbdirsmax, 2)
//  DOUBLE PRECISION, INTENT(OUT) :: k(3), eps, eps_dot, delta, sod, sov, my, beta
//  DOUBLE PRECISION :: kd(nbdirsmax, 3)
//	INTEGER :: nbdirs
void __FC_DECL__(hydrodynamic_plain_bearing_force_dv)(
						const doublereal& b,
						const doublereal& d, 
						const doublereal& Psi, 
						const doublereal& eta, 
						const doublereal omega[2],
						const doublereal omegad[2][NBDIRSMAX],
						const doublereal e[2], 
						const doublereal ed[2][NBDIRSMAX],
						const doublereal e_dot[2], 
						const doublereal e_dotd[2][NBDIRSMAX],
						doublereal k[3],
						doublereal kd[3][NBDIRSMAX], 
						doublereal& eps, 
						doublereal& eps_dot, 
						doublereal& delta, 
						doublereal& SoD, 
						doublereal& SoV, 
						doublereal& my, 
						doublereal& beta,
						const int& nbdirs);


// computes the contact force between shaft and bearing on the bearing side	


//~ SUBROUTINE PLAIN_BEARING_CONTACT_FORCE( d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot, e, e_dot, k)
//~ DOUBLE PRECISION, INTENT(IN) :: d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot(2), e(2), e_dot(2)
//~ DOUBLE PRECISION, INTENT(OUT) :: k(3)
void __FC_DECL__(plain_bearing_contact_force)(
					const doublereal& d,
					const doublereal& Psi,
					const doublereal& sP,
					const doublereal& DL,
					const doublereal& m,
					const doublereal& abs_FPrs1,
					const doublereal& myP,
					const doublereal& signum_delta_omega,
					const doublereal Phi_dot[2],
					const doublereal e[2],
					const doublereal e_dot[2],
					doublereal k[3]);

// computes the contact stiffness matrix

//SUBROUTINE PLAIN_BEARING_CONTACT_FORCE_DV(d, psi, sp, dl, m, abs_fprs1, myp, signum_delta_omega, phi_dot, phi_dotd, e, ed, e_dot, e_dotd, k, kd, nbdirs)
//  DOUBLE PRECISION, INTENT(IN) :: d, psi, sp, dl, m, abs_fprs1, myp, signum_delta_omega, phi_dot(2), e(2), e_dot(2)
//  DOUBLE PRECISION, INTENT(IN) :: phi_dotd(nbdirsmax, 2), ed(nbdirsmax, 2), e_dotd(nbdirsmax, 2)
//  DOUBLE PRECISION, INTENT(OUT) :: k(3)
//  DOUBLE PRECISION, INTENT(OUT) :: kd(nbdirsmax, 3)
//  INTEGER :: nbdirs

void __FC_DECL__(plain_bearing_contact_force_dv)(
					const doublereal& d,
					const doublereal& psi, 
					const doublereal& sp, 
					const doublereal& dl, 
					const doublereal& m, 
					const doublereal& abs_fprs1,  
					const doublereal& myp, 
					const doublereal& signum_delta_omega, 
					const doublereal phi_dot[2], 
					const doublereal phi_dotd[2][NBDIRSMAX],
					const doublereal e[2], 
					const doublereal ed[2][NBDIRSMAX], 
					const doublereal e_dot[2], 
					const doublereal e_dotd[2][NBDIRSMAX],
					doublereal k[3], 
					doublereal kd[3][NBDIRSMAX], 
					const int& nbdirs);					


#ifdef __cplusplus
} // extern "C"
#endif
												
#endif												
