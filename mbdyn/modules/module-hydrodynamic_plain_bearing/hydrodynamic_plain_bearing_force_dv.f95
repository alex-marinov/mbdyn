! $Header$
! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2012
! 
! Pierangelo Masarati	<masarati@aero.polimi.it>
! Paolo Mantegazza	<mantegazza@aero.polimi.it>
! 
! Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
! via La Masa, 34 - 20156 Milano, Italy
! http://www.aero.polimi.it
! 
! Changing this copyright notice is forbidden.
! 
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation (version 2 of the License).
! 
! 
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
MODULE MATH_CONST_DV
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979d0
END MODULE MATH_CONST_DV

!  Differentiation of hydrodynamic_plain_bearing_and_contact_force in forward (tangent) mode: (multi-directional mode)
!   variations   of useful results: k
!   with respect to varying inputs: e omega e_dot
!   RW status of diff variables: e:in omega:in k:out e_dot:in
SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_AND_CONTACT_FORCE_DV(b, d, psi, &
&  eta, omega, omegad, e, ed, e_dot, e_dotd, sp, dl, m, abs_fprs1, myp, &
&  signum_delta_omega, k, kd, eps, eps_dot, delta, sod, sov, my, beta, &
&  nbdirs)
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: b, d, psi, eta, omega(2), e(2), e_dot(&
&  2)
  DOUBLE PRECISION, INTENT(IN) :: omegad(nbdirsmax, 2), ed(nbdirsmax, 2)&
&  , e_dotd(nbdirsmax, 2)
  DOUBLE PRECISION, INTENT(IN) :: sp, dl, m, abs_fprs1, myp, &
&  signum_delta_omega
  DOUBLE PRECISION, INTENT(OUT) :: k(3), eps, eps_dot, delta, sod, sov, &
&  my, beta
  DOUBLE PRECISION, INTENT(OUT) :: kd(nbdirsmax, 3)
  DOUBLE PRECISION :: k_hydro(3), k_contact(3)
  DOUBLE PRECISION :: k_hydrod(nbdirsmax, 3), k_contactd(nbdirsmax, 3)
  EXTERNAL HYDRODYNAMIC_PLAIN_BEARING_FORCE, PLAIN_BEARING_CONTACT_FORCE
  EXTERNAL HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV, &
&      PLAIN_BEARING_CONTACT_FORCE_DV
  INTEGER :: nd
  INTEGER :: nbdirs
  CALL HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV(b, d, psi, eta, omega, omegad&
&                                     , e, ed, e_dot, e_dotd, k_hydro, &
&                                     k_hydrod, eps, eps_dot, delta, sod&
&                                     , sov, my, beta, nbdirs)
  CALL PLAIN_BEARING_CONTACT_FORCE_DV(d, psi, sp, dl, m, abs_fprs1, myp&
&                                , signum_delta_omega, omega, omegad, e, &
&                                ed, e_dot, e_dotd, k_contact, k_contactd&
&                                , nbdirs)
  DO nd=1,nbdirs
    kd(nd, :) = k_hydrod(nd, :) + k_contactd(nd, :)
  END DO
  k = k_hydro + k_contact
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_AND_CONTACT_FORCE_DV

!  Differentiation of hydrodynamic_plain_bearing_force in forward (tangent) mode: (multi-directional mode)
!   variations   of useful results: k
!   with respect to varying inputs: e omega e_dot
SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV(b, d, psi, eta, omega, &
&  omegad, e, ed, e_dot, e_dotd, k, kd, eps, eps_dot, delta, sod, sov, my&
&  , beta, nbdirs)
  USE MATH_CONST_DV
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!       PRINT OUT OF RESULTS
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!~ PRINT *,'HYDRODYNAMIC_PLAIN_BEARING_FORCE'
!~ PRINT *,'b=',b
!~ PRINT *,'d=',d
!~ PRINT *,'eta=',eta
!~ PRINT *,'Psi=',Psi
!~ PRINT *,'e=',e
!~ PRINT *,'e_dot=',e_dot
!~ PRINT *,'omega_w=',omega_w
!~ PRINT *,'delta=',delta
!~ PRINT *,'phi=',phi
!~ PRINT *,'kappa=',kappa
!~ PRINT *,'abs_e=',abs_e
!~ PRINT *,'abs_e_dot=',abs_e_dot
!~ PRINT *,'eps_dot=',eps_dot
!~ PRINT *,'eps=',eps
!~ PRINT *,'delta_dot=',delta_dot
!~ PRINT *,'a1=',a1
!~ PRINT *,'a2=',a2
!~ PRINT *,'a3=',a3
!~ PRINT *,'a4=',a4
!~ PRINT *,'a5=',a5
!~ PRINT *,'a6=',a6
!~ PRINT *,'a7=',a7
!~ PRINT *,'a8=',a8
!~ PRINT *,'a9=',a9
!~ PRINT *,'SoD=',SoD
!~ PRINT *,'SoV=',SoV
!~ PRINT *,'beta=',beta
!~ PRINT *,'my=',my
!~ PRINT *,'alpha=',alpha
!~ PRINT *,'omega_res=',omega_res
!~ PRINT *,'abs_FD=',abs_FD
!~ PRINT *,'abs_FV=',abs_FV
!~ PRINT *,'abs_MR=',abs_MR
!~ PRINT *, 'k=',k
!-------------------------------------------------------------------------------------------------------------------------
!       hydrodynamic plain bearing calculation according to Butenschoen's theory
!-------------------------------------------------------------------------------------------------------------------------
!       COORDINATE SYSTEM:
!       x ... axial direction
!       y, z ... radial direction
!-------------------------------------------------------------------------------------------------------------------------
!        INPUT PARAMETERS
!-------------------------------------------------------------------------------------------------------------------------
!       b   ... bearing width [m]
!       d   ... bearing diameter [m]
!       Psi ... relative radial clearance Psi = ( D - d ) / d [1]
!       eta ... dynamic oil viscosity [Pa*s]
!       omega(1) ... angular velocity of the shaft [rad/s]
!	    omega(2) ... angular velocity of the bearing [rad/s]
!       e ...  radial eccentricity of the shaft 
!               e(1) = ey  [m]
!               e(2) = ez  [m]
!       e_dot(1) ... velocity of the shaft relative to the bearing 
!               e_dot(1) = ey_dot [m/s]
!               e_dot(2) = ez_dot [m/s]
!
!--------------------------------------------------------------------------------------------------------------------------
!       OUTPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------
!       k ... force on the bearing
!               k(1) = Fy [N]
!               k(2) = Fz [N]
!               k(3) = Mx [Nm]
!       eps ...             relative eccentricity of the shaft [1]
!       eps_dot ...      time derivative of the relative eccentricity [1/s]
!       delta ...          angle of minimum clearance between shaft and bearing [rad]
!       SoD ....           Sommerfeld number for rotation [1]
!       SoV ...            Sommerfeld number for displacement [1]
!       my ...             friction coefficient [N/N]
!       beta  ...          angle between force for rotation and minimum clearance [rad]
  DOUBLE PRECISION, INTENT(IN) :: b, d, psi, eta, omega(2), e(2), e_dot(&
&  2)
  DOUBLE PRECISION, INTENT(IN) :: omegad(nbdirsmax, 2), ed(nbdirsmax, 2)&
&  , e_dotd(nbdirsmax, 2)
  DOUBLE PRECISION, INTENT(OUT) :: k(3), eps, eps_dot, delta, sod, sov, &
&  my, beta
  DOUBLE PRECISION :: kd(nbdirsmax, 3), epsd(nbdirsmax), &
&  eps_dotd(nbdirsmax), deltad(nbdirsmax), sodd(nbdirsmax), sovd(&
&  nbdirsmax), myd(nbdirsmax), betad(nbdirsmax)
  DOUBLE PRECISION :: abs_e, abs_e_dot, s
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: abs_ed, abs_e_dotd
  DOUBLE PRECISION :: delta_dot, alpha, kappa, phi
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: delta_dotd, alphad, kappad, &
&  phid
  DOUBLE PRECISION :: omega_res, abs_fd, abs_fv, abs_mr
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: omega_resd, abs_fdd, abs_fvd&
&  , abs_mrd
  DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9
  INTRINSIC DATAN2, DSQRT, SUM, DCOS, DSIN, DABS, DSIGN
  DOUBLE PRECISION, DIMENSION(2) :: arg1
  DOUBLE PRECISION, DIMENSION(nbdirsmax, 2) :: arg1d
  DOUBLE PRECISION :: arg2
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg2d
  DOUBLE PRECISION :: arg10
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg10d
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: result1d
  DOUBLE PRECISION :: pwx1
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwx1d
  DOUBLE PRECISION :: pwr1
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: pwr1d
  INTEGER :: nd
  INTEGER :: nbdirs
  DOUBLE PRECISION :: dabs1d(nbdirsmax)
  DOUBLE PRECISION :: dabs4d(nbdirsmax)
  DOUBLE PRECISION :: dabs7d(nbdirsmax)
  DOUBLE PRECISION :: dabs11d(nbdirsmax)
  DOUBLE PRECISION :: dabs3d(nbdirsmax)
  DOUBLE PRECISION :: dabs6d(nbdirsmax)
  DOUBLE PRECISION :: dabs11
  DOUBLE PRECISION :: dabs10
  DOUBLE PRECISION :: dabs9d(nbdirsmax)
  DOUBLE PRECISION :: dabs10d(nbdirsmax)
  DOUBLE PRECISION :: dabs9
  DOUBLE PRECISION :: dabs8
  DOUBLE PRECISION :: dabs7
  DOUBLE PRECISION :: dabs6
  DOUBLE PRECISION :: dabs2d(nbdirsmax)
  DOUBLE PRECISION :: dabs5
  INTRINSIC DACOS
  DOUBLE PRECISION :: dabs4
  DOUBLE PRECISION :: dabs3
  INTRINSIC ATAN
  DOUBLE PRECISION :: dabs2
  DOUBLE PRECISION :: dabs1
  DOUBLE PRECISION :: dabs8d(nbdirsmax)
! clearance of the shaft inside the bearing
  s = d*psi/(1-psi)
! absolute value of the eccentricity of the shaft inside the bearing
  arg1(:) = e(:)**2
  arg2 = SUM(arg1(:))
  DO nd=1,nbdirs
    deltad(nd) = (ed(nd, 2)*e(1)-ed(nd, 1)*e(2))/(e(2)**2+e(1)**2)
    phid(nd) = (e_dotd(nd, 2)*e_dot(1)-e_dotd(nd, 1)*e_dot(2))/(e_dot(2)&
&      **2+e_dot(1)**2)
    kappad(nd) = phid(nd) - deltad(nd)
    arg1d(nd, :) = 2*e(:)*ed(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
    IF (arg2 .EQ. 0.0) THEN
      abs_ed(nd) = 0.D0
    ELSE
      abs_ed(nd) = arg2d(nd)/(2.D0*DSQRT(arg2))
    END IF
    arg1d(nd, :) = 2*e_dot(:)*e_dotd(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
    epsd(nd) = 2d0*abs_ed(nd)/s
  END DO
! angle of the position with minimum clearance between shaft and bearing
  delta = DATAN2(e(2), e(1))
! angle of the velocity vector of the shaft relative to the bearing
  phi = DATAN2(e_dot(2), e_dot(1))
! angle between velocity vector and minimum clearance
  kappa = phi - delta
  abs_e = DSQRT(arg2)
! absolute value of the velocity of the shaft relative to the bearing
  arg1(:) = e_dot(:)**2
  arg2 = SUM(arg1(:))
  abs_e_dot = DSQRT(arg2)
  DO nd=1,nbdirs
    IF (arg2 .EQ. 0.0) THEN
      abs_e_dotd(nd) = 0.D0
    ELSE
      abs_e_dotd(nd) = arg2d(nd)/(2.D0*DSQRT(arg2))
    END IF
    eps_dotd(nd) = 2d0*(DCOS(kappa)*abs_e_dotd(nd)-kappad(nd)*DSIN(kappa&
&      )*abs_e_dot)/s
  END DO
! time derivative of the relative eccentricity of the shaft
  eps_dot = 2d0*DCOS(kappa)*abs_e_dot/s
! relative eccentricity of the shaft
  eps = 2d0*abs_e/s
  IF (eps_dot .NE. 0d0) THEN
    DO nd=1,nbdirs
      epsd(nd) = DSIGN(1d0, eps_dot)*epsd(nd)
    END DO
! eps is positive if it's time derivative is positive too
!       attention the signum function is zero if eps_dot is zero
!       but eps must not be zero in this case
    eps = DSIGN(1d0, eps_dot)*eps
  END IF
! time derivative of angle of minimum clearance
  IF (abs_e .EQ. 0d0) THEN
! avoid division by zero
    delta_dot = 0d0
    DO nd=1,nbdirs
      delta_dotd(nd) = 0.D0
    END DO
  ELSE
    DO nd=1,nbdirs
      delta_dotd(nd) = ((ed(nd, 1)*e_dot(2)+e(1)*e_dotd(nd, 2)-ed(nd, 2)&
&        *e_dot(1)-e(2)*e_dotd(nd, 1))*(e(2)**2+e(1)**2)-(e(1)*e_dot(2)-e&
&        (2)*e_dot(1))*(2*e(2)*ed(nd, 2)+2*e(1)*ed(nd, 1)))/(e(2)**2+e(1)&
&        **2)**2
    END DO
!~ delta_dot = DSIN(kappa) * abs_e_dot / abs_e;
    delta_dot = (e(1)*e_dot(2)-e(2)*e_dot(1))/(e(2)**2+e(1)**2)
  END IF
! coefficients for SoD for rotation according to Butenschoen
  a1 = 1.1642d0 - 1.9456d0*(b/d) + 7.1161d0*(b/d)**2 - 10.1073d0*(b/d)**&
&    3 + 5.0141d0*(b/d)**4
  a2 = -1.000026d0 - 0.023634d0*(b/d) - 0.4215d0*(b/d)**2 - 0.038817d0*(&
&    b/d)**3 - 0.090551d0*(b/d)**4
! coefficients for beta for rotation
  a3 = 1.152624d0 - 0.104565d0*(b/d)
  a4 = -2.5905d0 + 0.798745d0*(b/d)
  a5 = 8.73393d0 - 2.3291d0*(b/d)
  a6 = -13.3414d0 + 3.424337d0*(b/d)
  a7 = 6.6294d0 - 1.591732d0*(b/d)
! coefficients for displacement
  a8 = 0.70038d0 + 3.2415d0*(b/d) - 12.2486d0*(b/d)**2 + 18.895d0*(b/d)&
&    **3 - 9.3561d0*(b/d)**4
  a9 = -0.999935d0 + 0.0157434d0*(b/d) - 0.74224d0*(b/d)**2 + 0.42278d0*&
&    (b/d)**3 - 0.368928*(b/d)**4
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs1d(nd) = epsd(nd)
    END DO
    dabs1 = eps
  ELSE
    DO nd=1,nbdirs
      dabs1d(nd) = -epsd(nd)
    END DO
    dabs1 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs7d(nd) = epsd(nd)
    END DO
    dabs7 = eps
  ELSE
    DO nd=1,nbdirs
      dabs7d(nd) = -epsd(nd)
    END DO
    dabs7 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs10d(nd) = epsd(nd)
    END DO
    dabs10 = eps
  ELSE
    DO nd=1,nbdirs
      dabs10d(nd) = -epsd(nd)
    END DO
    dabs10 = -eps
  END IF
! Sommerfeld number for rotation according to Butenschoen 1976
  arg10 = pi**2*(1d0-eps**2) + 16d0*eps**2
  result1 = DSQRT(arg10)
! Sommerfeld number for displacement according to Butenschoen 1976
  pwx1 = 1d0 - eps**2
  DO nd=1,nbdirs
    arg10d(nd) = 16d0*2*eps*epsd(nd) - pi**2*2*eps*epsd(nd)
    IF (arg10 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
      result1d(nd) = arg10d(nd)/(2.D0*DSQRT(arg10))
    END IF
    sodd(nd) = (a1*((b**2*dabs1d(nd)*2d0*(1d0-eps**2)**2/d**2+b**2*dabs1&
&      *2d0*2**2*(1d0-eps**2)*eps*epsd(nd)/d**2)*result1*(dabs7-1d0)/(2d0&
&      **2*(1d0-eps**2)**4)+b**2*dabs1*(result1d(nd)*(dabs7-1d0)+result1*&
&      dabs7d(nd))/(d**2*2d0*(1d0-eps**2)**2))*(a2+dabs10)-b**2*dabs1*&
&      result1*a1*(dabs7-1d0)*dabs10d(nd)/(d**2*2d0*(1d0-eps**2)**2))/(a2&
&      +dabs10)**2
    pwx1d(nd) = -(2*eps*epsd(nd))
    IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. -(5d0/2d0) .EQ. INT(-(&
&        5d0/2d0)))) THEN
      pwr1d(nd) = -(5d0*pwx1**((-1)-5d0/2d0)*pwx1d(nd)/2d0)
    ELSE IF (pwx1 .EQ. 0.0 .AND. -(5d0/2d0) .EQ. 1.0) THEN
      pwr1d(nd) = pwx1d(nd)
    ELSE
      pwr1d(nd) = 0.0
    END IF
    arg10d(nd) = -(2*eps*epsd(nd))
  END DO
  sod = (b/d)**2*dabs1/(2d0*(1d0-eps**2)**2)*result1*a1*(dabs7-1d0)/(a2+&
&    dabs10)
  pwr1 = pwx1**(-(5d0/2d0))
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
  DO nd=1,nbdirs
    IF (arg10 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
      result1d(nd) = arg10d(nd)/(2.D0*DSQRT(arg10))
    END IF
    sovd(nd) = (4d0*b**2*a8*((pwr1d(nd)*(1d0-eps)-pwr1*epsd(nd))*((pi/&
&      2d0-1d0/2d0*DACOS(eps))*(1d0+2d0*eps**2)+3d0/2d0*eps*result1)+pwr1&
&      *(1d0-eps)*(epsd(nd)*(1d0+2d0*eps**2)/(2d0*SQRT(1.D0-eps**2))+(pi/&
&      2d0-1d0/2d0*DACOS(eps))*2d0*2*eps*epsd(nd)+3d0*(epsd(nd)*result1+&
&      eps*result1d(nd))/2d0))*(-a9-eps)/d**2+4d0*b**2*pwr1*((pi/2d0-1d0/&
&      2d0*DACOS(eps))*(1d0+2d0*eps**2)+3d0/2d0*eps*result1)*a8*(1d0-eps)&
&      *epsd(nd)/d**2)/(-a9-eps)**2
  END DO
  sov = 4d0*(b/d)**2*pwr1*((pi/2d0-1d0/2d0*DACOS(eps))*(1d0+2d0*eps**2)+&
&    3d0/2d0*eps*result1)*a8*(1d0-eps)/(-a9-eps)
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs2d(nd) = epsd(nd)
    END DO
    dabs2 = eps
  ELSE
    DO nd=1,nbdirs
      dabs2d(nd) = -epsd(nd)
    END DO
    dabs2 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs8d(nd) = epsd(nd)
    END DO
    dabs8 = eps
  ELSE
    DO nd=1,nbdirs
      dabs8d(nd) = -epsd(nd)
    END DO
    dabs8 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs11d(nd) = epsd(nd)
    END DO
    dabs11 = eps
  ELSE
    DO nd=1,nbdirs
      dabs11d(nd) = -epsd(nd)
    END DO
    dabs11 = -eps
  END IF
! angle between force for rotation and minimum clearance according to Butenschoen 1976
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
  arg2 = pi*result1/(2d0*dabs2)
  DO nd=1,nbdirs
    arg10d(nd) = -(2*eps*epsd(nd))
    IF (arg10 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
      result1d(nd) = arg10d(nd)/(2.D0*DSQRT(arg10))
    END IF
    arg2d(nd) = (pi*result1d(nd)*2d0*dabs2-pi*result1*2d0*dabs2d(nd))/(&
&      2d0*dabs2)**2
    betad(nd) = arg2d(nd)*(a3+a4*dabs8+a5*eps**2+a6*dabs11**3+a7*eps**4)&
&      /(1.0+arg2**2) + ATAN(arg2)*(a4*dabs8d(nd)+a5*2*eps*epsd(nd)+a6*3*&
&      dabs11**2*dabs11d(nd)+a7*4*eps**3*epsd(nd))
  END DO
  beta = ATAN(arg2)*(a3+a4*dabs8+a5*eps**2+a6*dabs11**3+a7*eps**4)
  IF ((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot) .GE. 0.) &
&  THEN
    DO nd=1,nbdirs
      dabs3d(nd) = ((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1)-2d0&
&        *delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(nd, 1)-2d0&
&        *delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot)**2
    END DO
    dabs3 = (omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot)
  ELSE
    DO nd=1,nbdirs
      dabs3d(nd) = -(((omegad(nd, 1)-omegad(nd, 2))*(omega(2)+omega(1)-&
&        2d0*delta_dot)-(omega(1)-omega(2))*(omegad(nd, 2)+omegad(nd, 1)-&
&        2d0*delta_dotd(nd)))/(omega(2)+omega(1)-2d0*delta_dot)**2)
    END DO
    dabs3 = -((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot))
  END IF
  IF (eps .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs9d(nd) = epsd(nd)
    END DO
    dabs9 = eps
  ELSE
    DO nd=1,nbdirs
      dabs9d(nd) = -epsd(nd)
    END DO
    dabs9 = -eps
  END IF
! friction coefficient according to Butenschoen
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
! effective hydrodynamic angular velocity according to Butenschoen 1976
  omega_res = omega(1) + omega(2) - 2d0*delta_dot
  DO nd=1,nbdirs
    arg10d(nd) = -(2*eps*epsd(nd))
    IF (arg10 .EQ. 0.0) THEN
      result1d(nd) = 0.D0
    ELSE
      result1d(nd) = arg10d(nd)/(2.D0*DSQRT(arg10))
    END IF
    myd(nd) = psi*((pi*dabs3d(nd)*result1*sod-dabs3*pi*(result1d(nd)*sod&
&      +result1*sodd(nd)))/(result1*sod)**2+(betad(nd)*DCOS(beta)*dabs9+&
&      DSIN(beta)*dabs9d(nd))/2d0)
    omega_resd(nd) = omegad(nd, 1) + omegad(nd, 2) - 2d0*delta_dotd(nd)
    alphad(nd) = deltad(nd) - DSIGN(1d0, omega_res)*betad(nd)
  END DO
  my = psi*(dabs3*pi/(result1*sod)+DSIN(beta)*dabs9/2d0)
! angle of the force for rotation
  alpha = delta - beta*DSIGN(1d0, omega_res)
  IF (omega_res .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs4d(nd) = omega_resd(nd)
    END DO
    dabs4 = omega_res
  ELSE
    DO nd=1,nbdirs
      dabs4d(nd) = -omega_resd(nd)
    END DO
    dabs4 = -omega_res
  END IF
  DO nd=1,nbdirs
    abs_fdd(nd) = b*d*eta*(sodd(nd)*dabs4+sod*dabs4d(nd))/psi**2
    abs_fvd(nd) = b*d*eta*(sovd(nd)*eps_dot+sov*eps_dotd(nd))/psi**2
  END DO
! absolute value of the force for rotation
  abs_fd = sod*(b*d*eta*dabs4)/psi**2
! absolute value of the force for displacement
  abs_fv = sov*(b*d*eta*eps_dot)/psi**2
  IF (eps .GE. 0.) THEN
    dabs5 = eps
  ELSE
    dabs5 = -eps
  END IF
  IF (dabs5 .LT. 1d-6) THEN
    IF (omega(1) - omega(2) .GE. 0.) THEN
      DO nd=1,nbdirs
        dabs6d(nd) = omegad(nd, 1) - omegad(nd, 2)
      END DO
      dabs6 = omega(1) - omega(2)
    ELSE
      DO nd=1,nbdirs
        dabs6d(nd) = -(omegad(nd, 1)-omegad(nd, 2))
      END DO
      dabs6 = -(omega(1)-omega(2))
    END IF
    DO nd=1,nbdirs
      abs_mrd(nd) = pi*b*d**2*eta*dabs6d(nd)/psi/2d0
    END DO
! avoid division infinite by infinite in case of zero relative eccentricity
! use analytical limit of abs_MR for eps going to zero
    abs_mr = pi*b*d**2*eta*dabs6/psi/2d0
  ELSE
    DO nd=1,nbdirs
      abs_mrd(nd) = d*(myd(nd)*abs_fd+my*abs_fdd(nd))/2d0
    END DO
    abs_mr = my*abs_fd*d/2d0
  END IF
! friction torque
  arg10 = omega(1) - omega(2)
  DO nd=1,nbdirs
    kd(nd, :) = 0.D0
    kd(nd, 1) = abs_fdd(nd)*DCOS(alpha) - abs_fd*alphad(nd)*DSIN(alpha) &
&      + abs_fvd(nd)*DCOS(delta) - abs_fv*deltad(nd)*DSIN(delta)
    kd(nd, 2) = abs_fdd(nd)*DSIN(alpha) + abs_fd*alphad(nd)*DCOS(alpha) &
&      + abs_fvd(nd)*DSIN(delta) + abs_fv*deltad(nd)*DCOS(delta)
    kd(nd, 3) = DSIGN(1d0, arg10)*abs_mrd(nd)
  END DO
! sum of force for rotation and force for displacement
  k(1) = abs_fd*DCOS(alpha) + abs_fv*DCOS(delta)
  k(2) = abs_fd*DSIN(alpha) + abs_fv*DSIN(delta)
  k(3) = abs_mr*DSIGN(1d0, arg10)
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_DV

SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_CD(b, d, psi, eta, omega, e&
&  , e_dot, k, eps, eps_dot, delta, sod, sov, my, beta)
  USE MATH_CONST_DV
  IMPLICIT NONE
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!       PRINT OUT OF RESULTS
!--------------------------------------------------------------------------------------------------------------------------------
!---------------
!~ PRINT *,'HYDRODYNAMIC_PLAIN_BEARING_FORCE'
!~ PRINT *,'b=',b
!~ PRINT *,'d=',d
!~ PRINT *,'eta=',eta
!~ PRINT *,'Psi=',Psi
!~ PRINT *,'e=',e
!~ PRINT *,'e_dot=',e_dot
!~ PRINT *,'omega_w=',omega_w
!~ PRINT *,'delta=',delta
!~ PRINT *,'phi=',phi
!~ PRINT *,'kappa=',kappa
!~ PRINT *,'abs_e=',abs_e
!~ PRINT *,'abs_e_dot=',abs_e_dot
!~ PRINT *,'eps_dot=',eps_dot
!~ PRINT *,'eps=',eps
!~ PRINT *,'delta_dot=',delta_dot
!~ PRINT *,'a1=',a1
!~ PRINT *,'a2=',a2
!~ PRINT *,'a3=',a3
!~ PRINT *,'a4=',a4
!~ PRINT *,'a5=',a5
!~ PRINT *,'a6=',a6
!~ PRINT *,'a7=',a7
!~ PRINT *,'a8=',a8
!~ PRINT *,'a9=',a9
!~ PRINT *,'SoD=',SoD
!~ PRINT *,'SoV=',SoV
!~ PRINT *,'beta=',beta
!~ PRINT *,'my=',my
!~ PRINT *,'alpha=',alpha
!~ PRINT *,'omega_res=',omega_res
!~ PRINT *,'abs_FD=',abs_FD
!~ PRINT *,'abs_FV=',abs_FV
!~ PRINT *,'abs_MR=',abs_MR
!~ PRINT *, 'k=',k
!-------------------------------------------------------------------------------------------------------------------------
!       hydrodynamic plain bearing calculation according to Butenschoen's theory
!-------------------------------------------------------------------------------------------------------------------------
!       COORDINATE SYSTEM:
!       x ... axial direction
!       y, z ... radial direction
!-------------------------------------------------------------------------------------------------------------------------
!        INPUT PARAMETERS
!-------------------------------------------------------------------------------------------------------------------------
!       b   ... bearing width [m]
!       d   ... bearing diameter [m]
!       Psi ... relative radial clearance Psi = ( D - d ) / d [1]
!       eta ... dynamic oil viscosity [Pa*s]
!       omega(1) ... angular velocity of the shaft [rad/s]
!	    omega(2) ... angular velocity of the bearing [rad/s]
!       e ...  radial eccentricity of the shaft 
!               e(1) = ey  [m]
!               e(2) = ez  [m]
!       e_dot(1) ... velocity of the shaft relative to the bearing 
!               e_dot(1) = ey_dot [m/s]
!               e_dot(2) = ez_dot [m/s]
!
!--------------------------------------------------------------------------------------------------------------------------
!       OUTPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------
!       k ... force on the bearing
!               k(1) = Fy [N]
!               k(2) = Fz [N]
!               k(3) = Mx [Nm]
!       eps ...             relative eccentricity of the shaft [1]
!       eps_dot ...      time derivative of the relative eccentricity [1/s]
!       delta ...          angle of minimum clearance between shaft and bearing [rad]
!       SoD ....           Sommerfeld number for rotation [1]
!       SoV ...            Sommerfeld number for displacement [1]
!       my ...             friction coefficient [N/N]
!       beta  ...          angle between force for rotation and minimum clearance [rad]
  DOUBLE PRECISION, INTENT(IN) :: b, d, psi, eta, omega(2), e(2), e_dot(&
&  2)
  DOUBLE PRECISION, INTENT(OUT) :: k(3), eps, eps_dot, delta, sod, sov, &
&  my, beta
  DOUBLE PRECISION :: abs_e, abs_e_dot, s
  DOUBLE PRECISION :: delta_dot, alpha, kappa, phi
  DOUBLE PRECISION :: omega_res, abs_fd, abs_fv, abs_mr
  DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9
  INTRINSIC DATAN2, DSQRT, SUM, DCOS, DSIN, DABS, DSIGN
  DOUBLE PRECISION, DIMENSION(2) :: arg1
  DOUBLE PRECISION :: arg2
  DOUBLE PRECISION :: arg10
  DOUBLE PRECISION :: result1
  DOUBLE PRECISION :: pwx1
  DOUBLE PRECISION :: pwr1
  DOUBLE PRECISION :: dabs11
  DOUBLE PRECISION :: dabs10
  DOUBLE PRECISION :: dabs9
  DOUBLE PRECISION :: dabs8
  DOUBLE PRECISION :: dabs7
  DOUBLE PRECISION :: dabs6
  DOUBLE PRECISION :: dabs5
  INTRINSIC DACOS
  DOUBLE PRECISION :: dabs4
  DOUBLE PRECISION :: dabs3
  INTRINSIC ATAN
  DOUBLE PRECISION :: dabs2
  DOUBLE PRECISION :: dabs1
! clearance of the shaft inside the bearing
  s = d*psi/(1-psi)
! angle of the position with minimum clearance between shaft and bearing
  delta = DATAN2(e(2), e(1))
! angle of the velocity vector of the shaft relative to the bearing
  phi = DATAN2(e_dot(2), e_dot(1))
! angle between velocity vector and minimum clearance
  kappa = phi - delta
! absolute value of the eccentricity of the shaft inside the bearing
  arg1(:) = e(:)**2
  arg2 = SUM(arg1(:))
  abs_e = DSQRT(arg2)
! absolute value of the velocity of the shaft relative to the bearing
  arg1(:) = e_dot(:)**2
  arg2 = SUM(arg1(:))
  abs_e_dot = DSQRT(arg2)
! time derivative of the relative eccentricity of the shaft
  eps_dot = 2d0*DCOS(kappa)*abs_e_dot/s
! relative eccentricity of the shaft
  eps = 2d0*abs_e/s
  IF (eps_dot .NE. 0d0) eps = DSIGN(1d0, eps_dot)*eps
! eps is positive if it's time derivative is positive too
!       attention the signum function is zero if eps_dot is zero
!       but eps must not be zero in this case
! time derivative of angle of minimum clearance
  IF (abs_e .EQ. 0d0) THEN
! avoid division by zero
    delta_dot = 0d0
  ELSE
!~ delta_dot = DSIN(kappa) * abs_e_dot / abs_e;
    delta_dot = (e(1)*e_dot(2)-e(2)*e_dot(1))/(e(2)**2+e(1)**2)
  END IF
! coefficients for SoD for rotation according to Butenschoen
  a1 = 1.1642d0 - 1.9456d0*(b/d) + 7.1161d0*(b/d)**2 - 10.1073d0*(b/d)**&
&    3 + 5.0141d0*(b/d)**4
  a2 = -1.000026d0 - 0.023634d0*(b/d) - 0.4215d0*(b/d)**2 - 0.038817d0*(&
&    b/d)**3 - 0.090551d0*(b/d)**4
! coefficients for beta for rotation
  a3 = 1.152624d0 - 0.104565d0*(b/d)
  a4 = -2.5905d0 + 0.798745d0*(b/d)
  a5 = 8.73393d0 - 2.3291d0*(b/d)
  a6 = -13.3414d0 + 3.424337d0*(b/d)
  a7 = 6.6294d0 - 1.591732d0*(b/d)
! coefficients for displacement
  a8 = 0.70038d0 + 3.2415d0*(b/d) - 12.2486d0*(b/d)**2 + 18.895d0*(b/d)&
&    **3 - 9.3561d0*(b/d)**4
  a9 = -0.999935d0 + 0.0157434d0*(b/d) - 0.74224d0*(b/d)**2 + 0.42278d0*&
&    (b/d)**3 - 0.368928*(b/d)**4
  IF (eps .GE. 0.) THEN
    dabs1 = eps
  ELSE
    dabs1 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    dabs7 = eps
  ELSE
    dabs7 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    dabs10 = eps
  ELSE
    dabs10 = -eps
  END IF
! Sommerfeld number for rotation according to Butenschoen 1976
  arg10 = pi**2*(1d0-eps**2) + 16d0*eps**2
  result1 = DSQRT(arg10)
  sod = (b/d)**2*dabs1/(2d0*(1d0-eps**2)**2)*result1*a1*(dabs7-1d0)/(a2+&
&    dabs10)
! Sommerfeld number for displacement according to Butenschoen 1976
  pwx1 = 1d0 - eps**2
  pwr1 = pwx1**(-(5d0/2d0))
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
  sov = 4d0*(b/d)**2*pwr1*((pi/2d0-1d0/2d0*DACOS(eps))*(1d0+2d0*eps**2)+&
&    3d0/2d0*eps*result1)*a8*(1d0-eps)/(-a9-eps)
  IF (eps .GE. 0.) THEN
    dabs2 = eps
  ELSE
    dabs2 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    dabs8 = eps
  ELSE
    dabs8 = -eps
  END IF
  IF (eps .GE. 0.) THEN
    dabs11 = eps
  ELSE
    dabs11 = -eps
  END IF
! angle between force for rotation and minimum clearance according to Butenschoen 1976
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
  arg2 = pi*result1/(2d0*dabs2)
  beta = ATAN(arg2)*(a3+a4*dabs8+a5*eps**2+a6*dabs11**3+a7*eps**4)
  IF ((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot) .GE. 0.) &
&  THEN
    dabs3 = (omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot)
  ELSE
    dabs3 = -((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot))
  END IF
  IF (eps .GE. 0.) THEN
    dabs9 = eps
  ELSE
    dabs9 = -eps
  END IF
! friction coefficient according to Butenschoen
  arg10 = 1d0 - eps**2
  result1 = DSQRT(arg10)
  my = psi*(dabs3*pi/(result1*sod)+DSIN(beta)*dabs9/2d0)
! effective hydrodynamic angular velocity according to Butenschoen 1976
  omega_res = omega(1) + omega(2) - 2d0*delta_dot
! angle of the force for rotation
  alpha = delta - beta*DSIGN(1d0, omega_res)
  IF (omega_res .GE. 0.) THEN
    dabs4 = omega_res
  ELSE
    dabs4 = -omega_res
  END IF
! absolute value of the force for rotation
  abs_fd = sod*(b*d*eta*dabs4)/psi**2
! absolute value of the force for displacement
  abs_fv = sov*(b*d*eta*eps_dot)/psi**2
  IF (eps .GE. 0.) THEN
    dabs5 = eps
  ELSE
    dabs5 = -eps
  END IF
  IF (dabs5 .LT. 1d-6) THEN
    IF (omega(1) - omega(2) .GE. 0.) THEN
      dabs6 = omega(1) - omega(2)
    ELSE
      dabs6 = -(omega(1)-omega(2))
    END IF
! avoid division infinite by infinite in case of zero relative eccentricity
! use analytical limit of abs_MR for eps going to zero
    abs_mr = pi*b*d**2*eta*dabs6/psi/2d0
  ELSE
    abs_mr = my*abs_fd*d/2d0
  END IF
! sum of force for rotation and force for displacement
  k(1) = abs_fd*DCOS(alpha) + abs_fv*DCOS(delta)
  k(2) = abs_fd*DSIN(alpha) + abs_fv*DSIN(delta)
! friction torque
  arg10 = omega(1) - omega(2)
  k(3) = abs_mr*DSIGN(1d0, arg10)
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE_CD

!  Differentiation of plain_bearing_contact_force in forward (tangent) mode: (multi-directional mode)
!   variations   of useful results: k
!   with respect to varying inputs: e e_dot phi_dot
SUBROUTINE PLAIN_BEARING_CONTACT_FORCE_DV(d, psi, sp, dl, m, abs_fprs1, &
&  myp, signum_delta_omega, phi_dot, phi_dotd, e, ed, e_dot, e_dotd, k, &
&  kd, nbdirs)
  USE MATH_CONST_DV
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
!    PRINT *,'e=',e
!    PRINT *,'e_dot=',e_dot
!    PRINT *,'eps=',eps
!    PRINT *,'eps_dot=',eps_dot
!    PRINT *,'abs_FPrs=',abs_FPrs
!    PRINT *,'abs_FPrd=',abs_FPrd
!    PRINT *,'abs_MP=',abs_MP
!    PRINT *,'FP23r=',FP23r
!    PRINT *,'FP23t=',FP23t
!    PRINT *,'k=',k
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!       phenomenological contact force calculation between shaft and bearing with penalty method 
!       phenomenological damping force with Lehr's damping law
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!       INPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!   d	                 shaft diameter       [m]
!   Psi                  relative clearance s = d * Psi / ( 1 - Psi ), D = d + s [m/m]
!   sP                   contact stiffness for wall contact ( epsilon == 1 ) 			[N/m]
!   DL		             Lehr damping ratio for wall contact					[1]
!   m		             rated mass for which DL is valid				[kg]
!   abs_FPrs1	         contact force for wall contact ( epsilon == 1 )	[N]
!   myP		             Coulomb friction coefficient ( abs_FPt == myP * abs_FPr )		[N/N]
!   signum_delta_omega   threshold for the angular velocity			[rad/s]
!   e                    eccentricity of the shaft                                                                       
!       e(1) = ey [m]
!       e(2) = ez [m]
!   e_dot                             velocity of the shaft
!       e_dot(1) = ey_dot [m/s]
!       e_dot(2) = ez_dot [m/s]
!   Phi_dot(1)                     angular velocity of the shaft [rad/s]
!   Phi_dot(2)                     angular velocity of the bearing [rad/s]
!
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------
!       OUTPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------
!   k(1)  = Fy [N]        force at the bearing in y direction
!   k(2)  = Fz [N]        force at the bearing in z direction
!   k(3)  = Mx [Nm]       torque at the bearing around the x axis
  DOUBLE PRECISION, INTENT(IN) :: d, psi, sp, dl, m, abs_fprs1, myp, &
&  signum_delta_omega, phi_dot(2), e(2), e_dot(2)
  DOUBLE PRECISION, INTENT(IN) :: phi_dotd(nbdirsmax, 2), ed(nbdirsmax, &
&  2), e_dotd(nbdirsmax, 2)
  DOUBLE PRECISION, INTENT(OUT) :: k(3)
  DOUBLE PRECISION, INTENT(OUT) :: kd(nbdirsmax, 3)
  DOUBLE PRECISION :: eps, eps_dot, delta, phi, kappa, abs_e, abs_e_dot&
&  , s
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: epsd, eps_dotd, deltad, phid&
&  , kappad, abs_ed, abs_e_dotd
  DOUBLE PRECISION :: abs_fprs, abs_fprd, abs_fpr, fp23r(2), &
&  delta_phi_dot, alpha, sign_mp
  DOUBLE PRECISION :: abs_fprsd(nbdirsmax), abs_fprdd(nbdirsmax), &
&  abs_fprd0(nbdirsmax), fp23rd(nbdirsmax, 2), delta_phi_dotd(nbdirsmax)&
&  , sign_mpd(nbdirsmax)
  DOUBLE PRECISION :: fp23t(2), abs_mp
  DOUBLE PRECISION :: fp23td(nbdirsmax, 2), abs_mpd(nbdirsmax)
  INTRINSIC SQRT, EXP, COS, SIN, ABS
  EXTERNAL SIGN_MOD
  EXTERNAL SIGN_MOD_DV
  DOUBLE PRECISION :: SIGN_MOD_CD
  DOUBLE PRECISION, DIMENSION(2) :: arg1
  DOUBLE PRECISION, DIMENSION(nbdirsmax, 2) :: arg1d
  DOUBLE PRECISION :: arg2
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg2d
  DOUBLE PRECISION :: arg10
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg10d
  DOUBLE PRECISION :: result1
  INTEGER :: nd
  INTEGER :: nbdirs
  INTRINSIC DCOS
  INTRINSIC DEXP
  DOUBLE PRECISION :: dabs1d(nbdirsmax)
  INTRINSIC DSIN
  INTRINSIC DATAN2
  INTRINSIC DABS
  INTRINSIC SUM
  DOUBLE PRECISION :: dabs2d(nbdirsmax)
  DOUBLE PRECISION :: dabs2
  DOUBLE PRECISION :: dabs1
  INTRINSIC DSQRT
! clearance between shaft and bearing
  s = d*psi/(1-psi)
! absolute value of the eccentricity of the shaft
  arg1(:) = e(:)**2
  arg2 = SUM(arg1(:))
  abs_e = DSQRT(arg2)
! relative eccentricity
  eps = 2d0*abs_e/s
! stiffness contact force
  arg10 = d*psi/2d0*sp/abs_fprs1*(eps-1d0)
  DO nd=1,nbdirs
    deltad(nd) = (ed(nd, 2)*e(1)-ed(nd, 1)*e(2))/(e(2)**2+e(1)**2)
    phid(nd) = (e_dotd(nd, 2)*e_dot(1)-e_dotd(nd, 1)*e_dot(2))/(e_dot(2)&
&      **2+e_dot(1)**2)
    kappad(nd) = phid(nd) - deltad(nd)
    arg1d(nd, :) = 2*e(:)*ed(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
    IF (arg2 .EQ. 0.0) THEN
      abs_ed(nd) = 0.D0
    ELSE
      abs_ed(nd) = arg2d(nd)/(2.D0*DSQRT(arg2))
    END IF
    arg1d(nd, :) = 2*e_dot(:)*e_dotd(nd, :)
    arg2d(nd) = SUM(arg1d(nd, :))
    epsd(nd) = 2d0*abs_ed(nd)/s
    arg10d(nd) = d*psi*sp*epsd(nd)/(2d0*abs_fprs1)
    abs_fprsd(nd) = abs_fprs1*arg10d(nd)*DEXP(arg10)
    fp23rd(nd, :) = 0.D0
    delta_phi_dotd(nd) = phi_dotd(nd, 1) - psi*epsd(nd)*phi_dot(2)/(1-&
&      psi) - (1d0+psi*eps/(1-psi))*phi_dotd(nd, 2)
  END DO
! angle at minimum clearance
  delta = DATAN2(e(2), e(1))
! moving direction of the shaft with respect to the bearing
  phi = DATAN2(e_dot(2), e_dot(1))
! angle between moving direction and minimum clearance
  kappa = phi - delta
! absolute value of the velocity of the shaft
  arg1(:) = e_dot(:)**2
  arg2 = SUM(arg1(:))
  abs_e_dot = DSQRT(arg2)
  DO nd=1,nbdirs
    IF (arg2 .EQ. 0.0) THEN
      abs_e_dotd(nd) = 0.D0
    ELSE
      abs_e_dotd(nd) = arg2d(nd)/(2.D0*DSQRT(arg2))
    END IF
    eps_dotd(nd) = 2d0*(DCOS(kappa)*abs_e_dotd(nd)-kappad(nd)*DSIN(kappa&
&      )*abs_e_dot)/s
    arg2d(nd) = d*psi*sp*epsd(nd)/(4d0*abs_fprs1)
  END DO
! derivative of the relative eccentricity versus time
  eps_dot = 2d0*DCOS(kappa)*abs_e_dot/s
  abs_fprs = abs_fprs1*DEXP(arg10)
! damping contact force
  arg10 = m*sp
  result1 = DSQRT(arg10)
  arg2 = d*psi/4d0*sp/abs_fprs1*(eps-1d0)
  abs_fprd = d*psi*result1*dl*eps_dot*DEXP(arg2)
! contact force
  abs_fpr = abs_fprs + abs_fprd
  DO nd=1,nbdirs
    abs_fprdd(nd) = d*psi*result1*dl*(eps_dotd(nd)*DEXP(arg2)+eps_dot*&
&      arg2d(nd)*DEXP(arg2))
    abs_fprd0(nd) = abs_fprsd(nd) + abs_fprdd(nd)
    fp23rd(nd, 1) = abs_fprd0(nd)*DCOS(delta) - abs_fpr*deltad(nd)*DSIN(&
&      delta)
    fp23rd(nd, 2) = abs_fprd0(nd)*DSIN(delta) + abs_fpr*deltad(nd)*DCOS(&
&      delta)
  END DO
! radial contact force
  fp23r(1) = abs_fpr*DCOS(delta)
  fp23r(2) = abs_fpr*DSIN(delta)
! The difference between the angular velocity of the shaft and the bearing is related to the circumference of the shaft.
! The term 1 + Psi * eps is needed because the sign of the torque
! must be related to the difference of the circumferential velocity.
  delta_phi_dot = phi_dot(1) - (1d0+psi*eps/(1-psi))*phi_dot(2)
! angle between radial and tangential component of the the contact force
! 	Attention: the signum function must not be used here since this would lead to a wrong direction of the force
! 	The direction of the tangential component of the contact force must be either pi/2 oder -pi/2
  IF (delta_phi_dot .GT. 0d0) THEN
    alpha = pi/2d0
  ELSE
    alpha = -(pi/2d0)
  END IF
! sign of the friction torque
! 	For numerical reasons the modified signum function must be used here
  CALL SIGN_MOD_DV(delta_phi_dot, delta_phi_dotd, signum_delta_omega, &
&             sign_mp, sign_mpd, nbdirs)
  IF (sign_mp .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs1d(nd) = sign_mpd(nd)
    END DO
    dabs1 = sign_mp
  ELSE
    DO nd=1,nbdirs
      dabs1d(nd) = -sign_mpd(nd)
    END DO
    dabs1 = -sign_mp
  END IF
! vector of the tangential component of the contact force
!	The term abs(sign_MP) is needed because the tangential force and the torque must be consistent
  arg10 = delta + alpha
  DO nd=1,nbdirs
    arg10d(nd) = deltad(nd)
    fp23td(nd, :) = 0.D0
    fp23td(nd, 1) = myp*((dabs1d(nd)*abs_fpr+dabs1*abs_fprd0(nd))*DCOS(&
&      arg10)-dabs1*abs_fpr*arg10d(nd)*DSIN(arg10))
  END DO
  fp23t(1) = dabs1*myp*abs_fpr*DCOS(arg10)
  IF (sign_mp .GE. 0.) THEN
    DO nd=1,nbdirs
      dabs2d(nd) = sign_mpd(nd)
    END DO
    dabs2 = sign_mp
  ELSE
    DO nd=1,nbdirs
      dabs2d(nd) = -sign_mpd(nd)
    END DO
    dabs2 = -sign_mp
  END IF
  arg10 = delta + alpha
! absolute value of the friction torque
  abs_mp = myp*abs_fpr*d/2d0
  DO nd=1,nbdirs
    arg10d(nd) = deltad(nd)
    fp23td(nd, 2) = myp*((dabs2d(nd)*abs_fpr+dabs2*abs_fprd0(nd))*DSIN(&
&      arg10)+dabs2*abs_fpr*arg10d(nd)*DCOS(arg10))
    abs_mpd(nd) = myp*d*abs_fprd0(nd)/2d0
    kd(nd, :) = 0.D0
    kd(nd, 1:2) = fp23rd(nd, :) + fp23td(nd, :)
    kd(nd, 3) = (sign_mpd(nd)*abs_mp+sign_mp*abs_mpd(nd))*(1d0+psi*eps/(&
&      1-psi)) + sign_mp*abs_mp*psi*epsd(nd)/(1-psi)
  END DO
  fp23t(2) = dabs2*myp*abs_fpr*DSIN(arg10)
! the radial and tangential component of the force are applied to the center of the bearing
  k(1:2) = fp23r + fp23t
! torque of the tangential component of the contact force around the x axis
  k(3) = sign_mp*abs_mp*(1d0+psi*eps/(1-psi))
END SUBROUTINE PLAIN_BEARING_CONTACT_FORCE_DV

!  Differentiation of sign_mod in forward (tangent) mode: (multi-directional mode)
!   variations   of useful results: sign_mod
!   with respect to varying inputs: x
SUBROUTINE SIGN_MOD_DV(x, xd, delta, sign_mod, sign_modd, nbdirs)
  USE MATH_CONST_DV
  USE DIFFSIZES
!  Hint: nbdirsmax should be the maximum number of differentiation directions
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x, delta
  DOUBLE PRECISION, DIMENSION(nbdirsmax), INTENT(IN) :: xd
  INTRINSIC DSIGN, DTANH
  DOUBLE PRECISION :: arg1
  DOUBLE PRECISION, DIMENSION(nbdirsmax) :: arg1d
  INTEGER :: nd
  DOUBLE PRECISION :: sign_mod
  DOUBLE PRECISION :: sign_modd(nbdirsmax)
  INTEGER :: nbdirs
  IF (delta .EQ. 0d0) THEN
    sign_mod = DSIGN(1d0, x)
    DO nd=1,nbdirs
      sign_modd(nd) = 0.D0
    END DO
  ELSE
    arg1 = 2*pi*x/delta
    DO nd=1,nbdirs
      arg1d(nd) = 2*pi*xd(nd)/delta
      sign_modd(nd) = arg1d(nd)*(1.D0-TANH(arg1)**2)
    END DO
    sign_mod = DTANH(arg1)
  END IF
END SUBROUTINE SIGN_MOD_DV

SUBROUTINE PLAIN_BEARING_CONTACT_FORCE_CD(d, psi, sp, dl, m, abs_fprs1, &
&  myp, signum_delta_omega, phi_dot, e, e_dot, k)
  USE MATH_CONST_DV
  IMPLICIT NONE
!    PRINT *,'e=',e
!    PRINT *,'e_dot=',e_dot
!    PRINT *,'eps=',eps
!    PRINT *,'eps_dot=',eps_dot
!    PRINT *,'abs_FPrs=',abs_FPrs
!    PRINT *,'abs_FPrd=',abs_FPrd
!    PRINT *,'abs_MP=',abs_MP
!    PRINT *,'FP23r=',FP23r
!    PRINT *,'FP23t=',FP23t
!    PRINT *,'k=',k
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!       phenomenological contact force calculation between shaft and bearing with penalty method 
!       phenomenological damping force with Lehr's damping law
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!       INPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------
!   d	                 shaft diameter       [m]
!   Psi                  relative clearance s = d * Psi / ( 1 - Psi ), D = d + s [m/m]
!   sP                   contact stiffness for wall contact ( epsilon == 1 ) 			[N/m]
!   DL		             Lehr damping ratio for wall contact					[1]
!   m		             rated mass for which DL is valid				[kg]
!   abs_FPrs1	         contact force for wall contact ( epsilon == 1 )	[N]
!   myP		             Coulomb friction coefficient ( abs_FPt == myP * abs_FPr )		[N/N]
!   signum_delta_omega   threshold for the angular velocity			[rad/s]
!   e                    eccentricity of the shaft                                                                       
!       e(1) = ey [m]
!       e(2) = ez [m]
!   e_dot                             velocity of the shaft
!       e_dot(1) = ey_dot [m/s]
!       e_dot(2) = ez_dot [m/s]
!   Phi_dot(1)                     angular velocity of the shaft [rad/s]
!   Phi_dot(2)                     angular velocity of the bearing [rad/s]
!
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------
!       OUTPUT PARAMETERS
!--------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------
!   k(1)  = Fy [N]        force at the bearing in y direction
!   k(2)  = Fz [N]        force at the bearing in z direction
!   k(3)  = Mx [Nm]       torque at the bearing around the x axis
  DOUBLE PRECISION, INTENT(IN) :: d, psi, sp, dl, m, abs_fprs1, myp, &
&  signum_delta_omega, phi_dot(2), e(2), e_dot(2)
  DOUBLE PRECISION, INTENT(OUT) :: k(3)
  DOUBLE PRECISION :: eps, eps_dot, delta, phi, kappa, abs_e, abs_e_dot&
&  , s
  DOUBLE PRECISION :: abs_fprs, abs_fprd, abs_fpr, fp23r(2), &
&  delta_phi_dot, alpha, sign_mp
  DOUBLE PRECISION :: fp23t(2), abs_mp
  INTRINSIC SQRT, EXP, COS, SIN, ABS
  EXTERNAL SIGN_MOD
  DOUBLE PRECISION :: SIGN_MOD_CD
  DOUBLE PRECISION, DIMENSION(2) :: arg1
  DOUBLE PRECISION :: arg2
  DOUBLE PRECISION :: arg10
  DOUBLE PRECISION :: result1
  INTRINSIC DCOS
  INTRINSIC DEXP
  INTRINSIC DSIN
  INTRINSIC DATAN2
  INTRINSIC DABS
  INTRINSIC SUM
  DOUBLE PRECISION :: dabs2
  DOUBLE PRECISION :: dabs1
  INTRINSIC DSQRT
! clearance between shaft and bearing
  s = d*psi/(1-psi)
! angle at minimum clearance
  delta = DATAN2(e(2), e(1))
! moving direction of the shaft with respect to the bearing
  phi = DATAN2(e_dot(2), e_dot(1))
! angle between moving direction and minimum clearance
  kappa = phi - delta
! absolute value of the eccentricity of the shaft
  arg1(:) = e(:)**2
  arg2 = SUM(arg1(:))
  abs_e = DSQRT(arg2)
! absolute value of the velocity of the shaft
  arg1(:) = e_dot(:)**2
  arg2 = SUM(arg1(:))
  abs_e_dot = DSQRT(arg2)
! derivative of the relative eccentricity versus time
  eps_dot = 2d0*DCOS(kappa)*abs_e_dot/s
! relative eccentricity
  eps = 2d0*abs_e/s
! stiffness contact force
  arg10 = d*psi/2d0*sp/abs_fprs1*(eps-1d0)
  abs_fprs = abs_fprs1*DEXP(arg10)
! damping contact force
  arg10 = m*sp
  result1 = DSQRT(arg10)
  arg2 = d*psi/4d0*sp/abs_fprs1*(eps-1d0)
  abs_fprd = d*psi*result1*dl*eps_dot*DEXP(arg2)
! contact force
  abs_fpr = abs_fprs + abs_fprd
! radial contact force
  fp23r(1) = abs_fpr*DCOS(delta)
  fp23r(2) = abs_fpr*DSIN(delta)
! The difference between the angular velocity of the shaft and the bearing is related to the circumference of the shaft.
! The term 1 + Psi * eps is needed because the sign of the torque
! must be related to the difference of the circumferential velocity.
  delta_phi_dot = phi_dot(1) - (1d0+psi*eps/(1-psi))*phi_dot(2)
! angle between radial and tangential component of the the contact force
! 	Attention: the signum function must not be used here since this would lead to a wrong direction of the force
! 	The direction of the tangential component of the contact force must be either pi/2 oder -pi/2
  IF (delta_phi_dot .GT. 0d0) THEN
    alpha = pi/2d0
  ELSE
    alpha = -(pi/2d0)
  END IF
! sign of the friction torque
! 	For numerical reasons the modified signum function must be used here
  sign_mp = SIGN_MOD_CD(delta_phi_dot, signum_delta_omega)
  IF (sign_mp .GE. 0.) THEN
    dabs1 = sign_mp
  ELSE
    dabs1 = -sign_mp
  END IF
! vector of the tangential component of the contact force
!	The term abs(sign_MP) is needed because the tangential force and the torque must be consistent
  arg10 = delta + alpha
  fp23t(1) = dabs1*myp*abs_fpr*DCOS(arg10)
  IF (sign_mp .GE. 0.) THEN
    dabs2 = sign_mp
  ELSE
    dabs2 = -sign_mp
  END IF
  arg10 = delta + alpha
  fp23t(2) = dabs2*myp*abs_fpr*DSIN(arg10)
! absolute value of the friction torque
  abs_mp = myp*abs_fpr*d/2d0
! the radial and tangential component of the force are applied to the center of the bearing
  k(1:2) = fp23r + fp23t
! torque of the tangential component of the contact force around the x axis
  k(3) = sign_mp*abs_mp*(1d0+psi*eps/(1-psi))
END SUBROUTINE PLAIN_BEARING_CONTACT_FORCE_CD

FUNCTION SIGN_MOD_CD(x, delta) RESULT (SIGN_MOD)
  USE MATH_CONST_DV
  IMPLICIT NONE
  DOUBLE PRECISION, INTENT(IN) :: x, delta
  INTRINSIC DSIGN, DTANH
  DOUBLE PRECISION :: arg1
  DOUBLE PRECISION :: sign_mod
  IF (delta .EQ. 0d0) THEN
    sign_mod = DSIGN(1d0, x)
  ELSE
    arg1 = 2*pi*x/delta
    sign_mod = DTANH(arg1)
  END IF
END FUNCTION SIGN_MOD_CD

