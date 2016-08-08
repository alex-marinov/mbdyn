! $Header$
! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2015
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

! AUTHOR: Reinhard Resch <R.RESCH@secop.com>
!        Copyright (C) 2011(-2016) all rights reserved.
!
!        The copyright of this code is transferred
!        to Pierangelo Masarati and Paolo Mantegazza
!        for use in the software MBDyn as described
!        in the GNU Public License version 2.1

MODULE HYDRODYNAMIC_PLAIN_BEARING
  USE ISO_C_BINDING
        IMPLICIT NONE
  DOUBLE PRECISION,PARAMETER :: pi = 4D0 * ATAN(1D0)
  TYPE,BIND(C) :: BEARING_DATA
     DOUBLE PRECISION :: b, d, Psi, eta, eps_max, s, a(9)
  END TYPE BEARING_DATA
CONTAINS
  SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_INIT(bdat) BIND(C)
    IMPLICIT NONE
    TYPE(BEARING_DATA),INTENT(INOUT) :: bdat
    ! absolute clearance of the shaft inside the bearing
    bdat%s = bdat%d * bdat%Psi / ( 1 - bdat%Psi )
    
    ! coefficients for SoD for rotation according to Butenschoen
    bdat%a(1) = 1.1642D0 - 1.9456D0 * ( bdat%b / bdat%d ) +7.1161D0 * ( bdat%b / bdat%d )**2 - 10.1073D0 * ( bdat%b / bdat%d )**3 &
         &   +5.0141D0 * ( bdat%b / bdat%d )**4

    bdat%a(2) = -1.000026D0 - 0.023634D0 * ( bdat%b / bdat%d ) - 0.4215D0 * ( bdat%b / bdat%d )**2 &
         &   - 0.038817D0 * ( bdat%b / bdat%d )**3 - 0.090551D0 * ( bdat%b / bdat%d )**4

    ! coefficients for beta for rotation
    bdat%a(3) = 1.152624D0 - 0.104565D0 * ( bdat%b / bdat%d )

    bdat%a(4) = -2.5905D0 + 0.798745D0 * ( bdat%b / bdat%d )

    bdat%a(5) = 8.73393D0 - 2.3291D0 * ( bdat%b / bdat%d )

    bdat%a(6) = -13.3414D0 + 3.424337D0 * ( bdat%b / bdat%d )

    bdat%a(7) = 6.6294D0 - 1.591732D0 * ( bdat%b / bdat%d )

    ! coefficients for displacement
    bdat%a(8) = 0.70038D0+3.2415D0 * ( bdat%b / bdat%d ) - 12.2486D0 * ( bdat%b / bdat%d )**2 + 18.895D0 * ( bdat%b / bdat%d )**3 &
         &   - 9.3561D0 * ( bdat%b / bdat%d )**4

    bdat%a(9) = -0.999935D0 + 0.0157434D0 * ( bdat%b / bdat%d ) - 0.74224D0 * ( bdat%b / bdat%d )**2 &
         &  + 0.42278D0 * ( bdat%b / bdat%d )**3 - 0.368928 * ( bdat%b / bdat%d )**4

  END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_INIT

  SUBROUTINE SOMMERFELD_NUMBERS(bdat, eps, omega, delta_dot, SoD, SoV, beta, mu)
    TYPE(BEARING_DATA),INTENT(IN) :: bdat
    DOUBLE PRECISION,INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION,INTENT(OUT) :: SoD, SoV, beta, mu
    DOUBLE PRECISION,PARAMETER :: eps_min = 1D-6
    
    ! Sommerfeld number for rotation according to Butenschoen 1976
    SoD = ( bdat%b / bdat%d )**2 * ABS(eps) / ( 2D0 * ( 1D0 - eps**2 )**2 ) & 
         & * SQRT( pi**2 * ( 1D0 - eps**2 ) + 16D0 * eps**2 ) * bdat%a(1) &
         & * ( ABS(eps) - 1D0 ) / ( bdat%a(2) + ABS(eps) )

    ! Sommerfeld number for displacement according to Butenschoen 1976
    SoV = 4D0 * ( bdat%b / bdat%d )**2 * ( 1D0 - eps**2 )**( -5D0 / 2D0 ) &
         & * ( ( pi / 2D0 - 1D0 / 2D0 * ACOS(eps) ) * ( 1D0 + 2D0 * eps**2 ) &
         & +  3D0 / 2D0 * eps * SQRT( 1D0 - eps**2 ) ) * bdat%a(8) * ( 1D0 - eps ) / ( -bdat%a(9) - eps )

    ! angle between force for rotation and minimum clearance according to Butenschoen 1976
    beta = ATAN2( pi * SQRT( 1D0 - eps**2 ), ( 2D0 * ABS(eps) ) ) &
         & * ( bdat%a(3) + bdat%a(4) * ABS(eps) + bdat%a(5) * eps**2 + bdat%a(6) * ABS(eps)**3 + bdat%a(7) * eps**4 )

    IF (ABS(eps) .LT. eps_min) THEN
       ! avoid division infinite by infinite in case of zero relative eccentricity
       ! use analytical limit of abs_MR for eps going to zero
       mu = HUGE(1D0)
    ELSE
       ! friction coefficient according to Butenschoen
       mu = bdat%Psi * ( ABS( ( omega(1) - omega(2) ) / ( omega(2) + omega(1) - 2D0 * delta_dot ) ) &
            & * pi / ( SQRT( 1D0 - eps**2 ) * SoD ) + SIN(beta) * ABS(eps) / 2D0 )
    ENDIF

  END SUBROUTINE SOMMERFELD_NUMBERS

  SUBROUTINE SOMMERFELD_NUMBERS_EXT(bdat, eps, omega, delta_dot, SoD, SoV, beta, mu)
    TYPE(BEARING_DATA),INTENT(IN) :: bdat
    DOUBLE PRECISION,INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION,INTENT(OUT) :: SoD, SoV, beta, mu
    DOUBLE PRECISION :: epsd, SoDd, SoVd, eps_max

    IF (ABS(eps) .LT. bdat%eps_max) THEN
       ! According to the thesis of Butenschoen those approximations are valid until eps = 0.999
       CALL SOMMERFELD_NUMBERS(bdat, eps, omega, delta_dot, SoD, SoV, beta, mu)
    ELSE
       ! Do a linear extrapolation above eps_max
       eps_max = SIGN(bdat%eps_max, eps)
       epsd = eps - eps_max
       CALL SOMMERFELD_NUMBERS_D(bdat, eps_max, epsd, omega, delta_dot, sod&
            &   , sodd, sov, sovd, beta, mu)
       SoD = SoD + SoDd
       SoV = SoV + SoVd
    ENDIF
  END SUBROUTINE SOMMERFELD_NUMBERS_EXT
  
  SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE(bdat, omega, e, e_dot, k, eps, eps_dot, delta, SoD, SoV, mu, beta) BIND(C)
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
!       d   ... shaft diameter [m]
!       Psi ... relative radial clearance Psi = ( D - d ) / D [1]
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
    !       mu ...             friction coefficient [N/N]
!       beta  ...          angle between force for rotation and minimum clearance [rad]
    TYPE(BEARING_DATA),INTENT(IN) :: bdat
    DOUBLE PRECISION,INTENT(IN) :: omega(2), e(2), e_dot(2)
    DOUBLE PRECISION,INTENT(OUT) :: k(3), eps, eps_dot, delta, SoD, SoV, mu, beta
    DOUBLE PRECISION :: abs_e, abs_e_dot
        DOUBLE PRECISION :: delta_dot, alpha, kappa, phi
        DOUBLE PRECISION :: omega_res, abs_FD, abs_FV, abs_MR
        
    INTRINSIC ATAN2, SQRT, SUM, COS, SIN, ABS, SIGN
        
        ! angle of the position with minimum clearance between shaft and bearing
    delta = ATAN2(e(2), e(1))
        
        ! angle of the velocity vector of the shaft relative to the bearing
    phi = ATAN2(e_dot(2), e_dot(1))
        
        ! angle between velocity vector and minimum clearance
        kappa = phi - delta
        
        ! absolute value of the eccentricity of the shaft inside the bearing
    abs_e = SQRT(SUM( e(:)**2 ))

        ! absolute value of the velocity of the shaft relative to the bearing
    abs_e_dot = SQRT( SUM( e_dot(:)**2 ) )
        
        ! time derivative of the relative eccentricity of the shaft
    eps_dot = 2D0 * COS(kappa) * abs_e_dot / bdat%s;

        ! relative eccentricity of the shaft
    eps = 2D0 * abs_e / bdat%s;
        
        IF ( eps_dot .NE. 0D0 ) THEN 
                ! eps is positive if it's time derivative is positive too
                !       attention the signum function is zero if eps_dot is zero
                !       but eps must not be zero in this case
       eps = SIGN(1D0, eps_dot) * eps
        endif

        ! time derivative of angle of minimum clearance
        if ( abs_e .EQ. 0D0 ) THEN
                ! avoid division by zero
                delta_dot = 0D0
        ELSE
                delta_dot = ( e(1) * e_dot(2) - e(2) * e_dot(1) ) / ( e(2)**2 + e(1)**2 )
        ENDIF
        
    CALL SOMMERFELD_NUMBERS_EXT(bdat, eps, omega, delta_dot, SoD, SoV, beta, mu)
        
        ! effective hydrodynamic angular velocity according to Butenschoen 1976
        omega_res = omega(1) + omega(2) - 2D0 * delta_dot

        ! angle of the force for rotation
    alpha = delta - beta * SIGN(1D0, omega_res)

        ! absolute value of the force for rotation
    abs_FD = SoD * ( bdat%b * bdat%d * bdat%eta * ABS(omega_res) ) / bdat%Psi**2

        ! absolute value of the force for displacement
    abs_FV = SoV * ( bdat%b * bdat%d * bdat%eta * eps_dot ) / bdat%Psi**2
        
    IF (mu .GE. HUGE(1D0)) THEN
       abs_MR = pi * bdat%b * bdat%d**2 * bdat%eta * ABS( omega(1) - omega(2) ) / bdat%Psi / 2D0
        ELSE
       abs_MR = mu * abs_FD * bdat%d / 2D0
        ENDIF

        ! sum of force for rotation and force for displacement
    k(1) = abs_FD * COS(alpha) + abs_FV * COS(delta)
    k(2) = abs_FD * SIN(alpha) + abs_FV * SIN(delta)
        
        ! friction torque
    k(3) = abs_MR * SIGN(1D0, omega(1) - omega(2) )
!-----------------------------------------------------------------------------------------------------------------------------------------------
!       PRINT OUT OF RESULTS
!-----------------------------------------------------------------------------------------------------------------------------------------------
        !~ PRINT *,'HYDRODYNAMIC_PLAIN_BEARING_FORCE'
    !~ PRINT *,'b=',bdat%b
    !~ PRINT *,'d=',bdat%d
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
    !~ PRINT *,'a(1)=',bdat%a(1)
    !~ PRINT *,'a(2)=',bdat%a(2)
    !~ PRINT *,'a(3)=',bdat%a(3)
    !~ PRINT *,'a(4)=',bdat%a(4)
    !~ PRINT *,'a(5)=',bdat%a(5)
    !~ PRINT *,'a(6)=',bdat%a(6)
    !~ PRINT *,'a(7)=',bdat%a(7)
    !~ PRINT *,'a(8)=',bdat%a(8)
    !~ PRINT *,'a(9)=',bdat%a(9)
        !~ PRINT *,'SoD=',SoD
        !~ PRINT *,'SoV=',SoV
        !~ PRINT *,'beta=',beta
    !~ PRINT *,'mu=',mu
        !~ PRINT *,'alpha=',alpha
        !~ PRINT *,'omega_res=',omega_res
        !~ PRINT *,'abs_FD=',abs_FD
        !~ PRINT *,'abs_FV=',abs_FV
        !~ PRINT *,'abs_MR=',abs_MR
        !~ PRINT *, 'k=',k
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE

!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.11 (r5902M) - 15 Dec 2015 09:00
!
        
!  Differentiation of sommerfeld_numbers in forward (tangent) mode:
!   variations   of useful results: sod beta sov mu
!   with respect to varying inputs: eps
!   RW status of diff variables: eps:in sod:out beta:out sov:out
!                mu:out
  SUBROUTINE SOMMERFELD_NUMBERS_D(bdat, eps, epsd, omega, delta_dot, sod&
&   , sodd, sov, sovd, beta, mu)
    TYPE(BEARING_DATA), INTENT(IN) :: bdat
    DOUBLE PRECISION, INTENT(IN) :: eps, delta_dot, omega(2)
    DOUBLE PRECISION, INTENT(IN) :: epsd
    DOUBLE PRECISION, INTENT(OUT) :: sod, sov, beta, mu
    DOUBLE PRECISION, INTENT(OUT) :: sodd, sovd
    DOUBLE PRECISION, PARAMETER :: eps_min = 1d-6
    INTRINSIC ABS
    INTRINSIC SQRT
    INTRINSIC ACOS
    INTRINSIC ATAN2
    INTRINSIC SIN
    INTRINSIC HUGE
    DOUBLE PRECISION :: arg1
    DOUBLE PRECISION :: arg1d
    DOUBLE PRECISION :: result1
    DOUBLE PRECISION :: result1d
    DOUBLE PRECISION :: pwx1
    DOUBLE PRECISION :: pwx1d
    DOUBLE PRECISION :: pwr1
    DOUBLE PRECISION :: pwr1d
    DOUBLE PRECISION :: result2
    DOUBLE PRECISION :: result2d
    DOUBLE PRECISION :: abs1d
    DOUBLE PRECISION :: abs4d
    DOUBLE PRECISION :: abs7d
    DOUBLE PRECISION :: abs0d
    DOUBLE PRECISION :: abs6d
    DOUBLE PRECISION :: abs8
    DOUBLE PRECISION :: abs7
    DOUBLE PRECISION :: abs6
    DOUBLE PRECISION :: abs5
    DOUBLE PRECISION :: abs4
    DOUBLE PRECISION :: abs3
    DOUBLE PRECISION :: abs2
    DOUBLE PRECISION :: abs1
    DOUBLE PRECISION :: abs0
    DOUBLE PRECISION :: abs5d
    DOUBLE PRECISION :: abs8d
    IF (eps .GE. 0.) THEN
      abs0d = epsd
      abs0 = eps
        ELSE
      abs0d = -epsd
      abs0 = -eps
        ENDIF
    IF (eps .GE. 0.) THEN
      abs4d = epsd
      abs4 = eps
        ELSE
      abs4d = -epsd
      abs4 = -eps
        ENDIF
    IF (eps .GE. 0.) THEN
      abs7d = epsd
      abs7 = eps
    ELSE
      abs7d = -epsd
      abs7 = -eps
    END IF
! Sommerfeld number for rotation according to Butenschoen 1976
    arg1d = 16d0*2*eps*epsd - pi**2*2*eps*epsd
    arg1 = pi**2*(1d0-eps**2) + 16d0*eps**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.D0
    ELSE
      result1d = arg1d/(2.0*SQRT(arg1))
    END IF
    result1 = SQRT(arg1)
    sodd = (bdat%a(1)*((bdat%b**2*abs0d*2d0*(1d0-eps**2)**2/bdat%d**2+&
&     bdat%b**2*abs0*2d0*2**2*(1d0-eps**2)*eps*epsd/bdat%d**2)*result1*(&
&     abs4-1d0)/(2d0**2*(1d0-eps**2)**4)+bdat%b**2*abs0*(result1d*(abs4-&
&     1d0)+result1*abs4d)/(bdat%d**2*2d0*(1d0-eps**2)**2))*(bdat%a(2)+&
&     abs7)-bdat%b**2*abs0*result1*bdat%a(1)*(abs4-1d0)*abs7d/(bdat%d**2&
&     *2d0*(1d0-eps**2)**2))/(bdat%a(2)+abs7)**2
    sod = (bdat%b/bdat%d)**2*abs0/(2d0*(1d0-eps**2)**2)*result1*bdat%a(1&
&     )*(abs4-1d0)/(bdat%a(2)+abs7)
! Sommerfeld number for displacement according to Butenschoen 1976
    pwx1d = -(2*eps*epsd)
    pwx1 = 1d0 - eps**2
    IF (pwx1 .GT. 0.0 .OR. (pwx1 .LT. 0.0 .AND. -(5d0/2d0) .EQ. INT(-(&
&       5d0/2d0)))) THEN
      pwr1d = -(5d0*pwx1**((-1)-5d0/2d0)*pwx1d/2d0)
    ELSE IF (pwx1 .EQ. 0.0 .AND. -(5d0/2d0) .EQ. 1.0) THEN
      pwr1d = pwx1d
    ELSE
      pwr1d = 0.0
    END IF
    pwr1 = pwx1**(-(5d0/2d0))
    IF (eps .EQ. 1.0 .OR. eps .EQ. (-1.0)) THEN
      result1d = 0.D0
    ELSE
      result1d = -(epsd/SQRT(1.0-eps**2))
    END IF
    result1 = ACOS(eps)
    arg1d = -(2*eps*epsd)
    arg1 = 1d0 - eps**2
    IF (arg1 .EQ. 0.0) THEN
      result2d = 0.D0
    ELSE
      result2d = arg1d/(2.0*SQRT(arg1))
    END IF
    result2 = SQRT(arg1)
    sovd = (4d0*bdat%b**2*bdat%a(8)*((pwr1d*(1d0-eps)-pwr1*epsd)*((pi/&
&     2d0-1d0/2d0*result1)*(1d0+2d0*eps**2)+3d0/2d0*eps*result2)+pwr1*(&
&     1d0-eps)*((pi/2d0-1d0/2d0*result1)*2d0*2*eps*epsd-result1d*(1d0+&
&     2d0*eps**2)/2d0+3d0*(epsd*result2+eps*result2d)/2d0))*(-bdat%a(9)-&
&     eps)/bdat%d**2+4d0*bdat%b**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+&
&     2d0*eps**2)+3d0/2d0*eps*result2)*bdat%a(8)*(1d0-eps)*epsd/bdat%d**&
&     2)/(-bdat%a(9)-eps)**2
    sov = 4d0*(bdat%b/bdat%d)**2*pwr1*((pi/2d0-1d0/2d0*result1)*(1d0+2d0&
&     *eps**2)+3d0/2d0*eps*result2)*bdat%a(8)*(1d0-eps)/(-bdat%a(9)-eps)
    IF (eps .GE. 0.) THEN
      abs1d = epsd
      abs1 = eps
    ELSE
      abs1d = -epsd
      abs1 = -eps
    END IF
    IF (eps .GE. 0.) THEN
      abs5d = epsd
      abs5 = eps
    ELSE
      abs5d = -epsd
      abs5 = -eps
    END IF
    IF (eps .GE. 0.) THEN
      abs8d = epsd
      abs8 = eps
    ELSE
      abs8d = -epsd
      abs8 = -eps
    END IF
! angle between force for rotation and minimum clearance according to Butenschoen 1976
    arg1d = -(2*eps*epsd)
    arg1 = 1d0 - eps**2
    IF (arg1 .EQ. 0.0) THEN
      result1d = 0.D0
    ELSE
      result1d = arg1d/(2.0*SQRT(arg1))
    END IF
    result1 = SQRT(arg1)
    beta = ATAN2(pi*result1, 2d0*abs1)*(bdat%a(3)+bdat%a(4)*abs5+bdat%a(&
&     5)*eps**2+bdat%a(6)*abs8**3+bdat%a(7)*eps**4)
    IF (eps .GE. 0.) THEN
      abs2 = eps
    ELSE
      abs2 = -eps
    END IF
    IF (abs2 .LT. eps_min) THEN
! avoid division infinite by infinite in case of zero relative eccentricity
! use analytical limit of abs_MR for eps going to zero
      mu = HUGE(1d0)
    ELSE
      IF ((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot) .GE. 0.&
&     ) THEN
        abs3 = (omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot)
      ELSE
        abs3 = -((omega(1)-omega(2))/(omega(2)+omega(1)-2d0*delta_dot))
      END IF
      IF (eps .GE. 0.) THEN
        abs6d = epsd
        abs6 = eps
      ELSE
        abs6d = -epsd
        abs6 = -eps
      END IF
! friction coefficient according to Butenschoen
      arg1d = -(2*eps*epsd)
      arg1 = 1d0 - eps**2
      IF (arg1 .EQ. 0.0) THEN
        result1d = 0.D0
      ELSE
        result1d = arg1d/(2.0*SQRT(arg1))
      END IF
      result1 = SQRT(arg1)
      mu = bdat%psi*(abs3*pi/(result1*sod)+SIN(beta)*abs6/2d0)
    END IF
  END SUBROUTINE SOMMERFELD_NUMBERS_D
END MODULE HYDRODYNAMIC_PLAIN_BEARING
