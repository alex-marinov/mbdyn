! $Header$
! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2013
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
MODULE MATH_CONST
        IMPLICIT NONE
        DOUBLE PRECISION,PARAMETER :: pi = 3.14159265358979D0
END MODULE MATH_CONST

SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE(b,d,Psi,eta,omega,e,e_dot,k,eps,eps_dot,delta,SoD,SoV,my,beta)
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

        use MATH_CONST
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN) :: b, d, Psi, eta, omega(2), e(2), e_dot(2)
        DOUBLE PRECISION,INTENT(OUT) :: k(3), eps, eps_dot, delta, SoD, SoV, my, beta
        DOUBLE PRECISION :: abs_e, abs_e_dot, s
        DOUBLE PRECISION :: delta_dot, alpha, kappa, phi
        DOUBLE PRECISION :: omega_res, abs_FD, abs_FV, abs_MR
        DOUBLE PRECISION :: a1, a2, a3, a4, a5, a6, a7, a8, a9
        INTRINSIC DATAN2, DSQRT, SUM, DCOS, DSIN, DABS, DSIGN
        
        ! clearance of the shaft inside the bearing
        s = d * Psi / ( 1 - Psi )
        
        ! angle of the position with minimum clearance between shaft and bearing
        delta = DATAN2(e(2),e(1))
        
        ! angle of the velocity vector of the shaft relative to the bearing
        phi = DATAN2(e_dot(2),e_dot(1))
        
        ! angle between velocity vector and minimum clearance
        kappa = phi - delta
        
        ! absolute value of the eccentricity of the shaft inside the bearing
        abs_e = DSQRT( SUM( e(:)**2 ) )

        ! absolute value of the velocity of the shaft relative to the bearing
        abs_e_dot = DSQRT( SUM( e_dot(:)**2 ) )
        
        ! time derivative of the relative eccentricity of the shaft
        eps_dot = 2D0 * DCOS(kappa) * abs_e_dot / s;

        ! relative eccentricity of the shaft
        eps= 2D0 * abs_e / s;
        
        IF ( eps_dot .NE. 0D0 ) THEN 
                ! eps is positive if it's time derivative is positive too
                !       attention the signum function is zero if eps_dot is zero
                !       but eps must not be zero in this case
                eps = DSIGN(1D0,eps_dot) * eps
        endif

        ! time derivative of angle of minimum clearance
        if ( abs_e .EQ. 0D0 ) THEN
                ! avoid division by zero
                delta_dot = 0D0
        ELSE
                !~ delta_dot = DSIN(kappa) * abs_e_dot / abs_e;
                delta_dot = ( e(1) * e_dot(2) - e(2) * e_dot(1) ) / ( e(2)**2 + e(1)**2 )
        ENDIF
        
        ! coefficients for SoD for rotation according to Butenschoen
        a1 = 1.1642D0 - 1.9456D0 * ( b / d ) +7.1161D0 * ( b / d )**2 - 10.1073D0 * ( b / d )**3 &
              &   +5.0141D0 * ( b / d )**4
              
        a2 = -1.000026D0 - 0.023634D0 * ( b / d ) - 0.4215D0 * ( b / d )**2 &
              &   - 0.038817D0 * ( b / d )**3 - 0.090551D0 * ( b / d )**4

        ! coefficients for beta for rotation
        a3 = 1.152624D0 - 0.104565D0 * ( b / d )
        
        a4 = -2.5905D0 + 0.798745D0 * ( b / d )
        
        a5 = 8.73393D0 - 2.3291D0 * ( b / d )
        
        a6 = -13.3414D0 + 3.424337D0 * ( b / d )
        
        a7 = 6.6294D0 - 1.591732D0 * ( b / d )
      
        ! coefficients for displacement
        a8 = 0.70038D0+3.2415D0 * ( b / d ) - 12.2486D0 * ( b / d )**2 + 18.895D0 * ( b / d )**3 &
                &   - 9.3561D0 * ( b / d )**4
                
        a9 = -0.999935D0 + 0.0157434D0 * ( b / d ) - 0.74224D0 * ( b / d )**2 &
                &  + 0.42278D0 * ( b / d )**3 - 0.368928 * ( b / d )**4
        
        ! Sommerfeld number for rotation according to Butenschoen 1976
        SoD = ( b / d )**2 * DABS(eps) / ( 2D0 * ( 1D0 - eps**2 )**2 ) & 
                & * DSQRT( pi**2 * ( 1D0 - eps**2 ) + 16D0 * eps**2 ) * a1 &
                & * ( DABS(eps) - 1D0 ) / ( a2 + DABS(eps) )

        ! Sommerfeld number for displacement according to Butenschoen 1976
        SoV= 4D0 * ( b / d )**2 * ( 1D0 - eps**2 )**( -5D0 / 2D0 ) &
                & * ( ( pi / 2D0 - 1D0 / 2D0 * DACOS(eps) ) * ( 1D0 + 2D0 * eps**2 ) &
                & +  3D0 / 2D0 * eps * DSQRT( 1D0 - eps**2 ) ) * a8 * ( 1D0 - eps ) / ( -a9 - eps )

        ! angle between force for rotation and minimum clearance according to Butenschoen 1976
        beta= ATAN( pi * DSQRT( 1D0 - eps**2 ) / ( 2D0 * DABS(eps) ) ) &
                & * ( a3 + a4 * DABS(eps) + a5 * eps**2 + a6 * DABS(eps)**3 + a7 * eps**4 )

        
        ! friction coefficient according to Butenschoen
        my = Psi * ( DABS( ( omega(1) - omega(2) ) / ( omega(2) + omega(1) - 2D0 * delta_dot ) ) &
                  & * pi / ( DSQRT( 1D0 - eps**2 ) * SoD ) + DSIN(beta) * DABS(eps) / 2D0 )
        
        ! effective hydrodynamic angular velocity according to Butenschoen 1976
        omega_res = omega(1) + omega(2) - 2D0 * delta_dot

        ! angle of the force for rotation
        alpha = delta - beta * DSIGN(1D0,omega_res)

        ! absolute value of the force for rotation
        abs_FD = SoD * ( b * d * eta * DABS(omega_res) ) / Psi**2

        ! absolute value of the force for displacement
        abs_FV = SoV * ( b * d * eta * eps_dot ) / Psi**2
        
        IF ( DABS(eps) .LT. 1D-6 ) THEN
                ! avoid division infinite by infinite in case of zero relative eccentricity
                ! use analytical limit of abs_MR for eps going to zero
                abs_MR = pi * b * d**2 * eta * DABS( omega(1) - omega(2) ) / Psi / 2D0
        ELSE
                abs_MR = my * abs_FD * d / 2D0
        ENDIF

        ! sum of force for rotation and force for displacement
        k(1) = abs_FD * DCOS(alpha) + abs_FV * DCOS(delta)
        k(2) = abs_FD * DSIN(alpha) + abs_FV * DSIN(delta)
        
        ! friction torque
        k(3) = abs_MR * DSIGN(1D0, omega(1) - omega(2) )
!-----------------------------------------------------------------------------------------------------------------------------------------------
!       PRINT OUT OF RESULTS
!-----------------------------------------------------------------------------------------------------------------------------------------------
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

        
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_FORCE

DOUBLE PRECISION FUNCTION SIGN_MOD(x,Delta)
        USE MATH_CONST
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN) :: x, Delta
        INTRINSIC DSIGN, DTANH
        
        IF ( Delta .EQ. 0D0 ) THEN
                SIGN_MOD = DSIGN(1D0,x)
        ELSE
                SIGN_MOD = DTANH( 2 * pi * x / Delta )
        ENDIF
END FUNCTION SIGN_MOD

SUBROUTINE PLAIN_BEARING_CONTACT_FORCE( &
& d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot, e, e_dot, k)
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
!       phenomenological contact force calculation between shaft and bearing with penalty method 
!       phenomenological damping force with Lehr's damping law
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
!       INPUT PARAMETERS
!---------------------------------------------------------------------------------------------------------------------------------------------------------------------
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
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------
!       OUTPUT PARAMETERS
!----------------------------------------------------------------------------------------------------------------------------------------------------------------------
!   k(1)  = Fy [N]        force at the bearing in y direction
!   k(2)  = Fz [N]        force at the bearing in z direction
!   k(3)  = Mx [Nm]       torque at the bearing around the x axis
        
        USE MATH_CONST
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(IN) :: d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega, Phi_dot(2), e(2), e_dot(2)
        DOUBLE PRECISION, INTENT(OUT) :: k(3)
        DOUBLE PRECISION :: eps, eps_dot, delta, phi, kappa, abs_e, abs_e_dot, s
        DOUBLE PRECISION :: abs_FPrs, abs_FPrd, abs_FPr, FP23r(2), delta_Phi_dot, alpha, sign_MP
        DOUBLE PRECISION :: FP23t(2), abs_MP
        INTRINSIC SQRT, EXP, COS, SIN, ABS
        EXTERNAL SIGN_MOD
        DOUBLE PRECISION SIGN_MOD

        ! clearance between shaft and bearing
        s  = d * Psi / ( 1 - Psi )

        ! angle at minimum clearance
        delta = DATAN2(e(2),e(1))
        
        ! moving direction of the shaft with respect to the bearing
        phi = DATAN2(e_dot(2),e_dot(1))
        
        ! angle between moving direction and minimum clearance
        kappa = phi - delta
        
        ! absolute value of the eccentricity of the shaft
        abs_e = DSQRT( SUM( e(:)**2 ) )

        ! absolute value of the velocity of the shaft
        abs_e_dot = DSQRT( SUM( e_dot(:)**2 ) )
        
        ! derivative of the relative eccentricity versus time
        eps_dot = 2D0 * DCOS(kappa) * abs_e_dot / s;

        ! relative eccentricity
        eps = 2D0 * abs_e / s;

        ! stiffness contact force
        abs_FPrs = abs_FPrs1 * DEXP( d * Psi / 2D0 * sP / abs_FPrs1 * ( eps - 1D0 ) )
        
        ! damping contact force
        abs_FPrd = d * Psi * DSQRT( m * sP ) * DL * eps_dot * DEXP( d * Psi / 4D0 * sP / abs_FPrs1 * ( eps - 1D0 ) )

        ! contact force
        abs_FPr = abs_FPrs + abs_FPrd 

        ! radial contact force
        FP23r(1) = abs_FPr * DCOS(delta) 
        FP23r(2) = abs_FPr * DSIN(delta)

        ! The difference between the angular velocity of the shaft and the bearing is related to the circumference of the shaft.
        ! The term 1 + Psi * eps is needed because the sign of the torque
        ! must be related to the difference of the circumferential velocity.
        delta_Phi_dot = Phi_dot(1) - ( 1D0 + Psi * eps  / ( 1 - Psi ) ) * Phi_dot(2)

        ! angle between radial and tangential component of the the contact force
        ! 	Attention: the signum function must not be used here since this would lead to a wrong direction of the force
        ! 	The direction of the tangential component of the contact force must be either pi/2 oder -pi/2
        IF ( delta_Phi_dot .GT. 0D0 ) THEN
                alpha = pi / 2D0
        ELSE
                alpha = -pi / 2D0
        ENDIF

        ! sign of the friction torque
        ! 	For numerical reasons the modified signum function must be used here
        sign_MP = SIGN_MOD( delta_Phi_dot, signum_delta_omega )

        ! vector of the tangential component of the contact force
        !	The term abs(sign_MP) is needed because the tangential force and the torque must be consistent
        FP23t(1) = DABS(sign_MP) * myP * abs_FPr * DCOS( delta + alpha )
        FP23t(2) = DABS(sign_MP) * myP * abs_FPr * DSIN( delta + alpha )
        
        ! absolute value of the friction torque
        abs_MP = myP * abs_FPr * d / 2D0

        ! the radial and tangential component of the force are applied to the center of the bearing
        k(1:2) = FP23r + FP23t
        
        ! torque of the tangential component of the contact force around the x axis
        k(3) = sign_MP * abs_MP * ( 1D0 + Psi * eps / ( 1 - Psi ) )
        

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
END SUBROUTINE PLAIN_BEARING_CONTACT_FORCE

SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_AND_CONTACT_FORCE( &
& b,d,Psi,eta,omega,e,e_dot,  sP, DL, m, abs_FPrs1, myP, signum_delta_omega, k,eps,eps_dot,delta,SoD,SoV,my,beta )
        IMPLICIT NONE
        DOUBLE PRECISION,INTENT(IN) :: b, d, Psi, eta, omega(2), e(2), e_dot(2)
        double PRECISION, INTENT(IN) :: sP, DL, m, abs_FPrs1, myP, signum_delta_omega
        DOUBLE PRECISION,INTENT(OUT) :: k(3), eps, eps_dot, delta, SoD, SoV, my, beta
        DOUBLE PRECISION :: k_hydro(3), k_contact(3)
        EXTERNAL HYDRODYNAMIC_PLAIN_BEARING_FORCE, PLAIN_BEARING_CONTACT_FORCE
        
        CALL HYDRODYNAMIC_PLAIN_BEARING_FORCE(b,d,Psi,eta,omega,e,e_dot,k_hydro,eps,eps_dot,delta,SoD,SoV,my,beta)
        CALL PLAIN_BEARING_CONTACT_FORCE( d, Psi, sP, DL, m, abs_FPrs1, myP, signum_delta_omega,omega,e,e_dot,k_contact)
        
        k = k_hydro + k_contact
END SUBROUTINE HYDRODYNAMIC_PLAIN_BEARING_AND_CONTACT_FORCE

