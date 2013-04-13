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
! Interface to NREL's AeroDyn library

SUBROUTINE MBDyn_init( Version, NBlades, RotRadius )

USE Identify
USE Blade
USE Wind  ! Use Global variable VX,VY,VZ.
USE Precision


IMPLICIT NONE

CHARACTER(26) Version
INTEGER(4) NBlades
REAL(ReKi) RotRadius

        DynProg = 'MBDyn '
        DynVer = Version
        Prog = ' '

        CALL SetProgName

        NB = NBlades

        B = NB    ! Need this value when we swith on "WAKE"
                  ! or "SWIRL" options
        R = RotRadius
        RETURN

END SUBROUTINE MBDyn_init

SUBROUTINE MBDyn_ad_inputgate( FileName, FileNameLen, ElemFileName, ElemFileNameLen, * )

IMPLICIT NONE

INTEGER(4) FileNameLen,ElemFileNameLen
CHARACTER(FileNameLen) FileName
CHARACTER(ElemFileNameLen) ElemFileName

        !print *,'### FileName=',FileName(1:FileNameLen),' FileNameLen=',FileNameLen

        !print *,'==> MBDyn_ad_inputgate(',FileName(1:FileNameLen),',',ElemFileName(1:ElemFileNameLen),')'

        if (FileNameLen .gt. 0) then
                CALL AD_inputgate(FileName)
        else
                CALL ADInputGate
        endif

        if (ElemFileNameLen .gt. 0) then
                !CALL ElemOpen('./aerodyn.elm')
                CALL ElemOpen(ElemFileName)
        endif

        !print *,'<== MBDyn_ad_inputgate(',FileName(1:FileNameLen),',',ElemFileName(1:ElemFileNameLen),')'

        RETURN 0

END SUBROUTINE MBDyn_ad_inputgate

SUBROUTINE MBDyn_true( val )

IMPLICIT NONE

LOGICAL val

        val = .TRUE.

        RETURN

END SUBROUTINE MBDyn_true

SUBROUTINE MBDyn_false( val )

IMPLICIT NONE

LOGICAL val

        val = .FALSE.

        RETURN

END SUBROUTINE MBDyn_false

! This subroutine is to make MBDyn-AeroDyn interface function
! can access the AeroDyn variables which was defined in the
! common module. By Fanzhong MENG 21 Feb. 2008

SUBROUTINE MBDyn_com_data( c_blade, c_elem )

USE Identify
USE Blade
USE Element

IMPLICIT NONE

INTEGER(4) c_blade
INTEGER(4) c_elem

        IBlade = c_blade
        JElem  = c_elem
        
        RETURN
END SUBROUTINE MBDyn_com_data

! This subroutine is to pass the current simulation time 
! of MBDyn to AeroDyn!
! c_time: current time

SUBROUTINE MBDyn_sim_time(c_time)

USE Identify
USE Blade
USE Element
USE Precision
USE AeroTime

IMPLICIT NONE

REAL(DbKi) c_time

        TIME = c_time
        
        RETURN
END SUBROUTINE MBDyn_sim_time

! This subroutine is to pass the current simulation time step 
! of MBDyn to AeroDyn!
! dt: time step

SUBROUTINE MBDyn_time_step(time_step)

USE Identify
USE Blade
USE Element
USE Precision
USE AeroTime

IMPLICIT NONE

REAL(ReKi) time_step

        DT = time_step
        
        RETURN
END SUBROUTINE MBDyn_time_step

!Tip loss constants are calculated and stored.

SUBROUTINE MBDyn_get_tl_const ( RLOCAL, Cur_elem )

! AeroDyn Modules:

USE               Blade
USE               Element
USE               Precision


IMPLICIT          NONE


! Passed Variables:

REAL(4)         RLOCAL
INTEGER(4)      Cur_elem


! Local Variables:

REAL(4)         DTIP
!!! INTEGER(4)      IELM

! Calculation of tip loss constants. 

! R = RLOCAL(NELM) + 0.5 * DR( NELM ) !* COS( PC )

! Calculate the tip-loss constant for each element
   DTIP           = R - RLOCAL
   TLCNST( Cur_elem ) = 0.5 * B * DTIP / RLOCAL

RETURN
END SUBROUTINE MBDyn_get_tl_const




!  Hub loss constants are calculated and stored.
SUBROUTINE MBDyn_get_hl_const( RLOCAL, Cur_elem, RHub )

! AeroDyn Modules:

USE               Blade
USE               Element
USE               Precision


IMPLICIT          NONE

! Passed variables:

REAL(4)        RLOCAL
INTEGER(4)     Cur_elem
REAL(4)        RHub

! Local Variables:

REAL(4)        DHUB
!!! INTEGER(4)     IELM


! Calculation of hub loss constants. 

! Calculate the hub-loss constant for each element
IF (RHub > 0.001) THEN
      DHUB           = RLOCAL - RHub
      HLCNST( Cur_elem ) = 0.5 * B * DHUB / RHub
ENDIF



RETURN
END SUBROUTINE MBDyn_get_hl_const

