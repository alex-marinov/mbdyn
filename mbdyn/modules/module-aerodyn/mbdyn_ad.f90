! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2008
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

SUBROUTINE MBDyn_init( Version, NBlades )

USE Identify
USE Blade
USE Wind  ! Use Global variable VX,VY,VZ.


IMPLICIT NONE

CHARACTER(26) Version
INTEGER(4) NBlades

        DynProg = 'MBDyn '
        DynVer = Version
        Prog = ' '

        CALL SetProgName

        NB = NBlades

        RETURN

END SUBROUTINE MBDyn_init

SUBROUTINE MBDyn_ad_inputgate( FileName, FileNameLen )

IMPLICIT NONE

INTEGER(4) FileNameLen
CHARACTER(FileNameLen) FileName

        print *,'### FileName=',FileName(1:FileNameLen),' FileNameLen=',FileNameLen

        CALL AD_inputgate(FileName)

        RETURN

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




