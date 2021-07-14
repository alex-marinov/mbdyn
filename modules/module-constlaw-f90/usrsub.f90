! $Header$
! MBDyn (C) is a multibody analysis code. 
! http://www.mbdyn.org
! 
! Copyright (C) 1996-2017
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
! generic helpers
      SUBROUTINE USINIT(D, N, V, ERR)
      INTEGER D, N, ERR
      DOUBLE PRECISION V(N)
      IF (N .NE. D) THEN
        ERR = 1
        RETURN
      ENDIF
      ERR = 0
      RETURN
      END

      SUBROUTINE USUPDT(D, N, V, E, EP, F, FDE, FDEP, ERR)
      INTEGER D, N, ERR
      DOUBLE PRECISION V(N), E(D), EP(D), F(D), FDE(D, D), FDEP(D, D)
      INTEGER R, C
      DO R = 1,D
        DO C = 1,D
           FDE(R, C) = 0.0
           FDEP(R, C) = 0.0
        ENDDO
        FDE(R, R) = V(R)
        F(R) = V(R)*E(R)
      ENDDO
      ERR = 0
      RETURN
      END

! 1D elastic constitutive law
      SUBROUTINE US1INIT(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      CALL USINIT(1, N, V, ERR)
      RETURN
      END

      SUBROUTINE US1DSTR(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      ERR = 0
      RETURN
      END

      SUBROUTINE US1UPDT(N, V, E, EP, F, FDE, FDEP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E, EP, F, FDE, FDEP
      CALL USUPDT(1, N, V, E, EP, F, FDE, FDEP, ERR)
      RETURN
      END

      SUBROUTINE US1AFTC(N, V, E, EP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E, EP
      ERR = 0
      RETURN
      END

! 3D elastic constitutive law
      SUBROUTINE US3INIT(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      CALL USINIT(3, N, V, ERR)
      RETURN
      END

      SUBROUTINE US3DSTR(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      ERR = 0
      RETURN
      END

      SUBROUTINE US3UPDT(N, V, E, EP, F, FDE, FDEP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E(3), EP(3), F(3), FDE(3, 3), FDEP(3, 3)
      CALL USUPDT(3, N, V, E, EP, F, FDE, FDEP, ERR)
      RETURN
      END

      SUBROUTINE US3AFTC(N, V, E, EP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E, EP
      ERR = 0
      RETURN
      END

! 6D elastic constitutive law
      SUBROUTINE US6INIT(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      CALL USINIT(6, N, V, ERR)
      RETURN
      END

      SUBROUTINE US6DSTR(N, V, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N)
      ERR = 0
      RETURN
      END

      SUBROUTINE US6UPDT(N, V, E, EP, F, FDE, FDEP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E(6), EP(6), F(6), FDE(6, 6), FDEP(6, 6)
      CALL USUPDT(6, N, V, E, EP, F, FDE, FDEP, ERR)
      RETURN
      END

      SUBROUTINE US6AFTC(N, V, E, EP, ERR)
      INTEGER N, ERR
      DOUBLE PRECISION V(N), E, EP
      ERR = 0
      RETURN
      END

