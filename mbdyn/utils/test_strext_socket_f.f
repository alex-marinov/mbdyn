C $Header$ */
C MBDyn (C) is a multibody analysis code. 
C http://www.mbdyn.org
C 
C Copyright (C) 1996-2014
C 
C Pierangelo Masarati	<masarati@aero.polimi.it>
C Paolo Mantegazza	<mantegazza@aero.polimi.it>
C 
C Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
C via La Masa, 34 - 20156 Milano, Italy
C http://www.aero.polimi.it
C 
C Changing this copyright notice is forbidden.
C 
C This program is free software; you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation (version 2 of the License).
C 
C 
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C 
C You should have received a copy of the GNU General Public License
C along with this program; if not, write to the Free Software
C Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C 
C NOTE: this is all we use from MBDyn.  It was intentionally designed
C to be configuration-independent.
C
C #include "mbc.h"

      SUBROUTINE MAINF

      IMPLICIT NONE

      INTEGER*4 REFNODE, NODES, ROT, ITERS, VERB
      INTEGER*4 STEPS, KEEPGOING, ITER, RC, I, J, N, CONVERGED
      REAL*4 RF(3), RM(3), NF(3, 100), NM(3, 100)
      REAL*4 RX(3), RR(3, 3), RTHETA(3), RXP(3), ROMEGA(3),
     & NX(3, 100), NR(3, 3, 100), NTHETA(3, 100),
     & NXP(3, 100), NOMEGA(3, 100)

      EQUIVALENCE(RR(1, 1), RTHETA(1))
      EQUIVALENCE(NR(1, 1, 1), NTHETA(1, 1))

      CALL TDATA(REFNODE, NODES, ROT, ITERS, VERB, RC)
      IF (NODES .GT. 100) THEN
        WRITE(*, *) 'NODES=',NODES,' exceeds max (100)'
        STOP
      ENDIF

      KEEPGOING = 1
      STEPS = 0
      DO WHILE (KEEPGOING .EQ. 1)
        CONVERGED = 0
        DO ITER = 1,ITERS
          CALL TRECV(RX, RR, RXP, ROMEGA, NX, NR, NXP, NOMEGA, RC)
          IF (RC .NE. 0) THEN
            WRITE (*, *) 'recv failed'
            STOP
          ENDIF

          IF (VERB .NE. 0) THEN
            IF (REFNODE .NE. 0) THEN
              WRITE (*, *) 'reference node:'
              WRITE (*, *) 'x=', (RX(I), I=1,3)
              IF (ROT .EQ. 0) THEN
                WRITE (*, *) 'R=', (RR(I,J), I,J = 1,3)
              ELSEIF (ROT .EQ. 1) THEN
                WRITE (*, *) 'THETA=', (RTHETA(I), I = 1,3)
              ELSEIF (ROT .EQ. 2) THEN
                WRITE (*, *) 'EULER123=', (RTHETA(I), I = 1,3)
              ENDIF
              WRITE (*, *) 'xp=', (RXP(I), I = 1,3)
              WRITE (*, *) 'omega=', (ROMEGA(I), I = 1,3)
            ENDIF

            IF (NODES .GT. 0) THEN
              DO N = 1,NODES
                WRITE (*, *) 'node', N, ':'
                WRITE (*, *) 'x=', (NX(I,N), I = 1,3)
                IF (ROT .EQ. 0) THEN
                  WRITE (*, *) 'R=', (NR(I,J,N), I,J = 1,3)
                ELSEIF (ROT .EQ. 1) THEN
                  WRITE (*, *) 'THETA=', (NTHETA(I, N), I = 1,3)
                ELSEIF (ROT .EQ. 2) THEN
                  WRITE (*, *) 'EULER123=', (NTHETA(I, N), I = 1,3)
                ENDIF
                WRITE (*, *) 'xp=', (NXP(I,N), I = 1,3)
                WRITE (*, *) 'omega=', (NOMEGA(I,N), I = 1,3)
              ENDDO
            ENDIF
          ENDIF

          CALL TFORCE(RF, RM, NF, NM)
          IF (ITER .EQ. ITERS) THEN
            CONVERGED = 1
          ENDIF
          CALL TSEND(RF, RM, NF, NM, CONVERGED, RC)
          IF (RC .NE. 0) THEN
            WRITE(*, *) 'send failed'
            STOP
          ENDIF
        ENDDO
      ENDDO
      END SUBROUTINE
