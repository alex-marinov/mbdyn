      PROGRAM WHEEL      

      IMPLICIT NONE

      INTEGER*4 MAXPTF
      PARAMETER(MAXPTF=10)

      INTEGER*4 IRCW, IPC, NPTF
      REAL*8 OME, RM, OMEP, BETA, BETAP, BLENG
      REAL*8 VM(6), CC(9), RES(6), TABFAT(2*MAXPTF), PP(6)

C
C inizializzo dati costanti
C
C raggio ruota
      CC(7) = 1.
C semilarghezze ruota
      CC(8) = -.5
      CC(9) = .5

      

C coseni direttori asse ruota
      CC(1) = 1.
      CC(2) = 0.
      CC(3) = 0.

C posizione riferimento su asse ruota
      CC(4) = 0.
      CC(5) = 0.
      CC(6) = 1.

      CALL WVEFR(VW, OME, CC, RES, RM, IRCW, IPC, OMEP, TABFAT,
     &  NPFT, BETA, BETAP, BLENG, PP)

      STOP
      END     
