C $Header$
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
C Copyright (C) 1985-2014 GRAALL
C
C Gian Luca Ghiringhelli	<ghiringhelli@aero.polimi.it>
C Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
C via La Masa, 34 - 20156 Milano, Italy
C http://www.aero.polimi.it
C 
C ======================================================================
C ======================================================================
C                   25-01-93                      SUBROUTINE DMPFR2
C ======================================================================
C ======================================================================
      SUBROUTINE DMPFR(DL,DLDT,V,
     &  TBDMR,TBDMA,TBCVR,TBCVA,VETVIS,
     &  NPDMR,NPDMA,NPCVR,NPCVA,NPVT,
     &  FTOT,FVISDL,FELDL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
      DOUBLE PRECISION DL,DLDT,V(*),
     &  TBDMR(*),TBDMA(*),TBCVR(*),TBCVA(*),VETVIS(*),
     &  FTOT,FVISDL,FELDL
      INTEGER NPDMR,NPDMA,NPCVR,NPCVA,NPVT
C
      INTEGER RLA,PNTR,P0,A0,RL0,PLTREX,TCOST,TLIN,TEMPORIF,
     &  VELRIF,CORONA,TRIFVUOTO,T,TRIFTEMP,RNURIF,ASTNT,RO,TIME
      PARAMETER(
     &  RLA=1,
     &  PNTR=2,
     &  P0=3,
     &  A0=4,
     &  RL0=5,
     &  PLTREX=6,
     &  TCOST=7,
     &  TLIN=8,
     &  TEMPORIF=9,
     &  VELRIF=10,
     &  CORONA=11,
     &  TRIFVUOTO=12,
     &  T=13,
     &  TRIFTEMP=14,
     &  RNURIF=15,
     &  ASTNT=16,
     &  RO=17,
     &  TIME=18)
C
c     RLA        lunghezza iniziale ammortizzatore
c     PNTR       coeff. di interazione corsa
c     P0         pressione di precarica 
c     A0         sezione gas
c     RL0        lunghezza iniziale camera gas
c	  PLTREX     esponente politropica
c     tcost,tlin legge lineare coeff attrito con corsa
c     temporif   tempo di riferimento per profilo attrito
c     velrif     velocita' di riferimento per profilo attrito    
c     corona     area corona per effetto risucchio
c     trifvuoto  tempo di riferimento per profilo attrito
c     T          temperatura olio in passaggio
c     triftemp   tempo di riferimento per effetto temperatura olio
c     RNURIF     cvis di riferimento
c     SWITCH     0 ammortizzatore passivo, >0 modo di controllo
c     ASTNT      sezione stantuffo olio
c	  RO         densita' olio

C     TIME    tempo dall'inizio della compressione dell'ammortizzatore

c     TBDMR   tabella sezione trafilamento ritorno  NPDMR punti
c     TBDMA   tabella sezione trafilamento andata   NPDMA punti
c     TBCVR   tabella coeff. efflusso ritorno       NPCVR punti
c     TBCVA   tabella coeff. efflusso ritorno       NPCVA punti
c     VETVIS  tabella viscosita'/temperatura        NPVT punti
      
C     ===================================
C     =    CALCOLA LA FORZA ELASTICA    =
C     ===================================
C
C
      FELAS=V(P0)*V(A0)*(V(RL0)/(V(RL0)+DL*V(PNTR)))**V(PLTREX)
C
C Rigidezza (parte principale)
      FELDL=-V(P0)*V(A0)*V(PNTR)*V(PLTREX)/V(RL0)*
     &  (V(RL0)/(V(RL0)+DL*V(PNTR)))**(V(PLTREX)+1.)

c      formula originale attrito
c      fatt=-FELAS*FRICO*TANH(DLDT/VREF))
c      attrito non profilato con tanh: problema all'inverisone
c      Fatt=-FELAS*FRICO*sign(1.,DLDT)

      FATT1=-FELAS*(V(TCOST)+V(TLIN)*DL)*    
     *             (1.-TANH((V(TIME))/V(TEMPORIF)))*
     *             TANH(DLDT/V(VELRIF))
      FATT2=V(P0)*V(CORONA)*(1.-TANH((V(TIME))/V(TRIFVUOTO)))
C ????
      FATT=FATT1+FATT2
C      fatt2=p0*corona*(1-tanh((TIME)/trifvuoto)) c fatt=fatt1+fatt2

C     =====================================================
C     = CALCOLA IL FATTORE CORRETTIVO DEL COEFFICIENTE DI =
C     =   EFFLUSSO DOVUTO ALL'EFFETTO DELLA TEMPERATURA   =
C     =====================================================
      COVERO=RINT(V(T),VETVIS,NPVT)/V(RNURIF)
      RIT=1./COVERO
      CORR=COVERO*(RIT-(RIT-1.)*TANH(2.*(V(TIME)-V(TRIFTEMP))))

C    =============================================================
C    = CALCOLA AREA PER CARRELLO PASSIVO (dldt<0 = compressione) =
C    =============================================================
      IF(ABS(DLDT).LT.1.E-10) DLDT=0.
cactv      IF(SWITCH.EQ.0) THEN
        IF(DLDT.GT.0.) THEN
          AREA=RINT(DL,TBDMR,NPDMR)
        ELSE IF (DLDT.LE.0.) THEN
          AREA=RINT(DL,TBDMA,NPDMA)
        ENDIF
cactv      ELSE
cactv        AREA=ATRAF
cactv      ENDIF
C    ========================================
C    = CALCOLA Ceff (dldt<0 = compressione) =
C    ========================================
      IF (DLDT.GT.0.) THEN
        CEFF=RINT(AREA,TBCVR,NPCVR)
      ELSE IF (DLDT.LE.0.) THEN
        CEFF=RINT(AREA,TBCVA,NPCVA)
      ENDIF
      CEFF=CEFF*CORR
C    =====================================
C    =      CALCOLA LA FORZA VISCOSA     =
C    =====================================
C
C Smorzamento (parte principale?)
      FVISDL=-(V(ASTNT)**3)*ABS(DLDT)*V(RO)/(AREA*AREA*CEFF*CEFF)
C
      FVIS=.5*DLDT*FVISDL
C     =====================================
C     =      CALCOLA LA FORZA TOTALE      =
C     =====================================
      FTOT=-(FELAS+FVIS+FATT)
      FELDL=-FELDL
      FVISDL=-FVISDL
      RETURN
      END


      FUNCTION RINT(P,X,NX)
      DOUBLE PRECISION RINT,P,X(*)
      INTEGER NX,I

      IF(P.LE.X(1)) THEN
        RINT=X(NX+1)
        RETURN
      END IF
      IF(P.GT.X(NX)) THEN
        RINT=X(NX+NX)
        RETURN
      END IF

      I=2
 2    IF(P.LT.X(I)) THEN
        RINT=((X(NX+I)-X(NX+I-1))*P+
     &    X(I)*X(NX+I-1)-X(I-1)*X(NX+I))/(X(I)-X(I-1))
        RETURN
      END IF
      I=I+1
      GOTO 2
      END
