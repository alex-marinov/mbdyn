C $Header$
C MBDyn (C) is a multibody analysis code. 
C http://www.mbdyn.org
C 
C Copyright (C) 1996-2017
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
C ============================================================================
C
C These routines have been supplied by
C Professor Massimiliano Lanz <massimiliano.lanz@polimi.it>
C
C ============================================================================
C 
C*************************          AEROD          *********************
c                        MODIFICATA NEL CALCOLO TNG
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE AEROD2(W, VAM, TNG, OUTA, INST, RSPEED, JPRO)
C
C   W: vettore lungo 6; velocita' del corpo aerodinamico nel sistema locale,
C      riferita ad un punto arbitrario
C   VAM:    vettore lungo 6, contiene dati del problema
C   TNG: vettore lungo 12, presumubilmente il termine noto
C        (ora l'ho ridotto a 6)
C   OUTA: vettore lungo 20, usato in i/o con subroutines di aerod
C   INST: flag di instazionario (>0)
C   RSPEED -> MODOMEGA: modulo della velocita' di rotazione
C   JPRO: tipo di profilo      


      IMPLICIT NONE

C Input/Output
      REAL*8 W(6), VAM(6), TNG(6), OUTA(*), RSPEED
      INTEGER*4 INST, JPRO

C Local
      REAL*8 DENS, CS, CORDA, B, D, VP2, VP, V, RK
      INTEGER*4 IGO, I

      REAL*8 CLIFT, CDRAG, CMOME, ASLOP, BSLOP,
     &  CSLOP, CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM
      
      REAL*8 DM(3,6), VCSTR(6), BM(6,4), AASTR(6,6),
     &  A1STR(4,3), AP(4,6), QASTR(6)
C
C     DEFINIZIONI VETTORE VAM
C
C VAM(1): densita' dell'aria
C VAM(2): celerita' del suono
C VAM(3): corda
C VAM(4): 1/4 corda
C VAM(5): 3/4 corda
C VAM(6): svergolamento (non usato; aggiunto esternamente alla rotaz. profilo)
C
      DENS  = VAM(1)
      CS    = VAM(2)
      CORDA = VAM(3)
      B     = VAM(4)
      D     = VAM(5)
C     SVER  = VAM(6)
C
C     COSTRUZIONE VC*
C
      CALL DZERO(DM, 18)
      DO 10 I = 1,3
        DM(I,I) = 1.D0
 10   CONTINUE
      DM(2,6) = D
CCC Ci vuole anche questo ?!?
      DM(3,5) = -D
C
C Ora DM ha la struttura seguente:
C   [ 1  0  0   0  0  0 ]
C   [ 0  1  0   0  0  D ]
C   [ 0  0  1   0 -D  0 ]
C
C e quindi se la velocita' W del corpo aerodinamico e' organizzata come
C W = [v1 v2 v3 w1 w2 w3 ]
C con:
C 1 direzione di avanzamento,
C 2 direzione normale,
C 3 direzione lungo l'apertura,
C DM*W da' la velocita' nel punto a 3/4 della corda      
C
C Moltiplica DM per W a dare VCSTR (velocita' nel sistema locale, punto 3/4 c)
      CALL DPROMV(DM, 3, 3, 6, W, VCSTR, 0)
C
C VP2 e' il modulo della velocita' nel piano del profilo al quadrato
C VP e' il modulo della velocita' nel piano del profilo
C V e' il modulo della velocita'
      VP2 = VCSTR(1)**2+VCSTR(2)**2
      VP = DSQRT(VP2)
      V = DSQRT(VP2+VCSTR(3)**2)
C
C Se il numero di Mach e' troppo piccolo, si ferma
      IF(V/CS.LT.1.D-6) THEN
        CALL DZERO(TNG, 6)
        RETURN
      ENDIF
C
C VP: modulo della velocita' nel piano del profilo
C V: modulo della velocita' totale
C
C     CALCOLO COEFFICIENTI AERODINAMICI
C
C INST e' il flag di instazionarieta';
C puo' valere 0, 1, 2 e forza l'utilizzo del metodo PK del relativo ordine
C FIXME: solo 0 e' valido, gli altri due sono probabilmente bacati ...
      IGO = INST+1
      GOTO(100, 200, 300), IGO
C
C Steady
 100  CALL COE0(VCSTR, OUTA, CS,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM, JPRO)
      GOTO 400
C
C Unsteady, HARRIS A.H.S. JULY 1970
 200  CALL COE1(VCSTR, OUTA, CS, CORDA, RSPEED, 
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM, JPRO)
      GOTO 400
C
C Unsteady, BIELAWA 31TH A.H.S. FORUM 1975
 300  CALL COE2(VCSTR, OUTA, CS, CORDA, RSPEED,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM, JPRO)

C
C     COSTRUZIONE A1* ?
C
 400  CONTINUE
      CALL DZERO(A1STR, 12)
      RK = -.5D0*DENS*CORDA
      A1STR(1,1) =  RK*(CRF*(VP-V)-CDRAG*VP)
      A1STR(1,2) = -RK*CLIFT*VP
      A1STR(2,1) = -A1STR(1, 2)
      A1STR(2,2) =  A1STR(1, 1)
      A1STR(3,3) = -RK*CRF*V
      A1STR(4,1) = -RK*CMOME*CORDA*VCSTR(1)
      A1STR(4,2) = -RK*CMOME*CORDA*VCSTR(2)
C
C La matrice BM ha la forma
C
C   [ 1  0  0    0 ]
C   [ 0  1  0    0 ]
C   [ 0  0  1    0 ]
C
C   [ 0  0  0    0 ]
C   [ 0  0  0    0 ]
C   [ 0  B  0    1 ]
C
C dove B e' la posizione del punto a 1/4 della corda rispetto al riferimento
C (punto di applicazione delle forze aerodinamiche)
C
      CALL DZERO(BM, 24)                                      
      DO 20 I = 1,3
        BM(I,I) = 1.D0
 20   CONTINUE
      BM(6,2) = B
      BM(6,4) = 1.D0
      CALL DPROMM(A1STR, 4, 4, 3, DM, 3, AP, 4, 6, 0, 0)
      CALL DPROMM(BM, 6, 6, 4, AP, 4, AASTR, 6, 6, 0, 0)
C
C     COSTRUZIONE Q* ?
C
CCCCC Copia gli ultimi 3 termini di W in VCSTR (elimina l'effetto della
CCCCC moltiplicazione iniziale per DM)
      CALL MOVE(VCSTR(4), W(4), 3)
CCCCC Quindi moltiplica la matrice AASTR per VCSTR a dare QASTR, le forze
      CALL DPROMV(AASTR, 6, 6, 6, VCSTR, QASTR, 0)
C     cambio segno alle forze perchè le calcola con segno opposto
      do I = 1,6
        TNG(I) = -QASTR(I)
      end do  
C
      END
      

      
C*************************          COE0          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE COE0(VCSTR, OUTA, CS,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DCMDM, DCRFDM, JPRO)


C Modificato da Pierangelo Masarati 19/11/97
      implicit none
C
C Input      
      real*8 VCSTR(*), CS
      integer*4 JPRO
C
C Output
      real*8 OUTA(*),
     &  CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,CRF,CFSLOP,DCPDM,DCRDRM,
     &  DCMDM,DCRFDM
C
C Local
      real*8 ALFA,VC1,GAM,COSGAM,RMACH,ASLRF,ASLOP0,CSLOP0
C
C Data      
      real*8 PG, DEGRAD
      DATA PG /3.14159265358979323846D0/, DEGRAD /.017453293D0/
C
C     CALCOLA I COEFFICIENTI AERODINAMICI CON TEORIA STAZIONARIA
C     (CON CORREZIONE PER FRECCIA SECONDO HARRIS A.H.S.JULY 1970)
C
      ALFA = DATAN2(-VCSTR(2),VCSTR(1))
      OUTA(2) = ALFA/DEGRAD
      VC1 = DABS(VCSTR(1))
      GAM = DATAN2(-VCSTR(3),VC1)
      OUTA(3) = GAM/DEGRAD
      IF(DABS(GAM).GT.PG/3.D0) GAM = PG/3.D0
      COSGAM = DCOS(GAM)
      RMACH = DSQRT(VCSTR(1)**2+VCSTR(2)**2+VCSTR(3)**2)/CS
      RMACH = RMACH*DSQRT(COSGAM)
      OUTA(4) = RMACH
      CALL CPCRCM(ALFA, RMACH,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DCMDM, DCRFDM,
     &  ASLOP0, CSLOP0, JPRO, 0)
      ASLRF = ASLOP0
      IF(DABS(ALFA).LT.1.D-6) GOTO 10
      ASLRF = CLIFT/(ALFA*COSGAM)
      IF(ASLRF.GT.ASLOP0) ASLRF = ASLOP0
   10 CLIFT = ASLRF*ALFA
      OUTA(5) = CLIFT
      OUTA(6) = CDRAG
      OUTA(7) = CMOME
C
      RETURN
C
C
      END
C*************************          COEPRD          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE COEPRD(DA, OUTA)
      IMPLICIT NONE
      REAL*8 DA, OUTA(*)
      REAL*8 ALF1, ALF2
      ALF1 = OUTA(9)/DA
      OUTA(9) = ALF1
      ALF2 = OUTA(10)/(DA*DA)
      OUTA(10) = ALF2
      RETURN
      END
C*************************          COE1          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE COE1(VCSTR, OUTA, CS, CORDA, RSPEED,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DCMDM, DCRFDM, JPRO)

      
      implicit none     
C      IMPLICIT REAL*8(A-H,O-Z)
C
C Input:
      real*8 VCSTR(*),CS,CORDA,RSPEED
      integer*4 JPRO
C
C Output:
      real*8 OUTA(*),
     &  CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,CRF,CFSLOP,
     &  DCPDM,DCRDRM,DCMDM,DCRFDM
C
C Functions:
      real*8 THF,THG
C
C Local:
      real*8 ABE,VC1,ALF1,ALF2,VP,GAM,COSGAM,RMACH,FR,RMM,DALF,
     &  SEGNO,AREF,ASLOP0,CSLOP0,ASLRF,RK,PA,ALF
C
C Data:
      real*8 PG,DEGRAD
      DATA PG /3.14159265358979323846D0/, DEGRAD /.017453293D0/
C
C     CALCOLA COEFFICIENTI AERODINAMICI CON TEORIA INSTAZIONARIA
C     SECONDO TEORIA HARRIS A.H.S. JULY 1970
C
      ABE = DATAN2(-VCSTR(2),VCSTR(1))
      OUTA(2) = ABE/DEGRAD
      VC1 = DABS(VCSTR(1))
C
C Questa operazione e' stata spostata all'esterno (COEPRD) per consentire
C l'uso di integratori impliciti, per i quali la correzione deve essere fatta
C solo a processo iterativo concluso
C 
C      DAA = DA*RSPEED
C      ALF1 = OUTA(9)/DAA
C      OUTA(9) = ALF1
C      ALF2 = OUTA(10)/(DAA*DAA)
C      OUTA(10) = ALF2
C
C Sono: ALF1 = theta/t
C       ALF2 = theta/tt
      ALF1 = OUTA(9)     
      ALF2 = OUTA(10)
C
      VP = DSQRT(VCSTR(1)**2+VCSTR(2)**2)
      GAM = DATAN2(-VCSTR(3),VC1)
      OUTA(3) = GAM/DEGRAD
      IF(DABS(GAM).GT.PG/3.D0) GAM = PG/3.D0
      COSGAM = DCOS(GAM)
      RMACH = DSQRT(VCSTR(1)**2+VCSTR(2)**2+VCSTR(3)**2)/CS
      RMACH = RMACH*DSQRT(COSGAM)
      OUTA(4) = RMACH
C Viene usato ALF1 per la prima volta: correzione instazionaria
      FR = .5D0*CORDA*ALF1*RSPEED/VP
      FR = DABS(FR)
      OUTA(11) = FR
C Mach effettivo usato solo tra .3 e .6
      RMM = .3D0
      IF(RMACH.GT..3D0) RMM = RMACH
      IF(RMACH.GT..6D0) RMM = .6D0
C
      DALF = 61.5D0*DLOG(.6D0/RMM)*DSQRT(FR)*PG/180.D0
      SEGNO = 1.D0
      IF(ALF1.LT.0.D0) SEGNO = -1.D0
      AREF = ABE-SEGNO*DALF
      OUTA(12) = AREF/DEGRAD
      CALL CPCRCM(AREF, RMACH,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DCMDM, DCRFDM,
     &  ASLOP0, CSLOP0, JPRO, 0)
      ASLRF = ASLOP0
C Se l'angolo di incidenza e' significativo corregge la slope di CL/alfa
      IF(DABS(AREF).LT.1.D-6) GOTO 10
      ASLRF = CLIFT/(AREF*COSGAM)
      IF(ASLRF.GT.ASLOP0) ASLRF = ASLOP0
 10   CONTINUE
C Frequenza ridotta (spero che con VP = 0 non si entri qui ...)
      RK = .5D0*CORDA*RSPEED/VP
      PA = .25D0
C Corregge l'angolo di incidenza attraverso le funzioni di Theodorsen
      ALF = THF(RK)*ABE+(.5D0*RK+THG(RK))*ALF1
      ALF = ALF+2.D0*(.75D0-PA)*THF(RK)*RK*ALF1
C Viene usato ALF2 per la prima volta
      ALF = ALF-RK*RK*(PA-.5D0)*ALF2
C Coefficiente di portanza corretto per Theodorsen (?)
      CLIFT = ASLRF*ALF
      OUTA(5) = CLIFT
      OUTA(6) = CDRAG
      OUTA(7) = CMOME
      OUTA(13) = ALF/DEGRAD
C
      RETURN
C
C
      END
C*************************          COE2          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE COE2(VCSTR, OUTA, CS, CORDA, RSPEED,
     &  CLIFT, CDRAG, CMOME, ASLOP, BSLOP, CSLOP,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DCMDM, DCRFDM, JPRO)

      implicit none
      integer*4 i
C      IMPLICIT REAL*8(A-H,O-Z)
C
C Input:
C
C VCSTR(*)	Speed
C OUTA(*)	Array with different parameters, used as I/O
C CS		?
C CORDA		Chord
C RSPEED	?
C
      real*8 VCSTR(*),CS,CORDA,RSPEED
      integer*4 JPRO
C
C Output:
C
C OUTA(*)	Array with different parameters, used as I/O
C CLIFT		?
C CDRAG		?
C CMOME		?
C ASLOP		?
C BSLOP		?
C CSLOP		?
C CRF		?
C CFSLOP	?
C DCPDM		?
C DCRDRM	?
C DCMDM		?
C DCRFDM	?
C DUM		?
C
      real*8 OUTA(*),CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,CRF,CFSLOP,
     &  DCPDM,DCRDRM,DCMDM,DCRFDM,DUM
C
C Local:
      real*8 ALFA,ALF1,ALF2,VP,VP2,A,B,VC1,GAM,COSGAM,RMACH,ETA,SEGNO,
     &  SEGNOE,
     &  A2,B2,ASN,ASM,SGN,SGM,SGMAX,S2,DAN,DCN,DAM,DCM,
     &  ASLOP0,CSLOP0,ASLRF,C1,
     &  ATMP
C
C Data:
      real*8 PG,DEGRAD,ASN0,ASM0,PN(14),QN(14),PM(14),QM(14)
      DATA PG /3.14159265358979323846D0/, DEGRAD /.017453293D0/
      DATA ASN0 /.22689D0/, ASM0 /.22689D0/
      DATA PN /-3.464003D-1,-1.549076D+0, 4.306330D+1,-5.397529D+1,
     &          5.781402D+0,-3.233003D+1,-2.162257D+1, 1.866347D+1,
     &          4.198390D+1, 3.295461D+2, 4*0.D0/
      DATA QN / 1.533717D+0, 6.977203D+0, 1.749010D+3, 1.694829D 3,
     &         -1.771899D+3,-3.291665D+4, 2.969051D+0,-3.632448D+1,
     &         -2.268578D+3, 6.601995D+3,-9.654208D+3, 8.533930D+4,
     &         -1.492624D+0, 1.163661D+1/
      DATA PM / 1.970065D+1,-6.751639D+1, 7.265269D+2, 4.865945D+4,
     &          2.086279D+4, 6.024672D+3, 1.446334D+2, 8.586896D+2,
     &         -7.550329D+2,-1.021613D+1, 2.247664D+1, 3*0.D0/
      DATA QM /-2.322808D+0,-1.322257D+0,-2.633891D+0,-2.180321D-1,
     &          4.580014D+0, 3.125497D-1,-2.828806D+1,-4.396734D+0,
     &          2.565870D+2,-1.204976D+1,-1.157802D+2, 8.612138D+0,
     &          2*0.D0/
C
C     CALCOLA I COEFFICIENTI AERODINAMICI CON TEORIA INSTAZIONARIA
C     CON DATI SINTETIZZATI DI BIELAWA 31TH A.H.S. FORUM 1975
C
      ALFA = DATAN2(-VCSTR(2),VCSTR(1))
      OUTA(2) = ALFA/DEGRAD
C
C Questa operazione e' stata spostata all'esterno (COEPRD) per consentire
C l'uso di integratori impliciti, per i quali la correzione deve essere fatta
C solo a processo iterativo concluso
C
C Coefficienti usati per il ritardo instazionario
      ALF1 = OUTA(9)
      ALF2 = OUTA(10)
C
      VP2 = VCSTR(1)**2+VCSTR(2)**2
      VP = DSQRT(VP2)

C	print *,'CORDA=',CORDA,', ALF1=',ALF1,', ALF2=',ALF2,
C     &	', VP=',VP,', VP2=',VP2
      
      A = .5D0*CORDA*ALF1/VP
      B = .25D0*CORDA*CORDA*ALF2/VP2
      VC1 = DABS(VCSTR(1))
      GAM = DATAN2(-VCSTR(3),VC1)
      OUTA(3) = GAM/DEGRAD
      IF(DABS(GAM).GT.PG/3.D0) GAM = PG/3.D0
      COSGAM = DCOS(GAM)
      RMACH = DSQRT((VP2+VCSTR(3)**2)*COSGAM)/CS
      OUTA(4) = RMACH
      IF(RMACH.GT..99D0) RMACH = .99D0
      ETA = DSQRT((A/.048D0)**2+(B/.016D0)**2)
      SEGNO = 1.D0
      IF(ALFA.LT.0D0) SEGNO = -1.D0
      SEGNOE = SEGNO
      IF(ETA.GT.1.D0) SEGNOE = SEGNOE/ETA
      A = SEGNOE*A
      B = SEGNOE*B
      A2 = A*A
      B2 = B*B
      ASN = ASN0*(1.D0-RMACH)
      ASM = ASM0*(1.D0-RMACH)
      SGN = DABS(ALFA/ASN)
      SGM = DABS(ALFA/ASM)
      SGMAX = 1.839D0-70.33D0*B
      IF(SGMAX.GT.1.86D0) SGMAX = 1.86D0
      IF(SGN.GT.SGMAX) SGN = SGMAX
      IF(SGM.GT.SGMAX) SGM = SGMAX

C	print *,'A=',A,', B=',B,', SGN=',SGN,', PN(1:6)=',(PN(i),i=1,6)
      
      DAN = A*(PN(1)+PN(5)*SGN)+B*(PN(2)+PN(6)*SGN)

C	print *,'1) DAN=',DAN

      DAN = DAN+DEXP(-1072.52D0*A2)*(A*(PN(3)+PN(7)*SGN)+
     &  A2*(PN(9)+PN(10)*SGN))

C	print *,'2) DAN=',DAN

      DAN = DAN+DEXP(-40316.42D0*B2)*B*(PN(4)+PN(8)*SGN)

C	print *,'3) DAN=',DAN

      DAN = DAN*ASN

C	print *,'4) DAN=',DAN

      DCN = A*(QN(1)+QN(3)*A2+SGN*(QN(7)+QN(9)*A2+QN(13)*SGN)+
     &  B2*(QN(5)+QN(11)*SGN))
      DCN = DCN+B*(QN(2)+QN(4)*A2+SGN*(QN(8)+QN(10)*A2+QN(14)*SGN)+
     &  B2*(QN(6)+QN(12)*SGN))
      DAM = A*(PM(1)+PM(3)*A2+PM(5)*B2+PM(10)*SGM+PM(7)*A)
      DAM = DAM+B*(PM(2)+PM(4)*B2+PM(6)*A2+PM(11)*SGM+PM(8)*B+PM(9)*A)
      DAM = DAM*ASM
      S2 = SGM*SGM
      DCM = A*(QM(2)+QM(8)*A+SGM*(QM(4)+QM(10)*A)+S2*(QM(6)+QM(12)*A))
      DCM = DCM+B*(QM(1)+QM(7)*B+SGM*(QM(3)+QM(9)*B)+
     &  S2*(QM(5)+QM(11)*B))
      OUTA(11) = DAN/DEGRAD
      OUTA(12) = DAM/DEGRAD
      OUTA(13) = DCN
      OUTA(14) = DCM

      DAN = SEGNO*DAN
      DCN = SEGNO*DCN
      DAM = SEGNO*DAM
      DCM = SEGNO*DCM

C
C Questa parte e' sbagliata, perche' per calcolare i coeff. di momento
C e le loro derivate vengono usati i coeff di portanza e le loro derivate
C ecc, quindi passando DUM come segnaposto, si ottengono risultati
C unpredictable.
C
C	print *,'ATMP = ALFA(=',ALFA,')-DAN(',DAN,')'

      ATMP = ALFA-DAN
      CALL CPCRCM(ATMP, RMACH,
     &  CLIFT, CDRAG, DUM, ASLOP, BSLOP, DUM,
     &  CRF, CFSLOP, DCPDM, DCRDRM, DUM, DCRFDM,
     &  ASLOP0, CSLOP0, JPRO, 0)

      ATMP = ALFA-DAM
      CALL CPCRCM(ATMP, RMACH,
     &  DUM, DUM, CMOME, DUM, DUM, CSLOP,
     &  DUM, DUM, DUM, DUM, DCMDM, DUM,
     &  ASLOP0, CSLOP0, JPRO, 1)
      
      CMOME = CMOME-.25D0*CLIFT
      CSLOP = CSLOP-.25D0*ASLOP
      DCMDM = DCMDM-.25D0*DCPDM      
      
      ASLRF = ASLOP0
      IF(DABS(ALFA-DAN).LT.1.D-6) GOTO 10
      ASLRF = CLIFT/((ALFA-DAN)*COSGAM)
      IF(ASLRF.GT.ASLOP0) ASLRF = ASLOP0
 10   CLIFT = ASLRF*(ALFA-DAN)
      C1 = .9457D0/DSQRT(1.D0-RMACH*RMACH)

C	print *,'COE2 ',CLIFT,ASLOP0*DAN,DCN*C1,CMOME,CSLOP0*DAM,DCM*C1
      
C La parte steady di CLIFT e' 0 quando alfa cambia segno!!!
C      CLIFT = CLIFT+SEGNO*(ASLOP0*DAN+DCN*C1)
C      CMOME = SEGNO*(CMOME+CSLOP0*DAM+DCM*C1)
      CLIFT = CLIFT+ASLOP0*DAN+DCN*C1
      CMOME = CMOME+CSLOP0*DAM+DCM*C1
      OUTA(5) = CLIFT
      OUTA(6) = CDRAG
      OUTA(7) = CMOME
C
      RETURN
C
C
      END


      
C*************************          CPCRCM          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE CPCRCM (APHIJ, EMIJ, CLIFT, CDRAG, CMOME, ASLOP,
     &  BSLOP,CSLOP, CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM,
     &  AS0, CS0, JPRO, IRETRN)
C
C APHIJ		ALFA
C EMIJ		RMACH
C CLIFT		CLIFT
C CDRAG		CDRAG
C CMOME		CMOME
C ASLOP		ASLOP
C BSLOP		BSLOP
C CSLOP		CSLOP
C CRF		CRF
C CFSLOP	CFSLOP
C DCPDM		DCPDM
C DCDRDM	DCDRDM
C DCMDM		DCMDM
C DCRFDM	DCRFDM
C AS0		ASLOP0
C CS0		CSLOP0
C JPRO		JPRO
C IRETRN	IRETRN
C
      IMPLICIT NONE
C
      REAL*8 APHIJ, EMIJ, CLIFT, CDRAG, CMOME, ASLOP,
     &  BSLOP,CSLOP, CRF, CFSLOP, DCPDM, DCDRDM, DCMDM, DCRFDM,
     &  AS0, CS0
      INTEGER*4 JPRO, IRETRN
C
      REAL*8 D1, D2, D3, CD0, B1, B2, B3,
     &  PG, SQT, C1, C2, C3, C4, C5, S, S2, S3, S4,
     &  P1, P2, P3, R0, R1, R2, R3, PEND0, PEND1, PEND2, PEND3,
     &  APHIJ2, APHIJ3
C
      INTEGER*4 NEG
C
C
      DIMENSION D1(2,5), D2(2,5), D3(2,5), CD0(2,5),
     &  B1(2,5), B2(2,5), B3(2,5)
      DATA D1 /.3D0,-.28532D2,.4D0,-.37392D2,.45D0,-.4597D2,.5D0,
     &         -.55819D2,.55D0,-.88102D2/
      DATA D2 /.3D0,.4826D1,.4D0,.70982D1,.45D0,.8914D1,.5D0,
     &         .10636D2,.55D0,.1993D2/
      DATA D3 /.3D0,.62085D1,.4D0,.64158D1,.45D0,.65334D1,.5D0,
     &         .66613D1,.55D0,.64065D1/
      DATA CD0 /.3D0,.60145D-2,.4D0,.35772D-2,.45D0,.2612D-2,.5D0,
     &         .46674D-2,.55D0,.29977D-2/
      DATA B1 /.3D0,.86845D-1,.4D0,.1915D0,.45D0,.24923D0,.5D0,
     &         .14568D0,.55D0,.25389D0/
      DATA B2 /.3D0,-.91795D0,.4D0,-.1966D1,.45D0,-.26795D1,.5D0,
     &         -.15943D1,.55D0,-.30426D1/
      DATA B3 /.3D0,.40357D1,.4D0,.67937D1,.45D0,.91332D1,.5D0,
     &         .62569D1,.55D0,.12049D2/
      DATA PG /3.14159265358979323846D0/
C**** APHIJ  = ANGOLO INCIDENZA
C**** EMIJ   = NUMERO DI MACH
C**** CLIFT  = COEFFICIENTE DI PORTANZA
C**** CDRAG  = COEFFICIENTE DI RESISTENZA TOTALE
C**** CMOME  = COEFFICIENTE DI MOMENTO (CONVENZIONE CABRANTE) RIFERITO
C****          AL CENTRO AERODINAMICO
C**** ASLOP  = PENDENZA CURVA DI PORTANZA
C**** BSLOP  = PENDENZA CURVA DI RESISTENZA TOTALE
C**** CSLOP  = PENDENZA CURVA DI MOMENTO
C**** CRF    = COEFFICIENTE DI RESISTENZA ATTRITO
C**** CFSLOP = PENDENZA CURVA DI RESISTENZA ATTRITO
C**** DCPDM  = DERIVATA DI CLIFT RISPETTO NUMERO DI MACH
C**** DCDRDM = DERIVATA DI CDRAG RISPETTO NUMERO DI MACH
C**** DCMDM  = DERIVATA DI CMOME RISPETTO NUMERO DI MACH
C**** DCRFDM = DERIVATA DI CRF RISPETTO NUMERO DI MACH
C**** AS0    = PENDENZA CURVA PORTANZA PER ALFA=0
C**** CS0    = PENDENZA CURVA MOMENTO PER ALFA=0
      CFSLOP = 0.D0
      CRF = .006D0
      DCRFDM = 0.D0
      NEG = 1
      IF(EMIJ.GT..99D0) EMIJ = .99D0
      SQT = DSQRT(1.D0-EMIJ*EMIJ)
      C1 = 1.D0-EMIJ
      C2 = .22689D0*C1
      C5 = EMIJ/(SQT*SQT)

C	print *,'APHIJ: ',APHIJ
      
97    IF(APHIJ.lt.0) then
         goto 181
      else 
         goto 182
      end if
181   APHIJ = -APHIJ
      NEG = -1*NEG
182   IF(APHIJ-PG.le.0) then
         goto 184
      else
          goto 183
      end if
183   APHIJ = APHIJ-PG*2.D0
      GOTO 97
184   IF(APHIJ-C2.lt.0) then
          goto 185
      else
          goto 187
      end if
185   ASLOP = 5.7296D0/SQT
      CLIFT = ASLOP*APHIJ
      CDRAG = .006D0+.13131D0*APHIJ*APHIJ
      BSLOP = .26262D0*APHIJ
      CMOME = 1.4324D0*APHIJ/SQT
      CSLOP = 1.4324D0/SQT
      DCPDM = CLIFT*C5
      DCDRDM = 0.D0
      DCMDM = CMOME*C5
      GOTO 250
187   IF(APHIJ-0.34906D0.lt.0) then
         goto 189
      else
         goto 191
      end if
189   CLIFT = .29269D0*C1+(1.3D0*EMIJ-.59D0)*APHIJ
      C2 = (.12217D0+.22689D0*EMIJ)*SQT
      CMOME = CLIFT/(4*C2)
      CSLOP = (1.3D0*EMIJ-.59D9)/(4.D0*C2)
      CLIFT = CLIFT/C2
      ASLOP = (1.3D0*EMIJ-.59D0)/C2
      DCPDM = (-.29269D0+1.3D0*APHIJ
     &  -CLIFT*(-(.12217D0+.22689D0*EMIJ)*EMIJ/SQT+
     &  .22689D0*SQT))/C2
      DCMDM = .25D0*DCPDM
      GOTO 210
  191 IF(APHIJ-2.7402D0.lt.0) then
          goto 193
      else
          goto 195
      end if
  193 S = DSIN(APHIJ)
      S2 = DSIN(2.*APHIJ)
      S3 = DSIN(3.*APHIJ)
      S4 = DSIN(4.*APHIJ)
      CLIFT = (.080373D0*S+1.04308D0*S2
     &  -.011059D0*S3+.023127D0*S4)/SQT
      CMOME = (-.02827D0*S+.14022D0*S2
     &  -.00622D0*S3+.01012D0*S4)/SQT
      C1 = DCOS(APHIJ)
      C2 = DCOS(2.D0*APHIJ)
      C3 = DCOS(3.D0*APHIJ)
      C4 = DCOS(4.D0*APHIJ)
CCCCC CLIFT = CLIFT/C2
      CSLOP = (-.02827D0*C1+.28044D0*C2
     &  -.01866D0*C3+.04048D0*C4)/SQT
      ASLOP = (.080373D0*C1+2.08616D0*C2
     &  -.033177D0*C3+.092508D0*C4)/SQT
      DCPDM = CLIFT*C5
      DCMDM = CMOME*C5
      GOTO 210
195   IF(APHIJ-3.0020D0.lt.0) then 
          goto 197
      else 
          goto 199
      endif
197   CLIFT = -(.4704D0+.10313D0*APHIJ)/SQT
      ASLOP = -.10313D0/SQT
      CMOME = -(.4786D0+.02578D0*APHIJ)/SQT
      CSLOP = .02578D0/SQT
      DCPDM = CLIFT*C5
      DCMDM = CMOME*C5
      GOTO 210
199   IF(APHIJ-PG.le.0) then
          goto 200
      else
          goto 260
      end if
200   CLIFT = (-17.55D0+5.5864D0*APHIJ)/SQT
      ASLOP = 5.5864D0/SQT
      CMOME = (-12.5109D0+3.9824D0*APHIJ)/SQT
      CSLOP = 3.9824D0/SQT
      DCPDM = CLIFT*C5
      DCMDM = CMOME*C5
210   CDRAG = 1.1233D0-.029894D0*DCOS(APHIJ)
     &  -1.00603D0*DCOS(2.D0*APHIJ)
      CDRAG = CDRAG+.003115D0*DCOS(3.D0*APHIJ)
     &  -.091487D0*DCOS(4.D0*APHIJ)
      CDRAG = CDRAG/SQT
      BSLOP = .029894D0*DSIN(APHIJ)+2.01206D0*DSIN(2.D0*APHIJ)
      BSLOP = (BSLOP+.009345D0*DSIN(3.D0*APHIJ)
     &  +.365948D0*DSIN(4.D0*APHIJ))/SQT
      DCDRDM = CDRAG*C5
      
 250  GOTO(251, 252), JPRO
C
C     CALCOLO COEFFICIENTI AERODINAMICI PER PROFILO 'RAE9671'
C
 252  IF(DABS(APHIJ).GT..21956242D0) GOTO 251
      IF(NEG.LT.0) APHIJ = -APHIJ
      CALL LININT(D1, 5, EMIJ, P1, PEND1)
      CALL LININT(D2, 5, EMIJ, P2, PEND2)
      CALL LININT(D3, 5, EMIJ, P3, PEND3)
      APHIJ2 = APHIJ*APHIJ
      APHIJ3 = APHIJ*APHIJ2
      CLIFT = P1*APHIJ3+P2*APHIJ2+P3*APHIJ
      ASLOP = 3.D0*P1*APHIJ2+2.D0*P2*APHIJ+P3
      AS0 = P3
      DCPDM = PEND1*APHIJ3+PEND2*APHIJ2+PEND3*APHIJ
C
      CALL LININT(CD0, 5, EMIJ, R0, PEND0)
      CRF = R0
      CFSLOP = 0.D0
      DCRFDM = PEND0
      CALL LININT(B1, 5, EMIJ, R1, PEND1)
      CALL LININT(B2, 5, EMIJ, R2, PEND2)
      CALL LININT(B3, 5, EMIJ, R3, PEND3)
      APHIJ = DABS(APHIJ)
      APHIJ3 = DABS(APHIJ3)
      CDRAG = R0+APHIJ*R1+APHIJ2*R2+APHIJ3*R3
      BSLOP = R1+2.D0*R2*APHIJ+3.D0*R3*APHIJ2
      DCDRDM = PEND0+APHIJ*PEND1+APHIJ2*PEND2+APHIJ3*PEND3
C
      IF(NEG.GT.0) GOTO 249
      CMOME = -CMOME
      CSLOP = -CSLOP
      DCMDM = -DCMDM
 249  CONTINUE
      CS0 = 0.D0
C
      GOTO 270
C
C      CMOME = CMOME-.25D0*CLIFT
C      CSLOP = CSLOP-.25D0*ASLOP
C      DCMDM = DCMDM-.25D0*DCPDM
C      RETURN
C
C
C     CALCOLO COEFFICIENTI AERODINAMICI PER PROFILO 'NACA0012'
C
 251  CONTINUE

      IF(NEG.le.0) then
          goto 255
      else
          goto 260
      end if
255   CLIFT = -CLIFT
      CMOME = -CMOME
      APHIJ = -APHIJ
      DCPDM = -DCPDM
      DCMDM = -DCMDM     
260   CONTINUE
C
      AS0 = 5.7296D0/SQT
C
C Qui arriva anche l'altro profilo per una parte di operazioni in comune
C
270   CONTINUE
      CS0 = 0.D0
      IF(IRETRN.EQ.1) RETURN
C
      CMOME = CMOME-.25D0*CLIFT
      CSLOP = CSLOP-.25D0*ASLOP
      DCMDM = DCMDM-.25D0*DCPDM
C
      RETURN
      END
C*************************          PSIROT          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE PSIROT(PSI, A, NRDA, NC, B, NRDB)
      IMPLICIT NONE

      INTEGER*4 NRDA, NC, NRDB
      REAL*8 PSI, A(NRDA,1), B(NRDB,1)
      
      REAL*8 SINP, COSP, TMP1, TMP2
      INTEGER*4 K

      
      SINP = DSIN(PSI)
      COSP = DCOS(PSI)
      DO 20 K = 1,NC
      TMP1 = +COSP*A(1,K)-SINP*A(2,K)
      TMP2 = +SINP*A(1,K)+COSP*A(2,K)
      B(1,K) = TMP1
      B(2,K) = TMP2
      B(3,K) = A(3,K)
  20  CONTINUE
  30  CONTINUE
      RETURN
      END
C*************************          THF F          *********************
C=    COMPILER (LINK=IBJ$)
      FUNCTION THF(K)
      IMPLICIT NONE
      REAL*8 THF, K, A, B, C, D
C ****  FUNZIONE F DI THEODORSEN ***
      A =  K**4-.102058D0*K**2+9.55732D-6
      B = -.761036D0*K**3+2.55067D-3*K
      C =  2.D0*K**4-.113928D0*K**2+9.55732D-6
      D = -1.064D0*K**3+2.6268D-3*K
      THF = (A*C+B*D)/(C*C+D*D)
      RETURN
      END
C*************************          THG F          *********************
C=    COMPILER (LINK=IBJ$)
      FUNCTION THG(K)
      IMPLICIT NONE
      REAL*8 THG, K, A, B, C, D
C ****  FUNZIONE G DI THEODORSEN ***
      A =  K**4-.102058D0*K**2+9.55732D-6
      B = -.761036D0*K**3+2.55067D-3*K
      C =  2.D0*K**4-.113928D0*K**2+9.55732D-6
      D = -1.064D0*K**3+2.6268D-3*K
      THG = (B*C-A*D)/(C*C+D*D)
      RETURN
      END
C*************************          UNST0          *********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE UNST0(ALF, PP, COE, PUNTI, NUPIA, NUPIR)
      IMPLICIT NONE

      INTEGER*4 NUPIA, NUPIR
      REAL*8 ALF(25,3), PP(5,5), COE(5,5), PUNTI(2,25),
     &  PERM(5), CD(5)

      REAL*8 DSCALL
      REAL*8 QPSI1, QPSI2,QPSI3
      INTEGER*4 I, K, L, INDER
      
      DO 10 I = 1,NUPIA
        PP(I,1) = 1.D0
  10  CONTINUE
      DO 20 K = 2,NUPIA
        DO 19 I = 1,NUPIA
          L = (I-1)*NUPIR+1
          PP(I,K) = PP(I,K-1)*PUNTI(2,L)
 19   CONTINUE
 20   CONTINUE
      CALL DRUFCT(PP, PERM, NUPIA, 5, INDER)
      DO 50 K = 1,NUPIR
        DO 30 I = 1,NUPIA
          L = (I-1)*NUPIR+K
          COE(I,1) = ALF(L,1)
  30    CONTINUE
        CALL DRUSOL(PP, COE, 5, NUPIA, PERM)
        DO 40 I = 1,NUPIA
          L = (I-1)*NUPIR+K
          QPSI1 = PUNTI(2,L)
          QPSI2 = QPSI1*QPSI1
          QPSI3 = QPSI2*QPSI1
          CD(1) = 0.D0
          CD(2) = 1.D0
          CD(3) = 2.D0*QPSI1
          CD(4) = 3.D0*QPSI2
          CD(5) = 4.D0*QPSI3
          ALF(L,2) = DSCALL(CD,COE,NUPIA)
C
          CD(1) = 0.D0
          CD(2) = 0.D0
          CD(3) = 2.D0
          CD(4) = 6.D0*QPSI1
          CD(5) = 12.D0*QPSI2
          ALF(L,3) = DSCALL(CD,COE,NUPIA)
   40   CONTINUE
   50 CONTINUE
      RETURN
      END























*********************************  STLIB  *******************************
C*************************          ADD          ***********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE ADD(A,NRDA,B,NRDB,NR,NC,PLUS)
      IMPLICIT NONE

      INTEGER*4 NRDA, NRDB, NR, NC
      REAL*8 A(NRDA,1),B(NRDB,1),PLUS

      INTEGER*4 I, K
      
      DO 10 I=1,NR
      DO 9 K=1,NC
      A(I,K)=A(I,K)+PLUS*B(I,K)
 9    end do
 10   end do
      RETURN
      END

C*************************          DABCT3          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DABCT3(A,B,NRDB,C)
C
C     CALCOLA IL PRODOTTO DELLE MATRICI (3*3) "A" * "B" * "C-TRASPOSTA"
C
      IMPLICIT NONE

      INTEGER*4 NRDB
      REAL*8 A,B,C,D
      
      DIMENSION A(3,3),B(NRDB,1),C(3,3),D(3,3)

      INTEGER*4 I, K, L
      
      DO 10 I=1,3
      DO 9 K=1,3
      D(I,K)=0.D0
      DO 8 L=1,3
      D(I,K)=D(I,K)+A(I,L)*B(L,K)
8     end do
9     end do
10    end do
      DO 20 I=1,3
      DO 19 K=1,3
      B(I,K)=0.D0
      DO 18 L=1,3
      B(I,K)=B(I,K)+D(I,L)*C(K,L)
18    end do
19    end do
20    end do
      RETURN
      END

C     ******************          DDEMA          ***********************
      SUBROUTINE DDEMA(RO,RD,DEMA,ND,FACT)
C
C     GIVEN THE ROTATION VECTOR RO AND THE DERIVATIVE OF THE ROTATION
C     VECTOR RD, CALCULATE THE CAPITAL GAMMA DOT.
C
      IMPLICIT NONE

      INTEGER*4 ND
      REAL*8 RO,RD,DEMA,FACT

      REAL*8 A,B,C,D,E,G,RDRO,ARG,ARG2,T
      INTEGER*4 I
      
      DIMENSION RO(3),RD(3),DEMA(ND,3),T(3)

      REAL*8 TRESH2
      DATA TRESH2 /1.D-10/
      
      ARG2=0.D0
      RDRO=0.D0
      DO 10 I=1,3
      T(I)=RO(I)*FACT
      ARG2=ARG2+T(I)**2
      RDRO=RDRO+T(I)*RD(I)
   10 CONTINUE
      IF (ARG2.LT.TRESH2) THEN
C
        A=(1.D0-ARG2/30.D0*(2.D0+ARG2/56.D0*(3.D0-ARG2/22.5D0)))/12.D0
        B=(1.D0-ARG2/42.D0*(2.D0-ARG2/72.D0*(3.D0-ARG2/27.5D0)))/60.D0
        D=(1.D0-ARG2/12.D0*(1.D0-ARG2/30.D0*(1.D0-ARG2/56.0D0)))/ 2.D0
        E=(1.D0-ARG2/20.D0*(1.D0-ARG2/42.D0*(1.D0-ARG2/72.0D0)))/ 6.D0
        A=-A*RDRO
        B=-B*RDRO
C
      ELSE
C
         ARG=DSQRT(ARG2)
         D=(1.D0-DCOS(ARG))/ARG2
         G=DSIN(ARG)/ARG
         E=(1.D0-G)/ARG2
         A=(G-2.D0*D)*RDRO/ARG2
         B=(D-3.D0*E)*RDRO/ARG2
C
      END IF
C
      C=-B*ARG2
      DEMA(1,1)=T(1)*T(1)*B+C
      DEMA(1,2)=T(1)*T(2)*B-T(3)*A
      DEMA(1,3)=T(1)*T(3)*B+T(2)*A
      DEMA(2,1)=T(2)*T(1)*B+T(3)*A
      DEMA(2,2)=T(2)*T(2)*B+C
      DEMA(2,3)=T(2)*T(3)*B-T(1)*A
      DEMA(3,1)=T(3)*T(1)*B-T(2)*A
      DEMA(3,2)=T(3)*T(2)*B+T(1)*A
      DEMA(3,3)=T(3)*T(3)*B+C
      DEMA(1,1)=DEMA(1,1)+ 2.D0*(RD(1)*T(1)-RDRO)*E
      DEMA(1,2)=DEMA(1,2)+(RD(1)*T(2)+T(1)*RD(2))*G-RD(3)*D
      DEMA(1,3)=DEMA(1,3)+(RD(1)*T(3)+T(1)*RD(3))*E+RD(2)*D
      DEMA(2,1)=DEMA(2,1)+(RD(2)*T(1)+T(2)*RD(1))*E+RD(3)*D
      DEMA(2,2)=DEMA(2,2)+ 2.D0*(RD(2)*T(2)-RDRO)*E
      DEMA(2,3)=DEMA(2,3)+(RD(2)*T(3)+T(2)*RD(3))*E-RD(1)*D
      DEMA(3,1)=DEMA(3,1)+(RD(3)*T(1)+T(3)*RD(1))*E-RD(2)*D
      DEMA(3,2)=DEMA(3,2)+(RD(3)*T(2)+T(3)*RD(2))*E+RD(1)*D
      DEMA(3,3)=DEMA(3,3)+ 2.D0*(RD(3)*T(3)-RDRO)*E
      RETURN
      END

C     ******************          DLGMA          ***********************
      SUBROUTINE DLGMA(RO,SPMA,ND,FACT,ID)
C
C     GIVEN THE ROTATION VECTOR RO, CALCULATE THE CAPITAL GAMMA MATRIX
C     SPMA.
C
      IMPLICIT NONE

      INTEGER*4 ND, ID
      REAL *8 RO,SPMA,FACT
      
      REAL*8 A,B,C,ARG,ARG2,T,TRESH2
      INTEGER*4 I
      
      DIMENSION RO(3),SPMA(ND,3),T(3)
      
      DATA TRESH2 /1.D-10/
      
      ARG2=0.D0
      DO 10 I=1,3
      T(I)=RO(I)*FACT
      ARG2=ARG2+T(I)**2
   10 CONTINUE
      IF (ARG2.LT.TRESH2) THEN
         A=(1.D0-ARG2/12.D0*(1-ARG2/30.D0*(1-ARG2/56.D0)))/2.D0
         B=(1.D0-ARG2/20.D0*(1-ARG2/42.D0*(1-ARG2/72.D0)))/6.D0
      ELSE
         ARG=DSQRT(ARG2)
         A=(1.D0-DCOS(ARG))/ARG2
         B=(1.D0-DSIN(ARG)/ARG)/ARG2
      END IF
      IF (ID.LT.0) THEN
          IF (ARG2.LT.TRESH2) THEN
             B=(1.D0+ARG2/60.D0*(1.D0+ARG2/420.D0))/12.D0
          ELSE
             B=(1.D0-.5D0*(1.D0-B*ARG2)/A)/ARG2
          END IF
        A=-.5D0
      END IF
         C=1.D0-B*ARG2
      SPMA(1,1)=T(1)*T(1)*B+C
      SPMA(1,2)=T(1)*T(2)*B-T(3)*A
      SPMA(1,3)=T(1)*T(3)*B+T(2)*A
      SPMA(2,1)=T(2)*T(1)*B+T(3)*A
      SPMA(2,2)=T(2)*T(2)*B+C
      SPMA(2,3)=T(2)*T(3)*B-T(1)*A
      SPMA(3,1)=T(3)*T(1)*B-T(2)*A
      SPMA(3,2)=T(3)*T(2)*B+T(1)*A
      SPMA(3,3)=T(3)*T(3)*B+C
      RETURN
      END

C*************************          DPRMMA          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DPRMMA(A,NRDA,NRA,NCA,B,NRDB,C,NRDC,NCC,IA,IB)
C
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B" E SOMMA
C     IL RISULTATO IN "C"
C     A    = PRIMA MATRICE
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"
C     NRA  = N.RO DI RIGHE UTILIZZATE PER IL PRODOTTO DI "A"
C     NCA  = N.RO DI COLONNE DI "A"
C     B    = SECONDA MATRICE
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"
C     C    = MATRICE PRODOTTO
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"
C     NCB  = N.RO DI COLONNE DI "C"
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO
C                                            1 = SI
C     IB   = INDICE DI TRASPOSIZIONE DI "B"
C
      IMPLICIT NONE

      INTEGER*4 NRDA,NRA,NCA,NRDB,NRDC,NCC,IA,IB
      REAL*8 A,B,C
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1)

      INTEGER*4 IGO,I,K,L
      
      IGO=2*IA+IB+1
      GOTO(10,20,30,40),IGO
10    DO 113 I=1,NRA
      DO 112 K=1,NCC
      DO 111 L=1,NCA
      C(I,K)=C(I,K)+A(I,L)*B(L,K)
111    CONTINUE
112    CONTINUE
113    CONTINUE
      RETURN
20    DO 120 I=1,NRA
      DO 121 K=1,NCC
      DO 122 L=1,NCA
      C(I,K)=C(I,K)+A(I,L)*B(K,L)
122   CONTINUE
121   CONTINUE
120   CONTINUE
      RETURN
30    DO 130 I=1,NCA
      DO 131 K=1,NCC
      DO 132 L=1,NRA
      C(I,K)=C(I,K)+A(L,I)*B(L,K)
132   CONTINUE
131   CONTINUE
130   CONTINUE
      RETURN
40    DO 140 I=1,NCA
      DO 141 K=1,NCC
      DO 142 L=1,NRA
      C(I,K)=C(I,K)+A(L,I)*B(K,L)
142   CONTINUE
141   CONTINUE
140   CONTINUE
      RETURN
      END

C*************************          DPROMM          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DPROMM(A,NRDA,NRA,NCA,B,NRDB,C,NRDC,NCC,IA,IB)
C
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B"
C     A    = PRIMA MATRICE
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"
C     NRA  = N.RO DI RIGHE UTILIZZATE PER IL PRODOTTO DI "A"
C     NCA  = N.RO DI COLONNE DI "A"
C     B    = SECONDA MATRICE
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"
C     C    = MATRICE PRODOTTO
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"
C     NCB  = N.RO DI COLONNE DI "C"
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO
C                                            1 = SI
C     IB   = INDICE DI TRASPOSIZIONE DI "B"
C
      IMPLICIT NONE

      INTEGER*4 NRDA,NRA,NCA,NRDB,NRDC,NCC,IA,IB
      REAL*8 A,B,C
      
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1)

      INTEGER*4 IGO,I,K,L
      
      IGO=2*IA+IB+1
      GOTO(10,20,30,40),IGO
10    DO 110 I=1,NRA
      DO 111 K=1,NCC
      C(I,K)=0.D0
      DO 112 L=1,NCA
      C(I,K)=C(I,K)+A(I,L)*B(L,K)
112   CONTINUE
111   CONTINUE
110   CONTINUE
      RETURN
20    DO 120 I=1,NRA
      DO 121 K=1,NCC
      C(I,K)=0.D0
      DO 122 L=1,NCA
      C(I,K)=C(I,K)+A(I,L)*B(K,L)
122   CONTINUE
121   CONTINUE
120   CONTINUE
      RETURN
30    DO 130 I=1,NCA
      DO 131 K=1,NCC
      C(I,K)=0.D0
      DO 132 L=1,NRA
      C(I,K)=C(I,K)+A(L,I)*B(L,K)
132   CONTINUE
131   CONTINUE
130   CONTINUE
      RETURN
40    DO 140 I=1,NCA
      DO 141 K=1,NCC
      C(I,K)=0.D0
      DO 142 L=1,NRA
      C(I,K)=C(I,K)+A(L,I)*B(K,L)
142   CONTINUE
141   CONTINUE
140   CONTINUE
      RETURN
      END


C*************************          DPROMV          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DPROMV(A,NRDA,NRA,NCA,B,C,IA)
C
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER IL VETTORE "B"
C
      IMPLICIT NONE

      INTEGER*4 NRDA,NRA,NCA,IA
      REAL*8 A,B,C

      INTEGER*4 IGO,I,K
      REAL*8 D
      
      DIMENSION A(NRDA,1),B(NCA),C(NRA)
      
      IGO=IA+1
      GOTO(100,200),IGO
100   DO 20 I=1,NRA
      D=0.D0
      DO 10 K=1,NCA
      D=D+A(I,K)*B(K)
10    end do
      C(I)=D
20    end do
      RETURN
200   DO 40 I=1,NRA
      D=0.D0
      DO 30 K=1,NCA
      D=D+A(K,I)*B(K)
30    end do
      C(I)=D
40    end do
      RETURN
      END
 
C*************************          DPRVAD          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DPRVAD(X,Y,Z)
C
C     ESEGUE IL PRODOTTO DEL VETTORE "X"
C     PER IL VETTORE "Y" E SOMMA IL RISULTATO IN "Z"
C
      IMPLICIT NONE
      
      REAL*8  X,Y,Z,W
      
      DIMENSION X(3),Y(3),Z(3),W(3)
      
      W(1)= X(2)*Y(3)-X(3)*Y(2)
      W(2)= X(3)*Y(1)-X(1)*Y(3)
      W(3)= X(1)*Y(2)-X(2)*Y(1)
      Z(1)=Z(1)+W(1)
      Z(2)=Z(2)+W(2)
      Z(3)=Z(3)+W(3)
      RETURN
      END

C*************************          DPRVVA          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DPRVVA(A,B,C,NRDC,NRA,NRB)

      IMPLICIT NONE

      INTEGER*4 NRDC,NRA,NRB
      REAL*8 A(1),B(1),C(NRDC,1)

      INTEGER*4 I,K
C
C     ESEGUE IL PRODOTTO DEL VETTORE A PER IL VETTORE B TRASPOSTO
C     LA MATRICE RISULTANTE VIENE SOMMATA IN C
C
      DO 20 I=1,NRA
      DO 21 K=1,NRB
      C(I,K)=C(I,K)+A(I)*B(K)
 21   CONTINUE
 20   CONTINUE
C
      RETURN
C
C
      END

C*************************          DRMOVE          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DRMOVE(A,B,N,H)
C
C     SOMMA AI TERMINI DEL VETTORE "A" QUELLI DI "B" MOLTIPLICATI PER H
C
      IMPLICIT NONE

      INTEGER*4 N
      REAL*8 A,B,H
      DIMENSION A(1),B(1)

      INTEGER*4 I
      
      DO 10 I=1,N
        A(I)=A(I)+B(I)*H
10    CONTINUE
      RETURN
      END

C*************************          DRUFCT          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DRUFCT(A,PERM,N,NR,INDER)      
      IMPLICIT NONE

      INTEGER*4 N,NR,INDER
      REAL*8 A,PERM
      
      DIMENSION A(NR,1),PERM(1)

      REAL*8 X,DP
      INTEGER*4 I,K,IM1,IP1,IPVT,J
      
      DO10I=1,N
      X=0.D0
      DO20K=1,N
      IF(DABS(A(I,K)).LT.X)GOTO20
      X=DABS(A(I,K))
20    CONTINUE
      IF(X.EQ.0.D0)GOTO110
      PERM(I)=1.D0/X
10    CONTINUE
      DO100 I=1,N
      IM1=I-1
      IP1=I+1
      IPVT=I
      X=0.D0
      DO50K=I,N
      DP=A(K,I)
      IF(I.EQ.1)GOTO40
      DO30J=1,IM1
      DP=DP-A(K,J)*A(J,I)
30    CONTINUE
      A(K,I)=DP
40    IF(X.GT.(DABS(DP)*PERM(K)))GOTO50
      IPVT=K
      X=DABS(DP)*PERM(K)
50    CONTINUE
      IF(X.LE.0.D0)GOTO110
      IF(IPVT.EQ.I)GOTO70
      DO60J=1,N
      X=A(IPVT,J)
      A(IPVT,J)=A(I,J)
      A(I,J)=X
60    CONTINUE
      PERM(IPVT)=PERM(I)
70    PERM(I)=IPVT
      IF(I.EQ.N)GOTO100
      X=A(I,I)
      DO90K=IP1,N
      A(K,I)=A(K,I)/X
      IF(I.EQ.1)GOTO90
      DP=A(I,K)
      DO80J=1,IM1
      DP=DP-A(I,J)*A(J,K)
80    CONTINUE
      A(I,K)=DP
90    CONTINUE
100   CONTINUE
      INDER=0
      RETURN
110   INDER=I
      RETURN
      END

C*************************          DRUSOL          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DRUSOL(A,B,NR,N,PERM)
      
      IMPLICIT NONE

      INTEGER*4 NR,N
      REAL*8 A,B,PERM
      
      DIMENSION A(NR,1),B(1),PERM(1)

      REAL*8 X,DP
      INTEGER*4 I,K,IM1,INF
      
      IF(N.GT.1)GOTO 2
      B(1)=B(1)/A(1,1)
      RETURN
2     DO10I=1,N
      K=PERM(I)
      IF(K.EQ.I)GOTO10
      X=B(K)
      B(K)=B(I)
      B(I)=X
10    CONTINUE
      DO20I=2,N
      IM1=I-1
      DP=B(I)
      DO40K=1,IM1
      DP=DP-A(I,K)*B(K)
40    CONTINUE
      B(I)=DP
20    CONTINUE
      B(N)=B(N)/A(N,N)
      DO60I=2,N
      IM1=N-I+1
      INF=IM1+1
      DP=B(IM1)
      DO70K=INF,N
      DP=DP-A(IM1,K)*B(K)
70    CONTINUE
      B(IM1)=DP/A(IM1,IM1)
60    CONTINUE
      RETURN
      END


C*************************          DSCALL         *********************
C=    COMPILER (LINK=IBJ$)
      REAL*8 FUNCTION DSCALL(X,Y,N)
C
C     ESEGUE IL PRODOTTO SCALARE DI X PER Y
C
      IMPLICIT NONE
      INTEGER*4 N
      REAL*8 X, Y

      INTEGER*4 I
      DIMENSION X(1),Y(1)
      DSCALL=0.D0
      DO 1 I=1,N
      DSCALL=DSCALL+X(I)*Y(I)
1     end do
      RETURN
      END

C*************************          DVET          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE DVET(W,WV,IVER)
C
C     CALCOLA L'OPERATORE VETTORE DI UN VETTORE
C
      IMPLICIT NONE

      INTEGER*4 IVER
      REAL*8 W,WV

      REAL*8 H
      
      DIMENSION W(3),WV(9)
      
      H=DFLOAT(IVER)
      WV(1)= 0.D0
      WV(2)= H*W(3)
      WV(3)=-H*W(2)
      WV(4)= -WV(2)
      WV(5)= 0.D0
      WV(6)= H*W(1)
      WV(7)= -WV(3)
      WV(8)= -WV(6)
      WV(9)= 0.D0
      RETURN
      END


C*************************          LININT          ********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE LININT(XY,N,X,Y,PEND)
      
      IMPLICIT NONE
      
      INTEGER*4 N
      REAL*8 XY(2,1),X,Y,PEND

      INTEGER*4 K1, K2, I
      IF(X.GT.XY(1,1)) GO TO 10
      K1=1
      GO TO 30
  10  DO 20 I=2,N
      K1=I-1
      IF(XY(1,I).GT.X) GO TO 30
  20  CONTINUE
      K1=N-1
  30  K2=K1+1
      PEND=(XY(2,K2)-XY(2,K1))/(XY(1,K2)-XY(1,K1))
      Y=XY(2,K1)+PEND*(X-XY(1,K1))
      RETURN
      END

C*************************          MOVE          **********************
C=    COMPILER (LINK=IBJ$)
      SUBROUTINE MOVE(A,B,N)
C
C     TRASFERISCE N TERMINI DEL VETTORE 'B' NEL VETTORE 'A'
C
      IMPLICIT NONE

      INTEGER*4 N
      REAL*8 A,B
      DIMENSION A(1),B(1)

      INTEGER*4 I
      
      DO 10 I=1,N
      A(I)=B(I)
10    CONTINUE
      RETURN
      END


C*************************          PRM3            ********************
      SUBROUTINE PRM3  (A,NRDA,B,NRDB,C,NRDC,IA,IB)
C
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B"
C     IL PRODOTTO E' DI 3 RIGHE PER 3 COLONNE
C     A    = PRIMA MATRICE
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"
C     NCA  = N.RO DI COLONNE DI "A"
C     B    = SECONDA MATRICE
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"
C     C    = MATRICE PRODOTTO
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO
C                                            1 = SI
C     IB   = INDICE DI TRASPOSIZIONE DI "B"
C
      IMPLICIT NONE

      INTEGER*4 NRDA,NRDB,NRDC,IA,IB
      REAL*8 A,B,C

      INTEGER*4 NRA, NCA, NCC, IGO, I, K, L
      REAL*8 D
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1),D(3,3)
      
      NRA=3
      NCA=3
      NCC=3
      IGO=2*IA+IB+1
      GOTO(10,20,30,40),IGO
10    DO 110 I=1,NRA
      DO 111 K=1,NCC
      D(I,K)=0.D0
      DO 112 L=1,NCA
      D(I,K)=D(I,K)+A(I,L)*B(L,K)
112   CONTINUE
111   CONTINUE
110   CONTINUE
      GO TO 50
20    DO 120 I=1,NRA
      DO 121 K=1,NCC
      D(I,K)=0.D0
      DO 122 L=1,NCA
      D(I,K)=D(I,K)+A(I,L)*B(K,L)
122   CONTINUE
121   CONTINUE
120   CONTINUE
      GO TO 50
30    DO 130 I=1,NCA
      DO 131 K=1,NCC
      D(I,K)=0.D0
      DO 132 L=1,NRA
      D(I,K)=D(I,K)+A(L,I)*B(L,K)
132   CONTINUE
131   CONTINUE
130   CONTINUE
      GO TO 50
40    DO 140 I=1,NCA
      DO 141 K=1,NCC
      D(I,K)=0.D0
      DO 142 L=1,NRA
      D(I,K)=D(I,K)+A(L,I)*B(K,L)
142   CONTINUE
141   CONTINUE
140   CONTINUE
   50 DO 60 I=1,3
      DO 61 K=1,3
      C(I,K)=D(I,K)
   61 CONTINUE
   60 CONTINUE
      RETURN
      END



C*************************          DZERO            ********************
      SUBROUTINE DZERO(A,N)

      IMPLICIT NONE
      
      INTEGER*4 N,I
      REAL*8 A(N)

      DO I=1,N
        A(I)=0.D0
      END DO
      RETURN
      END

C*************************          POLCOE           ********************
      SUBROUTINE POLCOE(X,Y,N,COF)

      IMPLICIT NONE
      
      REAL*8 X(1),Y(1),COF(1)
      INTEGER*4 N
C
      REAL*8 PHI,B,FF,S(10)
      INTEGER*4 I,J,K
C
C Fixme: check for N <= 10?
      DO I=1,N
        S(I)=0.D0
        COF(I)=0.D0
      ENDDO
      S(N)=-X(1)
      DO I=2,N
        DO J=N+1-I,N-1
          S(J)=S(J)-X(I)*S(J+1)
        ENDDO
        S(N)=S(N)-X(I)
      ENDDO
      DO J=1,N
        PHI=N
        DO K=N-1,1,-1
          PHI=K*S(K+1)+X(J)*PHI
        ENDDO
        FF=Y(J)/PHI
        B=1.D0
        DO K=N,1,-1
          COF(K)=COF(K)+B*FF
          B=S(K)+X(J)*B
        ENDDO
      ENDDO
      RETURN
      END

