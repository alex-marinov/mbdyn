C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE WVEFR
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C Driver per il calcolo delle forze scambiate tra pneumatico e terreno

      SUBROUTINE WVEFR(VW, OME, CC, RES, RM, IRCW, IPC, OMEP, TABFAT,
     &  NPFT, BETA, BETAP, BLENG, PP, WP, GP)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I, IRCW, IPC, NPFT
      DIMENSION VW(6), CC(9), WP(12), GP(3), RES(7),CDP(3, 4),
     &  ZCG(6,3), VEL(3), TABFAT(6), PP(6)

C      common /forces/res1,res2,res3
C   Inizializza variabili
          
CCCCCC          print *,'cc(1->3)',(cc(i),i=1,3)
CCCCCC          print *,'cc(4->6)',(cc(i),i=4,6)
CCCCCC          print *,'cc(7->9)',(cc(i),i=7,9)

C DESCRIZIONE AREE DI LAVORO LOCALI
C
C  CC(9) DEFINIZIONE GEOMETRICA RUOTA
C 1-3    COSENI DIRETTORI ASSE RUOTA
C 4-6    POSIZIONE ORIGINE ASSE RUOTA
C 7      RAGGIO DELLA RUOTA
C 8-9    SEMISPESSORI DELLA RUOTA NELLA DIREZIONE ASSE
C
C CDP(12) DEFINIZIONE DEL TERRENO
C 1-9     COSENI DIRETTORI PIANO TERRENO
C 10-12   PUNTO DI RIFERIMENTO TERRENO
C
C  PP PROPRIETA' DEL PIANO DI CONTATTO
C     1 PRESSIONE RIFERIMENTO
C     2 VOLUME DI RIFERIMENTO
C     3 VELOCITA' DI RIFERIMENTO ARATURA
C     4 COEFFICIENTE RIDUTTIVO ATTRITO
C     5 COEFFICIENTE RIDUTTIVO ARATURA
C     6 COEFFICIENTE RIDUTTIVO DERIVA
C
C  WP PROPRIETA' PNEUMATICO
C     1 PRESSIONE DI RIFERIMENTO       press
C     2 VOLUME DI RIFERIMENTO          calcolato=2.*PIGRE^2*RC^2*(RR-RC)*Cvol
C     3 VELOCITA' DI RIFERIMENTO PER SMORZAMENTO  wrv
C     4 COEFFICIENTE ATTRITO MASSIMO (METTERE A ZERO)
C     5 AREA DI RIFERIMENTO PER LA DERIVA   arrif
C     6 COEFFICIENTE DERIVA                 rmust
C     7 COEFFICIENTE ARRETRAMENTO           xch
C     8 COEFFICIENTE RIDUTTIVO DELLA PRESSIONE PER CAMBER  cracmb
C     9 COEFFICIENTE STRISCIAMENTO DI ATTRITO MASSIMO      smumax
C    10 COEFFICIENTE DI ATTRITO MASSIMO                    rmum
C    11 COEFFICIENTE DI ATTRITO MINIMO PER MASSIMO STRISCIAMENTO rmus
C
C  GP PROPRIETA' GLOBALI
C     1 COEFFICIENTE DELLA POLITROPICA PER PRESSIONE PNEUMATICO   espr
C     2 VELOCITA' LIMITE PER IL CALCOLO DELLE FORZE               trv
C     3 LIBERO PER FUTURI USI   
C
C     RES(1->6) Forze in output
C     RES(7)    Schiacciamento pneumatico

      RM     = 0.     
C      PRESS  = 465000.
      RC     = CC(9)
      RR     = CC(7)

C      print *,'rc=',rc,' rr=',rr

!      COEVOL = .3
!      COEVOL = 1.
C      COEVOL = .8
C  coevol:vol rif ruota?
C      WRV    = 40.
C      RMUST  = 0.1
C      XCH    = 0.
C      CRACMB = 1.
C      SMUMAX = .075
C      RMUM   = .7
C      RMUS   = .1
C      ESPR   = 1.25
C      TRV    = .005       
      ARRIF  = 3.7*DSQRT(RR**2.-(RR-RC)**2.)*RC
      WP(5) = ARRIF
      
C      RK     = 1
      RK = GP(3)
C      
C Verificare se e' il caso di passarla tra i dati
      PI     = ACOS(-1.0)
      alfarif= pi/15.
  
      DO I = 1, 6
        RES(I) = 0.
      END DO

C inizializza la matrice CDP con un sistema di riferimento
C "alquanto arbitrario" e con (X, Y) del punto passato
      CALL SETPL(CC(4), CDP)

      CALL WHTOR(CC, CDP, IRCW, ZCG, VOL, PEN, ARRIF)

CCCCCC      print *,'zcg=',zcg
CCCCCC      print *,'vol=',vol
CCCCCC      print *,'pen=',pen
CCCCCC      print *,'arrif=',arrif

CCCCCC      print *,'ircw',ircw,' arrif',arrif
      
      IF (IRCW .EQ. 0) THEN
        RETURN
      END IF
      
C ?
      CALL VERCP(CC(4), VW, ZCG(4, 3), VEL)

C      print *,'vel=',vel

C ?
      CALL TYFRC(VEL, OME, ZCG, VOL, AP, CDP, PP, WP, GP, CC, RES,
     &  RM, IRCW, IPC, OMEP, RK, TABFAT, NPFT, BETA, BETAP, BLENG,
     &  RMUST, ALFARIF)
      res1 = res(1)
      res2 = res(2)
      res3 = res(3)

      res(7) = pen

CCCCCC                 print *,'res: ',(res(i),i=1,6)
      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE WHTOR
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE WHTOR (CC, CDP, IRCW, ZCG, VOL, DL, ARRIF)
      
      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 IRCW, I, K

      DIMENSION CC(9), CDP(3, 4), ZCG(6, 3), CR(3), TVET(9),
     &  X(3), A(3), R(2)

C ?
      COMMON /PNEU/ DR

C ------------------------------
  
      DO  I = 1, 3
        DO  K = 1, 6
          ZCG(K, I) = 0.
        END DO
      END DO

C posizione del centro ruota data da posizione del punto
C di riferimento (CC(I), I = 4, 6) piu' il punto medio
C della ruota, .5*(CC(8)+CC(9)) proiettato lungo l'asse ruota,
C di cui (CC(I), I = 1, 3) sono i coseni direttori
      DO I = 1, 3
        CR(I) = CC(I+3)+CC(I)*.5*(CC(8)+CC(9))
      END DO

C copia la colonna 3 di CDP negli ultimi elementi di TVET
      CALL DMOVE(CDP(1, 3), TVET(7), 3)

C mette all'inizio di TVET il prodotto vettore normalizzato
C degli ultimi elementi di TVET (colonna 3 di CDP) con
C l'asse ruota (CC(I), I = 1, 3)
      CALL PRVN(TVET(7), CC(1), TVET(1), FF)

C mette a meta' di TVET il prodotto vettore normalizzato
C di (TVET(I), I = 7, 9) e (TVET(I), I = 1, 3)
      CALL PRVN(TVET(7), TVET(1), TVET(4), FF)

      D1 = 0.
      D2 = 0.
      A1 = 0.
      A2 = 0.

      DO I = 1, 3
        D1 = D1-TVET(6+I)*CDP(I, 4)
        D2 = D2-CC(I)*CR(I)
        A1 = A1+CC(I)*TVET(3+I)
        A2 = A2+CC(I)*TVET(6+I)
      END DO

      A3 = 0.

      DO I = 1, 3
        X(I) = -D1*TVET(6+I)-(D2+A2*D1)*TVET(3+I)/A1
        A3   = A3+TVET(I)*(CR(I)-X(I))
      END DO

      DO I = 1, 3
        ZCG(3+I, 1) = X(I)+A3*TVET(I)
        DO K = 2, 3
          ZCG(3+I, K) = ZCG(3+I, 1)
        END DO
      END DO

      DR = 0.

      DO I = 1,3
        DR = DR+(ZCG(3+I, 1)-CR(I))**2
      END DO

      DR = DSQRT(DR)

      IF (DR .LT. CC(7)) THEN
        GOTO 35
      END IF

      IRCW = 0
      DL   = 0.

      RETURN

  35  IRCW = 1
      DL   = CC(7)-DR
      A(3) = ARRIF*DL*2./DABS(CC(9)-CC(8))
      R(1) = .5*(DABS(CC(8)-CC(9)))
      R(2) = CC(7)

      DO I = 1, 2
        DRL  = R(I)-DL
        DRC  = DSQRT(R(I)*R(I)-DRL*DRL)
        ALFA = DATAN2(DRC, DRL)*2.
        A(I) = .5*R(I)*R(I)*(ALFA-DSIN(ALFA))
      END DO

      DO K = 1, 3
        DO I = 1, 3
          ZCG(I, K) = A(K)*TVET((K-1)*3+I)
        END DO
      END DO

      VOL = A(3)*DL*.5

C      print *,'dr=',dr,'dl=',dl,'a(3)=',a(3)

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   29- 9-83                      SUBROUTINE SETPL      
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C prepara la matrice CDP, che contiene una matrice di rotazione nelle   
C prime tre colonne, e la posizione (solo X e Y) nella quarta colonna   

      SUBROUTINE SETPL(XYZ, CDP)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I
      DIMENSION XYZ(3), CDP(3, 4)

C -------------------------------------

      ALFA = 0.

C azzera la matrice CDP
      DO I = 1, 4
        CALL DSET(CDP(1, I), 3, 0.D0)
      END DO

C ma ALFA e' 0. !
      CDP(1, 1) = DCOS(ALFA)
      CDP(3, 1) = DSIN(ALFA)
      CDP(1, 3) = DSIN(ALFA)
      CDP(3, 3) = -DCOS(ALFA)

C mette nella colonna 2 di CDP il prodotto vettore tra le colonne 3 e 1
C nota: il risultato deve essere 0., -1., 0.
      CALL PRVU(CDP(1, 3), CDP(1, 1), CDP(1, 2))

C mette la posizione (X, Y) nella colonna 4 di CDP
      CDP(1, 4) = XYZ(1)
      CDP(2, 4) = XYZ(2)
      CDP(3, 4) = 0.

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE TYFRC
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE TYFRC(VEL, OME, ZCG, VOL, AP, CDP, PP, WP, GP, CC,
     &  RES, RM, IRCW, IPC, OMEP, RK, TABFAT, NPFAT, BETA, BETAP,
     &  BLENG, RMUST, ALFARIF)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 IRCW, IPC, NPFAT, IPCAD, I
      DIMENSION VEL(3), ZCG(6, 3), CDP(3, 4), PP(6), WP(12), GP(3),
     &  CC(9), RES(6), VEP(3), B(3), RD(3), RS(3), RR(3), F(6),
     &  TABFAT(1), VASS(2)

      COMMON /VELANG/ ALFA, VASS1, VASS2, CCY, SSY, SR, PIPPO     

C -------------------------------------

      OMEZER = 0.
      IPCAD  = 0
      BI     = 0.
      SR     = 0.

C ----------------------------------------------------
C     ELASTIC,DAMPING,FRICTION AND PLOWING FORCES     
C ----------------------------------------------------

      IF (IPC .EQ. 0) THEN
        GOTO 5
      END IF

      IF (IPCAD .EQ. 0) THEN
        OMEZER = OME
      END IF

      IPCAD = 1

   5  CALL FORCES(CC, VEL, ZCG, VOL, CDP, PP, WP, WP(8), GP, RES,
     &  CC(4), AP, FBET, IRCW)

C      print *,'res=',res
C      print *,'wp=',(wp(i),i=5,6)
      
      DO I = 1, 6
        F(I) = 0.
      END DO  

      IF(WP(5) .LE. 0. .OR. WP(6) .LE. 0.) THEN
CCCCCC      WRITE(*,*) ' FORZE TANGENTI RUOTA NON CALCOLATE !!',wp(5),wp(6)
        RETURN
      END IF

C -----------------------------------------
C     VELOCITY IN EQUIVALENT PLANE FRAME   
C -----------------------------------------

      CALL TRSFV(CDP, VEL, VEP, -1)

      VA   = DSQRT(VEP(1)**2+VEP(2)**2)
      AREA = DSQRT(ZCG(1, 3)**2+ZCG(2, 3)**2+ZCG(3, 3)**2)

      DO I = 1, 3
        B(I) = ZCG(I+3,3)-CC(I+3)
      END DO

      CALL PRVN(B, CC, RR, A)

c      DO 16 I=1,2 
c        VREL(I) = OME*A*RR(I)+VEL(I)
C        VREL(I) = OME*A*RR(I)+A*VEL(I)/.252 

c  16  CONTINUE 


C -----------------------
C     FRICTION FORCES   
C -----------------------

      IF(A .EQ. 0.) THEN
        GOTO 20
      END IF

      TH = 1.

      IF (IPC .EQ. 1 .AND. OMEZER .EQ. 0.) THEN

C ------------------------------------------------
C     PROVA DI CADUTA SENZA PRE-SPIN DELLA RUOTA   
C ------------------------------------------------

        RAGGIO = 0.
        DO I = 1, 3
          RAGGIO = RAGGIO+(ZCG(I+3, 3)-CC(I+3))**2
        END DO
        RAGGIO = DSQRT(RAGGIO)
        if (DABS(BETA*BETAP) .GT. 1.E-7) THEN
          OMER   = BLENG/RAGGIO*BETAP*DSIN(BETA)
          AAA    = OME/OMER
        ELSE
          AAA = 1
        ENDIF
        TH     = DTANH(DABS((OME-OMER)/70.))
        SR     = 1.-AAA
        SIG    = -1.
        IF (OME .GT. OMER) THEN
          SIG = 1.
        END IF 
      ELSE IF (IPC .NE. 1 .OR. OMEZER .NE. 0.) THEN

C --------------------------
C     TUTTI GLI ALTRI CASI   
C --------------------------

        OMEGA = 0.0
        IF(SR .EQ. 0. .OR. OME .EQ. OMEZER) THEN
          GOTO 30
        END IF
        BI = TANGE(SR, TABFAT, NPFAT)
        IF (OME .NE. 0.0) THEN
          IF (OME-OMEZER .NE. 0.0) THEN
            OMEGA = OME-OMEP*AREA*AP*BI/(RK*(OME-OMEZER))
          ELSE
            OMEGA = 0.0
          END IF 
        ELSE
          OMEGA = OME
        ENDIF
  30    OMER = 0.
        DO I = 1, 3
          OMER = OMER+VEL(I)*RR(I)
        END DO
C        OMER = OMER/A
        IF (DABS(OMER-OMEGA*CC(7)).EQ.0.D0) THEN
          SR = 0.
        ELSE
          SR   = (OMER+OME*A)/(OMER-OMEGA*CC(7))
        END IF
        SIG  = SIGN(1.,SR)
        TH   = 1.
      END IF

      SRO = SR
      SR  = DABS(SR)
      
      IF (SR .GT. 1.) THEN
        SR = 1.
      END IF
         
C      CALL RONT(SR,TABFAT,NPFAT,RANT)

Calcolo cosangolo tra piano ruota e asse x del sistema di riferimento
CC   terreno

      CCY = CC(2)/DSIN(DACOS(CC(3)))        
c           print *,cc(1),cc(2),cc(3)
Calcolo sinangolo tra piano ruota e asse x del sistema di riferimento
CC   terreno

      SSY = CC(1)/DSIN(DACOS(CC(3)))

CC   Attenzione vale solo se asse non ä normale alla pista. in
CC   cosbeta (denominatore) non
CC   ce il segno perchä il suo valore ä positivo indipententemente
CC   dalla sua posizione

Calcolo angolo tra piano ruota e asse x del sistema di riferimento
CC   terreno

      AY = DATAN2(SSY, CCY)

CC   Non ci sono problenmi se l'asse mozzo non ä orientato come x

      VASS(1) = -VEL(1)*CCY+VEL(2)*SSY
      VASS(2) = -VEL(2)*CCY+VEL(1)*SSY
C           vassabs = sqrt(vass(1)**2+vass(2)**2)

Calcolo angolo di deriva

      VA1 = DABS(VASS(1))
      VASS1 = VASS(1)
      VA2 = DABS(VASS(2))
      VASS2 = VASS(2)
      RIF = 0.000001
      IF (VA1 .GE. RIF) THEN
        ALFA = DATAN2(VASS(2), VASS(1))
      ELSE IF (VA1 .LT. RIF .AND. VA2 .GE. RIF) THEN
        ALFA = DATAN2(VASS(1), VASS(2))
        ALFA = SIGN(1., VASS(2))*(PI/2.-ALFA)
      ELSE IF (VA1 .LT. RIF .AND. VA2 .LT. RIF) THEN
        ALFA = 0.
      END IF

C      print *,'prima di MUX:',va1,va2,sr,alfa
      
      CALL MUX(SR, ALFA, TABFAT, NPFAT, RANT)

      FRIC = TH*RANT      

      Fx = -SIG*PP(4)*FRIC*AREA*AP
      F(1) = FX*(-CCY)
      F(2) = FX*(SSY)
      F(3) = 0.
      
C      print *,'before tmmnt:',res    
C      print *,FX,CCY,SSY,FRIC,TH,RANT
      
      CALL TMMNT(ZCG(4, 3), F, CC(4), RES)
      
C      print *,'after tmmnt: ',res

  20  CONTINUE

C --------------------
C     DRIVING FORCE    
C --------------------

      IF (VA .GT. GP(2)) THEN
        IF(A .GT. 0.) THEN
          GOTO 45
        END IF

        DFCR = .7*WP(6)*PP(6)
        CD   = DFCR*AP*AREA
      
        DO I = 1, 3
          F(I) = -CD*(VEP(1)*CDP(I, 1)+VEP(2)*CDP(I, 2))/VA
        END DO
        
        GOTO 55

  45    CALL PRVN(CDP(1, 3), RR, RD, A)

        CALL PRVN(RD, CDP(1, 3), RS, A)

        CALL TRSFV(CDP, RD, RR, -1)


c      SINA  = (VREL(1)*RR(1)+VREL(2)*RR(2))/SQRT(VREL(1)**2+VREL(2)**2)     
        VEP1 = VEP(1)
        VEP2 = VEP(2)
      
        SINA = (VEP(1)*RR(1)+VEP(2)*RR(2))/VA
      
        COSA  = DSQRT(1.-SINA**2)
c      ANDER = ASIN(SINA)
        ANDER = DATAN2(SINA, COSA)
        BB    = 1.32*AREA/WP(5)+2.E-3
        BR    = .5*(BB+1.75-DSQRT((BB-1.75)**2+1.E-2))
        BR    = (2.-BR)*BR
        FD    = BR/(BB+1.E-4)
        ANDERVAL = ANDER
        
        ANDER = FD*ANDER

        DFCR  = WP(6)*PP(6)
      
C ------------------------------------------------------        
c      CD    = TANH(vep(2)/GP(2))*DFCR*AP*AREA*ABS(TANH(.7*ANDER/DFCR))
C ------------------------------------------------------

        CD = DFCR*AP*AREA*DTANH(.7*ANDER/DFCR)
        PIPPO = DTANH(ALFA/ALFARIF)
c      print *,'pvep1,vep2,va,sina,cosa', vep(1),vep(2),va,sina,cosa 
        CALL MUY(SR, ALFA, ALFARIF, RMUST, RUNT)

        FY = -RUNT*AP*AREA
        F(1) = FY*(-SSY)
        F(2) = FY*(-CCY)
        F(3) = 0.
c      print *,'Fy',fy,alfa,vass(2)

        A = 2.*CC(7)*WP(7)*DTANH(.7*DABS(ANDER))*COSA

        DO I = 1, 3
          B(I)=(ZCG(I+3,3)-RS(I)*A*SIG)*FBET
        END DO
        
  55    CALL TMMNT(B, F, CC(4), RES)
      ELSE
CCCCCC        write(*,*) ' Calcoli non effettuati',va,gp(2)
      END IF

C ----------------------------------
C     CALCOLO MOMENTO SULLA RUOTA 
C ----------------------------------

      DO I = 1, 3
	RM = RM+RES(I+3)*CC(I)
      END DO

      RETURN
      END


C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE FORCES
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE FORCES(CC, VEL, ZCG, VOL, CDP, PP, OP, RED, GP,
     &  RES, RP, AP, FBET, IRCW)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 IRCW, I, IGO, J, K
      DIMENSION CC(9), VEL(3), ZCG(6, 3), CDP(3, 4), PP(6), OP(4),
     &  GP(3), RES(6), VEP(3), F(6)

      DATA PGRQ/.785398D0/

C --------------------------------------

      SB = 0.

      DO I = 1, 3
        SB = SB+CC(I)*CDP(I, 3)
      END DO

      SB   = DABS(SB)
      CB   = DSQRT(1.-SB**2)
      BET  = DATAN2(SB, CB)
      FBET = 0.

      IF (BET .LT. PGRQ) THEN
        FBET = DCOS(2.*BET)
      END IF

      OPR = OP(1)*(RED+(1.-RED)*FBET)

C ------------------------
C     ACTUAL PRESSURE        
C ------------------------

      GM1 = 1./GP(1)
      RAP = (OPR/PP(1))**GM1
      RAP = RAP*(OP(2)/PP(2))
      X   = (VOL+RAP*PP(2)-OP(2))/(1.+RAP)

      IF (X .GT. 0. .AND. X .LT. VOL) THEN
        GOTO 30
      END IF

      IF (X .GT. 0.) THEN
        GOTO 20
      END IF

  10  IGO  = 2
      IRCW = 1
      X    = 0.
      AP   = OPR*(OP(2)/(OP(2)-VOL))**GP(1)

C      print *,opr,op(2),vol,gp(1)

      GOTO 40

  20  IGO  = 1
      IRCW = 3
      X    = VOL
      AP   = PP(1)*(PP(2)/(PP(2)-VOL))**GP(1)

      GOTO 40

  30  IGO  = 1
      IRCW = 2
      AP   = OPR*(OP(2)/(OP(2)-VOL+X))**GP(1)

  40  CONTINUE

C ---------------------------------
C     ZERO-ING RESULTANT FORCE        
C ---------------------------------

      CALL DSET(RES, 6, 0.D0)

      CALL DSET(F, 6, 0.D0)
C ----------------------
C     ELASTIC FORCE          
C ----------------------

      DO I = 1, 3
        F(I) = ZCG(I, 3)*AP
      END DO

      CALL TMMNT(ZCG(4, 3), F, RP, RES)

C ----------------------------------------------------
C     VELOCITY IN EQUIVALENT PLANE FRAME         
C ----------------------------------------------------

      CALL TRSFV(CDP, VEL, VEP, -1)

      VA = DSQRT(VEP(1)**2+VEP(2)**2)
      VM = DSQRT(VEP(1)**2+VEP(2)**2+VEP(3)**2)

C ---------------------------------
C             DAMPING FORCE           
C ---------------------------------

      DC = DTANH(VEP(3)/OP(3))

      DO I = 1, 3
        F(I) = -DC*ZCG(I, 3)*AP
      END DO

      CALL TMMNT(ZCG(4, 3), F, RP, RES)
C -------------------------
C     FRICTION FORCE          
C -------------------------

      IF (VA .LT. GP(2)) THEN
        GOTO 55
      END IF

      AN = DSQRT(ZCG(1, 3)**2+ZCG(2, 3)**2+ZCG(3, 3)**2)
      FC = OP(4)*PP(4)*AN*AP/VA

      DO I = 1, 3
        F(I) = -FC*(VEP(1)*CDP(I, 1)+VEP(2)*CDP(I, 2))
      END DO

      CALL TMMNT(ZCG(4, 3), F, RP, RES)

  55  CONTINUE

C ---------------------
C     PLOWING FORCE          
C ---------------------

      IF (IGO .EQ. 2) THEN
        GOTO 65
      END IF

      IF (VA .LT. GP(2)) THEN
        GOTO 65
      END IF

      IF (VOL .LE. 0.) THEN
        GOTO 65
      END IF

      PC = X*2.*AP*PP(5)*DTANH(VA/PP(3))/(VA*VOL)

      DO J = 1, 2
        A1 = 0.
        DO K = 1, 3
          A1 = A1+ZCG(K, J)*VEL(K)
        END DO
        DO I = 1, 3
          F(I) = -PC*A1*CDP(I, J)
        END DO
        CALL TMMNT(ZCG(4, J), F, RP, RES)
      END DO

  65  CONTINUE

      RETURN
      END

C -------------------------------------------------------------------
C -------------------------------------------------------------------
C                 30- 1-82                         FUNCTION TANGE
C -------------------------------------------------------------------
C -------------------------------------------------------------------

      REAL*8 FUNCTION TANGE(X, XV, NXV)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 NXV, I
      DIMENSION XV(1)

C -----------------------------------

      IF (NXV .GT. 1) THEN
        DO I = 1, NXV
          IF (X .LT. XV(I)) THEN
            GOTO 10
          END IF
        END DO
        I = NXV
  10    IF (I .EQ. 1) THEN
          I = I+1
        END IF
        TANGE = (XV(I+NXV)-XV(I+NXV-1))/(XV(I)-XV(I-1))
      ELSE
        TANGE = 0.
      END IF

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE VERCP
C ----------------------------------------------------------------------
C ---------------------------------------------------------------------

      SUBROUTINE VERCP(O, V, P, VC)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I
      DIMENSION O(3), P(3), V(6), VC(3), PMO(3)

C --------------------------------------------

      DO I = 1, 3
        PMO(I) = P(I)-O(I)
c        print *,i,Pmo(1),v(i),v(3+i)  
      END DO
      
      CALL PRVU(V(4), PMO, VC)

      DO I = 1, 3
        VC(I) = VC(I)+V(I)
      END DO

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE PRVN
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C ritorna in T il prodotto vettore tra X e Y normalizzato;              
C in D mette il modulo di X vettore Y                                   

      SUBROUTINE PRVN(X, Y, T, D)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I
      DIMENSION X(3), Y(3), T(3), Z(3)

C --------------------------------------

C mette in Z il prodotto vettore tra X e Y
      CALL PRVU(X, Y, Z)
C calcola la norma di Z
      D = Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3)

      IF (D .GT. 1.E-12) THEN
        GOTO 5
      END IF

      D = 0.

      RETURN

   5  D = DSQRT(D)

C se la norma di Z non e' nulla, la usa per normalizzare Z
      DO I = 1, 3
        T(I) = Z(I)/D
      END DO

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE DMOVE
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C copia N valori di RVI in RVO
      
      SUBROUTINE DMOVE(DVI, DVO, N)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I, N
      DIMENSION DVI(1), DVO(1)

C -----------------------------------------

      DO I = 1, N
        DVO(I) = DVI(I)
      END DO

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE RONT
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE RONT(X, XYT, NP, RANT)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I, NP, IM1 
      DIMENSION XYT(6)

C --------------------------------------

      IF (NP .GT. 1) THEN
        GOTO 10
      END IF

      RANT = XYT(2)

      RETURN

  10  CONTINUE

      DO I = 2, NP
        IF (X .LT. XYT(I)) THEN
          GOTO 30
        END IF
      END DO

      I    = NP
  30  IM1  = I-1
      RANT = (XYT(NP+I)-XYT(NP+IM1))/(XYT(I)-XYT(IM1))*(X-XYT(IM1))
      RANT = XYT(NP+IM1)+RANT

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE DSET
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE DSET(DV, N, DVAL)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I, N 
      DIMENSION DV(1)

C ---------------------------------------

      DO I = 1, N
        DV(I) = DVAL
      END DO

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE TMMNT
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE TMMNT(PA, FA, PB, FB)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I
      DIMENSION PA(3), PB(3), FA(6), FB(6), T(3)

C -------------------------------------

      DO I = 1, 3
        FB(I) = FB(I)+FA(I)
        T(I)  = PA(I)-PB(I)
      END DO
      CALL PRVU(T, FA, T)
      DO I = 1, 3
        FB(I+3) = FB(I+3)+T(I)
      END DO

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE TRSFV
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE TRSFV(TM, RV, RN, IVER)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 IVER, ITR
      DIMENSION TM(3, 1), RV(1), RN(1)

C ------------------------------------------

      ITR = IABS((IVER-1)/2)

      CALL MV3(TM, RV, RN, ITR)

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE PRVU
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C restituisce in T il prodotto vettore tra X e Y                        

      SUBROUTINE PRVU(X, Y, T)

      IMPLICIT REAL*8 (A-H, O-Z)

      DIMENSION X(3), Y(3), T(3)

C ------------------------------------------

      Z1   = X(2)*Y(3)-X(3)*Y(2)
      Z2   = X(3)*Y(1)-X(1)*Y(3)
      Z3   = X(1)*Y(2)-X(2)*Y(1)
      T(1) = Z1
      T(2) = Z2
      T(3) = Z3

      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE MV3
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------

      SUBROUTINE MV3(RM, VIN, VOUT, IT)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 IT
      DIMENSION RM(9), VIN(3), VOUT(3)

C ----------------------------------------------

      V1 = VIN(1)
      V2 = VIN(2)
      V3 = VIN(3)

      IF (IT .GT. 0) THEN
        GOTO 10
      END IF

      VOUT(1) = RM(1)*V1+RM(4)*V2+RM(7)*V3
      VOUT(2) = RM(2)*V1+RM(5)*V2+RM(8)*V3
      VOUT(3) = RM(3)*V1+RM(6)*V2+RM(9)*V3

      RETURN

  10  VOUT(1) = RM(1)*V1+RM(2)*V2+RM(3)*V3
      VOUT(2) = RM(4)*V1+RM(5)*V2+RM(6)*V3
      VOUT(3) = RM(7)*V1+RM(8)*V2+RM(9)*V3

      RETURN
      END
      
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                                              sub mux
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------


      SUBROUTINE MUX(X, ALFA, XYT, NP, RANT)

      IMPLICIT REAL*8 (A-H, O-Z)

      INTEGER*4 I, NP, IM1 
      DIMENSION XYT(6)

C --------------------------------------
      CRID = 1.-2./DACOS(-1.D0)*ALFA

      IF (NP .GT. 1) THEN
        GOTO 10
      END IF

      RANT = XYT(2)

      RETURN      

  10  CONTINUE

      DO I = 1, NP
        XYT(I+NP) = XYT(I+NP)*CRID
c      print *,'xyt',I,xyt(I+np)
      END DO      
      DO I = 2, NP
        IF (X .LT. XYT(I)) THEN
          GOTO 30
        END IF
      END DO

      I    = NP
  30  IM1  = I-1
      RANT = (XYT(NP+I)-XYT(NP+IM1))/(XYT(I)-XYT(IM1))*(X-XYT(IM1))
      RANT = XYT(NP+IM1)+RANT
c      print *,'rant',rant
      
      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                                              sub muy
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------


      SUBROUTINE MUY(X, ALFA, ALFARIF, RMUST, RUNT)

      IMPLICIT REAL*8 (A-H, O-Z)

C --------------------------------------
crod Dipende da SR

      
c      crod = 1.-x
      CROD1 = ALFA/(DACOS(-1.D0)/2.D0)*RMUST
      CROD0 = RMUST*DTANH(ALFA/ALFARIF)
      RUNT = CROD0+((CROD1-CROD0)/1.)*X
c      runt = 0.
c      print *,'runt',runt
      
      RETURN
      END
