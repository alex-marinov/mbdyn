C $Header$
C MBDyn (C) is a multibody analysis code.
C http://www.mbdyn.org
C
C Copyright (C) 1996-2012
C    
C Pierangelo Masarati  <masarati@aero.polimi.it>
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

      program ruota

      IMPLICIT REAL*8 (A-H,O-Z)
      
      DIMENSION RES(6),cc(9),tabfat(6),vw(6)
      read(*,'(a)') aa
      read(*,*) cc
      read(*,'(a)') aa
      read(*,*) (tabfat(io),io=1,6)
      read(*,'(a)') aa
      read(*,*) vw
      read(*,'(a)') aa
      read(*,*) ipc
      read(*,'(a)') aa
      read(*,*) ome, omep
      read(*,'(a)') aa
      read(*,*) bleng,beta,betap
      npft=3
      call WVEFR(VW,OME,CC,RES,RM,IRCW,IPC,OMEP,TABFAT,
     1                 NPFT,BETA,BETAP,BLENG)
      stop
      end
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE WVEFR
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE WVEFR(VW,OME,CC,RES,RM,IRCW,IPC,OMEP,TABFAT,
     &  NPFT,BETA,BETAP,BLENG)
      
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  VW VELOCITA' LINEARI E ANGOLARI PERNO IN COORDINATA SUOLO
C
C OME VELOCITA' ANGOLARE RUOTA (TETA PUNTO)
C
C  CC DEFINIZIONE GEOMETRICA RUOTA
C 1-3   COSENI DIRETTORI ASSE RUOTA
C 4-6   POSIZIONE ORIGINE ASSE RUOTA
C 7-8   SEMISPESSORI DELLA RUOTA NELLA DIREZIONE ASSE
C 9     RAGGIO DELLA RUOTA
C
C CDP DEFINIZIONE DEL TERRENO
C 1-9   COSENI DIRETTORI PIANO TERRENO
C 10-12 PUNTO DI RIFERIMENTO TERRENO
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
C     1 PRESSIONE DI RIFERIMENTO
C     2 VOLUME DI RIFERIMENTO
C     3 VELOCITA' DI RIFERIMENTO PER SMORZAMENTO
C     4 COEFFICIENTE ATTRITO MASSIMO (METTERE A ZERO)
C     5 AREA DI RIFERIMENTO PER LA DERIVA
C     6 COEFFICIENTE DERIVA
C     7 COEFFICIENTE ARRETRAMENTO
C     8 COEFFICIENTE RIDUTTIVO DELLA PRESSIONE PER CAMBER
C     9 COEFFICIENTE STRISCIAMENTO DI ATTRITO MASSIMO
C    10 COEFFICIENTE DI ATTRITO MASSIMO
C    11 COEFFICIENTE DI ATTRITO MINIMO PER MASSIMO STRISCIAMENTO
C
C  GP PROPRIETA' GLOBALI
C     1 COEFFICIENTE DELLA POLITROPICA PER PRESSIONE PNEUMATICO
C     2 VELOCITA' LIMITE PER IL CALCOLO DELLE FORZE
C     3 LIBERO PER FUTURI USI
C
C
cc      INCLUDE 'COMM2.FOR'
C
      DIMENSION VW(6),CC(9),WP(12),GP(3),RES(6),CDP(3,4),
     1          ZCG(6,3),VEL(3),TABFAT(1),pp(6)
      RM=0.
      
      PRESS =465000      
C      RC    =.087
      RC    = CC(9)
C      RR    =.252
      RR     = CC(7)
      COEVOL=.8
      WRV   =20
      RMUST =.075
      XCH   =0.
      CRACMB=1.
      SMUMAX=.075
      RMUM  =.7
      RMUS  =.25
      ESPR  =1.1
      TRV   =.005
C      ARRIF =0.032672
C 2.*M_PI*M_PI*dRaggioToro*dRaggioToro*(p->cc[6]-dRaggioToro)*coevol;
      ARRIF = 3.7*DSQRT(RR**2.-(RR-RC)**2.)*RC
      RK    =1

      WP(1)=PRESS
      WP(2)=2.*3.1415926d0*3.1415926d0*RC*RC*(RR-RC)
      WP(2)=WP(2)*COEVOL
      WP(3)=WRV
      WP(4)=0.
      WP(5)=ARRIF
      WP(6)=RMUST
      WP(7)=XCH
      WP(8)=CRACMB
      WP(9)=SMUMAX
      WP(10)=RMUM
      WP(11)=RMUS
      GP(1)=ESPR
      GP(2)=TRV  
      gp(3)=0.  
      PP(1)=1.e10
      PP(2)=.3
      PP(3)=.5
      PP(4)=1.
      PP(5)=.8
      PP(6)=1.
C      WRITE(*,*)'WVEFR',gp
      CALL RSET(RES,6,0.D0)
C     WRITE(IOUT,*)' DG0=',DG0
C     --------------------------------------------------
      CALL SETPL(CC(4),CDP)
C     --------------------------------------------------
C     WRITE(IOUT,*)((CDP(II,KK),KK=1,4),II=1,3)
C     ------------------------------------------------
      CALL WHTOR(CC,CDP,IRCW,ZCG,VOL,PEN,ARRIF)

C     ------------------------------------------------
      IF(IRCW.eq.0) then
      write(*,*) ' Assenza di contatto ',pen
      RETURN
      endif
      
C     ---------------------------------
      CALL VERCP(CC(4),VW,ZCG(4,3),VEL)
C     ---------------------------------
C     WRITE(IOUT,333)VEL
C333  FORMAT(' VEL ',3E10.4)
C     --------------------------------------------------------------
      CALL TYFRC(VEL,OME,ZCG,VOL,AP,CDP,PP,WP,GP,CC,RES,RM,IRCW,IPC,
     1           OMEP,RK,TABFAT,NPFT,BETA,BETAP,BLENG)
C     --------------------------------------------------------------
C     WRITE(IOUT,89)AP,VOL,RM
89    FORMAT('  PRESS.',E12.5,'  VOL.',E12.5,'  MOM.',E12.5)
C      WRITE(6,*)'RETURN WVEFR' 
      write(*,'(a,6E10.4)') ' res ',(RES(J),J=1,6)
      write(*,'(a,1E10.4)') ' pen ',(pen)
      write(*,'(a,1E10.4)') ' vol ',(vol)   
      write(*,'(a,i8)')     ' irc ',(ircw)      
      RETURN
      END
C     COMPILER (ARGCHK=OFF)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE WHTOR
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE WHTOR (CC,CDP,IRCW,ZCG,VOL,DL,ARRIF)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION CC(9),CDP(3,4),ZCG(6,3),CR(3),TVET(9),X(3),A(3),R(2)
C    
C     INGRESSI
C     CC DEFINIZIONE GEOMETRICA RUOTA
C     1-3   COSENI DIRETTORI ASSE RUOTA
C     4-6   POSIZIONE ORIGINE ASSE RUOTA
C     8-9   SEMISPESSORI DELLA RUOTA NELLA DIREZIONE ASSE
C     7     RAGGIO DELLA RUOTA
C
C     CDP DEFINIZIONE DEL TERRENO
C     1-9   COSENI DIRETTORI PIANO TERRENO
C     10-12 PUNTO DI RIFERIMENTO TERRENO
C     ARRIF Area ellissoide di riferimento per il pneumatico
C 
C     USCITE
C     IRCW  Flag di contatto (0 no, 1 si)
C     ZCG   ZCG(1-3,1-3) coseni direttori di .... ?
C           ZCG(4-6,i)   coordinate punto di contatto
C     VOL   Volume di contatto
C     DL    Entita' della compenetrazione
C

      WRITE(6,*)'ENTRO IN WHTOR'
      print'(A,3e10.4)',' cc(1-3) ',cc(1),cc(2),cc(3)
      print'(A,3e10.4)',' cc(4-6) ',cc(4),cc(5),cc(6)
      print'(A,3e10.4)',' cc(7-9) ',cc(7),cc(8),cc(9)
      print*,' cdp'          
      do l=1,3
      print*,(cdp(l,i),i=1,4)
      enddo
C
C     Calcolo posizione centro ruota
C
      DO 5 I=1,3
        CR(I)=CC(I+3)+CC(I)*.5*(CC(8)+CC(9))
        print*,i,cc(i+3),cc(i),.5*(CC(8)+CC(9))
5     CONTINUE
      print*,' cr ',cr
C
C     Calcolo versori terna terreno con asse y complanare con y ruota
C     (il terreno puo' essere comunque oriantato)
C
      CALL RMOVE(CDP(1,3),TVET(7),3)
      CALL PRVN(TVET(7),CC,TVET(1),FF)
      CALL PRVN(TVET(7),TVET(1),TVET(4),FF)
      print*,' tvet'
      print*,tvet(1),tvet(4),tvet(7)
      print*,tvet(2),tvet(5),tvet(8)
      print*,tvet(3),tvet(6),tvet(9)
C
C     D1=
C     D2=
C     A1=angolo relativo terreno ruota
C     A2=angolo relativo terreno ruota
C
      D1=0.
      D2=0.
      A1=0.
      A2=0.
      DO 10 I=1,3
        D1=D1-TVET(6+I)*CDP(I,4)
        D2=D2-CC(I)*CR(I)
        A1=A1+CC(I)*TVET(3+I)
        A2=A2+CC(I)*TVET(6+I)
10    CONTINUE                     
      print*,' d1 ',d1
      print*,' d2 ',d2
      print*,' a1 ',a1
      print*,' a2 ',a2
C
C     A3 e' ?
C     X  e' il vettore
C
      A3=0.
      DO 15 I=1,3
        X(I)=-D1*TVET(6+I)-(D2+A2*D1)*TVET(3+I)/A1
        A3=A3+TVET(I)*(CR(I)-X(I))
15    CONTINUE        
      print*,' x  ',x
      print*,' a3 ',a3
C
C     Calcolo posizione del punto di contatto a centro ruota
C     (le tre colonne contengono lo stesso vettore)
C
      DO 20 I=1,3
        ZCG(3+I,1)=X(I)+A3*TVET(I)
        DO 20 K=2,3
          ZCG(3+I,K)=ZCG(3+I,1)
20    CONTINUE
C
C     Calcolo intersezione 
C
      DR=0.
      DO 25 I=1,3
        DR=DR+(ZCG(3+I,1)-CR(I))**2
25    CONTINUE
      DR=DSQRT(DR)
      print*,' dr ',dr,'    <cc(7)',cc(7),'?'      
      IF(DR.ge.CC(7)) then
        IRCW=0
        DL=0.
        RETURN              
      else
        IRCW=1
        DL=CC(7)-DR                            
      print*,' dl ',dl
C
C     Calcolo aree di impronta
C     A(1) sezione frontale dell'intersezione toro/piano
C     A(2) sezione laterale dell'intersezione toro/piano
C     A(3) sezione impronta: calcolata come variabile linearmente da 0 a ARRIF
C
        A(3)=ARRIF*DL*2./DABS(CC(9)-CC(8))
        R(1)=.5*(DABS(CC(8)-CC(9)))
        R(2)=CC(7) 
      print*,' A(3) ',a(3)
      print*,' R    ',r      
C
C     le sezioni frontali e laterali sono intersezioni di cerchi con un piano
C     DL   e' la porzione di raggio intersecata dal piano
C     DRL  e' la porzione di raggio tra piano e centro del cerchio
C     DRC  e' la corda dell'intersezione
C     ALFA e' l'angolo sotteso dalla corda
C     A(I) e' l'area della calotta circolare I=1 sezione frontale
C                                            I=2 sezione laterale
C
        DO 30 I=1,2
          DRL=R(I)-DL
          DRC=DSQRT(R(I)*R(I)-DRL*DRL)
          ALFA=DATAN2(DRC,DRL)*2.
          A(I)=.5*R(I)*R(I)*(ALFA-DSIN(ALFA))
      print*,' drl ',drl,'     drc',drc
      print*,' alfa      ',alfa*180/3.1415
      print*,' i,A(i)    ',i,a(I)
   30   CONTINUE     
C
C     proiezione delle sezioni di contatto nelle direzioni coordinate
C     (il piano di contatto puo' essere comunque orientato)
C
        DO 35 K=1,3
          DO 35 I=1,3
            ZCG(I,K)=A(K)*TVET((K-1)*3+I)
   35   CONTINUE                            
C
C     Calcolo volume di intersezione: Vol=Aimpronta*Dl/2
C
        VOL=A(3)*DL*.5

        print'(A,e10.4)', ' arrif   ',arrif
        print'(A,e10.4)', ' A(3)    ',A(3)
        print'(A,2e10.4)',' dl,dr   ',dl,dr
        print*,' vol  ',vol
        do i=1,6
        print'(A,3e10.4)',' zcg ',(zcg(i,k),k=1,3)
        enddo
      endif
      WRITE(6,*)'RETURN WHTOR' 
      RETURN
      END
C     COMPILER (ARGCHK=OFF)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   29- 9-83                      SUBROUTINE SETPL
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE SETPL(XYZ,CDP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XYZ(3),CDP(3,4)
      ALFA= 0.*3.1415927d0/180.
120   CALL RSET(CDP,12,0.D0)
      CDP(1,1)=DCOS(ALFA)
      CDP(3,1)=DSIN(ALFA)
      CDP(1,3)=DSIN(ALFA)
      CDP(3,3)=-DCOS(ALFA)
      CALL PRVU(CDP(1,3),CDP(1,1),CDP(1,2))
      CDP(1,4)=XYZ(1)
      CDP(2,4)=XYZ(2)
      CDP(3,4)=0
      RETURN
      END
C     COMPILER (ARGCHK=OFF)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE TYFRC
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE TYFRC(VEL,OME,ZCG,VOL,AP,CDP,PP,WP,GP,CC,RES,RM,IRCW,
     1                 IPC,OMEP,RK,TABFAT,NPFAT,BETA,BETAP,
     2                 BLENG)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     VEL        VELOCITY OF THE CENTRE OF THE CONTACT AREA
C     OME        SPIN OF THE WHEEL
C     ZCG        ORIENTED AREA AND CENTRE OF GRAVITY
C     VOL        CONTACT VOLUME
C     CDP        MATRIX OF THE DIRECTION COSINES OF THE EQUIVALENT PLANE
C     PP         PLANE  PARAMETERS
C     WP         WHEEL  PARAMETERS
C     GP         GLOBAL PARAMETERS
C     CC         DEFINIZIONE GEOMETRICA RUOTA
C                1-3 : COSENI DIRETTORI ASSE RUOTA
C                4-6 : POSIZIONE ORIGINE ASSE RUOTA
C                7-8 : SEMISPESSORI DELLA RUOTA NELLA DIREZIONE ASSE
C                  9 : RAGGIO DELLA RUOTA
C     RES        1-3 : CONTACT REACTION FORCE
C                4-5 : CONTACT REACTION MOMENT
C     IRCW =1    DEFORMATION OF TYRE ONLY
C          =2    DEFORMATION OF TYRE AND PLANE
C          =3    DEFORMATION OF PLANE ONLY
C     IPC  =1    PROVA DI CADUTA CARRELLO
C          =0    PROVE CON Xcg NON NULLO
C
      
      DIMENSION VEL(3),ZCG(6,3),CDP(3,4),PP(6),WP(12),GP(3),CC(9),
     1          RES(6),VEP(3),B(3),RD(3),RS(3),RR(3),F(6),
     2          TABFAT(1)
      DATA OMEZER/0./,IPCAD/0/,BI/0./,SR/0./
C    ------------------------------------------------------
C    |    ELASTIC,DAMPING,FRICTION AND PLOWING FORCES     |
C    ------------------------------------------------------
      WRITE(6,*)'ENTRO IN TYFRC',ipc,ipcad,ome,omezer
      IF(IPC.EQ.0) GOTO 1
      IF(IPCAD.EQ.0) OMEZER=OME
      IPCAD=1
1     CALL FORCES(CC,VEL,ZCG,VOL,CDP,PP,WP,WP(8),GP,RES,CC(4),AP,FBET,
     +           IRCW)
      CALL RSET(F,6,0.D0)
      write(*,'(A,3e12.5)') ' FFORCES ',(RES(iol),iol=1,3)
      write(*,'(A,2e12.5)') ' wp(5,6) ',wp(5),wp(6)
      IF(WP(5).LE.0. .OR. WP(6).LE.0.) THEN
      WRITE(*,*) ' FORZE TANGENTI RUOTA NON CALCOLATE !! vedi dati:',
     *  ' wp(5 o 6)<=0 ',wp(5),wp(6)
      RETURN                                                      
      endif
C    -------------------------------------------
C    |    VELOCITY IN EQUIVALENT PLANE FRAME   |
C    -------------------------------------------
      CALL TRSFV(CDP,VEL,VEP,-1)
C     modulo velocita' impronta
      VA=DSQRT(VEP(1)**2+VEP(2)**2)
      do iolp=1,3
      WRITE(6,'(a,5e10.4)')' CDP ',(CDP(IOL,iolp),iol=1,3)
      enddo
      WRITE(6,'(a,5e10.4)')' VEL ',(VEL(IOL),iol=1,3)
      WRITE(6,'(a,5e10.4)')' VEP ',(VEP(IOL),iol=1,3)
      WRITE(6,'(a,5e10.4)')' VA  ',VA
      IF(VA.LT.GP(2)) then
      write(*,*) ' Calcoli non effettuati VA.lt.GP(2)',va,gp(2)
      RETURN
      endif
      AREA=DSQRT(ZCG(1,3)**2+ZCG(2,3)**2+ZCG(3,3)**2)
      DO 2 I=1,3
       B(I)=ZCG(I+3,3)-CC(I+3)
2     CONTINUE
      CALL PRVN(B,CC,RR,A)
      print*,' RR=B x CC,  A=|RR|'
      WRITE(6,'(a,5e10.4)')' B   ',(B(IOL),iol=1,3)
      WRITE(6,'(a,5e10.4)')' CC  ',(CC(IOL),iol=1,3)
      WRITE(6,'(a,5e10.4)')' RR  ',(RR(IOL),iol=1,3)
      WRITE(6,'(a,5e10.4)')' A   ',A
C     -------------------------
C     |    FRICTION FORCES    |
C     -------------------------
      IF(A.EQ.0.D0) GOTO 100
      TH=1.
      print*,ipc,omezer,ome
      IF (IPC.EQ.1.AND.OMEZER.EQ.0.D0) THEN
C          --------------------------------------------------
C          |   PROVA DI CADUTA SENZA PRE-SPIN DELLA RUOTA   |
C          --------------------------------------------------
           RAGGIO=0.
           DO 111 I=1,3
              RAGGIO=RAGGIO+(ZCG(I+3,3)-CC(I+3))**2
111        CONTINUE
           RAGGIO=DSQRT(RAGGIO)
C          WRITE(IOUT,777)RAGGIO,BETAP,BETA
777        FORMAT(' RAGGIO,BETAP,BETA=',3E12.5)
           OMER=BLENG/RAGGIO*BETAP*DSIN(BETA)
           AAA=OME/OMER
           TH=DTANH(DABS((OME-OMER)/70.))
           SR=1.-AAA
           SIG=-1.
           IF(OME.GT.OMER) SIG=1.
      ELSE IF(IPC.NE.1.OR.OMEZER.NE.0.)THEN
C          ----------------------------
C          |   TUTTI GLI ALTRI CASI   |
C          ----------------------------
           OMEGA=0.
           IF(SR.EQ.0..OR.OME.EQ.OMEZER) GOTO 31
cc           if(sr.eq.0.) goto 31
           BI=TANGE(SR,TABFAT,NPFAT)
           IF(OME.NE.0.0) THEN
            if(ome-omezer.ne.0.0)then
              OMEGA=OME-OMEP*AREA*AP*BI/(RK*(OME-OMEZER))
            else
              omega=0.0
            end if
           ELSE
C              OMEGA=0.
               omega=ome
           ENDIF
31         OMER=0.
           DO 30 I=1,3
             OMER=OMER+VEL(I)*RR(I)
30         CONTINUE
           OMER=OMER/A
           SR=(OMER-OMEGA)/(OMER-OMEZER)
      WRITE(6,'(a,5e10.4)')' VEL ',(VEL(IOL),iol=1,3)
      WRITE(6,'(5(a,e10.4,2x))')' OMER',OMER,'OMEGA',OMEGA
      WRITE(6,'(5(a,e10.4,2x))')' OMER',OMER,'OMEZER',OMEZER
      WRITE(6,'(2(a,e10.4,2x))')' SR  ',SR
2222      FORMAT(' SR',4E10.4)
           SIG=SIGN(1.,SR)
           TH=1.
      ENDIF
      SRO=SR
      SR=DABS(SR)
      IF(SR.GT.1.) SR=1.
      FRIC=TH*RINT(SR,TABFAT,NPFAT)
       WRITE(6,*)'TH',TH,'SR',SR,fric
      DO 40 I=1,3
	 F(I)=SIG*PP(4)*FRIC*AREA*AP*RR(I)
40    CONTINUE
      WRITE(6,*)' SIG   ',SIG,' PP4',PP(4)
      WRITE(6,*)' FRIC  ',FRIC
      WRITE(6,*)' AREA  ',AREA,' AP ',AP
      WRITE(6,*)' RR(I) ',(RR(IIO),IIO=1,3)
      WRITE(6,*)' TMMNT ',(RES(I),I=1,6)
      CALL TMMNT(ZCG(4,3),F,CC(4),RES)
      WRITE(6,*)' ZCG   ',(ZCG(3+I,3),I=1,3)
      WRITE(6,*)' F     ',(F(I),I=1,3),ipc,ipcad
      WRITE(6,*)' CC4   ',(CC(3+I),I=1,3)
      WRITE(6,*)' TMMNT ',(RES(I),I=1,6)
100   CONTINUE
       write(*,'(4(a,e12.5))') 'fric',fric,'th',th,'cat',
     *  fric/th,'sr',sr
       write(*,'(A,6e12.5)') ' F2 ',(RES(iol),iol=1,3)
     
C    ----------------------
C    |   DRIVING FORCE    |
C    ----------------------
      IF(VA.gt.GP(2)) then

      IF(A.GT.0.) GOTO 10
      DFCR=.7*WP(6)*PP(6)
      CD=DFCR*AP*AREA
      DO 5 I=1,3
5        F(I)=-CD*(VEP(1)*CDP(I,1)+VEP(2)*CDP(I,2))/VA
      GOTO 15
10    CALL PRVN(CDP(1,3),RR,RD,A)
      CALL PRVN(RD,CDP(1,3),RS,A)
      CALL TRSFV(CDP,RD,RR,-1)
      SINA=(VEP(1)*RR(1)+VEP(2)*RR(2))/VA
      COSA=DSQRT(1.-SINA**2)
      ANDER=DATAN2(SINA,COSA)
      BB=1.32*AREA/WP(5)+2.E-3
      BR=.5*(BB+1.75-DSQRT((BB-1.75)**2+1.E-2))
      BR=(2.-BR)*BR
      FD=BR/(BB+1.E-4)
      ANDER=FD*ANDER
      DFCR=WP(6)*PP(6)
      CD=DFCR*AP*AREA*DTANH(.7*ANDER/DFCR)
      DO 12 I=1,3
        F(I)=-CD*RD(I)
12    CONTINUE
      A=2.*CC(7)*WP(7)*DTANH(.7*DABS(ANDER))*COSA
      DO 14 I=1,3
        B(I)=(ZCG(I+3,3)-RS(I)*A*SIG)*FBET
14    CONTINUE
15    CALL TMMNT(B,F,CC(4),RES)
      else
      write(*,*) ' Calcoli Driving non effettuati VA.gt.GP(2)',va,gp(2)
      endif
      write(*,'(A,5e12.5)') ' F3 ',(RES(iol),iol=1,3)

C     -------------------------------
C     | CALCOLO MOMENTO SULLA RUOTA |
C     -------------------------------
      DO 16 I=1,3
	RM=RM+RES(I+3)*CC(I)
16    CONTINUE
      WRITE(6,*)'RETURN TYFRC' 
      RETURN
      END
C     COMPILER (ARGCHK=OFF)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-82                      SUBROUTINE FORCES
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE FORCES(CC,VEL,ZCG,VOL,CDP,PP,OP,RED,GP,RES,RP,AP,FBET,
     1                 IRCW)
      IMPLICIT REAL*8 (A-H,O-Z)
C       VEL      VELOCITY OF THE CENTRE OF CONTACT AREA
C       ZCG      ORIENTED AREA AND CENTRE OF GRAVITY
C       CDP      MATRIX OF THE DIRECTION COSINE OF THE EQUIVALENT PLANE
C       VOL      CONTACT VOLUME
C       PP       PLANE PARAMETERS
C       OP       DEFORMABLE OBSTACLE PARAMETERS
C       GP       GLOBAL PARAMETERS
C       RES 1-3  CONTACT REACTION FORCE
C       RES 4-6  CONTACT REACTION MOMENT
C       RP       MOMENT REDUCTION POINT
C       AP       ACTUAL PRESSURE
C       IRCW     RETURN CONDITION WORD
C            =1  DEFORMATION OF CILINDER ONLY
C            =2  DEFORMATION OF CILINDER AND PLANE
C            =3  DEFORMATION OF PLANE ONLY
C
      DIMENSION VEL(3),ZCG(6,3),CDP(3,4),PP(6),OP(4),GP(3),RES(6)
      DIMENSION VEP(3),F(6),CC(9)
      DATA PGRQ/.785398/
      SINB=0.
C      WRITE(6,*)'ENTRO IN FORCES'
      DO 1 I=1,3
        SINB=SINB+CC(I)*CDP(I,3)
1     CONTINUE
      SINB=DABS(SINB)
      COSB=DSQRT(1.-SINB**2)
      BET=DATAN2(SINB,COSB)
      FBET=0.
      IF(BET.LT.PGRQ)FBET=DCOS(2.*BET)
      OPR=OP(1)*(RED+(1.-RED)*FBET)
C    ---------------------------------
C    |        ACTUAL PRESSURE        |
C    ---------------------------------
2     GM1=1./GP(1)
      RAP=(OPR/PP(1))**GM1
      RAP=RAP*(OP(2)/PP(2))
      X=(VOL+RAP*PP(2)-OP(2))/(1.+RAP)
c      WRITE(6,'(4e10.4)')VOL,OPR
c      WRITE(6,8880)PP(1),PP(2),OP(1),OP(2),GP(1),VOL,X
c 8880 FORMAT(' PP1,PP2,OP1,OP2,GP,VOL,X'/7E12.5)
      IF(X.GT.0..AND.X.LT.VOL) GOTO 30
      IF(X.GT.0.) GOTO 20
10    IGO=2
      IRCW=1
      X=0.
      AP=OPR*(OP(2)/(OP(2)-VOL))**GP(1)

      print *,opr,op(2),vol,gp(1)
      
CC      WRITE(6,'(4e10.4)')VOL,AP
      GOTO 40
20    IGO=1
      IRCW=3
      X=VOL
      AP=PP(1)*(PP(2)/(PP(2)-VOL))**GP(1)
      GOTO 40
30    IGO=1
      IRCW=2
      AP=OPR*(OP(2)/(OP(2)-VOL+X))**GP(1)
40    CONTINUE
C    -------------------------------------------
C    |         ZERO-ING RESULTANT FORCE        |
C    -------------------------------------------
      CALL RSET(RES,6,0.D0)
      CALL RSET(F,6,0.D0)
C    ----------------------------------
C    |         ELASTIC FORCE          |
C    ----------------------------------
      DO 100 I=1,3
        F(I)=ZCG(I,3)*AP
100   CONTINUE
      CALL TMMNT(ZCG(4,3),F,RP,RES)
C     WRITE(6,8887)(F(IOO),IOO=1,3)
8887  FORMAT(' >>>F ELAS',3E12.5)
C    ------------------------------------------------------
C    |         VELOCITY IN EQUIVALENT PLANE FRAME         |
C    ------------------------------------------------------
C     WRITE(6,5010)VEL(1),VEL(3)
C5010 FORMAT(' VEL(1,3)=',3E15.5)
      CALL TRSFV(CDP,VEL,VEP,-1)
      VA=DSQRT(VEP(1)**2+VEP(2)**2)
      VM=DSQRT(VEP(1)**2+VEP(2)**2+VEP(3)**2)
C     WRITE(6,8886)VEP,VA,VM
C8886 FORMAT(' VEP,VA,VM',5E10.5)
C    -----------------------------------
C    |         DAMPING FORCE           |
C    -----------------------------------
      DC=DTANH(VEP(3)/OP(3))
      DO 200 I=1,3
        F(I)=-DC*ZCG(I,3)*AP
200   CONTINUE
C     WRITE(6,5887)(F(IOO),IOO=1,3)
5887  FORMAT(' >>>F DAMP',3E12.5)
      CALL TMMNT(ZCG(4,3),F,RP,RES)
C    -----------------------------------
C    |         FRICTION FORCE          |
C    -----------------------------------
      IF(VA.LT.GP(2)) GOTO 310
      AN=DSQRT(ZCG(1,3)**2+ZCG(2,3)**2+ZCG(3,3)**2)
      FC=OP(4)*PP(4)*AN*AP/VA
      DO 300 I=1,3
        F(I)=-FC*(VEP(1)*CDP(I,1)+VEP(2)*CDP(I,2))
300   CONTINUE
      CALL TMMNT(ZCG(4,3),F,RP,RES)
C     WRITE(6,8884)(F(IOO),IOO=1,3)
8884  FORMAT(' >>>F ATTR',3E12.5)
310   CONTINUE
C    ----------------------------------
C    |         PLOWING FORCE          |
C    ----------------------------------
      IF(IGO.EQ.2) GOTO 420
      IF(VA.LT.GP(2)) GOTO 420
      IF(VOL.LE.0.) GOTO 420
      PC=X*2.*AP*PP(5)*DTANH(VA/PP(3))/(VA*VOL)
      DO 410 J=1,2
        A1=0.
        DO 402 K=1,3
          A1=A1+ZCG(K,J)*VEL(K)
402     CONTINUE
        DO 406 I=1,3
          F(I)=-PC*A1*CDP(I,J)
406     CONTINUE
        CALL TMMNT(ZCG(4,J),F,RP,RES)
410   CONTINUE
420   CONTINUE
7777  RETURN
      END
C -------------------------------------------------------------------
C -------------------------------------------------------------------
C                 30- 1-82                         FUNCTION TANGE
C -------------------------------------------------------------------
C -------------------------------------------------------------------
      REAL*8 FUNCTION TANGE(X,XV,NXV)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XV(1)
      IF(NXV.GT.1)THEN
        DO 20 I=1,NXV
          IF(X.LT.XV(I)) GOTO 30
20      CONTINUE
        I=NXV
30      IF(I.EQ.1) I=I+1
        TANGE=(XV(I+NXV)-XV(I+NXV-1))/(XV(I)-XV(I-1))
      ELSE
        TANGE=0.
      END IF
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE VERCP
C ----------------------------------------------------------------------
C ---------------------------------------------------------------------
C@                 CALCOLA LA VELOCITA' DEL PUNTO DI CONTATTO CON IL
C@                 PIANO EQUIVALENTE DALLA VELOCITA' DEL CENTRO DEL
C@                 CILINDRO
C@
      SUBROUTINE VERCP(O,V,P,VC)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     ****
C     ****       VC=LV+AV^(P-0)
C     ****
C     LV         LINEAR  VELOCITY --- V(1),V(2),V(3)
C     AV         ANGULAR VELOCITY --- V(4),V(5),V(6)
C
C
C
      DIMENSION O(3),P(3),V(6),VC(3),PMO(3)
      DO 10 I=1,3
      PMO(I)=P(I)-O(I)
   10 CONTINUE
      CALL PRVU(V(4),PMO,VC)
      DO 20 I=1,3
      VC(I)=VC(I)+V(I)
   20 CONTINUE
      RETURN
      END

C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE PRVN
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C@                  CALCOLA IL PRODOTTO VETTORIALE DI DUE VETTORI E LO
C@                  NORMALIZZA
C@
      SUBROUTINE PRVN(X,Y,T,D)
      IMPLICIT REAL*8 (A-H,O-Z)
C     ****
C     ****       NORMALIZED VECTORIAL PRODUCT
C     ****       T=/X^Y/
C     ****
      DIMENSION X(3),Y(3),T(3),Z(3)
      CALL PRVU(X,Y,Z)
      D=Z(1)*Z(1)+Z(2)*Z(2)+Z(3)*Z(3)
      IF(D.GT.1.E-12) GOTO 1
      D=0.
      RETURN
    1 D=DSQRT(D)
      DO 10 I=1,3
      T(I)=Z(I)/D
   10 CONTINUE
      RETURN
      END
C     COMPILER (ARGCHK=OFF)
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE RMOVE
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE RMOVE(RVI,RVO,N)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RVI(1),RVO(1)
      DO 10 I=1,N
      RVO(I)=RVI(I)
   10 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      REAL FUNCTION RINT
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      REAL*8 FUNCTION RINT(X,XYT,NP)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XYT(1)
      IF(NP.GT.1)GOTO10
      RINT=XYT(2)
      RETURN
   10 CONTINUE
      DO 20 I=2,NP
      IF(X.LT.XYT(I))GOTO30
   20 CONTINUE
      I=NP
   30 IM1=I-1
      RINT=(XYT(NP+I)-XYT(NP+IM1))/(XYT(I)-XYT(IM1))*(X-XYT(IM1))
      RINT=XYT(NP+IM1)+RINT
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE RSET
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE RSET(RV,N,RVAL)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RV(1)
      DO 10 I=1,N
      RV(I)=RVAL
   10 CONTINUE
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE TMMNT
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE TMMNT(PA,FA,PB,FB)
C
C     ACCUMULA NELLE FORZE FB,APPLICATE NEL PUNTO PB,
C                 LE FORZE FA,APPLICATE NEL PUNTO PA
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION PA(3),PB(3),FA(6),FB(6)
      DIMENSION T(3)
      DO 10 I=1,3
      FB(I)=FB(I)+FA(I)
   10 T(I)=PA(I)-PB(I)
      CALL PRVU(T,FA,T)
      DO 20 I=1,3
   20 FB(I+3)=FB(I+3)+T(I)
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE TRSFV
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C@                  OPERA UNA TRASFORMAZIONE DI SISTEMI DI RIFERIMENTO
C@                  CARTESIANI ORTOGONALI SULLE COMPONENTI DI UN VETTORE
C@
      SUBROUTINE TRSFV(TM,RV,RN,IVER)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION TM(3,1),RV(1),RN(1)
      ITR=IABS((IVER-1)/2)
      CALL MV3(TM,RV,RN,ITR)
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE PRVU
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE PRVU(X,Y,T)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(3),Y(3),T(3)
      Z1=X(2)*Y(3)-X(3)*Y(2)
      Z2=X(3)*Y(1)-X(1)*Y(3)
      Z3=X(1)*Y(2)-X(2)*Y(1)
      T(1)=Z1
      T(2)=Z2
      T(3)=Z3
      RETURN
      END
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
C                   14- 1-83                      SUBROUTINE MV3
C ----------------------------------------------------------------------
C ----------------------------------------------------------------------
      SUBROUTINE MV3(RM,VIN,VOUT,IT)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION RM(1),VIN(1),VOUT(1)
      V1=VIN(1)
      V2=VIN(2)
      V3=VIN(3)
      IF(IT.GT.0)GOTO10
      VOUT(1)=RM(1)*V1+RM(4)*V2+RM(7)*V3
      VOUT(2)=RM(2)*V1+RM(5)*V2+RM(8)*V3
      VOUT(3)=RM(3)*V1+RM(6)*V2+RM(9)*V3
      RETURN
   10 VOUT(1)=RM(1)*V1+RM(2)*V2+RM(3)*V3
      VOUT(2)=RM(4)*V1+RM(5)*V2+RM(6)*V3
      VOUT(3)=RM(7)*V1+RM(8)*V2+RM(9)*V3
      RETURN
      END
