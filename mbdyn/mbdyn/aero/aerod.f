C $Header$
C*************************           F_AERO         ********************STE20560
C*************************      BLADEA MODIFIED     ********************STE20560



      SUBROUTINE F_AERO(IFLAG,TC,                            
     *                 nnodi,naero,nodaero,paraero,
     *                 nodid,gdl,pos,rotmat,vel,pointtype, umeanout,
     *                 alfaero,tempo)  
c
c   "assembla" nel residuo i termini relativi alle forze aerodinamiche 
c
c   INPUT: 
c
c   IFLAG      se =0 => no calcolo jacobiano aerodin
c   TC         trazione coppia
c   nnodi      numero totale nodi del modello 
c   naero      numero totale degli elementi aerodinamici di rotore
c   nodaero    1-3 nodi dei corpi cui è collegata la trave,
c              4 numero pti di integraz per ogni vol fin, 5 label dell'elemento
c   paraero    parametri descrittori dei tipi di elemento aerod:
c              1-9 offsets,10-12 corde nodali,13-15 sverg nodali,
c              16-18 dist c/4-var passo.19-21 dist 3c/4-var passo,
c              22-24 vettore nel sistema di rif di un corpo predef (rotore)
c              per calcolo posizione azimutale
c              25-33 matrice di rotazione sist sezione -> sistema locale
c              (si presuppone che sia la stessa per i tre nodi)
c   nodid      vettore con le labels numeriche dei nodi
c   gdl        vettore dei vettori delle posizioni dei gdl dei singoli nodi
c              all'interno di DYNSTIF
c   pos        posizioni dei nodi
c   rotmat     matrici di rotazione dei nodi 
c   vel        velocità dei nodi 
c   pointtype  tipo dei punti di equilibrio per vol fin:
c                0=pti medi  1=pti di gauss  
c   umeanout   velocità indotta media
c   alfaero    matrice contenente tempi e angoli incidenza ai nodi elemento
c   t          tempo


c      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4(I-N)
      implicit none

      integer*4 nnodi, naero,nodaero(5,*), nodid(*),
     *          gdl(9,*),pointtype,IFLAG
      real*8  pos(9,*),rotmat(9,*),vel(9,*), paraero(33,*), umeanout
      real*8  alfaero(2,3,3,1), tempo       

      integer*4 IRFI(2,12),ICFI(2,12),int                             
      DOUBLE PRECISION XQ(3,3),W_inerz(6),BLQ(3,3),BLQT(3,3),       
     *                 RQ(3,3),ALQ(3,3),W(6),TC(6),TC0(6),T,
     *                 VAM(6),VRHL(6),OUTA(20),vcstr(2)                 

      real*8 cospsi,sinpsi,psi,sins,coss,dist,dl,wght


      integer*4 IV(1)
      common IV
      real*8 DV(1)
      equivalence(IV,DV)

      integer*4 i, aero, nodi(3), SEARCH, 
     *          JPNT, jrp, jcp, jr,LM(18)
     
      

      COMMON/DATAERO/UMEAN,ro,cs,pesoU,rrot,AREA,inst,
     *                 imeth,jpro,mozzo,veliv
      real*8  UMEAN,ro,cs,pesoU,rrot,AREA
      integer*4  inst,imeth,jpro,mozzo,veliv

c
c     UMEAN      velocità indotta media 
c     ro         densità dell'aria
c     cs         celerità del suono
c     pesoU      peso per calcolo media vel indotta
c     rrot       raggio del rotore
c     AREA       area del rotore
c     INST       se=0 allora stazionario
c     IMETH      tipo di distribuzione vel indotta: 1=cost,2=Glauert,3=Mangler,4=din-inflow
c     jpro       tipo di profilo 
c     mozzo      numero della posizione del corpo mozzo nei vettori interni
c     veliv      numero della posizione del corpo velivolo nei vettori interni
c
      

      common/dinflow/TCvel,lamd,mud,umeand,roA,OM,Rtip
      real*8 TCvel(6),lamd,mud,umeand,roA,OM,Rtip

      COMMON/velind/MU,LAM,CHI,MODOMEGA
      real*8 MU,LAM,CHI,MODOMEGA

      real*8 residuo(18),braccio(3),mom(3)
      real*8 sa,sb,s(6),csi(6),ds_dcsi1,ds_dcsi2,ds_dcsi3,pesi(6)
      real*8 versV(3),versP(3),RT(3,3),Vmoz(3),Velnodi(3,3)   
      real*8 RVnodi(3,3),xnodi(3,3),of1(3),of2(3),of3(3),jaco(12,12)
      real*8 rot1(3,3),rot2(3,3),rot3(3,3),tng(12)
      real*8 MODVTRASL, dcosad, dsinad, umean0, umean1
c      real*8 Vmoz2(3)
	real*8 time(3),alf(3),coe(3),alf1(3),alf2(3),sver

      jrp = JPNT('ROWPOS  ')
      jcp = JPNT('COLPOS  ')
      jr  = JPNT('RESIDUO ')
      
C      open (102, file='azaerod.dat')
c      rewind 102
      
      call RESET_DP(residuo(1), 18, 1, 0.D0)

      call COPY_DP(TC, 6, 1, TC0)
      call RESET_DP(TC(1), 6, 1, 0.D0)
      umean0=umean

c  preparo dati mozzo per calcoli alfad e azimut
        call TRANSP(rotmat(1,mozzo),3,3,RT)
        CALL mult_0sab(1.0D0,RT,3,3,vel(1,mozzo),1,Vmoz)

*************************** modifica galleria del vento *************************
**********************************************************************
c        Vmoz2(1)=20.
c        CALL mult_0sab(1.0D0,RT,3,3,Vmoz2,1,Vmoz)
**********************************************************************
        if (abs(Vmoz(1))+abs(Vmoz(2)).GT.1D-9) then
          call VERSOR(Vmoz,versV)
        else
          versV(1)=0D0
          versV(2)=0D0
          versV(3)=0D0
        end if  
        versV(1)=-versV(1)
        versV(2)=-versV(2)
        versV(3)=0D0

c
c       calcolo alfad e modvtrasl
c
      MODVTRASL=DSQRT(Vmoz(1)**2+Vmoz(2)**2+Vmoz(3)**2)
      if (modvtrasl.LT.1D-9) then
        dcosad=0.D0
        dsinad=1.D0
      else  
        dcosad=-Vmoz(1)/MODVTRASL
        dsinad=Vmoz(3)/MODVTRASL
      end if  

c metto in VRHL la vel del mozzo
        call COPY_DP(pos(7,mozzo), 3, 1, VRHL(4))



c       calcolo trazione in sistema mozzo (con z concidente con z veliv) 
c       e TC in sistema velivolo
      call TRANSP(rotmat(1,veliv),3,3,RT)
      CALL mult_0sab(1.0D0,RT,3,3,TC0(1),1,TCvel(1))
      CALL mult_0sab(1.0D0,RT,3,3,TC0(4),1,TCvel(4))

c      TCvel(3)=TCvel(3)*4
      T=TCvel(3) 

c aggiorno il common/velind/
      IF(UMEAN.EQ.0.) UMEAN=.10D-3                                      
      CALL INDVEL(PSI,XQ,DIST, 5,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)
c calcolo velocità indotta media passo attuale
      if (abs(T).LT. 1D-9) then
        UMEAN1=0.
      else  
        UMEAN1=T/(2*ro*area*MODOMEGA*RROT*sqrt(mu**2+lam**2))
      end if


c      if (T.LT.0.) then
c        UMEAN=0.
c      else
c        a=2*ro*area
c      end if
c        UMEAN=(-a*modvtrasl*dcosad+dsqrt((a*modvtrasl*dcosad)**2+
c     *                                   4*a*T))/(2*a)


c   medio con peso velocità indotta media con quella passo precedente
      if (imeth.NE.4) UMEAN=(umean0*pesoU+umean1*(1-pesoU))
      IF(UMEAN.EQ.0.) UMEAN=.10D-3                                      
      umeanout= umean
      

c     write(*,*) 'Umean',umean
c     write(*,*) 'T=',TCvel(3)
c     write(*,*) TC0
      
c aggiorno il common/dinflow/
      if (imeth.EQ.4) then
        lamd=lam
        mud=mu
        umeand=umean
        OM=MODOMEGA
      end if
c preparo i parametri di integrazione per il tipo di punti dei volumi finiti
      if (pointtype.EQ.0) then
        sa=-0.5D0
        sb=-sa
c        s1=0.25*csi-0.75
c        s2=0.5*csi
c        s3=0.25*csi+0.75
        ds_dcsi1=.25D0
        ds_dcsi2=.5D0
        ds_dcsi3=.25D0
      else  
        sa=-.577350
        sb=-sa
c        s1=0.211324865*csi-0.78867513
c        s2=0.577350*csi
c        s3=0.211324865*csi+0.78867513
        ds_dcsi1=0.211324865D0
        ds_dcsi2=0.577350D0
        ds_dcsi3=0.211324865D0
      end if         
      

      VAM(1)=ro
      VAM(2)=cs


      
c inizio ciclo sugli elementi aerodinamici      
      do aero = 1,naero
        call RESET_DP(residuo(1), 18, 1, 0.D0)
        do i=1,3
          nodi(i) = SEARCH(nodaero(i,aero), nodid, nnodi)
        end do

        do i=1,6
          LM(i)=gdl(i,nodi(1))
          LM(6+i)=gdl(i,nodi(2))
          LM(12+i)=gdl(i,nodi(3))
        end do

c metto in VRHL la vel del mozzo
        call COPY_DP(vel(1,mozzo), 3, 1, VRHL(1))
        call COPY_DP(pos(7,mozzo), 3, 1, VRHL(4))
*************************** modifica galleria del vento *************************
**********************************************************************
c        VRHL(1)=20.
**********************************************************************


c  calcolo posizione dei nodi aerodinamici in xnodi(3,3)
c ***nodo1
        call COPY_DP(pos(1,nodi(1)), 3, 1,xnodi(1,1))
        CALL mult_0sab(1.D0,rotmat(1,nodi(1)),3,3,paraero(1,aero),1,of1)
        CALL madd(of1,xnodi(1,1),3,1,xnodi(1,1))
c ***nodo2
        call COPY_DP(pos(1,nodi(2)), 3, 1,xnodi(1,2))
        CALL mult_0sab(1.D0,rotmat(1,nodi(2)),3,3,paraero(4,aero),1,of2)
        CALL madd(of2,xnodi(1,2),3,1,xnodi(1,2))
c ***nodo3
        call COPY_DP(pos(1,nodi(3)), 3, 1,xnodi(1,3))
        CALL mult_0sab(1.D0,rotmat(1,nodi(3)),3,3,paraero(7,aero),1,of3)
        CALL madd(of3,xnodi(1,3),3,1,xnodi(1,3))

c  calcolo le matrici di rotazione sist sezione -> sistema globale
        CALL mult_0sab(1.D0,rotmat(1,nodi(1)),3,3,
     *                        paraero(25,aero),3,rot1)
        CALL mult_0sab(1.D0,rotmat(1,nodi(2)),3,3,
     *                        paraero(25,aero),3,rot2)
        CALL mult_0sab(1.D0,rotmat(1,nodi(3)),3,3,
     *                        paraero(25,aero),3,rot3)


c  metto in RVnodi i rotatinal vectors nodali
        call INV_ROT(rot1,RVnodi(1,1))
        call INV_ROT(rot2,RVnodi(1,2))
        call INV_ROT(rot3,RVnodi(1,3))

c  metto in Velnodi le velocità nodali
c ***nodo1
        call VECT(pos(7,nodi(1)),of1,Velnodi(1,1))
        CALL madd(vel(1,nodi(1)),Velnodi(1,1),3,1,Velnodi(1,1))
c ***nodo2
        call VECT(pos(7,nodi(2)),of2,Velnodi(1,2))
        CALL madd(vel(1,nodi(2)),Velnodi(1,2),3,1,Velnodi(1,2))
c ***nodo3
        call VECT(pos(7,nodi(3)),of3,Velnodi(1,3))
        CALL madd(vel(1,nodi(3)),Velnodi(1,3),3,1,Velnodi(1,3))



c       calcolo psi(rad) elemento (presuppone sistema mozzo con asse z perpend a suo piano)

        if (abs(versV(1))+abs(versV(2)) .GT. 0.01) then
          cospsi=versV(1)*paraero(22,aero)+versV(2)*paraero(23,aero)
          call VECT(versV,paraero(22,aero),versP)
          sinpsi=versP(3)
        
          PSI=dasin(sinpsi)
          if((sinpsi.GE.0).AND.(cospsi.LT.0)) PSI=3.14159265359-PSI
          if((sinpsi.LT.0).AND.(cospsi.LT.0)) PSI=-3.14159265359D0-PSI
          if(PSI.LT.0) PSI=PSI+6.2831853072D0
        else
          PSI=0D0
        end if
c	print *, 'Psi'        
c        Write(*,*) psi, dsinad  

c	calcolo angolo incidenza ai nodi
c	1' nodo
	sver=paraero(13,aero)
	sins=dsin(sver)
	coss=dcos(sver)
	call reset_dp(blq,3,3,0.d0)
	do i=1,3
		blq(i,1)=rot1(i,1)*coss+rot1(i,2)*sins
		blq(i,2)=-rot1(i,1)*sins+rot1(i,2)*coss
		blq(i,3)=rot1(i,3)
	enddo
        call MSUB(xnodi(1,1),pos(1,mozzo),3,1,XQ(1,1))
	call indvel(PSI,XQ(1,1),DIST,IMETH,dsinad,dcosad,MODVTRASL,
     1              VRHL(4),RROT,umean)
        do i=1,3
		W_inerz(i)=velnodi(i,1)+DISt*UMEAN*rotmat(6+i,mozzo)
	enddo
	call TRANSP(BLQ,3,3,BLQT)
	call mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
	call mult_0sab(1.0D0,BLQT,3,3,pos(7,nodi(1)),1,W(4))
	vcstr(1)=W(1)
	vcstr(2)=W(2)+W(6)*paraero(19,aero)
	time(1)=alfaero(1,1,1,aero)
	alf(1)=alfaero(2,1,1,aero)
	time(2)=alfaero(1,2,1,aero)
	alf(2)=alfaero(2,2,1,aero)
	time(3)=tempo
	alf(3)=datan2(-vcstr(2),vcstr(1))
	write( 14,*) 'nodo 1, aer',aero
	write(14,*) time, alf
	call polcoe(time,alf,3,coe)
	alf1(1)=coe(2)+2*coe(3)*tempo
	alf2(1)=2*coe(3)
c	
	alfaero(1,3,1,aero)=tempo
	alfaero(2,3,1,aero)=alf(3)

c	2' nodo
	sver=paraero(14,aero)
	sins=dsin(sver)
	coss=dcos(sver)
	call reset_dp(blq,3,3,0.d0)
	do i=1,3
		blq(i,1)=rot2(i,1)*coss+rot2(i,2)*sins
		blq(i,2)=-rot2(i,1)*sins+rot2(i,2)*coss
		blq(i,3)=rot2(i,3)
	enddo
        call MSUB(xnodi(1,2),pos(1,mozzo),3,1,XQ(1,1))
	call indvel(PSI,XQ(1,1),DIST,IMETH,dsinad,dcosad,MODVTRASL,
     1              VRHL(4),RROT,umean)
        do i=1,3
		W_inerz(i)=velnodi(i,2)+DISt*UMEAN*rotmat(6+i,mozzo)
	enddo
	call TRANSP(BLQ,3,3,BLQT)
	call mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
	call mult_0sab(1.0D0,BLQT,3,3,pos(7,nodi(2)),1,W(4))
	vcstr(1)=W(1)
	vcstr(2)=W(2)+W(6)*paraero(20,aero)
	time(1)=alfaero(1,1,2,aero)
	alf(1)=alfaero(2,1,2,aero)
	time(2)=alfaero(1,2,2,aero)
	alf(2)=alfaero(2,2,2,aero)
	time(3)=tempo
	alf(3)=datan2(-vcstr(2),vcstr(1))
	write(14,*) 'nodo 2'
	write(14,*) alf
	call polcoe(time,alf,3,coe)
	alf1(2)=coe(2)+2*coe(3)*tempo
	alf2(2)=2*coe(3)
c
	alfaero(1,3,2,aero)=tempo
	alfaero(2,3,2,aero)=alf(3)
	
c	3' nodo
	sver=paraero(15,aero)
	sins=dsin(sver)
	coss=dcos(sver)
	call reset_dp(blq,3,3,0.d0)
	do i=1,3
		blq(i,1)=rot3(i,1)*coss+rot3(i,2)*sins
		blq(i,2)=-rot3(i,1)*sins+rot3(i,2)*coss
		blq(i,3)=rot3(i,3)
	enddo
        call MSUB(xnodi(1,3),pos(1,mozzo),3,1,XQ(1,1))
	call indvel(PSI,XQ(1,1),DIST,IMETH,dsinad,dcosad,MODVTRASL,
     1              VRHL(4),RROT,umean)
        do i=1,3
		W_inerz(i)=velnodi(i,3)+DISt*UMEAN*rotmat(6+i,mozzo)
	enddo
	call TRANSP(BLQ,3,3,BLQT)
	call mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
	call mult_0sab(1.0D0,BLQT,3,3,pos(7,nodi(3)),1,W(4))
	vcstr(1)=W(1)
	vcstr(2)=W(2)+W(6)*paraero(21,aero)
	time(1)=alfaero(1,1,3,aero)
	alf(1)=alfaero(2,1,3,aero)
	time(2)=alfaero(1,2,3,aero)
	alf(2)=alfaero(2,2,3,aero)
	time(3)=tempo
	alf(3)=datan2(-vcstr(2),vcstr(1))
	write(14,*) ' nodo 3'
	write(14,*) alf
	call polcoe(time,alf,3,coe)
	alf1(3)=coe(2)+2*coe(3)*tempo
	alf2(3)=2*coe(3)
	write(14,*) 'alf1 & 2'
	write(14,*) alf1,alf2
c
	alfaero(1,3,3,aero)=tempo
	alfaero(2,3,3,aero)=alf(3)
	
	
c       INIZIALIZZO INTEGRAZIONI
        call DATI_INT(nodaero(4,aero),csi,pesi)

C                                                                       
C     INTEGRAZIONE VOLUME1                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.25*csi(i)-0.75
          else  
            s(i)=0.211324865*csi(i)-0.78867513
          end if         
        end do                                                          
c       inizio ciclo su pti integrazione
        DO 50 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)
c outa(9)
	 call INTERP (s(int),alf1(1),alf1(2),alf1(3),outa(9),1)
c outa(10)
	call INTERP(s(int),alf2(1),alf2(2),alf2(3),outa(10),1)

C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),pos(7,nodi(1)),pos(7,nodi(2)),
     *                pos(7,nodi(3)),W_inerz(4),3)

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************


c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)

c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))

C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       


          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if

          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)


          WGHT=PESI(INT)*dl*ds_dcsi1

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(1),tng(1),6,1,residuo(1))

          call MSUB(pos(1,nodi(1)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)

          call VECT(braccio,tng(1),mom)
          call MADD(residuo(4),mom,3,1,residuo(4))
          
          call MADD(TC(1),tng(7),6,1,TC(1))
          
 50     CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME1                                    
C                                                                       


C                                                                       
C     INTEGRAZIONE VOLUME2                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.5*csi(i)
          else  
            s(i)=0.577350*csi(i)
          end if         
        end do                                                          
c       inizio ciclo su pti integrazione
        DO 500 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)
c outa(9)
	 call INTERP (s(int),alf1(1),alf1(2),alf1(3),outa(9),1)
c outa(10)
	call INTERP(s(int),alf2(1),alf2(2),alf2(3),outa(10),1)

C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),pos(7,nodi(1)),pos(7,nodi(2)),
     *                pos(7,nodi(3)),W_inerz(4),3)

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************

c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)
c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))


C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       

          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if

          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)


          WGHT=PESI(INT)*dl*ds_dcsi2

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(7),tng(1),6,1,residuo(7))

          call MSUB(pos(1,nodi(2)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)

          call VECT(braccio,tng(1),mom)
          call MADD(residuo(10),mom,3,1,residuo(10))
          
          call MADD(TC(1),tng(7),6,1,TC(1))

          
 500    CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME2                                    
C                                                                       


C                                                                       
C     INTEGRAZIONE VOLUME3                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.25*csi(i)+0.75
          else  
            s(i)=0.211324865*csi(i)+0.78867513
          end if         
        end do                                                          


c       inizio ciclo su pti integrazione
        DO 5000 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)
c outa(9)
	 call INTERP (s(int),alf1(1),alf1(2),alf1(3),outa(9),1)
c outa(10)
	call INTERP(s(int),alf2(1),alf2(2),alf2(3),outa(10),1)

C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),pos(7,nodi(1)),pos(7,nodi(2)),
     *                pos(7,nodi(3)),W_inerz(4),3)

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************

c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)
c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))

C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       

          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if

          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)


          WGHT=PESI(INT)*dl*ds_dcsi3

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(13),tng(1),6,1,residuo(13))

          call MSUB(pos(1,nodi(3)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)
          
          call VECT(braccio,tng(1),mom)
          call MADD(residuo(16),mom,3,1,residuo(16))
          
          call MADD(TC(1),tng(7),6,1,TC(1))

          
5000    CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME3                                    
C                                                                       

      
      
c   assemblo in residuo
          
        CALL vsparsad(DV(jr),LM,residuo,18)
        
c        write(*,*)
c        write(*,*) 'Trazione',TC(3)
        
        
      end do
      
      RETURN                                                            
      END                                                               


































C*************************          ACC_AERO         ********************



      SUBROUTINE ACC_AERO(IFLAG,TC,                            
     *                 nnodi,naero,nodaero,paraero,
     *                 nodid,gdl,pos,rotmat,vel,pointtype)
       
c
c   "assembla" nel residuo i termini relativi alle forze aerodinamiche 
c
c   INPUT: 
c
c   IFLAG      se =0 => no calcolo jacobiano aerodin
c   TC         trazione coppia
c   nnodi      numero totale nodi del modello 
c   naero      numero totale degli elementi aerodinamici di rotore
c   nodaero    1-3 nodi dei corpi cui è collegata la trave,
c              4 numero pti di integraz per ogni vol fin, 5 label dell'elemento
c   paraero    parametri descrittori dei tipi di elemento aerod:
c              1-9 offsets,10-12 corde nodali,13-15 sverg nodali,
c              16-18 dist c/4-var passo.19-21 dist 3c/4-var passo,
c              22-24 vettore nel sistema di rif di un corpo predef (rotore)
c              per calcolo posizione azimutale
c              25-33 matrice di rotazione sist sezione -> sistema locale
c              (si presuppone che sia la stessa per i tre nodi)
c   nodid      vettore con le labels numeriche dei nodi
c   gdl        vettore dei vettori delle posizioni dei gdl dei singoli nodi
c              all'interno di DYNSTIF
c   pos        posizioni dei nodi
c   rotmat     matrici di rotazione dei nodi 
c   vel        velocità dei nodi 
c   pointtype  tipo dei punti di equilibrio per vol fin:
c                0=pti medi  1=pti di gauss  
c


c      IMPLICIT REAL*8 (A-H,O-Z), INTEGER*4(I-N)
      implicit none

      integer*4 nnodi, naero,nodaero(5,*), nodid(*),
     *          gdl(9,*),pointtype,IFLAG
      real*8  pos(9,*),rotmat(9,*),vel(9,*), paraero(33,*)
      




      integer*4 IRFI(2,12),ICFI(2,12),int                             
      DOUBLE PRECISION XQ(3,3),W_inerz(6),BLQ(3,3),BLQT(3,3),       
     *                 RQ(3,3),ALQ(3,3),W(6),TC(6),TC0(6),
     *                 VAM(6),VRHL(6),OUTA(20)                 

      real*8 cospsi,sinpsi,psi,sins,coss,dist,dl,wght


      integer*4 IV(1)
      common IV
      real*8 DV(1)
      equivalence(IV,DV)

      integer*4 i, aero, nodi(3), SEARCH, 
     *          JPNT, jrp, jcp, jr,LM(18)
     
      

      COMMON/DATAERO/UMEAN,ro,cs,pesoU,rrot,AREA,inst,
     *                 imeth,jpro,mozzo,veliv
      real*8  UMEAN,ro,cs,pesoU,rrot,AREA
      integer*4  inst,imeth,jpro,mozzo,veliv

c
c     UMEAN      velocità indotta media 
c     ro         densità dell'aria
c     cs         celerità del suono
c     pesoU      peso per calcolo media vel indotta
c     rrot       raggio del rotore
c     AREA       area del rotore
c     INST       se=0 allora stazionario
c     IMETH      tipo di distribuzione vel indotta: 1=cost,2=Glauert,3=Mangler
c     jpro       tipo di profilo 
c     mozzo      numero della posizione del corpo mozzo nei vettori interni
c
      


      common/dinflow/TCvel,lamd,mud,umeand,roA,OM,Rtip
      real*8 TCvel(6),lamd,mud,umeand,roA,OM,Rtip


      COMMON/velind/MU,LAM,CHI,MODOMEGA
      real*8 MU,LAM,CHI,MODOMEGA

      real*8 residuo(18),braccio(3),mom(3)
      real*8 sa,sb,s(6),csi(6),ds_dcsi1,ds_dcsi2,ds_dcsi3,pesi(6)
      real*8 versV(3),versP(3),RT(3,3),Vmoz(3),Velnodi(3,3)   
      real*8 RVnodi(3,3),xnodi(3,3),of1(3),of2(3),of3(3),jaco(12,12)
      real*8 rot1(3,3),rot2(3,3),rot3(3,3),tng(12)
      real*8 MODVTRASL, dcosad, dsinad
c      real*8 Vmoz2(3)


      jrp = JPNT('ROWPOS  ')
      jcp = JPNT('COLPOS  ')
      jr  = JPNT('RESIDUO ')
      
      
      call RESET_DP(residuo(1), 18, 1, 0.D0)

      call COPY_DP(TC, 6, 1, TC0)
      call RESET_DP(TC(1), 6, 1, 0.D0)
      if (imeth .EQ. 4) umean=DV(JPNT('genvar0 '))
          
      
c  preparo dati mozzo per calcoli alfad e azimut
        call TRANSP(rotmat(1,mozzo),3,3,RT)
        CALL mult_0sab(1.0D0,RT,3,3,vel(1,mozzo),1,Vmoz)
*************************** modifica galleria del vento *************************
**********************************************************************
c        Vmoz2(1)=20.
c        CALL mult_0sab(1.0D0,RT,3,3,Vmoz2,1,Vmoz)
**********************************************************************
        if (abs(Vmoz(1))+abs(Vmoz(2)).GT.1D-9) then
          call VERSOR(Vmoz,versV)
        else
          versV(1)=0D0
          versV(2)=0D0
          versV(3)=0D0
        end if  
        versV(1)=-versV(1)
        versV(2)=-versV(2)
        versV(3)=0D0

c
c       calcolo alfad e modvtrasl
c
      MODVTRASL=DSQRT(Vmoz(1)**2+Vmoz(2)**2+Vmoz(3)**2)
      if (modvtrasl.LT.1D-9) then
        dcosad=0.D0
        dsinad=1.D0
      else  
        dcosad=-Vmoz(1)/MODVTRASL
        dsinad=Vmoz(3)/MODVTRASL
      end if  

c metto in VRHL la vel del mozzo
        call COPY_DP(vel(4,mozzo), 3, 1, VRHL(4))


c aggiorno il common/velind/
      IF(UMEAN.EQ.0.) UMEAN=1.D-3                                      
      CALL INDVEL(PSI,XQ,DIST, 5,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)




c      if (T.LT.0.) then
c        UMEAN=0.
c      else
c        a=2*ro*area
c      end if
c        UMEAN=(-a*modvtrasl*dcosad+dsqrt((a*modvtrasl*dcosad)**2+
c     *                                   4*a*T))/(2*a)


c   medio con peso velocità indotta media con quella passo precedente
      IF(UMEAN.EQ.0.) UMEAN=.10D-3                                      

c      write(*,*) 'Umean',umean


c aggiorno il common/dinflow/
      if (imeth.EQ.4) then
        lamd=lam
        mud=mu
        umeand=umean
        roA=ro*area
        OM=MODOMEGA
        Rtip=rrot
      end if


c preparo i parametri di integrazione per il tipo di punti dei volumi finiti
      if (pointtype.EQ.0) then
        sa=-0.5D0
        sb=-sa
c        s1=0.25*csi-0.75
c        s2=0.5*csi
c        s3=0.25*csi+0.75
        ds_dcsi1=.25D0
        ds_dcsi2=.5D0
        ds_dcsi3=.25D0
      else  
        sa=-.577350
        sb=-sa
c        s1=0.211324865*csi-0.78867513
c        s2=0.577350*csi
c        s3=0.211324865*csi+0.78867513
        ds_dcsi1=0.211324865D0
        ds_dcsi2=0.577350D0
        ds_dcsi3=0.211324865D0
      end if         
      

      VAM(1)=ro
      VAM(2)=cs


      
c inizio ciclo sugli elementi aerodinamici      
      do aero = 1,naero
        call RESET_DP(residuo(1), 18, 1, 0.D0)
        do i=1,3
          nodi(i) = SEARCH(nodaero(i,aero), nodid, nnodi)
        end do

        do i=1,6
          LM(i)=gdl(i,nodi(1))
          LM(6+i)=gdl(i,nodi(2))
          LM(12+i)=gdl(i,nodi(3))
        end do

c metto in VRHL la vel del mozzo
        call COPY_DP(vel(1,mozzo), 3, 1, VRHL(1))
        call COPY_DP(vel(4,mozzo), 3, 1, VRHL(4))
*************************** modifica galleria del vento *************************
**********************************************************************
c        VRHL(1)=20.
**********************************************************************


c  calcolo posizione dei nodi aerodinamici in xnodi(3,3)
c ***nodo1
        call COPY_DP(pos(1,nodi(1)), 3, 1,xnodi(1,1))
        CALL mult_0sab(1.D0,rotmat(1,nodi(1)),3,3,paraero(1,aero),1,of1)
        CALL madd(of1,xnodi(1,1),3,1,xnodi(1,1))
c ***nodo2
        call COPY_DP(pos(1,nodi(2)), 3, 1,xnodi(1,2))
        CALL mult_0sab(1.D0,rotmat(1,nodi(2)),3,3,paraero(4,aero),1,of2)
        CALL madd(of2,xnodi(1,2),3,1,xnodi(1,2))
c ***nodo3
        call COPY_DP(pos(1,nodi(3)), 3, 1,xnodi(1,3))
        CALL mult_0sab(1.D0,rotmat(1,nodi(3)),3,3,paraero(7,aero),1,of3)
        CALL madd(of3,xnodi(1,3),3,1,xnodi(1,3))

c  calcolo le matrici di rotazione sist sezione -> sistema globale
        CALL mult_0sab(1.D0,rotmat(1,nodi(1)),3,3,
     *                        paraero(25,aero),3,rot1)
        CALL mult_0sab(1.D0,rotmat(1,nodi(2)),3,3,
     *                        paraero(25,aero),3,rot2)
        CALL mult_0sab(1.D0,rotmat(1,nodi(3)),3,3,
     *                        paraero(25,aero),3,rot3)


c  metto in RVnodi i rotatinal vectors nodali
        call INV_ROT(rot1,RVnodi(1,1))
        call INV_ROT(rot2,RVnodi(1,2))
        call INV_ROT(rot3,RVnodi(1,3))

c  metto in Velnodi le velocità nodali
c ***nodo1
        call VECT(vel(4,nodi(1)),of1,Velnodi(1,1))
        CALL madd(vel(1,nodi(1)),Velnodi(1,1),3,1,Velnodi(1,1))
c ***nodo2
        call VECT(vel(4,nodi(2)),of2,Velnodi(1,2))
        CALL madd(vel(1,nodi(2)),Velnodi(1,2),3,1,Velnodi(1,2))
c ***nodo3
        call VECT(vel(4,nodi(3)),of3,Velnodi(1,3))
        CALL madd(vel(1,nodi(3)),Velnodi(1,3),3,1,Velnodi(1,3))



c       calcolo psi(rad) elemento (presuppone sistema mozzo con asse z perpend a suo piano)

        if (abs(versV(1))+abs(versV(2)) .GT. 0.01) then
          cospsi=versV(1)*paraero(22,aero)+versV(2)*paraero(23,aero)
          call VECT(versV,paraero(22,aero),versP)
          sinpsi=versP(3)
        
          PSI=dasin(sinpsi)
          if((sinpsi.GE.0).AND.(cospsi.LT.0)) PSI=3.14159265359-PSI
          if((sinpsi.LT.0).AND.(cospsi.LT.0)) PSI=-3.14159265359D0-PSI
          if(PSI.LT.0) PSI=PSI+6.2831853072D0
        else
          PSI=0D0
        end if  

c       INIZIALIZZO INTEGRAZIONI
        call DATI_INT(nodaero(4,aero),csi,pesi)

C                                                                       
C     INTEGRAZIONE VOLUME1                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.25*csi(i)-0.75
          else  
            s(i)=0.211324865*csi(i)-0.78867513
          end if         
        end do                                                          
c       inizio ciclo su pti integrazione
        DO 50 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)


C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),vel(4,nodi(1)),vel(4,nodi(2)),
     *                vel(4,nodi(3)),W_inerz(4),3)

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************

c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)

c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))

C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       
          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if
          
          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)

          WGHT=PESI(INT)*dl*ds_dcsi1

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(1),tng(1),6,1,residuo(1))

          call MSUB(pos(1,nodi(1)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)
          call VECT(braccio,tng(1),mom)
          call MADD(residuo(4),mom,3,1,residuo(4))
          call MADD(TC(1),tng(7),6,1,TC(1))
          
 50     CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME1                                    
C                                                                       


C                                                                       
C     INTEGRAZIONE VOLUME2                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.5*csi(i)
          else  
            s(i)=0.577350*csi(i)
          end if         
        end do                                                          
c       inizio ciclo su pti integrazione
        DO 500 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)


C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),vel(4,nodi(1)),vel(4,nodi(2)),
     *                vel(4,nodi(3)),W_inerz(4),3)

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************


c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)
c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))


C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       

          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if

          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)

          WGHT=PESI(INT)*dl*ds_dcsi2

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(7),tng(1),6,1,residuo(7))
          call MSUB(pos(1,nodi(2)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)
          call VECT(braccio,tng(1),mom)
          call MADD(residuo(10),mom,3,1,residuo(10))
          call MADD(TC(1),tng(7),6,1,TC(1))

          
 500    CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME2                                    
C                                                                       


C                                                                       
C     INTEGRAZIONE VOLUME3                                   
C       
c       calcolo punti int espressi secondo coord curvilinea s
        do i=1,nodaero(4,aero)
          if (pointtype.EQ.0) then
            s(i)=0.25*csi(i)+0.75
          else  
            s(i)=0.211324865*csi(i)+0.78867513
          end if         
        end do                                                          


c       inizio ciclo su pti integrazione
        DO 5000 INT=1,nodaero(4,aero)
      
c         creo il vettore VAM
c  corda
          call INTERP(s(int),paraero(10,aero),paraero(11,aero),
     *                paraero(12,aero),VAM(3),1)
c  B=c/4
          call INTERP(s(int),paraero(16,aero),paraero(17,aero),
     *                paraero(18,aero),VAM(4),1)
c  D=3c/4
          call INTERP(s(int),paraero(19,aero),paraero(20,aero),
     *                paraero(21,aero),VAM(5),1)
c  sver
          call INTERP(s(int),paraero(13,aero),paraero(14,aero),
     *                paraero(15,aero),VAM(6),1)


C
C         calcolo   XQ RQ W_inerz BLQ          
C                                                                       
          call INTERP(s(int),xnodi(1,1),xnodi(1,2),
     *                xnodi(1,3),XQ(1,1),3)
          call MSUB(XQ(1,1),pos(1,mozzo),3,1,XQ(1,1))
          
          
          call INTERP(s(int),RVnodi(1,1),RVnodi(1,2),
     *                RVnodi(1,3),RQ(1,1),3)


c   ALQ : rotazione da sist sezione a sist globale
          call INTERP(s(int),rot1,rot2,rot3,ALQ(1,1),9)


c   BLQ : rotazione da aerodinamico a globale

          SINS=DSIN(VAM(6))                                                   
          COSS=DCOS(VAM(6))                                                   
          DO I=1,3                                                      
            BLQ(I,1)=+ALQ(I,1)*COSS+ALQ(I,2)*SINS                            
            BLQ(I,2)=-ALQ(I,1)*SINS+ALQ(I,2)*COSS                             
            BLQ(I,3)=+ALQ(I,3)                                                
          END DO                                                      

*************************** modifica galleria del vento *************************
**********************************************************************
c        W_inerz(1)=W_inerz(1)+20.
**********************************************************************


c calcolo vel assoluta pto nel sist inerziale
          call INTERP(s(int),Velnodi(1,1),Velnodi(1,2),
     *                Velnodi(1,3),W_inerz(1),3)
          call INTERP(s(int),vel(4,nodi(1)),vel(4,nodi(2)),
     *                vel(4,nodi(3)),W_inerz(4),3)


c calcolo velind in sistema mozzo 
          CALL INDVEL(PSI,XQ,DIST,IMETH,dsinad,dcosad,MODVTRASL,
     *                VRHL(4),RROT,umean)
c la sommo a W_inerz sempre nel sist inerziale per avere vel relativa all'aria
          do i=1,3
            W_inerz(i)=W_inerz(i)+DIST*UMEAN*rotmat(6+i,mozzo)
          end do            

c calcolo vel relativa all'aria del pto nel sist aerodinamico: W
          call TRANSP(BLQ,3,3,BLQT)
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(1),1,W(1))
          CALL mult_0sab(1.0D0,BLQT,3,3,W_inerz(4),1,W(4))

C                                                                       
C     CALCOLO FORZE AERODINAMICHE
C                                                                       
C                                                                       

          if ((W(1)**2+W(2)**2+W(3)**2).gt.0.1) then
            CALL AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,jaco,12,         
     *                tng,IFLAG,IRFI,ICFI,OUTA,INST,MODOMEGA,1.D0)         
          else
            call reset_dp(tng,12,1,0.D0)
          end if

          call CALC_DL(s(int),xnodi(1,1),xnodi(1,2),xnodi(1,3),dl)

          WGHT=PESI(INT)*dl*ds_dcsi3

          do i=1,12
            tng(i)=tng(i)*WGHT
          end do  

          call MADD(residuo(13),tng(1),6,1,residuo(13))
          call MSUB(pos(1,nodi(3)),pos(1,mozzo),3,1,braccio)
          call MSUB(XQ(1,1),braccio,3,1,braccio)
          call VECT(braccio,tng(1),mom)
          call MADD(residuo(16),mom,3,1,residuo(16))
          call MADD(TC(1),tng(7),6,1,TC(1))

          
5000    CONTINUE                                                          
C                                                                       
C     FINE LOOP PUNTI INTEGRAZIONE  VOLUME3                                    
C                                                                       

      
      
c   assemblo in residuo
          
        CALL vsparsad(DV(jr),LM,residuo,18)
        
c        write(*,*)
c        write(*,*) TC
        
      end do
      
      RETURN                                                            
      END                                                               




















































C*************************          AEROD          *********************STE17770
c                        MODIFICATA NEL CALCOLO TNG      
C=    COMPILER (LINK=IBJ$)                                              STE17780
      SUBROUTINE AEROD (XQ,RQ,BLQ,DIST,PSI,W,VRHL,VAM,G,NRDG,           STE17790
     *                 TNG,IFLAG,IRFI,ICFI,OUTA,INST,RSPEED,DA)         STE17800
      IMPLICIT REAL*8(A-H,O-Z)                                          STE17810
      DOUBLE PRECISION XQ(9),RQ(9),BLQ(3,3),W(6),VRHL(6),               STE17820
     *                 VAM(*),G(NRDG,*),TNG(*),OUTA(20)                 STE17830
      DIMENSION IRFI(2,1),ICFI(2,1)                                     STE17840
      DOUBLE PRECISION DM(3,6),VCSTR(6),BM(6,4),AASTR(6,6),A1STR(4,3),  STE17850
     *                 AP(4,6),QASTR(6),DSTR(6,6),D1STR(4,3),B0(3),     STE17860
     *                  OGV(3,3),VOV(3,3),QAV(3,3),XQV(3,3),            STE17870
     *                 GAMFIR(3,3),GAMDOT(3,3)                          STE17880
C                                                                       STE17890
C     DEFINIZIONI VETTORE VAM                                           STE17900
C                                                                       STE17910
      DENS =VAM(1)                                                      STE17920
      CS   =VAM(2)                                                      STE17930
      CORDA=VAM(3)                                                      STE17940
      B    =VAM(4)                                                      STE17950
      D    =VAM(5)                                                      STE17960
C     SVER =VAM(6)                                                      STE17970
C                                                                       STE17980
C                                                                       STE17990
C     COSTRUZIONE VC*                                                   STE18000
C                                                                       STE18010
      CALL DZERO(DM,18)                                                 STE18020
      DO 10 I=1,3                                                       STE18030
      DM(I,I)=1.D0                                                      STE18040
 10   CONTINUE                                                          STE18050
      DM(2,6)=D                                                         STE18060
      CALL DPROMV(DM,3,3,6,W,VCSTR,0)                                   STE18070
      VP=DSQRT(VCSTR(1)**2+VCSTR(2)**2)                                 STE18080
      V=DSQRT(VP*VP+VCSTR(3)**2)                                        STE18090
C                                                                       STE18100
C     CALCOLO COEFFICIENTI AERODINAMICI                                 STE18110
C                                                                       STE18120
      IGO=INST+1                                                        STE18130
      GOTO(100,200,300),IGO                                             STE18140
 100  CALL COE0(VCSTR,OUTA,CS,CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,      STE18150
     &          CRF,CFSLOP,DCPDM,DCDRDM,DCMDM,DCRFDM)                   STE18160
      GO TO 400                                                         STE18170
 200  CALL COE1(VCSTR,OUTA,CS,CORDA,RSPEED,DA,                          STE18180
     &          CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,                    STE18190
     &          CRF,CFSLOP,DCPDM,DCDRDM,DCMDM,DCRFDM)                   STE18200
      GO TO 400                                                         STE18210
 300  CALL COE2(VCSTR,OUTA,CS,CORDA,RSPEED,DA,                          STE18220
     &          CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,                    STE18230
     &          CRF,CFSLOP,DCPDM,DCDRDM,DCMDM,DCRFDM)                   STE18240
C                                                                       STE18250
C     COSTRUZIONE A1*                                                   STE18260
C                                                                       STE18270
 400  CONTINUE                                                          STE18280
	write(14,*) ' in aerod'
	write(14,*) outa
      CALL DZERO(A1STR,12)                                              STE18290
      RK=-.5*DENS*CORDA                                                 STE18300
      A1STR(1,1)= RK*(CRF*(VP-V)-CDRAG*VP)                              STE18310
      A1STR(1,2)=-RK*CLIFT*VP                                           STE18320
      A1STR(2,1)=-A1STR(1,2)                                            STE18330
      A1STR(2,2)= A1STR(1,1)                                            STE18340
      A1STR(3,3)=-RK*CRF*V                                              STE18350
      A1STR(4,1)=-RK*CMOME*CORDA*VCSTR(1)                               STE18360
      A1STR(4,2)=-RK*CMOME*CORDA*VCSTR(2)                               STE18370
C     CALL DWRITE(' A1STR',A1STR,4,4,3,0)                               STE18380
C                                                                       STE18390
C     COSTRUZIONE AA*                                                   STE18400
C                                                                       STE18410
      CALL DZERO(BM,24)                                                 STE18420
      DO 20 I=1,3                                                       STE18430
      BM(I,I)=1.D0                                                      STE18440
 20   CONTINUE                                                          STE18450
      BM(6,2)=B                                                         STE18460
      BM(6,4)=1.D0                                                      STE18470
      CALL DPROMM(A1STR,4,4,3,DM,3,AP,4,6,0,0)                          STE18480
      CALL DPROMM(BM,6,6,4,AP,4,AASTR,6,6,0,0)                          STE18490
C                                                                       STE18500
C     COSTRUZIONE Q*                                                    STE18510
C                                                                       STE18520
      CALL MOVE(VCSTR(4),W(4),3)                                        STE18530
      CALL DPROMV(AASTR,6,6,6,VCSTR,QASTR,0)
c     cambio segno alle forze perchè le calcola con segno opposto                                  
      do i=1,6
        QASTR(i)=-QASTR(i)
      end do  


c      paz= (XQ(1)**2+XQ(2)**2+XQ(3)**2)**.5
c      write(102,'(7(1Pe13.5))') paz,(QASTR(i),i=1,6)



C                                                                       STE18550
C     COSTRUZIONE TERMINE NOTO TNG                                      STE18560
C                                                                       STE18570
      CALL DZERO(TNG,12)                                                STE18580
      CALL DPROMV(BLQ,3,3,3,QASTR(1),TNG(1),0)                          STE18590
      CALL DPROMV(BLQ,3,3,3,QASTR(4),TNG(4),0)                          STE18600
C                                                                       STE18610
C                                                                       STE18620
      CALL DRMOVE(TNG(7),TNG(1),6,1.D0)                                 STE18630
      CALL DPRVAD(XQ(1),TNG(1),TNG(10))                                 STE18640
C                                                                       STE18650
C                                                                       STE18660
********************************************************************************
*    le prossime due righe sono commentate per avere TNG nel sist assoluto
*    dando la opportuna BLQ    (sett '95)
********************************************************************************
c      CALL PSIROT(PSI,TNG( 7),NRDG,1,TNG( 7),NRDG)                      STE18670
c      CALL PSIROT(PSI,TNG(10),NRDG,1,TNG(10),NRDG)                      STE18680
********************************************************************************

      DO 30 I=1,6                                                       STE18690
      IRFI(1,I)=1                                                       STE18700
      IRFI(2,I)=I                                                       STE18710
 30   CONTINUE                                                          STE18720
      IF(IFLAG.EQ.0) RETURN                                             STE18730
C     *********************                                             STE18740
C                                                                       STE18750
C     COSTRUZIONE D1*                                                   STE18760
C                                                                       STE18770
      CALL DZERO(A1STR,12)                                              STE18780
      A1STR(1,1)=-RK*((BSLOP-CFSLOP)*VP+CFSLOP*V)                       STE18790
      A1STR(1,2)=-RK*ASLOP*VP                                           STE18800
      A1STR(2,1)=-A1STR(1,2)                                            STE18810
      A1STR(2,2)= A1STR(1,1)                                            STE18820
      A1STR(3,3)=-RK*CFSLOP*V                                           STE18830
      A1STR(4,1)= RK*CSLOP*VCSTR(1)                                     STE18840
      A1STR(4,2)= RK*CSLOP*VCSTR(2)                                     STE18850
C     CALL DWRITE(' A1STR',A1STR,4,4,3,0)                               STE18860
      B0(1)=VCSTR(2)/(VP*VP)                                            STE18870
      B0(2)=-VCSTR(1)/(VP*VP)                                           STE18880
      B0(3)=0.D0                                                        STE18890
      CALL DPROMV(A1STR,4,4,3,VCSTR,AP,0)                               STE18900
      CALL DZERO(D1STR,12)                                              STE18910
      CALL DPRVVA(AP,B0,D1STR,4,4,3)                                    STE18920
C                                                                       STE18930
      A1STR(1,1)=-RK*(DCRFDM*V+(DCDRDM-DCRFDM)*VP)                      STE18940
      A1STR(1,2)=-RK*DCPDM*VP                                           STE18950
      A1STR(2,1)=-A1STR(1,2)                                            STE18960
      A1STR(2,2)= A1STR(1,1)                                            STE18970
      A1STR(3,3)=-RK*DCRFDM*V                                           STE18980
      A1STR(4,1)= RK*DCMDM*VCSTR(1)                                     STE18990
      A1STR(4,2)= RK*DCMDM*VCSTR(2)                                     STE19000
C     CALL DWRITE(' A1STR',A1STR,4,4,3,0)                               STE19010
      B0(1)=VCSTR(1)/(VP*CS)                                            STE19020
      B0(2)=VCSTR(2)/(VP*CS)                                            STE19030
      B0(3)=0.                                                          STE19040
      CALL DPROMV(A1STR,4,4,3,VCSTR,AP,0)                               STE19050
      CALL DPRVVA(AP,B0,D1STR,4,4,3)                                    STE19060
C                                                                       STE19070
      A1STR(1,1)=-RK*CRF                                                STE19080
      A1STR(1,2)=0.D0                                                   STE19090
      A1STR(2,1)=0.D0                                                   STE19100
      A1STR(2,2)= A1STR(1,1)                                            STE19110
      A1STR(3,3)=-RK*CRF                                                STE19120
      A1STR(4,1)=0.                                                     STE19130
      A1STR(4,2)=0.                                                     STE19140
C     CALL DWRITE(' A1STR',A1STR,4,4,3,0)                               STE19150
      B0(1)=VCSTR(1)/V                                                  STE19160
      B0(2)=VCSTR(2)/V                                                  STE19170
      B0(3)=VCSTR(3)/V                                                  STE19180
      CALL DPROMV(A1STR,4,4,3,VCSTR,AP,0)                               STE19190
      CALL DPRVVA(AP,B0,D1STR,4,4,3)                                    STE19200
C                                                                       STE19210
      A1STR(1,1)=-RK*(CDRAG-CRF)                                        STE19220
      A1STR(1,2)=-RK*CLIFT                                              STE19230
      A1STR(2,1)=-A1STR(1,2)                                            STE19240
      A1STR(2,2)= A1STR(1,1)                                            STE19250
      A1STR(3,3)= 0.D0                                                  STE19260
C     CALL DWRITE(' A1STR',A1STR,4,4,3,0)                               STE19270
      B0(1)=VCSTR(1)/VP                                                 STE19280
      B0(2)=VCSTR(2)/VP                                                 STE19290
      B0(3)=0.D0                                                        STE19300
      CALL DPROMV(A1STR,4,4,3,VCSTR,AP,0)                               STE19310
      CALL DPRVVA(AP,B0,D1STR,4,4,3)                                    STE19320
C                                                                       STE19330
      D1STR(4,1)=D1STR(4,1)+RK*CMOME*CORDA*VCSTR(1)                     STE19340
      D1STR(4,2)=D1STR(4,2)+RK*CMOME*CORDA*VCSTR(2)                     STE19350
C                                                                       STE19360
C     COSTRUZIONE D*                                                    STE19370
C                                                                       STE19380
      CALL DPROMM(D1STR,4,4,3,DM,3,AP,4,6,0,0)                          STE19390
      CALL DPROMM(BM,6,6,4,AP,4,DSTR,6,6,0,0)                           STE19400
C     CALL DWRITE(' AASTR',AASTR,6,6,6,0)                               STE19410
C     CALL DWRITE('  DSTR',DSTR ,6,6,6,0)                               STE19420
C     CALL DWRITE(' QASTR',QASTR,1,6,1,-1)                              STE19430
C     PRINT '(A)',' ***************************************'            STE19440
C                                                                       STE19450
C     COSTRUZIONE G*                                                    STE19460
C                                                                       STE19470
      DO 35 I=1,12                                                      STE19480
      DO 35 K=1,19                                                      STE19490
      G(I,K)=0.D0                                                       STE19500
 35   CONTINUE                                                          STE19510
      CALL ADD(G,NRDG,AASTR,6,6,6,1.D0)                                 STE19520
      CALL ADD(G,NRDG,DSTR,6,6,6,1.D0)                                  STE19530
      CALL DVET(W,VOV,1)                                                STE19540
      CALL DPROMM(G,NRDG,6,3,VOV,3,G(1,10),NRDG,3,0,0)                  STE19550
C     CALL DVET(W(4),VOV,1)                                             STE19560
C     CALL DPROMM(G(1,4),NRDG,6,3,VOV,3,G(1,10),NRDG,3,0,0)             STE19570
      CALL DVET(QASTR(1),QAV,1)                                         STE19580
      CALL ADD(G(1,10),NRDG,QAV,3,3,3,-1.D0)                            STE19590
      CALL DVET(QASTR(4),QAV,1)                                         STE19600
      CALL ADD(G(4,10),NRDG,QAV,3,3,3,-1.D0)                            STE19610
C                                                                       STE19620
C     COSTRUZIONE G                                                     STE19630
C                                                                       STE19640
      CALL DABCT3(BLQ,G( 1, 1),NRDG,BLQ)                                STE19650
      CALL DABCT3(BLQ,G( 1, 4),NRDG,BLQ)                                STE19660
      CALL DABCT3(BLQ,G( 1,10),NRDG,BLQ)                                STE19670
      CALL DABCT3(BLQ,G( 4, 1),NRDG,BLQ)                                STE19680
      CALL DABCT3(BLQ,G( 4, 4),NRDG,BLQ)                                STE19690
      CALL DABCT3(BLQ,G( 4,10),NRDG,BLQ)                                STE19700
C                                                                       STE19710
      CALL DVET(XQ(1),XQV,1)                                            STE19720
      CALL DVET(VRHL(4),OGV,1)                                          STE19730
C                                                                       STE19740
      CALL DPRMMA(G(1,1),NRDG,6,3,OGV,3,G(1,7),NRDG,3,0,0)              STE19750
      CALL DPRMMA(G(1,4),NRDG,6,3,OGV,3,G(1,10),NRDG,3,0,0)             STE19760
C                                                                       STE19770
      CALL ADD(G(1,13),NRDG,G(1,1),NRDG,6,6,1.D0)                       STE19780
      CALL DPRMMA(G(1,1),NRDG,6,3,XQV,3,G(1,16),NRDG,3,0,1)             STE19790
C                                                                       STE19800
C                                                                       STE19810
      CALL ADD(G(7,1),NRDG,G(1,1),NRDG,6,18,1.D0)                       STE19820
      CALL DPRMMA(XQV,3,3,3,G(1,1),NRDG,G(10,1),NRDG,18,0,0)            STE19830
      DO 40 I=1,6                                                       STE19840
      ICFI(1,I)=3                                                       STE19850
      ICFI(2,I)=I                                                       STE19860
      ICFI(1,I+6)=1                                                     STE19870
      ICFI(2,I+6)=I                                                     STE19880
 40   CONTINUE                                                          STE19890
C                                                                       STE19900
C                                                                       STE19910
      CALL PSIROT(PSI,G( 7,1),NRDG,18,G( 7,1),NRDG)                     STE19920
      CALL PSIROT(PSI,G(10,1),NRDG,18,G(10,1),NRDG)                     STE19930
C                                                                       STE19940
C                                                                       STE19950
      COSP=DCOS(PSI)                                                    STE19960
      SINP=DSIN(PSI)                                                    STE19970
      DO 50 I=1,4                                                       STE19980
      I1=(I-1)*3+1                                                      STE19990
      DO 50 K=1,2                                                       STE20000
      K1=(K-1)*3+13                                                     STE20010
      TMP1= G(I1,K1)*COSP-G(I1,K1+1)*SINP                               STE20020
      TMP2= G(I1,K1)*SINP+G(I1,K1+1)*COSP                               STE20030
      G(I1,K1  )=TMP1                                                   STE20040
      G(I1,K1+1)=TMP2                                                   STE20050
      TMP1= G(I1+1,K1)*COSP-G(I1+1,K1+1)*SINP                           STE20060
      TMP2= G(I1+1,K1)*SINP+G(I1+1,K1+1)*COSP                           STE20070
      G(I1+1,K1  )=TMP1                                                 STE20080
      G(I1+1,K1+1)=TMP2                                                 STE20090
  50  CONTINUE                                                          STE20100
      CALL DRMOVE(G(1,19),G(1,15),12,DIST)                              STE20110
C***********CORREZIONE DEL 10/1/90 PER ANALOGIA ALLA STIFF              STE20120
C************DELLA NUOVA FORMULAZIONE*********************              STE20130
C                                                                       STE20140
      CALL DDEMA(RQ,RQ(7),GAMDOT,3,1.D0)                                STE20150
      CALL DLGMA(RQ,GAMFIR,3,1.D0,-1)                                   STE20160
      CALL PRM3(GAMDOT,3,GAMFIR,3,GAMDOT,3,0,0)                         STE20170
      CALL DLGMA(RQ,GAMFIR,3,1.D0,0)                                    STE20180
C      CALL DWRITE('GAMDOT',GAMDOT,3,3,3,1)                             STE20190
C      CALL DWRITE('GAMFIR',GAMFIR,3,3,3,1)                             STE20200
C**************************************************************         STE20210
C    U N K N O W N    S U B S T I T U T I O N                           STE20220
C    ****************************************                           STE20230
      DO 36 JR=1,12,3                                                   STE20240
      CALL DPRMMA(G(JR,4),NRDG,3,3,GAMDOT,3,G(JR,10),NRDG,3,0,0)        STE20250
 36   CONTINUE                                                          STE20260
      DO 41 JR=1,12,3                                                   STE20270
         DO 41 JC=4,12,6                                                STE20280
         CALL PRM3(G(JR,JC),NRDG,GAMFIR,3,G(JR,JC),NRDG,0,0)            STE20290
 41   CONTINUE                                                          STE20300
C******************FINE CORREZIONE****************************          STE20310
      RETURN                                                            STE20320
      END                                                               STE20330
C*************************          COE0          **********************STE22920
C=    COMPILER (LINK=IBJ$)                                              STE22930
      SUBROUTINE COE0(VCSTR,OUTA,CS,                                    STE22940
     *                CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,              STE22950
     *                CRF,CFSLOP,DCPDM,DCRDRM,DCMDM,DCRFDM)             STE22960
      IMPLICIT REAL*8(A-H,O-Z)                                          STE22970
      DIMENSION VCSTR(6),OUTA(20)                                       STE22980
      DATA PG/3.1415926D0/,DEGRAD/.017453293D0/                         STE22990
C                                                                       STE23000
C     CALCOLA I COEFFICIENTI AERODINAMICI CON TEORIA STAZIONARIA        STE23010
C     (CON CORREZIONE PER FRECCIA SECONDO HARRIS A.H.S.JULY 1970)       STE23020
C                                                                       STE23030
      ALFA=DATAN2(-VCSTR(2),VCSTR(1))                                   STE23040
      OUTA(2)=ALFA/DEGRAD                                               STE23050
      VC1=DABS(VCSTR(1))                                                STE23060
      GAM=DATAN2(-VCSTR(3),VC1)                                         STE23070
      OUTA(3)=GAM/DEGRAD                                                STE23080
      IF(DABS(GAM).GT.PG/3.D0) GAM=PG/3.D0                              STE23090
      COSGAM=DCOS(GAM)                                                  STE23100
      RMACH=DSQRT(VCSTR(1)**2+VCSTR(2)**2+VCSTR(3)**2)/CS               STE23110
      RMACH=RMACH*DSQRT(COSGAM)                                         STE23120
      OUTA(4)=RMACH                                                     STE23130
      CALL CPCRCM(ALFA,RMACH,                                           STE23140
     *            CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,                  STE23150
     *            CRF,CFSLOP,DCPDM,DCRDRM,DCMDM,DCRFDM,                 STE23160
     *            ASLOP0,CSLOP0)                                        STE23170
      IF(DABS(ALFA*COSGAM).LT.1.D-6) GO TO 10                           STE23190
      ASLRF=CLIFT/(ALFA*COSGAM)                                         STE23200
      IF(ASLRF.GT.ASLOP0) CLIFT=ASLOP0*ALFA                             STE23210
 10   continue                                                          STE23220
      OUTA(5)=CLIFT                                                     STE23230
      OUTA(6)=CDRAG                                                     STE23240
      OUTA(7)=CMOME                                                     STE23250
C                                                                       STE23260
      RETURN                                                            STE23270
C                                                                       STE23280
C                                                                       STE23290
      END                                                               STE23300
C*************************          COE1          **********************STE23310
C=    COMPILER (LINK=IBJ$)                                              STE23320
      SUBROUTINE COE1(VCSTR,OUTA,CS,CORDA,RSPEED,DA,                    STE23330
     *                CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,              STE23340
     *                CRF,CFSLOP,DCPDM,DCRDRM,DCMDM,DCRFDM)             STE23350
      IMPLICIT REAL*8(A-H,O-Z)                                          STE23360
      DIMENSION VCSTR(6),OUTA(20)                                       STE23370
      DATA PG/3.1415926D0/,DEGRAD/.017453293D0/                         STE23380
C                                                                       STE23390
C     CALCOLA COEFFICIENTI AERODINAMICI CON TEORIA INSTAZIONARIA        STE23400
C     SECONDO TEORIA HARRIS A.H.S. JULY 1970                            STE23410
C                                                                       STE23420
      ABE=DATAN2(-VCSTR(2),VCSTR(1))                                    STE23430
      OUTA(2)=ABE/DEGRAD                                                STE23440
      VC1=DABS(VCSTR(1))                                                STE23450
      ALF1=OUTA(9)/RSPEED                                               STE23470
      OUTA(9)=ALF1                                                      STE23480
      ALF2=OUTA(10)/(RSPEED*RSPEED)                                     STE23490
      OUTA(10)=ALF2                                                     STE23500
      VP=DSQRT(VCSTR(1)**2+VCSTR(2)**2)                                 STE23510
      GAM=DATAN2(-VCSTR(3),VC1)                                         STE23520
      OUTA(3)=GAM/DEGRAD                                                STE23530
      IF(DABS(GAM).GT.PG/3.D0) GAM=PG/3.D0                              STE23540
      COSGAM=DCOS(GAM)                                                  STE23550
      RMACH=DSQRT(VCSTR(1)**2+VCSTR(2)**2+VCSTR(3)**2)/CS               STE23560
      RMACH=RMACH*DSQRT(COSGAM)                                         STE23570
      OUTA(4)=RMACH                                                     STE23580
      FR=.5D0*CORDA*RSPEED/VP                                           STE23590
      FR=DABS(FR)                                                       STE23600
      OUTA(11)=FR                                                       STE23610
      RMM=.3D0                                                          STE23620
      IF(RMACH.GT..3D0) RMM=RMACH                                       STE23630
      IF(RMACH.GT..6D0) RMM=.6D0                                        STE23640
      DALF=61.5D0*DLOG(.6D0/RMM)*DSQRT(FR)*PG/180.D0                    STE23650
      SEGNO=1.D0                                                        STE23660
      IF(ALF1.LT.0.D0) SEGNO=-1.D0                                      STE23670
       if (dabs(abe*cosgam).lt.(.22689*(1-rmach*rmach))) then
         AREF=abe*cosgam
       else	
	 amod=.2*dabs(abe*cosgam)
	 if (amod.lt.dalf) dalf=amod
         AREF=ABE-SEGNO*DALF                                            STE23680
       endif
      OUTA(12)=AREF/DEGRAD                                              STE23690
      CALL CPCRCM(AREF,RMACH,                                           STE23700
     *            CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,                  STE23710
     *            CRF,CFSLOP,DCPDM,DCRDRM,DCMDM,DCRFDM,                 STE23720
     *            ASLOP0,CSLOP0)                                        STE23730
      IF(DABS(AREF*COSGAM).LT.1.D-6) GO TO 10                           STE23750
      ASLRF=CLIFT/(AREF*COSGAM)                                         STE23760
      IF(ASLRF.GT.ASLOP0) CLIFT=ASLOP0*AREF
 10   CONTINUE                                                          STE23780
      RK=.5*CORDA*RSPEED/VP                                             STE23790
      PA=.25D0                                                          STE23800
      ALF=THF(RK)*ABE+(.5D0*RK+THG(RK))*ALF1                            STE23810
      ALF=ALF+2.D0*(.75D0-PA)*THF(RK)*RK*ALF1                           STE23820
      ALF=ALF-RK*RK*(PA-.5D0)*ALF2                                      STE23830
      if (ALF.gt.1.2*AREF) ALF=1.2*AREF
      if (ALF.lt.0.8*AREF) ALF=.8*AREF
      CLIFT=CLIFT*ALF/(AREF*COSGAM)
      OUTA(5)=CLIFT                                                     STE23850
      OUTA(6)=CDRAG                                                     STE23860
      OUTA(7)=CMOME                                                     STE23870
      OUTA(13)=ALF                                                      STE23880
      outa(14)=thf(rk)
      outa(15)=thg(rk)
      outa(16)=RSPEED
C                                                                       STE23890
      RETURN                                                            STE23900
C                                                                       STE23910
C                                                                       STE23920
      END                                                               STE23930
C*************************          COE2          **********************STE23940
C=    COMPILER (LINK=IBJ$)                                              STE23950
      SUBROUTINE COE2(VCSTR,OUTA,CS,CORDA,RSPEED,DA,                    STE23960
     *                CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,              STE23970
     *                CRF,CFSLOP,DCPDM,DCRDRM,DCMDM,DCRFDM)             STE23980
      IMPLICIT REAL*8(A-H,O-Z)                                          STE23990
      DIMENSION VCSTR(6),OUTA(20),PN(14),QN(14),PM(14),QM(14)           STE24000
      DATA PG/3.1415926D0/,DEGRAD/.017453293D0/                         STE24010
      DATA ASN0/.22689D0/,ASM0/.22689D0/                                STE24020
      DATA PN/-3.464003D-1,-1.549076D0,4.306330D+1,-5.397529D+1,        STE24030
     *        5.781402D0,-3.233003D+1,-2.162257D+1,1.866347D+1,         STE24040
     *        4.198390D+1,3.295461D+2,4*0.D0/                           STE24050
      DATA QN/1.533717D0,6.977203D+0,1.749010D+3,1.694829D 3,           STE24060
     *        -1.771899D+3,-3.291665D+4,2.969051D0,-3.632448D+1,        STE24070
     *        -2.268578D+3,6.601995D 3,-9.654208D+3,8.533930D+4,        STE24080
     *        -1.492624D0,1.163661D1/                                   STE24090
      DATA PM/1.970065D1,-6.751639D1,7.265269D2,4.865945D4,             STE24100
     *        2.086279D4,6.024672D3,1.446334D2,8.586896D2,              STE24110
     *        -7.550329D2,-1.021613D1,2.247664D1,3*0.D0/                STE24120
      DATA QM/-2.322808D0,-1.322257D0,-2.633891D0,-2.180321D-1,         STE24130
     *        4.580014D0,3.125497D-1,-2.828806D1,-4.396734D0,           STE24140
     *        2.565870D2,-1.204976D1,-1.157802D2,8.612138D0,            STE24150
     *        2*0.D0/                                                   STE24160
C                                                                       STE24170
C     CALCOLA I COEFFICIENTI AERODINAMICI CON TEORIA INSTAZIONARIA      STE24180
C     CON DATI SINTETIZZATI DI BIELAWA 31TH A.H.S. FORUM 1975           STE24190
C                                                                       STE24200
      ALFA=DATAN2(-VCSTR(2),VCSTR(1))                                   STE24210
      OUTA(2)=ALFA/DEGRAD                                               STE24220
      ALF1=OUTA(9)/DA                                                   STE24230
      OUTA(9)=ALF1                                                      STE24240
      ALF2=OUTA(10)/(DA*DA)                                             STE24250
      OUTA(10)=ALF2                                                     STE24260
      VP=DSQRT(VCSTR(1)**2+VCSTR(2)**2)                                 STE24270
      A=.5*CORDA*ALF1/VP                                                STE24280
      B=.25*CORDA*CORDA*ALF2/(VP*VP)                                    STE24290
      VC1=DABS(VCSTR(1))                                                STE24300
      GAM=DATAN2(-VCSTR(3),VC1)                                         STE24310
      OUTA(3)=GAM/DEGRAD                                                STE24320
      IF(DABS(GAM).GT.PG/3.D0) GAM=PG/3.D0                              STE24330
      COSGAM=DCOS(GAM)                                                  STE24340
      RMACH=DSQRT(VCSTR(1)**2+VCSTR(2)**2+VCSTR(3)**2)/CS               STE24350
      RMACH=RMACH*DSQRT(COSGAM)                                         STE24360
      OUTA(4)=RMACH                                                     STE24370
      IF(RMACH.GE.1.D0) RMACH=.99D0                                     STE24380
      ETA=DSQRT((A/.048)**2+(B/.016)**2)                                STE24390
      SEGNO=1.D0                                                        STE24400
      IF(ALFA.LT.0D0) SEGNO=-1.D0                                       STE24410
      IF(ETA.GT.1.D0) SEGNO=SEGNO/ETA                                   STE24420
      A=SEGNO*A                                                         STE24430
      B=SEGNO*B                                                         STE24440
      A2=A*A                                                            STE24450
      B2=B*B                                                            STE24460
      ASN=ASN0*(1.D0-RMACH)                                             STE24470
      ASM=ASM0*(1.D0-RMACH)                                             STE24480
      SIGN=DABS(ALFA/ASN)                                               STE24490
      SIGM=DABS(ALFA/ASM)                                               STE24500
      SIGMAX=1.839-70.33*B                                              STE24510
      IF(SIGMAX.GT.1.86D0) SIGMAX=1.86D0                                STE24520
      IF(SIGN.GT.SIGMAX) SIGN=SIGMAX                                    STE24530
      IF(SIGM.GT.SIGMAX) SIGM=SIGMAX                                    STE24540
      DAN=A*(PN(1)+PN(5)*SIGN)+B*(PN(2)+PN(6)*SIGN)                     STE24550
      DAN=DAN+DEXP(-1072.52*A2)*(A*(PN(3)+PN(7)*SIGN)+                  STE24560
     *      A2*(PN(9)+PN(10)*SIGN))                                     STE24570
      DAN=DAN+DEXP(-40316.42D0*B2)*B*(PN(4)+PN(8)*SIGN)                 STE24580
      DAN=DAN*ASN                                                       STE24590
      DCN=A*(QN(1)+QN(3)*A2+SIGN*(QN(7)+QN(9)*A2+QN(13)*SIGN)+          STE24600
     *      B2*(QN(5)+QN(11)*SIGN))                                     STE24610
      DCN=DCN+B*(QN(2)+QN(4)*A2+SIGN*(QN(8)+QN(10)*A2+QN(14)*SIGN)      STE24620
     *      +B2*(QN(6)+QN(12)*SIGN))                                    STE24630
      DAM=A*(PM(1)+PM(3)*A2+PM(5)*B2+PM(10)*SIGM+PM(7)*A)               STE24640
      DAM=DAM+B*(PM(2)+PM(4)*B2+PM(6)*A2+PM(11)*SIGM+PM(8)*B+PM(9)*A)   STE24650
      DAM=DAM*ASM                                                       STE24660
      S2=SIGM*SIGM                                                      STE24670
      DCM=A*(QM(2)+QM(8)*A+SIGM*(QM(4)+QM(10)*A)+S2*(QM(6)+QM(12)*A))   STE24680
      DCM=DCM+B*(QM(1)+QM(7)*B+SIGM*(QM(3)+QM(9)*B)+S2*(QM(5)+QM(11)*B))STE24690
      OUTA(11)=DAN/DEGRAD                                               STE24700
      OUTA(12)=DAM/DEGRAD                                               STE24710
      OUTA(13)=DCN                                                      STE24720
      OUTA(14)=DCM                                                      STE24730
      CALL CPCRCM(ALFA-DAN,RMACH,                                       STE24740
     *            CLIFT,CDRAG,DUM,ASLOP,BSLOP,DUM,                      STE24750
     *            CRF,CFSLOP,DCPDM,DCRDRM,DUM,DCRFDM,                   STE24760
     *            ASLOP0,CSLOP0)                                        STE24770
      CALL CPCRCM(ALFA-DAM,RMACH,                                       STE24780
     *            CL,DUM,CMOME,DUM,DUM,CSLOP,                           STE24790
     *            DUM,DUM,DUM,DUM,DCMDM,DUM,                            STE24800
     *            ASLOP0,CSLOP0)                                        STE24810
      IF(DABS((ALFA-DAN)*COSGAM).LT.1.D-6) GO TO 10                     STE24830
      ASLRF=CLIFT/((ALFA-DAN)*COSGAM)                                   STE24840
      IF(ASLRF.GT.ASLOP0) CLIFT=ASLOP0*(ALFA-DAN)                       STE24850
 10   CONTINUE                                                          STE24860
      C1=.9457/DSQRT(1.D0-RMACH*RMACH)                                  STE24870
      CLIFT=CLIFT+ASLOP0*DAN+DCN*C1                                     STE24880
      CMOME=CMOME+CSLOP0*DAM+DCM*C1                                     STE24890
      OUTA(5)=CLIFT                                                     STE24900
      OUTA(6)=CDRAG                                                     STE24910
      OUTA(7)=CMOME                                                     STE24920
C                                                                       STE24930
      RETURN                                                            STE24940
C                                                                       STE24950
C                                                                       STE24960
      END                                                               STE24970
C*************************          CPCRCM          ********************STE24980
C=    COMPILER (LINK=IBJ$)                                              STE24990
      SUBROUTINE CPCRCM (APHIJ,EMIJ,CLIFT,CDRAG,CMOME,ASLOP,BSLOP,CSLOP,STE25000
     1                  CRF,CFSLOP,DCPDM,DCDRDM,DCMDM,DCRFDM,AS0,CS0)   STE25010
      IMPLICIT REAL*8 (A-H,O-Z)                                         STE25020

c     il common serve per passare jpro
      COMMON/DATAERO/UMEAN,ro,cs,pesoU,rrot,AREA,inst,
     *                 imeth,jpro,mozzo,veliv
      real*8  UMEAN,ro,cs,pesoU,rrot,AREA
      integer*4  inst,imeth,jpro,mozzo,veliv


      DIMENSION D1(2,5),D2(2,5),D3(2,5),CD0(2,5),B1(2,5),B2(2,5),B3(2,5)STE25050
      DATA D1/.3,-.28532D2,.4,-.37392D2,.45,-.4597D2,.5,-.55819D2,      STE25060
     *        .55,-.88102D2/                                            STE25070
      DATA D2/.3,.4826D1,.4,.70982D1,.45,.8914D1,.5,.10636D2,           STE25080
     *        .55,.1993D2/                                              STE25090
      DATA D3/.3,.62085D1,.4,.64158D1,.45,.65334D1,.5,.66613D1,         STE25100
     *        .55,.64065D1/                                             STE25110
      DATA CD0/.3,.60145D-2,.4,.35772D-2,.45,.2612D-2,.5,.46674D-2,     STE25120
     *        .55,.29977D-2/                                            STE25130
      DATA B1/.3,.86845D-1,.4,.1915D0,.45,.24923D0,.5,.14568D0,         STE25140
     *        .55,.25389D0/                                             STE25150
      DATA B2/.3,-.91795D0,.4,-.1966D1,.45,-.26795D1,.5,-.15943D1,      STE25160
     *        .55,-.30426D1/                                            STE25170
      DATA B3/.3,.40357D1,.4,.67937D1,.45,.91332D1,.5,.62569D1,         STE25180
     *        .55,.12049D2/                                             STE25190
C**** APHIJ  = ANGOLO INCIDENZA                                         STE25200
C**** EMIJ   = NUMERO DI MACH                                           STE25210
C**** CLIFT  = COEFFICIENTE DI PORTANZA                                 STE25220
C**** CDRAG  = COEFFICIENTE DI RESISTENZA TOTALE                        STE25230
C**** CMOME  = COEFFICIENTE DI MOMENTO (CONVENZIONE CABRANTE) RIFERITO  STE25240
C****          AL CENTRO AERODINAMICO                                   STE25250
C**** ASLOP  = PENDENZA CURVA DI PORTANZA                               STE25260
C**** BSLOP  = PENDENZA CURVA DI RESISTENZA TOTALE                      STE25270
C**** CSLOP  = PENDENZA CURVA DI MOMENTO                                STE25280
C**** CRF    = COEFFICIENTE DI RESISTENZA ATTRITO                       STE25290
C**** CFSLOP = PENDENZA CURVA DI RESISTENZA ATTRITO                     STE25300
C**** DCPDM  = DERIVATA DI CLIFT RISPETTO NUMERO DI MACH                STE25310
C**** DCDRDM = DERIVATA DI CDRAG RISPETTO NUMERO DI MACH                STE25320
C**** DCMDM  = DERIVATA DI CMOME RISPETTO NUMERO DI MACH                STE25330
C**** DCRFDM = DERIVATA DI CRF RISPETTO NUMERO DI MACH                  STE25340
C**** AS0    = PENDENZA CURVA PORTANZA PER ALFA=0                       STE25350
C**** CS0    = PENDENZA CURVA MOMENTO PER ALFA=0                        STE25360
      PG=3.14159                                                        STE25370
      CFSLOP=0.                                                         STE25380
      CRF=.006                                                          STE25390
      DCRFDM=0.                                                         STE25400
      NEG=1                                                             STE25410
      IF(EMIJ.GT.1.D0)EMIJ=.99D0                                        STE25420
      SQT=DSQRT(1.-EMIJ*EMIJ)                                           STE25430
      C1=1.-EMIJ                                                        STE25440
      C2=.22689*C1                                                      STE25450
      C5=EMIJ/(SQT*SQT)                                                 STE25460
97    IF(APHIJ)181,182,182                                              STE25470
181   APHIJ=-APHIJ                                                      STE25480
      NEG=-1*NEG                                                        STE25490
182   IF(APHIJ-PG)184,184,183                                           STE25500
183   APHIJ=APHIJ-PG*2.                                                 STE25510
      GOTO97                                                            STE25520
184   IF(APHIJ-C2)185,187,187                                           STE25530
185   ASLOP=5.7296/SQT                                                  STE25540
      CLIFT=ASLOP*APHIJ                                                 STE25550
      CDRAG=.006+.13131*APHIJ*APHIJ                                     STE25560
      BSLOP=.26262*APHIJ                                                STE25570
      CMOME=1.4324*APHIJ/SQT                                            STE25580
      CSLOP=1.4324/SQT                                                  STE25590
      DCPDM=CLIFT*C5                                                    STE25600
      DCDRDM=0.                                                         STE25610
      DCMDM=CMOME*C5                                                    STE25620
      GOTO250                                                           STE25630
187   IF(APHIJ-.34906)189,191,191                                       STE25640
189   CLIFT=.29269*C1+(1.3*EMIJ-.59)*APHIJ                              STE25650
      C2=(.12217+.22689*EMIJ)*SQT                                       STE25660
      CMOME=CLIFT/(4*C2)                                                STE25670
      CSLOP=(1.3*EMIJ-.59)/(4*C2)                                       STE25680
      CLIFT=CLIFT/C2                                                    STE25690
      ASLOP=(1.3*EMIJ-.59)/C2                                           STE25700
      DCPDM=(-.29269+1.3*APHIJ-CLIFT*(-(.12217+.22689*EMIJ)*EMIJ/SQT+   STE25710
     1      .22689*SQT))/C2                                             STE25720
      DCMDM=.25*DCPDM                                                   STE25730
      GOTO210                                                           STE25740
191   IF(APHIJ-2.7402)193,195,195                                       STE25750
193   S=DSIN(APHIJ)                                                     STE25760
      S2=DSIN(2.*APHIJ)                                                 STE25770
      S3=DSIN(3.*APHIJ)                                                 STE25780
      S4=DSIN(4.*APHIJ)                                                 STE25790
      CLIFT=(.080373*S+1.04308*S2-.011059*S3+.023127*S4)/SQT            STE25800
      CMOME=(-.02827*S+.14022*S2-.00622*S3+.01012*S4)/SQT               STE25810
      C1=DCOS(APHIJ)                                                    STE25820
      C2=DCOS(2.*APHIJ)                                                 STE25830
      C3=DCOS(3.*APHIJ)                                                 STE25840
      C4=DCOS(4.*APHIJ)                                                 STE25850
      CSLOP=(-.02827*C1+.28044*C2-.01866*C3+.04048*C4)/SQT              STE25860
      ASLOP=(.080373*C1+2.08616*C2-.033177*C3+.092508*C4)/SQT           STE25870
      DCPDM=CLIFT*C5                                                    STE25880
      DCMDM=CMOME*C5                                                    STE25890
      GOTO210                                                           STE25900
195   IF(APHIJ-3.0020)197,199,199                                       STE25910
197   CLIFT=-(.4704+.10313*APHIJ)/SQT                                   STE25920
      ASLOP=-.10313/SQT                                                 STE25930
      CMOME=-(.4786+.02578*APHIJ)/SQT                                   STE25940
      CSLOP=.02578/SQT                                                  STE25950
      DCPDM=CLIFT*C5                                                    STE25960
      DCMDM=CMOME*C5                                                    STE25970
      GOTO210                                                           STE25980
199   IF(APHIJ-PG)200,200,260                                           STE25990
200   CLIFT=(-17.55+5.5864*APHIJ)/SQT                                   STE26000
      ASLOP=5.5864/SQT                                                  STE26010
      CMOME=(-12.5109+3.9824*APHIJ)/SQT                                 STE26020
      CSLOP=3.9824/SQT                                                  STE26030
      DCPDM=CLIFT*C5                                                    STE26040
      DCMDM=CMOME*C5                                                    STE26050
210   CDRAG=1.1233-.029894*DCOS(APHIJ)-1.00603*DCOS(2.*APHIJ)           STE26060
      CDRAG=CDRAG+.003115*DCOS(3.*APHIJ)-.091487*DCOS(4.*APHIJ)         STE26070
      CDRAG=CDRAG/SQT                                                   STE26080
      BSLOP=.029894*DSIN(APHIJ)+2.01206*DSIN(2.*APHIJ)                  STE26090
      BSLOP=(BSLOP+.009345*DSIN(3.*APHIJ)+.365948*DSIN(4.*APHIJ))/SQT   STE26100
      DCDRDM=CDRAG*C5                                                   STE26110
 250  GOTO(251,252),JPRO                                                STE26120
C                                                                       STE26130
C     CALCOLO COEFFICIENTI AERODINAMICI PER PROFILO 'RAE9671'           STE26140
C                                                                       STE26150
 252  IF(DABS(APHIJ).GT..21956242) GOTO 251                             STE26160
      IF(NEG.LT.0) APHIJ = -APHIJ                                       STE26170
      CALL LININT(D1,5,EMIJ,P1,PEND1)                                   STE26180
      CALL LININT(D2,5,EMIJ,P2,PEND2)                                   STE26190
      CALL LININT(D3,5,EMIJ,P3,PEND3)                                   STE26200
      APHIJ2=APHIJ*APHIJ                                                STE26210
      APHIJ3=APHIJ*APHIJ2                                               STE26220
      CLIFT=P1*APHIJ3+P2*APHIJ2+P3*APHIJ                                STE26230
      ASLOP=3.*P1*APHIJ2+2.*P2*APHIJ+P3                                 STE26240
      AS0=P3                                                            STE26250
      DCPDM=PEND1*APHIJ3+PEND2*APHIJ2+PEND3*APHIJ                       STE26260
C                                                                       STE26270
      CALL LININT(CD0,5,EMIJ,R0,PEND0)                                  STE26280
      CRF=R0                                                            STE26290
      CFSLOP=0.                                                         STE26300
      DCRFDM=PEND0                                                      STE26310
      CALL LININT(B1,5,EMIJ,R1,PEND1)                                   STE26320
      CALL LININT(B2,5,EMIJ,R2,PEND2)                                   STE26330
      CALL LININT(B3,5,EMIJ,R3,PEND3)                                   STE26340
      APHIJ=DABS(APHIJ)                                                 STE26350
      APHIJ3=DABS(APHIJ3)                                               STE26360
      CDRAG=R0+APHIJ*R1+APHIJ2*R2+APHIJ3*R3                             STE26370
      BSLOP=R1+2.*R2*APHIJ+3.*R3*APHIJ2                                 STE26380
      DCDRDM=PEND0+APHIJ*PEND1+APHIJ2*PEND2+APHIJ3*PEND3                STE26390
C                                                                       STE26400
      IF(NEG.GT.0) GOTO 249                                             STE26410
      CMOME=-CMOME                                                      STE26420
      CSLOP=-CSLOP                                                      STE26430
      DCMDM=-DCMDM                                                      STE26440
 249  CONTINUE                                                          STE26450
      CS0=0.                                                            STE26460
      CMOME=CMOME-.25*CLIFT                                             STE26470
      CSLOP=CSLOP-.25*ASLOP                                             STE26480
      DCMDM=DCMDM-.25*DCPDM                                             STE26490
      RETURN                                                            STE26500
 251  CONTINUE                                                          STE26510
      IF(NEG)255,255,260                                                STE26520
255   CLIFT=-CLIFT                                                      STE26530
      CMOME=-CMOME                                                      STE26540
      APHIJ=-APHIJ                                                      STE26550
      DCPDM=-DCPDM                                                      STE26560
      DCMDM=-DCMDM                                                      STE26570
260   CONTINUE                                                          STE26580
      CMOME=CMOME-.25*CLIFT                                             STE26590
      CSLOP=CSLOP-.25*ASLOP                                             STE26600
      DCMDM=DCMDM-.25*DCPDM                                             STE26610
      AS0=5.7296/SQT                                                    STE26620
      CS0=0.D0                                                          STE26630
      RETURN                                                            STE26640
      END                                                               STE26650
C*************************          PSIROT          ********************STE27750
C=    COMPILER (LINK=IBJ$)                                              STE27760
      SUBROUTINE PSIROT(PSI,A,NRDA,NC,B,NRDB)                           STE27770
      DOUBLE PRECISION PSI,A(NRDA,1),B(NRDB,1)                          STE27780
      DOUBLE PRECISION SINP,COSP,TMP1,TMP2                              STE27790
      SINP=DSIN(PSI)                                                    STE27800
      COSP=DCOS(PSI)                                                    STE27810
      DO 20 K=1,NC                                                      STE27820
      TMP1=+COSP*A(1,K)-SINP*A(2,K)                                     STE27830
      TMP2=+SINP*A(1,K)+COSP*A(2,K)                                     STE27840
      B(1,K)=TMP1                                                       STE27850
      B(2,K)=TMP2                                                       STE27860
      B(3,K)=A(3,K)                                                     STE27870
  20  CONTINUE                                                          STE27880
  30  CONTINUE                                                          STE27890
      RETURN                                                            STE27900
      END                                                               STE27910
C*************************          THF F          *********************STE27920
C=    COMPILER (LINK=IBJ$)                                              STE27930
      FUNCTION THF(K)                                                   STE27940
      REAL*8 THF,K,A,B,C,D                                              STE27950
C ****  FUNZIONE F DI THEODORSEN ***                                    STE27960
      A=K**4-.102058*K**2+9.55732E-6                                    STE27970
      B=-.761036*K**3+2.55067E-3*K                                      STE27980
      C= 2.*K**4-.113928*K**2+9.55732E-6                                STE27990
      D=-1.064*K**3+2.6268E-3*K                                         STE28000
      THF=(A*C+B*D)/(C*C+D*D)                                           STE28010
      RETURN                                                            STE28020
      END                                                               STE28030
C*************************          THG F          *********************STE28040
C=    COMPILER (LINK=IBJ$)                                              STE28050
      FUNCTION THG(K)                                                   STE28060
      REAL*8 THG,K,A,B,C,D                                              STE28070
C ****  FUNZIONE G DI THEODORSEN ***                                    STE28080
      A=K**4-.102058*K**2+9.55732E-6                                    STE28090
      B=-.761036*K**3+2.55067E-3*K                                      STE28100
      C= 2.*K**4-.113928*K**2+9.55732E-6                                STE28110
      D=-1.064*K**3+2.6268E-3*K                                         STE28120
      THG=(B*C-A*D)/(C*C+D*D)                                           STE28130
      RETURN                                                            STE28140
      END                                                               STE28150

	subroutine polcoe(x,y,n,cof)
	real*8 x(1),y(1),cof(1),s(10)
	real*8 phi,b,ff
	integer*4 n,i,j,k
	do i=1,n
	 s(i)=0.
	 cof(i)=0.
	enddo
	s(n)=-x(1)
	do i=2,n
	 do j=n+1-i,n-1
	  s(j)=s(j)-x(i)*s(j+1)
	 enddo
	 s(n)=s(n)-x(i)
	enddo
	do j=1,n
	 phi=n
	 do k=n-1,1,-1
	  phi=k*s(k+1)+x(j)*phi
	 enddo
	 ff=y(j)/phi
	 b=1.
	 do k=n,1,-1
	  cof(k)=cof(k)+b*ff
	  b=s(k)+x(j)*b
	 enddo
	enddo
	return
	end






















*********************************  STLIB  ***************************************
C*************************          ADD          ***********************STL00010
C=    COMPILER (LINK=IBJ$)                                              STL00020
      SUBROUTINE ADD(A,NRDA,B,NRDB,NR,NC,PLUS)                          STL00030
      DOUBLE PRECISION A(NRDA,1),B(NRDB,1),PLUS                         STL00040
      DO 10 I=1,NR                                                      STL00050
      DO 10 K=1,NC                                                      STL00060
 10   A(I,K)=A(I,K)+PLUS*B(I,K)                                         STL00070
      RETURN                                                            STL00080
      END                                                               STL00090

C*************************          DABCT3          ********************STL00280
C=    COMPILER (LINK=IBJ$)                                              STL00290
      SUBROUTINE DABCT3(A,B,NRDB,C)                                     STL00300
C                                                                       STL00310
C     CALCOLA IL PRODOTTO DELLE MATRICI (3*3) "A" * "B" * "C-TRASPOSTA" STL00320
C                                                                       STL00330
      DOUBLE PRECISION A,B,C,D                                          STL00340
      DIMENSION A(3,3),B(NRDB,1),C(3,3),D(3,3)                          STL00350
      DO 10 I=1,3                                                       STL00360
      DO 10 K=1,3                                                       STL00370
      D(I,K)=0.                                                         STL00380
      DO 10 L=1,3                                                       STL00390
10    D(I,K)=D(I,K)+A(I,L)*B(L,K)                                       STL00400
      DO 20 I=1,3                                                       STL00410
      DO 20 K=1,3                                                       STL00420
      B(I,K)=0.                                                         STL00430
      DO 20 L=1,3                                                       STL00440
20    B(I,K)=B(I,K)+D(I,L)*C(K,L)                                       STL00450
      RETURN                                                            STL00460
      END                                                               STL00470

C     ******************          DDEMA          ***********************STL10500
      SUBROUTINE DDEMA(RO,RD,DEMA,ND,FACT)                              STL10510
C                                                                       STL10520
C     GIVEN THE ROTATION VECTOR RO AND THE DERIVATIVE OF THE ROTATION   STL10530
C     VECTOR RD, CALCULATE THE CAPITAL GAMMA DOT.                       STL10540
C                                                                       STL10550
      IMPLICIT DOUBLE PRECISION  (A-H,O-Z)                              STL10560
      REAL*8 A,B,C,D,E,G,RDRO,ARG,ARG2                                  STL10570
      DIMENSION RO(3),RD(3),DEMA(ND,3),T(3)                             STL10580
      DATA TRESH2/1.D-10/                                               STL10590
      ARG2=0.D0                                                         STL10600
      RDRO=0.D0                                                         STL10610
      DO 10 I=1,3                                                       STL10620
      T(I)=RO(I)*FACT                                                   STL10630
      ARG2=ARG2+T(I)**2                                                 STL10640
      RDRO=RDRO+T(I)*RD(I)                                              STL10650
   10 CONTINUE                                                          STL10660
      IF (ARG2.LT.TRESH2) THEN                                          STL10670
C                                                                       STL10680
        A=(1.D0-ARG2/30.D0*(2.D0+ARG2/56.D0*(3.D0-ARG2/22.5D0)))/12.D0  STL10690
        B=(1.D0-ARG2/42.D0*(2.D0-ARG2/72.D0*(3.D0-ARG2/27.5D0)))/60.D0  STL10700
        D=(1.D0-ARG2/12.D0*(1.D0-ARG2/30.D0*(1.D0-ARG2/56.0D0)))/ 2.D0  STL10710
        E=(1.D0-ARG2/20.D0*(1.D0-ARG2/42.D0*(1.D0-ARG2/72.0D0)))/ 6.D0  STL10720
        A=-A*RDRO                                                       STL10730
        B=-B*RDRO                                                       STL10740
C                                                                       STL10750
      ELSE                                                              STL10760
C                                                                       STL10770
         ARG=DSQRT(ARG2)                                                STL10780
         D=(1.D0-DCOS(ARG))/ARG2                                        STL10790
         G=DSIN(ARG)/ARG                                                STL10800
         E=(1.D0-G)/ARG2                                                STL10810
         A=(G-2.D0*D)*RDRO/ARG2                                         STL10820
         B=(D-3.D0*E)*RDRO/ARG2                                         STL10830
C                                                                       STL10840
      END IF                                                            STL10850
C                                                                       STL10860
      C=-B*ARG2                                                         STL10870
      DEMA(1,1)=T(1)*T(1)*B+C                                           STL10880
      DEMA(1,2)=T(1)*T(2)*B-T(3)*A                                      STL10890
      DEMA(1,3)=T(1)*T(3)*B+T(2)*A                                      STL10900
      DEMA(2,1)=T(2)*T(1)*B+T(3)*A                                      STL10910
      DEMA(2,2)=T(2)*T(2)*B+C                                           STL10920
      DEMA(2,3)=T(2)*T(3)*B-T(1)*A                                      STL10930
      DEMA(3,1)=T(3)*T(1)*B-T(2)*A                                      STL10940
      DEMA(3,2)=T(3)*T(2)*B+T(1)*A                                      STL10950
      DEMA(3,3)=T(3)*T(3)*B+C                                           STL10960
      DEMA(1,1)=DEMA(1,1)+ 2.D0*(RD(1)*T(1)-RDRO)*E                     STL10970
      DEMA(1,2)=DEMA(1,2)+(RD(1)*T(2)+T(1)*RD(2))*G-RD(3)*D             STL10980
      DEMA(1,3)=DEMA(1,3)+(RD(1)*T(3)+T(1)*RD(3))*E+RD(2)*D             STL10990
      DEMA(2,1)=DEMA(2,1)+(RD(2)*T(1)+T(2)*RD(1))*E+RD(3)*D             STL11000
      DEMA(2,2)=DEMA(2,2)+ 2.D0*(RD(2)*T(2)-RDRO)*E                     STL11010
      DEMA(2,3)=DEMA(2,3)+(RD(2)*T(3)+T(2)*RD(3))*E-RD(1)*D             STL11020
      DEMA(3,1)=DEMA(3,1)+(RD(3)*T(1)+T(3)*RD(1))*E-RD(2)*D             STL11030
      DEMA(3,2)=DEMA(3,2)+(RD(3)*T(2)+T(3)*RD(2))*E+RD(1)*D             STL11040
      DEMA(3,3)=DEMA(3,3)+ 2.D0*(RD(3)*T(3)-RDRO)*E                     STL11050
      RETURN                                                            STL11060
      END                                                               STL11070

C     ******************          DLGMA          ***********************STL12570
      SUBROUTINE DLGMA(RO,SPMA,ND,FACT,ID)                              STL12580
C                                                                       STL12590
C     GIVEN THE ROTATION VECTOR RO, CALCULATE THE CAPITAL GAMMA MATRIX  STL12600
C     SPMA.                                                             STL12610
C                                                                       STL12620
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)                               STL12630
      REAL*8 A,B,C,ARG,ARG2                                             STL12640
      DIMENSION RO(3),SPMA(ND,3),T(3)                                   STL12650
      DATA TRESH2/1.D-10/                                               STL12660
      ARG2=0.D0                                                         STL12670
      DO 10 I=1,3                                                       STL12680
      T(I)=RO(I)*FACT                                                   STL12690
      ARG2=ARG2+T(I)**2                                                 STL12700
   10 CONTINUE                                                          STL12710
      IF (ARG2.LT.TRESH2) THEN                                          STL12720
         A=(1.D0-ARG2/12.D0*(1-ARG2/30.D0*(1-ARG2/56.D0)))/2.D0         STL12730
         B=(1.D0-ARG2/20.D0*(1-ARG2/42.D0*(1-ARG2/72.D0)))/6.D0         STL12740
      ELSE                                                              STL12750
         ARG=DSQRT(ARG2)                                                STL12760
         A=(1.D0-DCOS(ARG))/ARG2                                        STL12770
         B=(1.D0-DSIN(ARG)/ARG)/ARG2                                    STL12780
      END IF                                                            STL12790
      IF (ID.LT.0) THEN                                                 STL12800
          IF (ARG2.LT.TRESH2) THEN                                      STL12810
             B=(1.D0+ARG2/60.D0*(1.D0+ARG2/420.D0))/12.D0               STL12820
          ELSE                                                          STL12830
             B=(1.D0-.5D0*(1.D0-B*ARG2)/A)/ARG2                         STL12840
          END IF                                                        STL12850
        A=-.5D0                                                         STL12860
      END IF                                                            STL12870
         C=1.D0-B*ARG2                                                  STL12880
      SPMA(1,1)=T(1)*T(1)*B+C                                           STL12890
      SPMA(1,2)=T(1)*T(2)*B-T(3)*A                                      STL12900
      SPMA(1,3)=T(1)*T(3)*B+T(2)*A                                      STL12910
      SPMA(2,1)=T(2)*T(1)*B+T(3)*A                                      STL12920
      SPMA(2,2)=T(2)*T(2)*B+C                                           STL12930
      SPMA(2,3)=T(2)*T(3)*B-T(1)*A                                      STL12940
      SPMA(3,1)=T(3)*T(1)*B-T(2)*A                                      STL12950
      SPMA(3,2)=T(3)*T(2)*B+T(1)*A                                      STL12960
      SPMA(3,3)=T(3)*T(3)*B+C                                           STL12970
      RETURN                                                            STL12980
      END                                                               STL12990

C*************************          DPRMMA          ********************STL01230
C=    COMPILER (LINK=IBJ$)                                              STL01240
      SUBROUTINE DPRMMA(A,NRDA,NRA,NCA,B,NRDB,C,NRDC,NCC,IA,IB)         STL01250
C                                                                       STL01260
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B" E SOMMA              STL01270
C     IL RISULTATO IN "C"                                               STL01280
C     A    = PRIMA MATRICE                                              STL01290
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"                STL01300
C     NRA  = N.RO DI RIGHE UTILIZZATE PER IL PRODOTTO DI "A"            STL01310
C     NCA  = N.RO DI COLONNE DI "A"                                     STL01320
C     B    = SECONDA MATRICE                                            STL01330
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"                STL01340
C     C    = MATRICE PRODOTTO                                           STL01350
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"                STL01360
C     NCB  = N.RO DI COLONNE DI "C"                                     STL01370
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO                     STL01380
C                                            1 = SI                     STL01390
C     IB   = INDICE DI TRASPOSIZIONE DI "B"                             STL01400
C                                                                       STL01410
      DOUBLE PRECISION A,B,C                                            STL01420
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1)                           STL01430
      IGO=2*IA+IB+1                                                     STL01440
      GOTO(10,20,30,40),IGO                                             STL01450
10    DO 11 I=1,NRA                                                     STL01460
      DO 11 K=1,NCC                                                     STL01470
      DO 11 L=1,NCA                                                     STL01480
      C(I,K)=C(I,K)+A(I,L)*B(L,K)                                       STL01490
11    CONTINUE                                                          STL01500
      RETURN                                                            STL01510
20    DO 12 I=1,NRA                                                     STL01520
      DO 12 K=1,NCC                                                     STL01530
      DO 12 L=1,NCA                                                     STL01540
      C(I,K)=C(I,K)+A(I,L)*B(K,L)                                       STL01550
12    CONTINUE                                                          STL01560
      RETURN                                                            STL01570
30    DO 13 I=1,NCA                                                     STL01580
      DO 13 K=1,NCC                                                     STL01590
      DO 13 L=1,NRA                                                     STL01600
      C(I,K)=C(I,K)+A(L,I)*B(L,K)                                       STL01610
13    CONTINUE                                                          STL01620
      RETURN                                                            STL01630
40    DO 14 I=1,NCA                                                     STL01640
      DO 14 K=1,NCC                                                     STL01650
      DO 14 L=1,NRA                                                     STL01660
      C(I,K)=C(I,K)+A(L,I)*B(K,L)                                       STL01670
14    CONTINUE                                                          STL01680
      RETURN                                                            STL01690
      END                                                               STL01700

C*************************          DPROMM          ********************STL01970
C=    COMPILER (LINK=IBJ$)                                              STL01980
      SUBROUTINE DPROMM(A,NRDA,NRA,NCA,B,NRDB,C,NRDC,NCC,IA,IB)         STL01990
C                                                                       STL02000
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B"                      STL02010
C     A    = PRIMA MATRICE                                              STL02020
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"                STL02030
C     NRA  = N.RO DI RIGHE UTILIZZATE PER IL PRODOTTO DI "A"            STL02040
C     NCA  = N.RO DI COLONNE DI "A"                                     STL02050
C     B    = SECONDA MATRICE                                            STL02060
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"                STL02070
C     C    = MATRICE PRODOTTO                                           STL02080
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"                STL02090
C     NCB  = N.RO DI COLONNE DI "C"                                     STL02100
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO                     STL02110
C                                            1 = SI                     STL02120
C     IB   = INDICE DI TRASPOSIZIONE DI "B"                             STL02130
C                                                                       STL02140
      DOUBLE PRECISION A,B,C                                            STL02150
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1)                           STL02160
      IGO=2*IA+IB+1                                                     STL02170
      GOTO(10,20,30,40),IGO                                             STL02180
10    DO 11 I=1,NRA                                                     STL02190
      DO 11 K=1,NCC                                                     STL02200
      C(I,K)=0.                                                         STL02210
      DO 11 L=1,NCA                                                     STL02220
      C(I,K)=C(I,K)+A(I,L)*B(L,K)                                       STL02230
11    CONTINUE                                                          STL02240
      RETURN                                                            STL02250
20    DO 12 I=1,NRA                                                     STL02260
      DO 12 K=1,NCC                                                     STL02270
      C(I,K)=0.                                                         STL02280
      DO 12 L=1,NCA                                                     STL02290
      C(I,K)=C(I,K)+A(I,L)*B(K,L)                                       STL02300
12    CONTINUE                                                          STL02310
      RETURN                                                            STL02320
30    DO 13 I=1,NCA                                                     STL02330
      DO 13 K=1,NCC                                                     STL02340
      C(I,K)=0.                                                         STL02350
      DO 13 L=1,NRA                                                     STL02360
      C(I,K)=C(I,K)+A(L,I)*B(L,K)                                       STL02370
13    CONTINUE                                                          STL02380
      RETURN                                                            STL02390
40    DO 14 I=1,NCA                                                     STL02400
      DO 14 K=1,NCC                                                     STL02410
      C(I,K)=0.                                                         STL02420
      DO 14 L=1,NRA                                                     STL02430
      C(I,K)=C(I,K)+A(L,I)*B(K,L)                                       STL02440
14    CONTINUE                                                          STL02450
      RETURN                                                            STL02460
      END                                                               STL02470


C*************************          DPROMV          ********************STL02480
C=    COMPILER (LINK=IBJ$)                                              STL02490
      SUBROUTINE DPROMV(A,NRDA,NRA,NCA,B,C,IA)                          STL02500
C                                                                       STL02510
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER IL VETTORE "B"           STL02520
C                                                                       STL02530
      DOUBLE PRECISION A,B,C,                                           STL02540
     *                 D                                                STL02550
      DIMENSION A(NRDA,1),B(NCA),C(NRA)                                 STL02560
      IGO=IA+1                                                          STL02570
      GOTO(100,200),IGO                                                 STL02580
100   DO 20 I=1,NRA                                                     STL02590
      D=0.                                                              STL02600
      DO 10 K=1,NCA                                                     STL02610
10    D=D+A(I,K)*B(K)                                                   STL02620
20    C(I)=D                                                            STL02630
      RETURN                                                            STL02640
200   DO 40 I=1,NRA                                                     STL02650
      D=0.                                                              STL02660
      DO 30 K=1,NCA                                                     STL02670
30    D=D+A(K,I)*B(K)                                                   STL02680
40    C(I)=D                                                            STL02690
      RETURN                                                            STL02700
      END                                                               STL02710
 
C*************************          DPRVAD          ********************STL02890
C=    COMPILER (LINK=IBJ$)                                              STL02900
      SUBROUTINE DPRVAD(X,Y,Z)                                          STL02910
C                                                                       STL02920
C     ESEGUE IL PRODOTTO DEL VETTORE "X"                                STL02930
C     PER IL VETTORE "Y" E SOMMA IL RISULTATO IN "Z"                    STL02940
C                                                                       STL02950
      DOUBLE PRECISION X,Y,Z,W                                          STL02960
      DIMENSION X(3),Y(3),Z(3),W(3)                                     STL02970
      W(1)= X(2)*Y(3)-X(3)*Y(2)                                         STL02980
      W(2)= X(3)*Y(1)-X(1)*Y(3)                                         STL02990
      W(3)= X(1)*Y(2)-X(2)*Y(1)                                         STL03000
      Z(1)=Z(1)+W(1)                                                    STL03010
      Z(2)=Z(2)+W(2)                                                    STL03020
      Z(3)=Z(3)+W(3)                                                    STL03030
      RETURN                                                            STL03040
      END                                                               STL03050

C*************************          DPRVVA          ********************STL03060
C=    COMPILER (LINK=IBJ$)                                              STL03070
      SUBROUTINE DPRVVA(A,B,C,NRDC,NRA,NRB)                             STL03080
      DOUBLE PRECISION A(1),B(1),C(NRDC,1)                              STL03090
C                                                                       STL03100
C     ESEGUE IL PRODOTTO DEL VETTORE A PER IL VETTORE B TRASPOSTO       STL03110
C     LA MATRICE RISULTANTE VIENE SOMMATA IN C                          STL03120
C                                                                       STL03130
      DO 20 I=1,NRA                                                     STL03140
      DO 20 K=1,NRB                                                     STL03150
      C(I,K)=C(I,K)+A(I)*B(K)                                           STL03160
 20   CONTINUE                                                          STL03170
C                                                                       STL03180
      RETURN                                                            STL03190
C                                                                       STL03200
C                                                                       STL03210
      END                                                               STL03220

C*************************          DRMOVE          ********************STL03230
C=    COMPILER (LINK=IBJ$)                                              STL03240
      SUBROUTINE DRMOVE(A,B,N,H)                                        STL03250
C                                                                       STL03260
C     SOMMA AI TERMINI DEL VETTORE "A" QUELLI DI "B" MOLTIPLICATI PER H STL03270
C                                                                       STL03280
      DOUBLE PRECISION A,B,H                                            STL03290
      DIMENSION A(1),B(1)                                               STL03300
      DO 10 I=1,N                                                       STL03310
      A(I)=A(I)+B(I)*H                                                  STL03320
10    CONTINUE                                                          STL03330
      RETURN                                                            STL03340
      END                                                               STL03350

C*************************          DRUFCT          ********************STL03690
C=    COMPILER (LINK=IBJ$)                                              STL03700
      SUBROUTINE DRUFCT(A,PERM,N,NR,INDER)                              STL03710
      IMPLICIT REAL*8 (A-H,O-Z)                                         STL03720
      DIMENSION A(NR,1),PERM(1)                                         STL03730
      DO10I=1,N                                                         STL03740
      X=0.                                                              STL03750
      DO20K=1,N                                                         STL03760
      IF(DABS(A(I,K)).LT.X)GOTO20                                       STL03770
      X=DABS(A(I,K))                                                    STL03780
20    CONTINUE                                                          STL03790
      IF(X.EQ.0.)GOTO110                                                STL03800
      PERM(I)=1./X                                                      STL03810
10    CONTINUE                                                          STL03820
      DO100 I=1,N                                                       STL03830
      IM1=I-1                                                           STL03840
      IP1=I+1                                                           STL03850
      IPVT=I                                                            STL03860
      X=0.                                                              STL03870
      DO50K=I,N                                                         STL03880
      DP=A(K,I)                                                         STL03890
      IF(I.EQ.1)GOTO40                                                  STL03900
      DO30J=1,IM1                                                       STL03910
      DP=DP-A(K,J)*A(J,I)                                               STL03920
30    CONTINUE                                                          STL03930
      A(K,I)=DP                                                         STL03940
40    IF(X.GT.(DABS(DP)*PERM(K)))GOTO50                                 STL03950
      IPVT=K                                                            STL03960
      X=DABS(DP)*PERM(K)                                                STL03970
50    CONTINUE                                                          STL03980
      IF(X.LE.0.)GOTO110                                                STL03990
      IF(IPVT.EQ.I)GOTO70                                               STL04000
      DO60J=1,N                                                         STL04010
      X=A(IPVT,J)                                                       STL04020
      A(IPVT,J)=A(I,J)                                                  STL04030
      A(I,J)=X                                                          STL04040
60    CONTINUE                                                          STL04050
      PERM(IPVT)=PERM(I)                                                STL04060
70    PERM(I)=IPVT                                                      STL04070
      IF(I.EQ.N)GOTO100                                                 STL04080
      X=A(I,I)                                                          STL04090
      DO90K=IP1,N                                                       STL04100
      A(K,I)=A(K,I)/X                                                   STL04110
      IF(I.EQ.1)GOTO90                                                  STL04120
      DP=A(I,K)                                                         STL04130
      DO80J=1,IM1                                                       STL04140
      DP=DP-A(I,J)*A(J,K)                                               STL04150
80    CONTINUE                                                          STL04160
      A(I,K)=DP                                                         STL04170
90    CONTINUE                                                          STL04180
100   CONTINUE                                                          STL04190
      INDER=0                                                           STL04200
      RETURN                                                            STL04210
110   INDER=I                                                           STL04220
      RETURN                                                            STL04230
      END                                                               STL04240

C*************************          DRUSOL          ********************STL04250
C=    COMPILER (LINK=IBJ$)                                              STL04260
      SUBROUTINE DRUSOL(A,B,NR,N,PERM)                                  STL04270
      IMPLICIT REAL*8 (A-H,O-Z)                                         STL04280
      DIMENSION A(NR,1),B(1),PERM(1)                                    STL04290
      IF(N.GT.1)GOTO 2                                                  STL04300
      B(1)=B(1)/A(1,1)                                                  STL04310
      RETURN                                                            STL04320
2     DO10I=1,N                                                         STL04330
      K=PERM(I)                                                         STL04340
      IF(K.EQ.I)GOTO10                                                  STL04350
      X=B(K)                                                            STL04360
      B(K)=B(I)                                                         STL04370
      B(I)=X                                                            STL04380
10    CONTINUE                                                          STL04390
      DO20I=2,N                                                         STL04400
      IM1=I-1                                                           STL04410
      DP=B(I)                                                           STL04420
      DO40K=1,IM1                                                       STL04430
      DP=DP-A(I,K)*B(K)                                                 STL04440
40    CONTINUE                                                          STL04450
      B(I)=DP                                                           STL04460
20    CONTINUE                                                          STL04470
      B(N)=B(N)/A(N,N)                                                  STL04480
      DO60I=2,N                                                         STL04490
      IM1=N-I+1                                                         STL04500
      INF=IM1+1                                                         STL04510
      DP=B(IM1)                                                         STL04520
      DO70K=INF,N                                                       STL04530
      DP=DP-A(IM1,K)*B(K)                                               STL04540
70    CONTINUE                                                          STL04550
      B(IM1)=DP/A(IM1,IM1)                                              STL04560
60    CONTINUE                                                          STL04570
      RETURN                                                            STL04580
      END                                                               STL04590


C*************************          DSCAL          *********************STL04600
C=    COMPILER (LINK=IBJ$)                                              STL04610
      DOUBLE PRECISION FUNCTION DSCAL(X,Y,N)                            STL04620
C                                                                       STL04630
C     ESEGUE IL PRODOTTO SCALARE DI X PER Y                             STL04640
C                                                                       STL04650
      IMPLICIT REAL*8 (A-H,O-Z)                                         STL04660
      DIMENSION X(1),Y(1)                                               STL04670
      DSCAL=0.                                                          STL04680
      DO 1 I=1,N                                                        STL04690
1     DSCAL=DSCAL+X(I)*Y(I)                                             STL04700
      RETURN                                                            STL04710
      END                                                               STL04720

C*************************          DVET          **********************STL05420
C=    COMPILER (LINK=IBJ$)                                              STL05430
      SUBROUTINE DVET(W,WV,IVER)                                        STL05440
C                                                                       STL05450
C     CALCOLA L'OPERATORE VETTORE DI UN VETTORE                         STL05460
C                                                                       STL05470
      DOUBLE PRECISION W,WV                                             STL05480
      DIMENSION W(3),WV(9)                                              STL05490
      H=DFLOAT(IVER)                                                    STL05500
      WV(1)= 0.                                                         STL05510
      WV(2)= H*W(3)                                                     STL05520
      WV(3)=-H*W(2)                                                     STL05530
      WV(4)= -WV(2)                                                     STL05540
      WV(5)= 0.                                                         STL05550
      WV(6)= H*W(1)                                                     STL05560
      WV(7)= -WV(3)                                                     STL05570
      WV(8)= -WV(6)                                                     STL05580
      WV(9)= 0.                                                         STL05590
      RETURN                                                            STL05600
      END                                                               STL05610


C*************************          LININT          ********************STL08070
C=    COMPILER (LINK=IBJ$)                                              STL08080
      SUBROUTINE LININT(XY,N,X,Y,PEND)                                  STL08090
      DOUBLE PRECISION XY(2,1),X,Y,PEND                                 STL08100
      IF(X.GT.XY(1,1)) GO TO 10                                         STL08110
      K1=1                                                              STL08120
      GO TO 30                                                          STL08130
  10  DO 20 I=2,N                                                       STL08140
      K1=I-1                                                            STL08150
      IF(XY(1,I).GT.X) GO TO 30                                         STL08160
  20  CONTINUE                                                          STL08170
      K1=N-1                                                            STL08180
  30  K2=K1+1                                                           STL08190
      PEND=(XY(2,K2)-XY(2,K1))/(XY(1,K2)-XY(1,K1))                      STL08200
      Y=XY(2,K1)+PEND*(X-XY(1,K1))                                      STL08210
      RETURN                                                            STL08220
      END                                                               STL08230

C*************************          MOVE          **********************STL08370
C=    COMPILER (LINK=IBJ$)                                              STL08380
      SUBROUTINE MOVE(A,B,N)                                            STL08390
C                                                                       STL08400
C     TRASFERISCE N TERMINI DEL VETTORE 'B' NEL VETTORE 'A'             STL08410
C                                                                       STL08420
      DOUBLE PRECISION A,B                                              STL08430
      DIMENSION A(1),B(1)                                               STL08440
      DO 10 I=1,N                                                       STL08450
      A(I)=B(I)                                                         STL08460
10    CONTINUE                                                          STL08470
      RETURN                                                            STL08480
      END                                                               STL08490


C*************************          PRM3            ********************STL08500
      SUBROUTINE PRM3  (A,NRDA,B,NRDB,C,NRDC,IA,IB)                     STL08510
C                                                                       STL08520
C     ESEGUE IL PRODOTTO DELLA MATRICE "A" PER "B"                      STL08530
C     IL PRODOTTO E' DI 3 RIGHE PER 3 COLONNE                           STL08540
C     A    = PRIMA MATRICE                                              STL08550
C     NRDA = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "A"                STL08560
C     NCA  = N.RO DI COLONNE DI "A"                                     STL08570
C     B    = SECONDA MATRICE                                            STL08580
C     NRDB = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "B"                STL08590
C     C    = MATRICE PRODOTTO                                           STL08600
C     NRDC = N.RO DI RIGHE PER IL DIMENSIONAMENTO DI "C"                STL08610
C     IA   = INDICE DI TRASPOSIZIONE DI "A"  0 = NO                     STL08620
C                                            1 = SI                     STL08630
C     IB   = INDICE DI TRASPOSIZIONE DI "B"                             STL08640
C                                                                       STL08650
      DOUBLE PRECISION A,B,C,D                                          STL08660
      DIMENSION A(NRDA,1),B(NRDB,1),C(NRDC,1),D(3,3)                    STL08670
      NRA=3                                                             STL08680
      NCA=3                                                             STL08690
      NCC=3                                                             STL08700
      IGO=2*IA+IB+1                                                     STL08710
      GOTO(10,20,30,40),IGO                                             STL08720
10    DO 11 I=1,NRA                                                     STL08730
      DO 11 K=1,NCC                                                     STL08740
      D(I,K)=0.                                                         STL08750
      DO 11 L=1,NCA                                                     STL08760
      D(I,K)=D(I,K)+A(I,L)*B(L,K)                                       STL08770
11    CONTINUE                                                          STL08780
      GO TO 50                                                          STL08790
20    DO 12 I=1,NRA                                                     STL08800
      DO 12 K=1,NCC                                                     STL08810
      D(I,K)=0.                                                         STL08820
      DO 12 L=1,NCA                                                     STL08830
      D(I,K)=D(I,K)+A(I,L)*B(K,L)                                       STL08840
12    CONTINUE                                                          STL08850
      GO TO 50                                                          STL08860
30    DO 13 I=1,NCA                                                     STL08870
      DO 13 K=1,NCC                                                     STL08880
      D(I,K)=0.                                                         STL08890
      DO 13 L=1,NRA                                                     STL08900
      D(I,K)=D(I,K)+A(L,I)*B(L,K)                                       STL08910
13    CONTINUE                                                          STL08920
      GO TO 50                                                          STL08930
40    DO 14 I=1,NCA                                                     STL08940
      DO 14 K=1,NCC                                                     STL08950
      D(I,K)=0.                                                         STL08960
      DO 14 L=1,NRA                                                     STL08970
      D(I,K)=D(I,K)+A(L,I)*B(K,L)                                       STL08980
14    CONTINUE                                                          STL08990
   50 DO 60 I=1,3                                                       STL09000
      DO 60 K=1,3                                                       STL09010
      C(I,K)=D(I,K)                                                     STL09020
   60 CONTINUE                                                          STL09030
      RETURN                                                            STL09040
      END                                                               STL09050



C*************************          DZERO            ********************
      SUBROUTINE DZERO(A,n)                     
      
      integer*4 n,i
      real*8 A(n)
      
      do i=1,n
        A(i)=0D0
      end do
      return
      end  



 

