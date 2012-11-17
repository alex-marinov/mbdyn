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

      SUBROUTINE FEMGEN(MODNAM)
C
C     this routine reads OUTPUT4 and OUTPUT2 binary files, and generates
C     a "model.fem" ascii file that can be used as an input for the modal
C     element of MBDyn.
C    
C     Author: Giuseppe Quaranta: <quaranta@aero.polimi.it>
C     Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
C
C
      integer maxnod, maxmod
      parameter(maxnod=20000,maxmod=1000)
      character*72 modnam, cname*8,filein*9
      dimension lab(2),t(7),nam(2)
      integer outfil,in2f,in4f,sysout
      integer ncol, nrow, nmodes, err, iret 
      double precision db(6*maxnod)
      real b(4*maxnod)
      integer ib(4*maxnod)
      equivalence (b(1), ib(1))
      double precision m(maxmod,maxmod)
    
      integer j,i 
      integer nword,irtn

      data sysout /6/
      data outfil /10/
      data in2f   /11/
      data in4f   /12/

      print *, 'FEMGEN: translation of Binary NASTRAN modal data'
      print *, '        for use in the modal MBDyn element'
      print *, ' '
      print *, 'Files "mbdyn.mat" and "mbdyn.tab" must be present'
      print *, 'Maximum number of modes: ', maxmod
      print *, 'Maximum number of nodes: ', maxnod
      print *, 'Recompile with larger maxmod, maxnod as needed'
      print *, ' '   

      if (modnam(1:1).EQ.' ') then
        print*, 'Please enter the model name'
        read '(A72)', modnam          
      endif

      print*, 'Output to file ', modnam

C     Cancella gli spazi dopo il nome per aggiungere l'estensione .fem 
C     al file di uscita       
      do j=72,1,-1
         if (modnam(j:j).ne.' ') goto 6
      enddo     
 6    open(outfil,file=modnam(1:j)//'.fem')
C     apre i due file che contengono i dati di ingresso
C     mbdyn.mat contiene le matrici (OUTPUT4)
C     MHH KHH LUMPED-MASS 
      filein ='mbdyn.mat'
      open(in4f, file=filein, status='old',form='unformatted',
     *     err=900, iostat=iret)
C     mbdyn.tab contine le tabelle (OUTPUT2)
C     GPL BGPDT OUGV1
      filein = 'mbdyn.tab'
      open(in2f, file=filein, status='old',form='unformatted',
     *     err=900, iostat=iret)
      
C     legge la prima matrice in modo da conoscere il numero di modi      
      call GETIDS(m,ncol,nrow,cname,maxmod,maxmod,m,1,in4f,err)
      nmodes = nrow
      
      write(outfil,'(A24)') '** MBDyn MODAL DATA FILE'
      write(outfil,'(A18)') '** NODE SET "ALL" '
      write(outfil,*) ' '
      write(outfil,*) ' '

C     legge il primo record di GPL che contiene la lista dei nodi
      call IOPEN(in2f,sysout,lab)
      call IHEADR(in2f,sysout,nam,t)  
      write(sysout,'(A15,A4,A4)') ' Reading Table ', nam(1), nam(2) 
C     legge il trailer del blocco che indica il numero di elementi b(1) 
      call IREAD(in2f,sysout,ib,maxnod,0,nword,irtn)
      if (irtn .ne. 1) then
         print*, 'Errors in GPL table trailer. Aborting!'
         goto 902
      endif
      
      write(outfil,'(A25)') '** RECORD GROUP 1, HEADER'
      write(outfil,'(A31, A39)') '**   REVISION,  NODE,  NORMAL, ',
     *     'ATTACHMENT, CONSTRAINT, REJECTED MODES.'
      write(outfil,'(6X, A4, 5I10)') 'REV01', nword, nmodes, 0, 0, 0

      write(outfil,'(A2)') '**'
      write(outfil,'(A34,A9)') '** RECORD GROUP 2, FINITE ELEMENT ',
     *     'NODE LIST'
      do j=0,nword,6
	 if (j+6 .gt. nword) then
	    write(outfil,'(6I10)') (ib(i), i=j+1,nword)
         else
            write(outfil,'(6I10)') (ib(j+i), i=1,6) 
         endif
      enddo
C     il secondo record DI GRT non serve 
      irtn = 0
      do while (irtn .eq. 0)
         call IREAD(in2f,sysout,ib,maxnod,0,nword,irtn)
      enddo 
      irtn = 0
      do while (irtn .eq. 0)
         call IREAD(in2f,sysout,ib,maxnod,0,nword,irtn)
      enddo 

      write(outfil,'(A2)') '**'
      write(outfil,'(A33, A13)') '** RECORD GROUP 3, INITIAL MODAL ',
     *     'DISPLACEMENTS'
      write(outfil,'(500(1X,1PE17.10))') (0.0*i, i=1,nmodes)        
      write(outfil,'(A2)') '**'
      write(outfil,'(A32,A10)') '** RECORD GROUP 4, INITIAL MODAL ', 
     *     'VELOCITIES'
      write(outfil,'(500(1X,1PE17.10))') (0.0*i, i=1,nmodes)        
      
C     legge la matrice BGPDT che ha le posizioni dei nodi 
      call IHEADR(in2f,sysout,nam,t) 
      write(sysout,'(A15,A4,A4)') ' Reading Table ', nam(1), nam(2) 
C     legge il trailer del blocco  
C      call IREAD(in2f,sysout,ib,maxnod,0,nword,irtn)
C     if (irtn .ne. 1) then
C        print*, 'Errors in BGPDT table trailer. Aborting!'
C        goto 902
C     endif
C     legge le posizioni che possono essere al piu' 4000 da errore 
      call IREAD(in2f,sysout,b,4*maxnod,0,nword,irtn)
      if (irtn .ne. 1) then
         print*, 'Node saved by BGPDT are more than ',maxnod, 
     *        '. You must recompile with bigger arrays !!'
         goto 902
      endif
      write(outfil,'(A2)') '**'
      write(outfil,'(A38)') '** RECORD GROUP 5, NODAL X COORDINATES'   
      do j=2,nword,4
         write(outfil,'(E17.10)') b(j) 
      enddo
      write(outfil,'(A2)') '**'
      write(outfil,'(A38)') '** RECORD GROUP 6, NODAL Y COORDINATES'   
      do j=3,nword,4
         write(outfil,'(E17.10)') b(j) 
      enddo
      write(outfil,'(A2)') '**'
      write(outfil,'(A38)') '** RECORD GROUP 7, NODAL Z COORDINATES'   
      do j=4,nword,4
         write(outfil,'(E17.10)') b(j) 
      enddo

C     il secondo record DI BGPDT non serve 
      irtn = 0

      do while (irtn .eq. 0) 
         call IREAD(in2f,sysout,b,maxnod,0,nword,irtn)
      enddo 
      
C     legge la tabella delle forme dei modi OUGV1
      write(outfil,'(A2)') '**'
      write(outfil,'(A30)') '** RECORD GROUP 8, MODE SHAPES'    
      call tabrd(in2f, outfil, outfil,7)

C     stampa la matrice di massa
      write(outfil,'(A2)') '**'
      write(outfil,'(A36)') '** RECORD GROUP 9, MODAL MASS MATRIX'  
      do  j = 1,nrow 
         write(outfil,'(500(1X,1PE17.10))')(m(j,i),i=1,ncol) 	  
      enddo
C     legge la matrice di rigidezza     
      call GETIDS(m,ncol,nrow,cname,maxmod,maxmod,m,1,in4f,err)
C     stampa la matrice di rigidezza
      write(outfil,'(A2)') '**'
      write(outfil,'(A42)') '** RECORD GROUP 10, MODAL STIFFNESS MATRIX'  
      do  j = 1,nrow 
         write(outfil,'(500(1X,1PE17.10))')(m(j,i),i=1,ncol) 	  
      enddo      
C     legge il vettore lumped mass      
      call GETIDS(b,ncol,nrow,cname,6*maxnod,1,db,1,in4f,err)
C     stampa il vettore lumped mass
      write(outfil,'(A2)') '**'
      write(outfil,'(A38,A12)') '** RECORD GROUP 11, DIAGONAL OF LUMPED'
     *     ,' MASS MATRIX'
      do  j = 0,nrow,6 
         write(outfil,'(500(1X,1PE17.10))')(db(j+i),i=1,6) 	  
      enddo  
      
      stop
C     Errors      
 900  write(sysout,901) filein,iret
 901  format(/,'FEMGEN error while opening input file ', A9, 
     *     ',IRET= ',I4)
 902  stop 
      end 

C******************************************************************************
      SUBROUTINE GETIDS(B,NCOL,NROW,CNAME,NRB,NCB,DB,IDL,IU,ERR)
C     legge una matrice di double precision
C     Subroutine taken from mattest.f NASTRAN Version 70.7 utility program 
C******************************************************************************
C     
C     PURPOSE:
C     THIS ROUTINE INTERPRETS A MATRIX WRITTEN BY OUTPUT4
C     IT WILL HANDLE SINGLE OR DOUBLE PRECISION
C     IT WILL PROCESS EITHER THE SPARSE OR NON-SPARSE FORMS
C     
C     ARGUMENTS:
C     
C     B    -  SINGLE PRECISION ARRAY SUBSCRIPTED BY B(NRB,NCB)
C     THE MATRIX WILL BE RETURNED IN B IF IDL = 0
C     NCOL -  RETURNED AS THE ACTUAL NUMBER OF COLUMNS OF MATRIX
C     NROW -  RETURNED AS THE ACTUAL NUMBER OF ROWS OF MATRIX
C     CNAME -  RETURNED AS THE NAME OF THE INCOMMING MATRIX
C     NRB  -  MAXIMUM ROWS IN B AND DB ARRAYS
C     NCB  -  THE MAXIMUM NUMBER OF COLUMNS IN B AND DB ARRAYS
C     DB   -  DOUBLE PRECISION ARRAY LIKE B
C     IDL  -  0 ASKS FOR SINGLE RETURN  1  ASKS FOR D.P.
C     IU   -  FORTRAN UNIT TO READ
C     
      INTEGER SYSOUT, ERR
      DOUBLE PRECISION DB(NRB,NCB)
      DOUBLE PRECISION DA,DD(1)
      REAL B(NRB,NCB)
      CHARACTER*8 CNAME
      REAL DR(2)
C     
C     A,DA AND IA REPRESENT THE LARGEST COLUMN OF A MATRIX WHICH CAN BE
C     READ
C     
      DIMENSION A(360000),DA(180000),IA(360000)
      EQUIVALENCE (A(1),DA(1),IA(1))
      EQUIVALENCE (DD(1),DR(1))
C     
      DATA LENA/360000/
      DATA SYSOUT /6/
      
C     
C     READ MATRIX DESCRIPTORS
C     
      ERR = 0
      READ (IU) NCOL,NROW,NF,NTYPE,CNAME
      PRINT*, 'Reading Matrix ', CNAME,NCOL,NROW
C     
C     CHECK IF MATRIX IS TOO LARGE
C     
      IF(NCOL .GT. NCB .OR. NROW .GT. NRB) GO TO 50
      IF(NTYPE*NROW .GT. LENA) GO TO 50
      IF(NTYPE .GT. 2) GO TO 200
C     
C     ZERO B MATRIX
C     
      DO 7 I = 1,NCOL
         DO 66 J = 1,NROW
            IF(IDL .EQ. 1)  GO TO 67
            B(J,I) = 0.0
            GO TO 66
 67         DB(J,I) = 0.0D0
 66      CONTINUE
    7 CONTINUE
C     
C     FOR EACH COLUMN (ONLY NON-ZERO COLUMNS ON FILE)
C     
      DO 10 I = 1,NCOL
         READ(IU) ICOL,IROW,NW,(A(K),K=1,NW)
         IF(ICOL .GT. NCOL) GO TO 20
C     
C     TEST FOR SPARSE MATRIX OPTION
C     
         IF(IROW .EQ. 0) GO TO 30
C     
C     DENSE MATRIX FORMAT
C     
         IF(NTYPE .EQ. 2) GO TO 100
         DO 5 J=1,NW
            K = IROW+J-1
            IF(IDL .EQ. 1) GO TO 8
            B(K,ICOL) = A(J)
            GO TO 5
 8          DB(K,ICOL) = A(J)
 5       CONTINUE
         GO TO 10
C     
C     DOUBLE INCOMMING MATRIX
C     
 100     NW = NW/2
         DO 6 J = 1,NW
            K = IROW+J-1
            IF(IDL .EQ. 1) GO TO 11
            B(K,ICOL) = DA(J)
            GO TO 6
 11         DB(K,ICOL) = DA(J)
 6       CONTINUE
         GO TO 10
C     
C     SPARSE INCOMMING MATRIX
C     
 30      CONTINUE
         NTR = 1
 32      L = IA(NTR)/65536
         IROW = IA(NTR)-L*65536
         NTW = L-1
         IF(NTYPE .EQ. 2) GO TO 40
C     
C     SINGLE INCOMMING MATRIX
C     
         DO 31 J = 1,NTW
            K = IROW+J-1
            IF(IDL .EQ. 1) GO TO 34
            B(K,ICOL) = A(NTR+J)
            GO TO 31
 34         DB(K,ICOL) = A(NTR+J)
 31      CONTINUE
 33      NTR = NTR+L
         IF(NTR .GE. NW) GO TO 10
         GO TO 32
C     
C     DOUBLE INCOMMING MATRIX
C     
 40      CONTINUE
         DO 41 J = 1,NTW,2
            K = IROW+J/2
            DR(1) = A(NTR+J)
            DR(2) = A(NTR+J+1)
            IF(IDL .EQ. 1) GO TO 42
            B(K,ICOL) = DD(1)
            GO TO 41
 42         DB(K,ICOL) = DD(1)
 41      CONTINUE
         GO TO 33
 10   CONTINUE
C     
C     THERE IS DUMMY RECORD FOR LAST COLUMN +1
C     
      READ(IU) ICOL,IROW,NW,(A(K),K=1,NW)
 20   RETURN
 50   CONTINUE
      WRITE(SYSOUT,55)
 55   FORMAT(' MSC/NASTRAN MATRIX TOO LARGE.')
      STOP
 200  WRITE(SYSOUT,*)('Matrix is not real. Trying complex d.p....')
      ERR = 1
      REWIND(IU)
      RETURN 
      END
      
C***********************************************************************
C     
      SUBROUTINE TABRD(IUN,IOUT,ITOUT,TTYP)
C     legge alcune tabelle standard di NASTRAN
C     Subroutine taken from tabtest.f NASTRAN Version 70.7 utility program 
C     modified by Giuseppe Quaranta <quaranta@aero.polimi.it>
C     July 2001
C***********************************************************************
C     
C     Reads OUTPUT2 TABLES output
C     IUN    -  input file  unit
C     IOUT   -  output file unit
C     ITOUT  -  output table unit 
C     TTYP   -  table type:
C     7 modal displacements
C     3 modal spc forces
C     5 modal element stresses
C     4 modal element forces  
C     
      REAL BLOCK(200000)
      INTEGER IBLOCK(200000)
      EQUIVALENCE (BLOCK(1),IBLOCK(1))
C     
C     DATA BLOCK NAME (NAM) IS 2 WORDS, NASTRAN TRAILER (T) IS 7 WORDS,
C     AND THE OUTPUT2 LABEL (LAB) IS 2 WORDS.
C     
      INTEGER TTYP, SYSOUT
      INTEGER IRTN,ITYPE,NWDS,G,H,D,DM1,IEID,IGID,ISA,ISB
      INTEGER IMODE, COUNT, STARTS, DELTAS, STARTF, DELTAF, GRIDID
      INTEGER  I, K, J, IFLAG 
      DIMENSION NAM(2),T(7),LAB(2)
      CHARACTER*6 ELETYP(100)
      DIMENSION STARTS(100), DELTAS(100), STARTF(100), DELTAF(100)
      DIMENSION GRIDID(100)
      
      DATA SYSOUT /6/ 
      
C     THESE ARE THE ELEMENT NAMES.
      DATA ELETYP/' ROD  ',' BEAM ',' TUBE ','SHEAR ','      ',
     1     'TRIA1 ','TRBSC ','TRPLT ','TRMEM ','CONROD','ELAS1 ',
     2     'ELAS2 ','ELAS3 ','ELAS4 ','QDPLT ','QDMEM ','TRIA2 ',
     3     'QUAD2 ','QUAD1 ','DAMP1 ','DAMP2 ','DAMP3 ','DAMP4 ',
     4     'VISC  ','MASS1 ','MASS2 ','MASS3 ','MASS4 ','CONM1 ',
     5     'CONM2 ','PLOTEL','      ','QUAD4 ',' BAR  ','CONEAX',
     6     'TRIARG','TRAPRG',' GAP  ','TETRA ','WEDGE ','HEXA1 ',
     7     'HEXA2 ','FLUID2','FLUID3','FLUID4','FLMASS','AXIF2 ',
     8     'AXIF3 ','AXIF4 ','SLOT3 ','SLOT4 ',' HBDY ','TRIAX6',
     9     'TRIM6 ',7*'DUMMY ','QDMEM1','QDMEM2','QUAD8 ',' HEX8 ',
     A     'HEX20 ',' HEXA ','PENTA ',' BEND ',2*'      ','AER01 ',
     B     'FTUBE ','TRIA3 ','TRIA6 ',19*'      ','QUAD4 ','QUAD8 ',
     C     'TRIA3 ','TRIA6 ','      ',' BAR  '/
      
C     Starting point for element stress data
      DATA STARTS/  2,2,2,2,0,2,2,2,2,2,2,2,2,0,2,2,2,2,2,0,0,0,0,0,0,
     1     0,0,0,0,0,0,0,2,2,2,2,2,0,5,2,5,5,0,0,0,0,2,2,2,2,
     2     2,0,2,2,2,2,2,2,2,2,2,2,2,3,2,2,5,5,2,0,0,0,0,2,3,
     3     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3,3,3,3,0,2/
C     Number of output data for each element stress
      DATA DELTAS/  5,10,5,4,0,17,17,17,8,5,2,2,2,0,17,8,17,17,17,0,
     1     0,0,0,0,0,0,0,0,0,0,0,0,17,16,18,5,21,0,21,9,
     2     21,21,0,0,0,0,5,10,12,6,7,0,9,9,10,10,10,10,10,10,
     3     10,8,8,17,21,21,21,21,11,0,0,0,0,17,17,0,0,0,0,0,
     4     0,0,0,0,0,0,0,0,0,0,0,0,0,0,11,11,11,11,0,10/ 
C     Starting point for element force data  
      DATA STARTF/ 2,2,2,2,0,2,0,2,0,2,2,2,2,2,2,0,2,2,2,2,2,2,2,0,0,
     1     0,0,0,0,0,0,0,2,2,2,2,2,1,0,0,0,0,0,0,0,0,0,0,0,0,
     2     0,2,0,0,2,2,2,2,2,2,2,0,2,3,0,0,0,0,2,0,0,0,2,2,3,
     3     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,2/
C     Number of output data for each element stress
      DATA DELTAF/  2,9,2,17,0,6,0,6,0,2,2,2,2,2,6,0,6,6,6,2,
     1     2,2,2,0,0,
     2     0,0,0,0,0,0,0,9,9,7,10,13,8,0,0,0,0,0,0,0,
     3     0,0,0,0,0,
     4     0,8,0,0,10,10,10,10,10,10,10,0,17,10,0,0,0,0,7,0,
     5     0,0,2,9,10,
     6     0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     7     0,0,0,0, 8/ 
C     Grid id data(Y/N).
      DATA GRIDID/ 0,1,36*0,1,0,2*1,10*0,2*1,9*0,6*1,5*0,1,25*0/
      

C     READ THE HEADER RECORD AND START TO PROCESS
      CALL IHEADR(IUN,IOUT,NAM,T)
      write(sysout,'(A15,A4,A4)') ' Reading Table ', NAM(1), NAM(2) 

C     GET THE DATA BLOCK
      IRTN = 1
      COUNT = 1
      IMODE = 0              
      NWDS = 0
      ITYPE = 0
      IWR = 0
      IEID = 0
      IFLAG = 0
C     
C     ASK FOR 1000 WORDS TO BE TRANSMITTED PER CALL TO IREAD
C     RETURN UP TO 10000 WORDS IF NORMAL RETURN (IWR = WORDS RETURNED)
C     
C     legge l'header della tabella      
 30   CALL IREAD(IUN,IOUT,BLOCK,1000,0,IWR,IRTN)
      IF (IRTN .EQ. 0)IWR = 10000
      IF (IRTN .GE. 2) GO TO 999
      
C     WHICH DATA BLOCK IS IT?
      
      IF (IBLOCK(2) .EQ. TTYP) THEN
         
C     ... EVALUATING IN THE CONTROL BLOCK
	  IMODE = IMODE + 1
          WRITE(IOUT,'(A26,I2)')'**    NORMAL MODE SHAPE # ', IMODE  
         
C     istruzioni per i blocchi conteneti i dati relativi agli elementi  
         IF ((TTYP .GE. 4).AND.(TTYP .LE. 5)) THEN
            ITYPE = IBLOCK(3)
            IF (IMODE .EQ. 1) THEN
               WRITE(ITOUT,5600) ELETYP(ITYPE), ITYPE
 5600          FORMAT(A,I5)
            ENDIF
         ENDIF
         NWDS = IBLOCK(10)
         
C     ... GET THE NEXT RECORD BLOCK
C     ASK FOR 20000 WORDS TO BE TRANSMITTED PER CALL TO IREAD
C     legge le prime mille parole del record 

 32      CALL IREAD(IUN,IOUT,BLOCK,20000,IFLAG,IWR,IRTN)
C     RETURN UP TO 10000 WORDS IF NORMAL RETURN (IWR = WORDS RETURNED)
         IF (IRTN .EQ. 0) IWR = 20000
C     
C     ISA - THE FIRST OUTPUT FOR THAT ELEMENT AND GRID
C     ISB - THE LAST OUTPUT FOR THE ELEMENT AND GRID
         IF ((TTYP .EQ. 3).OR.(TTYP .EQ. 7)) THEN
            G = 3
            H = NWDS
            D = 7
            DM1 = D - 1
            DO 5100 I=1,IWR,NWDS
               ISB = G + D - 2
               WRITE(IOUT,5630) (BLOCK(K),K=G,ISB)
 5630          FORMAT(6(1X, 1PE17.10))
               COUNT = COUNT + D - 1
               G = G + D + 1
 5100       CONTINUE  
            
         ELSE IF (TTYP .EQ. 4) THEN
            IF (STARTF(ITYPE) .NE. 0) THEN
               G = STARTF(ITYPE)
               H = NWDS
               D = DELTAF(ITYPE)
               DM1 = D - 1
               DO 6100 I=1,IWR,NWDS
                  IEID=IBLOCK(I)/10
                  IF (GRIDID(ITYPE) .EQ. 1) THEN
                     DO 6150 J=G,H,D
                        IGID=IBLOCK(J)
                        ISA=J+1
                        ISB=J+DM1
                        IF ((IMODE .EQ. 1)) THEN 
                           WRITE(ITOUT,6620) IEID, COUNT, D-1, IGID
 6620                      FORMAT(1X, I14, I14, I14, I14)
                        ENDIF                      
                        WRITE(IOUT,6630) (BLOCK(K),K=ISA,ISB)
                        COUNT = COUNT + 1
 6630                   FORMAT(40(1X, 1PE17.10))
 6150                CONTINUE
                     G=G+NWDS
                     H=H+NWDS
                  ELSE 
                     ISB = G + D - 2
                     IF (IMODE .EQ. 1) THEN 
                        WRITE(ITOUT,6640) IEID, COUNT, D-1
 6640                   FORMAT(1X, I14, I14, I14)
                     ENDIF
                     WRITE(IOUT,6630) (BLOCK(K),K=G,ISB)
                     COUNT = COUNT + 1
                     G = G + D
                  ENDIF
 6100          CONTINUE  
            ELSE
               WRITE(SYSOUT,6800) ELETYP(ITYPE), ITYPE
 6800          FORMAT(/,
     1              '  ** ERROR ** IN OEF1 DATA BLOCK, ',A,' ELEMENT',
     2              ' -- TYPE =',I5,' IS NOT SUPPORTED.')
            ENDIF
            
            
         ELSE IF (TTYP .EQ. 5) THEN
            IF (STARTS(ITYPE) .NE. 0) THEN
               G = STARTS(ITYPE)
               H = NWDS
               D = DELTAS(ITYPE)
               DM1 = D - 1
               DO 7100 I=1,IWR,NWDS
                  IEID=IBLOCK(I)/10
                  IF (GRIDID(ITYPE) .EQ. 1) THEN
                     DO 7150 J=G,H,D
                        IGID=IBLOCK(J)
                        ISA=J+1
                        ISB=J+DM1
                        IF (IMODE .EQ. 1) THEN 
                           WRITE(ITOUT,7620) IEID, COUNT, D-1, IGID
 7620                      FORMAT(1X, I14, I14, I14, I14)
                        ENDIF                      
                        WRITE(IOUT,7630) (BLOCK(K),K=ISA,ISB)
                        COUNT = COUNT + 1
 7630                   FORMAT(40(1X, 1PE17.10))
 7150                CONTINUE
                     G=G+NWDS
                     H=H+NWDS
                  ELSE 
                     ISB = G + D - 2
                     IF (IMODE .LE. 1) THEN 
                        WRITE(ITOUT,7640) IEID, COUNT, D-1
 7640                   FORMAT(1X, I14, I14, I14)
                     ENDIF
                     WRITE(IOUT,7630) (BLOCK(K),K=G,ISB)
                     COUNT = COUNT + 1
                     G = G + D
                  ENDIF
 7100          CONTINUE  
            ELSE
               WRITE(SYSOUT,7800) ELETYP(ITYPE), ITYPE
 7800          FORMAT(/,
     1              '  ** ERROR ** IN OES1 DATA BLOCK, ',A,' ELEMENT',
     2              ' -- TYPE =',I5,' IS NOT SUPPORTED.')
            ENDIF
            
         ELSE
            WRITE(IOUT,800) IBLOCK(2)
 800        FORMAT(/,'  ** ERROR ** DATA BLOCK NOT SUPPORTED')
            IRTN = 3
            
         ENDIF
         
      ENDIF
      IF (IRTN .GE. 2) GO TO 999
      IF (IRTN .EQ. 0) GO TO 32
      GO TO 30
 999  CONTINUE
      IF (IRTN .GT. 2) GO TO 900
C     
      RETURN
C     ERRORS
C     
 900  CONTINUE
      WRITE(SYSOUT,910) IRTN
 910  FORMAT(20H1BAD PROGRAM IRTN = ,I14 )
      STOP
      END          
C     ======================================================================
      SUBROUTINE IOPEN(IUN,IOUT,L)
C     ======================================================================
C     
C     THIS ROUTINE OPENS (TO READ) AN INPUTT2 FILE
C     
      IMPLICIT INTEGER (A-Z)
C     
C     IT MUST BE CALLED ONCE AND ONLY ONCE PER UNIT
C     
C     IUN--INPUT--INTEGER--FORTRAN UNIT NUMBER (IE 20)
C     L(2)--OUTPUT--INTEGER--LABEL FROM DMAP (IE MYPROG)
C     
      INTEGER SYSIN,SYSOUT
      DIMENSION T(7),L(2)
      CHARACTER*72 INPT2,OUTPT
      DATA SYSIN /5/
      DATA SYSOUT /6/
      
      REWIND (IUN)
      READ (IUN) KEY
      IF(KEY .NE. 3) GO TO 9000
C     BRING IN DATE
      READ (IUN) M,D,Y
      READ (IUN) KEY
      IF(KEY .NE. 7) GO TO 9000
C     BRING IN TITLE
      READ (IUN) (T(I),I=1,7)
      READ (IUN) KEY
      IF(KEY .NE. 2) GO TO 9000
C     BRING IN LABEL
      READ (IUN) L
      READ (IUN) KEY
      IF(KEY .NE. -1) GO TO 9000
      READ (IUN) KEY
      IF(KEY .NE. 0) GO TO 9000
      RETURN
 9000 WRITE(IOUT,900) KEY
 900  FORMAT(1H1,15HIOPEN BAD KEY =,I14)
      STOP
      END
C     ======================================================================
      SUBROUTINE IHEADR(IUN,IOUT,NAM,T)
C     ======================================================================
C     
C     IHEADR IS CALLED ONCE PER DATA BLOCK ON UNIT IUN
C     
      IMPLICIT INTEGER (A-Z)
C     
C     IUN--INTEGER--INPUT--FORTRAN UNIT NUMBER (IE 20)
C     NAM--BCD--OUTPUT--DATA BLOCK NAME (2 WORDS) (IE EQEXIN)
C     T--OUTPUT--INTEGER--DATA BLOCK CONTROL BLOCK (7 WORDS IE
C     TRAILER(1)=T(2)
C     
      DIMENSION NAM(2),T(7),IH(2)
C     
      COMMON /BUFFER/IN,ICBP,IRL,IZ(200000)
C     
      READ (IUN) KEY
      IF(KEY .NE. 2) GO TO 9000
C     READ IN DATA BLOCK NAME
      READ (IUN) NAM
      READ (IUN) KEY
      IF(KEY .NE. -1) GO TO 9000
      READ (IUN) KEY
      IF (KEY .NE. 7) GO TO 9000
C     READ IN TRAILER
      READ (IUN) T
      READ (IUN) KEY
      IF(KEY .NE. -2) GO TO 9000
      READ (IUN) KEY
      IF(KEY .NE. 1) GO TO 9000
      READ (IUN) KEY
      IF(KEY .NE. 0) GO TO 9000
      READ (IUN) KEY
      IF (KEY .LT. 2) GO TO 9000
C     READ IN (SKIP) DATA FROM DATA BLOCK HEADER
      READ (IUN) IH
      READ (IUN) KEY
      IF(KEY .NE. -3) GO TO 9000
C     
C     START BUFFER ROUTINE A BEGINING
C     
      IN=0
      RETURN
 9000 WRITE(IOUT,900) KEY
 900  FORMAT (1H1,16HIHEADR BAD KEY = , I14)
      STOP
      END
C     ======================================================================
      SUBROUTINE IREAD(IUN,IOUT,ID,NW1,IFLAG,NWR,IRTN)
C     ======================================================================
C     
C     THIS ROUTINE READS IN INPUTT2 FILE WITH USER CONTROLS
C     
      IMPLICIT INTEGER (A-Z)
C     
C     IUN--INPUT--FORTRAN UNIT NUMBER (IE 20)
C     ID--OUTPUT--ARRAY CONTAINING RETURNED WORDS
C     NW1--INPUT--INTEGER--THE NUMBER OF WORDS TO RETURN
C     IF NW1 IS NEGATIVE THIS MANY WORDS WILL BE SKIPPED
C     IFLAG--INPUT--IF IFLAG IS 0 , A CONTINUATION READ WILL OCCUR
C     IF IFLAG IS 1 THE REMAINER OF THE RECORD WILL BE
C     SKIPPED
C     NWR--OUTPUT--IF THE RECORD DBLOCK NOT CONTAIN NW1 REMAINING
C     WORDS THEN NWR WILL INDICATE THE ACTUAL NUMBER
C     RETURNED
C     IRTN--OUTPUT--0 NORMAL RETURN
C     1 END OF RECORD (NWR WILL BE SET)
C     2 END OF DATA BLOCK
C     
      DIMENSION ID(1)
C     
      COMMON /BUFFER/IN,ICBP,IRL,IZ(200000)
C     
C     IN INDICATES PRESENCE IN IZ
C     ICBP IS THE CURRENT BUFFER POINTER
C     IRL IS THE CURRENT RECORD LENGTH
C     IZ CONTAINS THE DATA
C     
C     THIS ROUTINE CAN ONLY HANDLE ONE FORTRAN FILE AT A TIME
C     MORE CODE COULD EXTEND THIS CONCEPT
C     
      NW = NW1
      IRTN= 0
      IUAP=0
      IF(IN .NE. 0) GO TO 70
C     
C     READ INTO BUFFER ON INITIAL READ START OF RECORD FLAG
C     
      READ(IUN) KEY
      IF(KEY .NE. 1) GO TO 9000
C     READ RECORD TYPE (TABLES=0, ONLY TABLES ALLOWED)
      READ (IUN) KEY
      IF(KEY .NE. 0) GO TO 9000
 30   READ (IUN) KEY
C     END-OF-RECORD OR END-OF-FILE OR LENGTH-OF-RECORD INDICATION
      IF( KEY )160,200,40
 40   IF(KEY.LE.200000)GO TO 60
C     OUTPUT2 RECORD IS LARGER THAN 20000 DIMENSION OF IZ ARRAY
      WRITE(IOUT,50)
 50   FORMAT(50H0 ARRAY IZ DIMENSION EXCEEDED IN SUBROUTINE IREAD.)
      KEY = 200000
      WRITE(IOUT,51) KEY
 51   FORMAT(' READING ONLY THE FIRST ',I14)
 60   READ (IUN) (IZ(I),I=1,KEY)
      IRL = KEY
      ICBP = 0
      IN = 1
C     
C     DATA IS IN IZ
C     
 70   CONTINUE
      IF(NW) 180,170,80
C     
C     ARE ENOUGH WORDS IN BUFFER
C     
 80   CONTINUE
      IF(IRL -ICBP .LT. NW) GO TO 140
C     
C     SUPPLY WORDS
C     
      DO 90 I = 1,NW
         ID(IUAP+I) = IZ(ICBP+I)
 90   CONTINUE
      IUAP = IUAP + NW
C     
C     UPDATE BUFFER POINTER
C     
      ICBP=ICBP+NW
 100  IF(IFLAG .NE. 0) GO TO 120
 110  RETURN
C     
C     SKIP TO NEXT LOGICAL RECORD
C     
 120  CONTINUE
C     READ END OF RECORD FLAG IF NOT NEGATIVE IS MORE DATA FOR
C     CURRENT RECORD, READ DATA LENGTH KEY AND NEXT RECORD TO
C     SKIP THAT LENGTH OF DATA, AND LOOK FOR END-OF-RECORD FLAG
      READ (IUN) KEY
      IF(KEY .GT. 0) GO TO 130
      IN = 0
      GO TO 110
C     
C     SKIP MORE
C     
 130  READ(IUN) KEY
      GO TO 120
C     
C     MORE WORDS REQUESTED THAN AVAILABLE
C     
 140  CONTINUE
      K = IRL-ICBP
C     
C     SUPPLY REMAINDER OF CURRENT BUFFER
C     
      DO 150 I = 1,K
         ID(IUAP+I) = IZ(ICBP+I)
 150  CONTINUE
C     
C     BUMP USER POINTER AND DECREMENT WORDS
C     
      IUAP = IUAP +K
      NW = NW-K
C     
C     READ REMAINDER OF DATA RECORD
C     
      GO TO 30
C     
C     END OF RECORD
C     
 160  CONTINUE
      NWR = IUAP
      IRTN = 1
      IN = 0
      GO TO 110
C     
C     USER REQUESTED ZERO WORDS HONOR IFLAG AND QUIT
C     
 170  CONTINUE
      GO TO 100
C     
C     USER REQUESTED WORDS BE SKIPPED ( NW IS NEGATIVE HERE )
C     
 180  CONTINUE
      IF(ICBP-IRL .GT. NW) GO TO 190
C     
C     SKIP WITHIN CURRENT BUFFER
C     
      ICBP = ICBP - NW
      GO TO 100
C     
C     SKIP OVER MORE THAN CURRENT BUFFER
C     
 190  CONTINUE
C     REDUCE THE NEGATIVE COUNT, MAKE IT LESS NEGATIVE
      NW = NW+IRL-ICBP
C     
C     READ REMAINDER OF DATA RECORD
C     
      GO TO 30
C     
C     END OF FILE
C     
 200  CONTINUE
      IRTN = 2
      IN = 0
      GO TO 110
C     
C     ERRORS
C     
 9000 WRITE(IOUT,9001)  KEY,ICBP,IRL,NW,NW1,IFLAG
 9001 FORMAT(1H0,6I14)
      STOP
      END

      





