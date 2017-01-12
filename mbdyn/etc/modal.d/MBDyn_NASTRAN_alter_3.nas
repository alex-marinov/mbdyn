$ $Header$
$ Alters for MSC/NASTRAN v 70.7
$ ALTERS for SOL 103 MODAL ANALYSIS
$
$ MBDyn (C) is a multibody analysis code.
$ http://www.mbdyn.org
$ 
$ Copyright (C) 1996-2017 
$
$ Pierangelo Masarati     <masarati@aero.polimi.it>
$ Paolo Mantegazza        <mantegazza@aero.polimi.it>
$
$ Author: Giuseppe Quaranta <quaranta@aero.polimi.it>
$         Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
$
$  These alters may not produce the desired datablocks if SOL 200
$  auxiliary models are present or the new partitioned superelement
$  option is used.
$  In v69+ some of the datablocks have changed to accomodate the
$  64-bit precision for all data.  These Alter statements are for v69+
$  only; they will not work with v68.2.
$
$  ASSIGN FILES FOR OUTPUT2 AND OUTPUT4
$  TO WRITE THE DATA BLOCKS ON
$
ASSIGN OUTPUT2='mbdyn.tab' STATUS=UNKNOWN UNIT=11  $ tables
ASSIGN OUTPUT4='mbdyn.mat' STATUS=UNKNOWN UNIT=15  $ matrices
$
$
$$$ EXECUTIVE CONTROL DECK
$
TIME 500          $ MAXIMUM COMPUTATION TIME (MIN.)
SOL 103           $ SOLUTION SEQUENCE SPECIFICATION
$
$--------------------------------------------------------------------------
COMPILE     PHASE0 
ALTER       'IF ( NOGEOM1>0 OR NOGEOM2>0 ) GP1'  $ v69+ only
MESSAGE     //' WRITE GPL AND BGPDT DATABLOCKS'
$
$ BACKUP v69 TABLES
$
TYPE        PARM,,CHAR8,N,BGPDT,GPL  $
EQUIVX      BGPDT/BGPDT69/-1 $
EQUIVX      GPL/GPL69/-1  $
DELETE      /BGPDT,GPL,,,  $
MAKEOLD     BGPDT69,,,,/BGPDT,GPL,,,/
	    'BGPDT'/////'BGPDT'/'GPL'/ $
BGPDT='BGPDT68'
GPL='GPL68'
OUTPUT2      GPL,BGPDT,,,//-1/11  $
EQUIVX       BGPDT69/BGPDT/-1 $
EQUIVX       GPL69/GPL/-1 $
$
$--------------------------------------------------------------------------
$ PUT THE MASS AND THE STIFFNESS MATRICES IN ZUZR DB TO
$ TRANSMIT THEM TO  THE FOLLOWING ALTER 
$ 
COMPILE      PHASE1A
ALTER        'ENDIF $ SELR'  $
TYPE         DB ZUZR01 $
TYPE         DB ZUZR02 $
EQUIVX       MGG/ZUZR01/-1 $
EQUIVX       KGG/ZUZR02/-1 $
$
$---------------------------------------------------------------------------
$ COMPUTE THE NEW MASS AND STIFFNESS MATRICES
$ OUTPUT OF MODAL MASS AND MATRIX, OF LUMPED MASS AND OF MODAL FORMS
$
COMPILE     SEDRCVR
ALTER       'SDR2     CASEDR,CSTMS,MPTS,DIT,EQEXINS,,ETT,OL2,BGPDTN,' 
MESSAGE     //' WRITE MHH, KHH, UG, AND LUMPED MASS DATABLOCKS'
$
TYPE        DB ZUZR01
TYPE        DB ZUZR02
EQUIVX      ZUZR01/MGG/-1 $
EQUIVX      ZUZR02/KGG/-1 $
$
SMPYAD      UG,MGG,UG,,,/MHH/3////1////6 $
SMPYAD      UG,KGG,UG,,,/KHH/3////1////6 $
DIAGONAL    MGG/LUMPMS/'COLUMN'/1.  $
OUTPUT4     MHH,KHH,LUMPMS,,//-1/15 $
OUTPUT4     ,,,,//-2/15 $
OUTPUT2     OUGV1,,,,//0/11 $
OUTPUT2     ,,,,//-9/11 $
ENDALTER
