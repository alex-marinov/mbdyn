$ $Header$
$ Alters for MSC/NASTRAN v 70.7
$ ALTERS for SOL 101 STATIC MODES ANALYSIS
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
$  ASSIGN FILES FOR OUTPUT2 AND OUTPUT4
$  TO WRITE THE DATA BLOCKS ON
$
ASSIGN  OUTPUT4='mbdyn.stm' STATUS=UNKNOWN UNIT=19  $ static modes
$
$
$$$ EXECUTIVE CONTROL DECK
$
TIME 500          $ MAXIMUM COMPUTATION TIME (MIN.)
SOL 101           $ SOLUTION SEQUENCE SPECIFICATION
$
$-------------------------------------------------------------------------
COMPILE    SEDRCVR
ALTER      'SDR2     CASEDR,CSTMS,MPTS,DIT,EQEXINS,,ETT,OL2,BGPDTN,' $ 
MATGEN     EQEXINS/EXTINT/9//LUSETS $
MPYAD      EXTINT,UG,/UGVEXT/1 $ Displacements in External Sequence
OUTPUT4    UGVEXT//-1/19 $
ENDALTER

