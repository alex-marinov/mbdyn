%{
/*
 *
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 * 
 * Copyright (C) 1996-2007
 * 
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 * 
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 * 
 * Changing this copyright notice is forbidden.
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * 
 * 
 * ----------------------------------------------------------------------------
 *
 * ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
 * into raw MBDyn input files.
 * 
 * Copyright (C) 1999-2007
 * Leonardo Cassan		<lcassan@tiscalinet.it>
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include <ctype.h>
#include <b++.tab.h>
#include <string.h>

/* Generic declarations */

FILE* logfile=fopen("lex.log","w");
int ncount=1;
char ch;
char* header_buf= new char[255];
char* formula_text = new char[1024];
int buf_count=0;
int bcalling=0;
int Exp_State=0;

/* Prepare the variables for the expression parsing */
void LaunchParser (void);

/* Subs for formula */
void fconcat (char*);

%}


DIGIT      [0-9]
ALPHA      ([_])?[[:alpha:]*]
NUM        ([+-])?{DIGIT}*("."{DIGIT}+)?([eE]([+-])?{DIGIT}+)?
ASNUM      {DIGIT}*("."{DIGIT}+)?([eE]([+-])?{DIGIT}+)?
INVSLASH   [\\]
ID         {ALPHA}+[[:alnum:]*]*
IDN        [a-zA-Z]+[[:alnum:]*]*
SEP        [,]
TERM       [;]
NL	   [\n]
EQ	   [=]
EXC	   [!]
SLASH	   [/]
COLON      [:]
MUL        [*]

%s L ACCGRAV_START BEAM_START PART_START MARKER_START EQUILIBRIUM_START
%s GFORCE_START IC_START MATERIAL_START JOINT_START
%s JPRIM_START SPRINGDAMPER_START LSE_START RESULTS_START
%s REQUEST_START ADAMS_START UNITS_START SFORCE_START GRAPHICS_START
%s OUTPUT_START VARIABLE_START ARRAY_START

%x EXCL_START EXPRESSION_START

%%


^{EXC}  {
        fprintf (logfile,"%s",yytext);
	buf_count=0;
	BEGIN (EXCL_START);
	return EXCL;
}

{NUM}  { 
        fprintf (logfile,"%s",yytext); 
	yylval.value=atof(yytext);
	return NUM;
}


{NL}/{SEP} {
	ncount++;
	fprintf (logfile,"\n");
}

{SEP}{NL} {
	ncount++;
	fprintf (logfile,",\n");
	return SEP;
}


{SEP}   {
	fprintf (logfile,"%s",yytext);
	return SEP;
}

{TERM}  {
        fprintf (logfile,"%s",yytext);
        BEGIN(0);
        /* OPZIONE PUNTO E VIRGOLA DISATTIVATA (PER ATTIVARE return TERM;) */
        return NL;
}

{NL}    {
        BEGIN(0);
        fprintf (logfile,"\n");
        ncount++;
        return NL;
}

{EQ}    {
        fprintf (logfile,"%s",yytext);
        return EQ;
}

{EXC}	{
        fprintf (logfile,"%s",yytext);
        return EXCL;
}

{SLASH} {
        fprintf (logfile,"%s",yytext);
        return SLASH;
}

{COLON} {
        fprintf (logfile,"%s",yytext);
        return COLON;
}

{INVSLASH}  {
        fprintf(logfile,"%s",yytext);
	return INVSLASH;
}

([dD]?)    {
        fprintf (logfile,"%s",yytext);
        return D;
}

"("      {
        fprintf (logfile,"%s",yytext);
        Exp_State=Exp_State+1;
        return '(';
}

")"      {
        fprintf (logfile,"%s",yytext);
        Exp_State=Exp_State-1;
        return ')';
}



ACC	|
ACCG	|
ACCGR	|
ACCGRA	|
ACCGRAV	{
	fprintf (logfile,"%s",yytext);
	BEGIN(ACCGRAV_START);
	return ACCGRAV;
}




FU       |
FUN      |
FUNC     |
FUNCT    |
FUNCTI   |
FUNCTIO  |
FUNCTION {
         fprintf (logfile,"%s",yytext);
         LaunchParser();
	 return FUNCTION;
}



<ACCGRAV_START>	/* ACCGRAV START CONDITION */
{

I	|
IG	|
IGR	|
IGRA	|
IGRAV	{
	fprintf (logfile,"%s",yytext);
	return IGRAV;
}

J	|
JG	|
JGR	|
JGRA	|
JGRAV	{
	fprintf (logfile,"%s",yytext);
	return JGRAV;
}

K	|
KG	|
KGR	|
KGRA	|
KGRAV	{
	fprintf (logfile,"%s",yytext);
	return KGRAV;
}

}  /* END OF ACCGRAV START CONDITION */


ADAMS         {
              fprintf (logfile,"%s",yytext);
	      BEGIN(ADAMS_START);
	      return ADAMS;
}

<ADAMS_START> 
{

View          {
              fprintf (logfile,"%s",yytext);
	      return VIEW;
}

model         {
              fprintf (logfile,"%s",yytext);
	      return MODEL;
}

name          {
              fprintf (logfile,"%s",yytext);
	      BEGIN(0);
	      return NAME;
}

} /* END OF ADAMS HEADER TITLE */


BE	|
BEA	|
BEAM	{
	fprintf (logfile,"%s",yytext);
	BEGIN(BEAM_START);
	return BEAM;
}


<BEAM_START> /*ENTERING BEAM START CONDITION*/
{

I	{
	fprintf (logfile,"%s",yytext);
	return I;
}

J	{
	fprintf (logfile,"%s",yytext);
	return J;
}

L      |
LE     |
LEN    |
LENG   |
LENGT  |
LENGTH {
	fprintf (logfile,"%s",yytext);
	/* yylval.stringa="LENGTH"; */
	return LENGTH;
}

IX	|
IXX	{
	fprintf (logfile,"%s",yytext);
	return IXX;
}

IY	|
IYY	{
	fprintf (logfile,"%s",yytext);
	return IYY;
}

IZ	|
IZZ	{
	fprintf (logfile,"%s",yytext);
	return IZZ;
}

AREA	{
	fprintf (logfile,"%s",yytext);
	return AREA;
}

ASY	{
	fprintf (logfile,"%s",yytext);
	return ASY;
}

ASZ	{
	fprintf (logfile,"%s",yytext);
	return ASZ;
}

E               |
EM		|
EMO		|
EMOD		|
EMODU		|
EMODUL		|
EMODULU		|
EMODULUS	{
		fprintf (logfile,"%s",yytext);
		return EMODULUS;
}

G               |
GM		|
GMO		|
GMOD		|
GMODU		|
GMODUL		|
GMODULU		|
GMODULUS	{
		fprintf (logfile,"%s",yytext);
		return GMODULUS;
}

CM		|
CMA		|
CMAT		|
CMATR		|
CMATRI		|
CMATRIX		{
		fprintf (logfile,"%s",yytext);
		return CMATRIX;
}

CR		|
CRA		|
CRAT		|
CRATI		|
CRATIO		{
                fprintf (logfile,"%s",yytext);
		return CRATIO;
}

} /* END BEAM START CONDITION */


EQU		|
EQUI		|
EQUIL		|
EQUILI		|
EQUILIB		|
EQUILIBR	|
EQUILIBRI	|
EQUILIBRIU	|
EQUILIBRIUM	{
		fprintf (logfile,"%s",yytext);
		BEGIN (EQUILIBRIUM_START);
		return EQUILIBRIUM;
}


<EQUILIBRIUM_START> /* EQUILIBRIUM START CONDITION */
{

AL	|
ALI	|
ALIM	|
ALIMI	|
ALIMIT	{
	fprintf (logfile,"%s",yytext);
	return ALIMIT;
}

ERR	|
ERRO	|
ERROR	{
	fprintf (logfile,"%s",yytext);
	return ERROR;
}

IM		|
IMB		|
IMBA		|
IMBAL		|
IMBALA		|
IMBALAN		|
IMBALANC	|
IMBALANCE	{
		fprintf(logfile,"%s",yytext);
		return IMBALANCE;
}

MAX	|
MAXI	|
MAXIT	{
	fprintf (logfile,"%s",yytext);
	return MAXIT;
}

PAT	|
PATT	|
PATTE	|
PATTER	|
PATTERN	{
	fprintf (logfile,"%s",yytext);
	return PATTERN;
}

ST		|
STA		|
STAB		|
STABI		|
STABIL		|
STABILI		|
STABILIT	|
STABILITY	{
		fprintf (logfile,"%s",yytext);
		return STABILITY;
}

TL	|
TLI	|
TLIM	|
TLIMI	|
TLIMIT	{
	fprintf (logfile,"%s",yytext);
	return TLIMIT;
}

T       |
TRUE    {
        fprintf (logfile,"%s",yytext);
	return TRUE;
}

F       |
FALSE   {
        fprintf (logfile,"%s",yytext);
	return FALSE;
}


} /* END OF EQULIBRIUM START CONDITION */

<INITIAL>
{
M	|
MA	|
MAR	|
MARK	|
MARKE	|
MARKER	{
	fprintf (logfile,"%s",yytext);
	BEGIN(MARKER_START);
        return MARKER;
}
}

<MARKER_START> /* MARKER START CONDITION */
{

P	|
PA	|
PAR	|	
PART 	{
        fprintf (logfile,"%s",yytext);
	return PART;
}

FLO	|
FLOA	|
FLOAT	|
FLOATI	|
FLOATIN	|
FLOATING {
	fprintf (logfile,"%s",yytext);
	return FLOATING;
}

POI		|
POIN		|
POINT		|
POINT_		|
POINT_M		|
POINT_MA	|
POINT_MAS	|
POINT_MASS 	{
	fprintf (logfile,"%s",yytext);
	return POINT_MASS;
}

Q	|
QP	{
	fprintf (logfile,"%s",yytext);
	return QP;
}

Z	|
ZP	{
	fprintf (logfile,"%s",yytext);
	return ZP;
}

X	|
XP	{
	fprintf (logfile,"%s",yytext);
	return XP;
}

USEXP	{
	fprintf (logfile,"%s",yytext);
	return USEXP;
}

FLE		|
FLEX		|
FLEX_		|
FLEX_B		|
FLEX_BO		|
FLEX_BOD	|
FLEX_BODY 	{
	fprintf (logfile,"%s",yytext);
	return FLEX_BODY;
}

NOD	|
NODE	|
NODE_	|
NODE_I	|
NODE_ID	{
	fprintf (logfile,"%s",yytext);
	return NODE_ID;
}

R	|
RE	|
REU	|
REUL	|
REULE	|
REULER {
        fprintf (logfile,"%s",yytext);
	return REULER;
}

} /* END OF MARKER START CONDITION */

<INITIAL>
{
GF	|
GFO	|
GFOR	|
GFORC	|
GFORCE	{
	fprintf (logfile,"%s",yytext);
	BEGIN(GFORCE_START);
	return GFORCE;
}
}

<GFORCE_START>	/* START OF GFORCE START CONDITION */ 
{

I       {
        fprintf (logfile,"%s",yytext);
	return I;
}


JF	|
JFL	|
JFLO	|
JFLOA	|
JFLOAT	{
	fprintf (logfile,"%s",yytext);
	return JFLOAT;
}

RM	{
	fprintf (logfile,"%s",yytext);
	return RM;
}

FX	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
        return FX;
}

FY	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
	return FY;
}

FZ	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
	return FZ;
}

TX	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
	return TX;
}

TY	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
	return TY;
}

TZ	{
	fprintf (logfile,"%s",yytext);
	LaunchParser();
	return TZ;
}


} /* END OF GFORCE START CONDITION */

GRAPHICS {
         fprintf (logfile,"%s",yytext);
	 BEGIN (GRAPHICS_START);
         return GRAPHICS;
}

<GRAPHICS_START>
{

FORCE   {
        fprintf (logfile,"%s",yytext);
	return FORCE;
}

CM      {
        fprintf (logfile,"%s",yytext);
	return CM;
}

LENGTH  {
        fprintf (logfile,"%s",yytext);
	return LENGTH;
}

RADIUS  {
        fprintf (logfile,"%s",yytext);
	return RADIUS;
}

EMARKER {
        fprintf (logfile,"%s",yytext);
	return EMARKER;
}

EID     {
        fprintf (logfile,"%s",yytext);
        return EID;
}

ETYPE {
      fprintf(logfile,"%s",yytext);
      return ETYPE;
}

CYLINDER {
         fprintf (logfile,"%s",yytext);
	 return CYLINDER;
}

SFORCE   {
         fprintf (logfile,"%s",yytext);
	 return SFORCE;
	 }
	 
GFORCE   {
         fprintf (logfile,"%s",yytext);
	 return GFORCE;
	 }
	 
VFORCE   {
         fprintf (logfile,"%s",yytext);
	 return VFORCE;
}

} /* END OF GRAPHICS START CONDITION */


<INITIAL>IC      {
        fprintf (logfile,"%s",yytext);
        BEGIN (IC_START);
        return IC;
}

<IC_START>  /* IC START CONDITION */
{

T      |
TRUE   {
       fprintf (logfile,"%s",yytext);
       return TRUE;
}

F      |
FALSE  {
       fprintf (logfile,"%s",yytext);
       return FALSE;
}


AERR   |
AERRO  |
AERROR {
	fprintf (logfile,"%s",yytext);
        return AERROR;
}
    
AL     |
ALI    |
ALIM   |
ALIMI  |
ALIMIT {
       fprintf (logfile,"%s",yytext);
       return ALIMIT;
}

AMAX   |
AMAXI  |
AMAXIT {
       fprintf (logfile,"%s",yytext);
       return AMAXIT;
}
   
APAT     |
APATT    |
APATTE   |
APATTERN {
         fprintf (logfile,"%s",yytext);
         return APATTERN;
}

ERR      |
ERRO     |
ERROR    {
         fprintf (logfile,"%s",yytext);
         return ERROR;
}

MAX      |
MAXI     |
MAXIT    {
         fprintf (logfile,"%s",yytext);
         return MAXIT;
}

PAT      |
PATT     |
PATTE    |
PATTER   |
PATTERN  {
         fprintf (logfile,"%s",yytext);
         return PATTERN;
}

TL       |
TLI      |
TLIM     |
TLIMI    |
TLIMIT   {
         fprintf (logfile,"%s",yytext);
         return TLIMIT;
}
 
VER      |
VERR     |
VERRO    |
VERROR   {
         fprintf (logfile,"%s",yytext);
         return VERROR;
}

}  /* END OF IC START CONDITION */


<INITIAL>
{
J	|
JO	|
JOI	|
JOIN	|
JOINT	{
	fprintf (logfile,"%s",yytext);
	BEGIN (JOINT_START);
	return JOINT;
}
} /* obbligato a farlo, causa il fatto che esiste un duplicato in sintassi*/


<JOINT_START>   /* START OF JOIN CONDITIONS */
{

CONV	|
CONVE	|
CONVEL	{
	fprintf (logfile,"%s",yytext);
	return CONVEL;
}

C		|
CY		|
CYL		|
CYLI		|
CYLIN		|
CYLIND		|
CYLINDR		|
CYLINDRI	|
CYLINDRIC	|
CYLINDRICA	|
CYLINDRICAL	{
		fprintf (logfile,"%s",yytext);	
		return CYLINDRICAL;
}

FIX	|
FIXE	|
FIXED	{
	fprintf (logfile,"%s",yytext);
	return FIXED;
}

HO	|
HOO	|
HOOK	|
HOOKE	{
	fprintf (logfile,"%s",yytext);
	return HOOKE;
}

PL	|
PLA	|
PLAN	|
PLANA	|
PLANAR	{
	fprintf (logfile,"%s",yytext);
	return PLANAR;
}

RA	|
RAC	|
RACK	|
RACKP	|
RACKPI	|
RACKPIN	{
	fprintf (logfile,"%s",yytext);
	return RACKPIN;
}

RE		|
REV		|
REVO		|
REVOL		|
REVOLU		|
REVOLUT		|
REVOLUTE	{
		fprintf (logfile,"%s",yytext);
		return REVOLUTE;
}

SC	|
SCR	|
SCRE	|
SCREW	{
	fprintf (logfile,"%s",yytext);
	return SCREW;
}

PI	|
PIT	|
PITC	|
PITCH	{
	fprintf (logfile,"%s",yytext);
	return PITCH;
}

SP		|
SPH		|
SPHE		|
SPHER		|
SPHERI		|
SPHERIC		|
SPHERICA	|
SPHERICAL	{
		fprintf (logfile,"%s",yytext);
		return SPHERICAL;
}

T		|
TR		|
TRA		|
TRAN		|
TRANS		|
TRANSL		|
TRANSLA		|
TRANSLAT	|
TRANSLATI	|
TRANSLATIO	|
TRANSLATION	|
TRANSLATIONA	|
TRANSLATIONAL	{
		fprintf (logfile,"%s",yytext);
		return TRANSLATIONAL;
}

U		|
UN		|
UNI		|
UNIV		|
UNIVE		|
UNIVER		|
UNIVERS		|
UNIVERSA	|
UNIVERSAL	{
		fprintf (logfile,"%s",yytext);
		return UNIVERSAL;
}

ICT	|
ICTR	|
ICTRA	|
ICTRAN	{
	fprintf (logfile,"%s",yytext);
	return ICTRAN;
}

ICR	|
ICRO	|
ICROT	{
	fprintf (logfile,"%s",yytext);
	return ICROT;
}

DEL	|
DELT	|
DELTA	|
DELTA_	|
DELTA_V	{
	fprintf (logfile,"%s",yytext);
	return DELTA_V;
}

INN		|
INNE		|
INNER		|
INNER_		|
INNER_R		|
INNER_RA	|
INNER_RAD	|
INNER_RADI	|
INNER_RADIU	|
INNER_RADIUS	{
		fprintf (logfile,"%s",yytext);
		return INNER_RADIUS;
}

FRI		|
FRIC		|
FRICT		|
FRICTI		|
FRICTIO		|
FRICTION	{
		fprintf (logfile,"%s",yytext);
		return FRICTION;
}

MAXIMUM_D		|
MAXIMUM_DE		|
MAXIMUM_DEF		|
MAXIMUM_DEFO		|
MAXIMUM_DEFOR		|
MAXIMUM_DEFORM		|
MAXIMUM_DEFORMA		|
MAXIMUM_DEFORMAT	|
MAXIMUM_DEFORMATI	|
MAXIMUM_DEFORMATIO	|
MAXIMUM_DEFORMATION	{
			fprintf (logfile,"%s",yytext);
			return MAXIMUM_DEFORMATION;
}

MU_DYN_R	|
MU_DYN_RO	|
MU_DYN_ROT	{
		fprintf (logfile,"%s",yytext);
		return MU_DYN_ROT;
}

MU_STAT_R	|
MU_STAT_RO	|
MU_STAT_ROT	{
		fprintf (logfile,"%s",yytext);
		return MU_STAT_ROT;
}

OUT		|
OUTE		|
OUTER		|
OUTER_		|
OUTER_R		|
OUTER_RA	|
OUTER_RAD	|
OUTER_RADI	|
OUTER_RADIU	|
OUTER_RADIUS	{
		fprintf (logfile,"%s",yytext);
		return OUTER_RADIUS;
}

PRELOAD_R	|
PRELOAD_RA	|
PRELOAD_RAD	|
PRELOAD_RADI    |
PRELOAD_RADIA   |
PRELOAD_RADIAL  {
		fprintf (logfile,"%s",yytext);
		return PRELOAD_RADIAL;
}

PRELOAD_A	|
PRELOAD_AX	|
PRELOAD_AXI	|
PRELOAD_AXIA	|
PRELOAD_AXIAL	{
		fprintf (logfile,"%s",yytext);
		return PRELOAD_AXIAL;
}
  
HEI	|
HEIG	|
HEIGH	|
HEIGHT	{
	fprintf (logfile,"%s",yytext);
	return HEIGHT;
}

MU_DYN_T	|
MU_DYN_TR	|
MU_DYN_TRA	|
MU_DYN_TRAN	|
MU_DYN_TRANS	{
		fprintf (logfile,"%s",yytext);
		return MU_DYN_TRANS;
}

MU_STAT_T	|
MU_STAT_TR	|
MU_STAT_TRA	|
MU_STAT_TRAN	|
MU_STAT_TRANS	{
		fprintf (logfile,"%s",yytext);
		return MU_STAT_TRANS;
}

PRELOAD_X	{
		fprintf (logfile,"%s",yytext);
		return PRELOAD_X;
}

PRELOAD_Y	{
		fprintf (logfile,"%s",yytext);
		return PRELOAD_Y;
}

WID	|
WIDT	|
WIDTH	{
	fprintf (logfile,"%s",yytext);
	return WIDTH;
}

ON	{
	fprintf (logfile,"%s",yytext);
	return ON;
}

OFF	{
	fprintf (logfile,"%s",yytext);
	return OFF;
}

PRE		|
PREL		|
PRELO		|
PRELOA		|
PRELOAD		|
PRELOAD_	|
PRELOAD_O	|
PRELOAD_ON	|
PRELOAD_ONL	|
PRELOAD_ONLY	{
		fprintf (logfile,"%s",yytext);
		return PRELOAD_ONLY;
}

I	{
	fprintf (logfile,"%s",yytext);
	return I;
}

J	{
	fprintf (logfile,"%s",yytext);
	return J;
}

	
IC	{
	fprintf (logfile,"%s",yytext);
	return IC;
}

PD	{
	fprintf (logfile,"%s",yytext);
	return PD;
}

MAX_         |
MAX_F        |
MAX_FR       |
MAX_FRI      |
MAX_FRIC     |
MAX_FRIC_    |
MAX_FRIC_R   |
MAX_FRIC_RO  |
MAX_FRIC_ROT {
        fprintf (logfile,"%s",yytext);
	return MAX_FRIC_ROT;
}

}  /* END OF JOINT START CONDITION */


JP     |
JPR    |
JPRI   |
JPRIM  {
       fprintf (logfile,"%s",yytext);
       BEGIN(JPRIM_START);
       return JPRIM;
}

<JPRIM_START>
{

I     {
      fprintf (logfile,"%s",yytext);
      return I;
}

J     {
      fprintf (logfile,"%s",yytext);
      return J;
}

AT      |
ATP     |
ATPO    |
ATPOI   |
ATPOIN  |
ATPOINT {
        fprintf (logfile,"%s",yytext);
	return ATPOINT;
}

INL     |
INLI    |
INLIN   |
INLINE  {
        fprintf (logfile,"%s",yytext);
	return INLINE;
}

INP     |
INPL    |
INPLA   |
INPLAN  |
INPLANE {
        fprintf (logfile,"%s",yytext);
	return INPLANE;
}

OR          |
ORI         |
ORIE        |
ORIEN       |
ORIENT      |
ORIENTA     |
ORIENTAT    |
ORIENTATI   |
ORIENTATIO  |
ORIENTATION {
            fprintf (logfile,"%s",yytext);
	    return ORIENTATION;
}

PAR             |
PARA            |
PARAL           |
PARALL          |
PARALLE         |
PARALLEL        |
PARALLEL_       |
PARALLEL_A      |
PARALLEL_AX     |
PARALLEL_AXE    |
PARALLEL_AXES   {
                fprintf (logfile,"%s",yytext);
		return PARALLEL_AXES;
}

PERP            |
PERPE           |
PERPEN          |
PERPEND         |
PERPENDI        |
PERPENDIC       |
PERPENDICU      |
PERPENDICUL     |
PERPENDICULA    |
PERPENDICULAR   {
                fprintf (logfile,"%s",yytext);
		return PERPENDICULAR;
}

}  /* END OF JPRIM START CONDITION */


LSE {
    fprintf (logfile,"%s",yytext);
    BEGIN(LSE_START);
    return LSE;
}

<LSE_START>  /* BEGIN OF LSE START CONDITION */
{

X      {
       fprintf (logfile,"%s",yytext);
       return LSE_X;
}

A      {
       fprintf (logfile,"%s",yytext);
       return LSE_A;
}

U      {
       fprintf (logfile,"%s",yytext);
       return LSE_U;
}

C      {
       fprintf (logfile,"%s",yytext);
       return LSE_C;
}

Y      {
       fprintf (logfile,"%s",yytext);
       return LSE_Y;
}

B      {
       fprintf (logfile,"%s",yytext);
       return LSE_B;
}

IC     {
       fprintf (logfile,"%s",yytext);
       return IC;
}

STA            |
STAT           |
STATI          |
STATIC         |
STATIC_        |
STATIC_H       |
STATIC_HO      |
STATIC_HOL     |
STATIC_HOLD    {
               fprintf (logfile,"%s",yytext);
	       return STATIC_HOLD;
}

}  /* END OF LSE START CONDITION */

OUT     |
OUTP    |
OUTPU   |
OUTPUT  {
        fprintf (logfile,"%s",yytext);
	BEGIN(OUTPUT_START);
	return OUTPUT;
}

<OUTPUT_START>
{


REQ      |
REQS     |
REQSA    |
REQSAV   |
REQSAVE  {
         fprintf (logfile,"%s",yytext);
	 return REQSAVE;
}
   
GRS      |
GRSA     |
GRSAV    |
GRSAVE   {
         fprintf (logfile,"%s",yytext);
         return GRSAVE;
}
   
GR5      |
GR52     |
GR521    |
GR521S   |
GR521SA  |
GR521SAV {
         fprintf (logfile,"%s",yytext);
         return GR521SAV;
}
   
C        |
CH       |
CHA      |
CHAR     |
CHART    {
         fprintf (logfile,"%s",yytext);
         return CHART;
}

NOP      |
NOPR     |
NOPRI    |
NOPRIN   |
NOPRINT  {
         fprintf (logfile,"%s",yytext);
         return NOPRINT;
}
   
NOS           |
NOSE          |
NOSEP         |
NOSEPA        |
NOSEPAR       |
NOSEPARA      |
NOSEPARAT     |
NOSEPARATO    |
NOSEPARATOR   {
              fprintf (logfile,"%s",yytext);
              return NOSEPARATOR;
}
   
OS          |
OSF         |
OSFO        |
OSFOR       |
OSFORM      |
OSFORMA     |
OSFORMAT    {
            fprintf (logfile,"%s",yytext);
            return OSFORMAT;
}
   
FI          |
FIX         |
FIXE        |
FIXED       {
            fprintf (logfile,"%s",yytext);
            return FIXED;
}
   
TE          |
TEL         |
TELE        |
TELET       |
TELETY      |
TELETYP     |
TELETYPE    {
            fprintf (logfile,"%s",yytext);
            return TELETYPE;
}
   
V           |
VP          |
VPR         {
            fprintf (logfile,"%s",yytext);
            return VPR;
}
   
DS          |
DSC         |
DSCA        |
DSCAL       |
DSCALE      {
            fprintf (logfile,"%s",yytext);
            return DSCALE;
}

VS          |
VSC         |
VSCA        |
VSCAL       |
VSCALE      {
            fprintf (logfile,"%s",yytext);
            return VSCALE;
}

AS          |
ASC         |
ASCA        |
ASCAL       |
ASCALE      {
            fprintf (logfile,"%s",yytext);
            return ASCALE;
}

FS          |
FSC         |
FSCA        |
FSCAL       |
FSCALE      {
            fprintf (logfile,"%s",yytext);
            return FSCALE;
}

DZ          |
DZE         |
DZER        |
DZERO       {
            fprintf (logfile,"%s",yytext);
            return DZERO;
}

VZ          |
VZE         |
VZER        |
VZERO       {
            fprintf (logfile,"%s",yytext);
            return VZERO;
}

AZ          |
AZE         |
AZER        |
AZERO       {
            fprintf (logfile,"%s",yytext);
            return AZERO;
}

FZ          |
FZE         |
FZER        |
FZERO       {
            fprintf (logfile,"%s",yytext);
             return FZERO;
}

} /* END OF OUTPUT START CONDITION */


<INITIAL>
{
P	|
PA	|
PAR	|
PART	{
	fprintf (logfile,"%s",yytext);
        BEGIN (PART_START);
        return PART;
}
}

<PART_START> /* PART START CONDITION */
{

G	|
GR	|
GRO	|
GROU	|
GROUND	{
	fprintf (logfile,"%s",yytext);
	return GROUND;
}

M	|
MA	|
MAS	|
MASS	{
	fprintf (logfile,"%s",yytext);
	return MASS;
}

CM	{
	fprintf (logfile,"%s",yytext);
	return CM;
}

IM	{
	fprintf (logfile,"%s",yytext);
	return IM;
}

MATERIAL {
	fprintf (logfile,"%s",yytext);
	return MATERIAL;
}

IP	{
	fprintf (logfile,"%s",yytext);
	return IP;
}

Q	|
QG	{
	fprintf (logfile,"%s",yytext);
	return QG;
}

R	|
RE	|
REU	|
REUL	|
REULE	|
REULER	{
	fprintf (logfile,"%s",yytext);
	return REULER;
}

Z/{EQ}	|
ZG	{
	fprintf (logfile,"%s",yytext);
	return ZG;
}

X/{EQ}	|
XG	{
	fprintf (logfile,"%s",yytext);
        return XG;
}

EX	|
EXA	|
EXAC	|
EXACT	{
	fprintf (logfile,"%s",yytext);
	return EXACT;
}

X       {
	fprintf (logfile,"%s",yytext);
	return XCOORD;
}

Y       {
	fprintf (logfile,"%s",yytext);
	return YCOORD;
}

Z       {
	fprintf (logfile,"%s",yytext);
	return ZCOORD;
}


VM	{
	fprintf (logfile,"%s",yytext);
	return VM;
}

WM	{
	fprintf (logfile,"%s",yytext);
	return WM;
}


VX	{
	fprintf (logfile,"%s",yytext);
	return VX;
}

VY	{
	fprintf (logfile,"%s",yytext);
	return VY;
}

VZ	{
	fprintf (logfile,"%s",yytext);
	return VZ;
}

WX	{
	fprintf (logfile,"%s",yytext);
	return WX;
}

WY	{
	fprintf (logfile,"%s",yytext);
	return WY;
}

WZ	{
	fprintf (logfile,"%s",yytext);
	return WZ;
}


PSI	{
	fprintf (logfile,"%s",yytext);
	return PSI;
}

THETA	{
	fprintf (logfile,"%s",yytext);
	return THETA;
}

PHI	{
	fprintf (logfile,"%s",yytext);
	return PHI;
}

} /* END OF PART START CONDITION */


POI          |
POIN         |
POINT        |
POINT_       |
POINT_M      |
POINT_MA     |
POINT_MAS    |
POINT_MASS   {
             fprintf (logfile,"%s",yytext);
	     BEGIN(PART_START);
             return POINT_MASS;
}


R            |
RE           |
REQ          |
REQU         |
REQUE        |
REQUES       |
REQUEST      {
             fprintf (logfile,"%s",yytext);
	     BEGIN(REQUEST_START);
             return REQUEST;
}

<REQUEST_START>   /* BEGIN OF REQUEST START CONDITION */
{

F1      {
        fprintf (logfile,"%s",yytext);
	return F1;
}

F2      {
        fprintf (logfile,"%s",yytext);
	return F2;
}

F3      {
        fprintf (logfile,"%s",yytext);
	return F3;
}

F4      {
        fprintf (logfile,"%s",yytext);
	return F4;
}

F5      {
        fprintf (logfile,"%s",yytext);
	return F5;
}

F6      {
        fprintf (logfile,"%s",yytext);
	return F6;
}

F7      {
        fprintf (logfile,"%s",yytext);
	return F7;
}

F8      {
        fprintf (logfile,"%s",yytext);
	return F8;
}

DI             |
DIS            |
DISP           |
DISPL          |
DISPLA         |
DISPLAC        |
DISPLACE       |
DISPLACEM      |
DISPLACEME     |
DISPLACEMEN    |
DISPLACEMENT   {
               fprintf (logfile,"%s",yytext);
	       return DISPLACEMENT;
}

V              |
VE             |
VEL            |
VELO           |
VELOC          |
VELOCI         |
VELOCIT        |
VELOCITY       {
               fprintf (logfile,"%s",yytext);
	       return VELOCITY;
}

C        |
CO       |
COM      |
COMM     |
COMME    |
COMMEN   |
COMMENT  |
COMMENTS {
         fprintf (logfile,"%s",yytext);
         return COMMENT;
}

ACCE           |
ACCEL          |
ACCELE         |
ACCELER        |
ACCELERA       |
ACCELERAT      |
ACCELERATI     |
ACCELERATIO    |
ACCELERATION   {
               fprintf (logfile,"%s",yytext);
	       return ACCELERATION;
}

F              |
FO             |
FOR            |
FORC           |
FORCE          {
               fprintf (logfile,"%s",yytext);
	       return FORCE;
}

I  {
   fprintf (logfile,"%s",yytext);
   return I;
}

J  {
   fprintf (logfile,"%s",yytext);
   return J;
}

RM {
   fprintf (logfile,"%s",yytext);
   return RM;
}

T      |
TI     |
TIT    |
TITL   |
TITLE  {
       fprintf (logfile,"%s",yytext);
       return TITLE;
}

}  /* END OF REQUEST START CONDITION */



RES      |
RESU     |
RESUL    |
RESULT   |
RESULTS  {
         fprintf (logfile,"%s",yytext);
	 BEGIN(RESULTS_START);
         return RESULTS;
}

<RESULTS_START>  /* RESULT START CONDITION */
{

FORM       |
FORMA      |
FORMAT     |
FORMATT    |
FORMATTE   |
FORMATTED  {
           fprintf (logfile,"%s",yytext);
	   return FORMATTED;
}

NOACC            |
NOACCE           |
NOACCEL          |
NOACCELE         |
NOACCELER        |
NOACCELERA       | 
NOACCELERAT      |
NOACCELERATI     |
NOACCELLERATIO   |
NOACCELLERATIONS {
                 fprintf (logfile,"%s",yytext);
		 return NOACCELERATIONS;
}

NOAPP            |
NOAPPL           |
NOAPPLI          |
NOAPPLIE         |
NOAPPLIED        |
NOAPPLIEDF       |
NOAPPLIEDFO      |
NOAPPLIEDFOR     |
NOAPPLIEDFORC    |
NOAPPLIEDFORCE   |
NOAPPLIEDFORCES  {
                 fprintf (logfile,"%s",yytext);
		 return NOAPPLIEDFORCES;
}

NODAT               |
NODATA              |
NODATAS             |
NODATAST            |
NODATASTR           |
NODATASTRU          |
NODATASTRUC         |
NODATASTRUCT        |
NODATASTRUCTU       |
NODATASTRUCTUR      |
NODATASTRUCTURE     |
NODATASTRUCTURES    {
                    fprintf (logfile,"%s",yytext);
		    return NODATASTRUCTURES;
}

NODIS               |
NODISP              |
NODISPL             |
NODISPLA            |
NODISPLAC           |
NODISPLACE          |
NODISPLACEM         |
NODISPLACEME        |
NODISPLACEMEN       |
NODISPLACEMENT      |
NODISPLACEMENTS     {
                    fprintf (logfile,"%s",yytext);
		    return NODISPLACEMENTS;
}

NOFLO                  |
NOFLOA                 |
NOFLOAT                |
NOFLOATI               |
NOFLOATIN              |
NOFLOATING             |
NOFLOATINGM            |
NOFLOATINGMA           |
NOFLOATINGMAR          |
NOFLOATINGMARK         |
NOFLOATINGMARKE        |
NOFLOATINGMARKER       |
NOFLOATINGMARKERS      {
                       fprintf (logfile,"%s",yytext);
		       return NOFLOATINGMARKERS;
}

NOLIN        |
NOLINE       |
NOLINEA      |
NOLINEAR     {
             fprintf (logfile,"%s",yytext);
	     return NOLINEAR;
}

NOREA                |
NOREAC               |
NOREACT              |
NOREACTI             |
NOREACTIO            |
NOREACTION           |
NOREACTIONF          |
NOREACTIONFO         |
NOREACTIONFOR        |
NOREACTIONFORC       |
NOREACTIONFORCE      |
NOREACTIONFORCES     {
                     fprintf (logfile,"%s",yytext);
		     return NOREACTIONFORCES;
}

NOSYS           |       
NOSYST          |
NOSYSTE         |
NOSYSTEM        |
NOSYSTEME       |
NOSYSTEMEL      |
NOSYSTEMELE     |
NOSYSTEMELEM    |
NOSYSTEMELEME   |
NOSYSTEMELEMEN  |
NOSYSTEMELEMENT {
                fprintf (logfile,"%s",yytext);
		return NOSYSTEMELEMENTS;
}

NOTIR       |
NOTIRE      |
NOTIRES     {
            fprintf (logfile,"%s",yytext);
	    return NOTIRES;
}

NOVEL             |
NOVELO            |
NOVELOC           |
NOVELOCI          |
NOVELOCIT         |
NOVELOCITI        |
NOVELOCITIE       |
NOVELOCITIES      {
                  fprintf (logfile,"%s",yytext);
		  return NOVELOCITIES;
}

COM      |
COMM     |
COMME    |
COMMEN   |
COMMENT  {
         fprintf (logfile,"%s",yytext);
	 return COMMENT;
}

}

MATERIAL {
         fprintf (logfile,"%s",yytext);
	 BEGIN(MATERIAL_START);
         return MATERIAL;
}

<MATERIAL_START>
{

YOUNGS_MODULUS {
               fprintf (logfile,"%s",yytext);
	       return YOUNGS_MODULUS;
}

DENSITY        {
               fprintf (logfile,"%s",yytext);
	       return DENSITY;
}

NAME           {
               fprintf (logfile,"%s",yytext);
	       return NAME;
}

POISSONS_RATIO {
               fprintf (logfile,"%s",yytext);
	       return POISSONS_RATIO;
}

}  /* END OF MATERIAL START CONDITION */

<INITIAL>
SF      |
SFO     |
SFOR    |
SFORC   |
SFORCE  {
        fprintf (logfile,"%s",yytext);
	BEGIN (SFORCE_START);
	return SFORCE;
}

<SFORCE_START>
{

I       {
        fprintf (logfile,"%s",yytext);
	return I;
}

J       {
        fprintf (logfile,"%s",yytext);
	return J;
}

TR            |
TRA           |
TRAN          |
TRANS         |
TRANSL        |
TRANSLA       |
TRANSLAT      |
TRANSLATI     |
TRANSLATIO    |
TRANSLATION   |
TRANSLATIONA  |
TRANSLATIONAL {
             fprintf (logfile,"%s",yytext);
	     return TRANSLATION;
}

ROT          |
ROTA         |
ROTAT        |
ROTATI       |
ROTATIO      |
ROTATION     |
ROTATIONA    |
ROTATIONAL   {
             fprintf (logfile,"%s",yytext);
	     return ROTATION;
}

A            |
AC           |
ACT          |
ACTI         |
ACTIO        |
ACTION       |
ACTIONO      |
ACTIONON     |
ACTIONONL    |
ACTIONONLY   {
             fprintf (logfile,"%s",yytext);
	     return ACTIONONLY;
}


} /* END OF SFORCE START CONDITION */

SPR             |
SPRI            |
SPRIN           |
SPRING          |
SPRINGD         |
SPRINGDA        |
SPRINGDAM       |
SPRINGDAMP      |
SPRINGDAMPE     |
SPRINGDAMPER    {
                fprintf (logfile,"%s",yytext);
		BEGIN (SPRINGDAMPER_START);
                return SPRINGDAMPER;
}

<SPRINGDAMPER_START>  /* INITI OF SPRINGDAMPER START CONDITION */
{

I     {
      fprintf (logfile,"%s",yytext);
      return I;
}

J     {
      fprintf (logfile,"%s",yytext);
      return J;
}


TR              |
TRA             |
TRAN            |
TRANS           |
TRANSL          |
TRANSLA         |
TRANSLAT        |
TRANSLATI       |
TRANSLATIO      |
TRANSLATION     |
TRANSLATIONA    |
TRANSLATIONAL   {
                fprintf (logfile,"%s",yytext);
		return TRANSLATION;
}

RO              |
ROT             |
ROTA            |
ROTAT           |
ROTATI          |
ROTATIO         |
ROTATION        |
ROTATIONA       |
ROTATIONAL      {
                fprintf (logfile,"%s",yytext);
		return ROTATION;
}

C   {
    fprintf (logfile,"%s",yytext);
    return C0;
}

K   {
    fprintf (logfile,"%s",yytext);
    return K0;
}

F     |
FO    |
FOR   |
FORC  |
FORCE {
      fprintf (logfile,"%s",yytext);
      return FORCE;
}

L      |
LE     |
LEN    |
LENG   |
LENGT  |
LENGTH {
      fprintf (logfile,"%s",yytext);
      return LENGTH;
}

CT    {
      fprintf (logfile,"%s",yytext);
      return CT;
}

KT    {
      fprintf (logfile,"%s",yytext);
      return KT;
}

T
TO     |
TOR    |
TORQ   |
TORQU  |
TORQUE {
       fprintf (logfile,"%s",yytext);
       return TORQUE;
}

A      |
AN     |
ANG    |
ANGL   |
ANGLE  {
       fprintf (logfile,"%s",yytext);
       return ANGLE;
}

} /* END OF SPRINGDAMPER START CONDITION */


UNI     |
UNIT    |
UNITS   {
        fprintf (logfile,"%s",yytext);
	BEGIN (UNITS_START);
        return UNITS;
}

<UNITS_START>
{

FOR     |
FORC    |
FORCE   {
        fprintf (logfile,"%s",yytext);
        return FORCE;
}

MAS       |
MASS      {
          fprintf (logfile,"%s",yytext);
	  return MASS;
}

LEN       |
LENG      |
LENGT     |
LENGTH    {
        fprintf (logfile,"%s",yytext);
	return LENGTH;
}

TIM        |
TIME       {
           fprintf (logfile,"%s",yytext);
	   return TIME;
}

SYS       |
SYST      |
SYSTE     |
SYSTEM    {
          fprintf (logfile,"%s",yytext);
	  return SYSTEM;
}

UCF       {
          fprintf (logfile,"%s",yytext);
	  return UCF;
}

CGS       {
          fprintf (logfile,"%s",yytext);
	  return CGS;
}

FPS       {
          fprintf (logfile,"%s",yytext);
	  return FPS;
}

IPS       {
          fprintf (logfile,"%s",yytext);
	  return IPS;
}

MKS       {
          fprintf (logfile,"%s",yytext);
	  return MKS;
}

NON       |
NONE      {
          fprintf (logfile,"%s",yytext);
	  return NONE;
}

DYN     |
DYNE    {
        fprintf (logfile,"%s",yytext);
	return DYNE;
}

KILOGRAM_           |
KILOGRAM_F          |
KILOGRAM_FO         |
KILOGRAM_FOR        |
KILOGRAM_FORC       |
KILOGRAM_FORCE      {
         fprintf (logfile,"%s",yytext);
	 return KILOGRAM_FORCE;
}

KNEW     |
KNEWT    |
KNEWTO   |
KNEWTON  {
         fprintf (logfile,"%s",yytext);
	 return KNEWTON;
}

KPOUND_F     |
KPOUND_FO    |
KPOUND_FOR   |
KPOUND_FORC  |
KPOUND_FORCE {
             fprintf (logfile,"%s",yytext);
	     return KPOUND_FORCE;
}

NEW     |
NEWT    |
NEWTO   |
NEWTON  {
        fprintf (logfile,"%s",yytext);
	return NEWTON;
}

OUNCE_F      |
OUNCE_FO     |
OUNCE_FOR    |
OUNCE_FORC   |
OUNCE_FORCE  {
        fprintf (logfile,"%s",yytext);
	return OUNCE_FORCE;
}

POUND_F      |
POUND_FO     |
POUND_FOR    |
POUND_FORC   |
POUND_FORCE  {
        fprintf (logfile,"%s",yytext);
	return POUND_FORCE;
}

GRA       |
GRAM      {
          fprintf (logfile,"%s",yytext);
	  return GRAM;
}

KIL       |
KILO      |
KILOG     |
KILOGR    |
KILOGRA   |
KILOGRAM  {
          fprintf (logfile,"%s",yytext);
	  return KILOGRAM;
}

KPOUND_M    |
KPOUND_MA   |
KPOUND_MAS  |
KPOUND_MASS {
          fprintf (logfile,"%s",yytext);
	  return KPOUND_MASS;
}

OUNCE_M    |
OUNCE_MA   |
OUNCE_MAS  |
OUNCE_MASS {
           fprintf (logfile,"%s",yytext);
	   return OUNCE_MASS;
}

POUND_M    |
POUND_MA   |
POUND_MAS  |
POUND_MASS {
           fprintf (logfile,"%s",yytext);
	   return POUND_MASS;
}

SLUG       {
           fprintf (logfile,"%s",yytext);
	   return SLUG;
}

CENT         |
CENTI        |
CENTIM       |
CENTIME      |
CENTIMET     |
CENTIMETE    |
CENTIMETER   {
        fprintf (logfile,"%s",yytext);
	return CENTIMETER;
}

FOOT     {
         fprintf (logfile,"%s",yytext);
	 return FOOT;
}

KILOM       |
KILOME      |
KILOMET     |
KILOMETE    |
KILOMETER   {
         fprintf (logfile,"%s",yytext);
	 return KILOMETER;
}

INCH    {
        fprintf (logfile,"%s",yytext);
	return INCH;
}

MET     |
METE    |
METER   {
        fprintf (logfile,"%s",yytext);
	return METER;
}

MILLIM     |
MILLIME    |
MILLIMET   |
MILLIMETE  |
MILLIMETER {
           fprintf (logfile,"%s",yytext);
	   return MILLIMETER;
}

MILE       {
           fprintf (logfile,"%s",yytext);
	   return MILE;
}

HOU        |
HOUR       {
           fprintf (logfile,"%s",yytext);
	   return HOUR;
}

MILLIS      |
MILLISE     |
MILLISEC    |
MILLISECO   |
MILLISECOND {
            fprintf (logfile,"%s",yytext);
	    return MILLISECOND;
}

MIN       |
MINU      |
MINUT     |
MINUTE    {
          fprintf (logfile,"%s",yytext);
	  return MINUTE;
}

SEC       |
SECO      |
SECON     |
SECOND    {
          fprintf (logfile,"%s",yytext);
	  return SECOND;
}

} /* END OF TIME START CONDITION FOR UNITS */

VT       |
VTO      |
VTOR     |
VTORQ    |
VTORQU   |
VTORQUE  {
         fprintf (logfile,"%s",yytext);
	 BEGIN(GFORCE_START);
	 return VTORQUE;
}

<INITIAL>
VF       |
VFO      |
VFOR     |
VFORC    |
VFORCE   {
         fprintf (logfile,"%s",yytext);
	 BEGIN (GFORCE_START);
	 return VFORCE;
}

VAR      |
VARI     |
VARIA    |
VARIAB   |
VARIABL  |
VARIABLE {
         fprintf (logfile,"%s",yytext);
	 BEGIN (VARIABLE_START);
	 return VARIABLE;
}

AR       |
ARR      |
ARRA     |
ARRAY    {
         fprintf (logfile,"%s",yytext);
	 BEGIN (ARRAY_START);
         return _ARRAY_;
}

<ARRAY_START>
{

VARIABLES {
          fprintf (logfile,"%s",yytext);
	  return VARIABLE;
}

NUMBERS   {
          fprintf (logfile,"%s",yytext);
	  return NUMBERS;
}

S         |
SI        |
SIZ       |
SIZE      {
          fprintf (logfile,"%s",yytext);
	  return SIZE;
}

X         {
          fprintf (logfile,"%s",yytext);
	  return ARR_X;
}

Y         {
          fprintf (logfile,"%s",yytext);
	  return ARR_Y;
}

U         {
          fprintf (logfile,"%s",yytext);
	  return ARR_U;
}

I         |
IC        {
          fprintf (logfile,"%s",yytext);
	  return IC;
}

} /* END OF ARRAY START CONDITION */

<VARIABLE_START>
{

IC       {
         fprintf (logfile,"%s",yytext);
	 return IC;
}

} /* End of variable start */



<EXPRESSION_START>
{


{ASNUM} { 
      fprintf (logfile,"%s",yytext); 
      yylval.value=atof(yytext);
      fconcat (yytext);
      return NUM;
}

USER        {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return USER;
}

"+"         {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return '+';
}

"-"         {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return '-';
}

{MUL}{2}    {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return '@';
}

"*"          {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return '*';
}


"/"         {
      fprintf (logfile,"%s",yytext); 
      fconcat (yytext);
      return '/';
}

"^"         {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return '^';
}


"("      {
      fprintf (logfile,"%s",yytext);
      Exp_State=Exp_State+1;
      fconcat (yytext);
      return '(';
}

")"      {
      fprintf (logfile,"%s",yytext);
      Exp_State=Exp_State-1;
      fconcat (yytext);
      return ')';
}


"="         {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return EQ;
}

":"         {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return ':';
}

{NL}{SEP} {
      ncount++;
      fprintf (logfile,"\n");
}

{SEP}{NL} {
      ncount++;
      fprintf (logfile,",\n");
}


{SEP} {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      if (Exp_State==0) {
	 BEGIN (bcalling);
      }
      return SEP;
}

{TERM} {
      fprintf (logfile,"%s",yytext);
      BEGIN(0);
      /* OPZIONE PUNTO E VIRGOLA DISATTIVATA (PER ATTIVARE return TERM;) */
      return NL;
}

{NL}   {
      BEGIN(0);
      fprintf (logfile,"%s",yytext);
      ncount++;
      return NL;
}

{INVSLASH}    {
       /* Attenzione: l'utilizzo di \ disattiva il formula mode */
       fprintf(logfile,"%s",yytext);
       BEGIN(bcalling);
       return INVSLASH;
}

IF    {
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return IF;
}

{IDN}  {
      int l=strlen(yytext);
      yylval.stringa=new char [l];
      for (int i=0; i<l; i++) yylval.stringa[i]=yytext[i];
      yylval.stringa[l]='\0';
      fprintf (logfile,"%s",yytext);
      fconcat (yytext);
      return ID;
}



} /* EXPRESSION_END */



END	{
	fprintf (logfile,"%s",yytext);
	fprintf (logfile,"\n [NOTE: ADAMS WILL NOT READ FOLLOWING LINES.]");
	return END;
}

<EXCL_START>
{

.       {
        header_buf[buf_count++]=yytext[0];
	fprintf (logfile,"%s",yytext);
}

{NL}    {
	fprintf (logfile,"\n");
        ncount++;
        BEGIN(0);
}

} /* End of exclamation start */


{ID}  {
   int l=strlen(yytext);
   yylval.stringa=new char [l];
   for (int i=0; i<l; i++) yylval.stringa[i]=yytext[i];
   yylval.stringa[l]='\0';
   fprintf (logfile,"%s",yytext);
   return ID;
}


[^[:space:]] {
  fprintf (logfile,"%c",yytext);
  return NOT_CLASSIFIED;
}

%%


/* Auxiliary subroutines defined in the end of code cause declaration */

void LaunchParser (void)
{
   bcalling = YY_START;
   for (int i=0;i<1024;i++) formula_text[i]=' ';
   formula_text[0]='\0';
   BEGIN (EXPRESSION_START);
   Exp_State=0;
}

void fconcat (char* s)
{
   /* This sub concatenate the token at the formula */
   int j = strlen(formula_text);
   int q = strlen(s);
   for (int i=0;i<q;i++)
     formula_text[i+j]=s[i];
   formula_text[q+j]='\0';
   return;
}
