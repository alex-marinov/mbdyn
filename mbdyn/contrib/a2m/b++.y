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

// 31-7-2000 - DATA DELL'ULTIMO AGGIORNAMENTO


// INCLUDE FILES
#include <carddecl.h>
#include <translate.h>
#include <helpers.h>
#include <fstream.h>

// VARIABLE DECLARATIONS
// EXTERNAL VARIABLES COMING FROM THE LEXER
/* dichiarate EXTERN perché non è possibile dichiararle, essendo
 * create all'interno del programma dallo stesso flex */

extern int yylex();
extern int yyerror(char *s);
extern FILE* yyin;
extern FILE* logfile;
extern int ncount;
extern FILE* yyout;
extern char* formula_text;
extern void azzera_formula();

// DICHIARAZIONI DI VARIABILI DI ERRORE GLOBALE
/* Utilizzate dalle routine di errore */

int nerr;
int nwarnings;
int nerror;
FILE* errfile;
FILE* messagefile;
FILE* tablefile;
Boolean ES;

// DICHIARAZIONE DELLE STRUTTURE DI MEMORIZZAZIONE DELLE CARD ADAMS
// UN CARD DECK PER OGNI TIPOLOGIA DI CARD

card_deck beams;
card_deck markers;
card_deck parts;
card_deck joints;
card_deck materials;
card_deck accgravs;
card_deck equilibriums;
card_deck jprims;
card_deck ics;
card_deck springdampers;
card_deck pointmasses;
card_deck sforces;
card_deck gforces;
card_deck vtorques;
card_deck vforces;
card_deck variables;

// DICHIARAZIONE DELLA STRUTTURA DI MEMORIZZAZIONE
MTP_deck Marker_Table;
PTR_deck Part_Table;
PTR_deck Pmass_Table;
MTR_deck Reference_Table;
PTN_deck PartToNode_Table;
BTP_deck Body_Table;
BTR_deck Beam_reference_Table;
   
// DICHIARAZIONE DELLE STRUTTURE DI MEMORIZZAZIONE DELLE CARD MBDYN
MBDyn_deck MBReference;
MBDyn_deck MBBeams;
MBDyn_deck MBJoints;
MBDyn_deck MBBodies;
MBDyn_deck MBNodes;
MBDyn_deck MBConstLaw;
MBDyn_deck MBForces;
MBDyn_deck MBAccs;

// VARIABILI GLOBALI DI VARIA NATURA

Id GROUND_PART_ID=0;

/* CREA IL SET DI VARIABILI DI CONTROLLO A2M */

double A2M_initial_time = 0;
double A2M_final_time = 1;
double A2M_time_step = 1.e-3;
double A2M_rho = .6;
char A2M_method [] = "ms, rho, rho";
double A2M_tolerance = 1.e-6;
double A2M_max_iterations = 20;

double A2M_deriv_coef = 1.e-6;
double A2M_deriv_tol = 1.e-9;
double A2M_deriv_iter = 10;

// PER OGNI TIPO DI ELEMENTO CREA L'OPPORTUNO PUNTATORE
// Un puntatore deve essere sempre disponibile a contenere la card
// corrente, qualora venga riconosciuta

s_beam* curr_beam = new s_beam;
s_marker* curr_marker = new s_marker;
s_part* curr_part  = new s_part;
s_joint* curr_joint = new s_joint;
s_material* curr_material = new s_material;
s_accgrav* curr_accgrav = new s_accgrav;
s_equilibrium* curr_equilibrium= new s_equilibrium;
s_jprim* curr_jprim = new s_jprim;
s_ic* curr_ic = new s_ic;
s_springdamper* curr_springdamper= new s_springdamper;
s_pointmass* curr_pointmass = new s_pointmass;
s_sforce* curr_sforce = new s_sforce;
s_gforce* curr_gforce = new s_gforce;
s_vforce* curr_vforce = new s_vforce;
s_vtorque* curr_vtorque = new s_vtorque;
s_variable* curr_variable = new s_variable;


// DICHIARAZIONE DI VARIABILI PER IL TOKEN STORAGE
token_deck tokens;
p_token_entry tok_idx;

// VARIABILI DI VARIA UTILITA'
const double JotaMin=1e-12;
int ntokens;
int i_arg=0;
double null_param=0;

// MODI DI FUNZIONAMENTO DELL'INTERFACCIA - LINEA DI COMANDO
Boolean VERBOSE_MODE=N;
Boolean OVERWRITE_MODE=N;
Boolean THROW_MODE=Y;
Boolean SKIP_MODE=N;
Boolean EXTENDED_MATRIX_DISPLAY=N;
Boolean SPECIFY_FILE=N;
Boolean REMOVE_REMARK=N;

// VETTORI DI SERVIZIO PER LO SCAMBIO DI INFORMAZIONI
// Utilizzati all'interno del codice per il trasferimento dei dati
// all'interno delle strutture grammaticali

Vec3 Service_vector3; 
Vec2 Service_vector2;
Euler Angle_Vector;       
Angle Service_Angle(0,RADIANS);
double* coef_matrix=new double[21];
double Service_value=0;
coord_type* Coord_Vector=new coord_type[10];
Boolean* Bool_Vector=new Boolean[30];

// DICHIARAZIONI DELLE SUBROUTINE INTERNE ALLA PARTE BISON
void Init_Token(void);
inline const char* Find_Token(int);
void Verbose();
void Report_form();

%}

%union {
	double value;
	char carat;
	unsigned int unreal;
	char* stringa;
	char angle_spec;
}


%token <value>   NUM NUM_DEG
%token <carat>   SEP TERM NL EQ EXCL SLASH COLON INVSLASH
%token <stringa> ID
%token <angle_type> ANG_TYPE
%token D

%token TRUE FALSE

%token BEAM PART I J LENGTH AREA
%token IXX IYY IZZ ASY ASZ EMODULUS GMODULUS
%token CMATRIX CRATIO

%token PART GROUND MASS CM IM IP
%token QG REULER ZG XG VX VY VZ WX WY WZ VM WM EXACT IP
%token PSI THETA PHI XCOORD YCOORD ZCOORD 

%token MARKER USEXP POINT_MASS FLEX_BODY NODE_ID FLOATING
%token QP XP ZP 

%token ACCGRAV IGRAV JGRAV KGRAV
%token END EQUILIBRIUM ALIMIT ERROR IMBALANCE MAXIT PATTERN STABILITY TLIMIT

%token GFORCE JFLOAT RM FX FY FZ TX TY TZ FUNCTION USER

%token IC AERROR AMAXIT APATTERN VERROR

%token MATERIAL DENSITY NAME YOUNGS_MODULUS POISSONS_RATIO

%token JOINT
%token CONVEL CYLINDRICAL FIXED HOOKE PLANAR RACKPIN
%token REVOLUTE SCREW PITCH SPHERICAL TRANSLATIONAL UNIVERSAL
%token ICTRAN ICROT 
%token DELTA_V INNER_RADIUS FRICTION MAXIMUM_DEFORMATION
%token MU_DYN_ROT MU_STAT_ROT OUTER_RADIUS PRELOAD_RADIAL
%token PRELOAD_AXIAL
%token HEIGHT MU_DYN_TRANS MU_STAT_TRANS PRELOAD_X PRELOAD_Y
%token WIDTH PRELOAD_ONLY RACKPID ON OFF
%token IC PD MAX_FRIC_ROT

%token JPRIM ATPOINT INLINE INPLANE ORIENTATION PARALLEL_AXES
%token PERPENDICULAR

%token POINT_MASS

%token SPRINGDAMPER TRANSLATION ROTATION
%token CT KT TORQUE ANGLE FORCE C0 K0

%token LSE LSE_C LSE_X LSE_A STATIC_HOLD LSE_U LSE_B LSE_Y

%token RESULTS FORMATTED COMMENT
%token NOACCELERATIONS NOAPPLIEDFORCES NODATASTRUCTURES
%token NODISPLACEMENTS NOFLOATINGMARKERS
%token NOLINEAR NOREACTIONFORCES NOSYSTEMELEMENTS NOTIRES
%token NOVELOCITIES

%token REQUEST DISPLACEMENT VELOCITY ACCELERATION FORCE
%token COMMENTS TITLE F1 F2 F3 F4 F5 F6 F7 F8

%token UNITS FORCE DYNE KILOGRAM_FORCE KNEWTON KPOUND_FORCE
%token NEWTON OUNCE_FORCE POUND_FORCE
%token MASS GRAM KILOGRAM KPOUND_MASS OUNCE_MASS POUND_MASS SLUG
%token LENGTH CENTIMETER FOOT KILOMETER INCH METER MILLIMETER MILE
%token TIME HOUR MILLISECOND MINUTE SECOND
%token SYSTEM CGS FPS IPS MKS NONE
%token UCF
%token SFORCE ACTIONONLY
%token VTORQUE VFORCE

%token ADAMS VIEW MODEL NAME

%token REQUEST GRAPHICS OUTPUT
%token VARIABLE

%token _ARRAY_ NUMBERS VARIABLES SIZE ARR_X ARR_Y ARR_U
%token EMARKER EID CYLINDER ETYPE RADIUS

%token OUTPUT GRSAVE GR521SAV CHART NOPRINT NOSEPARATOR OSFORMAT
%token FIXED TELETYPE VPR DSCALE VSCALE ASCALE FSCALE
%token DZERO VZERO AZERO FZERO REQSAVE

// Definisce il token _SET_ per il passaggio di parametri indiretti
// e il token NON_CLASSIFIED per i caratteri non riconosciuti
%token _SET_
%token NOT_CLASSIFIED

%token IF

%left  '-' '+'
%left '*' '/'
%left '@'
%right '^'
%left NEG

%%

/* GENERIC INPUT PARSING */


input: /* empty */
       | input line
;

line:   NL		{}
        | statement     {}
	| remark {}
	| sim_title NL {}
	| error	NL {error_recovery(ncount,logfile); }
;

statement: BEAM SLASH NUM SEP beam_arglist { 
                Set_Card (curr_beam,(Id) $3, beams); } 
	|  PART SLASH NUM SEP part_arglist {
		Set_Card (curr_part,(Id) $3, parts); }
	|  MARKER SLASH NUM mark_arglist {
	        Set_Card (curr_marker, (Id) $3, markers); }
	|  ACCGRAV SLASH accgravlist {
	        Set_Card (curr_accgrav, 0, accgravs); }
	|  GFORCE SLASH NUM SEP gforce_arglist {
	        Set_Card (curr_gforce,(Id) $3, gforces); }
        |  JOINT SLASH NUM SEP joint_arglist{
	        Set_Card (curr_joint,(Id) $3, joints);}
	|  POINT_MASS SLASH NUM SEP pointmass_arglist {
	        Set_Card (curr_pointmass,(Id) $3, pointmasses); }
	|  SFORCE SLASH NUM SEP sforce_arglist {
	        Set_Card (curr_sforce, (Id) $3, sforces); }
        |  VFORCE SLASH NUM SEP vforce_arglist {
	        Set_Card (curr_vforce, (Id) $3, vforces); }
	|  VTORQUE SLASH NUM SEP vtorque_arglist {
	        Set_Card (curr_vtorque, (Id) $3, vtorques); }
	|  END {}
	|  JPRIM SLASH NUM SEP jprim_arglist {
	        Set_Card (curr_jprim,(Id) $3, jprims);}
	|  VARIABLE SLASH NUM SEP variable_arglist {
	        Set_Card (curr_variable, (Id) $3, variables); }
        |  EQUILIBRIUM SLASH equil_arglist {
	        Set_Card (curr_equilibrium, -2, equilibriums); }
        |  SPRINGDAMPER SLASH NUM SEP springdamper_arglist {
                Set_Card (curr_springdamper,(Id) $3, springdampers); }
	|  IC SLASH ic_arglist {
	        Set_Card (curr_ic, -3, ics); }
	|  MATERIAL SLASH NUM SEP material_arglist{
	        Set_Card (curr_material,(Id) $3, materials);}
        |  LSE SLASH NUM SEP lse_arglist {
	   out_warning ("CARD LSE NOT IMPLEMENTED YET"); }
	|  RESULTS SLASH results_arglist {
	   out_warning ("CARD RESULTS NOT IMPLEMENTED YET"); }
	|  UNITS SLASH units_arglist {
	   out_warning ("CARD UNITS NOT IMPLEMENTED YET"); }
	|  REQUEST SLASH NUM SEP request_arglist {
	   out_warning ("CARD REQUEST NOT IMPLEMENTED YET"); }
	|  GRAPHICS SLASH NUM SEP graphics_arglist {
	   out_warning ("CARD GRAPHICS NOT IMPLEMENTED YET"); }
	|  OUTPUT SLASH output_arglist {
	   out_warning ("CARD OUTPUT NOT IMPLEMENTED YET"); }
	|  _ARRAY_ SLASH NUM SEP array_arglist {
	   out_warning ("CARD ARRAY NOT IMPLEMENTED YET"); }
;

remark: EXCL {}
;

sim_title: ADAMS SLASH VIEW MODEL NAME COLON idlist {}
;

/* END OF GENERIC INPUT PARSING */


//------------------------------ ACC GRAV START------------------------------

accgravlist:	  accgrav_arg NL {}
               |  accgrav_arg SEP accgravlist {}
	       |  error SEP accgravlist {}
	       |  error NL {}
;

accgrav_arg:      IGRAV EQ NUM { curr_accgrav->Set_Part(IGRAV,$3); }
               |  JGRAV EQ NUM { curr_accgrav->Set_Part(JGRAV,$3); }
	       |  KGRAV EQ NUM { curr_accgrav->Set_Part(KGRAV,$3); }
;


//------------------------------ ARRAY START --------------------------------
/* warning : strict grammatical form: terms "variables" & "numbers" must be
   in the end of the statement line */
   

array_arglist:   par SEP par SEP varlist NL {}
               | par SEP varlist NL {}
	       | par SEP par SEP numblist NL {}
	       | par SEP numblist NL {}
	       | par NL {}
	       | par SEP par NL {}
;
 
par:            ARR_X       {}
              | ARR_Y       {}
	      | ARR_U       {}
	      | IC          {}
	      | SIZE EQ NUM {}
	      | error       {}
;

num_list:       NUM {}
              | NUM SEP num_list {}
;

varlist:       VARIABLE EQ num_list {}
;

numblist:      NUMBERS EQ num_list {}
;

//------------------------------- BEAM START --------------------------------

beam_arglist:    beam_argument NL {}
               | beam_argument SEP beam_arglist {}
	       | error SEP beam_arglist { yyclearin; i_arg=0;}
	       | error NL { i_arg=0; }
;

beam_argument:	 I EQ NUM { curr_beam->Set_Part(I,(Id)$3); }
               | J EQ NUM { curr_beam->Set_Part(J,(Id)$3); }
	       | LENGTH EQ NUM { curr_beam->Set_Part(LENGTH,$3); }
	       | IXX EQ NUM { curr_beam->Set_Part(IXX, $3); }
	       | IYY EQ NUM { curr_beam->Set_Part(IYY, $3); }
	       | IZZ EQ NUM { curr_beam->Set_Part(IZZ, $3); }
	       | AREA EQ NUM { curr_beam->Set_Part (AREA, $3); }
	       | ASZ EQ NUM { curr_beam->Set_Part( ASZ, $3); }
	       | ASY EQ NUM { curr_beam->Set_Part (ASY, $3); }
	       | EMODULUS EQ NUM { curr_beam->Set_Part (EMODULUS,$3); }
	       | GMODULUS EQ NUM { curr_beam->Set_Part (GMODULUS,$3); }
	       | cmatrix { 
	            curr_beam->Set_Part(CMATRIX,coef_matrix);
		    i_arg=0; }
	       | CRATIO EQ NUM { curr_beam->Set_Part(CRATIO,$3); }
;

cmatrix:       CMATRIX EQ block SEP block SEP block {}
;

block:         NUM SEP NUM SEP NUM SEP NUM SEP NUM SEP NUM SEP NUM {
               coef_matrix[i_arg++]=$1;
	       coef_matrix[i_arg++]=$3;
	       coef_matrix[i_arg++]=$5;
	       coef_matrix[i_arg++]=$7;
	       coef_matrix[i_arg++]=$9;
	       coef_matrix[i_arg++]=$11;
	       coef_matrix[i_arg++]=$13; }	       
;

// ----------------------------- EQULIBRIUM START ----------------------------

equil_arglist:   equi_argument NL
               | equi_argument SEP equil_arglist {}
               | error SEP equil_arglist { yyclearin; i_arg=0; }
               | error NL { i_arg=0; }
;

equi_argument:  alimit { curr_equilibrium->Set_Part(ALIMIT,Service_Angle);}
               |_error { curr_equilibrium->Set_Part(ERROR,Service_value);}
	       |equi_imbalance {}
	       |maxit {curr_equilibrium->Set_Part(MAXIT,(Id) Service_value);}
	       |pattern { curr_equilibrium->ncoord_pattern=i_arg;
	                  curr_equilibrium->Set_Part(PATTERN,Bool_Vector);
		          i_arg=0;}
	       |equi_stab {}
	       |tlimit {curr_equilibrium->Set_Part(TLIMIT,Service_value);}
;

equi_stab:      STABILITY EQ NUM { curr_equilibrium->Set_Part(STABILITY,$3);}
;
equi_imbalance:	IMBALANCE EQ NUM { curr_equilibrium->Set_Part(IMBALANCE,$3);}
;

//------------------------------ GFORCE START -------------------------------

gforce_arglist:	 gforce_arg NL {}
               | gforce_arg SEP gforce_arglist {}
	       | error SEP gforce_arglist { yyclearin; }
	       | error NL {}
;

gforce_arg:    I EQ NUM { curr_gforce->Set_Part(I,(Id) $3); }
             | JFLOAT EQ NUM { curr_gforce->Set_Part(JFLOAT,(Id) $3); }
	     | RM EQ NUM { curr_gforce->Set_Part(RM,(Id) $3); }
	     | function_user { curr_gforce->Store_Formula("FUNCTION",formula_text);
	                       curr_gforce->Set_Part(FUNCTION,null_param); }
	     | gforce_magnitude {}
;

gforce_magnitude:  gforcem_arg {}
                 | gforcem_arg INVSLASH gforce_magnitude {}
;

gforcem_arg :  FX EQ magnitude { curr_gforce->Store_Formula("FX",formula_text);
                                 curr_gforce->Set_Part(FX,Service_value); }
             | FY EQ magnitude { curr_gforce->Store_Formula("FY",formula_text);
	                         curr_gforce->Set_Part(FY,Service_value); }
	     | FZ EQ magnitude { curr_gforce->Store_Formula("FZ",formula_text);
	                         curr_gforce->Set_Part(FZ,Service_value); }
	     | TX EQ magnitude { curr_gforce->Store_Formula("TX",formula_text);
	                         curr_gforce->Set_Part(TX,Service_value); }
	     | TY EQ magnitude { curr_gforce->Store_Formula("TY",formula_text);
	                         curr_gforce->Set_Part(TY,Service_value); }
	     | TZ EQ magnitude { curr_gforce->Store_Formula("TZ",formula_text);
	                         curr_gforce->Set_Part(TZ,Service_value); }
;

//----------------------------- GRAPHICS START -----------------------------//

graphics_arglist: graphics_argument {}
                | graphics_argument SEP graphics_arglist {}
		| error SEP graphics_arglist { yyclearin; i_arg=0; }
		| error NL { i_arg=0; }
;

graphics_argument:   FORCE {}
                   | CYLINDER {}
		   | ETYPE EQ force_type {}
		   | EID EQ NUM {}
		   | EMARKER EQ NUM {}
		   | LENGTH EQ NUM {}
		   | RADIUS EQ NUM {}
		   | CM EQ NUM {}
;

force_type:   SFORCE {}
            | GFORCE {}
	    | VFORCE {}
	    | VTORQUE {}
;

//-------------------------------- IC START ---------------------------------

ic_arglist:     ic_argument NL {}
		| ic_argument SEP ic_arglist {}
                | error SEP ic_arglist { yyclearin; i_arg=0; }
                | error NL { i_arg=0; }
;

ic_argument:	ic_aerror {} 
		| alimit {}
		| ic_amaxit {}	
		| ic_apattern {}
		| _error { curr_ic->Set_Part(ERROR,Service_value);}
		| maxit { curr_ic->Set_Part(MAXIT,(int) Service_value);}
		| pattern { curr_ic->ncoord_pattern=i_arg; 
			curr_ic->Set_Part(PATTERN,Bool_Vector);
			i_arg=0;}
		| tlimit { curr_ic->Set_Part(TLIMIT,Service_value);}
		| ic_verror {}
;

ic_aerror:	AERROR EQ NUM { curr_ic->Set_Part(AERROR,$3); }
;

ic_amaxit:	AMAXIT EQ NUM { curr_ic->Set_Part(AMAXIT,(int) $3); }
;

ic_apattern:	APATTERN EQ J_coordlist { curr_ic->ncoord_apattern=i_arg;
		  curr_ic->Set_Part(APATTERN,Bool_Vector);
		  i_arg=0;}
;

ic_verror:	VERROR EQ NUM { curr_ic->Set_Part(VERROR,$3);}
;

// ------------------------------ JOINT START --------------------------------

joint_arglist: 	joint_argument NL {}
               	| joint_argument SEP joint_arglist{}
                | error SEP joint_arglist { yyclearin;}
                | error NL {}
;

joint_argument: joint_type {}
		| nodeid {}
	      	| generic_param {}
	      	| joint_def {}
	      	| joint_deltav {}
	      	| jtl_paramspec {}
	      	| joint_friction {}
;

joint_type:       CONVEL { curr_joint->Set_Part(_SET_,_CONVEL); }
               	| CYLINDRICAL { curr_joint->Set_Part(_SET_,_CYLINDRICAL);}
	       	| FIXED { curr_joint->Set_Part(_SET_,_FIXED); }
	       	| HOOKE { curr_joint->Set_Part(_SET_,_HOOKE); }
	       	| PLANAR { curr_joint->Set_Part(_SET_,_PLANAR); }
	       	| RACKPIN { curr_joint->Set_Part(_SET_,_RACKPIN); }
	       	| REVOLUTE { curr_joint->Set_Part(_SET_,_REVOLUTE);}
	       	| SCREW { curr_joint->Set_Part(_SET_,_SCREW); }
	       	| SPHERICAL { curr_joint->Set_Part (_SET_,_SPHERICAL); }
	       	| TRANSLATIONAL {curr_joint->Set_Part(_SET_,_TRANSLATIONAL);}
	       	| UNIVERSAL { curr_joint->Set_Part (_SET_,_UNIVERSAL); }
;

nodeid:          I EQ NUM { curr_joint->Set_Part(I,(Id)$3); }
		| J EQ NUM { curr_joint->Set_Part(J,(Id)$3); }
;

generic_param:    ICTRAN EQ NUM SEP NUM {
	          Vec2 Aux($3,$5); 
		  curr_joint->Set_Part(ICTRAN,Aux); }
		| ICROT EQ NUM SEP NUM {
	          Vec2 Aux($3,$5);
		  curr_joint->Set_Part(ICROT,Aux); }
		| IC EQ NUM SEP NUM {
	          Vec2 Aux($3,$5);
		  curr_joint->Set_Part(IC,Aux); }
               	| INNER_RADIUS EQ NUM {curr_joint->Set_Part(INNER_RADIUS,$3);}
	       	| MU_DYN_ROT EQ NUM {curr_joint->Set_Part(MU_DYN_ROT,$3); }
	       	| MU_STAT_ROT EQ NUM  {curr_joint->Set_Part(MU_STAT_ROT,$3); }
	       	| OUTER_RADIUS EQ NUM {curr_joint->Set_Part(OUTER_RADIUS,$3);}
	       	| PRELOAD_RADIAL EQ NUM {curr_joint->Set_Part(PRELOAD_RADIAL,$3);}
	       	| PRELOAD_AXIAL EQ NUM {curr_joint->Set_Part(PRELOAD_AXIAL,$3);}
	       	| PITCH EQ NUM {curr_joint->Set_Part(PITCH,$3); }
	       	| PD EQ NUM {curr_joint->Set_Part(PD,$3); }
		| MAX_FRIC_ROT EQ NUM {curr_joint->Set_Part(MAX_FRIC_ROT,$3); }
;


friction_flag:   ON { curr_joint->Set_Part(FRICTION,_ON);}
               | OFF { curr_joint->Set_Part(FRICTION,_OFF);}
	       | PRELOAD_ONLY {curr_joint->Set_Part(FRICTION,_PRELOAD_ONLY);}
;

jtl_paramspec :  HEIGHT EQ NUM { curr_joint->Set_Part(HEIGHT,$3);}
               | PRELOAD_X EQ NUM {curr_joint->Set_Part(PRELOAD_X,$3);}
	       | PRELOAD_Y EQ NUM {curr_joint->Set_Part(PRELOAD_Y,$3);}
	       | WIDTH EQ NUM { curr_joint->Set_Part(WIDTH,$3);}
	       | MU_DYN_TRANS EQ NUM { curr_joint->Set_Part(MU_DYN_TRANS,$3);}
	       | MU_STAT_TRANS EQ NUM {curr_joint->Set_Part(MU_STAT_TRANS,$3);}
;

joint_friction:  FRICTION EQ friction_flag {}
;
joint_deltav:    DELTA_V EQ NUM { curr_joint->Set_Part(DELTA_V,$3);}
;
joint_def:       MAXIMUM_DEFORMATION EQ NUM {
                 curr_joint->Set_Part(MAXIMUM_DEFORMATION,$3);}
;

//----------------------------- JPRIM START ---------------------------------

jprim_arglist:   jprim_argument NL {}
               | jprim_argument SEP jprim_arglist {}
               | error SEP jprim_arglist { yyclearin; }
               | error NL {}
;

jprim_argument:  ATPOINT { curr_jprim->Set_Part(_SET_,_ATPOINT); }
	       | INLINE { curr_jprim->Set_Part(_SET_,_INLINE); }
               | INPLANE { curr_jprim->Set_Part (_SET_,_INPLANE); }
	       | ORIENTATION { curr_jprim->Set_Part (_SET_,_ORIENTATION);}
               | PARALLEL_AXES { curr_jprim->Set_Part (_SET_,_PARALLEL_AXES);}
	       | PERPENDICULAR { curr_jprim->Set_Part (_SET_,_PERPENDICULAR);}
               | I EQ NUM {curr_jprim->Set_Part (I, (Id) $3);}
               | J EQ NUM {curr_jprim->Set_Part (J, (Id) $3);}
;

//------------------------------- LSE START ---------------------------------

lse_arglist:   	lse_argument NL {}
              | lse_argument SEP lse_arglist {}
	      | error SEP lse_arglist { yyclearin; }
	      | error NL {}
;

lse_argument:   LSE_X EQ NUM {}
              | LSE_A EQ NUM {}
	      | STATIC_HOLD {}
	      | IC EQ NUM {}
	      | LSE_U EQ NUM {}
	      | LSE_B EQ NUM {}
	      | LSE_Y EQ NUM {}
	      | D EQ NUM {}
	      | LSE_C EQ NUM {}
;

//----------------------------- MARKER START --------------------------------

mark_arglist:     SEP marker_arglist {}
                | NL {}
;

marker_arglist:   marker_argument NL {}
                | marker_argument SEP marker_arglist {}
                | error SEP marker_arglist { yyclearin; }
                | error NL {}

marker_argument:  PART EQ NUM { curr_marker->Set_Part(PART,(Id) $3); }
		| POINT_MASS EQ NUM { 
		       curr_marker->Set_Part(POINT_MASS,(Id)$3); }
		| FLEX_BODY EQ NUM { curr_marker->Set_Part(FLEX_BODY,(Id) $3);}
		| NODE_ID EQ NUM { curr_marker->Set_Part(NODE_ID,(Id) $3); }
		| qp { curr_marker->Set_Part(QP,Service_vector3); }
		| zp { curr_marker->Set_Part(ZP,Service_vector3); }
		| xp { curr_marker->Set_Part(XP,Service_vector3); }
		| reuler { curr_marker->Set_Part(REULER,Angle_Vector); }
		| FLOATING { curr_marker->_Floating=Y; }
		| USEXP { curr_marker->_UseXP=Y; }
;

//---------------------------- MATERIAL START -------------------------------

material_arglist:   mat_argument NL {}
                  | mat_argument SEP material_arglist {}
                  | error SEP material_arglist { yyclearin; }
                  | error NL {}
;

mat_argument:       YOUNGS_MODULUS EQ NUM {
                    curr_material->Set_Part(YOUNGS_MODULUS,$3);}
                  | DENSITY EQ NUM {
		    curr_material->Set_Part(DENSITY,$3);}
		  | NAME EQ ID {
		    curr_material->Set_Part(NAME,(char*) $3);}
		  | POISSONS_RATIO EQ NUM {
		    curr_material->Set_Part(POISSONS_RATIO,$3); }
;

//----------------------------- OUTPUT ARGLIST -----------------------------//

output_arglist:     output_argument NL {}
                  | output_argument SEP output_arglist {}
		  | error SEP output_arglist {}
		  | error NL {}
;

output_argument:    REQSAVE {}
                  | GRSAVE {}
		  | GR521SAV {}
		  | CHART {}
		  | NOPRINT {}
		  | NOSEPARATOR {}
		  | OSFORMAT {}
		  | FIXED {}
		  | TELETYPE {}
		  | VPR {}
		  | DSCALE EQ NUM {}
                  | NUM {}
		  | VSCALE EQ NUM {}
		  | ASCALE EQ NUM {}
		  | FSCALE EQ NUM {}
		  | DZERO EQ NUM {}
		  | VZERO EQ NUM {}
		  | AZERO EQ NUM {}
		  | FZERO EQ NUM {}
;



//------------------------------- PART START --------------------------------

ground:          GROUND { curr_part->Set_Part(GROUND,$<stringa>1); }
;
material:        MATERIAL EQ ID {curr_part->Set_Part(MATERIAL,$<stringa>3); }
;

part_arglist:   part_argument NL {}
               	| part_argument SEP part_arglist {}
                | error SEP part_arglist { yyclearin; i_arg=0; }
                | error NL { i_arg=0; }
;
part_argument:   MASS EQ NUM { curr_part->Set_Part(MASS,$3); }
		| CM EQ NUM   { curr_part->Set_Part(CM,(Id) $3);}
	      	| IM EQ NUM   { curr_part->Set_Part(IM,(Id) $3);}
               	| IP EQ NUM SEP NUM SEP NUM {
	       		Vec3 S($3,$5,$7); 
			curr_part->Set_Part(IP,S); }
               	| material {}
		| ground {}
	       	| VM EQ NUM   { curr_part->Set_Part(VM,(Id) $3);}
	       	| WM EQ NUM   { curr_part->Set_Part(WM,(Id) $3);}
	       	| VX EQ NUM   { curr_part->Set_Part(VX,$3);}
	       	| VY EQ NUM   { curr_part->Set_Part(VY,$3);}
	       	| VZ EQ NUM   { curr_part->Set_Part(VZ,$3);}
	       	| WX EQ NUM   { curr_part->Set_Part(WX,$3); }
	       	| WY EQ NUM   { curr_part->Set_Part(WY,$3); }
	       	| WZ EQ NUM   { curr_part->Set_Part(WZ,$3); }
	       	| qg { curr_part->Set_Part(QG,Service_vector3); }
	       	| zg { curr_part->Set_Part(ZG,Service_vector3); }
	       	| xg { curr_part->Set_Part(XG,Service_vector3); }
	       	| reuler { curr_part->Set_Part(REULER,Angle_Vector); }
	       	| add_ip_list { }
	       	| EXACT EQ exact_list { curr_part->ncoord=i_arg;
	           curr_part->Set_Part(EXACT,Coord_Vector);
		   i_arg=0; }
;

add_ip_list:     NUM SEP NUM SEP NUM { 
                   Vec3 S($1,$3,$5); 
		   Set_Param(curr_part->AddIp,S,curr_part->_AddIp); }
;

//---------------------------- POINT_MASS START ------------------------------

pointmass_arglist: pointmass_arg NL {}
		| pointmass_arg SEP pointmass_arglist {}
                | error SEP pointmass_arglist { yyclearin; i_arg=0; }
                | error NL { i_arg=0; }
;

pointmass_arg: 	MASS EQ NUM { curr_pointmass->Set_Part(MASS,$3);}
        	| CM EQ NUM   { curr_pointmass->Set_Part(CM,(Id) $3);}
	       	| VX EQ NUM   { curr_pointmass->Set_Part(VX,$3);}
	       	| VY EQ NUM   { curr_pointmass->Set_Part(VY,$3);}
	       	| VZ EQ NUM   { curr_pointmass->Set_Part(VZ,$3);}
	       	| qg { curr_pointmass->Set_Part(QG,Service_vector3); }
	       	| zg { curr_pointmass->Set_Part(ZG,Service_vector3); }
	       	| xg { curr_pointmass->Set_Part(XG,Service_vector3); }
	       	| reuler { curr_pointmass->Set_Part(REULER,Angle_Vector); }
	       	| EXACT EQ exact_list {
		   curr_pointmass->ncoord=i_arg;
	           curr_pointmass->Set_Part(EXACT,Coord_Vector);
                   i_arg=0; }
;

//------------------------------ RESULTS START -------------------------------

results_arglist:  error SEP results_arglist { yyclearin; }
                | results_argument NL {}
                | results_argument SEP results_arglist {}
		| error NL {}
                | NL {}
;
results_argument: NOACCELERATIONS {}
		| NOAPPLIEDFORCES {}
		| NODATASTRUCTURES {}
		| NODISPLACEMENTS {}
		| NOFLOATINGMARKERS {}
		| NOLINEAR {}
		| NOREACTIONFORCES {}
		| NOSYSTEMELEMENTS {}
		| NOTIRES {}
		| NOVELOCITIES {}
		| COMMENT EQ ID {}
		| FORMATTED {}
;

//-----------------------------  REQUEST START --------------------------------

request_arglist: request_argument NL {}
               | request_argument SEP request_arglist {}
	       | error SEP request_arglist { yyclearin; }
	       | error NL {}
;

request_argument: request_type {}
                | request_value {}
		| request_flist {}
		| request_function {}
		| request_comment {}
;

request_type:  DISPLACEMENT {}
               | VELOCITY {}
	       | ACCELERATION {}
	       | FORCE {}
;

request_value:  I EQ NUM {}
              | J EQ NUM {}
              | RM EQ NUM {}
;

request_flist:  request_farg {}
              | request_farg request_flist {}
;

request_farg:   F1 EQ ID INVSLASH {}
              | F2 EQ ID INVSLASH {}
	      | F3 EQ ID INVSLASH {}
	      | F4 EQ ID INVSLASH {}
	      | F5 EQ ID INVSLASH {}
	      | F6 EQ ID INVSLASH {}
	      | F7 EQ ID INVSLASH {}
	      | F8 EQ ID INVSLASH {}
;
request_function: FUNCTION EQ USER ID {} INVSLASH
;

request_comment: TITLE EQ result_clist {}
                 | COMMENT EQ idlist {}
;

result_clist: Rc {}
            | Rc COLON result_clist {}
;

Rc:          ID {}
;

//------------------------------ SFORCE START --------------------------------

sforce_arglist : sforce_arg NL {}
               | sforce_arg SEP sforce_arglist {}
	       | error SEP sforce_arglist { yyclearin; }
	       | error NL {}
;

sforce_arg:    I EQ NUM { curr_sforce->Set_Part(I,(Id) $3);}
             | J EQ NUM { curr_sforce->Set_Part(J,(Id) $3);}
	     | TRANSLATION { curr_sforce->Set_Part(_SET_,_TRANSLATION); }
	     | ROTATION { curr_sforce->Set_Part(_SET_,_ROTATION); }
	     | ACTIONONLY {curr_sforce->Set_Part(ACTIONONLY,Y);}
	     | function_user { curr_sforce->Store_Formula("FUNCTION",formula_text);
		               curr_sforce->Set_Part(FUNCTION,Service_value); }
;

//--------------------------- SPRINGDAMPER START -----------------------------

springdamper_arglist:  springdamper_arg NL {}
                     | springdamper_arg SEP springdamper_arglist {}
                     | error SEP springdamper_arglist { yyclearin; }
                     | error NL {}
;

springdamper_arg :     I EQ NUM { curr_springdamper->Set_Part(I, (Id) $3);}
                     | J EQ NUM { curr_springdamper->Set_Part(J, (Id) $3);}
                     | damper_mode {}
                     | stl_argument {}
                     | srl_argument {}
;

damper_mode : TRANSLATION { curr_springdamper->Set_Part(_SET_,_TRANSLATION);}
            | ROTATION { curr_springdamper->Set_Part(_SET_,_ROTATION);}
;

stl_argument: C0 EQ NUM { curr_springdamper->Set_Part(C0,$3);}
            | K0 EQ NUM { curr_springdamper->Set_Part(K0,$3);}
            | FORCE EQ NUM {curr_springdamper->Set_Part(FORCE,$3);}
            | LENGTH EQ NUM  {curr_springdamper->Set_Part(LENGTH,$3);}
;
srl_argument: CT EQ NUM { curr_springdamper->Set_Part (CT,$3);}
            | KT EQ NUM { curr_springdamper->Set_Part (KT,$3);}
            | TORQUE EQ NUM { curr_springdamper->Set_Part(TORQUE,$3);}
            | ANGLE EQ ang { curr_springdamper->Set_Part(ANGLE,Service_Angle);}
;


// ----------------------------- UNITS START ---------------------------------


units_arglist: units_argument NL {}
               | units_argument SEP units_arglist {}
	       | error SEP units_arglist {}
	       | error NL {}
;

units_argument: FORCE EQ force_type {}
              | MASS EQ mass_type {}
	      | LENGTH EQ length_type {}
	      | TIME EQ time_type {}
	      | SYSTEM EQ system_type {}
	      | UCF EQ NUM {}
;

force_type:   DYNE {}
            | KILOGRAM_FORCE {}
	    | KNEWTON {}
	    | KPOUND_FORCE {}
	    | NEWTON {}
	    | OUNCE_FORCE {}
	    | POUND_FORCE {}
;

mass_type:    GRAM {}
            | KILOGRAM {}
	    | KPOUND_MASS {}
	    | OUNCE_MASS {}
	    | POUND_MASS {}
	    | SLUG {}
;

length_type:  CENTIMETER {}
            | FOOT {}
	    | KILOMETER {}
	    | INCH {}
	    | METER {}
            | MILLIMETER {}
	    | MILE {}
;

time_type:   HOUR {}
           | MILLISECOND {}
	   | MINUTE {}
	   | SECOND {}
;

system_type: CGS {}
           | FPS {}
	   | IPS {}
	   | MKS {}
	   | NONE {}
;

// ----------------------------- VARIABLE START -------------------------------

variable_arglist :  variable_arg NL {}
                 |  variable_arg SEP variable_arglist {}
		 |  error SEP variable_arglist { yyclearin; }
		 |  error NL {}
;

variable_arg:    IC EQ NUM { curr_variable->Set_Part (IC,$3); }
               | function_user { curr_variable->Set_Part(FUNCTION,Service_value); }
;

// ------------------------------ VFORCE START -------------------------------

vforce_arglist:   vforce_arg NL {}
                | vforce_arg SEP vforce_arglist {}
		| error SEP vforce_arglist { yyclearin; }
		| error NL {}
;

vforce_arg:    I EQ NUM { curr_vforce->Set_Part(I,(Id) $3); }
             | JFLOAT EQ NUM { curr_vforce->Set_Part(JFLOAT,(Id) $3); }
	     | RM EQ NUM { curr_vforce->Set_Part(RM,(Id) $3); }
	     | function_user { curr_vforce->Store_Formula("FUNCTION",formula_text);
	                       curr_vforce->Set_Part(FUNCTION,Service_value);}
	     | vforce_magnitude {}
;

vforce_magnitude:  vforcem_arg {}
                 | vforcem_arg INVSLASH vforce_magnitude {}
;

vforcem_arg:   FX EQ magnitude { curr_vforce->Store_Formula("FX",formula_text);
                                 curr_vforce->Set_Part(FX,Service_value); }
             | FY EQ magnitude { curr_vforce->Store_Formula("FY",formula_text);
	                         curr_vforce->Set_Part(FY,Service_value); }
	     | FZ EQ magnitude { curr_vforce->Store_Formula("FZ",formula_text);
	                         curr_vforce->Set_Part(FZ,Service_value); }
;

// ---------------------------- VTORQUE START ------------------------------

vtorque_arglist:  vtorque_arg NL {}
                | vtorque_arg SEP vtorque_arglist {}
		| error SEP vtorque_arglist { yyclearin; }
		| error NL {}
;

vtorque_arg: I EQ NUM { curr_vtorque->Set_Part(I,(Id) $3); }
             | JFLOAT EQ NUM { curr_vtorque->Set_Part(JFLOAT,(Id) $3); }
	     | RM EQ NUM { curr_vtorque->Set_Part(RM,(Id) $3); }
	     | function_user { curr_vtorque->Store_Formula("FUNCTION",formula_text);
	                       curr_vtorque->Set_Part(FUNCTION,Service_value); }
	     | vtorque_magnitude {}
;

vtorque_magnitude: vtorquem_arg {}
                 | vtorquem_arg INVSLASH vtorque_magnitude {}
;

vtorquem_arg:  TX EQ magnitude { curr_vtorque->Store_Formula("TX",formula_text);
                                 curr_vtorque->Set_Part(TX,Service_value); }
	     | TY EQ magnitude { curr_vtorque->Store_Formula("TY",formula_text);
	                         curr_vtorque->Set_Part(TY,Service_value); }
	     | TZ EQ magnitude { curr_vtorque->Store_Formula("TZ",formula_text);
	                         curr_vtorque->Set_Part(TZ,Service_value); }
;

// * * * * * * * * * * GENERAL POURPOUSE GRAMMAR FORMS * * * * * * * * * * * //


idlist:         ID {}
              | ID idlist {}
;

reuler:		REULER EQ ang1 SEP ang2 SEP ang3 {}
;

ang1:           ang { Angle_Vector.Xangle=Service_Angle; }
;
ang2:           ang { Angle_Vector.Yangle=Service_Angle; }
;
ang3:           ang { Angle_Vector.Zangle=Service_Angle; }
;

ang:		NUM { Service_Angle.value=$<value>1;
                      Service_Angle.ref=RADIANS; }
		| NUM D { Service_Angle.value=$<value>1;
		                   Service_Angle.ref=DEGREE;}
;

qp:	        QP EQ NUM SEP NUM SEP NUM {
                      Service_vector3.Set ($3,$5,$7); }
;
zp:	        ZP EQ NUM SEP NUM SEP NUM {
                      Service_vector3.Set ($3,$5,$7); }
;
xp:             XP EQ NUM SEP NUM SEP NUM {
                      Service_vector3.Set ($3,$5,$7); }
;

qg:		QG EQ NUM SEP NUM SEP NUM  { 
                      Service_vector3.Set ($3,$5,$7); }
;
xg:		XG EQ NUM SEP NUM SEP NUM  { 
                      Service_vector3.Set ($3,$5,$7); }
;
zg:		ZG EQ NUM SEP NUM SEP NUM  { 
                      Service_vector3.Set ($3,$5,$7); }
;

alimit:    	ALIMIT EQ ang { }
;
_error:     	ERROR EQ NUM { Service_value=$3; }
;
maxit:	   	MAXIT EQ NUM { Service_value=$3; }
;
pattern:   	PATTERN EQ J_coordlist {}
;
tlimit:	   	TLIMIT EQ NUM { Service_value=$3; }
;

/* SISTEMA DI COORDINATE PER POINT_MASS E PART */

exact_list:	exact_arg {}
		| exact_arg COLON exact_list {}
;

exact_arg:	XCOORD { Coord_Vector[i_arg++]=_XC; }
		| YCOORD { Coord_Vector[i_arg++]=_YC; }
		| ZCOORD { Coord_Vector[i_arg++]=_ZC; }
		| PSI { Coord_Vector[i_arg++]=_PSIC; }
		| THETA { Coord_Vector[i_arg++]=_THETAC; }
		| PHI { Coord_Vector[i_arg++]=_PHIC; }
;

/* SISTEMA DI COORDINATE PER ITERAZIONE SULLO JACOBIANO*/

J_coordlist:   J_argument {}
               | J_argument COLON J_coordlist  {}
;
J_argument:     TRUE  { Bool_Vector[i_arg++]=Y;}
	      | FALSE { Bool_Vector[i_arg++]=N;}
;

/* VARIABILI PER LA FUNZIONE USER - FUNCTION , E PER LE FORZE */


function_user:  FUNCTION EQ USER '(' explist ')' { Service_value=0; }
              | FUNCTION EQ exp { Service_value=0; }
	      | FUNCTION EQ NUM { Service_value=$3; formula_text[0]='\0';}
;

magnitude:      exp { Service_value=0; }
              | NUM { Service_value=$1; formula_text[0]='\0';}
;

// * * * * * * * * * * * * *  EXPRESSION PARSING * * * * * * * * * * * * * * //


/* Per ora ne viene fatta solo la lettura, non vengono interpretate 
   per essere tradotte - see warnings.. */

exp:            function {}
	      | variable {}
	      | if_stmt {}
              | expl '+' expl {  }
	      | expl '-' expl {  }
	      | expl '*' expl {  }
	      | expl '/' expl {  }
	      | expl '^' expl {}
	      | expl '@' expl {}
              | '-' expl %prec NEG { }
	      | '(' expl ')' { }
	      | error {}
;

/* Con questo stratagemma ogni qualvolta viene identificato un numero qls
   non singolo dopo l'uguale viene passato al parser delle espressioni */
   
expl:           exp {}
              | NUM { }
;

if_stmt:      IF '(' expl ':' explist ')' {}
; 

function:     ID '(' explist ')' { }
;

variable:     ID { }
;

explist:     expl {}
           | expl SEP explist {}
;

// * * * * * * * * * * * * * * END OF EXPR PARSER * * * * * * * * * * * * *  //

%%

// - - - - - - - - - - - - FUNZIONI STATICHE X TOKENS - - - - - - - - - - -  //


// QUESTE FUNZIONI DEVONO RIMANERE QUA PERCHE' LE VARIABILI
// YYNTOKENS, YYTNAME, YYTOKNUM SONO DEFINITE SOLO ALL'INTERNO DI B++.CC
// ESSENDO STATIC
    
const char** Token_Table;

inline const char* Find_Token(int P)
{
     tok_idx=tokens.find(P);
     int idx= (*tok_idx).second;
     return Token_Table[idx];
}

void Init_Token(void)
{
   ntokens=YYNTOKENS;
   Token_Table=new const char*[ntokens];
   for (int i=0; i< YYNTOKENS; i++) {
      Token_Table[i]=new char[strlen(yytname[i])];
      tokens.insert(token_entry(yytoknum[i],i));
      if (Token_Table[i] != NULL) Token_Table[i]=yytname[i];
   }
   return;
}

// - - - - - - - - - - - -  FUNZIONE DI VISUALIZ - - - - - - - - - - - - - - //

void Verbose ()
{
   Display_deck(beams,curr_beam,cout);
   Display_deck(markers,curr_marker,cout);
   Display_deck(joints,curr_joint,cout);
   Display_deck(parts,curr_part,cout);
   Display_deck(materials,curr_material,cout);
   Display_deck(accgravs,curr_accgrav,cout);
   Display_deck(equilibriums,curr_equilibrium,cout);
   Display_deck(jprims,curr_jprim,cout);
   Display_deck(ics,curr_ic,cout);
   Display_deck(pointmasses,curr_pointmass,cout);
   Display_deck(springdampers,curr_springdamper,cout);
   Display_deck(sforces,curr_sforce,cout);
   Display_deck(gforces,curr_gforce,cout);
   Display_deck(vforces,curr_vforce,cout);
   Display_deck(vtorques,curr_vtorque,cout);
   Display_deck(variables,curr_variable,cout);
}

void Report_form (char* file_to_write)
{
   fclose (errfile);
   fclose (messagefile);
   fclose (tablefile);
   fclose (logfile);
   if (VERBOSE_MODE==Y) Verbose();
   cout << "END OF RECOGNITION" << endl << "LINES: ";
   cout << ncount << endl;
   cout << "ERRORS FOUND:" << nerr << endl;
   cout << "WARNINGS:" << nwarnings << endl;
   if (SPECIFY_FILE==Y) {
      cout << "FILE TO OUTPUT ? ";
      cin >> file_to_write;
   }
   cout << endl << "WRITING MBDYN CODE .." << endl;
   ofstream outfile(file_to_write,ios::out);
   if (!outfile) { 
     cout << "FILE CANNOT BE OPENED - FATAL ERROR" << endl;
     exit(-1);
   } else Translate(outfile);
   cout << "COMPLETED" << endl;
}


// ---------------------------- MAIN CODE ------------------------------------


int main (int argn, char* argv[])
{
   FILE* temp;
   Boolean To_process;
   temp=tmpfile();
   char *reportfile, *infile, *outfile;
   Init_Token();
   stdout=temp;
   reportfile = CheckInputFile(argn, argv, errfile, messagefile, tablefile);
   
   if (reportfile) 
     {
   	nerr=0; 
  	ncount=1;
	
	infile = new char[strlen(reportfile)+4+1];
	strcpy(infile, reportfile);
	strcat(infile, ".adm");

	FILE *test;
	test = fopen (infile,"r");
	if (test==NULL) {
		cout << endl << "FILE " << infile << " DOES NOT EXIST" << endl;
		exit(EXIT_FAILURE);
	}
	outfile = new char[strlen(reportfile)+4+1];
	strcpy(outfile, reportfile);
	strcat(outfile, ".mbd");
        cout << endl << "RECOGNITION START" << endl;
	cout << "FILE OPENED FOR SCANNING.." << endl;
   	yyin=fopen(infile,"r");
	yyparse();
     }
   else
     {
        cout << endl << "RECOGNITION START" << endl;
	cout << "CONSOLE OPENED FOR READING" << endl;
        nerr=0; ncount=1;
	yyin=stdin;
	errfile=fopen("console.err","w");
	messagefile=fopen("console.msg","w");
	tablefile=fopen("console.ref","w");
        yyparse();
	outfile="console.mbd";
     }
     Report_form (outfile);
     fclose(temp);
}
