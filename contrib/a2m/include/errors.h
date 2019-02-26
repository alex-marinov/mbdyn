/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


#ifndef ERRORS_H
#define ERRORS_H

char* Error_code[] = {
   "UNRECOGNIZED CARD",
   "DUPLICATE CARD ENTRY",
   "UNKNOWN ENTRY (SYNTAX ERROR)",
   "INTERNAL ERROR, UNABLE TO SET",
   "INTERNAL ERROR, DATA MISMATCH",
   "DUPLICATE LABEL (ALREADY USED)",
   "BEAM, A MARKER WAS NOT DEFINED",
   "BEAM, MISSING LENGTH",
   "BEAM, MISSING INERTIA PROPERTY",
   "BEAM, MISSING MATERIAL PROPERTY",
   "PART, MARKER CM UNDEFINED BUT MASS NOT NULL",
   "PART, ENTRY NOT CONGRUENT WITH NULL MASS",
   "MARKER, ENTRY NOT CONGRUENT WITH NODE_ID METHOD",
   "MARKER, ENTRY NOT CONGRUENT WITH PART/P.MASS METHOD",
   "MARKER, DEFINITION IN FLOATING MODE DOES NOT NEED",
   "CONFUSED ABOUT BEAM DAMPING MODE: BOTH CMATRIX AND CRATIO DEFINED",
   "MARKER, EULER ANGLES METHOD DOES NOT NEED",
   "REFERENCE USED DOES NOT EXIST",
   "MARKER, MISSING PARAMETER FOR NODE_ID METHOD",
   "JOINT, A MARKER ID WAS NOT DEFINED",
   "JOINT, EXTENSIONS NOT ALLOWED FOR THIS TYPE OF JOINT",
   "JOINT, CYLINDRICAL JOINT DOES NOT NEED",
   "JOINT, RACKPINT JOINT DOES NOT NEED",
   "JOINT, REVOLUTE JOINT DOES NOT NEED",
   "JOINT, SCREW JOINT DOES NOT NEED",
   "JOINT, TRANSLATIONAL JOINT DOES NOT NEED",
   "JOINT, EXTENSIONS MUST BE SPECIFIED FOR THIS JOINT",
   "JPRIM, A MARKER ID WAS NOT DEFINED",
   "JPRIM, MISSING PRIMITIVE TYPE",
   "POINTMASS, REULER METHOD DOES NOT NEED",
   "GFORCE, FUNCTION METHOD DOES NOT NEED",
   "GFORCE, MISSING PARAMETER",
   "SFORCE, A MARKER WAS NOT DEFINED",
   "SFORCE, UNDEFINED MODE (TRANSLATION OR ROTATION?)",
   "SPRINGDAMPER, A MARKER WAS NOT DEFINED",
   "SPRINGDAMPER, ENTRY NOT CONGRUENT WITH TRANSLATION MODE",
   "SPRINGDAMPER, ENTRY NOT CONGRUENT WITH ROTATION MODE",
   "VFORCE, MARKER RM IS NEEDED BUT IT'S NOT DEFINED!",
   "VFORCE, MARKER I IS NEEDED BUT IT'S NOT DEFINED!",
   "VFORCE, MARKER J(FLOAT) IS NEEDED BUT IT'S NOT DEFINED!",
   "VFORCE, USE OF FORCE EXPRESSION AND FUNCTION IS AMBIGOUS",
   "VFORCE, AT LEAST A EXPRESSION OR A FUNCTION IS NEEDED",
   "VTORQUE, MARKER RM IS NEEDED BUT IT'S NOT DEFINED!",
   "VTORQUE, MARKER I IS NEEDED BUT IT'S NOT DEFINED!",
   "VTORQUE, MARKER J(FLOAT) IS NEEDED BUT IT'S NOT DEFINED!",
   "VTORQUE, USE OF FORCE EXPRESSION AND FUNCTION IS AMBIGOUS",
   "VTORQUE, AT LEAST A EXPRESSION OR A FUNCTION IS NEEDED",
   "GROUND PART ALREADY DEFINED"
};

#endif
