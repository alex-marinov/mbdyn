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


// TRANSLATE


// Funzioni di conversione da Adams a Mbdyn
// GESTORE DELLE CHIAMATE DELLE ROUTINE DI TRASFORMAZIONE

#include <translate.h>
#include <storage.h>
#include <carddecl.h>

p_card_entry i;
p_MTP_entry iMTP;
MTP Tempbuffer;
extern MTP_deck Marker_Table;
extern PTR_deck Part_Table;
extern BTR_deck Beam_reference_Table;

extern MBDyn_deck MBReference;
extern MBDyn_deck MBBeams;
extern MBDyn_deck MBBodies;
extern MBDyn_deck MBNodes;
extern MBDyn_deck MBForces;
extern MBDyn_deck MBAccs;
extern MBDyn_deck MBJoints;

MBDyn_reference* buffer;
MBDyn_beam* Bbuffer;
MBDyn_body* BBody;
MBDyn_node* BNode;
MBDyn_node* BPart;
MBDyn_force* BForce;
MBDyn_joint* BJoint;
MBDyn_gravity* BGravity;

p_MBDyn_entry pbuffer;

extern card_deck markers;
extern card_deck beams;
extern card_deck joints;
extern card_deck parts;
extern card_deck pointmasses;
extern card_deck sforces;
extern card_deck gforces;
extern card_deck vtorques;
extern card_deck vforces;
extern card_deck accgravs;

extern s_beam* curr_beam;
extern s_marker* curr_marker;
extern s_joint* curr_joint;
extern s_part* curr_part;
extern s_pointmass* curr_pointmass;
extern s_sforce* curr_sforce;
extern s_gforce* curr_gforce;
extern s_vtorque* curr_vtorque;
extern s_vforce* curr_vforce;
extern s_accgrav* curr_accgrav;

extern double A2M_initial_time;
extern double A2M_final_time;
extern double A2M_time_step;
extern double A2M_rho;
extern char A2M_method [];
extern double A2M_tolerance;
extern double A2M_max_iterations;
extern double A2M_deriv_coef;
extern double A2M_deriv_tol;
extern double A2M_deriv_iter;

int NO_STRUCTURAL_NODES, NO_FORCES, NO_JOINTS, NO_RIGIDBODIES, NO_BEAMS;
double IS,IT;

void Translate (ofstream& out)
{
   //
   /*  T R A N S L A T I O N    O F   F I L E */
   //

   /* A.  - Marker, creating queue */
   for (i=markers.begin();i!=markers.end();i++) {
      curr_marker=(s_marker*) (*i).second;
      curr_marker->queue();
   }
   /* end of queueing operation */

   /* 03. - Part, pointmassess translation section */
   Translate_deck ("Translating parts...",parts,curr_part,out);
   Translate_deck ("Translating pointmasses...",pointmasses,curr_pointmass,out);

   /* B.  - Dequeueing markers . */
   for (i=markers.begin();i!=markers.end();i++) {
      curr_marker=(s_marker*) (*i).second;
      curr_marker->dequeue();
   }
   /* End of dequeueing .. */

   /* 05. - Marker, translation section */
   Translate_deck ("Translating markers...",markers,curr_marker,out);
   /* 06. - Forces, translation section */
   
   Translate_deck ("Translating sforces...",sforces,curr_sforce,out);
   Translate_deck ("Translating gforce...",gforces,curr_gforce,out);
   Translate_deck ("Translating vtorques...",vtorques,curr_vtorque,out);
   Translate_deck ("Translating vforces...",vforces,curr_vforce,out);
   
   /* 07. - Gravity, translation section */
   Translate_deck ("Translating accgrav...",accgravs,curr_accgrav,out);
   /* 08. - Joints, translation section */
   Translate_deck ("Translating joints...",joints,curr_joint,out);
   /* Set the total amount of elements for each type */
   /* 09. - Beam, translation section */
   Translate_deck ("Translating beams...",beams,curr_beam,out);

   NO_STRUCTURAL_NODES=MBNodes.size();
   NO_FORCES=MBForces.size();
   NO_JOINTS=MBJoints.size();
   NO_RIGIDBODIES=MBBodies.size();
   NO_BEAMS=MBBeams.size();
   cout << "Processing output file ... " << endl;
   //
   /* 01. - HEADER OF THE TRANSLATED FILE */
   //
   out << endl;
   out << "# File generato da Adams2MBDyn - versione 1.0 Beta" << endl
     << endl;
   /* 02. - Integrator type , data section */
   out << "begin: data;" << endl
     << "  integrator: multistep;" << endl
     << "end: data;" << endl << endl;
   /* 02a.- Multistep definitions */
   out << "begin: multistep;" << endl
     << "  initial time: " << A2M_initial_time << ";" << endl
     << "  final time: " << A2M_final_time << ";" << endl
     << "  time step: " << A2M_time_step << ";" << endl
     << endl
     << "  set: real rho = " << A2M_rho << ";" << endl
     << "  method: " << A2M_method << ";" << endl
     << "  tolerance: " << A2M_tolerance << ";" << endl
     << "  max iterations: " << A2M_max_iterations << "; "<< endl
     << endl
     << "  derivatives coefficient: " << A2M_deriv_coef << ";" << endl
     << "  derivatives tolerance: " << A2M_deriv_tol << "; " << endl
     << "  derivatives max iterations: " << A2M_deriv_iter << ";" << endl
     << "end: multistep;" << endl << endl;
   /* 03. - Control data section */
   out << "begin: control data;" << endl
     << "  structural nodes: " << NO_STRUCTURAL_NODES<< ";" << endl
     << "  rigid bodies: " << NO_RIGIDBODIES << ";" << endl
     << "  joints: " << NO_JOINTS << ";" << endl
     << "  forces: " << NO_FORCES << ";" << endl
     << "  beams: " << NO_BEAMS << ";" << endl;
   
   /* Se è presente la gravità la inserisce nella lista */
   if (accgravs.size()!=0) out << "  gravity;" << endl;
   
   out << endl
     << "  initial stiffness: " << IS << ";" << endl
     << "  initial tolerance: " << IT << ";" << endl
     << endl
     << "#  skip initial joint assembly;" << endl
     << "#  use: rigid bodies, in assembly;" << endl
     << "#  make restart file;" << endl
     << "end: control data;" << endl << endl;
   //
   /*  R E S T A R T   S U L  F I L E   D I   O U T P U T  */
   //
   /* REFERENCE */
   out << "# input related card - Reference section" << endl << endl;
   
   /* Riordina i reference affinché i BCS compaiono prima degli altri marker */
   Id idx;
   Id counter=1;
   MBDyn_reference* TEMPREF;
   MBDyn_deck MBMarkers;
   p_MBDyn_entry index;
   p_PTR_entry indexPTR;
   for (indexPTR=Part_Table.begin();indexPTR!=Part_Table.end();indexPTR++)
     {
	idx = (*indexPTR).second;
	TEMPREF = (MBDyn_reference*) Find_MBCard (idx,MBReference);
	MBMarkers.insert (MBDyn_entry(counter++,(MBDyn_card*) TEMPREF));
     }
   p_BTR_entry indexBTR;
   for (indexBTR=Beam_reference_Table.begin();
	indexBTR!=Beam_reference_Table.end();
	indexBTR++)
     {
	idx = (*indexBTR).second;
	TEMPREF = (MBDyn_reference*) Find_MBCard (idx,MBReference);
	MBMarkers.insert (MBDyn_entry(counter++,(MBDyn_card*) TEMPREF));
     }
   p_MTP_entry indexMTP;
   for (indexMTP=Marker_Table.begin();indexMTP!=Marker_Table.end();indexMTP++)
     {
	idx = (*indexMTP).first;
	TEMPREF = (MBDyn_reference*) Find_MBCard (idx,MBReference);
	MBMarkers.insert (MBDyn_entry(counter++,(MBDyn_card*) TEMPREF));
     }
   if (MBReference.size()!=MBMarkers.size())
     cout << "[DEBUG] : REORDERED MARKERS DOES NOT LINK WITH ORIGINAL DECK" << endl;
   /* End of reordering */
   
   Restart_MBDYNDeck (MBMarkers,out);
   /* NODE SECTION */
   out << "begin: nodes;" << endl;
   Restart_MBDYNDeck (MBNodes,out);
   out << "end: nodes;" << endl << endl;
   /* ELEMENT SECTION */
   out << "begin: elements;" << endl;
   Restart_MBDYNDeck (MBBeams,out);
   Restart_MBDYNDeck (MBBodies,out);
   Restart_MBDYNDeck (MBForces,out);
   Restart_MBDYNDeck (MBAccs,out);
   Restart_MBDYNDeck (MBJoints,out);
   out << endl;
   out << "end: elements;" << endl;
   out.close();
}
