
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

int NO_STRUCTURAL_NODES, NO_FORCES, NO_JOINTS, NO_RIGIDBODIES;
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
   /* 04. - Beam, translation section */
   Translate_deck ("Translating beams...",beams,curr_beam,out);

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
   NO_STRUCTURAL_NODES=MBNodes.size();
   NO_FORCES=MBForces.size();
   NO_JOINTS=MBJoints.size();
   NO_RIGIDBODIES=MBBodies.size();
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
     << "  initial time: 0.;" << endl
     << "  final time: 1.;" << endl
     << "  time step: .1;" << endl
     << endl
     << "  set: real rho = .6;" << endl
     << "  method: ms, rho, rho;" << endl
     << "  tolerance: 1.e-6;" << endl
     << "  max iterations: 20;" << endl
     << endl
     << "  derivatives coefficient: 1.;" << endl
     << "  derivatives tolerance: 1.e-9;" << endl
     << "  derivatives max iterations: 10;" << endl
     << "end: multistep;" << endl << endl;
   /* 03. - Control data section */
   out << "begin: control data;" << endl
     << "  structural nodes: " << NO_STRUCTURAL_NODES<< ";" << endl
     << "  rigid bodies: " << NO_RIGIDBODIES << ";" << endl
     << "  joints: " << NO_JOINTS << ";" << endl
     << "  forces: " << NO_FORCES << ";" << endl
     << endl
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
   Restart_MBDYNDeck (MBReference,out);
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
