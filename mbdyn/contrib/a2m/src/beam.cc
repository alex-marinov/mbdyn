//                                 BEAM.CC

//Definizioni per la card Beam

#include <beam.h>

extern MBDyn_deck MBBeams;
extern MBDyn_deck MBReference;
extern MBDyn_deck MBNodes;
extern MBDyn_deck MBBodies;
extern MBDyn_deck MBConstLaw;

s_beam::s_beam(void)	: Asy(0), Asz(0), Cratio(0), Ixx(0),Iyy(0),Izz(0),
		        Area(0),Length(0),Midnode(0),Emodulus(0),Gmodulus(0),
			Cmatrix (new double[21]),
                    	_Node1(N), _Node2(N), _Asy(N), _Asz(N), _Cratio(N),
                        _Ixx(N), _Iyy(N), _Izz(N), _Area(N), _Length(N), 
                        _Midnode(N), _Emodulus(N), _Gmodulus(N),
                        _Cmatrix(N) { Node[0]=0; Node[1]=0; }


s_beam::~s_beam(void)     {}


Boolean s_beam::Test()
{
   const int err_before=nerr;
   // CONTROLLA CHE SIANO PRESENTI I PARAMETRI NECESSARI
   if (_Node1==N) { out_error(6,"OF ID I"); }
   if (_Node2==N) { out_error(6,"OF ID J"); }
   if (_Length==N) { out_error(7,""); }
   if (_Ixx==N) { out_error(8,"IXX"); }
   if (_Iyy==N) { out_error(8,"IYY"); }
   if (_Izz==N) { out_error(8,"IZZ"); }
   if (_Area==N) { out_error(8,"AREA"); }
   if (_Emodulus==N) { out_error (9,"YOUNG MODULUS"); }
   if (_Gmodulus==N) { out_error (9,"SHEAR MODULUS"); }
   if ((_Cmatrix==Y) & (_Cratio==Y)) { out_error (15,""); }
   if (err_before != nerr) return Y; else return N;
}

inline const char* const s_beam::Gettype(void) const
{
   return "BEAM";
}

ostream& s_beam::Print (ostream& out) const
{
   out << endl;
   out << "BEAM:" << label << endl ;
   out << "     " << "Node 1 [" << _Node1 << "] = "  << Node[0] << endl
       << "     " << "Node 2 [" << _Node2 << "] = "  << Node[2] << endl;
   out << "     " << "Midnode: " << Midnode << endl;
   out << "     " << "Length [" << _Length << "] = " << Length << endl;
   out << "     " << "Area [" << _Area << "] = " << Area << endl;
   out << "     " << "Asy [" << _Asy << "] = " << Asy << endl
       << "     " << "Asz [" << _Asz << "] = " << Asz << endl;
   out << "     " << "Ixx [" << _Ixx << "] = " << Ixx << endl
       << "     " << "Iyy [" << _Iyy << "] = " << Iyy << endl
       << "     " << "Izz [" << _Izz << "] = " << Izz << endl;
   out << "     " << "Emodulus [" << _Emodulus << "] = " << Emodulus << endl
       << "     " << "Gmodulus [" << _Gmodulus << "] = " << Gmodulus << endl
       << "     " << "Cratio [" << _Cratio << "] = " << Cratio << endl
       << "     " << "Cmatrix [" << _Cmatrix << "]" << endl;
   if (_Cmatrix==Y) {
      out << "     ";
      for (int i=0; i< 8; i++) out << Cmatrix[i] << " ";
      out << endl << "     " ;
      for (int i=8; i< 16; i++) out << Cmatrix[i] << " " ;
      out << endl << "     ";
      for (int i=16; i<21; i++) out << Cmatrix[i] << " " ;
      out << endl;
   }
   out << endl;
   return out;
}

void s_beam::Translate(ostream& out)
{
   char* comment = new char [160];
   // Traduzione : TRAVE ADAMS = STRUCT NODE + BEAM MBDYN
   pCL6D pK1=new ConstitutiveLaw6D(Emodulus,Gmodulus);
   pCL6D pK2=new ConstitutiveLaw6D(Emodulus,Gmodulus);
   /* Dati da inserire nella formula della BEAM */
   Vec3 F1,F2,F3;
   Mat3x3 R1,R2;
   
   /* INSERIMENTO DEL NODO STRUTTURALE A SINISTRA */
   Id Node1 = GetFreeLabel (MBNodes,Node[0]);
   MBDyn_node_structural* LEFTNODE = new MBDyn_node_structural (Node1);
   /* Node[0] contiene l'ID del Marker I, che in MBDYN, data la corrispon-
    * denza imposta è uguale al reference I */
   NodeSet (LEFTNODE,Node[0]);
   MBNodes.insert(MBDyn_entry(Node1,(MBDyn_card*) LEFTNODE));

   for (int i=0;i<160;i++) comment[i]=' ';
   sprintf (comment,"\n# Structural Node [Left] %i relative to Adams BEAM %i",
	    Node1,label);
   LEFTNODE->Remark(comment);
   
   /* INSERIMENTO DEL NODO STRUTTURALE A DESTRA */
   Id Node3 = GetFreeLabel (MBNodes,Node[2]);
   MBDyn_node_structural* RIGHTNODE = new MBDyn_node_structural (Node3);
   /* Node[2] contiene l'ID del Marker J, che in MBDYN, data la corrispon-
    * denza imposta è uguale al reference J */
   NodeSet (RIGHTNODE,Node[2]);
   MBNodes.insert(MBDyn_entry(Node3,(MBDyn_card*) RIGHTNODE));

   for (int i=0;i<160;i++) comment[i]=' ';
   sprintf (comment,"\n# Structural Node [Right] %i relative to Adams BEAM %i",
	    Node3,label);
   RIGHTNODE->Remark(comment);

   /* INSERIMENTO DEL NODO STRUTTURALE CENTRALE */
   Id Marker_to_Node2 = GetFreeLabel (MBReference);
   /* I marker I e J esistono in quanto definiti come references. */
   MBDyn_reference* TMPREF2 = new MBDyn_reference(Marker_to_Node2);
   /* Interpolazione tra i marker I e J per ottenere il nodo centrale */
   Ref_Interp (TMPREF2,Node[0],Node[2]);
   Id Node2 = GetFreeLabel (MBNodes);
   MBDyn_node_structural* CENTRALNODE = new MBDyn_node_structural (Node2);
   NodeSet (CENTRALNODE,TMPREF2);
   MBNodes.insert(MBDyn_entry(Node2,(MBDyn_card*) CENTRALNODE));

   for (int i=0;i<160;i++) comment[i]=' ';
   sprintf (comment,"\n# Structural Node [Center] %i relative to Adams BEAM %i",
	    Node2,label);
   CENTRALNODE->Remark(comment);
   
   /* Distruzione del marker temporaneo */
   delete TMPREF2;
   // Creazione della beam attaccata al nodo centrale.
   // Al solito verifica che ci sia una label libera per la Beam.
   Id IDBEAM = GetFreeLabel (MBBeams,label);
   MBDyn_beam* CENTRALBEAM = new
     MBDyn_beam (IDBEAM,Node1,Node2,Node3,F1,F2,F3,R1,R2,pK1,pK2);
   MBBeams.insert(MBDyn_entry(IDBEAM,(MBDyn_card*) CENTRALBEAM));

   for (int i=0;i<160;i++) comment[i]=' ';
   sprintf (comment,"\n# MBDyn Beam %i relative to Adams BEAM %i",
	    IDBEAM,label);
   CENTRALBEAM->Remark(comment);

   // CREAZIONE DEI CORPI RIGIDI - Strutture dati che contengono informazioni
   // 
   // 3 Rigid body MBDyn con l'equivalente massa
   // per la scelta dei punti interpoli a 1radice di 3 (punti di Gauss)
   // e calcoli la corrispondente massa
   //
   /* QUESTA PARTE DI CODICE HA SENSO SOLO SE LA BEAM DI ADAMS HA MASSA */
   /* MA ALLO STATO ATTUALE TALE ELEMENTO NE E' PRIVO DI CONSEGUENZA NON*/
   /* OCCORRE SPECIFICARE ALCUNCHE' DI AGGIUNTIVO                       */
   /*
   Id LB[3];
   double M[3];
   Vec3 XGC[3];
   Mat3x3 J0[3];
   // La beam di MBDYN non ha massa - neanche quella di ADAMS
   // di conseguenza:
   for (int i=0;i<3;i++) M[i]=0;
   // occorre stabilire la posizione dei centri di massa delle parti
   LB[0]=GetFreeLabel (MBBodies);
   LB[1]=GetFreeLabel (MBBodies);
   LB[2]=GetFreeLabel (MBBodies);
   MBDyn_body* MB1=new MBDyn_body (LB[0],Node1,M[0],XGC[0],J0[0]);
   MBDyn_body* MB2=new MBDyn_body (LB[1],Node2,M[1],XGC[1],J0[1]);
   MBDyn_body* MB3=new MBDyn_body (LB[2],Node3,M[2],XGC[2],J0[2]);
   MBBodies.insert(MBDyn_entry(LB[0],(MBDyn_card*) MB1));
   MBBodies.insert(MBDyn_entry(LB[1],(MBDyn_card*) MB2));
   MBBodies.insert(MBDyn_entry(LB[2],(MBDyn_card*) MB3));
   */
   return;
}

