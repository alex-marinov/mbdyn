/* Elemento aerodinamico modale */

/* Elemento aerodinamico modale by Felice Felippone(C) - Marzo 2000 
 * L'elemento aerodinamico modale e' collegato al modal joint.
 * La posizione dei punti aeodinamici e' definita dalle deformate
 * modali ( + il moto rigido) della struttura. L'utente deve specificare
 * il nodo modale (rigido) collegato e l'orientamento del corpo aerodinamico 
 * rispetto a questo. Per adesso si assume che dall'esterno qualcuno si
 * preoccupi di fornire all'elemento aerodinamico una tabella con le
 * 'zone di influenza' , vale a dire i nodi FEM che l'elemento aerodinamico
 * utilizzera' per ricavare la posizione di ogni punto aerodinamico.
 * NOTA: IN FASE DI SVILUPPO!!!
*/

#include "aerodyn.h"
#include "rotor.h"
#include "beam.h"
#include "modal.h"

#include "gauss.h"
#include "aerod2.h"
#include "shape.h"

/* per l'output !?! */
#define AEROD_OUT_STD 0
#define AEROD_OUT_PGAUSS 1
#define AEROD_OUT_NODE 2

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
typedef struct {
   Vec3 f;
   doublereal alpha;
} Aero_output;
#endif /* AEROD_OUTPUT */



/* AerodynamicModal - begin */


class AerodynamicModal 
: virtual public Elem, public AerodynamicElem, public InitialAssemblyElem, 
public DriveOwner {
 protected:
   const StructNode* pModalNode; /* Nodo modale per il moto rigido */
   const Modal* pModalJoint;     /* puntatore all'elemento modale di riferimento */
   Rotor* pRotor;
   
   const Mat3x3 Ra;             /* Rotaz. del sistema aerodinamico al nodo */
   const Vec3 Ra3;              /* Terza colonna della matrice Ra */
   
   const ShapeOwner Chord;         /* corda */
   const ShapeOwner ForcePoint;    /* punto di applicazione della forza (1/4) */
   const ShapeOwner VelocityPoint; /* punto di applicazione b.c. (3/4) */
   const ShapeOwner Twist;         /* svergolamento */

   double dL;                   /* lunghezza dell'elemento */
   
   integer iInst;               /* switch di instazionarieta': 0/1/2 */  
   integer iProfile;            /* Tipo di profilo */
   
   integer NModes;              /* Numero di modi */
   integer NFemNodes;           /* Numero di nodi FEM */
   integer NAeroElems;          /* Numero di elementi aerodinamici */
   unsigned int* iAeroTable;    /* Tabella con le connessioni nodi FEM <-> aero */
   //integer Zaxis;               
   
   Mat3xN* pPHIt;               /* Autovettori (saranno sostituiti dalla matrice int. H) */
   Mat3xN* pPHIr;
   Mat3xN* pFemNodesPosition;   /* Posizione dei nodi FEM */

   FullMatrixHandler* pGTKG;     /* Matrice di interpolazione tra nodi aerodinamici e strutturali */    
   FullMatrixHandler* pH;        /* Matrice che interpola gli spostamenti aerodinamici con
                                   quelli modali */  
   /* */

   GaussDataIterator GDI;       /* Iteratore sui punti di Gauss */
   doublereal* pdOuta;          /* Dati privati */
   doublereal** pvdOuta;        /* */
   
   
   Vec3 F;                      /* Forza */
   Vec3 M;                      /* Momento */

   VecN a;                      /* deformate modali */
   VecN aPrime;

   Mat3x3** ppR1tot;            

#if AEROD_OUTPUT == AEROD_OUT_PGAUSS
   /* temporaneo, per output */
   Aero_output* pOutput;
#endif /* AEROD_OUTPUT */

   /* overload della funzione di ToBeOutput();
    * serve per allocare il vettore dei dati di output se il flag
    * viene settato dopo la costruzione */
   virtual void SetOutputFlag(flag f = flag(1));   
   
   /* Assemblaggio residuo */
   void AssVec(SubVectorHandler& WorkVec);
   
 public:
   AerodynamicModal(unsigned int uLabel, 
		   const StructNode* pN, const Modal* pModalJoint, Rotor* pR,
		   const Mat3x3& RaTmp,
		   const Shape* pC, const Shape* pF, 
		   const Shape* pV, const Shape* pT,
		   integer iI, integer iN, integer iP, 
                   integer iM, integer NFN,
		   integer iAP,
                   Mat3xN* pModeShapest, Mat3xN* pModeShapesr,
                   FullMatrixHandler* pH, FullMatrixHandler* pGTKG, double dL,
                   Mat3xN* pFemNodesPosition, //integer Zaxis,
		   const DriveCaller* pDC, flag fOut);
   virtual ~AerodynamicModal(void);
   
   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const;
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const {
      return Elem::AERODYNAMIC;
   };   
   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 6;
      *piNumCols = 1;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ ,
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	DEBUGCOUTFNAME("AerodynamicModal::AssJac");
	WorkMat.SetNullMatrix();
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal dCoef,
				    const VectorHandler& XCurr,
				    const VectorHandler& XPrimeCurr);

   /* output; si assume che ogni tipo di elemento sappia, attraverso
    * l'OutputHandler, dove scrivere il proprio output */
   virtual void Output(OutputHandler& OH) const;

   /* Numero di GDL iniziali */
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 0;
   };
   
   /* Dimensioni del workspace */
   virtual void InitialWorkSpaceDim(integer* piNumRows,
				    integer* piNumCols) const {
      *piNumRows = 6;
      *piNumCols = 1;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,	  
	    const VectorHandler& /* XCurr */) {
	DEBUGCOUTFNAME("AerodynamicModal::InitialAssJac");
	WorkMat.SetNullMatrix();
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
					   const VectorHandler& XCurr);
   
   
   /* Tipo di elemento aerodinamico */
   virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
      return AerodynamicElem::AERODYNAMICMODAL;
   };
   
};

/* Funzioni di interpolazione */
extern doublereal ShapeFunc2N(doublereal d, integer iNode);
extern doublereal DxDcsi2N(const Vec3& X1, const Vec3& X2);


/* AerodynamicModal - end */



class DataManager;
class MBDynParser;

extern Elem* ReadAerodynamicModal(DataManager* pDM,
				  MBDynParser& HP,
				  unsigned int uLabel);
