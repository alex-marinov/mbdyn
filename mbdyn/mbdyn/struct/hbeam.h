/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
 */

/* Trave a volumi finiti, con predisposizione per forze di inerzia consistenti
 * e legame cositutivo piezoelettico */


#ifndef HBEAM_H
#define HBEAM_H

#include "myassert.h"
#include "except.h"

#include "strnode.h"
#include "elem.h"

#include "constltp.h"

/* HBeam - begin */

class HBeam 
: virtual public Elem, public InitialAssemblyElem {
    friend class AerodynamicBeam;

  public:
    class ErrGeneric {};
   
  protected:   
    enum NodeName { NODE1 = 0, NODE2 = 1, NUMNODES = 2 };
    enum Deformations { STRAIN = 0, CURVAT = 1, NUMDEFORM = 2 };
   
    /* Puntatori ai nodi */
    const StructNode* pNode[NUMNODES];
   
    /* Offset dei nodi */
    Vec3 f[NUMNODES];
    Vec3 fRef[NUMNODES];
   
    /* Matrice di rotazione delle sezioni - non sono const perche' vengono
     * aggiornate ad ogni iterazione */
    Mat3x3 R;
    Mat3x3 RRef;

    /* Constitutive laws*/
    ConstitutiveLaw6DOwner* pD;
   
    /* Reference constitutive laws */
    Mat6x6 DRef;
   
    /* Velocita' angolare delle sezioni */
    Vec3 Omega;
    Vec3 OmegaRef;
   
    /* Dati temporanei che vengono passati da AssRes ad AssJac */
    Vec6 Az;
    Vec6 AzRef;
    Vec6 AzLoc;
    Vec6 AzLocRef;
    Vec6 DefLoc;
    Vec6 DefLocRef;

    Vec3 p;   
    Vec3 g;   
    Vec3 L0;   
    Vec3 L;   
       
    Vec3 LRef;   
   
    doublereal dsdxi;
   
    /* Is first res? */
    flag fFirstRes;
         
    virtual Vec3 
    InterpState(const Vec3& v1,
                const Vec3& v2);
    virtual Vec3 
    InterpDeriv(const Vec3& v1,
                const Vec3& v2);
   
    /* Funzioni di calcolo delle matrici */
    virtual void 
    AssStiffnessMat(FullSubMatrixHandler& WMA,
                    FullSubMatrixHandler& WMB,
		    doublereal dCoef, 
		    const VectorHandler& XCurr,	    
		    const VectorHandler& XPrimeCurr);
   
    virtual void 
    AssStiffnessVec(SubVectorHandler& WorkVec,
                    doublereal dCoef,
		    const VectorHandler& XCurr, 
		    const VectorHandler& XPrimeCurr);

    /* Per le beam che aggiungono qualcosa alle az. interne */
    virtual void 
    AddInternalForces(Vec6& /* AzLoc */) {
        NO_OP;
    };   
   
    /*
     * Inizializza le derivate delle funzioni di forma 
     * e calcola le deformazioni geometriche iniziali
     */
    virtual void DsDxi(void);

    /*
     * Calcola la velocita' angolare delle sezioni a partire da quelle
     * dei nodi in modo consistente
     */
    virtual void Omega0(void);
   
    /* Funzione interna di restart */
    virtual ostream& Restart_(ostream& out) const;
   
  public:
    /* Costruttore normale */
    HBeam(unsigned int uL, 
	 const StructNode* pN1, const StructNode* pN2,
	 const Vec3& F1, const Vec3& F2,
	 const Mat3x3& r,
	 const ConstitutiveLaw6D* pD_I,
	 flag fOut);
   
    /* Distruttore banale */
    virtual ~HBeam(void);
   
    virtual inline void* pGet(void) const { 
        return (void*)this;
    };

    /* Tipo di trave */
    virtual Beam::Type GetBeamType(void) const {
        return Beam::ELASTIC; 
    };   
   
    /* Tipo di elemento */
    virtual Elem::Type GetElemType(void) const {
        return Elem::BEAM; 
    };
   
    /* Contributo al file di restart */
    virtual ostream& Restart(ostream& out) const;
   
    /* funzioni proprie */
   
    /* Dimensioni del workspace */
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
	*piNumRows = 12;
        *piNumCols = 12; 
    };
   
    /* Settings iniziali, prima della prima soluzione */
    void SetValue(VectorHandler& /* X */ , VectorHandler& /* XP */ ) const;
   
    /* Prepara i parametri di riferimento dopo la predizione */
    virtual void 
    AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ );
   
    /* assemblaggio residuo */
    virtual SubVectorHandler& 
    AssRes(SubVectorHandler& WorkVec,
           doublereal dCoef,
	   const VectorHandler& XCurr, 
	   const VectorHandler& XPrimeCurr);
   
    /* assemblaggio jacobiano */
    virtual VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat,
	   doublereal dCoef, 
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr);

    /* output; si assume che ogni tipo di elemento sappia, attraverso
     * l'OutputHandler, dove scrivere il proprio output */
    virtual void Output(OutputHandler& OH) const;   

    /* Funzioni proprie tipiche dei vincoli, usate durante l'assemblaggio 
     * iniziale. Le travi non sono vincoli (o meglio, non vengono considerate
     * tali), ma richiedono comunque assemblaggio iniziale. Innanzitutto 
     * le coordinate dei nodi che sono fornite vengono usate per calcolare
     * la lunghezza della trave e per la metrica della trave stessa. 
     * Inoltre durante l'assemblaggio la trave esplica la propria 
     * deformabilita' quando i vincoli tentano di muovere i nodi. */
   
    /* Numero di gradi di liberta' definiti durante l'assemblaggio iniziale
     * e' dato dai gradi di liberta' soliti piu' le loro derivate necessarie; 
     * tipicamente per un vincolo di posizione il numero di dof raddoppia, in
     * quanto vengono derivate tutte le equazioni, mentre per un vincolo di 
     * velocita' rimane inalterato. Sono possibili casi intermedi per vincoli
     * misti di posizione e velocita' */
    virtual unsigned int iGetInitialNumDof(void) const { 
        return 0;
    };
   
    /* Dimensione del workspace durante l'assemblaggio iniziale. Occorre tener
     * conto del numero di dof che l'elemento definisce in questa fase e dei
     * dof dei nodi che vengono utilizzati. Sono considerati dof indipendenti
     * la posizione e la velocita' dei nodi */
    virtual void 
    InitialWorkSpaceDim(integer* piNumRows, 
			integer* piNumCols) const {
        *piNumRows = 12;
        *piNumCols = 12; 
    };
   
    /* Setta il valore iniziale delle proprie variabili */
    virtual void SetInitialValue(VectorHandler& /* X */ ) const {
        NO_OP; 
    };
   
    /* Contributo allo jacobiano durante l'assemblaggio iniziale */
    virtual VariableSubMatrixHandler& 
    InitialAssJac(VariableSubMatrixHandler& WorkMat, 
                  const VectorHandler& XCurr);
   
    /* Contributo al residuo durante l'assemblaggio iniziale */   
    virtual SubVectorHandler& 
    InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);

    /* Accesso ai dati privati */
    virtual unsigned int iGetNumPrivData(void) const {
        return 6;
    };   
      
    virtual doublereal dGetPrivData(unsigned int i) const {
        ASSERT(i > 0 && i <= 6);
        switch (i) {
        case 1:
        case 4:
        case 5:
        case 6:
	    return DefLoc.dGet(i);
        case 2:
        case 3:
	    cerr << "Beam " << GetLabel() 
	        << ": not allowed to return shear strain" << endl;
	    THROW(ErrGeneric());
        default:
	    cerr << "Beam " << GetLabel() << ": illegal private data " 
	       << i << endl;
	    THROW(ErrGeneric());
        }
#ifndef USE_EXCEPTIONS
        return 0.;
#endif /* USE_EXCEPTIONS */
    };
   
   
    /* Accesso ai nodi */
    virtual const StructNode* pGetNode(unsigned int i) const;

    /******** PER IL SOLUTORE PARALLELO *********/        
    /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
     * utile per l'assemblaggio della matrice di connessione fra i dofs */
    virtual void 
    GetConnectedNodes(int& NumNodes, 
                      Node::Type* NdTyps, 
		      unsigned int* NdLabels) {
        NumNodes = NUMNODES;
        for (int i = 0; i <= NUMNODES-1; i++) {
            NdTyps[i] = pNode[i]->GetNodeType();
            NdLabels[i] = pNode[i]->GetLabel();
        }
    };
    /**************************************************/
   
    /* Adams output stuff */
    unsigned int iGetNumAdamsDummyParts(void) const {
        return 1;
    };
    void GetAdamsDummyPart(unsigned int part, Vec3& x, Mat3x3& R) const;
    ostream& 
    WriteAdamsDummyPartCmd(ostream& out,
                           unsigned int part, 
			   unsigned int firstId) const;
};

/* Beam - end */


class DataManager;
class MBDynParser;

extern Elem* 
ReadHBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* HBEAM_H */

