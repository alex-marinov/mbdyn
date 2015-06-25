/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
 * the Free Software Foundation (version 2 of the License).
 * 
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

#include <matvecexp.h>

#include "strnode.h"
#include "elem.h"
#include "gravity.h"

#include "constltp.h"

// #define VISCOELASTIC_BEAM /* uncomment when ViscoElasticBeam is available */
#undef VISCOELASTIC_BEAM
// #define PIEZO_BEAM /* uncomment when PiezoBeam is available */
#undef PIEZO_BEAM

/* Beam - begin */

class HBeam
: virtual public Elem, public ElemGravityOwner, public InitialAssemblyElem {
    friend class AerodynamicBeam;

  public:
    class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};
   
  private:
    Beam::Type BeamT;
   
  protected:   
    enum NodeName { NODE1 = 0, NODE2 = 1, NUMNODES = 2 };
    enum Deformations { STRAIN = 0, CURVAT = 1, NUMDEFORM = 2 };
   
    /* Puntatori ai nodi */
    const StructNode* pNode[NUMNODES];
   
    /* Offset dei nodi */
    Vec3 f[NUMNODES];
    Vec3 fRef[NUMNODES];

    /* Matarici di rotazione iniziali dei nodi */
    Mat3x3 Rn[NUMNODES];
   
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
    Vec6 DefLoc;
    Vec6 DefLocRef;

    Vec3 p;
    Vec3 g;
    
    Vec3 L0;
    Vec3 L;
    Vec3 LRef;   
   
    doublereal xi;
    doublereal dsdxi;
    doublereal dxids;
   
    /* Is first res? */
    bool bFirstRes;

    /* Temporanei per interpolazione elicoidale */
    /* posiz. punto medio:		p */
    /* orientazione punto medio:	R */
    /* F (dp/ds)			L */
    /* assiale di dR/ds*R^T		Rho */
    Vec3 Rho;
    Vec3 Rho0;
    /* derivate di p e L rispetto a p e PhiDelta,
     * di R e Rho rispetto a PhiDelta */
    Mat3x3 pdp[NUMNODES];
    Mat3x3 pdP[NUMNODES];
    Mat3x3 Ldp[NUMNODES];
    Mat3x3 LdP[NUMNODES];
    Mat3x3 RdP[NUMNODES];
    Mat3x3 RhodP[NUMNODES];

    /*
     * Nota: Rho0 e L0 non sono derivate "iniziali" nel sistema iniziale,
     * ma nel sistema del materiale: "L0" = R0^T L0, "Rho0" = R0^T Rho0
     */

    /* pesi e derivate */

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
   
    /* Inizializza le derivate delle funzioni di forma 
     * e calcola le deformazioni geometriche iniziali */
    virtual void DsDxi(void);

    /* Calcola la velocita' angolare delle sezioni a partire da quelle
     * dei nodi; per ora lo fa in modo brutale, ma si puo' fare anche 
     * in modo consistente */
    virtual void Omega0(void);
  
    /* cambia il tipo (?) */
    void SetBeamType(Beam::Type T) { 
        BeamT = T;
    };

    /* Funzione interna di restart */
    virtual std::ostream& Restart_(std::ostream& out) const;
   
  public:
    /* Costruttore normale */
    HBeam(unsigned int uL, 
	 const StructNode* pN1, const StructNode* pN2,
	 const Vec3& F1, const Vec3& F2,
	 const Mat3x3& R1, const Mat3x3& R2,
	 const ConstitutiveLaw6D* pd,
	 flag fOut);

    /* Distruttore banale */
    virtual ~HBeam(void);
   
    /* Tipo di trave */
    virtual Beam::Type GetBeamType(void) const {
        return Beam::ELASTIC; 
    };   
   
    /* Tipo di elemento */
    virtual Elem::Type GetElemType(void) const {
        return Elem::BEAM; 
    };
   
    /* Contributo al file di restart */
    virtual std::ostream& Restart(std::ostream& out) const;
   
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

    /* funzioni proprie */
   
    /* Dimensioni del workspace; sono 36 righe perche' se genera anche le 
     * forze d'inerzia consistenti deve avere accesso alle righe di definizione
     * della quantita' di moto */
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
	*piNumRows = 12;
        *piNumCols = 12; 
    };
   
    /* Settings iniziali, prima della prima soluzione */
    void SetValue(DataManager *pDM,
		    VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		    SimulationEntity::Hints *ph = 0);
   
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
    virtual void SetInitialValue(VectorHandler& /* X */ ) {
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
    virtual unsigned int iGetNumPrivData(void) const;
    virtual unsigned int iGetPrivDataIdx(const char *s) const;
    virtual doublereal dGetPrivData(unsigned int i) const;
   
   
    /* Accesso ai nodi */
    virtual const StructNode* pGetNode(unsigned int i) const;

    /******** PER IL SOLUTORE PARALLELO *********/        
    /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
     * utile per l'assemblaggio della matrice di connessione fra i dofs */
    virtual void 
    GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
        connectedNodes.resize(NUMNODES);
        for (int i = 0; i < NUMNODES; i++) {
            connectedNodes[i] = pNode[i];
        }
    };
    /**************************************************/

    /* Adams output stuff */
    unsigned int iGetNumDummyParts(void) const {
        return 1;
    };
    void GetDummyPartPos(unsigned int part, Vec3& x, Mat3x3& R) const;
    void GetDummyPartVel(unsigned int part, Vec3& v, Vec3& w) const;
#ifdef USE_ADAMS 
    std::ostream& 
    WriteAdamsDummyPartCmd(std::ostream& out,
                           unsigned int part, 
			   unsigned int firstId) const;
#endif /* USE_ADAMS */
};

/* Beam - end */


#ifdef VISCOELASTIC_BEAM
/* ViscoElasticBeam - begin */

class ViscoElasticHBeam : virtual public Elem, public HBeam {
  protected:
   
    /* Derivate di deformazioni e curvature */
    Vec3 LPrime; 
    Vec3 gPrime;  
      
    Vec3 LPrimeRef;
   
    Vec6 DefPrimeLoc;
    Vec6 DefPrimeLocRef;

    Mat6x6 ERef;
      
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
   
  public:
    /* Costruttore normale */
    ViscoElasticHBeam(unsigned int uL, 
	             const StructNode* pN1, 
		     const StructNode* pN2, 
	             const Vec3& F1, 
		     const Vec3& F2, 
	             const Mat3x3& r, 
	             const ConstitutiveLaw6D* pd, 
		     flag fOut);
   
    /* Distruttore banale */
    virtual ~ViscoElasticHBeam(void) { 
        NO_OP;
    };

    /* Tipo di trave */
    virtual Beam::Type GetBeamType(void) const {
        return Beam::VISCOELASTIC; 
    };

    /* Settings iniziali, prima della prima soluzione */
    void SetValue(DataManager *pDM,
		    VectorHandler& /* X */ , VectorHandler& /* XP */ ,
		    SimulationEntity::Hints *ph = 0);
   
    /* Prepara i parametri di riferimento dopo la predizione */
    virtual void 
    AfterPredict(VectorHandler& /* X */ , VectorHandler& /* XP */ );   
};

/* ViscoElasticBeam - end */
#endif /* VISCOELASTIC_BEAM */

class DataManager;
class MBDynParser;

extern Elem* 
ReadHBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* BEAM2_H */

