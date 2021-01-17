/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
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


#ifndef BEAM_H
#define BEAM_H

#include "myassert.h"
#include "except.h"

#include "strnode.h"
#include "elem.h"
#include "gravity.h"

#include "constltp.h"

#ifdef USE_SPARSE_AUTODIFF
#include <array>
#include "sp_gradient.h"
#include "sp_matrix_base.h"
#include "sp_matvecass.h"
#endif

extern const char* psBeamNames[];

/* ... */
class DataManager;
class MBDynParser;

/* Beam - begin */

class Beam
: virtual public Elem, public ElemGravityOwner, public InitialAssemblyElem {
    friend class AerodynamicBeam;
    friend Elem* ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);
    friend class Beam2;

  public:
    /* Tipi di travi */
    enum Type {
        UNKNOWN = -1,
	ELASTIC = 0,
	VISCOELASTIC,
	PIEZOELECTRICELASTIC,
	PIEZOELECTRICVISCOELASTIC,
	FULLYCOUPLEDPIEZOELECTRICELASTIC,

	LASTBEAMTYPE
    };

	// output mask
	enum {
		OUTPUT_NONE = 0x000U,

		OUTPUT_EP_X = (ToBeOutput::OUTPUT_PRIVATE << 0),
		OUTPUT_EP_R = (ToBeOutput::OUTPUT_PRIVATE << 1),
		OUTPUT_EP_CONFIGURATION = (OUTPUT_EP_X | OUTPUT_EP_R),

		OUTPUT_EP_F = (ToBeOutput::OUTPUT_PRIVATE << 2),
		OUTPUT_EP_M = (ToBeOutput::OUTPUT_PRIVATE << 3),
		OUTPUT_EP_FORCES = (OUTPUT_EP_F | OUTPUT_EP_M ),

		OUTPUT_EP_NU = (ToBeOutput::OUTPUT_PRIVATE << 4),
		OUTPUT_EP_K = (ToBeOutput::OUTPUT_PRIVATE << 5),
		OUTPUT_EP_STRAINS = (OUTPUT_EP_NU | OUTPUT_EP_K),

		OUTPUT_EP_NUP = (ToBeOutput::OUTPUT_PRIVATE << 6),
		OUTPUT_EP_KP = (ToBeOutput::OUTPUT_PRIVATE << 7),
		OUTPUT_EP_STRAIN_RATES = (OUTPUT_EP_NUP | OUTPUT_EP_KP),

		OUTPUT_DEFAULT = (OUTPUT_EP_F | OUTPUT_EP_M),

		OUTPUT_EP_ALL = (ToBeOutput::OUTPUT_PRIVATE_MASK)
	};

protected:
 	static const unsigned int iNumPrivData =
		+3		//  0 ( 1 ->  3) - strain
		+3		//  3 ( 4 ->  6) - curvature
		+3		//  6 ( 7 ->  9) - force
		+3		//  9 (10 -> 12) - moment
		+3		// 12 (13 -> 15) - position
		+3		// 15 (16 -> 18) - orientation vector
		+3		// 18 (19 -> 21) - angular velocity
		+3		// 21 (22 -> 24) - strain rate
		+3		// 24 (25 -> 27) - curvature rate
	;

    static unsigned int iGetPrivDataIdx_int(const char *s,
	ConstLawType::Type type);

	// output flags
	OrientationDescription od;

#ifdef USE_NETCDF
	MBDynNcVar Var_X[2],
		Var_Phi[2],
		Var_F[2],
		Var_M[2],
		Var_Nu[2],
		Var_K[2],
		Var_NuP[2],
		Var_KP[2];
#endif /* USE_NETCDF */

  public:
    class ErrGeneric : public MBDynErrBase {
	public:
		ErrGeneric(MBDYN_EXCEPT_ARGS_DECL) : MBDynErrBase(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
	};

  public:
    enum Section { S_I = 0, SII = 1, NUMSEZ = 2 };
    enum NodeName { NODE1 = 0, NODE2 = 1, NODE3 = 2, NUMNODES = 3 };
    enum Deformations { STRAIN = 0, CURVAT = 1, NUMDEFORM = 2 };

  protected:
    /* Puntatori ai nodi */
    const StructNode* pNode[NUMNODES];

    /* Offset dei nodi */
    const Vec3 f[NUMNODES];
    Vec3 fRef[NUMNODES];
    const Mat3x3 RNode[NUMNODES];

    /* Matrice di rotazione delle sezioni - non sono const perche' vengono
     * aggiornate ad ogni iterazione */
    Mat3x3 R[NUMSEZ];
    Mat3x3 RRef[NUMSEZ];
    Mat3x3 RPrev[NUMSEZ];

    /* Constitutive laws*/
    ConstitutiveLaw6DOwner* pD[NUMSEZ];

    /* Reference constitutive laws */
    Mat6x6 DRef[NUMSEZ];

    /* Per forze d'inerzia consistenti: */
    const bool bConsistentInertia;

    const doublereal dMass_I;
    const Vec3 S0_I;
    const Mat3x3 J0_I;

    const doublereal dMassII;
    const Vec3 S0II;
    const Mat3x3 J0II;

    /* Velocita' angolare delle sezioni */
    Vec3 Omega[NUMSEZ];
    Vec3 OmegaRef[NUMSEZ];

    /* Dati temporanei che vengono passati da AssRes ad AssJac */
    Vec6 Az[NUMSEZ];
    Vec6 AzRef[NUMSEZ];
    Vec6 AzLoc[NUMSEZ];
    Vec6 DefLoc[NUMSEZ];
    Vec6 DefLocRef[NUMSEZ];
    Vec6 DefLocPrev[NUMSEZ];

    // NOTE: Moved to Beam from ViscoElasticBeam for output purposes
    Vec6 DefPrimeLoc[NUMSEZ];

    Vec3 p[NUMSEZ];
    Vec3 g[NUMSEZ];
    Vec3 L0[NUMSEZ];
    Vec3 L[NUMSEZ];

    Vec3 LRef[NUMSEZ];

    doublereal dsdxi[NUMSEZ];

    /* Is first res? */
    bool bFirstRes;

    /* Funzioni di servizio */
    static Vec3
    InterpState(const Vec3& v1,
                const Vec3& v2,
		const Vec3& v3,
		enum Section Sec);
    Vec3
    InterpDeriv(const Vec3& v1,
                const Vec3& v2,
		const Vec3& v3,
		enum Section Sec);

#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     static sp_grad::SpColVector<T, 3>
     InterpState(const sp_grad::SpColVector<T, 3>& v1,
                 const sp_grad::SpColVector<T, 3>& v2,
                 const sp_grad::SpColVector<T, 3>& v3,
                 Section Sec);

     template <typename T>
     sp_grad::SpColVector<T, 3>
     InterpDeriv(const sp_grad::SpColVector<T, 3>& v1,
                 const sp_grad::SpColVector<T, 3>& v2,
                 const sp_grad::SpColVector<T, 3>& v3,
                 Section Sec);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLoc);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& AzLoc);

    virtual void
    AddInternalForces(sp_grad::SpColVector<doublereal, 6>& AzLoc, unsigned int iSez) {
        NO_OP;
    }
     
    virtual void
    AddInternalForces(sp_grad::SpColVector<sp_grad::SpGradient, 6>& AzLoc, unsigned int iSez) {
        NO_OP;
    }     
#endif
    /* Funzioni di calcolo delle matrici */
    virtual void
    AssStiffnessMat(FullSubMatrixHandler& WMA,
                    FullSubMatrixHandler& WMB,
		    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr);

    /* Per le beam che aggiungono qualcosa alle az. interne */
    virtual void
    AddInternalForces(Vec6& /* AzLoc */ , unsigned int /* iSez */ ) {
        NO_OP;
    };
     
    virtual void
    AssStiffnessVec(SubVectorHandler& WorkVec,
                    doublereal dCoef,
		    const VectorHandler& XCurr,
		    const VectorHandler& XPrimeCurr);
     
    virtual void
    AssInertiaMat(FullSubMatrixHandler& /* WMA */ ,
                  FullSubMatrixHandler& /* WMB */ ,
		  doublereal /* dCoef */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) {
        NO_OP;
    };

    virtual void
    AssInertiaVec(SubVectorHandler& /* WorkVec */ ,
                  doublereal /* dCoef */ ,
		  const VectorHandler& /* XCurr */ ,
		  const VectorHandler& /* XPrimeCurr */ ) {
        NO_OP;
    };

    /* Inizializza le derivate delle funzioni di forma
     * e calcola le deformazioni geometriche iniziali */
    virtual void DsDxi(void);

    /* Calcola la velocita' angolare delle sezioni a partire da quelle
     * dei nodi; per ora lo fa in modo brutale, ma si puo' fare anche
     * in modo consistente */
    virtual void Omega0(void);

    /* Funzione interna di restart */
    virtual std::ostream& Restart_(std::ostream& out) const;

    /* Inizializza i dati */
    void Init(void);

  public:
    /* Costruttore normale */
    Beam(unsigned int uL,
	 const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
	 const Vec3& F1, const Vec3& F2, const Vec3& F3,
	 const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3,
	 const Mat3x3& r_I, const Mat3x3& rII,
	 const ConstitutiveLaw6D* pD_I, const ConstitutiveLaw6D* pDII,
	 OrientationDescription ood,
	 flag fOut);

    /* Costruttore per la trave con forze d'inerzia consistenti */
    Beam(unsigned int uL,
	 const StructNode* pN1, const StructNode* pN2, const StructNode* pN3,
	 const Vec3& F1, const Vec3& F2, const Vec3& F3,
	 const Mat3x3& R1, const Mat3x3& R2, const Mat3x3& R3,
	 const Mat3x3& r_I, const Mat3x3& rII,
	 const ConstitutiveLaw6D* pD_I, const ConstitutiveLaw6D* pDII,
	 doublereal dM_I,
	 const Vec3& s0_I, const Mat3x3& j0_I,
	 doublereal dMII,
	 const Vec3& s0II, const Mat3x3& j0II,
	 OrientationDescription ood,
	 flag fOut);

    /* Distruttore banale */
    virtual ~Beam(void);

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

    /* Inverse Dynamics */
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP);

    /* funzioni proprie */

    /* Dimensioni del workspace; sono 36 righe perche' se genera anche le
     * forze d'inerzia consistenti deve avere accesso alle righe di definizione
     * della quantita' di moto */
    virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
        if (bConsistentInertia) {
	    *piNumRows = 36;
        } else {
	    *piNumRows = 18;
        }

        *piNumCols = 18;
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

#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    enum sp_grad::SpFunctionCall func);
#endif
     
    /* Inverse Dynamics: */
    virtual SubVectorHandler&
    AssRes(SubVectorHandler& WorkVec,
	   const VectorHandler&  XCurr ,
	   const VectorHandler&  XPrimeCurr ,
	   const VectorHandler&  XPrimePrimeCurr ,
	   InverseDynamics::Order iOrder = InverseDynamics::INVERSE_DYNAMICS);

    /* inverse dynamics capable element */
    virtual bool bInverseDynamics(void) const;

    /* assemblaggio jacobiano */
    virtual VariableSubMatrixHandler&
    AssJac(VariableSubMatrixHandler& WorkMat,
	   doublereal dCoef,
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr);

    /* assemblaggio matrici per autovalori */
    void
    AssMats(VariableSubMatrixHandler& WorkMatA,
	   VariableSubMatrixHandler& WorkMatB,
	   const VectorHandler& XCurr,
	   const VectorHandler& XPrimeCurr);

    virtual void OutputPrepare(OutputHandler &OH);

    /* output; si assume che ogni tipo di elemento sappia, attraverso
     * l'OutputHandler, dove scrivere il proprio output */
    virtual void Output(OutputHandler& OH) const;

#if 0
    /* Output di un modello NASTRAN equivalente nella configurazione corrente */
    virtual void Output_pch(std::ostream& out) const;
#endif

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
        *piNumRows = 18;
        *piNumCols = 18;
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

};

/* Beam - end */


/* ViscoElasticBeam - begin */

class ViscoElasticBeam : virtual public Elem, public Beam {
  protected:

    /* Derivate di deformazioni e curvature */
    Vec3 LPrime[NUMSEZ];
    Vec3 gPrime[NUMSEZ];

    Vec3 LPrimeRef[NUMSEZ];

    // NOTE: Moved to Beam from ViscoElasticBeam for output purposes
    // Vec6 DefPrimeLoc[NUMSEZ];

    Vec6 DefPrimeLocRef[NUMSEZ];

    Mat6x6 ERef[NUMSEZ];

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

    /* Inizializza i dati */
    void Init(void);

#ifdef USE_SPARSE_AUTODIFF
     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<doublereal, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& gPrime,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& Omega,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<doublereal, 3>, NUMSEZ>& LPrime,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& DefPrimeLoc,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<doublereal, 6>, NUMSEZ>& AzLoc);

     inline void
     UpdateState(const std::array<sp_grad::SpMatrixA<sp_grad::SpGradient, 3, 3>, NUMSEZ>& R,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& p,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& g,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& gPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& Omega,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& L,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 3>, NUMSEZ>& LPrime,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& DefPrimeLoc,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& Az,
                 const std::array<sp_grad::SpColVectorA<sp_grad::SpGradient, 6>, NUMSEZ>& AzLoc);
#endif
  public:
    /* Costruttore normale */
    ViscoElasticBeam(unsigned int uL,
	             const StructNode* pN1,
		     const StructNode* pN2,
		     const StructNode* pN3,
	             const Vec3& F1,
		     const Vec3& F2,
		     const Vec3& F3,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     const Mat3x3& R3,
	             const Mat3x3& r_I,
		     const Mat3x3& rII,
	             const ConstitutiveLaw6D* pD_I,
		     const ConstitutiveLaw6D* pDII,
		     OrientationDescription ood,
		     flag fOut);

    /* Costruttore per la trave con forze d'inerzia consistenti */
    ViscoElasticBeam(unsigned int uL,
                     const StructNode* pN1,
		     const StructNode* pN2,
		     const StructNode* pN3,
		     const Vec3& F1,
		     const Vec3& F2,
		     const Vec3& F3,
		     const Mat3x3& R1,
		     const Mat3x3& R2,
		     const Mat3x3& R3,
		     const Mat3x3& r_I,
		     const Mat3x3& rII,
		     const ConstitutiveLaw6D* pD_I,
		     const ConstitutiveLaw6D* pDII,
		     doublereal dM_I,
		     const Vec3& s0_I,
		     const Mat3x3& j0_I,
		     doublereal dMII,
		     const Vec3& s0II,
		     const Mat3x3& j0II,
		     OrientationDescription ood,
		     flag fOut);

    /* Distruttore banale */
    virtual ~ViscoElasticBeam(void) {
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

    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

    /* Inverse Dynamics */
    virtual void
    AfterConvergence(const VectorHandler& X, const VectorHandler& XP,
    		const VectorHandler& XPP);

#ifdef USE_SPARSE_AUTODIFF
     template <typename T>
     inline void
     AssRes(sp_grad::SpGradientAssVec<T>& WorkVec,
	    doublereal dCoef,
	    const sp_grad::SpGradientVectorHandler<T>& XCurr,
	    const sp_grad::SpGradientVectorHandler<T>& XPrimeCurr,
	    enum sp_grad::SpFunctionCall func);

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
#endif
     
    virtual doublereal dGetPrivData(unsigned int i) const;
};

/* ViscoElasticBeam - end */

extern void
ReadBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od);

extern void
ReadOptionalBeamCustomOutput(DataManager* pDM, MBDynParser& HP, unsigned int uLabel,
	Beam::Type BT, unsigned& uFlags, OrientationDescription& od);

extern Elem*
ReadBeam(DataManager* pDM, MBDynParser& HP, unsigned int uLabel);

#endif /* BEAM_H */
