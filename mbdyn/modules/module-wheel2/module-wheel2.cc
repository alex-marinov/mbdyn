#include <mbconfig.h>

#include <loadable.h>

/*
 *	2 corpi: 
 *	-	Wheel
 *	-	Ground
 *
 *	Note: 
 *	-	Axle e Wheel sono collegati da un giunto che consente
 *		solo rotazione relativa attorno ad un asse (axle)
 *	-	Si assume che il centro della ruota coincida con
 *		la posizione del corpo ruota
 *	-	Ground supporta un piano definito dalla normale e da
 *		un punto (potra' essere reso deformabile e "non piano"
 *		in futuro)
 *	-	Le forze sono applicate nel "punto di contatto", che viene
 *		calcolato in base a considerazioni geometriche sulla
 *		posizione ed orientazione relativa tra Wheel e Ground.
 */

struct module_wheel {
	/*
	 * Ruota
	 */
	StructNode *pWheel;

	/*
	 * Direzione asse ruota
	 */
	Vec3 WheelAxle;
	
	/*
	 * Terreno
	 */
	StructNode *pGround;

	/*
	 * Posizione e orientazione del terreno
	 */
	Vec3 GroundPosition;
	Vec3 GroundDirection;

	/*
	 * Geometria ruota
	 */
	doublereal dRadius;
	doublereal dInternalRadius;
	doublereal dVolCoef;

	doublereal dRefArea;
	doublereal dRNP;		/* R+nG'*pG */
	doublereal dV0;

	/*
	 * Proprieta' pneumatico
	 */
	doublereal dP0;
	doublereal dGamma;
	doublereal dHystVRef;

	/*
	 * Output
	 */
	Vec3 F;
	Vec3 M;
	doublereal dInstRadius;
	doublereal dDeltaL;
	doublereal dVn;
	doublereal dSr;
};

/* funzioni di default */
static void*
read(LoadableElem* pEl,
		DataManager* pDM,
		MBDynParser& HP,
		const DriveHandler* pDH)
{
	DEBUGCOUTFNAME("read");
	
	/* allocation of user-defined struct */
	module_wheel* p = NULL;
	SAFENEW(p, module_wheel, EMmm);

	/*
	 * leggo la ruota
	 */
	p->pWheel = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);

	/*
	 * leggo l'orientazione dell'asse ruota nel sistema locale
	 */
	ReferenceFrame RF = ReferenceFrame(p->pWheel);
	p->WheelAxle = HP.GetVecRel(RF);
	
	/*
	 * leggo il terreno
	 */
	p->pGround = (StructNode *)pDM->ReadNode(HP, Node::STRUCTURAL);
	
	/*
	 * leggo posizione ed orientazione del terreno nel sistema del nodo
	 */
	RF = ReferenceFrame(p->pGround);
	p->GroundPosition = HP.GetPosRel(RF);
	p->GroundDirection = HP.GetVecRel(RF);

	/*
	 * normalizzo l'orientazione del terreno
	 */
	doublereal d = p->GroundDirection.Dot();
	if (d <= DBL_EPSILON) {
		cerr << "null direction at line " << HP.GetLineData() << endl;
		THROW(DataManager::ErrGeneric());
	}
	p->GroundDirection /= sqrt(d);

	/*
	 * Geometria ruota
	 */
	p->dRadius = HP.GetReal();
	p->dInternalRadius = HP.GetReal();
	p->dVolCoef = HP.GetReal();
	
	/*
	 * Area di riferimento
	 */
	p->dRefArea = 3.7*p->dInternalRadius*sqrt(
			p->dInternalRadius*(2.*p->dRadius-p->dInternalRadius)
			);

	/*
	 * termine per il calcolo di Delta L
	 */
	p->dRNP = p->dRadius+p->GroundPosition*p->GroundDirection;

	/*
	 * Volume di riferimento
	 */
	p->dV0 = 2.*M_PI*(p->dRadius-p->dInternalRadius)
		*M_PI*p->dInternalRadius*p->dInternalRadius*p->dVolCoef;

	/*
	 * Proprieta' pneumatico
	 */
	p->dP0 = HP.GetReal();
	p->dGamma = HP.GetReal();
	p->dHystVRef = HP.GetReal();
	
	return (void *)p;
}

unsigned int
i_get_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_dof");
	return 0;
}

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
	DEBUGCOUTFNAME("output");
	
	if (pEl->fToBeOutput()) {
		module_wheel* p = (module_wheel *)pEl->pGetData();      
		ostream& out = OH.Loadable();

		out << setw(8) << pEl->GetLabel() << " ",
			p->F.Write(out, " ") << " ",
			p->M.Write(out, " ") << " "
			<< p->dInstRadius << " "
			<< p->dDeltaL << " "
			<< p->dVn << " " 
			<< p->dSr << endl;
	}
}

ostream&
restart(const LoadableElem* pEl, ostream& out)
{
	DEBUGCOUTFNAME("restart");
	return out << "not implemented yet;" << endl;
}

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
	DEBUGCOUTFNAME("work_space_dim");
	*piNumRows = 12;
	*piNumCols = 12;
}

static VariableSubMatrixHandler& 
ass_jac(LoadableElem* pEl, 
	VariableSubMatrixHandler& WorkMat,
	doublereal dCoef, 
	const VectorHandler& XCurr,
	const VectorHandler& XPrimeCurr)
{  
	DEBUGCOUTFNAME("ass_jac");   
	WorkMat.SetNullMatrix();
	
	return WorkMat;
}

static SubVectorHandler& 
ass_res(LoadableElem* pEl, 
	SubVectorHandler& WorkVec,
	doublereal dCoef,
	const VectorHandler& XCurr, 
	const VectorHandler& XPrimeCurr)
{
	DEBUGCOUTFNAME("ass_res");
	
	module_wheel* p = (module_wheel *)pEl->pGetData();

	/*
	 * Orientazione del terreno nel sistema assoluto
	 */
	Vec3 n = p->pGround->GetRCurr()*p->GroundDirection;
	
	/*
	 * Distanza Wheel Ground nel sistema assoluto
	 */
	Vec3 x = p->pWheel->GetXCurr()-p->pGround->GetXCurr();

	/*
	 * se dDeltaL > 0 c'e' contatto, altrimenti no
	 */
	doublereal dDeltaL = p->dRNP - x*n;

	/*
	 * Reset dati per output
	 */
	p->dDeltaL = dDeltaL;
	p->dInstRadius = p->dRadius-dDeltaL;
	p->F = Zero3;
	p->M = Zero3;
	
	if (dDeltaL < 0.) {
		/*
		 * Non assemblo neppure il vettore ;)
		 */
		WorkVec.Resize(0);

		return WorkVec;
	}
	
	/*
	 * Dimensiono il vettore
	 */
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
	WorkVec.Resize(iNumRows);
	WorkVec.Reset(0.);

	integer iGroundFirstMomIndex = p->pGround->iGetFirstMomentumIndex();
	integer iWheelFirstMomIndex = p->pWheel->iGetFirstMomentumIndex();

	/* Indici equazioni */
	for (int iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.fPutRowIndex(iCnt, iGroundFirstMomIndex+iCnt);
		WorkVec.fPutRowIndex(6+iCnt, iWheelFirstMomIndex+iCnt);
	}

	/*
	 * Velocita' tra Wheel (nell'asse)
	 * e Ground nel sistema assoluto
	 */
	Vec3 va = p->pWheel->GetVCurr()
		-p->pGround->GetVCurr()-(p->pGround->GetWCurr()).Cross(
			p->pGround->GetRCurr()*p->GroundPosition
			);
	/*
	 * Velocita' tra Wheel (nel punto di contatto) 
	 * e Ground nel sistema assoluto
	 */
	Vec3 v = va-(p->pWheel->GetWCurr()).Cross(n*p->dInstRadius);
	
	/*
	 * Componente normale al terreno della velocita'
	 * (positiva se la ruota si allontana dal terreno)
	 */
	p->dVn = n*v;

	/*
	 * Direzione di "avanzamento": asse ruota cross normale al ground
	 */
	Vec3 fwd = (p->pWheel->GetRCurr()*p->WheelAxle).Cross(n);
	doublereal d = fwd.Dot();
	if (d < DBL_EPSILON) {
		cerr << "wheel axle is orthogonal to the ground" << endl;
		THROW(DataManager::ErrGeneric());
	}
	fwd /= sqrt(d);

	/*
	 * Slip ratio
	 */
	p->dSr = 0.;
	doublereal fwdva = fwd.Dot(va);
	if (fabs(fwdva) > DBL_EPSILON) {
		p->dSr = fwd.Dot(v)/fwdva;
		if (fabs(p->dSr) > 1.) {
			p->dSr = copysign(1., p->dSr);
		}
	}
	
	/*
	 * Stima dell'area di contatto
	 */
	doublereal dA = p->dRefArea*(dDeltaL/p->dRadius);
	
	/*
	 * Stima del volume di compenetrazione tra ruota e terreno
	 */
	doublereal dDeltaV = .5*dA*dDeltaL;

	/*
	 * Stima della pressione del pneumatico (adiabatica)
	 */
	doublereal dP = p->dP0*pow(p->dV0/(p->dV0-dDeltaV), p->dGamma);

	/*
	 * Il punto di applicazione della forza e' xW - R * nG ;
	 * per la forza normale si puo' usare anche la posizione
	 * della ruota, nell'ipotesi che il punto di contatto
	 * stia nell'intersezione della normale al ground passante
	 * per l'asse ruota
	 */
	Vec3 pc = p->pWheel->GetXCurr()-(n*p->dRadius);

	/*
	 * Forza
	 */
	p->F = n*(dA*dP*(1.+tanh(-p->dVn/p->dHystVRef)));
	p->M = (pc-p->pWheel->GetXCurr()).Cross(p->F);
	
	WorkVec.Sub(1, p->F);
	WorkVec.Sub(4, (pc-p->pGround->GetXCurr()).Cross(p->F));
	WorkVec.Add(7, p->F);
	WorkVec.Add(10, p->M);

	return WorkVec;
}

void
before_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP,
		VectorHandler& XPrev,
		VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

void
after_predict(const LoadableElem* pEl, 
		VectorHandler& X,
		VectorHandler& XP)
{
	DEBUGCOUTFNAME("after_predict");
}

void
update(LoadableElem* pEl, 
		const VectorHandler& X,
		const VectorHandler& XP)
{
	DEBUGCOUTFNAME("update");
}

unsigned int
i_get_initial_num_dof(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_initial_num_dof");
	return 0;
}

void
initial_work_space_dim(const LoadableElem* pEl, 
		integer* piNumRows, 
		integer* piNumCols)
{
	DEBUGCOUTFNAME("initial_work_space_dim");
	*piNumRows = 0;
	*piNumCols = 0;   
}

VariableSubMatrixHandler& 
initial_ass_jac(LoadableElem* pEl, 
		VariableSubMatrixHandler& WorkMat, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_jac");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->InitialWorkSpaceDim(&iNumRows, &iNumCols);
   
	FullSubMatrixHandler& WM = WorkMat.SetFull();
	WM.ResizeInit(iNumRows, iNumCols, 0.);
	
	module_wheel* p = (module_wheel *)pEl->pGetData();
   
	/* set sub-matrix indices and coefs */

	return WorkMat;
}

SubVectorHandler& 
initial_ass_res(LoadableElem* pEl, 
		SubVectorHandler& WorkVec, 
		const VectorHandler& XCurr)
{
	DEBUGCOUTFNAME("initial_ass_res");
	integer iNumRows = 0;
	integer iNumCols = 0;
	pEl->WorkSpaceDim(&iNumRows, &iNumCols);
	
	WorkVec.Resize(iNumRows);
	
	module_wheel* p = (module_wheel *)pEl->pGetData(); 
	
	/* set sub-vector indices and coefs */
   
	return WorkVec;
}

void
set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
	DEBUGCOUTFNAME("set_value");
}
   
void
set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
	DEBUGCOUTFNAME("set_initial_value");
}

unsigned int
i_get_num_priv_data(const LoadableElem* pEl)
{
	DEBUGCOUTFNAME("i_get_num_priv_data");
	return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
	DEBUGCOUTFNAME("d_get_priv_data");
	ASSERT(pEl->iGetNumPrivData() > 0);
	if (i > pEl->iGetNumPrivData()) {
		cerr << "Module-template Elem: illegal private data index "
			<< i << endl;      
		THROW(ErrGeneric());
	}
	
	/* return i-th priv data */

	return 0.;
}

static void
destroy(LoadableElem* pEl)
{
	DEBUGCOUTFNAME("destroy");
	module_wheel* p = (module_wheel *)pEl->pGetData();
	
	SAFEDELETE(p, EMmm);
}

static struct
LoadableCalls lc = {
	read, /* */
	NULL /* i_get_num_dof */ ,
	NULL /* set_dof */ ,
	output, /* */
	NULL /* restart */ ,
	work_space_dim, /* */
	ass_jac, /* */
	NULL /* ass_eig */ ,
	ass_res, /* */
	NULL /* before_predict */ ,
	NULL /* after_predict */ ,
	NULL /* update */ ,
	NULL /* i_get_initial_num_dof */ ,
	NULL /* initial_work_space_dim */ ,
	NULL /* initial_ass_jac */ ,
	NULL /* initial_ass_res */ ,
	NULL /* set_value */ ,
	NULL /* set_initial_value */ ,
	NULL /* i_get_num_priv_data */ ,
	NULL /* d_get_priv_data */ ,
	destroy
};

extern "C" {
void *calls = &lc;
}

