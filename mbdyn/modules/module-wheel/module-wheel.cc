#include <mbconfig.h>

#include <loadable.h>
#include <myf2c.h>
#include <wheel.h>

const int iNumOutputs = 7;

struct module_wheel {
   
   StructNode* pWheel;
   StructNode* pStrut;
   
   Vec3 WheelOffset;
   Mat3x3 WheelR;
   
   /*
    C DESCRIZIONE AREE DI LAVORO LOCALI
    C
    C  CC(9) DEFINIZIONE GEOMETRICA RUOTA
    C 1-3    COSENI DIRETTORI ASSE RUOTA
    C 4-6    POSIZIONE ORIGINE ASSE RUOTA
    C 7      RAGGIO DELLA RUOTA
    C 8-9    SEMISPESSORI DELLA RUOTA NELLA DIREZIONE ASSE
    C
    C CDP(12) DEFINIZIONE DEL TERRENO
    C 1-9     COSENI DIRETTORI PIANO TERRENO
    C 10-12   PUNTO DI RIFERIMENTO TERRENO
    C
    C  PP PROPRIETA' DEL PIANO DI CONTATTO
    C     1 PRESSIONE RIFERIMENTO
    C     2 VOLUME DI RIFERIMENTO
    C     3 VELOCITA' DI RIFERIMENTO ARATURA
    C     4 COEFFICIENTE RIDUTTIVO ATTRITO
    C     5 COEFFICIENTE RIDUTTIVO ARATURA
    C     6 COEFFICIENTE RIDUTTIVO DERIVA
    C
    C  WP PROPRIETA' PNEUMATICO
    C     1 PRESSIONE DI RIFERIMENTO       press
    C     2 VOLUME DI RIFERIMENTO          calcolato=2.*PIGRE^2*RC^2*(RR-RC)*Cvol
    C     3 VELOCITA' DI RIFERIMENTO PER SMORZAMENTO  wrv
    C     4 COEFFICIENTE ATTRITO MASSIMO (METTERE A ZERO)
    C     5 AREA DI RIFERIMENTO PER LA DERIVA   arrif
    C     6 COEFFICIENTE DERIVA                 rmust
    C     7 COEFFICIENTE ARRETRAMENTO           xch
    C     8 COEFFICIENTE RIDUTTIVO DELLA PRESSIONE PER CAMBER  cracmb
    C     9 COEFFICIENTE STRISCIAMENTO DI ATTRITO MASSIMO      smumax
    C    10 COEFFICIENTE DI ATTRITO MASSIMO                    rmum
    C    11 COEFFICIENTE DI ATTRITO MINIMO PER MASSIMO STRISCIAMENTO rmus
    C
    C  GP PROPRIETA' GLOBALI
    C     1 COEFFICIENTE DELLA POLITROPICA PER PRESSIONE PNEUMATICO   espr
    C     2 VELOCITA' LIMITE PER IL CALCOLO DELLE FORZE               trv
    C     3 LIBERO PER FUTURI USI   
    */
   
   // user-defined struct
   doublereal cc[9];
   doublereal cdp[12];
   doublereal wp[12];     /* abbiamo aggiunto ARRIF alla fine */
   
   doublereal pp[6];
   doublereal gp[3];      /* nota: gp[3] era libero per futuri usi;
			   * abbiamo aggiunto RK */

   doublereal vw[6];      /* velocita', input */
   doublereal res[iNumOutputs];     /* forze e mom., output */
   doublereal rm;         /* momento asse ruota, output */
   doublereal ome;        /* vel. angolare relativa ruota */
   doublereal omep;       /* acc. angolare ruota, input (non supportato) */
   
   integer ircw;          /* flag di contatto, output */
   integer ipc;           /* tipo di simulazione, input */
   
   integer npft;          /* numero punti tabella coeff. attrito */
   doublereal* tabfat;    /* tabella coeff. attrito */
   
   doublereal beta;       /* geometria eventuale braccio oscillante (non supportato) */
   doublereal betap;
   doublereal bleng;
};

// #ifdef DEBUG
static void
out_array(ostream& out, const char *name, doublereal *pd, integer size)
{
   out << name << "=[";
   for (int i = 0; i < size-1; i++) {
      out << pd[i] << ",";
   }
   out << pd[size-1] << "]" << endl;
}
// #endif

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
   
   /* initialisation */
   for (int i = 0; i < 9; i++) {
      p->cc[i] = 0.;
   }
   for (int i = 0; i < 12; i++) {
      p->cdp[i] = 0.;
   }
   for (int i = 0; i < 12; i++) {
      p->wp[i] = 0.;
   }
   for (int i = 0; i < 6; i++) {
      p->pp[i] = 0.;
   }
   for (int i = 0; i < 3; i++) {
      p->gp[i] = 0.;
   }
   for (int i = 0; i < 6; i++) {
      p->vw[i] = 0.;
   }
   for (int i = 0; i < iNumOutputs; i++) {
      p->res[i] = 0.;
   }
   
   p->rm = 0.;
   p->ome = 0.;
   p->omep = 0.; /* non supportato */
   
   p->ircw = 0;
   p->ipc = 0;
   
   p->npft = 0;
   p->tabfat = NULL;
   
   p->beta = 0.; /* non supportati */
   p->betap = 0.;
   p->bleng = 0.;
    
   
   /* nodo ruota */
   unsigned int uNode = (unsigned int)HP.GetInt();
   DEBUGCOUT("Wheel Node " << uNode << endl);
       
   /* verifica di esistenza del nodo */  
   if ((p->pWheel = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << "line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;    
      THROW(DataManager::ErrGeneric());
   }                  
   
   ReferenceFrame RF(p->pWheel);
   p->WheelOffset = HP.GetPosRel(RF);
   p->WheelR = HP.GetRotRel(RF);
   
   /* nodo strut */
   uNode = (unsigned int)HP.GetInt();
   DEBUGCOUT("Strut Node " << uNode << endl);
       
   /* verifica di esistenza del nodo */  
   if ((p->pStrut = pDM->pFindStructNode(uNode)) == NULL) {
      cerr << "line " << HP.GetLineData() 
	<< ": structural node " << uNode
	<< " not defined" << endl;    
      THROW(DataManager::ErrGeneric());
   }                  

   /* legge dati ruota */
   const char* filename = HP.GetFileName();
   
   /* apre il file di ingresso in formato graal */
   ifstream in(filename);
   
   /* salta le prime 3 righe */
   char* s = NULL;
   
   /* legge alcuni parametri */
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s 
	      << "\")" << endl);
   /* legge alcuni parametri */
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s
	      << "\")" << endl);
   /* legge alcuni parametri */
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s 
	      << "\")" << endl);
   
   /* raggio ruota */
   in >> p->cc[6];
   
   /* raggio toro */
   doublereal dRaggioToro;
   in >> dRaggioToro;
   p->cc[7] = -dRaggioToro;
   p->cc[8] = dRaggioToro;
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading R (=" << p->cc[6]
	      << ") RT (=" << dRaggioToro
	      << ", remark: \"" << s << "\")" << endl);

   in >> p->wp[0];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading P (=" << p->wp[0]
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->gp[0];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading GP(1) (polytropic exp.=" << p->gp[0]
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->gp[2];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading RK (=" << p->gp[2]
	      << ", remark: \"" << s << "\")" << endl);
   
#warning "VERIFICARE!"   
   doublereal csi;
   in >> csi;
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading CSI (=" << csi
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->wp[5];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading RMUST (=" << p->wp[5]
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->wp[6];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading XCH (=" << p->wp[6]
	      << ", remark: \"" << s << "\")" << endl);
   
   doublereal coevol;
   in >> coevol;
   p->wp[1] = 2.*M_PI*M_PI*p->cc[8]*p->cc[8]*(p->cc[6]-p->cc[8])*coevol;
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading CVOL (=" << coevol
	      << ", VRef=" << p->wp[1]
	      << ", remark: \"" << s << "\")" << endl);
 
   in >> p->wp[2];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading Vrif1 (=" << p->wp[2]
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->gp[1];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading Vrif2 (=" << p->gp[1]
	      << ", remark: \"" << s << "\")" << endl);
   
   in >> p->wp[7];
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading CRACMB (=" << p->wp[7]
	      << ", remark: \"" << s << "\")" << endl);
   
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s
	      << "\")" << endl);   
   
   /* tabella coeff. di attrito */
   in >> p->npft;
   in.gets(&s);
   DEBUGLCOUT(MYDEBUG_INPUT, "Reading NPFT (=" << p->npft
	      << ", remark: \"" << s << "\")" << endl);
   if (p->npft <= 0) {
      cerr << "Loadable::read(" << pEl->GetLabel() << "): illegal NPFT " 
	<< p->npft << endl;
      THROW(ErrGeneric());
   }

   SAFENEWARR(p->tabfat, doublereal, 2*p->npft, EMmm);
   
   for (int i = 0; i < p->npft; i++) {
      in >> p->tabfat[i] >> p->tabfat[p->npft+i];
      in.gets(&s);
      DEBUGLCOUT(MYDEBUG_INPUT, "Reading s(" << i+1 << ")=" << p->tabfat[i]
		 << ", f(" << i+1 << ")=" << p->tabfat[p->npft+i]
		 << " (remark: \"" << s << "\")" << endl);
   }
   
#warning "VERIFICARE (dati terreno, da rendere inseribili)"   
   p->pp[0] = 1.e10;
   p->pp[1] = .3;
   p->pp[2] = .5;
   p->pp[3] = 1.;
   p->pp[4] = .8;
   p->pp[5] = 1.;
   
   return (void *)p;
}

#if 0
unsigned int i_get_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_dof");
   return 0;
}
#endif

static void
output(const LoadableElem* pEl, OutputHandler& OH)
{
   DEBUGCOUTFNAME("output");
   
   if (pEl->fToBeOutput()) {
      module_wheel* p = (module_wheel *)pEl->pGetData();      
      ostream& out = OH.Loadable();
   
      out << setw(8) << pEl->GetLabel();
      for (int i = 0; i < iNumOutputs; i++) {
	 out << " " << p->res[i];
      }
      out << endl;
   }
}

#if 0
ostream& restart(const LoadableElem* pEl, ostream& out)
{
   DEBUGCOUTFNAME("restart");
   return out << "not implemented yet;" << endl;
}
#endif

static void
work_space_dim(const LoadableElem* pEl, 
		    integer* piNumRows, 
		    integer* piNumCols)
{
   DEBUGCOUTFNAME("work_space_dim");
   *piNumRows = 6;
   *piNumCols = 1;
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
   integer iNumRows = 0;
   integer iNumCols = 0;
   pEl->WorkSpaceDim(&iNumRows, &iNumCols);
   
   WorkVec.Resize(iNumRows);
   WorkVec.Reset(0.);
   
   module_wheel* p = (module_wheel *)pEl->pGetData();

   /* set sub-vector indices and coefs */
   integer iIndex = p->pWheel->iGetFirstMomentumIndex();
   
   Mat3x3 RStrut(p->pStrut->GetRCurr());      /* mat. R wheel */
   Mat3x3 RWheel(p->pWheel->GetRCurr());      /* mat. R wheel */
   Mat3x3 R(RWheel*p->WheelR);                /* mat. R asse ruota (asse 3 locale) */
   Vec3 WheelOffset(RStrut*p->WheelOffset);   /* offset wheel sist. globale */
   Vec3 X(p->pStrut->GetXCurr()+WheelOffset); /* pos. asse ruota sist. glob. */
   Vec3 OmegaStrut(p->pStrut->GetWCurr());    /* omega strut sist. globale */
   Vec3 OmegaWheel(p->pWheel->GetWCurr());    /* omega ruota sist. globale */
   Vec3 V(p->pStrut->GetVCurr()+OmegaStrut.Cross(WheelOffset)); /* vel. asse ruota sist. glob. */
   Vec3 R3(R.GetVec(3));                      /* direz. asse ruota */
   
   p->ome = R3.Dot(OmegaStrut-OmegaWheel);    /* omega relativa */
   
   p->cc[0] = R3.dGet(1);
   p->cc[1] = R3.dGet(2);
   p->cc[2] = R3.dGet(3);
   
   p->cc[3] = X.dGet(1);
   p->cc[4] = X.dGet(2);
   p->cc[5] = X.dGet(3);

   p->vw[0] = V.dGet(1);
   p->vw[1] = V.dGet(2);
   p->vw[2] = V.dGet(3);
 
   p->vw[3] = OmegaStrut.dGet(1);
   p->vw[4] = OmegaStrut.dGet(2);
   p->vw[5] = OmegaStrut.dGet(3);

//#ifdef DEBUG
   cout << "in:" << endl;
   out_array(cout, "\tvw", p->vw, 6);
   out_array(cout, "\tome", &(p->ome), 1);
   out_array(cout, "\tcc", p->cc, 9);
   out_array(cout, "\trm", &p->rm, 1);
   cout << "\tircw=" << p->ircw << endl;
   cout << "\tipc=" << p->ipc << endl;
   out_array(cout, "\tomep", &p->omep, 1);
   out_array(cout, "\ttabfat", p->tabfat, 2*p->npft);
   cout << "\tnpft=" << p->npft << endl;
   out_array(cout, "\tbeta", &p->beta, 1);
   out_array(cout, "\tbetap", &p->betap, 1);
   out_array(cout, "\tbleng", &p->bleng, 1);
   out_array(cout, "\tpp", p->pp, 6);
   out_array(cout, "\twp", p->wp, 12);
   out_array(cout, "\tgp", p->gp, 3);
   out_array(cout, "\tres", p->res, 6);
//#endif
   
   /* chiama funzione fortran */
   __FC_DECL__(wvefr) (p->vw,
		       &p->ome,
		       p->cc,
		       p->res,
		       &p->rm,
		       &p->ircw,
		       &p->ipc,
		       &p->omep,
		       p->tabfat,
		       &p->npft,
		       &p->beta,
		       &p->betap,
		       &p->bleng,
		       p->pp,
		       p->wp,
		       p->gp);

//#ifdef DEBUG
   cout << "out:" << endl;
   out_array(cout, "\tvw", p->vw, 6);
   out_array(cout, "\tome", &(p->ome), 1);
   out_array(cout, "\tcc", p->cc, 9);
   out_array(cout, "\trm", &p->rm, 1);
   cout << "\tircw=" << p->ircw << endl;
   cout << "\tipc=" << p->ipc << endl;
   out_array(cout, "\tomep", &p->omep, 1);
   out_array(cout, "\ttabfat", p->tabfat, 2*p->npft);
   cout << "\tnpft=" << p->npft << endl;
   out_array(cout, "\tbeta", &p->beta, 1);
   out_array(cout, "\tbetap", &p->betap, 1);
   out_array(cout, "\tbleng", &p->bleng, 1);
   out_array(cout, "\tpp", p->pp, 6);
   out_array(cout, "\twp", p->wp, 12);
   out_array(cout, "\tgp", p->gp, 3);
   out_array(cout, "\tres", p->res, 6+1);
//#endif
   
   Vec3 F(p->res);
   Vec3 M(p->res+3);
   // M += WheelOffset.Cross(F);	  

   for (int i = 1; i <= 6; i++) {
      WorkVec.fPutRowIndex(i, iIndex+i);
   }
   
   cout << "R=" << RStrut << endl;
   cout << "\tF={" << F << "}, M={" << M << "}" << endl;
   cout << "\tF={" << RStrut*F << "}, M={" << RStrut*M << "}" << endl;
   
   WorkVec.Sub(1, RStrut*F);
   WorkVec.Sub(4, RStrut*M);
   
   return WorkVec;
}

#if 0
void before_predict(const LoadableElem* pEl, 
		    VectorHandler& X,
		    VectorHandler& XP,
		    VectorHandler& XPrev,
		    VectorHandler& XPPrev)
{
   DEBUGCOUTFNAME("before_predict");
}

void after_predict(const LoadableElem* pEl, 
		   VectorHandler& X,
		   VectorHandler& XP)
{
   DEBUGCOUTFNAME("after_predict");
}

void update(LoadableElem* pEl, 
	    const VectorHandler& X,
	    const VectorHandler& XP)
{
   DEBUGCOUTFNAME("update");
}

unsigned int i_get_initial_num_dof(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_initial_num_dof");
   return 0;
}

void initial_work_space_dim(const LoadableElem* pEl, 
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
   
   // set sub-matrix indices and coefs
   
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

   // set sub-vector indices and coefs
   
   return WorkVec;
}

void set_value(const LoadableElem* pEl, VectorHandler& X, VectorHandler& XP)
{
   DEBUGCOUTFNAME("set_value");
}
   
void set_initial_value(const LoadableElem* pEl, VectorHandler& X)
{
   DEBUGCOUTFNAME("set_initial_value");
}

unsigned int i_get_num_priv_data(const LoadableElem* pEl)
{
   DEBUGCOUTFNAME("i_get_num_priv_data");
   return 0;
}

doublereal d_get_priv_data(const LoadableElem* pEl, unsigned int i)
{
   DEBUGCOUTFNAME("d_get_priv_data");
   ASSERT(pEl->iGetNumPrivData() > 0);
   if (i > pEl->iGetNumPrivData()) {
      cerr << "Module-template Elem: illegal private data index " << i << endl;      
      THROW(ErrGeneric());
   }
   
   // return i-th priv data
   return 0.;
}
#endif

static void
destroy(LoadableElem* pEl)
{
   DEBUGCOUTFNAME("destroy");
   module_wheel* p = (module_wheel *)pEl->pGetData();
   
   if (p->npft > 0 && p->tabfat != NULL) {
      SAFEDELETEARR(p->tabfat, EMmm);
   }
   
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

extern "C" void *calls = &lc;

