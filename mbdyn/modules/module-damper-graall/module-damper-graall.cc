/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>
#include <cfloat>
#include <fstream>

#include "dataman.h"
#include "drive_.h"
#include "constltp.h"

/* GRAALLDamperConstitutiveLaw - begin */

extern "C" {
int __FC_DECL__(dmpfr) (doublereal* DL,
	doublereal* DLDT,
	doublereal* V,
	doublereal* TBDMR,
	doublereal* TBDMA,
	doublereal* TBCVR,
	doublereal* TBCVA,
	doublereal* VETVIS,
	integer* NPDMR,
	integer* NPDMA,
	integer* NPCVR,
	integer* NPCVA,
	integer* NPVT,
	doublereal* FTOT,
	doublereal* FVISDL,
	doublereal* FELDL);
}

class GRAALLDamperConstitutiveLaw
: public ConstitutiveLaw<doublereal, doublereal> {
private:
	enum {
		RLA = 1,
		PNTR = 2,
		P0 = 3,
		A0 = 4,
		RL0 = 5,
		PLTREX = 6,
		TCOST = 7,
		TLIN = 8,
		TEMPORIF = 9,
		VELRIF = 10,
		CORONA = 11,
		TRIFVUOTO = 12,
		T = 13,
		TRIFTEMP = 14,
		RNURIF = 15,
		ASTNT = 16,
		R0 = 17,
		TIME = 18
	};

	doublereal v[TIME];

	/*
	      RLA = 0,     // lunghezza iniziale ammortizzatore
	      PNTR;    // coeff. di interazione corsa
	      P0;      // pressione di precarica
	      A0;      // sezione gas
	      RL0;     // lunghezza iniziale camera gas
	      PLTREX   // esponente politropica
	      tcost;
	      tlin;    // legge lineare coeff attrito con corsa
	C     temporif   tempo di riferimento per profilo attrito
	C     velrif     velocita' di riferimento per profilo attrito
	C     corona     area corona per effetto risucchio
	C     trifvuoto  tempo di riferimento per profilo attrito
	C     T          temperatura olio in passaggio
	C     triftemp   tempo di riferimento per effetto temperatura olio
	C     RNURIF     cvis di riferimento
	C     SWITCH     0 ammortizzatore passivo, >0 modo di controllo
	C     ASTNT      sezione stantuffo olio
	C	  RO         densita' olio

	C     TIME    tempo dall'inizio della compressione dell'ammortizzatore

	C     TBDMR   tabella sezione trafilamento ritorno  NPDMR punti
	C     TBDMA   tabella sezione trafilamento andata   NPDMA punti
	C     TBCVR   tabella coeff. efflusso ritorno       NPCVR punti
	C     TBCVA   tabella coeff. efflusso ritorno       NPCVA punti
	C     VETVIS  tabella viscosita'/temperatura        NPVT punti
	 */

	integer NPDMR,NPDMA,NPCVR,NPCVA,NPVT;
	doublereal* TBDMR;
	doublereal* TBDMA;
	doublereal* TBCVR;
	doublereal* TBCVA;
	doublereal* VETVIS;

	doublereal dT0;
	DriveOwner Time;
	flag fWorking;

public:
	GRAALLDamperConstitutiveLaw(const DriveCaller* pDC,
		const doublereal& rla,
		const char* const filename)
	: ConstitutiveLaw<doublereal, doublereal>(),
	NPDMR(0), NPDMA(0), NPCVR(0), NPCVA(0), NPVT(0),
	TBDMR(NULL), TBDMA(NULL), TBCVR(NULL), TBCVA(NULL), VETVIS(NULL),
	dT0(0.), Time(pDC), fWorking(0)
	{
		char s[1024];

		/* Si assume che il drive che viene passato sia in realta'
		 * un drive di tipo Time, che verra' usato per conoscere
		 * il tempo della simulazione */

		/* apre il file di ingresso in formato graal */
		std::ifstream in(filename);

		/* legge alcuni parametri */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s
			<< "\")" << std::endl);

		in >> v[P0-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading P0 (=" << v[P0-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[A0-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading A0 (=" << v[A0-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[PNTR-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading PNTR (=" << v[PNTR-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[RL0-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading RL0 (=" << v[RL0-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[PLTREX-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading PLTREX (=" << v[PLTREX-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Skipping ... (remark: \"" << s
			<< "\")" << std::endl);

		in >> v[ASTNT-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading ASTNT (=" << v[ASTNT-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[R0-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading R0 (=" << v[R0-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[TCOST-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TCOST (=" << v[TCOST-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[TLIN-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TLIN (=" << v[TLIN-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[TEMPORIF-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading A0 (=" << v[TEMPORIF-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[VELRIF-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading VELRIF (=" << v[VELRIF-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[CORONA-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading CORONA (=" << v[CORONA-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[TRIFVUOTO-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TRIFVUOTO (=" << v[TRIFVUOTO-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[RNURIF-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading RNURIF (=" << v[RNURIF-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		in >> v[TRIFTEMP-1];
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TRIFTEMP (=" << v[TRIFTEMP-1]
			<< ", remark: \"" << s << "\")" << std::endl);

		/* legge le tabelle */

		/* legge TBDMR */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBDMR (remark: \"" << s
			<< "\")" << std::endl);

		in >> NPDMR;
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "NPDMR=" << NPDMR
			<< " (remark: \"" << s << "\")" << std::endl);

		SAFENEWARR(TBDMR, doublereal, 2*NPDMR);

		for (integer i = 0; i < NPDMR; i++) {
			in >> TBDMR[i];
			in >> TBDMR[NPDMR+i];
			in.getline(s, sizeof(s));
			DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBDMR[" << i+1
				<< "]=(" << TBDMR[i]
				<< "," << TBDMR[NPDMR+i]
				<< ") (remark: \"" << s << "\")" << std::endl);
		}

		/* legge TBDMA */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBDMA (remark: \"" << s
			<< "\")" << std::endl);

		in >> NPDMA;
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "NPDMA=" << NPDMA
			<< " (remark: \"" << s << "\")" << std::endl);

		SAFENEWARR(TBDMA, doublereal, 2*NPDMA);

		for (integer i = 0; i < NPDMA; i++) {
			in >> TBDMA[i];
			in >> TBDMA[NPDMA+i];
			in.getline(s, sizeof(s));
			DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBDMA[" << i+1
				<< "]=(" << TBDMA[i]
				<< "," << TBDMA[NPDMA+i]
				<< ") (remark: \"" << s << "\")" << std::endl);
		}

		/* legge TBCVR */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBCVR (remark: \"" << s
			<< "\")" << std::endl);

		in >> NPCVR;
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "NPCVR=" << NPCVR
			<< " (remark: \"" << s << "\")" << std::endl);

		SAFENEWARR(TBCVR, doublereal, 2*NPCVR);

		for (integer i = 0; i < NPCVR; i++) {
			in >> TBCVR[i];
			in >> TBCVR[NPCVR+i];
			in.getline(s, sizeof(s));
			DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBCVR[" << i+1
				<< "]=(" << TBCVR[i]
				<< "," << TBCVR[NPCVR+i]
				<< ") (remark: \"" << s << "\")" << std::endl);
		}

		/* legge TBCVA */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBCVA (remark: \"" << s
			<< "\")" << std::endl);

		in >> NPCVA;
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "NPCVA=" << NPCVA
			<< " (remark: \"" << s << "\")" << std::endl);

		SAFENEWARR(TBCVA, doublereal, 2*NPCVA);

		for (integer i = 0; i < NPCVA; i++) {
			in >> TBCVA[i];
			in >> TBCVA[NPCVA+i];
			in.getline(s, sizeof(s));
			DEBUGLCOUT(MYDEBUG_INPUT, "Reading TBCVA[" << i+1
				<< "]=(" << TBCVA[i]
				<< "," << TBCVA[NPCVA+i]
				<< ") (remark: \"" << s << "\")" << std::endl);
		}

		/* legge VETVIS */
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "Reading VETVIS (remark: \"" << s
			<< "\")" << std::endl);

		in >> NPVT;
		in.getline(s, sizeof(s));
		DEBUGLCOUT(MYDEBUG_INPUT, "NPVT=" << NPVT
			<< " (remark: \"" << s << "\")" << std::endl);

		SAFENEWARR(VETVIS, doublereal, 2*NPVT);

		for (integer i = 0; i < NPVT; i++) {
			in >> VETVIS[i];
			in >> VETVIS[NPVT+i];
			in.getline(s, sizeof(s));
			DEBUGLCOUT(MYDEBUG_INPUT, "Reading VETVIS[" << i+1
				<< "]=(" << VETVIS[i]
				<< "," << VETVIS[NPVT+i]
				<< ") (remark: \"" << s << "\")" << std::endl);
		}

		/* fine */
		in.close();

		/* Inizializza altri dati */
		v[RLA-1] = rla;   // passato come argomento, per ora
		v[T-1] = 273.16;  // T assoluta ?!?
		v[TIME-1] = Time.dGet();
	};

	virtual ~GRAALLDamperConstitutiveLaw(void) {
		SAFEDELETEARR(VETVIS);
		SAFEDELETEARR(TBCVA);
		SAFEDELETEARR(TBCVR);
		SAFEDELETEARR(TBDMA);
		SAFEDELETEARR(TBDMR);
	};

	ConstLawType::Type GetConstLawType(void) const {
		return ConstLawType::VISCOELASTIC;
	};

	virtual ConstitutiveLaw<doublereal, doublereal>* pCopy(void) const {
		return 0;
	};

	virtual std::ostream& Restart(std::ostream& out) const {
		return out << "GRAALLDamperConstitutiveLaw not implemented yet!" << std::endl;
	};

	virtual void Update(const doublereal& Eps, const doublereal& EpsPrime = 0.) {
		Epsilon = Eps;
		EpsilonPrime = EpsPrime;

		/* se non e' attivo, esegue un test sulla velocita' relativa;
		 * se si sta muovendo, allora setta il flag di stato su attivo,
		 * e calcola il tempo di inizio */
		if (fWorking == 0) {
			if (fabs(EpsPrime) < 1.e-6) {
				F = 0.;
				FDE = 0.;
				FDEPrime = 0.;
				return;
			} else {
				fWorking = 1;
				dT0 = Time.dGet();
			}
		}

		v[TIME-1] = Time.dGet() - dT0;

		doublereal E = Epsilon*v[RLA-1];
		doublereal EP = EpsilonPrime*v[RLA-1];

		__FC_DECL__(dmpfr) (
			&E,
			&EP,
			v,
			TBDMR,
			TBDMA,
			TBCVR,
			TBCVA,
			VETVIS,
			&NPDMR,
			&NPDMA,
			&NPCVR,
			&NPCVA,
			&NPVT,
			&F,
			&FDE,
			&FDEPrime);
	};
};

/* GRAALLDamperConstitutiveLaw - end */

struct GRAALLDamperCLR : public ConstitutiveLawRead<doublereal, doublereal> {
	virtual ConstitutiveLaw<doublereal, doublereal> *
	Read(const DataManager* pDM, MBDynParser& HP, ConstLawType::Type& CLType) {
		ConstitutiveLaw<doublereal, doublereal>* pCL = 0;

		CLType = ConstLawType::VISCOELASTIC;

		const char* filename = HP.GetFileName();
		DEBUGCOUT("Graall damper input file: \""
				<< filename << "\"" << std::endl);

		doublereal rla = HP.GetReal();
		DEBUGCOUT("Reference length: " << rla << std::endl);

		DriveCaller* pDC = NULL;
		SAFENEWWITHCONSTRUCTOR(pDC,
				TimeDriveCaller,
				TimeDriveCaller(pDM->pGetDrvHdl()));

		typedef GRAALLDamperConstitutiveLaw L;
		SAFENEWWITHCONSTRUCTOR(pCL,
				L,
				L(pDC, rla, filename));

		return pCL;
	};
};

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif

	ConstitutiveLawRead<doublereal, doublereal> *rf1D
		= new GRAALLDamperCLR;
	if (!SetCL1D("GRAALL" "damper", rf1D)) {
		delete rf1D;

		silent_cerr("GRAALLDamperConstitutiveLaw: "
			"module_init(" << module_name << ") "
			"failed" << std::endl);

		return -1;
	}

	return 0;
}

