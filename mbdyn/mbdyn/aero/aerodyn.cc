/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2006
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

/* Elementi aerodinamici
 * 
 * - Proprieta' dell'aria:
 * elemento unico, contiene:
 *   - direzione e modulo del vento relativo. Il modulo e' associato
 *     ad un driver che ne consente la variazione, per consentire la 
 *     simulazione di un transitorio di galleria.
 *     Verra' aggiunta una variazione della direzione per simulare la raffica.
 *   - densita' dell'aria
 *   - celerita' del suono
 *   - ...
 * 
 * - Elemento di rotore:
 * classe virtuale che possiede alcuni dati topologici e geometrici di rotore
 * ed i metodi per calcolare grandezze utili agli elementi aerodinamici che
 * fanno parte di un rotore
 * 
 * - Elemento aerodinamico:
 * stazionario o quasi-stazionario, con metodo p-k di ordine 0, 1 o 2,
 * associato ad un corpo rigido, basato sulla strip theory.
 * 
 * - Elemento aerodinamico:
 * analogo al precedente, ma associato alla trave a Volumi Finiti a tre nodi
 * 
 * - Elemento aerodinamico instazionario:
 * in fase di sviluppo, modella dinamicamente alcuni stati a dare
 * il comportamento instazionario di una superficie aerodinamica
 * modellata con la strip theory
 * 
 */

#ifdef HAVE_CONFIG_H
#include <mbconfig.h>           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <aerodyn_.h>
#include <dataman.h>
#include <drive_.h>
#include <tpldrive_.h>

/*
 * Gust
 */
Gust::~Gust(void)
{
	NO_OP;
}

Gust1D::Gust1D(const Vec3& f, const Vec3& g, const doublereal& v,
		DriveCaller *pT, DriveCaller *pG)
: FrontDir(f),
GustDir(g),
dVRef(v),
Time(pT),
GustProfile(pG)
{
	ASSERT(pT != NULL);
	ASSERT(pG != NULL);
}

Gust1D::~Gust1D(void)
{
	NO_OP;
}

std::ostream&
Gust1D::Restart(std::ostream& out) const
{
	out << "front 1D, ",
		FrontDir.Write(out, ", ")
		<< ", ", GustDir.Write(out, ", ")
		<< ", " << dVRef
		<< ", ", GustProfile.pGetDriveCaller()->Restart(out);
	return out;
}

Vec3
Gust1D::GetVelocity(const Vec3& X) const
{
	Vec3 V;
	GetVelocity(X, V);
	return V;
}

void
Gust1D::GetVelocity(const Vec3& X, Vec3& V) const
{
	doublereal x = FrontDir*X + dVRef*Time.dGet();
	doublereal v = GustProfile.dGet(x);
	V = GustDir*v;
}

/* AirProperties - begin */

AirProperties::AirProperties(const TplDriveCaller<Vec3>* pDC,
		Gust *pG, flag fOut)
: Elem(1, Elem::AIRPROPERTIES, fOut),
InitialAssemblyElem(1, Elem::AIRPROPERTIES, fOut),
TplDriveOwner<Vec3>(pDC),
Velocity(0.),
pGust(pG)
{
	NO_OP;
}
   
AirProperties::~AirProperties(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
AirProperties::Restart(std::ostream& out) const
{
	pGetDriveCaller()->Restart(out);
	if (pGust) {
		out << ", ", pGust->Restart(out);
	}
	return out << ';' << std::endl;	   
}
   
/* Dimensioni del workspace */
void
AirProperties::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
AirProperties::AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal /* dCoef */ , 
		const VectorHandler& /* XCurr */ ,
		const VectorHandler& /* XPrimeCurr */ )
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}
   
/* assemblaggio residuo */
SubVectorHandler&
AirProperties::AssRes(SubVectorHandler& WorkVec,
		doublereal /* dCoef */ ,
		const VectorHandler& /* XCurr */ , 
		const VectorHandler& /* XPrimeCurr */ )
{
  	WorkVec.Resize(0);
  	
	/* Approfitto del fatto che AirProperties viene aggiornato prima 
	 * degli altri elementi (vedi l'enum Elem::Type e la sequenza di
	 * assemblaggio) per fargli calcolare Velocity una volta per tutte.
	 * Quindi, quando viene chiamata GetVelocity(void), 
	 * questa restituisce un reference alla velocita' con il
	 * minimo overhead */
	Velocity = Get();
  	
	return WorkVec;
}
   
   /* Numero di GDL iniziali */
unsigned int
AirProperties::iGetInitialNumDof(void) const
{
	return 0;
}
     
/* Dimensioni initiali del workspace */
void
AirProperties::InitialWorkSpaceDim(integer* piNumRows,
		integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}
   
/* assemblaggio jacobiano */
VariableSubMatrixHandler& 
AirProperties::InitialAssJac(VariableSubMatrixHandler& WorkMat,	  
		const VectorHandler& /* XCurr */ )
{
	WorkMat.SetNullMatrix();
	return WorkMat;
}
   
/* assemblaggio residuo */
SubVectorHandler&
AirProperties::InitialAssRes(SubVectorHandler& WorkVec,
		const VectorHandler& XCurr)
{
	/* Chiama AssRes, che si limita ad aggiornare la velocita' */
	return AssRes(WorkVec, 1, XCurr, XCurr);
}

void
AirProperties::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		OH.AirProps() << std::setw(8) << 0 << " "
			<< dGetAirDensity(Zero3) << " "
			<< dGetSoundSpeed(Zero3) << " "
			<< GetVelocity(Zero3) << std::endl;
	}
}
 
Vec3
AirProperties::GetVelocity(const Vec3& X) const
{
	Vec3 V(0.);

	if (pGust) {
		pGust->GetVelocity(X, V);
	}

	return V += Velocity;
}

void
AirProperties::GetVelocity(const Vec3& X, Vec3& V) const
{
	V = GetVelocity(X);
}

/* AirProperties - end */


/* BasicAirProperties - begin */

BasicAirProperties::BasicAirProperties(const TplDriveCaller<Vec3>* pDC,
		DriveCaller *pRho, doublereal dSS, Gust *pG, flag fOut)
: Elem(1, Elem::AIRPROPERTIES, fOut),
AirProperties(pDC, pG, fOut),
pAirDensity(pRho),
dSoundSpeed(dSS)
{
	NO_OP;
}
   
BasicAirProperties::~BasicAirProperties(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
BasicAirProperties::Restart(std::ostream& out) const
{
	out << "  air properties: ",
		pAirDensity->Restart(out) << ", " << dSoundSpeed << ", ";
	return AirProperties::Restart(out);
}
   
doublereal
BasicAirProperties::dGetAirDensity(const Vec3& /* X */ ) const
{
	return pAirDensity->dGet();
}
   
doublereal
BasicAirProperties::dGetAirPressure(const Vec3& /* X */ ) const
{
	return -1.;
}

doublereal
BasicAirProperties::dGetAirTemperature(const Vec3& /* X */ ) const
{
	return -1.;
}

doublereal
BasicAirProperties::dGetSoundSpeed(const Vec3& /* X */ ) const
{
	return dSoundSpeed;
}

bool
BasicAirProperties::GetAirProps(const Vec3& X, doublereal& rho,
		doublereal& c, doublereal& p, doublereal& T) const
{
	/* FIXME */
	rho = dGetAirDensity(X);
	c = dGetSoundSpeed(X);
	p = -1.;
	T = -1.;

	return true;
}

/* BasicAirProperties - end */


/* StdAirProperties - begin */

StdAirProperties::StdAirProperties(const TplDriveCaller<Vec3>* pDC,
		 doublereal PRef_, DriveCaller *RhoRef_, doublereal TRef_,
		 doublereal a_, doublereal R_, doublereal g0_,
		 doublereal z0_, doublereal z1_, doublereal z2_,
		 Gust *pG, flag fOut)
: Elem(1, Elem::AIRPROPERTIES, fOut),
AirProperties(pDC, pG, fOut),
PRef(PRef_),
RhoRef(RhoRef_),
TRef(TRef_),
a(a_),
R(R_),
g0(g0_),
z0(z0_),
z1(z1_),
z2(z2_)
{
	ASSERT(PRef > 0.);
	ASSERT(RhoRef != NULL);
	ASSERT(TRef > 0.);
	ASSERT(R > 0.);
	ASSERT(g0 > 0.);
	ASSERT(z1 > 0.);
	ASSERT(z2 > z1);
}
   
StdAirProperties::~StdAirProperties(void)
{
	NO_OP;
}

/* Scrive il contributo dell'elemento al file di restart */
std::ostream&
StdAirProperties::Restart(std::ostream& out) const
{
	out << "  air properties: std, "
		<< PRef << ", ",
		RhoRef->Restart(out) << ", "
		<< TRef << ", "
		<< a << ", "
		<< R << ", "
		<< g0 << ", "
		<< z1 << ", "
		<< z2 << ", ";
	if (z0 != 0.) {
		out << "reference altitude, " << z0 << ", ";
	}
	return AirProperties::Restart(out);
}
   
doublereal
StdAirProperties::dGetAirDensity(const Vec3& X) const
{
	doublereal rho, c, p, T;

	GetAirProps(X, rho, c, p, T);

	return rho;
}
   
doublereal
StdAirProperties::dGetAirPressure(const Vec3& X) const
{
	doublereal rho, c, p, T;

	GetAirProps(X, rho, c, p, T);

	return p;
}

doublereal
StdAirProperties::dGetAirTemperature(const Vec3& X) const
{
	doublereal rho, c, p, T;

	GetAirProps(X, rho, c, p, T);

	return T;
}

doublereal
StdAirProperties::dGetSoundSpeed(const Vec3& X) const
{
	doublereal rho, c, p, T;

	GetAirProps(X, rho, c, p, T);

	return c;
}

bool
StdAirProperties::GetAirProps(const Vec3& X, doublereal& rho,
		doublereal& c, doublereal& p, doublereal& T) const
{
	/* FIXME */
	doublereal z = X.dGet(3) + z0;

	if (z < z1) {
		/* regione del gradiente (troposfera) */
		T = TRef + a * z;
		p = PRef * pow(T / TRef, -g0 / (a * R));
		doublereal rhoRef = RhoRef->dGet();
		if (rhoRef < 0.) {
			silent_cerr("illegal reference density "
				<< rhoRef << std::endl);
			throw ErrGeneric();
		}
		rho = rhoRef * pow(T / TRef, -(g0 / (a * R) + 1));
	} else {
		/* regione a T = const. (stratosfera) */
		T = TRef + a * z1;
		doublereal p1 = PRef * pow(T / TRef, -g0 / (a * R));
		p = p1*exp(-g0 / (R * T) * (z - z1));
		doublereal rhoRef = RhoRef->dGet();
		if (rhoRef < 0.) {
			silent_cerr("illegal reference density "
				<< rhoRef << std::endl);
			throw ErrGeneric();
		}
		doublereal rho1 = rhoRef * pow(T / TRef, -(g0 / (a * R) + 1));
		rho = rho1 * exp(-(g0 / (R*T) * (z - z1)));
	}
	c = sqrt(1.4 * R * T);

	return true;
}

/* StdAirProperties - end */

static void
ReadAirstreamData(DataManager *pDM, MBDynParser& HP,
		TplDriveCaller<Vec3>*& pDC, Gust*& pG)
{
	ASSERT(pDC == NULL);
	ASSERT(pG == NULL);

	/* Driver multiplo */
     	pDC = ReadTplDrive(pDM, HP, Vec3(0.));
	if (HP.IsKeyWord("gust")) {
		if (HP.IsKeyWord("front" "1d")) {
			/* front direction */
			Vec3 f = HP.GetVecAbs(AbsRefFrame);

			/* gust velocity direction */
			Vec3 g = HP.GetVecAbs(AbsRefFrame);

			/* reference velocity */
			doublereal v = HP.GetReal();

			/* time drive caller
			 * FIXME: not needed if v = 0 */
			DriveCaller *pT = NULL;
			SAFENEWWITHCONSTRUCTOR(pT, TimeDriveCaller,
					TimeDriveCaller(pDM->pGetDrvHdl()));

			/* gust profile drive caller */
			DriveCaller *pP = HP.GetDriveCaller();

			/* gust */
			SAFENEWWITHCONSTRUCTOR(pG, Gust1D,
					Gust1D(f, g, v, pT, pP));

		} else {
			silent_cerr("unknown gust type at line "
				<< HP.GetLineData() << std::endl);
			throw ErrGeneric();
		}
	}
}

Elem *
ReadAirProperties(DataManager* pDM, MBDynParser& HP)
{
	Elem *pEl = NULL;

	if (HP.IsKeyWord("std")) {
		doublereal PRef(0.);
		doublereal rhoRef(0.);
		DriveCaller *RhoRef(NULL);
		doublereal TRef(0.);
		doublereal a(0.);
		doublereal R(0.);
		doublereal g0(0.);
		doublereal z0(0.);
		doublereal z1(0.);
		doublereal z2(0.);

		bool Std = false;

		if (HP.IsKeyWord("SI")) {
			Std = true;

			PRef = 101325.;		/* Pa */
			rhoRef = 1.2250;	/* kg/m^3 */
			TRef = 288.16;		/* K */
			a = -6.5e-3;		/* K/m */
			R = 287.;		/* J/kgK */
			g0 = 9.81; 		/* m/s^2 */
			z1 = 11000.; 		/* m */
			z2 = 25000.; 		/* m */

		} else if (HP.IsKeyWord("british")) {
			Std = true;
			PRef = 2116.2; 		/* lb/ft^2 */
			rhoRef = 0.002377;	/* slug/ft3 */
			TRef = 518.69;		/* R */
			a = -3.566e-3;		/* R/ft */
			R = 1716;		/* ft lb/slug R */
			g0 = 32.17;		/* ft/s^2 */
			z1 = 36089;		/* ft */
			z2 = 82021;		/* ft */

		} else {
			PRef = HP.GetReal();
			if (PRef <= 0.) {
				silent_cerr("illegal reference "
					"pressure" << PRef 
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			RhoRef = HP.GetDriveCaller();
			/* FIXME: we need to do runtime checks ... */

			TRef = HP.GetReal();
			if (TRef <= 0.) {
				silent_cerr("illegal reference "
					"temperature " << TRef 
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			a = HP.GetReal();
			if (a >= 0.) {
				/* FIXME: should we leave this free? */
				silent_cerr("illegal temperature gradient "
					<< a << " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			R = HP.GetReal();
			if (R <= 0.) {
				silent_cerr("illegal gas constant " << R
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			g0 = HP.GetReal();
			if (g0 <= 0.) {
				silent_cerr("illegal reference "
					"gravity acceleration " << g0
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			z1 = HP.GetReal();
			if (z1 <= 0.) {
				silent_cerr("illegal troposphere altitude "
					<< z1
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}

			z2 = HP.GetReal();
			if (z2 <= z1) {
				silent_cerr("illegal stratosphere altitude "
					<< z2
					<< " at line " << HP.GetLineData()
					<< std::endl);
				throw ErrGeneric();
			}
		}
	
		if (Std) {
			if (HP.IsKeyWord("temperature" "deviation")) {
				doublereal T = HP.GetReal();

				if (TRef + T <= 0.) {
					silent_cerr("illegal "
						"temperature deviation " << T
						<< " at line " 
						<< HP.GetLineData()
						<< std::endl);
					throw ErrGeneric();
				}

				/*
				 * trasformazione isobara applicata
				 * all'equazione
				 * di stato: rho * R * T = cost
				 *
				 * rho = rho_0 T_0 / T
				 */
				rhoRef *= TRef / (TRef + T);
				TRef += T;
			}
			
			SAFENEWWITHCONSTRUCTOR(RhoRef, ConstDriveCaller,
					ConstDriveCaller(rhoRef));
		}

		if (HP.IsKeyWord("reference" "altitude")) {
			z0 = HP.GetReal();
		}

	     	/* Driver multiplo */	   
	     	TplDriveCaller<Vec3>* pDC = NULL;
		Gust *pG = NULL;
		ReadAirstreamData(pDM, HP, pDC, pG);
	     	flag fOut = pDM->fReadOutput(HP, Elem::AIRPROPERTIES);
	     
	     	SAFENEWWITHCONSTRUCTOR(pEl, 
				StdAirProperties,
				StdAirProperties(pDC, 
					PRef, RhoRef, TRef, a, R, g0, z0, z1, z2,
					pG, fOut));
	} else {
		/* Legacy: density and sound celerity at one altitude;
		 * no altitude dependency */
		
		DriveCaller *pRho = HP.GetDriveCaller();

	     	doublereal dSS = HP.GetReal();
	     	DEBUGLCOUT(MYDEBUG_INPUT, "Sound speed: " << dSS << std::endl);
	     	if (dSS <= 0.) {
			silent_cerr("illegal null or negative sound speed "
				"at line " << HP.GetLineData() << std::endl);
		
			throw DataManager::ErrGeneric();
	     	}	      
	     
	     	/* Driver multiplo */	   
	     	TplDriveCaller<Vec3>* pDC = NULL;
		Gust *pG = NULL;
		ReadAirstreamData(pDM, HP, pDC, pG);
	     	flag fOut = pDM->fReadOutput(HP, Elem::AIRPROPERTIES);
	     
	     	SAFENEWWITHCONSTRUCTOR(pEl, 
				BasicAirProperties,
				BasicAirProperties(pDC, pRho, dSS,
					pG, fOut));
	}

	return pEl;
}


/* AirPropOwner - begin */

AirPropOwner::AirPropOwner(void)
: pAirProperties(NULL)
{
	NO_OP;
}

AirPropOwner::~AirPropOwner(void)
{
	NO_OP;
}
   
void
AirPropOwner::PutAirProperties(const AirProperties* pAP)
{
	ASSERT(pAirProperties == NULL);

	(AirProperties*&)pAirProperties = (AirProperties*)pAP;
}
   
flag
AirPropOwner::fGetAirVelocity(Vec3& Velocity, const Vec3& X) const
{
	if (pAirProperties == NULL) {
		return 0;
	}

	Velocity = pAirProperties->GetVelocity(X);
	return 1;
}
   
doublereal
AirPropOwner::dGetAirDensity(const Vec3& X) const
{
	return pAirProperties->dGetAirDensity(X);
}
   
doublereal
AirPropOwner::dGetAirPressure(const Vec3& X) const
{
	return pAirProperties->dGetAirPressure(X); 
}
   
doublereal
AirPropOwner::dGetAirTemperature(const Vec3& X) const
{
	return pAirProperties->dGetAirTemperature(X); 
}
   
doublereal
AirPropOwner::dGetSoundSpeed(const Vec3& X) const
{
	return pAirProperties->dGetSoundSpeed(X); 
}

bool
AirPropOwner::GetAirProps(const Vec3& X, doublereal& rho,
		doublereal& c, doublereal& p, doublereal& T) const
{
	return pAirProperties->GetAirProps(X, rho, c, p, T);
}

/* AirPropOwner - end */


/* AerodynamicElem - begin */

AerodynamicElem::AerodynamicElem(unsigned int uL,
		AerodynamicElem::Type T, flag fOut)
: Elem(uL, Elem::AERODYNAMIC, fOut),
AeroT(T)
{
	NO_OP; 
}
   
AerodynamicElem::~AerodynamicElem(void)
{
	NO_OP; 
}
   
bool
AerodynamicElem::NeedsAirProperties(void) const
{
	return true;
}

const Rotor *
AerodynamicElem::pGetRotor(void) const
{
	return NULL;
}
   
/* AerodynamicElem - end */

