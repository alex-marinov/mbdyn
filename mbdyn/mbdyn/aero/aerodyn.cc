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

#include <aerodyn.h>

/* AirProperties - begin */

AirProperties::AirProperties(const TplDriveCaller<Vec3>* pDC,
		DriveCaller *pRho, doublereal dSS, flag fOut)
: Elem(1, Elem::AIRPROPERTIES, fOut),
InitialAssemblyElem(1, Elem::AIRPROPERTIES, fOut),
TplDriveOwner<Vec3>(pDC),    
Velocity(0.),
pAirDensity(pRho),
dSoundSpeed(dSS)
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
	return out << "  air properties: ",
		pAirDensity->Restart(out) << ", " << dSoundSpeed << ", ",
		pGetDriveCaller()->Restart(out) << ';' << std::endl;	   
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
     
const Vec3&
AirProperties::GetVelocity(const Vec3& /* X */ ) const
{
	return Velocity;
}
   
doublereal
AirProperties::dGetAirDensity(const Vec3& /* X */ ) const
{
	return pAirDensity->dGet();
}
   
doublereal
AirProperties::dGetAirPressure(const Vec3& /* X */ ) const
{
	return 0.;
}

doublereal
AirProperties::dGetAirTemperature(const Vec3& /* X */ ) const
{
	return 0.;
}

doublereal
AirProperties::dGetSoundSpeed(const Vec3& /* X */ ) const
{
	return dSoundSpeed;
}

/* AirProperties - end */


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

