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

#ifndef AERODYN_H
#define AERODYN_H

#include <elem.h>
#include <tpldrive.h>

extern const char* psAeroNames[];


/* AirProperties - begin */

class AirProperties 
: virtual public Elem, public InitialAssemblyElem, public TplDriveOwner<Vec3> {
 protected:
   mutable Vec3 Velocity;
   
 public:
   AirProperties(const TplDriveCaller<Vec3>* pDC, flag fOut);
   
   virtual ~AirProperties(void);

   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual Elem::Type GetElemType(void) const { 
      return Elem::AIRPROPERTIES; 
   };   
      
   /* funzioni di servizio */

   /* Il metodo iGetNumDof() serve a ritornare il numero di gradi di liberta'
    * propri che l'elemento definisce. Non e' virtuale in quanto serve a 
    * ritornare 0 per gli elementi che non possiedono gradi di liberta'.
    * Viene usato nella costruzione dei DofOwner e quindi deve essere 
    * indipendente da essi. In genere non comporta overhead in quanto il 
    * numero di dof aggiunti da un tipo e' una costante e non richede dati 
    * propri.
    * Il metodo pGetDofOwner() ritorna il puntatore al DofOwner dell'oggetto.
    * E' usato da tutti quelli che agiscono direttamente sui DofOwner.
    * Non e' virtuale in quanto ritorna NULL per tutti i tipi che non hanno
    * dof propri.
    * Il metodo SetDof() ritorna, per ogni dof dell'elemento, l'ordine.
    * E' usato per completare i singoli Dof relativi all'elemento.
    */
   
   /* funzioni proprie */
   
   /* Dimensioni del workspace */
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ , 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ );
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ , 
				    const VectorHandler& /* XPrimeCurr */ );
   
   /* Numero di GDL iniziali */
   virtual unsigned int iGetInitialNumDof(void) const;
     
   /* Dimensioni initiali del workspace */
   virtual void InitialWorkSpaceDim(integer* piNumRows,
		   integer* piNumCols) const;
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,	  
	    const VectorHandler& /* XCurr */ );
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
					   const VectorHandler& XCurr);

   virtual void Output(OutputHandler&) const;

   /*
    * Deprecated; use GetAirProps instead
    */
   virtual const Vec3& GetVelocity(const Vec3& /* X */ ) const = 0;
   
   virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const = 0;
   
   virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const = 0;

   virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const = 0;

   virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const = 0;

   /*
    * End of deprecated; use GetAirProps instead
    */

   virtual bool GetAirProps(const Vec3& X, Vec3& V, doublereal& rho,
		   doublereal& c, doublereal& p, doublereal& T) const = 0;
};

/* AirProperties - end */


/* BasicAirProperties - begin */

class BasicAirProperties 
: virtual public Elem, public AirProperties {
 protected:
   DriveCaller *pAirDensity;
   doublereal dSoundSpeed;
   
 public:
   BasicAirProperties(const TplDriveCaller<Vec3>* pDC,
		 DriveCaller *pRho, doublereal dSS, flag fOut);
   
   virtual ~BasicAirProperties(void);

   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /*
    * Deprecated; use GetAirProps instead
    */
   virtual const Vec3& GetVelocity(const Vec3& /* X */ ) const;
   
   virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const;
   
   virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const;

   virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const;

   virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const;

   /*
    * End of deprecated; use GetAirProps instead
    */

   virtual bool GetAirProps(const Vec3& X, Vec3& V, doublereal& rho,
		   doublereal& c, doublereal& p, doublereal& T) const;
};

/* BasicAirProperties - end */


/* StdAirProperties - begin */

class StdAirProperties 
: virtual public Elem, public AirProperties {
 protected:
   doublereal PRef;
   DriveCaller *RhoRef;
   doublereal TRef;
   doublereal a;
   doublereal R;
   doublereal g0;
   doublereal z1;
   doublereal z2;
   
 public:
   StdAirProperties(const TplDriveCaller<Vec3>* pDC,
		 doublereal PRef_, DriveCaller *RhoRef_, doublereal TRef_,
		 doublereal a_, doublereal R_, doublereal g0_,
		 doublereal z1_, doublereal z2_, flag fOut);
   
   virtual ~StdAirProperties(void);

   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual std::ostream& Restart(std::ostream& out) const;
   
   /*
    * Deprecated; use GetAirProps instead
    */
   virtual const Vec3& GetVelocity(const Vec3& /* X */ ) const;
   
   virtual doublereal dGetAirDensity(const Vec3& /* X */ ) const;
   
   virtual doublereal dGetAirPressure(const Vec3& /* X */ ) const;

   virtual doublereal dGetAirTemperature(const Vec3& /* X */ ) const;

   virtual doublereal dGetSoundSpeed(const Vec3& /* X */ ) const;

   /*
    * End of deprecated; use GetAirProps instead
    */

   virtual bool GetAirProps(const Vec3& X, Vec3& V, doublereal& rho,
		   doublereal& c, doublereal& p, doublereal& T) const;
};

/* StdAirProperties - end */


/* AirPropOwner - begin */

class AirPropOwner {
 protected:
   const AirProperties* pAirProperties;
 public:
   AirPropOwner(void);
   virtual ~AirPropOwner(void);
   
   virtual void PutAirProperties(const AirProperties* pAP);
   
   /*
    * Deprecated; use GetAirProps instead
    */
   virtual flag fGetAirVelocity(Vec3& Velocity, const Vec3& X) const;
   
   virtual doublereal dGetAirDensity(const Vec3& X) const;
   
   virtual doublereal dGetAirPressure(const Vec3& X) const;
   
   virtual doublereal dGetAirTemperature(const Vec3& X) const;
   
   virtual doublereal dGetSoundSpeed(const Vec3& X) const;

   /*
    * End of deprecated; use GetAirProps instead
    */

   virtual bool GetAirProps(const Vec3& X, Vec3& V, doublereal& rho,
		   doublereal& c, doublereal& p, doublereal& T) const;
};

/* AirPropOwner - end */


/* AerodynamicElem - begin */

class AerodynamicElem : virtual public Elem, public AirPropOwner {
 public:
   /* Tipi di elementi aerodinamici */
   enum Type {
      UNKNOWN = -1,

	ROTOR = 0,
	AEROMODAL,
	AERODYNAMICBODY,
	AERODYNAMICBEAM,

	AERODYNAMICLOADABLE,
	
	LASTAEROTYPE
   };
 private:
   AerodynamicElem::Type AeroT;
   
 protected:
   
 public:
   AerodynamicElem(unsigned int uL, AerodynamicElem::Type T, flag fOut);
   
   virtual ~AerodynamicElem(void);
   
   /* Tipo di elemento aerodinamico */
   virtual AerodynamicElem::Type GetAerodynamicElemType(void) const {
	   return AeroT;
   };

   /* Consente di effettuare un casting sicuro da Elem* a AerodynamicElem* */
   virtual AerodynamicElem* pGetAerodynamicElem(void) const {
	   return (AerodynamicElem *)this;
   };

   virtual bool NeedsAirProperties(void) const;

   virtual const Rotor *pGetRotor(void) const;
};

/* AerodynamicElem - end */

#endif /* AERODYN_H */

