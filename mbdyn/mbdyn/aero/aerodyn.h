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

#include "elem.h"
#include "tpldrive.h"
#include "memmans.h"


/* Tipi di elementi aerodinamici */
class AeroType {
 public:
   enum Type {
      UNKNOWN = -1,

	ROTOR = 0,
	AERODYNAMICBODY,
	AERODYNAMICBEAM,
        AERODYNAMICMODAL,
	
	LASTAEROTYPE
   };
};

extern const char* psAeroNames[];


/* AirProperties - begin */

class AirProperties 
: virtual public Elem, public InitialAssemblyElem, public TplDriveOwner<Vec3> {
 protected:
   Vec3 Velocity;
   doublereal dAirDensity;
   doublereal dSoundSpeed;
   
 public:
   AirProperties(const TplDriveCaller<Vec3>* pDC,
		 doublereal dRho, doublereal dSS, flag fOut)
     : Elem(1, ElemType::AIRPROPERTIES, fOut),
     InitialAssemblyElem(1, ElemType::AIRPROPERTIES, fOut),
     TplDriveOwner<Vec3>(pDC),    
     Velocity(0.), dAirDensity(dRho), dSoundSpeed(dSS) {
	NO_OP;
     };
   
   virtual ~AirProperties(void) {
      NO_OP;
   };

   virtual inline void* pGet(void) const { 
      return (void*)this;
   };
   
   /* Scrive il contributo dell'elemento al file di restart */
   virtual ostream& Restart(ostream& out) const {
      return out << "  air properties: " 
	<< dAirDensity << ", " << dSoundSpeed << ", ",
	pGetDriveCaller()->Restart(out) << ';' << endl;	   
   };
   
   /* Tipo dell'elemento (usato per debug ecc.) */
   virtual ElemType::Type GetElemType(void) const { 
      return ElemType::AIRPROPERTIES; 
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
   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 0;
      *piNumCols = 0;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     AssJac(VariableSubMatrixHandler& WorkMat,
	    doublereal /* dCoef */ , 
	    const VectorHandler& /* XCurr */ ,
	    const VectorHandler& /* XPrimeCurr */ ) {
	WorkMat.SetNullMatrix();
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
				    doublereal /* dCoef */ ,
				    const VectorHandler& /* XCurr */ , 
				    const VectorHandler& /* XPrimeCurr */ ) {
      WorkVec.Resize(0);
      
      /* Approfitto del fatto che AirProperties viene aggiornato prima 
       * degli altri elementi (vedi l'enum ElemType e la sequenza di
       * assemblaggio) per fargli calcolare Velocity una volta per tutte.
       * Quindi, quando viene chiamata GetVelocity(void), 
       * questa restituisce un reference all'accelerazione con il
       * minimo overhead */
      Velocity = Get();
      
      return WorkVec;
   };
   
   
   /* Numero di GDL iniziali */
   virtual unsigned int iGetInitialNumDof(void) const { 
      return 0;
   };
     
   /* Dimensioni initiali del workspace */
   virtual void InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
      *piNumRows = 0;
      *piNumCols = 0;
   };
   
   /* assemblaggio jacobiano */
   virtual VariableSubMatrixHandler& 
     InitialAssJac(VariableSubMatrixHandler& WorkMat,	  
	    const VectorHandler& /* XCurr */ ) {
	WorkMat.SetNullMatrix();
	return WorkMat;
     };
   
   /* assemblaggio residuo */
   virtual SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec,
					   const VectorHandler& XCurr) {
      /* Chiama AssRes, che si limita ad aggiornare la velocita' */
      return AssRes(WorkVec, 1, XCurr, XCurr);
   };
     
     
     
     
    
   /* output; si assume che ogni tipo di elemento sappia, attraverso
    * l'OutputHandler, dove scrivere il proprio output */
   virtual void Output(OutputHandler& /* OH */ ) const { 
      NO_OP; 
   };
   
   virtual const Vec3& GetVelocity(const Vec3& /* X */ ) const { 
      return Velocity;
   };
   
   virtual doublereal dGetAirDensity(void) const {
      return dAirDensity;
   };
   
   virtual doublereal dGetSoundSpeed(void) const { 
      return dSoundSpeed;
   };

#ifdef DEBUG
   virtual const char* sClassName(void) const { return "AirProperties"; };
#endif   
};

/* AirProperties - end */


/* AirPropOwner - begin */

class AirPropOwner {
 protected:
   const AirProperties* pAirProperties;
 public:
   AirPropOwner(void) : pAirProperties(NULL) { NO_OP; };
   virtual ~AirPropOwner(void) { NO_OP; };
   
   virtual void PutAirProperties(const AirProperties* pAP) {
      ASSERT(pAirProperties == NULL);
      (AirProperties*&)pAirProperties = (AirProperties*)pAP;
   };
   
   virtual flag fGetAirVelocity(Vec3& Velocity, const Vec3& X) const {
      if (pAirProperties == NULL) {
	 return 0;
      }
      
      Velocity = pAirProperties->GetVelocity(X);
      return 1;
   };
   
   virtual doublereal dGetSoundSpeed(void) const { 
      return pAirProperties->dGetSoundSpeed(); 
   };
   
   virtual doublereal dGetAirDensity(void) const { 
      return pAirProperties->dGetAirDensity();
   };
};

/* AirPropOwner - end */


/* AerodynamicElem - begin */

class AerodynamicElem : virtual public Elem, public AirPropOwner {
 private:
   AeroType::Type AeroT;
   
 protected:
   
 public:
   AerodynamicElem(unsigned int uL, AeroType::Type T, flag fOut)
     : Elem(uL, ElemType::AERODYNAMIC, fOut), AeroT(T) { 
	NO_OP; 
     };
   
   virtual ~AerodynamicElem(void) { 
      NO_OP; 
   };
   
   /* Tipo di elemento aerodinamico */
   virtual AeroType::Type GetAeroType(void) const { 
      return AeroT;
   };


   /* Consente di effettuare un casting sicuro da Elem* a AerodynamicElem* */
   virtual AerodynamicElem* pGetAerodynamicElem(void) const {
      return (AerodynamicElem*)this; 
   };
   
#ifdef DEBUG
   virtual flag fIsAerodynamicElem(void) const { 
      return flag(1);
   };
#endif   
};

/* AerodynamicElem - end */

#endif
