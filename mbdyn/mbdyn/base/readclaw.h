/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2003
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

#ifndef READCLAW_H
#define READCLAW_H

#include <drive_.h>		/* per TimeDriveCaller */
#include <tpldrive.h>
#include <tpldrive_.h>
#include <constltp_.h>
#include <symcltp.h>
#ifdef USE_GRAALLDAMPER
#include <damper.h>
#endif /* USE_GRAALLDAMPER */
#include <shockabsorber.h>


template <class T>
void GetPreStress(MBDynParser& HP, T& PreStress)
{
   if (HP.IsKeyWord("prestress")) {
      PreStress = HP.Get(PreStress);
   }
}


template <class T>
TplDriveCaller<T>* GetPreStrain(DataManager* pDM, 
				MBDynParser& HP,
				DriveHandler* pDH,
				T& PreStrain)
{
   if (HP.IsKeyWord("prestrain")) {
      return ReadTplDrive(pDM, HP, pDH, PreStrain);
   } else {
      DriveCaller* pDC = NULL;
      SAFENEWWITHCONSTRUCTOR(pDC,
			     NullDriveCaller,
			     NullDriveCaller(pDH));

      T t(0.);
      TplDriveCaller<T>* pTplDC = NULL;
      SAFENEWWITHCONSTRUCTOR(pTplDC,
			     SingleTplDriveCaller<T>,
			     SingleTplDriveCaller<T>(pDC, t));
      return pTplDC;
   }
}


template <class T, class Tder>
ConstitutiveLaw<T, Tder>* ReadConstLaw(DataManager* pDM,
				       MBDynParser& HP,
				       DriveHandler* pDH,
				       ConstLawType::Type& CLType,
				       ConstitutiveLaw<T, Tder>*)
{
   DEBUGCOUT("Entering ReadConstLaw" << std::endl);
   
   const char* sKeyWords[] = {
      "linear" "elastic",
      "linear" "elastic" "isotropic",
	"linear" "elastic" "generic",
	"linear" "elastic" "generic" "axial" "torsion" "coupling",
	"linear" "elastic" "bistop",
	"log" "elastic",
	"double" "linear" "elastic",
	"isotropic" "hardening" "elastic",
	"contact" "elastic",
	"symbolic" "elastic" "isotropic",
	"linear" "viscous",
	"linear" "viscous" "isotropic",
	"linear" "viscous" "generic",
	"linear" "viscoelastic",
	"linear" "viscoelastic" "isotropic",
	"linear" "viscoelastic" "generic",
	"doublelinear" "viscoelastic",
	"turbulent" "viscoelastic",
	"linear" "viscoelastic" "bistop",
	"graall" "damper",
	"shock" "absorber",
	NULL
   };
   
   /* enum delle parole chiave */
   enum KeyWords {
      UNKNOWN = -1,
	
	LINEARELASTIC = 0,
	LINEARELASTICISOTROPIC,
	LINEARELASTICGENERIC,
	LINEARELASTICGENERICAXIALTORSIONCOUPLING,
	LINEARELASTICBISTOP,
	LOGELASTIC,
	DOUBLELINEARELASTIC,
	ISOTROPICHARDENINGELASTIC,
	CONTACTELASTIC,
	SYMBOLICELASTICISOTROPIC,
	LINEARVISCOUS,
	LINEARVISCOUSISOTROPIC,
	LINEARVISCOUSGENERIC,
	LINEARVISCOELASTIC,
	LINEARVISCOELASTICISOTROPIC,
	LINEARVISCOELASTICGENERIC,
	DOUBLELINEARVISCOELASTIC,
	TURBULENTVISCOELASTIC,
	LINEARVISCOELASTICBISTOP,
	GRAALLDAMPER,
	SHOCKABSORBER,
	
	LASTKEYWORD
   };
   int CurrKW;
   
   /* tabella delle parole chiave */
   KeyTable K(HP, sKeyWords);
   
   ConstitutiveLaw<T, Tder>* pCL = NULL;   
   CurrKW = HP.GetWord();
   switch (CurrKW) {
    case LINEARELASTIC:
    case LINEARELASTICISOTROPIC: {
       CLType = ConstLawType::ELASTIC;
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("Linear Elastic Isotropic Constitutive Law, stiffness = "
		 << dS << std::endl);
       
       if (dS <= 0.) {		      
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }

       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);

       typedef LinearElasticIsotropicConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS));
       
       break;
    }

    case LINEARELASTICGENERIC: {
       CLType = ConstLawType::ELASTIC;
             
       DEBUGCOUT("Linear Elastic Generic Constitutive Law" << std::endl);
       Tder S(0.);
       S = HP.Get(S);
              
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef LinearElasticGenericConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S));
       
       break;
    }

    case LINEARELASTICGENERICAXIALTORSIONCOUPLING: {
       CLType = ConstLawType::ELASTIC;
             
       DEBUGCOUT("Linear Elastic Generic Constitutive Law with Axial-Torsion Coupling" << std::endl);
       Tder S(0.);
       S = HP.Get(S);
       
       // coefficiente di accoppiamento
       doublereal dCoupl = HP.GetReal();
       DEBUGCOUT("coupling coefficient: " << dCoupl << std::endl);
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef LinearElasticGenericAxialTorsionCouplingConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S, dCoupl));
       
       break;
    }

   case LOGELASTIC: {
       CLType = ConstLawType::ELASTIC;
             
       DEBUGCOUT("Logaritmic Elastic Constitutive Law" << std::endl);
       doublereal dS = HP.GetReal();
       if (dS <= 0.) {		      
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef LogConstitutiveLaw<T, Tder> L;
       // typedef LogConstitutiveLaw L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS));
       
       break;
    }      
      
    case DOUBLELINEARELASTIC: {
       CLType = ConstLawType::ELASTIC;
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("stiffness = " << dS << std::endl);
       
       if (dS <= 0.) {
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dUpp = HP.GetReal();
       if (dUpp <= 0.) {
	  std::cerr << "warning, null or negative upper limit strain at line "
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dLow = HP.GetReal();
       if (dLow >= 0.) {
	  std::cerr << "warning, null or positive lower limit strain at line "
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dSecondS = HP.GetReal();
       if (dSecondS <= 0.) {		      
	  std::cerr << "warning, null or negative second stiffness at line "
	    << HP.GetLineData() << std::endl;
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef DoubleLinearElasticConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, 
			      L, 
			      L(pTplDC, PreStress, dS, dUpp, dLow, dSecondS));
       
       break;
    }
      
    case ISOTROPICHARDENINGELASTIC: {
       CLType = ConstLawType::ELASTIC;
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("Stiffness = " << dS << std::endl);
       
       if (dS <= 0.) {		      
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dE = HP.GetReal();
       DEBUGCOUT("Reference strain = " << dE << std::endl);
       
       if (dE <= 0.) {		      
	  std::cerr << "error, null or negative reference strain at line " 
	    << HP.GetLineData() << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef IsotropicHardeningConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dE));

       break;
    }
      
    case CONTACTELASTIC: {
       CLType = ConstLawType::ELASTIC;
       
       doublereal dK = HP.GetReal();
       DEBUGCOUT("Stiffness = " << dK << std::endl);
       
       if (dK <= 0.) {		      
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dGamma = HP.GetReal();
       DEBUGCOUT("Exponent = " << dGamma << std::endl);
       
       if (dGamma < 1.) {
	  std::cerr << "error, exponent < 1. at line "
	    << HP.GetLineData() << std::endl;	  
	  THROW(DataManager::ErrGeneric());
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef ContactConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dK, dGamma));

       break;
    }

    case SYMBOLICELASTICISOTROPIC: {
       CLType = ConstLawType::ELASTIC;

       const char *epsilon = 0;
       if (HP.IsKeyWord("epsilon")) {
	       const char *tmp = HP.GetStringWithDelims();
       	       if (tmp == 0) {
		       std::cerr << "unable to get \"epsilon\" symbol at line " << HP.GetLineData() << std::endl;
		       THROW(DataManager::ErrGeneric());
	       }
	       SAFESTRDUP(epsilon, tmp);
       } else {
	       epsilon = "epsilon";
       }

       if (!HP.IsKeyWord("expression")) {
	       std::cerr << "keyword \"expression\" expected at line " << HP.GetLineData() << std::endl;
	       THROW(DataManager::ErrGeneric());
       }
       const char *tmp = HP.GetStringWithDelims();
       if (tmp == 0) {
	       std::cerr << "unable to get \"expression\" at line " << HP.GetLineData() << std::endl;
	       THROW(DataManager::ErrGeneric());
       }
       const char *expression = 0;
       SAFESTRDUP(expression, tmp);

       const char **symbols = 0;
       if (HP.IsKeyWord("symbols")) {
	       /* FIXME! */
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef SymbolicElasticIsotropicConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, epsilon, expression, symbols));

       break;
    }
      
    case LINEARVISCOUS:
    case LINEARVISCOUSISOTROPIC: {
       CLType = ConstLawType::VISCOUS;
       
       doublereal dSP = HP.GetReal();
       DEBUGCOUT("stiffness prime = " << dSP << std::endl);
       
       if (dSP <= 0.) {
	  std::cerr << "warning, null or negative stiffness prime at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       /* Prestress (no prestrain) */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       
       typedef LinearViscousIsotropicConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(NULL, PreStress, dSP));
       
       break;
    }
      
    case LINEARVISCOUSGENERIC: {
       CLType = ConstLawType::VISCOUS;
       
       Tder SP(0.);
       SP = HP.Get(SP);
       
       /* Prestress (no prestrain) */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       
       typedef LinearViscousGenericConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(NULL, PreStress, SP));
       
       break;
    }
      
    case LINEARVISCOELASTIC:
    case LINEARVISCOELASTICISOTROPIC: {
       CLType = ConstLawType::VISCOELASTIC;     
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("Stiffness = " << dS << std::endl);
       
       if (dS <= 0.) {
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dSP = 0.;
       if (HP.IsKeyWord("proportional")) {
	  doublereal k = HP.GetReal();
	  dSP = k*dS;
       } else {	 
	  dSP = HP.GetReal();
       }       
       DEBUGCOUT("stiffness prime = " << dSP << std::endl);
       
       if (dSP <= 0.) {
	  std::cerr << "warning, null or negative stiffness prime at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef LinearViscoElasticIsotropicConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dSP));
       
       break;
    }
      
    case LINEARVISCOELASTICGENERIC: {
       CLType = ConstLawType::VISCOELASTIC;     
       
       Tder S(0.);
       S = HP.Get(S);      

       Tder SP(0.);
       if (HP.IsKeyWord("proportional")) {
	  doublereal k = HP.GetReal();
	  SP = S*k;
       } else {	 
	  SP = HP.Get(SP);
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef LinearViscoElasticGenericConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, S, SP));
       
       break;
    }
      
    case DOUBLELINEARVISCOELASTIC: {
       CLType = ConstLawType::VISCOELASTIC;
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("stiffness = " << dS << std::endl);
       
       if (dS <= 0.) {
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dUpp = HP.GetReal();
       if (dUpp <= 0.) {
	  std::cerr << "warning, null or negative upper limit strain at line "
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dLow = HP.GetReal();
       if (dLow >= 0.) {
	  std::cerr << "warning, null or positive lower limit strain at line "
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dSecondS = HP.GetReal();
       if (dSecondS <= 0.) {		      
	  std::cerr << "warning, null or negative second stiffness at line "
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dSP = HP.GetReal();
       DEBUGCOUT("stiffness prime = " << dSP << std::endl);
       
       if (dSP <= 0.) {
	  std::cerr << "warning, null or negative stiffness prime at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       typedef DoubleLinearViscoElasticConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, 
			      L, 
			      L(pTplDC, PreStress, 
				dS, dUpp, dLow, dSecondS, dSP));
       
       break;
    }
      
    case TURBULENTVISCOELASTIC:	{
       CLType = ConstLawType::VISCOELASTIC;
       
       doublereal dS = HP.GetReal();
       DEBUGCOUT("Visco-Elastic Turbulent Rod Joint, stiffness = " 
		 << dS << std::endl);
       
       if (dS <= 0.) {
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dParabStiff = HP.GetReal();
       DEBUGCOUT("stiffness prime = " << dParabStiff << std::endl);
       
       if (dParabStiff <= 0.) {
	  std::cerr << "warning, null or negative derivative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }
       
       doublereal dTreshold = 0.;
       if (HP.IsArg()) {
	  dTreshold = HP.GetReal(dTreshold);
	  
	  // Il legame costitutivo ha la forma seguente:
	  //    F = Kp*e + Kd*(de/dt)
	  // con Kp costante e Kd dato dalla seguente legge:
	  //    Kd = cost2                per fabs(de/dt) < Treshold
	  //    Kd = 2*cost1*fabs(de/dt)  per fabs(de/dt) > Treshold
	  // se non viene inserito il valore di treshold, lo si
	  // assume = 0. e quindi il legame e' sempre del secondo
	  // tipo. Altrimenti, se non viene inserita la seconda
	  // costante cost2, si assume che vi sia raccordo tra 
	  // i due tipi di legge, ovvero cost2 = cost1*Treshold
	  // altrimenti e' possibile avere un comportamento,
	  // che in prima approssimazione e' valido 
	  // per numerosi fluidi, in cui vi e' un salto tra i due 
	  // tipi di legge costitutiva.	 
       }
       
       doublereal dSP = dTreshold*dParabStiff;		      
       if (HP.IsArg()) {			 
	  dSP = HP.GetReal(dSP);
       }
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
	
       typedef TurbulentViscoElasticConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, 
			      L, 
			      L(pTplDC, PreStress, 
				dS, dSP, dTreshold, dParabStiff));
       
       break;
    }	      	     	      
        
				
    case LINEARELASTICBISTOP:
    case LINEARVISCOELASTICBISTOP: {
       typedef LinearViscoElasticBiStopConstitutiveLaw<T, Tder> L;
       CLType = ConstLawType::VISCOELASTIC;
             
       DEBUGCOUT("Linear Viscoelastic Bi Stop Constitutive Law" << std::endl);
       doublereal dS = HP.GetReal();
       if (dS <= 0.) {
	  std::cerr << "warning, null or negative stiffness at line " 
	    << HP.GetLineData() << std::endl;
       }

       doublereal dSp = 0.;
       if (CurrKW == LINEARVISCOELASTICBISTOP) {
	  dSp = HP.GetReal();
	  if (dSp <= 0.) {
	     std::cerr << "warning, null or negative stiffness prime at line " 
		     << HP.GetLineData() << std::endl;
	  }
       }

       typedef typename LinearViscoElasticBiStopConstitutiveLaw<T, Tder>::Status LS;
       LS s = L::INACTIVE;
       if (HP.IsKeyWord("initialstatus")) {
          if (HP.IsKeyWord("active")) {
	     s = L::ACTIVE;
	  } else if (HP.IsKeyWord("inactive")) {
	     s = L::INACTIVE;
	  } else {
	     std::cerr << "unknown initial status at line " 
		     << HP.GetLineData() << std::endl;
	     THROW(ErrGeneric());
	  }
       }

       const DriveCaller *pA = ReadDriveData(pDM, HP, pDH);
       const DriveCaller *pD = ReadDriveData(pDM, HP, pDH);
       
       /* Prestress and prestrain */
       T PreStress(0.);
       GetPreStress(HP, PreStress);
       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
       
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pTplDC, PreStress, dS, dSp, s, pA, pD));
       
       break;
    }      
	   
    case GRAALLDAMPER: {
#ifdef USE_GRAALLDAMPER    
       CLType = ConstLawType::VISCOELASTIC;

       const char* filename = HP.GetFileName();
       DEBUGCOUT("Graall damper input file: \"" 
		       << filename << "\"" << std::endl);
       
       doublereal rla = HP.GetReal();
       DEBUGCOUT("Reference length: " << rla << std::endl);
       
       DriveCaller* pDC = NULL;
       SAFENEWWITHCONSTRUCTOR(pDC,
			      TimeDriveCaller,
			      TimeDriveCaller(pDH));
       
       T t(1.);
       TplDriveCaller<T>* pTplDC = NULL;
       SAFENEWWITHCONSTRUCTOR(pTplDC,
                             SingleTplDriveCaller<T>,
                             SingleTplDriveCaller<T>(pDC, t));
       
       typedef GRAALLDamperConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, 
			      L, 
			      L(pTplDC, rla, filename));
              
       break;
#else /* USE_GRAALLDAMPER */
       std::cerr << "can't use GRAALL Damper" << std::endl;
       THROW(ErrGeneric());
#endif /* USE_GRAALLDAMPER */
    }

    /*
     * Shock absorber per Stefy:
     *
     * ``Riprogettazione dell'ammortizzatore del carrello anteriore
     * di un velivolo di aviazione generale'', 
     * S. Carlucci e S. Gualdi,
     * A.A. 1997-98
     */
    case SHOCKABSORBER: {
       CLType = ConstLawType::VISCOELASTIC;

       T PreStrain(0.);
       TplDriveCaller<T>* pTplDC = GetPreStrain(pDM, HP, pDH, PreStrain);
	      
       typedef ShockAbsorberConstitutiveLaw<T, Tder> L;
       SAFENEWWITHCONSTRUCTOR(pCL, L, L(pDM, pTplDC, HP));

       break;
    }

	   
    // Aggiungere altri rods
      
    default: {
       std::cerr << "Unknown constitutive law type at line " 
	 << HP.GetLineData() << std::endl;
       
       THROW(ErrGeneric());
    }		       
   }

   ASSERT(pCL != NULL);
   return pCL;
} /* End of ReadConstLaw */


/* Legge un legame costitutivo tridimensionale

ConstitutiveLaw3D* DataManager::ReadConstLaw3D(MBDynParser& HP, 
					       ConstLawType::Type& T)
{
   const char sFuncName[] = "DataManager::ReadConstLaw3D()";
   DEBUGCOUT("Entering " << sFuncName << std::endl);
   
   const char* sKeyWords[] = {
      "linear" "elastic" "isotropic",        // Lineare elastico
	"linear" "elastic" "orthotropic",    // Generalmente anisotropa
	"doublelinear" "elastic",            // Spezzata lineare elastica
	"linear" "viscous" "isotropic",      // Lineare viscoso
	"turbulent" "viscous" "isotropic",   // Quadratico viscoso
	"velocity" "damper",                 // smorzatore di velocita' ang.
	"linear" "viscoelastic" "isotropic"  // Lineare viscoelastico
   };
   
   // enum delle parole chiave
   enum KeyWords {
      UNKNOWN = -1,
	LINEARELASTICISOTROPIC = 0,
	LINEARELASTICORTHOTROPIC, 
	DOUBLELINEARELASTIC,
	LINEARVISCOUSISOTROPIC,
	TURBULENTVISCOUSISOTROPIC,
	VELOCITYDAMPER,
	LINEARVISCOELASTICISOTROPIC,
	LASTKEYWORD
   };
   
   // tabella delle parole chiave
   KeyTable K(HP, sKeyWords);
   
   ConstitutiveLaw3D* pCL = NULL;   
   switch(HP.GetWord()) {
      
    case LINEARELASTICISOTROPIC:	
	{	   
	   T = ConstLawType::ELASTIC;
		   
	   doublereal dS = HP.GetReal();
	   DEBUGCOUT("stiffness = " << dS << std::endl);
	   
	   if (dS <= 0.) {		      
	      std::cerr << "warning, null or negative stiffness at line " 
		<< HP.GetLineData() << std::endl;
	   }
	   	   
	   Vec3 PreStress(0.);
	   if(HP.IsArg()) {
	      PreStress = HP.GetVec3();
	   }
	   
	   Vec3 PreStrain(0.);
	   if(HP.IsArg()) {
	      PreStrain = HP.GetVec3();
	   }		   
	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL, 
				  LinearElasticIsotropicConstitutiveLaw3D,
				  LinearElasticIsotropicConstitutiveLaw3D(dS, 
									  PreStress,
									  PreStrain));
	   break;
	}
	      
    case LINEARELASTICORTHOTROPIC:
	{	   
	   T = ConstLawType::ELASTIC;
		   
	   Mat3x3 K(HP.GetMat3x3());
	   
	   Vec3 PreStress(0.);
	   if(HP.IsArg()) {
	      PreStress = HP.GetVec3();
	   }
	   
	   Vec3 PreStrain(0.);
	   if(HP.IsArg()) {
	      PreStrain = HP.GetVec3();
	   }		   
	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL, 
				  LinearElasticOrthotropicConstitutiveLaw3D,
				  LinearElasticOrthotropicConstitutiveLaw3D(K, 
									    PreStress,
									    PreStrain));
	   break;
	}
	      
    case DOUBLELINEARELASTIC:	
	{	   
	   T = ConstLawType::ELASTIC;
		   
	   doublereal dS = HP.GetReal();
	   DEBUGCOUT("stiffness = " << dS << std::endl);
	   
	   if (dS <= 0.) {		      
	      std::cerr << "warning, null or negative stiffness at line " 
		<< HP.GetLineData() << std::endl;
	   }
	   
	   doublereal dUpp = HP.GetReal();
	   if (dUpp <= 0.) {
	      std::cerr << "warning, null or negative upper limit strain at line "
		<< HP.GetLineData() << std::endl;
	   }
	   
	   doublereal dLow = HP.GetReal();
	   if(dLow >= 0.) {
	      std::cerr << "warning, null or positive lower limit strain at line "
		<< HP.GetLineData() << std::endl;
	   }
	   
	   doublereal dSecondS = HP.GetReal();
	   if(dSecondS <= 0.) {		      
	      std::cerr << "warning, null or negative second stiffness at line "
		<< HP.GetLineData() << std::endl;
	   }
	   	   
	   Vec3 PreStress(0.);
	   if(HP.IsArg()) {
	      PreStress = HP.GetVec3();
	   }
	   
	   Vec3 PreStrain(0.);
	   if(HP.IsArg()) {
	      PreStrain = HP.GetVec3();
	   }		   
	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL,
				  DoubleLinearElasticConstitutiveLaw3D,
				  DoubleLinearElasticConstitutiveLaw3D(dS,
								       dUpp, 
								       dLow,
								       dSecondS,
								       PreStress,
								       PreStrain));
	   break;
	}
	      
    case LINEARVISCOUSISOTROPIC:
	{	   
	   T = ConstLawType::VISCOUS;
		   
	   doublereal dS = HP.GetReal();
	   DEBUGCOUT("stiffness prime = " << dS << std::endl);
	   
	   if (dS <= 0.) {		      
	      std::cerr << "warning, null or negative stiffness prime at line "
		<< HP.GetLineData() << std::endl;
	   }
	   	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL,
				  LinearViscousIsotropicConstitutiveLaw3D,
				  LinearViscousIsotropicConstitutiveLaw3D(dS));
	   break;
	}
	 
    case TURBULENTVISCOUSISOTROPIC:           
	{
	   T = ConstLawType::VISCOUS;	   
	   
	   doublereal dParabStiff = HP.GetReal();
	   DEBUGCOUT("Visco-Elastic Turbulent Rod Joint, stiffness = " 
		     << dParabStiff << std::endl);
	   
	   if (dParabStiff <= 0.) {
	      std::cerr << "warning, null or negative stiffness at line " 
		<< HP.GetLineData() << std::endl;
	   }
	   
	   doublereal dTreshold = 0.;
	   if(HP.IsArg()) {
	      dTreshold = HP.GetReal(dTreshold);
	      
	      // Il legame costitutivo ha la forma seguente:
	      //    F = Kp*e + Kd*(de/dt)
	      // con Kp costante e Kd dato dalla seguente legge:
	      //    Kd = cost2                per fabs(de/dt) < Treshold
	      //    Kd = 2*cost1*fabs(de/dt)  per fabs(de/dt) > Treshold
	      // se non viene inserito il valore di treshold, lo si
	      // assume = 0. e quindi il legame e' sempre del secondo
	      // tipo. Altrimenti, se non viene inserita la seconda
	      // costante cost2, si assume che vi sia raccordo tra 
	      // i due tipi di legge, ovvero cost2 = cost1*Treshold
	      // altrimenti e' possibile avere un comportamento,
	      // che in prima approssimazione e' valido 
	      // per numerosi fluidi, in cui vi e' un salto tra i due 
	      // tipi di legge costitutiva.	     
	   }
	   
	   doublereal dSP = dTreshold*dParabStiff;
	   if(HP.IsArg()) {			 
	      dSP = HP.GetReal(dSP);
	   }
	   	   
	   // Legame costitutivo
	   if(dTreshold > 0.) {
	      SAFENEWWITHCONSTRUCTOR(pCL,
				     LinearTurbulentViscousIsotropicConstitutiveLaw3D,
				     LinearTurbulentViscousIsotropicConstitutiveLaw3D(dSP,
										      dTreshold,
										      dParabStiff));
	   } else {
	      SAFENEWWITHCONSTRUCTOR(pCL, 
				     TurbulentViscousIsotropicConstitutiveLaw3D,
				     TurbulentViscousIsotropicConstitutiveLaw3D(dParabStiff));
	   }	   
	   
	   break;
	}	      	     	            
      
    case VELOCITYDAMPER:
	{
	   T = ConstLawType::VISCOUS;
	
	   doublereal dR = HP.GetReal();	   	   
	   DEBUGCOUT("Reference value: " << dR << std::endl);
	   
	   doublereal dS = HP.GetReal();
	   DEBUGCOUT("stiffness prime = " << dS << std::endl);
	   
	   if (dS <= 0.) {		      
	      std::cerr << "warning, null or negative stiffness prime at line "
		<< HP.GetLineData() << std::endl;
	   }
	   	   
	   doublereal dV = HP.GetReal();
	   DEBUGCOUT("Reference velocity around axis 3: " << dV << std::endl);
	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL,
				  VelocityDamperConstitutiveLaw3D,
				  VelocityDamperConstitutiveLaw3D(dR, dS, dV));
	   break;
	}
            
    case LINEARVISCOELASTICISOTROPIC:
	{
	   T = ConstLawType::VISCOELASTIC;	   
	   
	   doublereal dS = HP.GetReal();
	   DEBUGCOUT("Visco-Elastic Linear Rod Joint, stiffness = " 
		     << dS << std::endl);
	   
	   if (dS <= 0.) {
	      std::cerr << "warning, null or negative stiffness at line " 
		<< HP.GetLineData() << std::endl;
	   }
	   
	   doublereal dSP = HP.GetReal();
	   DEBUGCOUT("stiffness prime = " << dSP << std::endl);
	   
	   if (dSP <= 0.) {
	      std::cerr << "warning, null or negative stiffness prime at line " 
		<< HP.GetLineData() << std::endl;
	   }
	   
	   Vec3 PreStress(0.);
	   if(HP.IsArg()) {
	      PreStress = HP.GetVec3();
	   }
	   
	   Vec3 PreStrain(0.);
	   if(HP.IsArg()) {
	      PreStrain = HP.GetVec3();
	   }		   
	   
	   // Legame costitutivo
	   SAFENEWWITHCONSTRUCTOR(pCL,
				  LinearViscoElasticIsotropicConstitutiveLaw3D,
				  LinearViscoElasticIsotropicConstitutiveLaw3D(dS, dSP,
									       PreStress,
									       PreStrain));
	   
	   break;
	}      
      
      // Aggiungere altri rods
    default:
	{
	   std::cerr << "Unknown constitutive law type at line " 
	     << HP.GetLineData() << std::endl;

	   THROW(DataManager::ErrGeneric());
	}
   }
   
   return pCL;
} End of ReadConstLaw3D */

#endif // READCLAW_H
