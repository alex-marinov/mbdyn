/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2004
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

/* gestore dell'output */


#ifndef OUTPUT_H
#define OUTPUT_H

/* se #define DEBUG_COUT l'output avviene su cout anziche' nei files */

#include <ac/iostream>
#include <ac/iomanip>
#include <ac/fstream>

#include <myassert.h>
#include <except.h>
#include <solman.h>
#include <filename.h>

/* OutputHandler - begin */

class OutputHandler : public FileName {
 public:
   enum OutFiles {
      UNKNOWN = -1,
      OUTPUT      = 0,
      STRNODES,
      ELECTRIC,
      ABSTRACT,      
      INERTIA,
      JOINTS,
      FORCES, 
      BEAMS,
      ROTORS,
      RESTART,
      AERODYNAMIC,	/* 10 */
      HYDRAULIC,
      PRESNODES,
      LOADABLE,
      GENELS,
      PARTITION,
      ADAMSRES,
      ADAMSCMD,
      AEROMODALS,
      REFERENCEFRAMES,
      LOG,		/* 20 */
      AIRPROPS,
      PARAMETERS,
      EXTERNALS,

      LASTFILE		/* 24 */
   };   
   
 private:
   
   /* Aggiungere qui i files che si desidera avere a disposizione */
   struct {
       std::ofstream* pof;
      bool UseDefaultPrecision;
      bool UseScientific;
      bool IsOpen;
   } OutData[LASTFILE];
   
   std::ofstream ofOutput;      /*  0 */
   std::ofstream ofStrNodes;
   std::ofstream ofElectric;
   std::ofstream ofAbstract;
   std::ofstream ofInertia;
   std::ofstream ofJoints;      /*  5 */
   std::ofstream ofForces;
   std::ofstream ofBeams;
   std::ofstream ofRotors;
   std::ofstream ofRestart;
   std::ofstream ofAerodynamic; /* 10 */
   std::ofstream ofHydraulic;
   std::ofstream ofPresNodes;
   std::ofstream ofLoadable;
   std::ofstream ofGenels;
   std::ofstream ofPartition;   /* 15 */
   std::ofstream ofAdamsRes;
   std::ofstream ofAdamsCmd;
   std::ofstream ofAeroModals;
   std::ofstream ofReferenceFrames;
   std::ofstream ofLog;
   std::ofstream ofAirProps;
   std::ofstream ofParameters;
   std::ofstream ofExternals;
   
   int iCurrWidth;
   int iCurrPrecision;
   
 public:
   OutputHandler(void);
   
   OutputHandler(const char* sFName, int iExtNum = -1);
   
   void Init(const char* sFName, int iExtNum = -1);
      
   ~OutputHandler(void);
   
   
   /* Aggiungere qui le funzioni che aprono i singoli stream */
   bool Open(const OutputHandler::OutFiles out);

   bool IsOpen(const OutputHandler::OutFiles out) const;
   
   bool Close(const OutputHandler::OutFiles out);
   
   bool OutputOpen(void);
   bool RestartOpen(void);

   bool PartitionOpen(void);
   bool AdamsResOpen(void);
   bool AdamsCmdOpen(void);
   bool LogOpen(void);

   /* Aggiungere qui le funzioni che ritornano gli stream desiderati */
   inline std::ostream& Get(const OutputHandler::OutFiles f) {
      ASSERT(f > -1 && f < LASTFILE);
      ASSERT(OutData[f].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(OutData[f].pof->is_open());
#endif /* HAVE_ISOPEN */
      return *(OutData[f].pof);
   };

   inline std::ostream& Output(void) const {
#ifdef DEBUG_COUT
      return (std::ostream&)cout;
#else
      ASSERT(OutData[OUTPUT].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofOutput.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofOutput;
#endif	
   };   
   
   inline std::ostream& StrNodes(void) const {
      ASSERT(OutData[STRNODES].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofStrNodes.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofStrNodes;	  
   };   
   
   inline std::ostream& Electric(void) const {
      ASSERT(OutData[ELECTRIC].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofElectric.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofElectric;	  
   };
   
   inline std::ostream& Abstract(void) const {
      ASSERT(OutData[ABSTRACT].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAbstract.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAbstract;
   };   
   
   inline std::ostream& Inertia(void) const {
      ASSERT(OutData[INERTIA].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofInertia.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofInertia;	  
   };   
   
   inline std::ostream& Joints(void) const {
      ASSERT(OutData[JOINTS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofJoints.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofJoints;	  
   };   
   
   inline std::ostream& Forces(void) const {
      ASSERT(OutData[FORCES].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofForces.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofForces;	  
   };   
   
   inline std::ostream& Beams(void) const {
      ASSERT(OutData[BEAMS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofBeams.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofBeams;	  
   };   
   
   inline std::ostream& Rotors(void) const {
      ASSERT(OutData[ROTORS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofRotors.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofRotors;	  
   };   
   
   inline std::ostream& Restart(void) const {
      ASSERT(OutData[RESTART].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofRestart.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofRestart;	  
   };   
   
   inline std::ostream& Aerodynamic(void) const {
      ASSERT(OutData[AERODYNAMIC].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAerodynamic.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAerodynamic;	  
   };   
   
   inline std::ostream& Hydraulic(void) const {
      ASSERT(OutData[HYDRAULIC].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofHydraulic.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofHydraulic;	  
   };   

   inline std::ostream& PresNodes(void) const {
      ASSERT(OutData[PRESNODES].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofPresNodes.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofPresNodes;	  
   };   

   inline std::ostream& Loadable(void) const {
      ASSERT(OutData[LOADABLE].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofLoadable.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofLoadable;	  
   };   

   inline std::ostream& Genels(void) const {
      ASSERT(OutData[GENELS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofGenels.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofGenels;	  
   };   

   inline std::ostream& Partition(void) const {
      ASSERT(OutData[PARTITION].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofPartition.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofPartition;	  
   };
   
   inline std::ostream& AdamsRes(void) const {
      ASSERT(OutData[ADAMSRES].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAdamsRes.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAdamsRes;	  
   };

   inline std::ostream& AdamsCmd(void) const {
      ASSERT(OutData[ADAMSCMD].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAdamsCmd.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAdamsCmd;
   };

   inline std::ostream& AeroModals(void) const {
      ASSERT(OutData[AEROMODALS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAeroModals.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAeroModals;	  
   };   

   inline std::ostream& ReferenceFrames(void) const {
      ASSERT(OutData[REFERENCEFRAMES].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofReferenceFrames.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofReferenceFrames;	  
   };   

   inline std::ostream& Log(void) const {
#ifdef DEBUG_COUT
      return (std::ostream&)cout;
#else
      ASSERT(OutData[LOG].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofLog.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofLog;
#endif	
   };   
   
   inline std::ostream& AirProps(void) const {
      ASSERT(OutData[AIRPROPS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofAirProps.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofAirProps;	  
   };

   inline std::ostream& Parameters(void) const {
      ASSERT(OutData[PARAMETERS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofParameters.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofParameters;	  
   };

   inline std::ostream& Externals(void) const {
      ASSERT(OutData[EXTERNALS].IsOpen);
#ifdef HAVE_ISOPEN
      ASSERT(ofExternal.is_open());
#endif /* HAVE_ISOPEN */
      return (std::ostream&)ofExternals;	  
   };

   inline int iW(void) const { 
      return iCurrWidth; 
   };
   
   inline int iP(void) const { 
      return iCurrPrecision; 
   };
   
   void SetWidth(int iNewWidth);
   
   void SetPrecision(int iNewPrecision);
};

extern OutputHandler OutHdl;

/* OutputHandler - end */


/* ToBeOutput - begin */

const flag fDefaultOut = 1;

class ToBeOutput {
 protected:
   flag fOutput;
   
 public:
   ToBeOutput(flag fOut = fDefaultOut);
   virtual ~ToBeOutput(void);
   
   /* Regular output */
   virtual void Output(OutputHandler& OH) const;

   /* Output of perturbed solution (modes ...) */
   virtual void Output(OutputHandler& OH,
		   const VectorHandler& X, const VectorHandler& XP) const;

   /* Output of modes in NASTRAN's pch/f06 format */
#define __NASTRAN_FORMAT_FIXED__	1
#define __NASTRAN_FORMAT_FIXED16__	2
#define __NASTRAN_FORMAT_FREE__		3
   virtual void Output_pch(std::ostream &pch) const;
   virtual void Output_f06(std::ostream &f06, const VectorHandler& X) const;
   virtual void Output_f06(std::ostream &f06, 
		   const VectorHandler& Xr, const VectorHandler& Xi) const;

   /* virtual void AdamsOutput(void) const; */
      
   virtual flag fToBeOutput(void) const;
   virtual void SetOutputFlag(flag f = flag(1));
};

/* ToBeOutput - end */

#endif /* OUTPUT_H */

