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

/* gestore dell'output */


#ifndef OUTPUT_H
#define OUTPUT_H

/* se #define DEBUG_COUT l'output avviene su cout anziche' nei files */

#include <myassert.h>
#include <mystddef.h>
#include <except.h>
#include <solman.h>


/* include di servizio */
#include "filename.h"

/* include generali */
#include <iostream.h>
#include <fstream.h>
#include <iomanip.h>

/* OutputHandler - begin */

class OutputHandler : public FileName {
 public:
   enum OutFiles {
      UNKNOWN = -1,
      OUTPUT      = 0,
      STRNODES    = 1,
      ELECTRIC    = 2,
      ABSTRACT    = 3,      
      INERTIA     = 4,
      JOINTS      = 5,
      FORCES      = 6, 
      BEAMS       = 7,
      ROTORS      = 8,
      RESTART     = 9,
      AERODYNAMIC = 10,
      HYDRAULIC   = 11,
      PRESNODES   = 12,
      LOADABLE    = 13,
      GENELS      = 14,
      PARTITION   = 15,
      ADAMSRES    = 16,
      ADAMSCMD    = 17,
      
      LASTFILE    /* = 18 */
   };   
   
 private:
   
   /* Aggiungere qui i files che si desidera avere a disposizione */
   struct {
      ofstream* pof;
      flag fToUseDefaultPrecision;
      flag fToUseScientific;
   } OutData[LASTFILE];
   
   ofstream ofOutput;      /*  0 */
   ofstream ofStrNodes;
   ofstream ofElectric;
   ofstream ofAbstract;
   ofstream ofInertia;
   ofstream ofJoints;      /*  5 */
   ofstream ofForces;
   ofstream ofBeams;
   ofstream ofRotors;
   ofstream ofRestart;
   ofstream ofAerodynamic; /* 10 */
   ofstream ofHydraulic;
   ofstream ofPresNodes;
   ofstream ofLoadable;
   ofstream ofGenels;
   ofstream ofPartition;   /* 15 */
   ofstream ofAdamsRes;
   ofstream ofAdamsCmd;
   
   int iCurrWidth;
   int iCurrPrecision;
   
 public:
   OutputHandler(void);
   
   OutputHandler(const char* sFName, int iExtNum = -1);
   
   void Init(const char* sFName, int iExtNum = -1);
      
   ~OutputHandler(void);
   
   
   /* Aggiungere qui le funzioni che aprono i singoli stream */
   int Open(OutFiles out);
   
   
   int OutputOpen(void);

#if 0
   int StrNodesOpen(void);
   int ElectricOpen(void);
   int AbstractOpen(void);
   int InertiaOpen(void);
   int JointsOpen(void);
   int ForcesOpen(void);
   int BeamsOpen(void);
   int RotorsOpen(void);
#endif /* 0 */
   
   int RestartOpen(void);

#if 0
   int AerodynamicOpen(void);
   int HydraulicOpen(void);
   int PresNodesOpen(void);
   int LoadableOpen(void);
   int GenelsOpen(void);
#endif /* 0 */

   int PartitionOpen(void);
   int AdamsResOpen(void);
   int AdamsCmdOpen(void);

   /* Aggiungere qui le funzioni che ritornano gli stream desiderati */
   inline ostream& Get(const OutputHandler::OutFiles f) {
      ASSERT(f > -1 && f < LASTFILE);
#if HAVE_ISOPEN
      ASSERT(OutData[f].pof->is_open());
#endif /* HAVE_ISOPEN */
      return *(OutData[f].pof);
   };

   inline ostream& Output(void) const {
#ifdef DEBUG_COUT
      return (ostream&)cout;
#else
#if HAVE_ISOPEN
      ASSERT(ofOutput.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofOutput;
#endif	
   };   
   
   inline ostream& StrNodes(void) const {
#if HAVE_ISOPEN
      ASSERT(ofStrNodes.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofStrNodes;	  
   };   
   
   inline ostream& Electric(void) const {
#if HAVE_ISOPEN
      ASSERT(ofElectric.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofElectric;	  
   };
   
   inline ostream& Abstract(void) const {
#if HAVE_ISOPEN
      ASSERT(ofAbstract.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofAbstract;
   };   
   
   inline ostream& Inertia(void) const {
#if HAVE_ISOPEN
      ASSERT(ofInertia.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofInertia;	  
   };   
   
   inline ostream& Joints(void) const {
#if HAVE_ISOPEN
      ASSERT(ofJoints.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofJoints;	  
   };   
   
   inline ostream& Forces(void) const {
#if HAVE_ISOPEN
      ASSERT(ofForces.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofForces;	  
   };   
   
   inline ostream& Beams(void) const {
#if HAVE_ISOPEN
      ASSERT(ofBeams.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofBeams;	  
   };   
   
   inline ostream& Rotors(void) const {
#if HAVE_ISOPEN
      ASSERT(ofRotors.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofRotors;	  
   };   
   
   inline ostream& Restart(void) const {
#if HAVE_ISOPEN
      ASSERT(ofRestart.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofRestart;	  
   };   
   
   inline ostream& Aerodynamic(void) const {
#if HAVE_ISOPEN
      ASSERT(ofAerodynamic.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofAerodynamic;	  
   };   
   
   inline ostream& Hydraulic(void) const {
#if HAVE_ISOPEN
      ASSERT(ofHydraulic.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofHydraulic;	  
   };   

   inline ostream& PresNodes(void) const {
#if HAVE_ISOPEN
      ASSERT(ofPresNodes.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofPresNodes;	  
   };   

   inline ostream& Loadable(void) const {
#if HAVE_ISOPEN
      ASSERT(ofLoadable.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofLoadable;	  
   };   

   inline ostream& Genels(void) const {
#if HAVE_ISOPEN
      ASSERT(ofGenels.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofGenels;	  
   };   

   inline ostream& Partition(void) const {
#if HAVE_ISOPEN
      ASSERT(ofPartition.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofPartition;	  
   };
   
   inline ostream& AdamsRes(void) const {
#if HAVE_ISOPEN
      ASSERT(ofAdamsRes.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofAdamsRes;	  
   };

   inline ostream& AdamsCmd(void) const {
#if HAVE_ISOPEN
      ASSERT(ofAdamsCmd.is_open());
#endif /* HAVE_ISOPEN */
      return (ostream&)ofAdamsCmd;
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
   ToBeOutput(flag fOut = fDefaultOut) : fOutput(fOut) { NO_OP; };
   virtual ~ToBeOutput(void) { NO_OP; };
   
   /* Regular output */
   virtual void Output(OutputHandler& OH) const = 0;

   /* Output of perturbed solution (modes ...) */
   virtual void Output(OutputHandler& OH,
		   const VectorHandler& X, const VectorHandler& XP) const {
	   NO_OP;
   };

   /* Output of modes in NASTRAN's pch/f06 format */
   virtual void Output_pch(ostream &pch) const {
	   NO_OP;
   };
   virtual void Output_f06(ostream &f06, const VectorHandler& X) const {
	   NO_OP;
   };

   /* virtual void AdamsOutput(void) const; */
      
   virtual flag fToBeOutput(void) const {
      return fOutput;
   };
   
   virtual void SetOutputFlag(flag f = flag(1)) {
      fOutput = f;
   };
};

/* ToBeOutput - end */

#endif
