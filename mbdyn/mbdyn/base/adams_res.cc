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

/* Continua il DataManager */

#include <mbconfig.h>

#include <dataman.h>
#include <dataman_.h>

flag
DataManager::fAdamsOutput(void) const
{
	return fAdamsResOutput;
}

void 
DataManager::AdamsResOutputInit(void)
{
   OutHdl.AdamsResOpen();
   ostream& out = OutHdl.AdamsRes();
   time_t t = time(NULL);

   unsigned int nStrNodes = NodeData[NodeType::STRUCTURAL].iNum;
   iAdamsOutputParts = nStrNodes;
   
   Elem* p = NULL;
   if (ElemIter.fGetFirst(p)) {
      do {
	 ASSERT(p != NULL);
	 if (p != NULL) {
	    iAdamsOutputParts += p->iGetNumAdamsDummyParts();
	 }
      } while (ElemIter.fGetNext(p));
   }
   
   
   
   /* Cmd Init */
   OutHdl.AdamsCmdOpen();
   ostream& cmd = OutHdl.AdamsCmd();
   
   cmd << sAdamsModelName << endl;

#if 0
   cmd 
     << "!" << endl
     << "model create  &" << endl
     << "   model_name = " << sAdamsModelName << endl
     << "!" << endl
     << "defaults model  &" << endl
     << "   part_name = ground" << endl
     << "!" << endl
     << "defaults coordinate_system  &" << endl
     << "   default_coordinate_system = ." << sAdamsModelName << ".ground" << endl;   
#endif
   
   out
     
     /* Header Block: File Type and Source Sub-Block */
     << "FILE TYPE AND SOURCE                           6" << endl
     << "ADAMS RESULTS FILE       M_KGS_N_RACA 0.100000000000000E+01" << endl
     << "ADAMS VERSION - 10.0                    A.RE7" << endl
     << "Mechanical Dynamics, Inc." << endl
     << "Unknown" << endl
     << ctime(&t)
     << "!" << endl
     
     /* Header Block: Data Set Title Sub-Block */
     << "DATA SET TITLE                                 3" << endl
     << "ADAMS/View model name: " << (sSimulationTitle ? sSimulationTitle : "") << endl
     << "" << endl
     << "!" << endl

     /* Header Block: Results Set Sub-Block */
     << "RESULTS TITLE                                  3" << endl
     << "" << endl
     << "" << endl
     << "!" << endl
     
     /* Header Block: Number of Analysis Output Blocs Sub-Block */
     << "NUMBER OF ANALYSIS BLOCKS                      2" << endl
     << setw(8) << 2+10000 << setw(8) << 123+3*iAdamsOutputParts << endl
     << "!" << endl
     
     /* Rigid Parts Map Block */
     << "RIGID PARTS MAP                         " << setw(8) << 2+iAdamsOutputParts << endl
     << setw(8) << iAdamsOutputParts << "YESNO NO" << endl;

   unsigned int i;
   for (i = 0; i < nStrNodes; i++) {
      out 
	<< setw(20) << 2+i  /* 2 perche' 1 e' il "ground" in Adams! */
	<< setw(8) << 5+i
	<< setw(8) << 5+i+iAdamsOutputParts
	<< setw(8) << 5+i+iAdamsOutputParts << endl;
   }
   
   if (ElemIter.fGetFirst(p)) {
      do {
	 ASSERT(p != NULL);
	 if (p != NULL) {
	    for (unsigned int part = 1; part <= p->iGetNumAdamsDummyParts(); part++) {
	       out
		 << setw(20) << 2+i
		 << setw(8) << 5+i
		 << setw(8) << 5+i+iAdamsOutputParts
		 << setw(8) << 5+i+iAdamsOutputParts << endl;
	       i++;
	    }
	 }
      } while (ElemIter.fGetNext(p));
   }
  
   
   out << "!" << endl

     /* Point Mass Map Block */
     << "POINT MASS MAP                                 2" << endl
     << setw(8) << 0 << "YESNO NO" << endl
     << "!" << endl
     
     /* Flexible Body Map Block */
     << "FLEXIBLE BODY MAP                              2" << endl
     << setw(8) << 0 << "YESNO NO" << endl
     << "!" << endl
     
     /* User Differential Equations Map Block */
     << "USER DIFFERENTIAL EQUATIONS MAP                2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Joints Map Block */
     << "JOINTS MAP                                     2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Joint Primitives Map Block */
     << "JOINT PRIMITIVES MAP                           2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Motion Generators Map Block */
     << "MOTION GENERATORS MAP                          2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Gears Map Block */
     << "GEARS MAP                                      2" << endl
     << setw(8) << 0 << "NO" << endl     
     << "!" << endl
     
     /* Couplers Map Block */
     << "COUPLERS MAP                                   2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* User Constraints Map Block */
     << "USER CONSTRAINTS MAP                           2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Sforces Map Block */
     << "SFORCES MAP                                    2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Spring Damper Forces Map Block */
     << "SPRINGDAMPER FORCES MAP                        2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Bushings Map Block */
     << "BUSHINGS MAP                                   2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Beams Map Block */
     << "BEAMS MAP                                      2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Fields Map Block */
     << "FIELDS MAP                                     2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* User Define Results Map Block */
     << "USER DEFINED RESULTS MAP                       2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Variable map Block */
     << "VARIABLE MAP                                   2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Pinput Map Block */
     << "PINPUT MAP                                     2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Poutput Map Block */
     << "POUTPUT MAP                                    2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Transfer Function SISO Map Block */
     << "TFSISO MAP                                     2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* LSE Map Block */
     << "LSE MAP                                        2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* GSE Map Block */
     << "GSE MAP                                        2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* PTCV Map Block */
     << "PTCV MAP                                       2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* CVCV Map Block */
     << "CVCV MAP                                       2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Nforces Map Block */
     << "NFORCES MAP                                    2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Vforces Map Block */
     << "VFORCES MAP                                    2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Vtorques Map Block */
     << "VTORQUES MAP                                   2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Gforces map Block */
     << "GFORCES MAP                                    2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Frictions Map Block */
     << "FRICTIONS MAP                                  2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Floating Markers Map Block */
     << "FLOATING MARKERS MAP                           2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Mode Shapes Map Block */
     << "MODE SHAPES MAP                                8" << endl
     << setw(8) << 12 << setw(8) << 1 << endl
     << setw(8) << 0 << setw(8) << 13 << endl
     << setw(8) << 0 << setw(8) << 13 << endl
     << setw(8) << 0 << setw(8) << 13 << endl
     << setw(8) << 0 << setw(8) << 13 << endl
     << setw(8) << 0 << setw(8) << 13 << endl
     << setw(8) << 0 << setw(8) << 13 << endl     
     << "!" << endl
     
     /* Tires Map Block */
     << "TIRES MAP                                      2" << endl
     << setw(8) << 0 << "NO" << endl
     << "!" << endl
     
     /* Name Map Block */
     << "NAME MAP                                " << setw(8) << 2+2*iAdamsOutputParts << endl
     << setw(8) << iAdamsOutputParts << "YES" << endl;
   for (i = 0; i < nStrNodes; i++) {
      StructNode *p = (StructNode *)NodeData[NodeType::STRUCTURAL].ppFirstNode[i];
      unsigned int l = p->GetLabel();
      out 
	<< setw(20) << 2+i << "RIGID PART                 1" << endl
	<< "." << sAdamsModelName << ".PART_" << l << endl;
      
      Vec3 x(p->GetXCurr());
      Vec3 e(EulerAngles(p->GetRCurr()));

      /* Commento con il nome del nodo */
      const char *sName = p->GetName();
      if (sName != NULL) {
	 cmd << "! --- " << sName << endl;
      } else {
	 cmd << "! --- Structural Node " << l << endl;
      }
      cmd 
	<< "PART_" << l << endl
	<< 2+i << " " << x << " " << e << endl;
     
#if 0 
      cmd
	<< "!" << endl
	<< "part create rigid_body name_and_position  &" << endl
	<< "   part_name = ." << sAdamsModelName << ".PART_" << 1+i << "  &" << endl
	<< "   adams_id = " << 2+i << "  &" << endl
	<< "   location = ", x.Write(cmd, ", ") << "  &" << endl
	<< "   orientation = ", e.Write(cmd, ", ") << endl
	<< "!" << endl
	<< "marker create &" << endl
	<< "   marker_name = ." << sAdamsModelName << ".PART_" << 1+i << ".MAR_1  &" << endl
	<< "   adams_id = " << 2+i << "  &" << endl
	<< "   location = ", x.Write(cmd, ", ") << "  &" << endl
	<< "   orientation = ", e.Write(cmd, ", ") << endl
	<< "!" << endl
	<< "geometry create shape ellipsoid  &" << endl
	<< "   ellipsoid_name = ." << sAdamsModelName << ".PART_" << 1+i << ".SPHERE_1  &" << endl
	<< "   center_marker = ." << sAdamsModelName << ".PART_" << 1+i << ".MAR_1  &" << endl
	<< "   x_scale_factor = 0.1414213562  &" << endl
	<< "   y_scale_factor = 0.1414213562  &" << endl
	<< "   z_scale_factor = 0.1414213562" << endl
	<< "!" << endl
	<< "part attributes  &" << endl
	<< "   part_name = ." << sAdamsModelName << ".PART_" << 1+i << "  &" << endl
	<< "   color = MAIZE  &" << endl
	<< "   name_visibility = off" << endl;
#endif
      
   }
   
   if (ElemIter.fGetFirst(p)) {
      do {
	 ASSERT(p != NULL);
	 if (p != NULL) {
	    unsigned int nParts = p->iGetNumAdamsDummyParts();
	    if (nParts > 0) {
	       unsigned int l = p->GetLabel();
	       const char *s = psAdamsElemCode[p->GetElemType()];	     

	       /* Nomi e dati delle parti */
	       for (unsigned int part = 1; part <= p->iGetNumAdamsDummyParts(); part++, i++) {
		  /* file intermedio per generazione .cmd */
		  const char *sName = p->GetName();
		  if (sName != NULL) {
		     cmd << "! --- " << sName << " (part " << part << ")" << endl;
		  } else {
		     cmd << "! --- " << psElemNames[p->GetElemType()] << " " << l << " (part " << part << ")" << endl;
		  }
		  p->WriteAdamsDummyPartCmd(cmd, part, 2+i);
		  
		  /* nome parte */
		  out 
		    << setw(20) << 2+i << "RIGID PART                 1" << endl
		    << "." << sAdamsModelName << "." << s << "_" << l << "_" << part << endl;
	       }
	    }
	 }
      } while (ElemIter.fGetNext(p));
   }
   
   out 
     << "!" << endl;
}

void 
DataManager::AdamsResOutput(integer iBlock, const char *type, const char *id) const
{
   ostream& out = OutHdl.AdamsRes();

   unsigned int nStrNodes = NodeData[NodeType::STRUCTURAL].iNum;
   
   out 
     << "ANALYSIS OUTPUT BLOCK                   " << setw(8) << 5+iAdamsOutputParts << endl
     << setw(8) << iAdamsOutputBlock;
#if HAVE_FORM_IN_OSTREAM
   out.form("%-20s", type);
#endif /* HAVE_FORM_IN_OSTREAM */
   out
     << id << endl;
#if HAVE_FORM_IN_OSTREAM
   out.form("%12.5e\n", pTime->GetVal().GetReal());
#endif /* HAVE_FORM_IN_OSTREAM */
   out
     << endl
     << endl;

   for (unsigned int i = 0; i < nStrNodes; i++) {
      Vec3 x;
      doublereal e0;
      Vec3 e;
      
      StructNode *pNode = 
	(StructNode *)NodeData[NodeType::STRUCTURAL].ppFirstNode[i];
      
      x = pNode->GetXCurr();
      EulerParams(pNode->GetRCurr(), e0, e);
#if HAVE_FORM_IN_OSTREAM
      out.form("%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n",
	       x.dGet(1), x.dGet(2), x.dGet(3),
	       e0, e.dGet(1), e.dGet(2), e.dGet(3));
#endif /* HAVE_FORM_IN_OSTREAM */
   }

   Elem *p = NULL;
   if (ElemIter.fGetFirst(p)) {
      do {
	 ASSERT(p != NULL);
	 if (p != NULL) {
	    for (unsigned int part = 1; part <= p->iGetNumAdamsDummyParts(); part++) {
	       Vec3 x;
	       Mat3x3 R;
	       doublereal e0;
	       Vec3 e;
	       p->GetAdamsDummyPart(part, x, R);
	       EulerParams(R, e0, e);
#if HAVE_FORM_IN_OSTREAM
	       out.form("%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e%12.5e\n",
			x.dGet(1), x.dGet(2), x.dGet(3),
			e0, e.dGet(1), e.dGet(2), e.dGet(3));
#endif /* HAVE_FORM_IN_OSTREAM */
	    }
	 }
      } while (ElemIter.fGetNext(p));
   }
   
   out << "!" << endl;
}
