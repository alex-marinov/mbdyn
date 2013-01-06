/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/* 
 * These routines are part of DataManager, but are especially meant
 * to prepare textual output for Altair's MotionView
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_MOTIONVIEW

#include "dataman.h"
#include "dataman_.h"

#define USE_MOTION_VIEW_C_API

#ifdef USE_MOTION_VIEW_C_API
#include <mbs_result_c_api.h>
#define MRF(f) (f ## _c)
#else /* USE_MOTION_VIEW_C_API */
#include <mbs_result_api.h>
#define MRF(f) (f)
#endif /* USE_MOTION_VIEW_C_API */

enum {
	MRF_RIGID_BODY = 0,
	MRF_LAST_TYPE
};
struct MRFDataType_s {
	int	index;
	char	*name;
} MRFDataType[] = {
	{ MRF_RIGID_BODY,	"Rigid Body" },
	{ MRF_LAST_TYPE,	NULL }
};

bool
DataManager::bMotionViewOutput(void) const
{
	return (ResMode & RES_MOTIONVIEW) ? true : false;
}

void 
DataManager::MotionViewResOutputInit(const char *sOutputFileName)
{
	/* FIXME: need initial time, final time and number of steps :( */
#if MRFOPENRESULT_NEEDS_ARGS
	MRF(mrfOpenResult)(0., 1., 1);
#else /* MRFOPENRESULT_NEEDS_ARGS */
	MRF(mrfOpenResult)();
#endif /* MRFOPENRESULT_NEEDS_ARGS */

	/* FIXME: need a file name */
	char mrf_buf[1024];
	char abf_buf[1024];
	char tab_buf[1024];
	snprintf(mrf_buf, sizeof(mrf_buf), "%s.mrf", sOutputFileName);
	snprintf(abf_buf, sizeof(abf_buf), "%s.abf", sOutputFileName);
	snprintf(tab_buf, sizeof(tab_buf), "%s.tab", sOutputFileName);
	MRF(mrfOpenResultHeader)(mrf_buf, abf_buf, tab_buf);
	
	/* Rigid bodies associated to structural nodes */
	MRF(mrfCreateDataType)(MRFDataType[MRF_RIGID_BODY].name);

	/* Rigid bodies associated to structural nodes */
	MRF(mrfOpenDataTypeByName)(MRFDataType[MRF_RIGID_BODY].name);
	MRF(mrfRegisterForRigidBodyAnimation)();
	MRF(mrfRegisterNoGroundBody)();
	MRF(mrfRegisterForTabulatedOutput)();

	/* Mandatory: position and orientation thru Euler params ... */
	MRF(mrfCreateComponent)("x");
	MRF(mrfCreateComponent)("y");
	MRF(mrfCreateComponent)("z");
	MRF(mrfCreateComponent)("e0");
	MRF(mrfCreateComponent)("e1");
	MRF(mrfCreateComponent)("e2");
	MRF(mrfCreateComponent)("e3");

	/* ... plus velocity and angular velocity */
	MRF(mrfCreateComponent)("vx");
	MRF(mrfCreateComponent)("vy");
	MRF(mrfCreateComponent)("vz");
	MRF(mrfCreateComponent)("wx");
	MRF(mrfCreateComponent)("wy");
	MRF(mrfCreateComponent)("wz");

	MRF(mrfCloseDataType)();

	MRF(mrfOpenDataTypeByIndex)(MRFDataType[MRF_RIGID_BODY].index);
	for (unsigned int i = 0; i < NodeData[Node::STRUCTURAL].iNum; i++) {
		StructNode *pStr = (StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

		/* Skip relative frame dummy nodes */
		if (pStr->GetStructNodeType() == StructNode::DUMMY) {
			DummyStructNode *pDmy = (DummyStructNode *)pStr;

			if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
				continue;
			}
		}

		char namebuf[1024];

		const char *sName = pStr->GetName();
		if (sName == NULL) {
			snprintf(namebuf, sizeof(namebuf), "Node %u", pStr->GetLabel());
			sName = namebuf;
		}
#ifdef USE_MOTION_VIEW_C_API
		mrfCreateElement_c_id(sName, pStr->GetLabel());
#else /* USE_MOTION_VIEW_C_API */
		mrfCreateElement(sName, pStr->GetLabel());
#endif /* USE_MOTION_VIEW_C_API */
	}
	MRF(mrfCloseDataType)();

	MRF(mrfCloseResultHeader)();
}

void 
DataManager::MotionViewResOutput(integer iBlock, const char *type,
		const char *id) const
{
	MRF(mrfOpenTimeStep)((float)dGetTime());

	MRF(mrfOpenDataTypeByIndex)(MRFDataType[MRF_RIGID_BODY].index);
	for (unsigned int i = 0; i < NodeData[Node::STRUCTURAL].iNum; i++) {
		StructNode *pStr = (StructNode *)NodeData[Node::STRUCTURAL].ppFirstNode[i];

		/* Skip relative frame dummy nodes */
		if (pStr->GetStructNodeType() == StructNode::DUMMY) {
			DummyStructNode *pDmy = (DummyStructNode *)pStr;

			if (pDmy->GetDummyType() == DummyStructNode::RELATIVEFRAME) {
				continue;
			}
		}

		MRF(mrfMoveToElementByIndex)(i);
		float q[13];

		Vec3 X(pStr->GetXCurr());
		
		q[0] = X.dGet(1);
		q[1] = X.dGet(2);
		q[2] = X.dGet(3);

		doublereal e0;
		Vec3 e;
		MatR2EulerParams(pStr->GetRCurr(), e0, e);

		q[3] = e0;
		q[4] = e.dGet(1);
		q[5] = e.dGet(2);
		q[6] = e.dGet(3);

		Vec3 V(pStr->GetVCurr());

		q[7] = V.dGet(1);
		q[8] = V.dGet(2);
		q[9] = V.dGet(3);

		Vec3 W(pStr->GetWCurr());

		q[10] = W.dGet(1);
		q[11] = W.dGet(2);
		q[12] = W.dGet(3);

		MRF(mrfPutComponentData)(q);
	}
	MRF(mrfCloseDataType)();

	MRF(mrfCloseTimeStep)();
}

void
DataManager::MotionViewResOutputFini(void) const
{
#ifdef USE_MOTION_VIEW_C_API
	mrfCloseResult_c_0();
#else /* USE_MOTION_VIEW_C_API */
	mrfCloseResult();
#endif /* USE_MOTION_VIEW_C_API */
}

/*
 * Emulate API calls while waiting for Linux-compatible version ...
 */
//#define MIMIC_MOTIONVIEW
#ifdef MIMIC_MOTIONVIEW

int
#if MRFOPENRESULT_NEEDS_ARGS
mrfOpenResult(const float start_time, const float end_time, int num_time_steps)
#else /* !MRFOPENRESULT_NEEDS_ARGS */
mrfOpenResult(void)
#endif /* !MRFOPENRESULT_NEEDS_ARGS */

{
	return 0;
}

int
mrfCloseResult(int flag)
{
	return 0;
}

int
mrfOpenResultHeader(const char *mrf_fliename, const char *abf_filename,
		const char *tab_filename)
{
	return 0;
}

int
mrfCloseResultHeader(void)
{
	return 0;
}

int
mrfCreateDataType(const char *type_name)
{
	return 0;
}

int
mrfOpenDataTypeByIndex(const int type_idx)
{
	return 0;
}

int
mrfOpenDataTypeByName(const char* type_name)
{
	return 0;
}

int
mrfRegisterForTabulatedOutput(void)
{
	return 0;
}

int
mrfRegisterForRigidBodyAnimation(void)
{
	return 0;
}

/*	   x,y,z,e0,e1,e2,e3	

*/
int
mrfRegisterNoGroundBody(void)
{
	return 0;
}
/* by default, the api assumes that the first rigid body is the ground body
   and its time history is discarded because they are time-invariant.
   If in your model there is no ground body or you do not want to make the
   ground body the first one in the rigid body element list, then you have
   to call this function. 
*/

int
mrfRegisterForModalFlexBodyAnimation(const int num_modal_coords)
{
	return 0;
}
/*	   x,y,z,e0,e1,e2,e3	
	   all followed by modal coordinates q1, q2, ....
*/

int
mrfRegisterForAbsNodalFlexBodyAnimation(const int num_nodes)
{
	return 0;
}
/*     x1,x2,...,xn,y1,y2,...,yn,z1,z2,...,zn where n=num_nodes  */

int
mrfCloseDataType(void)
{
	return 0;
}

int
mrfCreateComponent(const char *component_name)
{
	return 0;
}

int mrfCreateElement(const char *element_name)
{
	return 0;
}

int
mrfCreateElement(const char *element_name, const int element_id)
{
	return 0;
}

int
mrfOpenTimeStep(const float time)
{
	return 0;
}

int
mrfCloseTimeStep(void)
{
	return 0;
}

int
mrfMoveToElementByIndex(const int element_idx)
{
	return 0;
}

int
mrfMoveToElementByName(const char *element_name)
{
	return 0;
}

int
mrfPutComponentData(const float *component_data)
{
	return 0;
}

#endif /* MIMIC_MOTIONVIEW */

#endif /* USE_MOTIONVIEW */

