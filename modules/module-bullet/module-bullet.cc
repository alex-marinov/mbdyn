/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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

#include <iostream>
#include <cfloat>

#include "dataman.h"
#include "userelem.h"

//include common Bullet Collision Detection headerfiles
#include <btBulletCollisionCommon.h>


// not required
/*
#include "LinearMath/btIDebugDraw.h"
#include "BulletDynamics/Dynamics/btDynamicsWorld.h"

#include "BulletDynamics/ConstraintSolver/btPoint2PointConstraint.h"//picking
#include "BulletDynamics/ConstraintSolver/btGeneric6DofConstraint.h"//picking

#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btCompoundShape.h"
#include "BulletCollision/CollisionShapes/btUniformScalingShape.h"
#include "BulletDynamics/ConstraintSolver/btConstraintSolver.h"*/
#include "matvec3.h"
#include "strnode.h"

// ModuleBullet: begin

class ModuleBullet: virtual public Elem, public UserDefinedElem
{private:
	// add private data
    const StructNode *pNodeA;
    const StructNode *pNodeB;
    
    double scene_size;
    unsigned int max_objects;
    
    btCollisionConfiguration* bt_collision_configuration;
    btCollisionDispatcher* bt_dispatcher;
    btBroadphaseInterface* bt_broadphase;
    btCollisionWorld* bt_collision_world;
    
    Vec3 m_ptA;
    Vec3 m_ptB;
    
    btCollisionObject	objects[];
   
    
public:
	ModuleBullet(unsigned uLabel, const DofOwner *pDO,DataManager* pDM, MBDynParser& HP);
    
    
	virtual ~ModuleBullet(void);
    
	virtual void Output(OutputHandler& OH) const;
    
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
    
	VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,doublereal dCoef,const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
    
	SubVectorHandler& AssRes(SubVectorHandler& WorkVec,doublereal dCoef,const VectorHandler& XCurr,const VectorHandler& XPrimeCurr);
    
	unsigned int iGetNumPrivData(void) const;
    
	unsigned int iGetPrivDataIdx(const char *s) const;
	
    doublereal dGetPrivData(unsigned int i) const;
	
    int iGetNumConnectedNodes(void) const;
	
    void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	
    void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,SimulationEntity::Hints *ph);
	
    std::ostream& Restart(std::ostream& out) const;
	
    virtual unsigned int iGetInitialNumDof(void) const;
	
    virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	
    VariableSubMatrixHandler& InitialAssJac(VariableSubMatrixHandler& WorkMat,const VectorHandler& XCurr);
   	
    SubVectorHandler& InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
};

ModuleBullet::ModuleBullet(unsigned uLabel, const DofOwner *pDO,DataManager* pDM, MBDynParser& HP): Elem(uLabel, flag(0)),UserDefinedElem(uLabel, pDO)
{
    
	// help
	if (HP.IsKeyWord("help")) {
		silent_cout(
                    "									\n"
                    "Module: 	bullet							\n"
                    "Author: 	Vivek Kumar <vivekkumar0893@gmail.com>			\n"
                    "		Pierangelo Masarati <pierangelo.masarati@polimi.it>	\n"
                    "Organization:	Dipartimento di Scienze e Tecnologie Aerospaziali	\n"
                    "		Politecnico di Milano					\n"
                    "		http://www.aero.polimi.it/				\n"
                    "									\n"
                    "All rights reserved							\n"
                    "									\n"
                    "Syntax:								\n"
                    "	user defined : <label> , bullet , ... ;				\n"
                    << std::endl);
        
		if (!HP.IsArg()) {
			/*
			 * Exit quietly if nothing else is provided
			 */
			throw NoErr(MBDYN_EXCEPT_ARGS);
		}
	}
    
    pNodeA = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
    pNodeB = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);
    
	scene_size=10000;
	max_objects=2;
    btScalar sscene_size = (btScalar) scene_size;
    btVector3 worldAabbMin(-sscene_size, -sscene_size, -sscene_size);
    btVector3 worldAabbMax(sscene_size, sscene_size, sscene_size);
    
    bt_collision_configuration = new btDefaultCollisionConfiguration();
    bt_dispatcher = new btCollisionDispatcher(bt_collision_configuration);
    
    bt_broadphase = new bt32BitAxisSweep3(worldAabbMin, worldAabbMax, max_objects, 0, true);
    
    bt_collision_world = new btCollisionWorld(bt_dispatcher, bt_broadphase, bt_collision_configuration);
    
    //btCollisionObject* sphere_A = new btCollisionObject();
    //btCollisionObject* sphere_B = new btCollisionObject();
    
    
	//btAxisSweep3*	broadphase = new btAxisSweep3(worldAabbMin,worldAabbMax);
    
    
    ////////////////////////////////////////////////////////////////
    // A box shape has been used in this implementation.
    // But the connection of object[0] with Node A and of object[1] with Node B is still assumed.
    /////////////////////////////////////////////////////////////
    
    
    
    btMatrix3x3 basisA;
	basisA.setIdentity();
    
	btMatrix3x3 basisB;
	basisB.setIdentity();
    
	objects[0].getWorldTransform().setBasis(basisA);
	objects[1].getWorldTransform().setBasis(basisB);
    
	btBoxShape* boxA = new btBoxShape(btVector3(1,1,1));
	boxA->setMargin(0.f);
    
	btBoxShape* boxB = new btBoxShape(btVector3(0.5,0.5,0.5));
	boxB->setMargin(0.f);
    
    
	objects[0].setCollisionShape(boxA);//&hullA;
	objects[1].setCollisionShape(boxB);//&hullB;
    
    bt_collision_world->addCollisionObject(&objects[0]);
	bt_collision_world->addCollisionObject(&objects[1]);
    
	SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
}

ModuleBullet::~ModuleBullet(void)
{
	NO_OP;
}

void
ModuleBullet::Output(OutputHandler& OH) const
{
	if (fToBeOutput()) {
		 std::ostream& out = OH.Loadable();
        
        out << std::setw(8) << GetLabel()	// 1:	label
        << " " << m_ptA		// A
        << " " << m_ptB	// B
<< std::endl;

	}
}

Vec3
BTtoMBDyn(const btVector3& in)
{
	// check!
	return Vec3(in.getX(), in.getY(), in.getZ());
}

btVector3
MBDyntoBT(const Vec3& in)
{
	return btVector3(in(1), in(2), in(3));
}



Mat3x3
BTtoMBDynMat(const btMatrix3x3& in)
{
	// check!
    btVector3 one= in.getRow(0);
    btVector3 two= in.getRow(1);
    btVector3 three = in.getRow(2);
    
    Vec3 one1 = BTtoMBDyn(one);
    Vec3 two1 = BTtoMBDyn(two);
    Vec3 three1 = BTtoMBDyn(three);
    
	return Mat3x3(one1, two1 , three1);
}

btMatrix3x3
MBDyntoBTMat(const Mat3x3 in)
{
	return btMatrix3x3(in(1,1), in(1,2), in(1,3),in(2,1), in(2,2), in(2,3),in(3,1), in(3,2), in(3,3));
}





void
ModuleBullet::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	// set *piNumRows to 6 * number_of_nodes
	// set *piNumCols to 1 for residual only
	*piNumRows = 12;
	*piNumCols = 12;
}

VariableSubMatrixHandler& ModuleBullet::AssJac(VariableSubMatrixHandler& WorkMat,doublereal dCoef,const VectorHandler& XCurr,const VectorHandler& XPrimeCurr)
{
	WorkMat.SetNullMatrix();
    
	return WorkMat;
}

SubVectorHandler&
ModuleBullet::AssRes(SubVectorHandler& WorkVec,
                     doublereal dCoef,
                     const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr)
{
	WorkVec.ResizeReset(0);
    
    
    
    
   // btDefaultCollisionConfiguration* collisionConfiguration = new btDefaultCollisionConfiguration();
//	btCollisionDispatcher* dispatcher = new btCollisionDispatcher(collisionConfiguration);

    
//	collisionWorld = new btCollisionWorld(dispatcher,broadphase,collisionConfiguration);
//	collisionWorld->setDebugDrawer(&debugDrawer);
    
//#ifdef TEST_NOT_ADDING_OBJECTS_TO_WORLD
    //	collisionWorld->addCollisionObject(&objects[0]);
//	collisionWorld->addCollisionObject(&objects[1]);
//#endif //TEST_NOT_ADDING_OBJECTS_TO_WORLD
    
  
    
    
    //Translate.
    objects[0].getWorldTransform().setOrigin(MBDyntoBT(pNodeA->GetXCurr()));
    objects[1].getWorldTransform().setOrigin(MBDyntoBT(pNodeB->GetXCurr()));
    
   
    //Orientation
    objects[0].getWorldTransform().setBasis(MBDyntoBTMat(pNodeA->GetRCurr()));
    objects[1].getWorldTransform().setBasis(MBDyntoBTMat(pNodeB->GetRCurr()));
    
    
    bt_collision_world->performDiscreteCollisionDetection();
    
    int numManifolds = bt_collision_world->getDispatcher()->getNumManifolds();
    
    //For each contact manifold
    for (int i = 0; i < numManifolds; i++) {
        btPersistentManifold* contactManifold = bt_collision_world->getDispatcher()->getManifoldByIndexInternal(i);
        btCollisionObject* obA = (btCollisionObject*)(contactManifold->getBody0());
        btCollisionObject* obB = (btCollisionObject*)(contactManifold->getBody1());
        contactManifold->refreshContactPoints(obA->getWorldTransform(), obB->getWorldTransform());
        int numContacts = contactManifold->getNumContacts();
        
        //For each contact point in that manifold
        for (int j = 0; j < numContacts; j++) {
            
            //Get the contact information
            btManifoldPoint& pt = contactManifold->getContactPoint(j);
            btVector3 ptA = pt.getPositionWorldOnA();
            btVector3 ptB = pt.getPositionWorldOnB();
            double ptdist = pt.getDistance();
            
            Vec3 m_ptA = BTtoMBDyn(ptA);
            Vec3 m_ptB = BTtoMBDyn(ptB);
            
            
        }
    }

    
    // Only contact points detected by now but not forces.
    // No contribution to residual is made as no change in WorkVec.
    
    
	// compute contact forces and write them in contribution to residual
    
    
    
	return WorkVec;
}

unsigned int
ModuleBullet::iGetNumPrivData(void) const
{
	// return number of private data
	return 0;
}

unsigned int
ModuleBullet::iGetPrivDataIdx(const char *s) const
{
    
	// parse string and compute index of requested private data
    
	// shouldn't get here until private data are defined
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

doublereal
ModuleBullet::dGetPrivData(unsigned int i) const
{
	ASSERT(i > 1 && i <= iGetNumPrivData());
    
	// compute requested private data
    
	// shouldn't get here until private data are defined
	throw ErrGeneric(MBDYN_EXCEPT_ARGS);
}

int
ModuleBullet::iGetNumConnectedNodes(void) const
{
    
    
	// return number of connected nodes
	//return 0;
    return 2;
}

void
ModuleBullet::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
	// resize to number of connected nodes
	// fill array with connected nodes
	//connectedNodes.resize(0);
	// connectedNodes[0] = m_pNode;
    connectedNodes.resize(2);
    
	connectedNodes[0] = pNodeA;
	connectedNodes[1] = pNodeB;
    
}

void
ModuleBullet::SetValue(DataManager *pDM,
                       VectorHandler& X, VectorHandler& XP,
                       SimulationEntity::Hints *ph)
{
	NO_OP;
}

std::ostream&
ModuleBullet::Restart(std::ostream& out) const
{
	return out << "# ModuleBullet: not implemented" << std::endl;
}

unsigned int
ModuleBullet::iGetInitialNumDof(void) const
{
	return 0;
}

void
ModuleBullet::InitialWorkSpaceDim(
                                  integer* piNumRows,
                                  integer* piNumCols) const
{
	*piNumRows = 0;
	*piNumCols = 0;
}

VariableSubMatrixHandler&
ModuleBullet::InitialAssJac(
                            VariableSubMatrixHandler& WorkMat,
                            const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);
    
	WorkMat.SetNullMatrix();
    
	return WorkMat;
}

SubVectorHandler& 
ModuleBullet::InitialAssRes(
                            SubVectorHandler& WorkVec,
                            const VectorHandler& XCurr)
{
	// should not be called, since initial workspace is empty
	ASSERT(0);
    
	WorkVec.ResizeReset(0);
    
	return WorkVec;
}

// ModuleBullet: end

extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
	UserDefinedElemRead *rf1 = new UDERead<ModuleBullet>;
    
	if (!SetUDE("bullet", rf1)) {
		delete rf1;
        
		silent_cerr("ModuleBullet: "
                    "module_init(" << module_name << ") "
                    "failed" << std::endl);
        
		return -1;
	}
    
	return 0;
}

