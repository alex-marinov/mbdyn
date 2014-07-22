/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

#ifndef POSREL_H
#define POSREL_H

#include <strnode.h>
#include <submat.h>
#include <drive.h>

/* Posizione relativa
 * contiene un riferimento ad un nodo strutturale 
 * ed una posizione relativa rispetto a questo, piu' gli strumenti
 * per passare alla posizione ed alal velocita' nel sistema assoluto
 * (rotazione) ed in quello globale (traslazione del nodo).
 * In pratica un punto nel sistema locale di un nodo.
 * E' una classe virtuale pura, in quanto la posizione locale e' definita
 * attraverso un metodo virtuale puro che ritorna un vettore 3x1.
 * Sono definite contestualmente le classi derivate di:
 * 
 * - posizione costante 
 * - direzione fissa (usa un drive per la distanza dal nodo)
 * 
 * Sviluppi possibili:
 * - TplDrive<Vec3>
 * - ...
 * 
 * Occorre modificare gli elementi ed i reference per usare queste distanze
 * in modo uniforme e consistente.
 */


/* PosRel - begin */

class PosRel {
 protected:
   StructNode* pNode;
   
 public:
   PosRel(StructNode* pNode);
   virtual ~PosRel(void);   
   
   inline const StructNode* pGetNode(void) const;
   
   virtual const Vec3& GetPosRel(void) const = 0;
   
   virtual inline const Vec3& GetPosAbs(void) const;
   virtual inline const Vec3& GetPosGlob(void) const;   
   virtual inline const Vec3& GetVelAbs(void) const;
   virtual inline const Vec3& GetVelGlob(void) const;   
};

inline const StructNode* PosRel::pGetNode(void) const
{
   return pNode;
}

inline const StructNode* PosRel::GetPosAbs(void) const
{
   return pNode->GetRCurr()*GetPosRel();
}

inline const StructNode* PosRel::GetPosGlob(void) const
{
   return pNode->GetXCurr()+pNode->GetRCurr()*GetPosRel();
}

inline const StructNode* PosRel::GetVelAbs(void) const
{
   return pNode->GetWCurr().Cross(pNode->GetRCurr()*GetPosRel());
}

inline const StructNode* PosRel::GetVelGlob(void) const
{
   return pNode->GetVCurr()
     +pNode->GetWCurr().Cross(pNode->GetRCurr()*GetPosRel());
}

/* PosRel - end */


/* ConstPosRel - begin */

class ConstPosRel : public PosRel {
 protected:
   Vec3 d;
 public:
   ConstPosRel(StructNode* pNode, const Vec3& d);
   virtual ~ConstPosRel(void);
   
   virtual inline const Vec3& GetPosRel(void) const;
};

inline const Vec3& ConstPosRel::GetPosRel(void) const
{
   return d;
}

/* ConstPosRel - end */


/* FixedDirPosRel - begin */

class FixedDirPosRel : public PosRel {
 protected:
   Vec3 v;
   DriveOwner d;

 public:
   FixedDirPosRel(StructNode* pNode, const Vec3& v, const DriveCaller* pDC);
   virtual ~FixedDirPosRel(void);
   
   virtual inline const Vec3& GetPosRel(void) const;
};

inline const Vec3& FixedDirPosRel::GetPosRel(void) const
{
   return v*d.dGet();
}

/* FixedDirPosRel - end */

#endif /* POSREL_H */
