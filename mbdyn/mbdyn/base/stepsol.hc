/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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
 
#ifndef STEPSOL_HC
#define STEPSOL_HC

template<class T> void
StepIntegrator::UpdateLoop(
		const T* const t,
		void (T::* pUpd)(const int DCount,
			const DofOrder::Order Order,
			const VectorHandler* const pSol) const,
		const VectorHandler* const pSol
	) const
{
	DEBUGCOUTFNAME("StepIntegrator::UpdateLoop");
	ASSERT(pDM != NULL);
	
#ifdef USE_SCHUR
   	Dof CurrDof;
	SchurDataManager* pSDM;
	if ((pSDM = dynamic_cast<SchurDataManager*> (pDM)) != 0) {
		
		Dof* pDofs = pSDM->pGetDofsList();
		
		integer iNumLocDofs = pSDM->HowManyDofs(SchurDataManager::LOCAL);
		integer* pLocDofs = pSDM->GetDofsList(SchurDataManager::LOCAL);
		integer iNumIntDofs = pSDM->HowManyDofs(SchurDataManager::INTERNAL);
		integer* pIntDofs = pSDM->GetDofsList(SchurDataManager::INTERNAL);
		/* local dofs */
		int DCount = 0;
		for (int iCntp1 = 0; iCntp1 < iNumLocDofs; iCntp1++) {
			DCount = pLocDofs[iCntp1];
			CurrDof = pDofs[DCount-1];
			(t->*pUpd)(DCount,CurrDof.Order,pSol);
		}

		/* local interface dofs */
		DCount = 0;
		for (int iCntp1 = 0; iCntp1 < iNumIntDofs; iCntp1++) {
			DCount = pIntDofs[iCntp1];
			CurrDof = pDofs[DCount-1];
			(t->*pUpd)(DCount,CurrDof.Order,pSol);
		}

	} else
#endif // USE_SCHUR
	{
		DataManager::DofIterator_const CurrDof = pDofs->begin();
		integer iNumDofs = pDofs->size();

	   	for (int iCntp1 = 1; iCntp1 <= iNumDofs;
				iCntp1++, ++CurrDof)
		{
			ASSERT(CurrDof != pDofs->end());
			(t->*pUpd)(iCntp1, CurrDof->Order, pSol);
   		}
	}
	return;
}

#endif /* STEPSOL_HC */

