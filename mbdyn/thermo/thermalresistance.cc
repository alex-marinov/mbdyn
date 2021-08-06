
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "therm.h"
#include "thermalresistance.h"

ThermalResistance::ThermalResistance(unsigned int uL, 
   	const DofOwner* pDO,
	const ThermalNode* p1, 
	const ThermalNode* p2, 
	doublereal r, flag fOut) : 
Elem(uL, fOut), 
Thermal(uL, pDO, fOut), 
pNode1(p1), 
pNode2(p2), 
thermalresistance(r) {
	NO_OP;
};
   
ThermalResistance::~ThermalResistance(void) {};
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
Thermal::Type ThermalResistance::GetThermalType(void) const {
	return THERMALRESISTANCE;
};

void ThermalResistance::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
	*piNumRows = 2;
	*piNumCols = 2;
};

      
VariableSubMatrixHandler& ThermalResistance::AssJac(VariableSubMatrixHandler& WorkMat,
				doublereal dCoef,
				const VectorHandler& XCurr, 
				const VectorHandler& XPrimeCurr) {
	/* Setta la sottomatrice come piena (e' un po' dispersivo, ma lo jacobiano 
	 * e' complicato */					
	FullSubMatrixHandler& WM = WorkMat.SetFull();
   
	/* Ridimensiona la sottomatrice in base alle esigenze */
	integer iNumRows = 0;
	integer iNumCols = 0;

	WorkSpaceDim(&iNumRows, &iNumCols);
	WM.ResizeReset(iNumRows, iNumCols);

	integer iNode1RowIndex = pNode1->iGetFirstRowIndex() + 1;
	integer iNode2RowIndex = pNode2->iGetFirstRowIndex() + 1;
	integer iNode1ColIndex = pNode1->iGetFirstColIndex() + 1;
	integer iNode2ColIndex = pNode2->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNode1RowIndex);
	WM.PutColIndex(1, iNode1ColIndex);
	WM.PutRowIndex(2, iNode2RowIndex);
	WM.PutColIndex(2, iNode2ColIndex);
	
	doublereal ri = 1. / thermalresistance * dCoef;

	WM.IncCoef(1, 1, ri);
	WM.DecCoef(1, 2, ri);
	WM.DecCoef(2, 1, ri);
	WM.IncCoef(2, 2, ri);

	return WorkMat;
};
   
SubVectorHandler& ThermalResistance::AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr) {
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);
 
	/* Indici */
	integer iNode1RowIndex = pNode1->iGetFirstRowIndex() + 1;
	integer iNode2RowIndex = pNode2->iGetFirstRowIndex() + 1;
	
	WorkVec.PutRowIndex(1, iNode1RowIndex);
	WorkVec.PutRowIndex(2, iNode2RowIndex);

	double t1 = pNode1->dGetX();
	double t2 = pNode2->dGetX();
	
	doublereal q21 = (t2 - t1) / thermalresistance;
	
	WorkVec.IncCoef(1, q21);
	WorkVec.DecCoef(2, q21);
	
	return WorkVec;
};
   
//    virtual void AfterConvergence(const VectorHandler& X, 
// 		   const VectorHandler& XP);
//    virtual void Output(OutputHandler& OH) const;
//    
//    virtual void SetValue(DataManager *pDM,
// 		   VectorHandler& X, VectorHandler& XP,
// 		   SimulationEntity::Hints *ph = 0);
// 
//    /* *******PER IL SOLUTORE PARALLELO******** */        
//    /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
//       utile per l'assemblaggio della matrice di connessione fra i dofs */
void ThermalResistance::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(2);
     connectedNodes[0] = pNode1;
     connectedNodes[1] = pNode2;
};
//    /* ************************************************ */
// };

const OutputHandler::Dimensions \
ThermalResistance::GetEquationDimension(integer index) const {
	// DOF == 0
	return OutputHandler::Dimensions::UnknownDimension;
}