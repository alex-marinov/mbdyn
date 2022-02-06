
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#include "thermalcapacitance.h"

ThermalCapacitance::ThermalCapacitance(unsigned int uL, 
   	const DofOwner* pDO,
	const ThermalNode* p1, 
	doublereal c, flag fOut) : 
Elem(uL, fOut), 
Thermal(uL, pDO, fOut), 
pNode1(p1), 
thermalcapacitance(c) {
	NO_OP;
};
   
ThermalCapacitance::~ThermalCapacitance(void) {};
   
/* Tipo di elemento idraulico (usato solo per debug ecc.) */
Thermal::Type ThermalCapacitance::GetThermalType(void) const {
	return THERMALCAPACITANCE;
};

void ThermalCapacitance::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const {
	*piNumRows = 1;
	*piNumCols = 1;
};

      
VariableSubMatrixHandler& ThermalCapacitance::AssJac(VariableSubMatrixHandler& WorkMat,
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
	integer iNode1ColIndex = pNode1->iGetFirstColIndex() + 1;

	WM.PutRowIndex(1, iNode1RowIndex);
	WM.PutColIndex(1, iNode1ColIndex);

	WM.IncCoef(1, 1, thermalcapacitance);
	return WorkMat;
};
   
SubVectorHandler& ThermalCapacitance::AssRes(SubVectorHandler& WorkVec,
			doublereal dCoef,
			const VectorHandler& XCurr, 
			const VectorHandler& XPrimeCurr) {
	integer iNumRows = 0;
	integer iNumCols = 0;
	WorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);
 
	/* Indici */
	integer iNode1RowIndex = pNode1->iGetFirstRowIndex() + 1;
	
	WorkVec.PutRowIndex(1, iNode1RowIndex);

	double tp = pNode1->dGetXPrime();
	WorkVec.IncCoef(1, -thermalcapacitance * tp);
	
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
void ThermalCapacitance::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const {
     connectedNodes.resize(1);
     connectedNodes[0] = pNode1;
};
//    /* ************************************************ */
// };

const OutputHandler::Dimensions 
ThermalCapacitance::GetEquationDimension(integer index) const {
	// DOF == 0
	return OutputHandler::Dimensions::UnknownDimension;
}

std::ostream&
ThermalCapacitance::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{

	out
		<< "It does not have any DOF" << std::endl;

	return out;
}