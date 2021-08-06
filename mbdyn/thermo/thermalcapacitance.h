#ifndef THERMALCAPACITANCE_H
#define THERMALCAPACITANCE_H

#include "therm.h"
#include "thermalnode.h"

class ThermalCapacitance : virtual public Thermal {
 private:
   const ThermalNode* pNode1;
   doublereal thermalcapacitance;
 public:
   ThermalCapacitance(unsigned int uL, 
   	const DofOwner* pDO,
   	const ThermalNode* p1, 
	doublereal c, 
	flag fOut);
   
   ~ThermalCapacitance(void);
   
   /* Tipo di elemento termico (usato solo per debug ecc.) */
   virtual Thermal::Type GetThermalType(void) const;

   virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
      
   VariableSubMatrixHandler& AssJac(VariableSubMatrixHandler& WorkMat,
				    doublereal dCoef,
				    const VectorHandler& XCurr, 
				    const VectorHandler& XPrimeCurr);
   
   SubVectorHandler& AssRes(SubVectorHandler& WorkVec,
			    doublereal dCoef,
			    const VectorHandler& XCurr, 
			    const VectorHandler& XPrimeCurr);
   
//    virtual void AfterConvergence(const VectorHandler& X, 
// 		   const VectorHandler& XP);
//    virtual void Output(OutputHandler& OH) const;
   
//    virtual void SetValue(DataManager *pDM,
// 		   VectorHandler& X, VectorHandler& XP,
// 		   SimulationEntity::Hints *ph = 0);

   /* *******PER IL SOLUTORE PARALLELO******** */        
   /* Fornisce il tipo e la label dei nodi che sono connessi all'elemento
      utile per l'assemblaggio della matrice di connessione fra i dofs */
   virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
   /* ************************************************ */

   /* returns the dimension of the component */
	const virtual OutputHandler::Dimensions GetEquationDimension(integer index) const;
};

#endif  /* THERMALCAPACITANCE_H */
