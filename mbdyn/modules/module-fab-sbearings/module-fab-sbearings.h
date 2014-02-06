/**
    Library of bearings for "digital fabrication" machines (alpha version) [2013]
    Eduardo Okabe (okabe@unicamp.br)
    Postdoc CNPq at Aero/Polimi
*/

// Hydrodynamic Bearing Model 01 v0 - September 12th, 2013

class HydrodynamicBearing01
: virtual public Elem, public UserDefinedElem {
private:
   const StructNode* pNode1;
   const StructNode* pNode2;
   Mat3x3 R1tilde;
   Mat3x3 R2tilde;
   Vec3 X1tilde;
   Vec3 X2tilde;
   Vec3 XRel;
   Vec3 XPRel;
   Vec3 HForce;
   Vec3 HMoment;
   doublereal cr;
   doublereal d0;
   doublereal l0;
   doublereal wr_min;
   const HydraulicFluid* hFluid;
   enum bModels {
      OCVIRK,
      CAPONE91,
      SOMMERFELD,
   };
   bModels bearing_model;

public:
	HydrodynamicBearing01(unsigned uLabel, const DofOwner *pDO,
		DataManager* pDM, MBDynParser& HP);
	virtual ~HydrodynamicBearing01(void);

	virtual void Output(OutputHandler& OH) const;
	virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
	VariableSubMatrixHandler&
	AssJac(VariableSubMatrixHandler& WorkMat,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	void JacNum(doublereal Xi[], doublereal JacMat[], integer n, integer m);
	SubVectorHandler&
	AssRes(SubVectorHandler& WorkVec,
		doublereal dCoef,
		const VectorHandler& XCurr,
		const VectorHandler& XPrimeCurr);
	void FForce(doublereal Xi[], doublereal Fi[]) const;
	doublereal hi(doublereal cr, doublereal x, doublereal y, doublereal v) const;
    doublereal dhi(doublereal dx, doublereal dy, doublereal v) const;
	unsigned int iGetNumPrivData(void) const;
	int iGetNumConnectedNodes(void) const;
	void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
	void SetValue(DataManager *pDM, VectorHandler& X, VectorHandler& XP,
		SimulationEntity::Hints *ph);
	std::ostream& Restart(std::ostream& out) const;
	virtual unsigned int iGetInitialNumDof(void) const;
	virtual void
	InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
   	VariableSubMatrixHandler&
	InitialAssJac(VariableSubMatrixHandler& WorkMat,
		      const VectorHandler& XCurr);
   	SubVectorHandler&
	InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr);
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
};


