#include "ScalarFunctions.h"
#include "simentity.h"

class BasicFriction : public SimulationEntity{
public:
/*
 * 	unsigned int iGetNumDof(void) const;
 * 	DofOrder::Order GetDofType(unsigned int i) const;
 * 	DofOrder::Order GetEqType (unsigned int i) const;
 * 	void SetValue(VectorHandler&X, VectorHandler&XP) const;
 * 	void BeforePredict(VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&,
 * 		VectorHandler&) const;
 * 	void AfterPredict(VectorHandler&X, VectorHandler&XP);
 * 	void Update(const VectorHandler&XCurr, const VectorHandler&XPrimeCurr);
 * 	void AfterConvergence(VectorHandler&X, VectorHandler&XP);
 */
	virtual doublereal fc(void) const = 0;
	virtual void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP);
};

class ModLugreFriction : public BasicFriction {
private:
	const doublereal sigma0;
	const doublereal sigma1;
	const doublereal sigma2;
	const doublereal kappa;
	doublereal f;
	const DifferentiableScalarFunction & fss;
	doublereal alpha(const doublereal z,
		const doublereal x1) const;
	doublereal alphad_x1(const doublereal z,
		const doublereal x1) const;
	doublereal alphad_z(const doublereal z,
		const doublereal x1) const;
public:
	ModLugreFriction(
		const doublereal sigma0,
		const doublereal sigma1,
		const doublereal sigma2,
		const doublereal kappa,
		const BasicScalarFunction *const f);
	unsigned int iGetNumDof(void) const;
	DofOrder::Order GetDofType(unsigned int i) const;
	DofOrder::Order GetEqType (unsigned int i) const;
	doublereal fc(void) const;
	void AssRes(
		SubVectorHandler& WorkVec,
		const unsigned int startdof,
		const doublereal F,
		const doublereal v,
		const VectorHandler& X,
		const VectorHandler& XP);
};
