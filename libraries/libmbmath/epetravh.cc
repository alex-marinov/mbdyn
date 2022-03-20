#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#ifdef USE_TRILINOS
#include "matvec3.h"
#include "submat.h"
#undef HAVE_BLAS // FIXME: avoid conflicting declaration
#include "epetravh.h"

EpetraVectorHandler::EpetraVectorHandler(integer iSize, const Epetra_Comm& oComm)
     :oEPV(Epetra_Map(iSize, 1, oComm), true)
{
}

EpetraVectorHandler::~EpetraVectorHandler()
{
}

#ifdef DEBUG
void EpetraVectorHandler::IsValid(void) const
{

}
#endif

doublereal* EpetraVectorHandler::pdGetVec(void) const
{
#ifdef DEBUG
     IsValid();
#endif

     return oEPV.Values();
}

integer EpetraVectorHandler::iGetSize(void) const
{
#ifdef DEBUG
     IsValid();
#endif

     return oEPV.GlobalLength();
}

void EpetraVectorHandler::Reset(void)
{
#ifdef DEBUG
     IsValid();
#endif

     oEPV.PutScalar(0.);
}

void EpetraVectorHandler::Resize(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif

     Epetra_Map oNewMap(iNewSize, 1, oEPV.Comm());
     Epetra_Vector oNewEPV(oNewMap, false);

     const integer iSizeCopy = std::min(iNewSize, iGetSize());

     for (integer i = 0; i < iSizeCopy; ++i) {
          oNewEPV[i] = oEPV[i];
     }

     for (integer i = iSizeCopy + 1; i < iNewSize; ++i) {
          oNewEPV[i] = 0.;
     }

     oEPV = oNewEPV;

#ifdef DEBUG
     IsValid();
#endif
}

void EpetraVectorHandler::ResizeReset(integer iNewSize)
{
#ifdef DEBUG
     IsValid();
#endif

     Epetra_Map oMap(iNewSize, 1, oEPV.Comm());
     oEPV = Epetra_Vector(oMap, true);

#ifdef DEBUG
     IsValid();
#endif
}

void EpetraVectorHandler::PutCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     oEPV[iRow - 1] = dCoef;

#ifdef DEBUG
     IsValid();
#endif
}

void EpetraVectorHandler::IncCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     oEPV[iRow - 1] += dCoef;

#ifdef DEBUG
     IsValid();
#endif
}

void EpetraVectorHandler::DecCoef(integer iRow, const doublereal& dCoef)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     oEPV[iRow - 1] -= dCoef;

#ifdef DEBUG
     IsValid();
#endif
}

const doublereal& EpetraVectorHandler::dGetCoef(integer iRow) const
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     return oEPV[iRow - 1];
}

const doublereal& EpetraVectorHandler::operator () (integer iRow) const
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     return oEPV[iRow - 1];
}

doublereal& EpetraVectorHandler::operator () (integer iRow)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow <= iGetSize());

     return oEPV[iRow - 1];
}

void EpetraVectorHandler::Add(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     for (integer i = 1; i <= 3; ++i) {
          oEPV[iRow + i - 1] += v(i);
     }
}

void EpetraVectorHandler::Sub(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     for (integer i = 1; i <= 3; ++i) {
          oEPV[iRow + i - 1] -= v(i);
     }
}

void EpetraVectorHandler::Put(integer iRow, const Vec3& v)
{
#ifdef DEBUG
     IsValid();
#endif

     ASSERT(iRow >= 1);
     ASSERT(iRow + 2 <= iGetSize());

     for (integer i = 1; i <= 3; ++i) {
          oEPV[iRow + i - 1] = v(i);
     }
}

VectorHandler&
EpetraVectorHandler::ScalarAddMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     const integer iSize = iGetSize();

     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] += d * VH(i);
     }

#ifdef DEBUG
     IsValid();
#endif
     return *this;
}

VectorHandler&
EpetraVectorHandler::ScalarAddMul(const VectorHandler& VH, const VectorHandler& VH1,
                                  const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
     VH1.IsValid();
     ASSERT(iGetSize() == VH1.iGetSize());
#endif
     const integer iSize = iGetSize();

     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] = VH(i) + d * VH1(i);
     }

#ifdef DEBUG
     IsValid();
#endif

     return *this;
}

VectorHandler&
EpetraVectorHandler::ScalarMul(const VectorHandler& VH, const doublereal& d)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     
     const integer iSize = iGetSize();
     
     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] = d * VH(i);
     }
     
#ifdef DEBUG
     IsValid();
#endif
     
     return *this;
}

VectorHandler& EpetraVectorHandler::operator+=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     
     const integer iSize = iGetSize();
     
     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] += VH(i);
     }
     
#ifdef DEBUG
     IsValid();
#endif     
     return *this;
}

VectorHandler& EpetraVectorHandler::operator+=(const SubVectorHandler& SubVH)
{
#ifdef DEBUG
     IsValid();
     SubVH.IsValid();
#endif

     SubVH.AddTo(*this);
     
     return *this;
}

VectorHandler& EpetraVectorHandler::operator-=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
     ASSERT(iGetSize() == VH.iGetSize());
#endif
     const integer iSize = iGetSize();
     
     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] -= VH(i);
     }
     
#ifdef DEBUG
     IsValid();
#endif     
     return *this;
}

VectorHandler& EpetraVectorHandler::operator*=(const doublereal &d)
{
#ifdef DEBUG
     IsValid();
#endif
     const integer iSize = iGetSize();
     
     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] *= d;
     }
     
#ifdef DEBUG
     IsValid();
#endif          
     return *this;
}

VectorHandler& EpetraVectorHandler::operator=(const VectorHandler& VH)
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
#endif
     
     ResizeReset(VH.iGetSize());

     const integer iSize = iGetSize();

     for (integer i = 1; i <= iSize; ++i) {
          oEPV[i - 1] = VH(i);
     }
     
#ifdef DEBUG
     IsValid();
#endif
     return *this;
}

doublereal EpetraVectorHandler::Dot(void) const
{
#ifdef DEBUG
     IsValid();
#endif
     
     const integer iSize = iGetSize();
     doublereal dProd = 0.;

     for (integer i = 1; i <= iSize; ++i) {
          const doublereal Xi = oEPV[i - 1];
          dProd += Xi * Xi;
     }

     return dProd;
}

doublereal EpetraVectorHandler::Norm(void) const
{
#ifdef DEBUG
     IsValid();
#endif
     
     return sqrt(Dot());
}

doublereal EpetraVectorHandler::InnerProd(const VectorHandler& VH) const
{
#ifdef DEBUG
     IsValid();
     VH.IsValid();
#endif
     
     ASSERT(iGetSize() == VH.iGetSize());

     const integer iSize = iGetSize();

     doublereal dProd = 0.;

     for (integer i = 1; i <= iSize; ++i) {
          dProd += oEPV[i - 1] * VH(i);
     }

     return dProd;
}
#endif
