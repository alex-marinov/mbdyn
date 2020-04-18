/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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
  AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
  Copyright (C) 2013(-2019) all rights reserved.

  The copyright of this code is transferred
  to Pierangelo Masarati and Paolo Mantegazza
  for use in the software MBDyn as described
  in the GNU Public License version 2.1
*/

// #define GRADIENT_DEBUG 1
// #define MATVEC_DEBUG 1
// #define HYDRO_DEBUG 1

#ifndef CREATE_PROFILE
#define CREATE_PROFILE 0
#endif

#ifndef HYDRO_DEBUG
#define HYDRO_DEBUG DEBUG
#endif

#include <algorithm>
#include <array>
#include <cassert>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <vector>

#ifdef HAVE_CONFIG_H
#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */
#endif /* HAVE_CONFIG_H */

#include <ac/lapack.h>
#include <clock_time.h>
#include <dataman.h>
#include <gauss.h>
#include <gradient.h>
#include <hfluid.h>
#include <matvec.h>
#include <matvecass.h>
#include <presnode.h>
#include <../thermo/thermalnode.h>
#include <userelem.h>
#include <modal.h>

#include "module-hydrodynamic_plain_bearing2.h"

#if defined(USE_AUTODIFF) && __cplusplus >= 201103L
namespace {
    using namespace grad;

#if !defined(DEBUG) && HYDRO_DEBUG > 0
#define HYDRO_ASSERT(expr) assert(expr)
#else
#define HYDRO_ASSERT(expr) ASSERT(expr)
#endif

#if HYDRO_DEBUG >= 2
#define HYDRO_TRACE(expr) static_cast<void>(std::cerr << expr)
#else
#define HYDRO_TRACE(expr) static_cast<void>(0)
#endif

#if CREATE_PROFILE == 1 || HYDRO_DEBUG > 0
    static inline index_type
    iGetGradientVectorSize(doublereal) {
        return 0;
    }

    template <index_type N_SIZE>
    static inline index_type
    iGetGradientVectorSize(const Gradient<N_SIZE>& g) {
        return g.iGetLocalSize();
    }
#endif

#if HYDRO_DEBUG > 0
    template <typename ElementType>
    bool bCheckNumColsWorkSpace(const ElementType* pElem, grad::FunctionCall eFunc, doublereal g, index_type iRowIndex)
    {
        return true;
    }

    template <typename ElementType, index_type N_SIZE>
    bool bCheckNumColsWorkSpace(const ElementType* pElem, grad::FunctionCall eFunc, const Gradient<N_SIZE>& g, index_type iRowIndex);

    template <typename ElementType, typename G, index_type N_ROWS>
    bool bCheckNumColsWorkSpace(const ElementType* pElem, grad::FunctionCall eFunc, const Vector<G, N_ROWS>& v, index_type iRowIndex)
    {
        for (index_type i = 1; i <= v.iGetNumRows(); ++i) {
            if (!bCheckNumColsWorkSpace(pElem, eFunc, v(i), iRowIndex + i - 1)) {
                return false;
            }
        }

        return true;
    }

#define CHECK_NUM_COLS_WORK_SPACE(pElem, eFunc, g, iRowIndex) HYDRO_ASSERT(bCheckNumColsWorkSpace(pElem, eFunc, g, iRowIndex))
#else
#define CHECK_NUM_COLS_WORK_SPACE(pElem, eFunc, g, iRowIndex) static_cast<void>(0)
#endif

    inline void SetDofMap(doublereal& g, LocalDofMap*)
    {
        g = 0;
    }

    template <index_type N_SIZE>
    inline void SetDofMap(Gradient<N_SIZE>& g, LocalDofMap* pDofMap)
    {
        g = Gradient<N_SIZE>(0., pDofMap);
    }

    namespace util {
        template <typename T>
        class variable_ptr
        {
        public:
            explicit variable_ptr(T* p = nullptr)
                :pPointer(p), bOwns(true) {
            }

            variable_ptr(variable_ptr&& oPtr) {
                steal(std::move(oPtr));
            }

            variable_ptr(const variable_ptr&) = delete;

            ~variable_ptr() {
                free();
            }

            variable_ptr& operator=(variable_ptr&& oPtr) {
                free();
                steal(std::move(oPtr));
                return *this;
            }

            void reset(T* p) {
                free();
                pPointer = p;
                bOwns = true;
            }

            variable_ptr& operator=(const variable_ptr&) = delete;

            T* release() {
                bOwns = false;
                return pPointer;
            }

            bool owner() const { return bOwns; }
            T* get() const { return pPointer; }
            operator T*() const { return pPointer; }

            T* operator->() const {
                HYDRO_ASSERT(pPointer != nullptr);
                return pPointer;
            }

            T& operator*() const {
                HYDRO_ASSERT(pPointer != nullptr);
                return *pPointer;
            }

        private:
            void steal(variable_ptr&& oPtr) {
                pPointer = oPtr.pPointer;
                bOwns = oPtr.bOwns;
                oPtr.bOwns = false;
            }

            void free() {
                if (bOwns) {
                    delete pPointer;
                    pPointer = nullptr;
                    bOwns = false;
                }
            }

            T* pPointer;
            bool bOwns;
        };
    }

    class HydroRootBase {
    public:
        enum FrictionLossType {
            FLUID_FRICTION = 0,
            CONTACT_FRICTION = 1
        };

        enum PrivateDataType {
            PD_CLEARANCE,
            PD_TOTAL_DEFORMATION,
            PD_PRESSURE,
            PD_CONT_PRESSURE,
            PD_DENSITY,
            PD_TEMPERATURE,
            PD_F1x,
            PD_F1y,
            PD_F1z,
            PD_M1x,
            PD_M1y,
            PD_M1z,
            PD_F2x,
            PD_F2y,
            PD_F2z,
            PD_M2x,
            PD_M2y,
            PD_M2z
        };
    };

    class Node2D;
    class FluxNode;
    class ThermoHydrNode;
    class HydroNode;
    class HydroUpdatedNode;
    class HydroElement;
    class HydroMesh;
    class HydroRootElement;
    class BearingGeometry;
    class FluidStateBoundaryCond;
    class LubricationGroove;
    class LubricationGrooveMaster;
    class LubricationGrooveSlave;
    class ContactModel;
    class FrictionModel;
    class ComplianceModel;

    class Geometry2D {
    protected:
        explicit Geometry2D(const Vector<doublereal, 2>& x);
    public:
        virtual ~Geometry2D();
        virtual std::unique_ptr<Geometry2D> Clone(const Vector<doublereal, 2>& x) const=0;
        virtual bool bPointIsInside(const Vector<doublereal, 2>& p1) const=0;
        static std::unique_ptr<Geometry2D> Read(HydroRootElement* pRoot, MBDynParser& HP);
        const Vector<doublereal, 2>& GetPosition() const {
            return x;
        }
    protected:
        const Vector<doublereal, 2> x; // position
    };

    class Circle2D: public Geometry2D {
    public:
        Circle2D(const Vector<doublereal, 2>& x, doublereal r);
        virtual std::unique_ptr<Geometry2D> Clone(const Vector<doublereal, 2>& x) const;
        virtual bool bPointIsInside(const Vector<doublereal, 2>& p1) const;

    private:
        const doublereal r;
    };

    class Rectangle2D: public Geometry2D {
    public:
        Rectangle2D(const Vector<doublereal, 2>& x, doublereal w, doublereal h);
        virtual std::unique_ptr<Geometry2D> Clone(const Vector<doublereal, 2>& x) const;
        virtual bool bPointIsInside(const Vector<doublereal, 2>& p1) const;

    private:
        const doublereal w, h;
    };

    class CompleteSurface2D: public Geometry2D {
    public:
        explicit CompleteSurface2D(const Vector<doublereal, 2>& x);
        virtual std::unique_ptr<Geometry2D> Clone(const Vector<doublereal, 2>& x) const;
        virtual bool bPointIsInside(const Vector<doublereal, 2>& p1) const;
    };

    class SurfaceGrid2D: public Geometry2D {
    public:
        explicit SurfaceGrid2D(const Vector<doublereal, 2>& xc,
                               const Vector<doublereal, DYNAMIC_SIZE>& x,
                               const Vector<doublereal, DYNAMIC_SIZE>& z,
                               doublereal tolx,
                               doublereal tolz,
                               const Matrix<int8_t, DYNAMIC_SIZE, DYNAMIC_SIZE>& status);
        virtual std::unique_ptr<Geometry2D> Clone(const Vector<doublereal, 2>& xc) const;
        virtual bool bPointIsInside(const Vector<doublereal, 2>& p1) const;

    private:
        const doublereal tolx, tolz;
        const Vector<doublereal, DYNAMIC_SIZE> x, z;
        const Matrix<int8_t, DYNAMIC_SIZE, DYNAMIC_SIZE> status;
    };

    class LubricationGroove {
    public:
        enum Type { FIXED, MOVING };

    protected:
        LubricationGroove(integer iLabel, std::unique_ptr<Geometry2D>&& pGeometry);

    public:
        virtual ~LubricationGroove();

        integer iGetLabel() const { return iLabel; }
        Geometry2D* pGetGeometry() const { return pGeometry.get(); }

        virtual FluidStateBoundaryCond* pGetBoundaryCond() const=0;
        virtual void AddNode(Node2D*)=0;
        virtual integer iGetNumNodes() const=0;
        virtual FluidStateBoundaryCond* pReleaseBoundaryCond()=0;
        virtual Type GetType() const=0;
        static std::unique_ptr<LubricationGrooveMaster>
        Read(integer iLabel, HydroRootElement* pRoot, BearingGeometry* pGeometry, MBDynParser& HP);

    private:
        integer iLabel;
        std::unique_ptr<Geometry2D> pGeometry;
    };

    class LubricationGrooveMaster: public LubricationGroove {
    public:
        LubricationGrooveMaster(integer iLabel,
                                std::unique_ptr<Geometry2D>&& pGeometry,
                                FluidStateBoundaryCond* pBoundaryCond,
                                enum Type type);

        virtual ~LubricationGrooveMaster();

        virtual void AddNode(Node2D*);

        virtual integer iGetNumNodes() const;

        virtual FluidStateBoundaryCond* pReleaseBoundaryCond();

        virtual Type GetType() const;

        virtual FluidStateBoundaryCond* pGetBoundaryCond() const;

    private:
        enum Type eType;
        integer iNumNodes;
        util::variable_ptr<FluidStateBoundaryCond> pBoundaryCond;
    };

    class LubricationGrooveSlave: public LubricationGroove {
    public:
        LubricationGrooveSlave(LubricationGrooveMaster* pMaster, const Vector<doublereal, 2>& x);

    public:
        virtual ~LubricationGrooveSlave();

        virtual FluidStateBoundaryCond* pGetBoundaryCond() const;

        virtual void AddNode(Node2D* pNode);

        virtual integer iGetNumNodes() const;

        virtual FluidStateBoundaryCond* pReleaseBoundaryCond();

        virtual Type GetType() const;

    private:
        LubricationGrooveMaster* pMaster;
    };

    class Pocket
    {
    public:
        explicit Pocket(std::unique_ptr<Geometry2D>&& pGeometry);
        virtual ~Pocket();
        static std::unique_ptr<Pocket> Read(HydroRootElement* pRoot, MBDynParser& HP, const class CylindricalBearing* pParent);
        virtual void GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const=0;
        virtual void GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const=0;
        virtual void GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const=0;
        virtual void GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const=0;
        virtual void GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const=0;
        virtual void GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const=0;
        virtual std::unique_ptr<Pocket> Clone(const Vector<doublereal, 2>& x) const=0;
        const Geometry2D* pGetGeometry() const { return pGeometry.get(); }

    private:
        std::unique_ptr<const Geometry2D> pGeometry;
    };

    class ConstHeightPocket: public Pocket
    {
    public:
        ConstHeightPocket(std::unique_ptr<Geometry2D>&& pGeometry, doublereal dy);
        virtual ~ConstHeightPocket();

        virtual void GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const;
        virtual void GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const;
        virtual void GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const;
        virtual void GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const;
        virtual void GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const;
        virtual void GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const;
        virtual std::unique_ptr<Pocket> Clone(const Vector<doublereal, 2>& x) const;

    private:
        doublereal dy;
    };

    class RectangularPocket: public Pocket
    {
    public:
        RectangularPocket(std::unique_ptr<Geometry2D>&& pGeometry,
                          const Vector<doublereal, 2>& x,
                          const Vector<doublereal, 2>& z,
                          const Matrix<doublereal, 2, 2>& Deltay);
        virtual ~RectangularPocket();

        virtual void GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const;
        virtual void GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const;
        virtual void GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const;
        virtual void GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const;
        virtual void GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const;
        virtual void GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const;
        virtual std::unique_ptr<Pocket> Clone(const Vector<doublereal, 2>& x) const;

    private:
        template <typename T> inline
        void GetHeightTpl(const Vector<T, 2>& x, T& Deltay) const;

        template <typename T> inline
        void GetHeightDerXTpl(const Vector<T, 2>& x, T& dDeltay_dx) const;

        template <typename T> inline
        void GetHeightDerZTpl(const Vector<T, 2>& x, T& dDeltay_dz) const;

    private:
        const Vector<doublereal, 2> x, z;
        const Matrix<doublereal, 2, 2> f;
        doublereal dfi1_dx, dfi2_dx;
    };

    class SurfaceGrid: public Pocket
    {
    public:
        explicit SurfaceGrid(std::unique_ptr<Geometry2D>&& pGeometry,
                             const Vector<doublereal, DYNAMIC_SIZE>& x,
                             const Vector<doublereal, DYNAMIC_SIZE>& z,
                             const Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE>& f);
        virtual ~SurfaceGrid();
        virtual void GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const;
        virtual void GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const;
        virtual void GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const;
        virtual void GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const;
        virtual void GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const;
        virtual void GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const;
        virtual std::unique_ptr<Pocket> Clone(const Vector<doublereal, 2>& x) const;

    private:
        template <typename T> inline
        index_type iFindGridX(const Vector<T, 2>& xci) const;

        template <typename T> inline
        index_type iFindGridZ(const Vector<T, 2>& zci) const;

        template <typename T> inline
        void GetHeightTpl(const Vector<T, 2>& xci, T& Deltay) const;

        template <typename T> inline
        void GetHeightDerXTpl(const Vector<T, 2>& xci, T& dDeltay_dx) const;

        template <typename T> inline
        void GetHeightDerZTpl(const Vector<T, 2>& xci, T& dDeltay_dz) const;

        const Vector<doublereal, DYNAMIC_SIZE> x, z;
        const Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> f;
    };

    class HelicalGroove: public Pocket
    {
    public:
        explicit HelicalGroove(std::unique_ptr<Geometry2D>&& pGeometry,
                               std::array<std::unique_ptr<DriveCaller>, 2>&& rgProfile,
                               const Matrix<doublereal, 2, 2>& R0,
                               const Vector<doublereal, 2>& x0,
                               doublereal P);
        virtual ~HelicalGroove();
        virtual void GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const;
        virtual void GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const;
        virtual void GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const;
        virtual void GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const;
        virtual void GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const;
        virtual void GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const;
        virtual std::unique_ptr<Pocket> Clone(const Vector<doublereal, 2>& x) const;

    private:
        template <typename T> inline
        void RelativePosition(const Vector<T, 2>& xci, Vector<T, 2>& x) const;

        template <typename T> inline
        void GetHeightTpl(const Vector<T, 2>& xci, T& Deltay) const;

        template <typename T> inline
        void GetHeightDerTpl(const Vector<T, 2>& xci, T& dDeltay_dx, index_type iDirection) const;

        enum ProfileIndex {
            PROFILE_X = 0,
            PROFILE_Z
        };

        std::array<std::unique_ptr<DriveCaller>, 2> rgProfile;
        const Matrix<doublereal, 2, 2> R0;
        const Vector<doublereal, 2> x0;
        const doublereal P;
    };

    template <typename T>
    class KinematicsBoundaryCond {
    public:
        inline
        KinematicsBoundaryCond();

        inline void
        Update(HydroUpdatedNode* pNode,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        inline void GetClearance(T& h, LocalDofMap* pDofMap) const;
        inline void GetClearanceDerTime(T& dh_dt, LocalDofMap* pDofMap) const;
        inline void GetVelocity(Vector<T, 2>& U1, Vector<T, 2>& U2, LocalDofMap* pDofMap) const;
        inline void GetHydraulicVelocity(Vector<T, 2>& U, LocalDofMap* pDofMap) const;
        inline bool GetContactPressure(T& pasp, LocalDofMap* pDofMap) const;
        inline bool GetContactStress(Vector<T, 2>& tauc_0, LocalDofMap* pDofMap) const;
        inline bool GetContactFrictionLossDens(T& Pfc, LocalDofMap* pDofMap) const;

    private:
        bool bContact;
        T h;                    // clearance
        T dh_dt;                // derivative of clearance versus time
        T pasp;                 // asperity contact pressure
        Vector<T, 2> U1;        // velocity at the shaft
        Vector<T, 2> U2;        // velocity at the bearing
        Vector<T, 2> U;         // effective hydraulic velocity
        Vector<T, 2> tauc_0;    // asperity contact shear stress at y=0
        T Pfc;                  // asperity contact power losses per unit area
    };

    class ThermWallBoundCond {
    public:
        enum WallTempIndex {
            TW_SHAFT = 0,
            TW_BEARING = 1
        };

        ThermWallBoundCond();
        virtual ~ThermWallBoundCond();

        void ParseInput(DataManager* pDM, MBDynParser& HP, const HydroRootElement* pParent);

        template <typename G>
        void GetWallTemperature(WallTempIndex eWall, G& T, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const {
            rgNodes[eWall]->GetX(T, dCoef, func, pDofMap);
        }

        template <typename G>
        void AddHeatFlux(WallTempIndex eWall, const G& Qdot, GradientAssVec<G>& WorkVec) const {
            WorkVec.AddItem(rgNodes[eWall]->iGetFirstRowIndex() + 1, Qdot);
        }

    private:
        std::array<ThermalNode*, 2> rgNodes;
    };

    class HydroFluidBase {
    public:
        enum HydraulicType {
            INCOMPRESSIBLE,
            COMPRESSIBLE,
        };

        enum ThermalType {
            ISOTHERMAL,
            NON_ISOTHERMAL,
        };

        enum HeatCapacityType {
            SPEC_HEAT_TRUE = 0,
            SPEC_HEAT_AVERAGED = 1
        };

        enum CavitationState {
            FULL_FILM_REGION,
            CAVITATION_REGION
        };
    };

    class ThermalFluidModel: public HydroFluidBase {
    public:
        explicit ThermalFluidModel(doublereal T0,
                                   doublereal rho0,
                                   doublereal eta0,
                                   doublereal beta);

        void ParseInput(DataManager* pDM, MBDynParser& HP, const HydroRootElement* pParent);

        template <typename U>
        inline U GetViscosityLiquid(const U& T) const;

        template <typename U>
        inline U GetSpecHeatPerVolume(const U& p, const U& T, HeatCapacityType eType) const;

        template <typename U>
        inline U GetDensityLiquid(const U& T, U* drho_dT = nullptr) const;

        template <typename U>
        inline U GetThermalConductivityLiquid(const U& T) const;

        template <typename U>
        U GetThermalConductivityMixture(const U& T, const U& rho) const;

        template <typename U>
        inline U GetSpecificHeatLiquid(const U& p, const U& T, HeatCapacityType eType) const;

        doublereal dGetRefDensity() const { return rho0; }

        doublereal dGetRefViscosity() const { return eta0; }

        doublereal dGetRefTemperature() const { return T0; }

        ThermalType GetThermalType() const { return eType; }

    private:
        bool bValid() const;

        ThermalType eType;
        doublereal T0;
        doublereal cp0;
        const doublereal rho0;
        const doublereal eta0;
        doublereal beta;
        doublereal lambda0;
        doublereal alphalambda;
        doublereal Aeta2_Aeta3;
        doublereal Aeta3;
        doublereal Ac1;
        doublereal Ac2;
        doublereal Ac3;
        doublereal Ac4;
        doublereal Ac5;
        doublereal Alambda;
    };

    class HydroFluid: public HydroFluidBase {
    public:
        static const index_type iNumDof = 2;

        class ErrNotUnique : public ErrNotImplementedYet {
        public:
            ErrNotUnique(MBDYN_EXCEPT_ARGS_DECL) : ErrNotImplementedYet(MBDYN_EXCEPT_ARGS_PASSTHRU) {};
        };

        HydroFluid(doublereal pc, const ThermalFluidModel& oThermModel);
        virtual ~HydroFluid();
        virtual void GetDensity(const doublereal& p,
                                const doublereal& T,
                                doublereal& rho,
                                doublereal* drho_dp = nullptr,
                                doublereal* drho_dT = nullptr) const=0;

        virtual void GetDensity(const Gradient<0>& p,
                                const Gradient<0>& T,
                                Gradient<0>& rho,
                                Gradient<0>* drho_dp = nullptr,
                                Gradient<0>* drho_dT = nullptr) const=0;

        virtual void GetPressure(const doublereal& rho,
                                 const doublereal& T,
                                 doublereal& p,
                                 doublereal* dp_drho = nullptr,
                                 doublereal* dp_dT = nullptr) const=0;

        virtual void GetPressure(const Gradient<0>& rho,
                                 const Gradient<0>& T,
                                 Gradient<0>& p,
                                 Gradient<0>* dp_drho = nullptr,
                                 Gradient<0>* dp_dT = nullptr) const=0;

        virtual void GetViscosity(const doublereal& rho,
                                  const doublereal& T,
                                  doublereal& eta) const=0;

        virtual void GetViscosity(const Gradient<0>& rho,
                                  const Gradient<0>& T,
                                  Gradient<0>& eta) const=0;

        template <typename U>
        void GetSpecificHeat(const U& p, const U& T, const U& rho, U& cp, HeatCapacityType eType) const {
            cp = oThermModel.GetSpecificHeatLiquid(p, T, eType);
        }

        template <typename U>
        void GetThermalConductivity(const U& T, const U& rho, U& lambda) const {
            lambda = oThermModel.GetThermalConductivityMixture(T, rho);
        }

        virtual void
        ThetaToPhysical(const std::array<doublereal, iNumDof>& Theta,
                        const std::array<doublereal, iNumDof>& dTheta_dt,
                        const doublereal& T,
                        const doublereal& dT_dt,
                        doublereal& p,
                        doublereal& dp_dt,
                        doublereal& rho,
                        doublereal& drho_dt) const=0;
        virtual void
        ThetaToPhysical(const std::array<Gradient<0>, iNumDof>& Theta,
                        const std::array<Gradient<0>, iNumDof>& dTheta_dt,
                        const Gradient<0>& T,
                        const Gradient<0>& dT_dt,
                        Gradient<0>& p,
                        Gradient<0>& dp_dt,
                        Gradient<0>& rho,
                        Gradient<0>& drho_dt) const=0;

        virtual doublereal GetTheta0(index_type iDofIndex) const=0;
        virtual CavitationState Cavitation(doublereal& p, doublereal* dp_dt=0) const=0;
        virtual CavitationState Cavitation(Gradient<0>& p, Gradient<0>* dp_dt=0) const=0;
        virtual doublereal dGetRefPressure() const;
        virtual doublereal dGetRefDensity() const;

        doublereal dGetRefTemperature() const {
            return oThermModel.dGetRefTemperature();
        }

        doublereal dGetRefViscosity() const {
            return oThermModel.dGetRefViscosity();
        }

        virtual HydraulicType GetHydraulicType() const=0;

        ThermalType GetThermalType() const {
            return oThermModel.GetThermalType();
        }

    protected:
        const doublereal pc;
        ThermalFluidModel oThermModel;
    };

    class HydroIncompressibleFluid: public HydroFluid {
    public:
        HydroIncompressibleFluid(doublereal pc, const ThermalFluidModel& oThermModel);
        virtual ~HydroIncompressibleFluid();

        virtual void GetDensity(const doublereal& p,
                                const doublereal& T,
                                doublereal& rho,
                                doublereal* drho_dp = nullptr,
                                doublereal* drho_dT = nullptr) const;

        virtual void GetDensity(const Gradient<0>& p,
                                const Gradient<0>& T,
                                Gradient<0>& rho,
                                Gradient<0>* drho_dp = nullptr,
                                Gradient<0>* drho_dT = nullptr) const;

        virtual void GetPressure(const doublereal& rho,
                                 const doublereal& T,
                                 doublereal& p,
                                 doublereal* dp_drho = nullptr,
                                 doublereal* dp_dT = nullptr) const;

        virtual void GetPressure(const Gradient<0>& rho,
                                 const Gradient<0>& T,
                                 Gradient<0>& p,
                                 Gradient<0>* dp_drho = nullptr,
                                 Gradient<0>* dp_dT = nullptr) const;

        virtual void GetViscosity(const doublereal& rho,
                                  const doublereal& T,
                                  doublereal& eta) const;

        virtual void GetViscosity(const Gradient<0>& rho,
                                  const Gradient<0>& T,
                                  Gradient<0>& eta) const;

        virtual void
        ThetaToPhysical(const std::array<doublereal, iNumDof>& Theta,
                        const std::array<doublereal, iNumDof>& dTheta_dt,
                        const doublereal& T,
                        const doublereal& dT_dt,
                        doublereal& p,
                        doublereal& dp_dt,
                        doublereal& rho,
                        doublereal& drho_dt) const;

        virtual void
        ThetaToPhysical(const std::array<Gradient<0>, iNumDof>& Theta,
                        const std::array<Gradient<0>, iNumDof>& dTheta_dt,
                        const Gradient<0>& T,
                        const Gradient<0>& dT_dt,
                        Gradient<0>& p,
                        Gradient<0>& dp_dt,
                        Gradient<0>& rho,
                        Gradient<0>& drho_dt) const;

        virtual doublereal GetTheta0(index_type iDofIndex) const;
        virtual CavitationState Cavitation(doublereal& p, doublereal* dp_dt=0) const;
        virtual CavitationState Cavitation(Gradient<0>& p, Gradient<0>* dp_dt=0) const;
        virtual HydraulicType GetHydraulicType() const;

    private:
        template <typename G>
        inline void
        GetDensityTpl(const G& p,
                      const G& T,
                      G& rho,
                      G* drho_dp,
                      G* drho_dT) const;

        template <typename G>
        inline void
        GetPressureTpl(const G& rho,
                       const G& T,
                       G& p,
                       G* dp_drho,
                       G* dp_dT) const;

        template <typename G>
        inline CavitationState
        CavitationTpl(G& p, G* dp_dt) const;

        template <typename G>
        inline void
        GetViscosityTpl(const G& rho, const G& T, G& eta) const;

        template <typename G>
        inline void
        ThetaToPhysicalTpl(const G& Theta,
                           const G& dTheta_dt,
                           const G& T,
                           const G& dT_dt,
                           G& p,
                           G& dp_dt,
                           G& rho,
                           G& drho_dt) const;
    };

    class LinearCompressibleFluid: public HydroFluid {
    public:
        LinearCompressibleFluid(doublereal etavap_etaliq,
                                doublereal beta,
                                const doublereal pc,
                                HydraulicType type,
                                const ThermalFluidModel& oThermModel);
        virtual ~LinearCompressibleFluid();

        virtual void GetDensity(const doublereal& p,
                                const doublereal& T,
                                doublereal& rho,
                                doublereal* drho_dp = nullptr,
                                doublereal* drho_dT = nullptr) const;

        virtual void GetDensity(const Gradient<0>& p,
                                const Gradient<0>& T,
                                Gradient<0>& rho,
                                Gradient<0>* drho_dp = nullptr,
                                Gradient<0>* drho_dT = nullptr) const;

        virtual void GetPressure(const doublereal& rho,
                                 const doublereal& T,
                                 doublereal& p,
                                 doublereal* dp_drho = nullptr,
                                 doublereal* dp_dT = nullptr) const;

        virtual void GetPressure(const Gradient<0>& rho,
                                 const Gradient<0>& T,
                                 Gradient<0>& p,
                                 Gradient<0>* dp_drho = nullptr,
                                 Gradient<0>* dp_dT = nullptr) const;

        virtual void GetViscosity(const doublereal& rho,
                                  const doublereal& T,
                                  doublereal& eta) const;

        virtual void GetViscosity(const Gradient<0>& rho,
                                  const Gradient<0>& T,
                                  Gradient<0>& eta) const;

        virtual void
        ThetaToPhysical(const std::array<doublereal, iNumDof>& Theta,
                        const std::array<doublereal, iNumDof>& dTheta_dt,
                        const doublereal& T,
                        const doublereal& dT_dt,
                        doublereal& p,
                        doublereal& dp_dt,
                        doublereal& rho,
                        doublereal& drho_dt) const;
        virtual void
        ThetaToPhysical(const std::array<Gradient<0>, iNumDof>& Theta,
                        const std::array<Gradient<0>, iNumDof>& dTheta_dt,
                        const Gradient<0>& T,
                        const Gradient<0>& dT_dt,
                        Gradient<0>& p,
                        Gradient<0>& dp_dt,
                        Gradient<0>& rho,
                        Gradient<0>& drho_dt) const;

        virtual doublereal GetTheta0(integer iDofIndex) const;
        virtual CavitationState Cavitation(doublereal& p, doublereal* dp_dt = nullptr) const;
        virtual CavitationState Cavitation(Gradient<0>& p, Gradient<0>* dp_dt = nullptr) const;
        virtual HydraulicType GetHydraulicType() const;

    private:
        template <typename G>
        inline void
        GetDensityTpl(const G& p,
                      const G& T,
                      G& rho,
                      G* drho_dp,
                      G* drho_dT) const;

        template <typename G>
        inline void
        GetPressureTpl(const G& rho,
                       const G& T,
                       G& p,
                       G* dp_drho,
                       G* dp_dT) const;

        template <typename G>
        inline CavitationState
        CavitationTpl(G& p, G* dp_dt) const;

        template <typename G>
        inline void
        GetViscosityTpl(const G& rho, const G& T, G& eta) const;

        template <typename G>
        inline void
        ThetaToPhysicalTpl(const std::array<G, iNumDof>& Theta,
                           const std::array<G, iNumDof>& dTheta_dt,
                           const G& T,
                           const G& dT_dt,
                           G& p,
                           G& dp_dt,
                           G& rho,
                           G& drho_dt) const;

    private:
        const doublereal etavap_etaliq;
        const doublereal beta;
        const HydraulicType type;
    };

    class Node2D {
    public:
        enum NodeType {
            HYDRAULIC_NODE      = 0x01000,
            THERMAL_NODE        = 0x02000,
            FLUX_NODE_X         = 0x04000,
            FLUX_NODE_Z         = 0x08000,
            CORNER_NODE         = 0x10000,
            CENTRAL_NODE        = 0x20000,
            ACTIVE_NODE         = 0x01,
            PASSIVE_NODE        = 0x02,
            COUPLED_NODE        = 0x04,
            COMPUTED_NODE       = 0x08,
            MASTER_NODE         = 0x10,
            SLAVE_NODE          = 0x20,
            UPDATED_NODE        = 0x40,
            INCOMPRESSIBLE_NODE = 0x80,
            COMPRESSIBLE_NODE   = 0x100
        };

        enum NodeMask {
            PHYSICS_MASK      = HYDRAULIC_NODE | THERMAL_NODE | FLUX_NODE_X | FLUX_NODE_Z,
            LOCATION_MASK     = CORNER_NODE | CENTRAL_NODE,
            ACTIVE_MASK       = ACTIVE_NODE | PASSIVE_NODE | COUPLED_NODE | COMPUTED_NODE,
            MASTER_SLAVE_MASK = MASTER_NODE | SLAVE_NODE,
            COMPRESSIBLE_MASK = COMPRESSIBLE_NODE | INCOMPRESSIBLE_NODE
        };

        Node2D(integer iNodeNo,
               const Vector<doublereal, 2>& x,
               HydroMesh* pMesh,
               integer iNodeFlags);

        virtual ~Node2D();

        const Vector<doublereal, 2>&
        GetPosition2D() const {
            return x;
        }

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const=0;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const=0;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap)=0;

        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void DofUpdate(VectorHandler& X, VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X,
                                      const VectorHandler& XP);

        virtual bool
        bGetPrivateData(HydroRootBase::PrivateDataType eType,
                        doublereal& dPrivData) const=0;

        virtual void
        Output(std::ostream& os, unsigned uOutputFlags) const=0;

        bool bIsNodeType(NodeType eNodeMask) const {
            return (iNodeFlags & eNodeMask) != 0;
        }

        integer iGetNodeFlags() const {
            return iNodeFlags;
        }

        NodeType GetNodePhysics() const {
            return static_cast<NodeType>(iNodeFlags & PHYSICS_MASK);
        }

        HydroMesh* pGetMesh() const {
            return pMesh;
        }

        const HydroFluid* pGetFluid() const {
            return pFluid;
        }

        integer iGetNodeNumber() const {
            return iNodeNo;
        }

    private:
        const integer iNodeNo;
        const integer iNodeFlags;
        const Vector<doublereal, 2> x;
        HydroMesh* const pMesh;
        const HydroFluid* const pFluid;
    };

    class FluidStateBoundaryCond {
    public:
        enum Type {
            BC_PRESSURE = 1,
            BC_DENSITY = 2,
            BC_FILLING_RATIO = 3
        };

        enum ExtrapMethod {
            EX_INLET = 0,
            EX_OUTLET = 1,
        };

        explicit FluidStateBoundaryCond(const HydroFluid* pFluid,
                                        Type eType,
                                        ExtrapMethod eExtrapMethod,
                                        unsigned uNodeMask,
                                        std::unique_ptr<DriveCaller>&& pTemp);
        virtual ~FluidStateBoundaryCond();
        virtual doublereal dGetPressure(doublereal h=0.) const=0;
        virtual doublereal dGetPressureDerTime(doublereal h=0., doublereal dh_dt=0.) const=0;
        virtual doublereal dGetDensity(doublereal h=0.) const=0;
        virtual doublereal dGetDensityDerTime(doublereal h=0., doublereal dh_dt=0.) const=0;
        inline doublereal dGetTemperature() const;
        inline doublereal dGetTemperatureDerTime() const;
        virtual void Update();
        Type GetType() const { return eType; }
        ExtrapMethod GetExtrapMethod() const { return eExtrapMethod; }

        bool bIncludeNode(Node2D::NodeType eNodeType) const {
            HYDRO_ASSERT((eNodeType == Node2D::HYDRAULIC_NODE) || (eNodeType == Node2D::THERMAL_NODE));
            return uNodeMask & eNodeType;
        }

        const HydroFluid* pGetFluid() const { return pFluid; }
        inline bool bNeedClearance() const;
        static std::unique_ptr<FluidStateBoundaryCond>
        Read(MBDynParser& HP, const HydroRootElement* pRoot);

    private:
        const HydroFluid* const pFluid;
        const Type eType;
        const ExtrapMethod eExtrapMethod;
        const unsigned uNodeMask;
        std::unique_ptr<DriveCaller> pTempDrv;
        doublereal T, dT_dt;
    };

    class FluidStateFunction: public FluidStateBoundaryCond {
    public:
        explicit FluidStateFunction(const HydroFluid* pFluid,
                                    Type eType,
                                    ExtrapMethod eExtrapMethod,
                                    unsigned uNodeMask,
                                    std::unique_ptr<DriveCaller>&& pPressDens,
                                    std::unique_ptr<DriveCaller>&& pTemp);
        virtual ~FluidStateFunction();
        virtual doublereal dGetPressure(doublereal h) const;
        virtual doublereal dGetPressureDerTime(doublereal h, doublereal dh_dt) const;
        virtual doublereal dGetDensity(doublereal h) const;
        virtual doublereal dGetDensityDerTime(doublereal h, doublereal dh_dt) const;
        virtual void Update();

    private:
        std::unique_ptr<DriveCaller> pPressDensDrv;
        doublereal p;
        doublereal dp_dt;
        doublereal rho;
        doublereal drho_dt;
    };

    class FillingRatioFunction: public FluidStateBoundaryCond {
    public:
        explicit FillingRatioFunction(const HydroFluid* pFluid,
                                      Type eType,
                                      ExtrapMethod eExtrapMethod,
                                      unsigned uNodeMask,
                                      std::unique_ptr<DriveCaller>&& pFillRatio,
                                      std::unique_ptr<DriveCaller>&& pTemp);
        virtual ~FillingRatioFunction();
        virtual doublereal dGetPressure(doublereal h) const;
        virtual doublereal dGetPressureDerTime(doublereal h, doublereal dh_dt) const;
        virtual doublereal dGetDensity(doublereal h) const;
        virtual doublereal dGetDensityDerTime(doublereal h, doublereal dh_dt) const;
        virtual void Update();

    private:
        std::unique_ptr<DriveCaller> pFillRatioDrv;
        doublereal h0;
        doublereal dh0_dt;
        doublereal rho;
        doublereal drho_dt;
    };

    class PressureCouplingSlave;
    class PressureCouplingMaster;

    class PressureCouplingCond {
    protected:
        PressureCouplingCond(integer iLabel, std::unique_ptr<Geometry2D>&& pGeometry);

    public:
        virtual ~PressureCouplingCond();

        const Geometry2D* pGetGeometry() const {
            return pGeometry.get();
        }

        virtual PressureNode* pGetNode() const=0;
        virtual ThermalNode* pGetThermalNode() const=0;

        static std::unique_ptr<PressureCouplingMaster>
        Read(integer iLabel, HydroRootElement* pRoot, DataManager* pDM, MBDynParser& HP);

        virtual void AddNode(HydroNode* pNode)=0;
        virtual integer iGetNumNodes() const=0;
        integer iGetLabel() const { return iLabel; }

    private:
        const integer iLabel;
        std::unique_ptr<Geometry2D> pGeometry;
    };

    class PressureCouplingMaster: public PressureCouplingCond {
    public:
        PressureCouplingMaster(integer iLabel,
                               PressureNode* pHydroNode,
                               ThermalNode* pThermalNode,
                               std::unique_ptr<Geometry2D>&& pGeometry);

        virtual ~PressureCouplingMaster();

        virtual PressureNode* pGetNode() const;
        virtual ThermalNode* pGetThermalNode() const;
        virtual void AddNode(HydroNode* pNode);
        virtual integer iGetNumNodes() const;
        std::unique_ptr<PressureCouplingSlave> Clone(const Vector<doublereal, 2>& x);

    private:
        PressureNode* const pHydroNode;
        ThermalNode* const pThermalNode;
        integer iNumNodes;
    };

    class PressureCouplingSlave: public PressureCouplingCond {
    public:
        PressureCouplingSlave(PressureCouplingMaster* pMaster,
                              std::unique_ptr<Geometry2D>&& pGeometry);
        virtual PressureNode* pGetNode() const;
        virtual ThermalNode* pGetThermalNode() const;
        virtual void AddNode(HydroNode* pNode);
        virtual integer iGetNumNodes() const;

    private:
        PressureCouplingMaster* pMaster;
    };

    class HydroDofOwner {
    public:
        enum OffsetIndexType {
            UNKNOWN_OFFSET = -1
        };

        explicit HydroDofOwner();
        virtual ~HydroDofOwner();
        inline integer iGetOffsetIndex(grad::FunctionCall eFunc) const;
        inline void SetOffsetIndex(integer iOffset, grad::FunctionCall eFunc);
        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)=0;
        virtual unsigned int iGetNumDof(void) const=0;
        virtual unsigned int iGetInitialNumDof(void) const=0;
        virtual DofOrder::Order GetDofType(unsigned int i) const=0;
        virtual DofOrder::Order GetEqType(unsigned int i) const=0;
        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const=0;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const=0;

    private:
        inline index_type iFuncCallToIndex(grad::FunctionCall eFunc) const {
            HYDRO_ASSERT((eFunc & grad::REGULAR_FLAG) || (eFunc & grad::INITIAL_ASS_FLAG));
            HYDRO_ASSERT(!((eFunc & grad::REGULAR_FLAG) && (eFunc & grad::INITIAL_ASS_FLAG)));

            return (eFunc & grad::REGULAR_FLAG) ? 1 : 0;
        }

        std::array<integer, 2> rgOffsetIndex;
    };

    class ThermoHydrNode: public Node2D {
    public:
        ThermoHydrNode(integer iNodeNo,
                       const Vector<doublereal, 2>& x,
                       HydroMesh* pMesh,
                       integer iNodeFlags);

        virtual ~ThermoHydrNode();
        virtual void GetTemperature(doublereal& T, LocalDofMap* = nullptr, doublereal = 0.) const=0;
        virtual void GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        virtual void GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* = nullptr, doublereal=0.) const=0;
        virtual void GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        virtual bool bGetPrivateData(HydroRootBase::PrivateDataType eType, doublereal& dPrivData) const;
        virtual void Output(std::ostream& os, unsigned uOutputFlags) const;
    };

    class ThermalActiveNode: public ThermoHydrNode, public HydroDofOwner {
    public:
        ThermalActiveNode(integer iNodeNo,
                          const Vector<doublereal, 2>& x,
                          HydroMesh* pParent,
                          doublereal T0,
                          bool bDoInitAss);
        virtual ~ThermalActiveNode();

        virtual void GetTemperature(doublereal& T, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);
        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual unsigned int iGetNumDof(void) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;

        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

    private:
        grad::FunctionCall eCurrFunc;
        doublereal T, dT_dt;
        const doublereal s;
        const bool bDoInitAss;
    };

    class ThermalCoupledNode: public ThermoHydrNode {
    public:
        ThermalCoupledNode(integer iNodeNo,
                           const Vector<doublereal, 2>& x,
                           HydroMesh* pMesh,
                           ThermalNode* pExtThermNode);

        virtual ~ThermalCoupledNode();
        virtual void GetTemperature(doublereal& T, LocalDofMap* = nullptr, doublereal = 0.) const;
        virtual void GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

    private:
        ThermalNode* const pExtThermNode;
        grad::FunctionCall eCurrFunc;
    };

    class ThermalInletNode: public ThermalActiveNode {
    public:
        ThermalInletNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         HydroMesh* pParent,
                         ThermalNode* pExtThermNode,
                         bool bDoInitAss);

        virtual ~ThermalInletNode();

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual integer
        iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        const ThermalCoupledNode* pGetInletNode() const {
            return &oInletNode;
        }
    private:
        ThermalCoupledNode oInletNode;
    };

    class ThermalPassiveNode: public ThermoHydrNode {
    public:
        ThermalPassiveNode(integer iNodeNo,
                           const Vector<doublereal, 2>& x,
                           HydroMesh* pParent,
                           const FluidStateBoundaryCond* pBoundCond);
        virtual ~ThermalPassiveNode();

        virtual void GetTemperature(doublereal& T, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        const FluidStateBoundaryCond* pGetFluidBoundCond() const { return pBoundCond; }

    private:
        const FluidStateBoundaryCond* const pBoundCond;
    };

    class ThermalSlaveNode: public ThermoHydrNode {
    public:
        ThermalSlaveNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         ThermoHydrNode* pMasterNode);

        virtual ~ThermalSlaveNode();
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;
        virtual void GetTemperature(doublereal& T, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

    private:
        ThermoHydrNode* const pMasterNode;
    };

    class FluxNode: public Node2D {
    public:
        static const index_type iNodeDown = 0;
        static const index_type iNodeUp = 1;
        static const index_type iNumNodes = 2;

        enum PressureSource {
            PRESSURE_FROM_NODE = 0,
            PRESSURE_FROM_MESH = 1,
        };

        enum NodeDataReq {
            ND_NONE = 0x0,
            ND_HYDRAULIC = 0x1,
            ND_THERMAL = ND_HYDRAULIC | 0x2,
            ND_THERMAL_WALL = ND_THERMAL | 0x4
        };

        static const index_type iNumPressSources = PRESSURE_FROM_MESH + 1;

        FluxNode(integer iNodeNo,
                 HydroMesh* pMesh,
                 const std::array<const HydroNode*, iNumNodes>& rgNodes,
                 PressureSource ePressSrc,
                 NodeDataReq eNodeDataReq);

        virtual ~FluxNode();
        inline void GetEnergyBalance(doublereal& Qu, LocalDofMap* = nullptr) const;
        inline void GetEnergyBalance(Gradient<0>& Qu, LocalDofMap* pDofMap) const;
        inline void GetDissipationFactors(doublereal& A0, doublereal& Ah, doublereal& Ac, LocalDofMap* = nullptr) const;
        inline void GetDissipationFactors(Gradient<0>& A0, Gradient<0>& Ah, Gradient<0>& Ac, LocalDofMap* pDofMap) const;
        inline void GetVolumeFluxDens(doublereal& qu, LocalDofMap* = nullptr, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;
        inline void GetVolumeFluxDens(Gradient<0>& qu, LocalDofMap* pDofMap, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;
        inline void GetVelocityAvg(doublereal& wu, LocalDofMap* = nullptr, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;
        inline void GetVelocityAvg(Gradient<0>& wu, LocalDofMap* pDofMap, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;
        inline void GetMassFluxDens(doublereal& mdotu, LocalDofMap* = nullptr, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;
        inline void GetMassFluxDens(Gradient<0>& mdotu, LocalDofMap* pDofMap, PressureSource ePressSrc = PRESSURE_FROM_NODE) const;

        inline void RequestPressureSource(PressureSource ePressSrc);
        inline void RequestNodeData(NodeDataReq eFlag);

        const HydroNode* pGetNode(index_type iNode) const {
            HYDRO_ASSERT(iNode >= 0);
            HYDRO_ASSERT(iNode < iNumNodes);
            return rgNodes[iNode];
        }

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;

        virtual bool
        bGetPrivateData(HydroRootBase::PrivateDataType eType,
                        doublereal& dPrivData) const;

        virtual void
        Output(std::ostream& os, unsigned uOutputFlags) const;

    private:
        static index_type iDirectionFromNodes(const std::array<const HydroNode*, iNumNodes>& rgNodes);

        template <typename T>
        struct NodeData {
            NodeData()
                :Qu{0.},
                 A0{0.},
                 Ah{0.},
                 Ac{0.} {
                 }

            T Qu, A0, Ah, Ac;
        };

        template <typename T>
        struct FluxData {
            FluxData()
                :wu{0.},
                 qu{0.},
                 mdotu{0.} {
                 }

            T wu;
            T qu;
            T mdotu;
        };

        template <typename T>
        inline void UpdateTpl(NodeData<T>& oNode,
                              std::array<FluxData<T>, iNumPressSources>& rgFlux,
                              doublereal dCoef,
                              FunctionCall func,
                              LocalDofMap* pDofMap) const;

        index_type iDirection;
        doublereal du;
        const std::array<const HydroNode*, iNumNodes> rgNodes;
        PressureSource ePressSource;
        unsigned uNodeDataReq;
        NodeData<doublereal> oNode;
        NodeData<Gradient<0> > oNode_grad;
        std::array<FluxData<doublereal>, iNumPressSources> rgFlux;
        std::array<FluxData<Gradient<0> >, iNumPressSources> rgFlux_grad;
        LocalDofMap oDofMap;
    };

    class HydroNode: public Node2D {
    public:
        HydroNode(integer iNodeNo,
                  const Vector<doublereal, 2>& x,
                  HydroMesh* pMesh,
                  integer iNodeFlags);
        virtual ~HydroNode();

        virtual const Vector<doublereal, 3>&
        GetPosition3D() const=0;

        virtual const Matrix<doublereal, 3, 3>&
        GetTangentCoordSys() const=0;
        virtual index_type iGetComplianceIndex() const=0;
        virtual void GetPressure(doublereal& p, LocalDofMap* = nullptr, doublereal=0.) const=0;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const=0;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        virtual void GetDensity(doublereal& rho, LocalDofMap* = nullptr, doublereal=0.) const=0;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        virtual HydroFluid::CavitationState GetCavitationState() const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* = nullptr, doublereal=0.) const=0;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const=0;
        template <typename G>
        inline void GetTemperature(G& T, LocalDofMap* pDofMap = nullptr, doublereal dCoef = 0.) const;
        template <typename G>
        inline void GetTemperatureDerTime(G& dT_dt, LocalDofMap* pDofMap = nullptr, doublereal dCoef = 0.) const;
        template <typename G>
        inline void GetViscosity(G& eta, LocalDofMap* pDofMap = nullptr, doublereal dCoef = 0.) const;
        virtual void GetStress(doublereal& tau_xy_0, doublereal& tau_yz_0, doublereal& tau_xy_h, doublereal& tau_yz_h) const=0;
        virtual void SetStress(doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h)=0;
        virtual bool GetContactPressure(doublereal& pasp, LocalDofMap* = nullptr) const=0;
        virtual bool GetContactPressure(Gradient<0>& pasp, LocalDofMap* pDofMap) const=0;
        virtual void GetContactStress(Vector<doublereal, 2>& tauc_0, LocalDofMap* = nullptr) const=0;
        virtual void GetContactStress(Vector<Gradient<0>, 2>& tauc_0, LocalDofMap* pDofMap) const=0;
        virtual void GetContactFrictionLossDens(doublereal& Pfc, LocalDofMap* = nullptr) const=0;
        virtual void GetContactFrictionLossDens(Gradient<0>& Pfc, LocalDofMap* pDofMap) const=0;
        virtual void GetClearance(doublereal& h, LocalDofMap* = nullptr) const=0;
        virtual void GetClearance(Gradient<0>& h, LocalDofMap* pDofMap) const=0;
        virtual void GetClearanceDerTime(doublereal& dh_dt, LocalDofMap* =nullptr) const=0;
        virtual void GetClearanceDerTime(Gradient<0>& dh_dt, LocalDofMap* pDofMap) const=0;
        virtual void GetRadialDeformation(doublereal& w, doublereal& dw_dt, doublereal dCoef=0, FunctionCall func=REGULAR_RES, LocalDofMap* pDofMap=nullptr) const=0;
        virtual void GetRadialDeformation(Gradient<0>& w, Gradient<0>& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const=0;
        virtual void GetRadialDeformation1(doublereal& w1) const=0;
        virtual void GetRadialDeformation2(doublereal& w2) const=0;
        virtual void GetVelocity(Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, LocalDofMap* =nullptr) const=0;
        virtual void GetVelocity(Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, LocalDofMap* pDofMap) const=0;
        virtual void GetHydraulicVelocity(Vector<doublereal, 2>& U, LocalDofMap* = nullptr) const=0;
        virtual void GetHydraulicVelocity(Vector<Gradient<0>, 2>& U, LocalDofMap* pDofMap) const=0;
        virtual const FluidStateBoundaryCond* pGetMovingPressBoundCond() const=0;
        virtual void SetMovingPressBoundCond(const FluidStateBoundaryCond* pBoundCond)=0;
        doublereal dGetClearance(const FluidStateBoundaryCond* pBoundCond, doublereal* dh_dt=nullptr) const;
        virtual bool bGetPrivateData(HydroRootBase::PrivateDataType eType, doublereal& dPrivData) const;
        virtual void Output(std::ostream& os, unsigned uOutputFlags) const;

        const ThermoHydrNode* pGetThermalNode() const {
            return pThermalNode;
        }

        ThermoHydrNode* pGetThermalNode() {
            return pThermalNode;
        }

        void SetThermalNode(ThermoHydrNode* pNode) {
            HYDRO_ASSERT(pThermalNode == nullptr);
            HYDRO_ASSERT(pNode->GetPosition2D()(1) == GetPosition2D()(1));
            HYDRO_ASSERT(pNode->GetPosition2D()(2) == GetPosition2D()(2));

            pThermalNode = pNode;
        }

    private:
        ThermoHydrNode* pThermalNode;
    };

    class HydroSlaveNode: public HydroNode {
    public:
        HydroSlaveNode(integer iNodeNo,
                       const Vector<doublereal, 2>& x,
                       HydroMesh* pMesh,
                       HydroNode* pMasterNode);
        virtual ~HydroSlaveNode();

        virtual const Vector<doublereal, 3>&
        GetPosition3D() const;

        virtual const Matrix<doublereal, 3, 3>&
        GetTangentCoordSys() const;

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;
        virtual index_type iGetComplianceIndex() const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetStress(doublereal& tau_xy_0, doublereal& tau_yz_0, doublereal& tau_xy_h, doublereal& tau_yz_h) const;
        virtual void SetStress(doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h);
        virtual bool GetContactPressure(doublereal& pasp, LocalDofMap*) const;
        virtual bool GetContactPressure(Gradient<0>& pasp, LocalDofMap* pDofMap) const;
        virtual void GetContactStress(Vector<doublereal, 2>& tauc_0, LocalDofMap* = nullptr) const;
        virtual void GetContactStress(Vector<Gradient<0>, 2>& tauc_0, LocalDofMap* = nullptr) const;
        virtual void GetContactFrictionLossDens(doublereal& Pfc, LocalDofMap* = nullptr) const;
        virtual void GetContactFrictionLossDens(Gradient<0>& Pfc, LocalDofMap* pDofMap) const;
        virtual void GetClearance(doublereal& h, LocalDofMap*) const;
        virtual void GetClearance(Gradient<0>& h, LocalDofMap* pDofMap) const;
        virtual void GetClearanceDerTime(doublereal& dh_dt, LocalDofMap*) const;
        virtual void GetClearanceDerTime(Gradient<0>& dh_dt, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation(doublereal& w, doublereal& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation(Gradient<0>& w, Gradient<0>& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation1(doublereal& w1) const;
        virtual void GetRadialDeformation2(doublereal& w2) const;
        virtual void GetVelocity(Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, LocalDofMap*) const;
        virtual void GetVelocity(Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, LocalDofMap* pDofMap) const;
        virtual void GetHydraulicVelocity(Vector<doublereal, 2>& U, LocalDofMap* =0) const;
        virtual void GetHydraulicVelocity(Vector<Gradient<0>, 2>& U, LocalDofMap* pDofMap) const;
        virtual const FluidStateBoundaryCond* pGetMovingPressBoundCond() const;
        virtual void SetMovingPressBoundCond(const FluidStateBoundaryCond* pBoundCond);

    private:
        HydroNode* pMasterNode;
    };

    class HydroMasterNode: public HydroNode
    {
    public:
        HydroMasterNode(integer iNodeNo,
                        const Vector<doublereal, 2>& x,
                        HydroMesh* pMesh,
                        integer iNodeFlags);
        virtual ~HydroMasterNode();

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual const FluidStateBoundaryCond* pGetMovingPressBoundCond() const;
        virtual void SetMovingPressBoundCond(const FluidStateBoundaryCond* pBoundCond);

    private:
        const FluidStateBoundaryCond* pMovingBoundCond;
    };

    class HydroUpdatedNode: public HydroMasterNode {
    public:
        HydroUpdatedNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         HydroMesh* pMesh,
                         ContactModel* pContactModel,
                         std::unique_ptr<FrictionModel>&& pFrictionModel,
                         integer iNodeFlags);
        virtual ~HydroUpdatedNode();

        virtual const Vector<doublereal, 3>&
        GetPosition3D() const;

        virtual const Matrix<doublereal, 3, 3>&
        GetTangentCoordSys() const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);
        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
        virtual void GetStress(doublereal& tau_xy_0, doublereal& tau_yz_0, doublereal& tau_xy_h, doublereal& tau_yz_h) const;
        virtual void SetStress(doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h);
        const ContactModel* pGetContactModel() const { return pContactModel; }
        const FrictionModel* pGetFrictionModel() const { return pFrictionModel.get(); }
        FrictionModel* pGetFrictionModel(){ return pFrictionModel.get(); }
        virtual bool GetContactPressure(doublereal& pasp, LocalDofMap* =0) const;
        virtual bool GetContactPressure(Gradient<0>& pasp, LocalDofMap* pDofMap) const;
        virtual void GetContactStress(Vector<doublereal, 2>& tauc_0, LocalDofMap* = nullptr) const;
        virtual void GetContactStress(Vector<Gradient<0>, 2>& tauc_0, LocalDofMap* pDofMap) const;
        virtual void GetContactFrictionLossDens(doublereal& Pfc, LocalDofMap* = nullptr) const;
        virtual void GetContactFrictionLossDens(Gradient<0>& Pfc, LocalDofMap* pDofMap) const;
        virtual void GetClearance(doublereal& h, LocalDofMap* pDofMap=0) const;
        virtual void GetClearance(Gradient<0>& h, LocalDofMap* pDofMap) const;
        virtual void GetClearanceDerTime(doublereal& dh_dt, LocalDofMap* pDofMap=0) const;
        virtual void GetClearanceDerTime(Gradient<0>& dh_dt, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation(doublereal& w, doublereal& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation(Gradient<0>& w, Gradient<0>& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const;
        virtual void GetRadialDeformation1(doublereal& w1) const;
        virtual void GetRadialDeformation2(doublereal& w2) const;
        virtual void GetVelocity(Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, LocalDofMap*) const;
        virtual void GetVelocity(Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, LocalDofMap* pDofMap) const;
        virtual void GetHydraulicVelocity(Vector<doublereal, 2>& U, LocalDofMap* =0) const;
        virtual void GetHydraulicVelocity(Vector<Gradient<0>, 2>& U, LocalDofMap* pDofMap) const;
        index_type iGetComplianceIndex() const;

    private:
        Vector<doublereal, 3> v;
        Matrix<doublereal, 3, 3> Rt;
        doublereal sum_tau_xy_0, sum_tau_yz_0, sum_tau_xy_h, sum_tau_yz_h;
        integer iNumStressEval;
        KinematicsBoundaryCond<doublereal> oBoundary;
        KinematicsBoundaryCond<Gradient<0> > oBoundary_grad;
        ContactModel* pContactModel;
        std::unique_ptr<FrictionModel> pFrictionModel;
        ComplianceModel* pComplianceModel;
        index_type iComplianceIndex;
    };

    class HydroIncompressibleNode: public HydroUpdatedNode
    {
    public:
        HydroIncompressibleNode(integer iNodeNo,
                                const Vector<doublereal, 2>& x,
                                HydroMesh* pParent,
                                ContactModel* pContactModel,
                                std::unique_ptr<FrictionModel>&& pFrictionModel,
                                integer iNodeFlags);
        virtual ~HydroIncompressibleNode();
        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

    private:
        template <typename G>
        struct FluidState {
            FluidState()
                :rho(0.),
                 drho_dt(0.) {
            }

            G rho;
            G drho_dt;
        };

        template <typename G>
        inline void
        UpdateState(FluidState<G>& oState, LocalDofMap* pDofMap, doublereal dCoef) const;

        FluidState<doublereal> oState;
        FluidState<Gradient<0> > oState_grad;
    };

    class HydroActiveNode: public HydroIncompressibleNode, public HydroDofOwner {
    public:
        HydroActiveNode(integer iNodeNo,
                        const Vector<doublereal, 2>& x,
                        HydroMesh* pParent,
                        ContactModel* pContactModel,
                        std::unique_ptr<FrictionModel>&& pFrictionModel);
        virtual ~HydroActiveNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);
        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;

        virtual unsigned int iGetNumDof(void) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;

        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

    private:
        doublereal p, dp_dt;
        const doublereal s;
        FunctionCall eCurrFunc;
    };

    class HydroPassiveNode: public HydroIncompressibleNode {
    public:
        HydroPassiveNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         HydroMesh* pParent,
                         ContactModel* pContactModel,
                         std::unique_ptr<FrictionModel>&& pFrictionModel,
                         const FluidStateBoundaryCond* pBoundaryCond);
        virtual ~HydroPassiveNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;

        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;

    private:
        inline const FluidStateBoundaryCond* pGetBoundCond() const;
        const FluidStateBoundaryCond* pBoundaryCond;
    };

    class HydroCoupledNode: public HydroIncompressibleNode {
    public:
        HydroCoupledNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         HydroMesh* pParent,
                         ContactModel* pContactModel,
                         std::unique_ptr<FrictionModel>&& pFrictionModel,
                         const PressureNode* pNode);
        virtual ~HydroCoupledNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;
        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);
    private:
        const PressureNode* const pExtNode;
        grad::FunctionCall eCurrFunc;
    };

    class HydroCompressibleNode: public HydroUpdatedNode {
    public:
        HydroCompressibleNode(integer iNodeNo,
                              const Vector<doublereal, 2>& x,
                              HydroMesh* pParent,
                              ContactModel* pContactModel,
                              std::unique_ptr<FrictionModel>&& pFrictionModel,
                              integer iNodeFlags);
        virtual ~HydroCompressibleNode();
    };

    class HydroActiveComprNode: public HydroCompressibleNode, public HydroDofOwner {
    public:
        HydroActiveComprNode(integer iNodeNo,
                             const Vector<doublereal, 2>& x,
                             HydroMesh* pParent,
                             ContactModel* pContactModel,
                             std::unique_ptr<FrictionModel>&& pFrictionModel);
        virtual ~HydroActiveComprNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual void
        AfterPredict(VectorHandler& X, VectorHandler& XP);

        virtual void
        DofUpdate(VectorHandler& X, VectorHandler& XP);

        virtual void
        AfterConvergence(const VectorHandler& X,
                         const VectorHandler& XP);

        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual HydroFluid::CavitationState GetCavitationState() const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual unsigned int iGetNumDof(void) const;

        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;

        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

        static const index_type iNumDofMax = 2;

    private:
        inline void UpdateTheta(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr);
        inline void UpdateCavitationState();
        inline void ResolveCavitationState(VectorHandler& X, VectorHandler& XP);

        template <typename G>
        struct FluidState {
            FluidState()
                :T{0.},
                 dT_dt{0.},
                 p{0.},
                 rho{0.},
                 drho_dt{0.} {
                     std::fill(Theta.begin(), Theta.end(), 0.); // backward compatible with g++-4.8
                     std::fill(dTheta_dt.begin(), dTheta_dt.end(), 0.);
                 }
            std::array<G, iNumDofMax> Theta;
            std::array<G, iNumDofMax> dTheta_dt;
            G T;
            G dT_dt;
            G p;
            G dp_dt;
            G rho;
            G drho_dt;
        };

        inline void GetTheta(std::array<doublereal, iNumDofMax>& Theta, LocalDofMap* pDofMap, doublereal dCoef) const;
        inline void GetTheta(std::array<Gradient<0>, iNumDofMax>& Theta, LocalDofMap* pDofMap, doublereal dCoef) const;
        inline void GetThetaDerTime(std::array<doublereal, iNumDofMax>& dTheta_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        inline void GetThetaDerTime(std::array<Gradient<0>, iNumDofMax>& dTheta_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

        template <typename G>
        inline void UpdateState(FluidState<G>& oState, LocalDofMap* pDofMap = nullptr, doublereal dCoef = 0.) const;

        FluidState<doublereal> oState;
        FluidState<Gradient<0> > oState_grad;

        struct State {
            doublereal t;
            HydroFluid::CavitationState eCavitationState;
            std::array<doublereal, iNumDofMax> Theta;
            std::array<doublereal, iNumDofMax> dTheta_dt;
        };

        State oRefState, oIncState; // Updated at the first iteration after Solve
        State oCurrState; // Updated in each iteration

        std::array<doublereal, iNumDofMax> s;
        FunctionCall eCurrFunc;
        LocalDofMap oDofMap;
    };

    class HydroPassiveComprNode: public HydroCompressibleNode {
    public:
        HydroPassiveComprNode(integer iNodeNo,
                              const Vector<doublereal, 2>& x,
                              HydroMesh* pParent,
                              ContactModel* pContactModel,
                              std::unique_ptr<FrictionModel>&& pFrictionModel,
                              const FluidStateBoundaryCond* pBoundaryCond);
        virtual ~HydroPassiveComprNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;

        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;

        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual HydroFluid::CavitationState GetCavitationState() const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

        inline const FluidStateBoundaryCond* pGetFluidBoundCond() const;
    private:
        const FluidStateBoundaryCond* pBoundaryCond;
    };

    class HydroComprOutletNode: public HydroPassiveComprNode {
    public:
        HydroComprOutletNode(integer iNodeNo,
                             const Vector<doublereal, 2>& x,
                             HydroMesh* pParent,
                             ContactModel* pContactModel,
                             std::unique_ptr<FrictionModel>&& pFrictionModel,
                             const FluidStateBoundaryCond* pBoundaryCond,
                             const HydroMasterNode* pMasterNode);
        virtual ~HydroComprOutletNode();

        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;

    private:
        const HydroMasterNode* const pMasterNode;
    };

    class HydroCoupledComprNode: public HydroCompressibleNode {
    public:
        HydroCoupledComprNode(integer iNodeNo,
                              const Vector<doublereal, 2>& x,
                              HydroMesh* pParent,
                              ContactModel* pContactModel,
                              std::unique_ptr<FrictionModel>&& pFrictionModel,
                              PressureNode* pExtNode);
        virtual ~HydroCoupledComprNode();

        virtual integer iGetFirstEquationIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetFirstDofIndex(grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;
        virtual void GetPressure(doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void GetPressureDerTime(doublereal& dp_dt, LocalDofMap* = nullptr, doublereal=0.) const;
        virtual void GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef=0.) const;

        virtual void GetDensity(doublereal& rho, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const;
        virtual HydroFluid::CavitationState GetCavitationState() const;
        virtual void GetDensityDerTime(doublereal& drho_dt, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        void Update(const VectorHandler& XCurr,
                    const VectorHandler& XPrimeCurr,
                    doublereal dCoef,
                    enum FunctionCall func,
                    LocalDofMap* pDofMap);

    private:
        template <typename G>
        struct FluidState {
            FluidState()
                :eCavitationState(HydroFluid::FULL_FILM_REGION),
                 p(0.),
                 dp_dt(0.),
                 rho(0.),
                 drho_dt(0.) {
            }

            HydroFluid::CavitationState eCavitationState;
            G p;
            G dp_dt;
            G rho;
            G drho_dt;
        };

        template <typename G>
        inline void UpdateState(FluidState<G>& oStateCurr, LocalDofMap* pDofMap, doublereal dCoef, FunctionCall func) const;
        void GetExtPressure(doublereal& p, doublereal& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        void GetExtPressure(Gradient<0>& p, Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const;
        const PressureNode* const pExtNode;
        doublereal pext;
        doublereal dpext_dt;
        FunctionCall eCurrFunc;

        FluidState<doublereal> oState;
        FluidState<Gradient<0> > oState_grad;
        LocalDofMap oDofMap;
    };

    class BearingGeometry {
    public:
        enum Type {
            CYLINDRICAL_MESH_UNKNOWN = 0x0,
            // Needed for post-processing
            CYLINDRICAL_MESH = 0x1,
            CYLINDRICAL_MESH_AT_SHAFT = CYLINDRICAL_MESH | 0x10,
            CYLINDRICAL_MESH_AT_BEARING = CYLINDRICAL_MESH | 0x20
        };

        explicit
        BearingGeometry(HydroRootElement* pParent);

        virtual
        ~BearingGeometry();

        virtual void
        ParseInput(DataManager* pDM, MBDynParser& HP)=0;

        virtual void Initialize() = 0;

        virtual void
        GetClosestDistance2D(const Vector<doublereal, 2>& x1,
                             const Vector<doublereal, 2>& x2,
                             Vector<doublereal, 2>& dx) const=0;

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              doublereal& h,
                              doublereal& dh_dt,
                              Vector<doublereal, 2>& U1,
                              Vector<doublereal, 2>& U2,
                              Vector<doublereal, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const=0;

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              Gradient<0>& h,
                              Gradient<0>& dh_dt,
                              Vector<Gradient<0>, 2>& U1,
                              Vector<Gradient<0>, 2>& U2,
                              Vector<Gradient<0>, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const=0;

        virtual void
        GetNonNegativeClearance(const doublereal& h,
                                doublereal& hn,
                                const doublereal* dh_dt = nullptr,
                                doublereal* dhn_dt = nullptr) const;
        virtual void
        GetNonNegativeClearance(const Gradient<0>& h,
                                Gradient<0>& hn,
                                const Gradient<0>* dh_dt = nullptr,
                                Gradient<0>* dhn_dt = nullptr) const;

        virtual void
        GetPosition3D(const Vector<doublereal, 2>& x,
                      Vector<doublereal, 3>& v) const=0;

        virtual void
        GetTangentCoordSys(const Vector<doublereal, 2>& x,
                           Matrix<doublereal, 3, 3>& Rt) const=0;

        doublereal
        dGetNodeDistance2D(const Node2D* pNode1,
                           const Node2D* pNode2,
                           index_type iDirection) const;

        virtual void
        GetStructNodeOffset(const HydroNode* pHydroNode, Vector<doublereal, 3>& v) const=0;

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<doublereal, 3>& dF_0_Rt,
                         const Vector<doublereal, 3>& dF_h_Rt,
                         const Vector<doublereal, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap)=0;

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<Gradient<0>, 3>& dF_0_Rt,
                         const Vector<Gradient<0>, 3>& dF_h_Rt,
                         const Vector<Gradient<0>, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap)=0;

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const=0;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const=0;

        virtual void
        Update(doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap)=0;

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode)=0;

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode)=0;

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode)=0;

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode)=0;

        virtual int iGetNumConnectedNodes(void) const=0;
        virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const=0;
        HydroRootElement* pGetParent() const { return pParent; }
        virtual std::ostream& PrintLogFile(std::ostream& os) const=0;
        virtual std::ostream& Output(std::ostream& os) const=0;
        virtual enum Type GetType() const=0;

        virtual enum LubricationGroove::Type
        ReadLubricationGrooveType(MBDynParser& HP) const=0;

        inline const FluidStateBoundaryCond* pGetMovingPressBoundCond(const HydroNode* pNode) const {
            return pNode->pGetMovingPressBoundCond();
        }

        virtual void AddMovingLubrGroove(LubricationGroove* pGroove);
        virtual void ReserveMovingLubrGrooves(size_t n);

        const LubricationGroove*
        pFindMovingLubricationGroove(const Vector<doublereal, 2>& x, Node2D::NodeType eNodeType) const;

        virtual doublereal dGetReferenceClearance() const=0;

        virtual bool
        bGetPrivateData(HydroRootBase::PrivateDataType eType,
                        doublereal& dPrivData) const=0;

        virtual void GetMovingMeshOffset(Vector<doublereal, 2>& x, LocalDofMap* = nullptr) const = 0;
        virtual void GetMovingMeshOffset(Vector<Gradient<0>, 2>& x, LocalDofMap* pDofMap) const = 0;

    protected:
        virtual doublereal dGetMinClearance() const=0;
        virtual doublereal dGetPocketHeightMesh(const Vector<doublereal, 2>& x) const=0;

#if CREATE_PROFILE == 1
        enum { PROF_RES = 0, PROF_JAC = 1};
        struct {
            doublereal dtAddForce[2];
            doublereal dtOperatorPlus[2];
        } profile;
#endif

    private:
        template <typename T> inline void
        GetNonNegativeClearanceTpl(const T& h,
                                   T& hn,
                                   const T* dh_dt,
                                   T* dhn_dt) const;

        HydroRootElement* pParent;

        typedef std::vector<LubricationGroove*> LubrGrooveVectorType;
        typedef LubrGrooveVectorType::const_iterator ConstLubrGrooveIterator;

        LubrGrooveVectorType rgMovingGrooves;
    };

    class RigidBodyBearing: public BearingGeometry {
    public:
        RigidBodyBearing(HydroRootElement* pParent);
        virtual ~RigidBodyBearing();
        virtual void ParseInput(DataManager* pDM, MBDynParser& HP);
        const StructNode* pGetNode1() const { return pNode1; }
        const StructNode* pGetNode2() const { return pNode2; }
        const Vector<doublereal, 3>&
        GetOffsetNode1() const { return o1_R1; }
        const Matrix<doublereal, 3, 3>&
        GetOrientationNode1() const { return Rb1; }
        const Vector<doublereal, 3>&
        GetOffsetNode2() const { return o2_R2; }
        const Matrix<doublereal, 3, 3>&
        GetOrientationNode2() const { return Rb2; }
        virtual int iGetNumConnectedNodes(void) const;
        virtual void GetConnectedNodes(std::vector<const Node *>& connectedNodes) const;
        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc) const;
        virtual std::ostream& PrintLogFile(std::ostream& os) const;
        virtual std::ostream& Output(std::ostream& os) const;

    protected:
        inline void
        SaveReactionForce(const Vector<doublereal, 3>& F1,
                          const Vector<doublereal, 3>& M1,
                          const Vector<doublereal, 3>& F2,
                          const Vector<doublereal, 3>& M2);

        inline void
        SaveReactionForce(const Vector<Gradient<0>, 3>& F1,
                          const Vector<Gradient<0>, 3>& M1,
                          const Vector<Gradient<0>, 3>& F2,
                          const Vector<Gradient<0>, 3>& M2);

    private:
        const StructNode* pNode1;   // shaft node
        const StructNode* pNode2;   // bearing node
        Vector<doublereal, 3> o1_R1;
        Vector<doublereal, 3> o2_R2;
        Matrix<doublereal, 3, 3> Rb1;
        Matrix<doublereal, 3, 3> Rb2;
        Vector<doublereal, 3> F1, M1, F2, M2;
    };

    class CylindricalBearing: public RigidBodyBearing {
    public:
        explicit CylindricalBearing(HydroRootElement* pParent);
        virtual ~CylindricalBearing();
        virtual void ParseInput(DataManager* pDM, MBDynParser& HP);
        doublereal dGetBearingWidth() const { return b; }
        doublereal dGetShaftRadius() const { return r; }
        doublereal dGetBearingRadius() const { return R; }

        virtual void
        GetClosestDistance2D(const Vector<doublereal, 2>& x1,
                             const Vector<doublereal, 2>& x2,
                             Vector<doublereal, 2>& dx) const;

        virtual void
        GetPosition3D(const Vector<doublereal, 2>& x,
                      Vector<doublereal, 3>& v) const;

        virtual void
        GetPosition3D(const Vector<Gradient<0>, 2>& x,
                      Vector<Gradient<0>, 3>& v) const;

        virtual void
        GetTangentCoordSys(const Vector<doublereal, 2>& x,
                           Matrix<doublereal, 3, 3>& Rt) const;

        virtual void
        GetTangentCoordSys(const Vector<Gradient<0>, 2>& x,
                           Matrix<Gradient<0>, 3, 3>& Rt) const;
        virtual void
        GetStructNodeOffset(const HydroNode* pHydroNode, Vector<doublereal, 3>& v) const;

        virtual void Update(doublereal dCoef, enum FunctionCall func, LocalDofMap* pDofMap);
        virtual std::ostream& PrintLogFile(std::ostream& os) const;

        virtual doublereal dGetMeshRadius() const=0; // radius at the mesh side

        virtual doublereal dGetMinClearance() const;
        virtual doublereal dGetReferenceClearance() const;

    protected:
        typedef std::unique_ptr<Pocket> PocketPtr;
        typedef std::vector<PocketPtr, GradientAllocator<PocketPtr> > PocketVector;
        typedef PocketVector::iterator PocketIterator;
        typedef PocketVector::const_iterator ConstPocketIterator;

        const Pocket*
        pFindBearingPocket(const Vector<doublereal, 2>& x) const;

        const Pocket*
        pFindShaftPocket(const Vector<doublereal, 2>& x) const;

        virtual const Pocket*
        pFindMeshPocket(const Vector<doublereal, 2>& x) const=0;

        virtual const Matrix<doublereal, 3, 3>&
        GetOrientationMeshNode() const=0;

        virtual doublereal
        dGetPocketHeightMesh(const Vector<doublereal, 2>& x) const;

    protected:
        LocalDofMap oDofMap;

    private:
        static const Pocket* pFindPocket(const Vector<doublereal, 2>& x, const PocketVector& rgPockets);
        void ReadPockets(MBDynParser& HP, PocketVector& rgPockets);

        template <typename T>
        inline void
        GetTangentCoordSysTpl(const Vector<T, 2>& x,
                              Matrix<T, 3, 3>& Rt) const;

        template <typename T>
        inline void
        GetPosition3DTpl(const Vector<T, 2>& x,
                         Vector<T, 3>& v) const;
    private:
        doublereal b; // bearing width
        doublereal r; // shaft radius
        doublereal R; // bearing radius
        doublereal hmin;
        PocketVector rgPocketsShaft;
        PocketVector rgPocketsBearing;
    };

    class CylindricalMeshAtShaft: public CylindricalBearing {
    public:
        CylindricalMeshAtShaft(HydroRootElement* pParent);
        virtual ~CylindricalMeshAtShaft();
        virtual void Initialize();
        virtual Type GetType() const;

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              doublereal& h,
                              doublereal& dh_dt,
                              Vector<doublereal, 2>& U1,
                              Vector<doublereal, 2>& U2,
                              Vector<doublereal, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const;

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              Gradient<0>& h,
                              Gradient<0>& dh_dt,
                              Vector<Gradient<0>, 2>& U1,
                              Vector<Gradient<0>, 2>& U2,
                              Vector<Gradient<0>, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const;

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<doublereal, 3>& dF_0_Rt,
                         const Vector<doublereal, 3>& dF_h_Rt,
                         const Vector<doublereal, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap);

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<Gradient<0>, 3>& dF_0_Rt,
                         const Vector<Gradient<0>, 3>& dF_h_Rt,
                         const Vector<Gradient<0>, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap);

        virtual void
        Update(doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkMat,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkMat,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual enum LubricationGroove::Type
        ReadLubricationGrooveType(MBDynParser& HP) const;

        virtual bool
        bGetPrivateData(HydroRootBase::PrivateDataType eType,
                        doublereal& dPrivData) const;

        virtual void GetMovingMeshOffset(Vector<doublereal, 2>& x, LocalDofMap* = nullptr) const;
        virtual void GetMovingMeshOffset(Vector<Gradient<0>, 2>& x, LocalDofMap* pDofMap) const;

    protected:
        virtual doublereal dGetMeshRadius() const;

        virtual const Matrix<doublereal, 3, 3>&
        GetOrientationMeshNode() const;

        virtual const Pocket*
        pFindMeshPocket(const Vector<doublereal, 2>& x) const;

    private:
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkMat,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        template <typename T>
        class Boundary {
        public:
            inline
            Boundary(const CylindricalMeshAtShaft& rParent);

            inline void Initialize();

            inline void
            Update(doublereal dCoef,
                   enum FunctionCall func,
                   LocalDofMap* pDofMap);

            inline void
            GetBoundaryConditions(HydroNode* pNode,
                                  T& h,
                                  T& dh_dt,
                                  Vector<T, 2>& U1,
                                  Vector<T, 2>& U2,
                                  Vector<T, 2>& U,
                                  doublereal dCoef,
                                  enum FunctionCall func,
                                  LocalDofMap* pDofMap) const;

            void GetMovingMeshOffset(Vector<T, 2>& x, LocalDofMap* pDofMap) const;

        private:
            typedef Vector<doublereal, 2> CVec2;
            typedef Vector<doublereal, 3> CVec3;
            typedef Matrix<doublereal, 3, 3> CMat3x3;
            typedef Vector<T, 2> VVec2;
            typedef Vector<T, 3> VVec3;
            typedef Matrix<T, 3, 3> VMat3x3;

            const CylindricalMeshAtShaft& rParent;
            CVec3 Rb2T_o2;
            VVec3 X1, X2, X1P, X2P, omega1, omega2;
            VMat3x3 R1, R2;
            VMat3x3 Rb2T_R2T;
        };

        template <typename T>
        struct ReactionForce {
            Vector<T, 3> F1_R1;
            Vector<T, 3> M1_R1;
            Vector<T, 3> F2_R1;
            Vector<T, 3> M2_R1;

            inline void Reset();
        };

        template <typename T>
        inline void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<T, 3>& dF_0_Rt,
                         const Vector<T, 3>& dF_h_Rt,
                         const Vector<T, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap,
                         ReactionForce<T>& oReact);

        inline const
        ReactionForce<doublereal>& GetReactionForce(const doublereal& dummy) const;

        inline const
        ReactionForce<Gradient<0> >& GetReactionForce(const Gradient<0>& dummy) const;

        Boundary<doublereal> oBound;
        Boundary<Gradient<0> > oBound_grad;
        ReactionForce<doublereal> oReaction;
        ReactionForce<Gradient<0> > oReaction_grad;
    };

    class CylindricalMeshAtBearing: public CylindricalBearing {
    public:
        CylindricalMeshAtBearing(HydroRootElement* pParent);
        virtual ~CylindricalMeshAtBearing();
        virtual Type GetType() const;
        virtual void Initialize();

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              doublereal& h,
                              doublereal& dh_dt,
                              Vector<doublereal, 2>& U1,
                              Vector<doublereal, 2>& U2,
                              Vector<doublereal, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const;

        virtual void
        GetBoundaryConditions(HydroNode* pNode,
                              Gradient<0>& h,
                              Gradient<0>& dh_dt,
                              Vector<Gradient<0>, 2>& U1,
                              Vector<Gradient<0>, 2>& U2,
                              Vector<Gradient<0>, 2>& U,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap) const;

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<doublereal, 3>& dF_0_Rt,
                         const Vector<doublereal, 3>& dF_h_Rt,
                         const Vector<doublereal, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap);

        virtual void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<Gradient<0>, 3>& dF_0_Rt,
                         const Vector<Gradient<0>, 3>& dF_h_Rt,
                         const Vector<Gradient<0>, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap);

        virtual void
        Update(doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkMat,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkMat,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual enum LubricationGroove::Type
        ReadLubricationGrooveType(MBDynParser& HP) const;

        bool bGetPrivateData(HydroRootBase::PrivateDataType eType,
                             doublereal& dPrivData) const;

        virtual void GetMovingMeshOffset(Vector<doublereal, 2>& x, LocalDofMap* = nullptr) const;
        virtual void GetMovingMeshOffset(Vector<Gradient<0>, 2>& x, LocalDofMap* pDofMap) const;

    protected:
        virtual doublereal dGetMeshRadius() const;

        virtual const Matrix<doublereal, 3, 3>&
        GetOrientationMeshNode() const;

        virtual const Pocket*
        pFindMeshPocket(const Vector<doublereal, 2>& x) const;

    private:
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkMat,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        template <typename T>
        struct Boundary {
            inline
            Boundary(const CylindricalMeshAtBearing& rParent);

            inline void Initialize();

            inline void
            Update(doublereal dCoef,
                   enum FunctionCall func,
                   LocalDofMap* pDofMap);

            inline void
            GetBoundaryConditions(HydroNode* pNode,
                                  T& h,
                                  T& dh_dt,
                                  Vector<T, 2>& U1,
                                  Vector<T, 2>& U2,
                                  Vector<T, 2>& U,
                                  doublereal dCoef,
                                  enum FunctionCall func,
                                  LocalDofMap* pDofMap) const;

            void GetMovingMeshOffset(Vector<T, 2>& x, LocalDofMap* pDofMap) const;

            typedef Vector<doublereal, 2> CVec2;
            typedef Vector<doublereal, 3> CVec3;
            typedef Matrix<doublereal, 3, 3> CMat3x3;
            typedef Vector<T, 2> VVec2;
            typedef Vector<T, 3> VVec3;
            typedef Matrix<T, 3, 3> VMat3x3;

            const CylindricalMeshAtBearing& rParent;
            CVec3 Rb1T_o1;
            VVec3 X1, X2, X1P, X2P, omega1, omega2;
            VMat3x3 R1, R2;
            VMat3x3 Rb1T_R1T;
        };

        template <typename T>
        struct ReactionForce {
            Vector<T, 3> F1_R2;
            Vector<T, 3> M1_R2;
            Vector<T, 3> F2_R2;
            Vector<T, 3> M2_R2;

            inline void Reset();
        };

        template <typename T>
        inline void
        AddReactionForce(const Vector<doublereal, 2>& x,
                         const Vector<doublereal, 3>& v,
                         const Matrix<doublereal, 3, 3>& Rt,
                         const Vector<T, 3>& dF_0_Rt,
                         const Vector<T, 3>& dF_h_Rt,
                         const Vector<T, 2>& dM_h_Rt,
                         doublereal dCoef,
                         enum FunctionCall func,
                         LocalDofMap* pDofMap,
                         ReactionForce<T>& oReact);

        inline const ReactionForce<doublereal>&
        GetReactionForce(const doublereal& dummy) const;

        inline const ReactionForce<Gradient<0> >&
        GetReactionForce(const Gradient<0>& dummy) const;

        Boundary<doublereal> oBound;
        Boundary<Gradient<0> > oBound_grad;
        ReactionForce<doublereal> oReaction;
        ReactionForce<Gradient<0> > oReaction_grad;
    };

    class Material
    {
    public:
        Material()
            :E(0.), nu(0.) {

        }

        // Only Hook materials are supported right now since the theory of the elastic half space is based on them
        Material(doublereal E, doublereal nu)
            :E(E), nu(nu) {
        }

        void ParseInput(integer iIndex, MBDynParser& HP, const HydroRootElement* pRoot);

        doublereal dGetReducedModulus(const Material& oMat2) const {
            return 1. / ((1. - nu * nu) / E + (1. - oMat2.nu * oMat2.nu) / oMat2.E);
        }

        doublereal dGetReducedModulus() const {
            return E / (1. - nu * nu);
        }

        static Material Rigid() {
            return Material(std::numeric_limits<doublereal>::max(), 0.);
        }

    private:
        doublereal E;
        doublereal nu;
    };

    class FrictionModel
    {
    protected:
        explicit FrictionModel(HydroMesh* pMesh);

    public:
        virtual ~FrictionModel();
        virtual void ParseInput(MBDynParser& HP)=0;
        virtual void GetFrictionForce(const doublereal h, const Vector<doublereal, 2>& U, doublereal p, Vector<doublereal, 2>& tau)=0;
        virtual void GetFrictionForce(const Gradient<0>& h, const Vector<Gradient<0>, 2>& U, const Gradient<0>& p, Vector<Gradient<0>, 2>& tau)=0;
        virtual void AfterPredict(const VectorHandler& X, const VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);
        virtual std::unique_ptr<FrictionModel> Clone() const=0;
        HydroMesh* pGetMesh() const { return pMesh; }

    private:
        HydroMesh* const pMesh;
    };


    class CoulombFriction: public FrictionModel
    {
    public:
        explicit CoulombFriction(HydroMesh* pMesh);
        virtual ~CoulombFriction();
        virtual void ParseInput(MBDynParser& HP);
        virtual void GetFrictionForce(const doublereal h, const Vector<doublereal, 2>& U, doublereal p, Vector<doublereal, 2>& tau);
        virtual void GetFrictionForce(const Gradient<0>& h, const Vector<Gradient<0>, 2>& U, const Gradient<0>& p, Vector<Gradient<0>, 2>& tau);
        virtual std::unique_ptr<FrictionModel> Clone() const;

    private:
        template <typename T>
        void GetFrictionForceTpl(const T& h, const Vector<T, 2>& U, const T& p, Vector<T, 2>& tau);

        doublereal mu;
        doublereal signumDeltaU;
    };

    class LugreFriction: public FrictionModel
    {
    public:
        explicit LugreFriction(HydroMesh* pMesh);
        virtual ~LugreFriction();
        virtual void ParseInput(MBDynParser& HP);
        virtual void GetFrictionForce(const doublereal h, const Vector<doublereal, 2>& U, doublereal p, Vector<doublereal, 2>& tau);
        virtual void GetFrictionForce(const Gradient<0>& h, const Vector<Gradient<0>, 2>& U, const Gradient<0>& p, Vector<Gradient<0>, 2>& tau);
        virtual std::unique_ptr<FrictionModel> Clone() const;
        virtual void AfterConvergence(const VectorHandler& X, const VectorHandler& XP);

    private:
        template <typename T>
        void GetFrictionForceTpl(const T& h, const Vector<T, 2>& U, const T& p, Vector<T, 2>& tau);

        void SaveStictionState(const Vector<doublereal, 2>& z, const Vector<doublereal, 2>& zP);
        void SaveStictionState(const Vector<Gradient<0>, 2>& z, const Vector<Gradient<0>, 2>& zP);

        Matrix<doublereal, 2, 2> Mk, Mk2, invMk2_sigma0, Ms, Ms2, sigma0, sigma1;
        doublereal beta, vs, gamma;
        Vector<doublereal, 2> zPrev, zCurr, zPPrev, zPCurr;
        doublereal tPrev, tCurr;
    };

    class ContactModel
    {
    protected:
        explicit ContactModel(HydroMesh* pMesh);

    public:
        virtual ~ContactModel();

        virtual void ParseInput(MBDynParser& HP);
        virtual bool GetContactPressure(const doublereal h, doublereal& pasp) const=0;
        virtual bool GetContactPressure(const Gradient<0>& h, Gradient<0>& pasp) const=0;
        HydroMesh* pGetMesh() const { return pMesh; }
        const Material& GetMaterial(index_type i) const {
            HYDRO_ASSERT(i >= 0);
            HYDRO_ASSERT(size_t(i) < rgMaterials.size());
            return rgMaterials[i];
        }

    private:
        HydroMesh* const pMesh;
        static const index_type iNumMaterials = 2;
        std::array<Material, iNumMaterials> rgMaterials;
    };

    class GreenwoodTrippCM: public ContactModel {
    public:
        explicit GreenwoodTrippCM(HydroMesh* pMesh);
        virtual ~GreenwoodTrippCM();
        virtual void ParseInput(MBDynParser& HP);
        virtual bool GetContactPressure(const doublereal h, doublereal& pasp) const;
        virtual bool GetContactPressure(const Gradient<0>& h, Gradient<0>& pasp) const;

    private:
        template <typename T>
        bool ContactPressureTpl(const T& h, T& pasp) const;

        doublereal k, sigmaDelta;
        doublereal H0, Hoffset, a0, a1;
    };

    class PenaltyCM: public ContactModel {
    public:
        explicit PenaltyCM(HydroMesh* pMesh, doublereal href);
        virtual ~PenaltyCM();
        virtual void ParseInput(MBDynParser& HP);
        virtual bool GetContactPressure(const doublereal h, doublereal& pasp) const;
        virtual bool GetContactPressure(const Gradient<0>& h, Gradient<0>& pasp) const;

    private:
        template <typename T>
        bool ContactPressureTpl(const T& h, T& pasp) const;

        doublereal a, b, c, h0, h1, href;
    };

    class HydroElement {
    public:
        enum ElementType {
            REYNOLDS_ELEM,
            FRICTION_ELEM,
            COUPLING_ELEM,
            COMPLIANCE_ELEM,
            THERMAL_ELEM
        };

        explicit HydroElement(HydroMesh* pMesh, ElementType eType);
        virtual ~HydroElement();

        ElementType GetElementType() const { return eType; }
        virtual int iGetNumNodes() const=0;
        virtual void SetNode(int iNode, HydroNode* pNode)=0;
        virtual HydroNode* pGetNode(int iNode) const=0;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const;

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode)=0;

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode)=0;

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode)=0;

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode)=0;

        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X,
                                      const VectorHandler& XP);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const=0;
        virtual void Initialize()=0;
        HydroMesh* pGetMesh() const { return pMesh; }
        const HydroFluid* pGetFluid() const { return pFluid; }

        friend inline std::ostream& operator<<(std::ostream& os, const HydroElement& oElem) {
            oElem.Print(os);
            return os;
        }

    protected:
        virtual void Print(std::ostream& os) const;

    private:
        HydroMesh* const pMesh;
        const HydroFluid* const pFluid;
        const ElementType eType;
    };

    class LinFD5Elem: public HydroElement {
    public:

        /*
         *              Node layout
         *
         *           North
         *             |
         *             |
         * West --- Center --- East
         *             |
         *             |
         *           South
         *
         *                      ^ z
         *                      |
         *                      +---> x
         */

        static const int iNodeCenter = 0;
        static const int iNodeWest = 1;
        static const int iNodeEast = 2;
        static const int iNodeSouth = 3;
        static const int iNodeNorth = 4;

        static const int iNodeFlxWest = 0;
        static const int iNodeFlxEast = 1;
        static const int iNodeFlzSouth = 2;
        static const int iNodeFlzNorth = 3;

        explicit LinFD5Elem(HydroMesh* pMesh, ElementType eType);
        virtual ~LinFD5Elem();
        virtual int iGetNumNodes() const;
        virtual void SetNode(int iNode, HydroNode* pNode);
        virtual void SetFluxNode(int iNode, FluxNode* pFluxNode);
        virtual HydroNode* pGetNode(int iNode) const;
        virtual void Initialize();

    protected:
        static const int iNumNodes = 5;
        static const int iNumFluxNodes = 4;
        std::array<HydroNode*, iNumNodes> rgHydroNodes;
        std::array<FluxNode*, iNumFluxNodes> rgFluxNodes;
        std::array<Vector<doublereal, 2>, iNumNodes> x;
        doublereal dx;
        doublereal dz;
        doublereal dA;
    };

    class LinFD4Elem: public HydroElement {
    public:
        explicit LinFD4Elem(HydroMesh* pMesh, ElementType eType);
        virtual ~LinFD4Elem();
        virtual int iGetNumNodes() const;
        virtual void SetNode(int iNode, HydroNode* pNode);
        virtual HydroNode* pGetNode(int iNode) const;
        virtual void Initialize();

        /*
         *              Node layout
         *
         * Node2NW ------- Node1NE
         *    |         ^ z       |
         *    |         |         |
         *    |         +---> x   |
         * Node3SW ------- Node4SE
         */

        static const int iNode1NE = 0; /**< north east */
        static const int iNode2NW = 1; /**< north west */
        static const int iNode3SW = 2; /**< south west */
        static const int iNode4SE = 3; /**< south east */

    protected:
        static const int iNumNodes = 4;
        std::array<HydroNode*, iNumNodes> rgHydroNodes;
        doublereal dx, dz, dA;
    };

    class LinFD5ReynoldsElem: public LinFD5Elem {
    public:
        explicit LinFD5ReynoldsElem(HydroMesh* pMesh);
        virtual ~LinFD5ReynoldsElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

    private:
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);
    private:
        LocalDofMap oDofMap;

#if CREATE_PROFILE == 1
        enum { PROF_RES = 0, PROF_JAC = 1 };
        static struct ProfileData {
            doublereal dtAss[2];
            HydroElement* pLastElem;
        } profile;
#endif
    };

#if CREATE_PROFILE == 1
    LinFD5ReynoldsElem::ProfileData LinFD5ReynoldsElem::profile;
#endif

    class LinFD5CouplingElem: public LinFD5Elem {
    public:
        explicit LinFD5CouplingElem(HydroMesh* pMesh);
        virtual ~LinFD5CouplingElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

    private:
        LocalDofMap oDofMap;
    };

    class LinFD4FrictionElem: public LinFD4Elem {
    public:
        explicit LinFD4FrictionElem(HydroMesh* pMesh);
        virtual ~LinFD4FrictionElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T> inline
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void Initialize();

    private:
        template <typename T> inline
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type, const std::array<Vector<doublereal, 2>, iNumNodes>& Ui, doublereal dTau_xy, doublereal dTau_yz) const;
        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type, const std::array<Vector<Gradient<0>, 2>, iNumNodes>& Ui, const Gradient<0>& dTau_xy, const Gradient<0>& dTau_yz) const;
        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type, doublereal dPf) const;
        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type, const Gradient<0>& dPf) const;
        static inline void SetStress(HydroNode* pNode, doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h);
        static inline void SetStress(HydroNode*, const Gradient<0>&, const Gradient<0>&, const Gradient<0>&, const Gradient<0>&);

        LocalDofMap oDofMap;
        Vector<doublereal, 2> xc;
        Vector<doublereal, 3> vc;
        Matrix<doublereal, 3, 3> Rtc;
        doublereal dScaleEnergy;

#if CREATE_PROFILE == 1
        enum { PROF_RES = 0, PROF_JAC = 1 };
        static struct ProfileData {
            doublereal dtAss[2];
            HydroElement* pLastElem;
        } profile;
#endif
    };

#if CREATE_PROFILE == 1
    LinFD4FrictionElem::ProfileData LinFD4FrictionElem::profile;
#endif

    class LinFD4MassFlowZ: public LinFD4Elem {
    public:
        explicit LinFD4MassFlowZ(HydroMesh* pMesh);
        virtual ~LinFD4MassFlowZ();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T> inline
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        virtual void Initialize();

        void SetFluxNode(int iNode, FluxNode* pFluxNode);

        static const index_type iFNodeWest = 0;
        static const index_type iFNodeEast = 1;

    private:
        static const index_type iNumFluxNodes = 2;
        std::array<FluxNode*, iNumFluxNodes> rgFluxNodes;
        LocalDofMap oDofMap;
    };

    class LinFD5ComprReynoldsElem: public LinFD5Elem {
    public:
        explicit LinFD5ComprReynoldsElem(HydroMesh* pMesh);
        virtual ~LinFD5ComprReynoldsElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);
    private:
        inline void SetMaxTimeStep(const std::array<doublereal, iNumFluxNodes>& w) const;
        inline void SetMaxTimeStep(const std::array<Gradient<0>, iNumFluxNodes>&) const;

        static const index_type iNumDofMax = HydroActiveComprNode::iNumDofMax;
        LocalDofMap oDofMap;
    };

    class LinFD5ThermalElem: public LinFD5Elem {
    public:
        explicit LinFD5ThermalElem(HydroMesh* pMesh);
        virtual ~LinFD5ThermalElem();

        virtual void SetNode(int iNode, HydroNode* pNode);
        virtual void SetFluxNode(int iNode, FluxNode* pFluxNode);

        virtual int iGetNumNodesTherm() const;
        virtual ThermoHydrNode* pGetNodeTherm(int iNode) const;
        virtual void Initialize();

    protected:
        template <typename G>
        inline void EnergyBalance(G& f, G& Qdot0, G& Qdoth, LocalDofMap* pDofMap, doublereal dCoef, grad::FunctionCall func) const;
        inline void SetMaxTimeStep(const std::array<doublereal, 2>& wxi, const std::array<doublereal, 2>& wzi) const;
        static inline void UpwindWeight(const std::array<doublereal, 2>& q, std::array<doublereal, 2>& alpha);

        std::array<ThermoHydrNode*, iNumNodes> rgThermNodes;
        doublereal dScale;
        LocalDofMap oDofMap;
    };

    class LinFD5ThermalElemImp: public LinFD5ThermalElem {
    public:
        explicit LinFD5ThermalElemImp(HydroMesh* pMesh, bool bDoInitAss);
        virtual ~LinFD5ThermalElemImp();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

    private:
        template <typename G>
        void UnivAssRes(GradientAssVec<G>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<G>& XCurr,
                        enum FunctionCall func);
        const bool bDoInitAss;
    };

    class LinFD5ThermalCouplingElem: public LinFD5ThermalElem {
    public:
        explicit LinFD5ThermalCouplingElem(HydroMesh* pMesh,
                                           bool bDoInitAss);
        virtual ~LinFD5ThermalCouplingElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows,
                                  integer* piNumCols,
                                  grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void Initialize();

    private:
        LocalDofMap oDofMap;
        const ThermalInletNode* pInletNode;
        const bool bDoInitAss;
    };

    class QuadFeIso9Elem: public HydroElement {
    public:
        class IntegrationRule {
        public:
            inline explicit
            IntegrationRule(ElementType eElemType,
                            index_type iGaussMin = 2,
                            index_type iGaussMax = 2,
                            index_type iGaussStep = 1);

            void ParseInput(MBDynParser& HP, const HydroMesh* pMesh);

            index_type iGetGaussFirst() const {
                HYDRO_ASSERT(bInvariant());
                return iGaussMin;
            }
            index_type iGetGaussLast() const {
                HYDRO_ASSERT(bInvariant());
                return iGaussMax;
            }
            index_type iGetGaussNext(index_type iGaussCurr) const {
                HYDRO_ASSERT(bInvariant());
                return iGaussCurr + iGaussStep;
            }
            index_type iGetGaussCount() const {
                HYDRO_ASSERT(bInvariant());
                return (iGaussMax - iGaussMin) / iGaussStep + 1;
            }
        private:
#if HYDRO_DEBUG > 0
            bool bInvariant() const;
#endif
            index_type iGaussMin;
            index_type iGaussMax;
            index_type iGaussStep;
            const ElementType eElemType;
        };

        /** Node layout
         *
         *   2----5----1
         *   |    |    |
         *   6----9----8
         *   |    |    |
         *   3----7----4
         *
         *      ^ s
         *      |
         *      +---> r
         */
        explicit QuadFeIso9Elem(HydroMesh* pMesh,
                                const IntegrationRule& oIntegRule,
                                ElementType eType);
        virtual ~QuadFeIso9Elem();
        virtual int iGetNumNodes() const;
        virtual void SetNode(int iNode, HydroNode* pNode);
        virtual HydroNode* pGetNode(int iNode) const;
        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X,
                                      const VectorHandler& XP);
        static const index_type iNumNodes = 9;

    protected:
        struct GaussPointData {
            Vector<doublereal, iNumNodes> N;
            Matrix<doublereal, 2, iNumNodes> B;
            doublereal detJ;
        };

        inline index_type iGetNumGaussPoints(index_type iIntegRule) const;

        inline index_type iGetNumIntegrationRules() const;

        index_type iSelectIntegrationRule(const Vector<doublereal, iNumNodes>& pe,
                                          const Vector<doublereal, iNumNodes>* paspe = nullptr);

        index_type iSelectIntegrationRule(const Vector<Gradient<0>, iNumNodes>& pe,
                                          const Vector<Gradient<0>, iNumNodes>* paspe = nullptr);

        index_type iGetCurrentIntegrationRule() const {
            HYDRO_ASSERT(iCurrIntegRule >= 0);
            HYDRO_ASSERT(size_t(iCurrIntegRule) < rgGauss.size());

            return iCurrIntegRule;
        }

        inline index_type iGetGaussPointSize() const;
        inline index_type iGetGaussPointSize1D() const;
        inline index_type iGetGaussPointIndex(index_type iGaussR, index_type iGaussS, index_type iIntegRule) const;
        inline index_type iGetGaussPointIndex1D(index_type iGaussR, index_type iIntegRule) const;
        inline doublereal dGetGaussWeight(index_type iGaussPoint, index_type iIntegRule) const;

        inline doublereal dGetGaussPos(index_type iGaussPoint, index_type iIntegRule) const;

        void NodePositionMatrix(Matrix<doublereal, iNumNodes, 2>& xe) const;

        void PressureInterpolMatrix(Vector<doublereal, iNumNodes>& N,
                                    doublereal r,
                                    doublereal s) const;

        void PressureInterpolMatrixDer(Matrix<doublereal, 2, iNumNodes>& dN_drs,
                                       doublereal r,
                                       doublereal s) const;

        void PressureGradInterpolMatrix(Matrix<doublereal, 2, iNumNodes>& B,
                                        doublereal& detJ,
                                        const Matrix<doublereal, iNumNodes, 2>& xe,
                                        doublereal r,
                                        doublereal s) const;

        std::array<HydroNode*, iNumNodes> rgNodes;

        static const index_type iNumNodesPerGroup = 4;
        static const index_type iNumNodeGroups = 4;

        struct NodeGroup {
            doublereal r, s;
            index_type nodes[iNumNodesPerGroup];
        };

        static const NodeGroup* pFindNodeGroup(doublereal r, doublereal s);

    private:
        static const NodeGroup rgNodeGroups[iNumNodeGroups];
        std::vector<GaussData> rgGauss;
        integer iCurrIntegRule;
    };

    class QuadFeIso9ReynoldsElem: public QuadFeIso9Elem {
    public:
        explicit QuadFeIso9ReynoldsElem(HydroMesh* pMesh,
                                        const IntegrationRule& oIntegRule);
        virtual ~QuadFeIso9ReynoldsElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkMat,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void Initialize();

    private:
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        std::vector<GaussPointData> rgGaussPntDat;
        LocalDofMap oDofMap;
        doublereal dA;
    };

    class QuadFeIso9FrictionElem: public QuadFeIso9Elem {
    public:
        explicit QuadFeIso9FrictionElem(HydroMesh* pMesh,
                                        const IntegrationRule& oIntegRule);

        virtual ~QuadFeIso9FrictionElem();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkMat,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkMat,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void Initialize();

    protected:
        struct GaussPointDataF: GaussPointData {
            Vector<doublereal, 2> xc;
            Vector<doublereal, 3> vc;
            Matrix<doublereal, 3, 3> Rtc;
        };

    private:
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkMat,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type,
                                    const Vector<doublereal, 2>& U,
                                    doublereal dTau_xy,
                                    doublereal dTau_yz) const;

        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type,
                                    const Vector<Gradient<0>, 2>& U,
                                    const Gradient<0>& dTau_xy,
                                    const Gradient<0>& dTau_yz) const;

        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type,
                                    doublereal dPf) const;
        inline void AddFrictionLoss(HydroRootBase::FrictionLossType type,
                                    const Gradient<0>& dPf) const;

        void SetStress(doublereal r,
                       doublereal s,
                       doublereal tau_xy_0,
                       doublereal tau_yz_0,
                       doublereal tau_xy_h,
                       doublereal tau_yz_h) const;

        void SetStress(doublereal r,
                       doublereal s,
                       const Gradient<0>& tau_xy_0,
                       const Gradient<0>& tau_yz_0,
                       const Gradient<0>& tau_xy_h,
                       const Gradient<0>& tau_yz_h) const {}

        std::vector<GaussPointDataF> rgGaussPntDat;
        LocalDofMap oDofMap;
    };

    class QuadFeIso9MassFlowZ: public QuadFeIso9Elem {
    public:
        explicit QuadFeIso9MassFlowZ(HydroMesh* pMesh,
                                     const IntegrationRule& oIntegRule,
                                     doublereal sref);
        virtual ~QuadFeIso9MassFlowZ();

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      enum GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkMat,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void Initialize();

    private:
        struct GaussPointDataMassFlow: GaussPointData {
            GaussPointDataMassFlow(): detJr{0.} {};
            doublereal detJr;
        };

        std::vector<GaussPointDataMassFlow> rgGaussPntDat;
        LocalDofMap oDofMap;
        const doublereal sref;
        static const index_type iNumNodesOutletBound = 3;
        const index_type* const pNodeIndexOutletBound;
        static const index_type rgNodeIndexOutletBound[2][iNumNodesOutletBound];
    };

    class HydroMesh {
    public:
        explicit HydroMesh(HydroRootElement* pParent);
        virtual ~HydroMesh();
        /**
         *  @param p returns the pressure at the node by taking into account
         *               the specific boundary conditions of the selected elements
         *               (e.g. boundary conditions according to Guembel)
         */
        virtual void GetPressure(const HydroNode* pNode, doublereal& p, LocalDofMap* =0, doublereal=0.) const;
        virtual void GetPressure(const HydroNode* pNode, Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef=0.) const;
        virtual void ParseInput(DataManager* pDM, MBDynParser& HP)=0;
        virtual integer iGetNumNodes() const=0;
        virtual integer iGetNumElements() const=0;
        virtual integer iGetNumBounaryConditions() const;
        virtual void Generate()=0;
        virtual CylindricalBearing* pGetGeometry()const;
        virtual ComplianceModel* pGetComplianceModel() const;
        virtual std::ostream& PrintLogFile(std::ostream& os) const;
        virtual std::ostream& Output(std::ostream& os) const=0;
        HydroRootElement* pGetParent() const { return pParent; }
        const ThermWallBoundCond* pGetThermWallBoundCond() const {
            return pThermWallBoundCond.get();
        }
        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr, doublereal dCoef,
               enum FunctionCall func,
               LocalDofMap* pDofMap);
        inline doublereal dGetMaxPressureGradient() const;

    protected:
        typedef std::vector<std::unique_ptr<FluidStateBoundaryCond> > BoundaryCondVector;
        typedef std::vector<std::unique_ptr<LubricationGroove> > LubrGrooveVector;
        typedef std::vector<std::unique_ptr<PressureCouplingCond> > CouplingVector;

        LubricationGroove* pFindGroove(const Vector<doublereal, 2>& x, Node2D::NodeType eNodeType, integer iNodeId);
        PressureCouplingCond* pFindCouplingCond(const Vector<doublereal, 2>& x, integer iNodeId);

        void ParseGeometry(DataManager* pDM, MBDynParser& HP);
        void ParseBoundaryCond(DataManager* pDM, MBDynParser& HP);
        void ParsePressureCouplingCond(DataManager* pDM, MBDynParser& HP);
        void ParseContactModel(DataManager* pDM, MBDynParser& HP);
        void ParseLubricationGrooves(DataManager* pDM, MBDynParser& HP);
        void ParseThermWallBoundCond(DataManager* pDM, MBDynParser& HP);
        void ParseComplianceModel(DataManager* pDM, MBDynParser& HP);
        void SetThermalModel(const HydroFluid* pFluid);
        void GenerateBoundaryConditions();
        void GenerateComplianceModel();
        void GenerateMovingLubricationGrooves();
        BoundaryCondVector rgBoundaryCond;
        LubrGrooveVector rgGrooves;
        CouplingVector rgCouplingCond;
        std::unique_ptr<CylindricalBearing> pGeometry;
        std::unique_ptr<ContactModel> pContact;
        std::unique_ptr<FrictionModel> pFriction;
        std::unique_ptr<ThermWallBoundCond> pThermWallBoundCond;
        util::variable_ptr<ComplianceModel> pCompliance;
        std::string ParseFileName(MBDynParser& HP);

        struct CouplingCondition {
            CouplingCondition(PressureNode* pHydroNode = nullptr,
                              ThermalNode* pThermalNode = nullptr)
                :pExtHydroNode(pHydroNode),
                 pExtThermNode(pThermalNode) {
            }
            PressureNode* pExtHydroNode;
            ThermalNode* pExtThermNode;
        };

        static const index_type iNumCouplingAxial = 2;
        std::array<CouplingCondition, iNumCouplingAxial> rgOutletAxial;
        bool bUseOutletAxial;
        const bool bThermalModel;

    private:
        HydroRootElement* pParent;
    };

/**
 *      \brief Rectangular mesh with linear finite differences
 */
    class LinFDMesh: public HydroMesh {
    public:
        explicit LinFDMesh(HydroRootElement* pParent);
        virtual ~LinFDMesh();
        virtual void ParseInput(DataManager* pDM, MBDynParser& HP);
        virtual integer iGetNumNodes() const;
        virtual integer iGetNumElements() const;
        virtual void Generate();
        virtual std::ostream& Output(std::ostream& os) const;

    private:
        inline Vector<doublereal, 2> GetNodePosition(integer i, integer j) const;
        inline integer iGetNodeIndexHydro(integer i, integer j) const;
        inline integer iGetNodeIndexTherm(integer i, integer j) const;
        inline integer iGetNodeIndexFluxX(integer i, integer j) const;
        inline integer iGetNodeIndexFluxZ(integer i, integer j) const;

    private:
        integer M, N;
        std::vector<doublereal> x, z;
        enum ElementType {
            CENT_DIFF_5
        } eElemType;
    };

    class QuadFeIso9Mesh: public HydroMesh {
    public:
        explicit QuadFeIso9Mesh(HydroRootElement* pParent);
        virtual ~QuadFeIso9Mesh();
        virtual void ParseInput(DataManager* pDM, MBDynParser& HP);
        virtual integer iGetNumNodes() const;
        virtual integer iGetNumElements() const;
        virtual void Generate();
        virtual std::ostream& Output(std::ostream& os) const;

    private:
        enum MeshGeometry
        {
            PLAIN_BEARING,
            HELICAL_GROOVE_PUMP
        };

        index_type iGetNodeIndex(index_type i, index_type j) const;
        Vector<doublereal, 2> GetNodePosition(index_type i, index_type j) const;
        const FluidStateBoundaryCond* pFindBoundaryCond(index_type i, index_type j) const;
        Vector<doublereal, DYNAMIC_SIZE> x, z;
        QuadFeIso9Elem::IntegrationRule oIntegRuleReynolds, oIntegRuleFriction;
        doublereal dSkewMesh;
        MeshGeometry eMeshGeometry;
    };

    class PressureElement {
    public:
        enum NodePosition {
            NODE_NE = 0,
            NODE_NW = 1,
            NODE_SW = 2,
            NODE_SE = 3
        };

        PressureElement() {
            std::fill(std::begin(rgNodes), std::end(rgNodes), nullptr);
        }

        index_type iGetNumNodes() const {
            return rgNodes.size();
        }

        void SetNode(index_type iNode, const HydroNode* pNode) {
            HYDRO_ASSERT(iNode >= 0);
            HYDRO_ASSERT(iNode < iGetNumNodes());

            rgNodes[iNode] = pNode;
        }

        const HydroNode* pGetNode(index_type iNode) const {
            HYDRO_ASSERT(iNode >= 0);
            HYDRO_ASSERT(iNode < iGetNumNodes());
            HYDRO_ASSERT(rgNodes[iNode] != nullptr);

            return rgNodes[iNode];
        }

    private:
        std::array<const HydroNode*, 4> rgNodes;
    };

    struct ComplianceMatrixCommon {
        static const grad::index_type MAX_NUM_MATRICES = 3;
        typedef Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> MatrixType;

        enum MatrixKind {
            MAT_FULL,
            MAT_MODAL,
            MAT_INVALID
        };

        enum MeshLocation {
            LOC_MESH_FIXED,
            LOC_MESH_MOVING
        };

        class MatrixArray {
        public:
            MatrixArray(MatrixType* pC, MatrixType* pD, MatrixType* pE)
                :pC(pC), pD(pD), pE(pE), eMatType(ComplianceMatrixCommon::MAT_FULL) {
            }

            MatrixArray(MatrixType* pRPhiK, MatrixType* pPhin)
                :pRPhiK(pRPhiK), pPhin(pPhin), eMatType(ComplianceMatrixCommon::MAT_MODAL) {
                rgMat[2] = nullptr;
            }

            MatrixKind GetMatrixType() const {
                return eMatType;
            }

            MatrixType* C() const {
                HYDRO_ASSERT(pC != nullptr);
                HYDRO_ASSERT(eMatType == ComplianceMatrixCommon::MAT_FULL);

                return pC;
            }

            MatrixType* D() const {
                HYDRO_ASSERT(pD != nullptr);
                HYDRO_ASSERT(eMatType == ComplianceMatrixCommon::MAT_FULL);

                return pD;
            }

            MatrixType* E() const {
                HYDRO_ASSERT(pE != nullptr);
                HYDRO_ASSERT(eMatType == ComplianceMatrixCommon::MAT_FULL);

                return pE;
            }

            MatrixType* RPhiK() const {
                HYDRO_ASSERT(pRPhiK != nullptr);
                HYDRO_ASSERT(rgMat[2] == nullptr);
                HYDRO_ASSERT(eMatType == ComplianceMatrixCommon::MAT_MODAL);

                return pRPhiK;
            }

            MatrixType* Phin() const {
                HYDRO_ASSERT(pPhin != nullptr);
                HYDRO_ASSERT(rgMat[2] == nullptr);
                HYDRO_ASSERT(eMatType == ComplianceMatrixCommon::MAT_MODAL);

                return pPhin;
            }

            MatrixType* operator[] (index_type i) const {
                HYDRO_ASSERT(i >= 0);
                HYDRO_ASSERT(i < MAX_NUM_MATRICES);
                HYDRO_ASSERT(rgMat[i] != nullptr);

                return rgMat[i];
            }

        private:
            union {
                struct {
                    MatrixType* pC;
                    MatrixType* pD;
                    MatrixType* pE;
                };

                struct {
                    MatrixType* pRPhiK;
                    MatrixType* pPhin;
                };

                MatrixType* rgMat[MAX_NUM_MATRICES];
            };

            ComplianceMatrixCommon::MatrixKind eMatType;
        };

        struct GridIndex {
            GridIndex(index_type ix, index_type iz)
                :ix(ix), iz(iz) {
            }

            bool operator==(const GridIndex& oGridIndex) const {
                return ix == oGridIndex.ix && iz == oGridIndex.iz;
            }

            index_type ix;
            index_type iz;
        };

        class MatrixData {
        public:
            template <typename... ARGS>
            MatrixData(MeshLocation eMeshLocation,
                       const Modal* pModalJoint,
                       ARGS... pMatrix)
                :eMeshLocation(eMeshLocation),
                 pModalJoint(pModalJoint),
                 rgMatrices(pMatrix...) {
            }

            const MeshLocation eMeshLocation;
            const Modal* const pModalJoint;
            MatrixArray rgMatrices;
            std::vector<doublereal> rgGridX, rgGridZ;
            std::vector<index_type> rgMatIdx;
            std::vector<GridIndex> rgActGridIdx;
            doublereal dBearingDiameter = 0.;
            doublereal dBearingWidth = 0.;
            doublereal dRefPressure = 0.;            
        };
    };

    class ComplianceMatrixFileParser: public ComplianceMatrixCommon {
        typedef std::vector<HydroUpdatedNode*> NodesContainer;
        typedef std::vector<PressureElement> ElementContainer;

    public:
        explicit ComplianceMatrixFileParser(MatrixData& oMatData)
            :oMatData(oMatData) {
        }

        void Parse(const std::string& strFileName,
                   const HydroMesh* pMesh,
                   const ElementContainer& rgElements,
                   const NodesContainer& rgNodes,
                   doublereal dPressScale);

    private:
        struct NodeRec {
            NodeRec(index_type iComplianceIndex, doublereal x, doublereal z)
                :iComplianceIndex(iComplianceIndex), x(x), z(z) {
            }

            index_type iComplianceIndex;
            doublereal x, z;

            bool operator <(const NodeRec& oRec) const {
                if (x < oRec.x) {
                    return true;
                } else if (x == oRec.x) {
                    return z < oRec.z;
                } else {
                    return false;
                }
            }
        };

        struct InpNodeRec {
            index_type iNode;
            index_type iDof;
            index_type iPosX;
            index_type iPosZ;
        };

        MatrixData& oMatData;
    };

    class ComplianceMatrix: public ComplianceMatrixCommon {
    public:
        typedef std::vector<PressureElement> ElementContainer;
        typedef std::vector<HydroUpdatedNode*> NodesContainer;

    protected:
        ComplianceMatrix();

    public:
        virtual ~ComplianceMatrix();
        virtual void AddCompliance(MatrixData& oMatData,
                                   const HydroMesh* pMesh,
                                   const ElementContainer& rgElements,
                                   const NodesContainer& rgNodes,
                                   doublereal dPressScale) const=0;
    };

    class ElasticHalfSpace: public ComplianceMatrix {
    public:
        explicit ElasticHalfSpace(const Material& oMaterial1,
                                  const Material& oMaterial2);
        virtual ~ElasticHalfSpace();
        virtual void AddCompliance(MatrixData& oMatData,
                                   const HydroMesh* pMesh,
                                   const ElementContainer& rgElements,
                                   const NodesContainer& rgNodes,
                                   doublereal dPressScale) const;
    private:
        doublereal Ered;
    };

    class ComplianceFromFile: public ComplianceMatrix {
    public:
        explicit ComplianceFromFile(const std::string& strFileName);
        virtual ~ComplianceFromFile();
        virtual void AddCompliance(MatrixData& oMatData,
                                   const HydroMesh* pMesh,
                                   const ElementContainer& rgElements,
                                   const NodesContainer& rgNodes,
                                   doublereal dPressScale) const;
    private:
        std::string strFileName;
    };

    class ComplianceModel: public HydroElement, public HydroDofOwner {
    public:
        enum Type {
            COMP_MOD_UNKNOWN = -1,
            COMP_MOD_MODAL = 0,
            COMP_MOD_NODAL,
            COMP_MOD_NODAL_DOUBLE
        };

        explicit ComplianceModel(HydroMesh* pMesh,
                                 doublereal dDefScale,
                                 doublereal dPressScale);
        virtual ~ComplianceModel();
        virtual int iGetNumNodes() const;
        virtual void SetNode(int iNode, HydroNode* pNode);
        virtual HydroNode* pGetNode(int iNode) const;
        virtual integer iGetFirstIndex(grad::FunctionCall eFunc) const;
        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void AfterConvergence(const VectorHandler& X,
                                      const VectorHandler& XP);
        virtual void Initialize();

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               FunctionCall func)=0;

        virtual void
        GetRadialDeformation(doublereal& w,
                             doublereal& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const=0;

        virtual void
        GetRadialDeformation(Gradient<0>& w,
                             Gradient<0>& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const=0;

        virtual void
        GetRadialDeformation1(doublereal& w1,
                              const HydroUpdatedNode* pNode) const;
        virtual void
        GetRadialDeformation2(doublereal& w2,
                              const HydroUpdatedNode* pNode) const;

        void EnableInitAss(bool bEnable) {
            bDoInitAss = bEnable;
        }
    protected:
        typedef std::vector<PressureElement> ElementContainer;
        typedef std::vector<HydroUpdatedNode*> NodesContainer;

        ElementContainer rgElements;
        NodesContainer rgNodes;
        const doublereal dDefScale, dPressScale;
        bool bDoInitAss;
    };

    class ComplianceModelNodal: public ComplianceModel {
    public:
        typedef std::unique_ptr<ComplianceMatrix> ComplianceMatrixPtr;
        typedef std::array<ComplianceMatrixPtr, 2> ComplianceMatrixArray;

        explicit ComplianceModelNodal(HydroMesh* pMesh,
                                      const Modal* pModalJoint,
                                      doublereal dDefScale,
                                      doublereal dPressScale,
                                      ComplianceMatrixArray&& rgMatrices);
        virtual ~ComplianceModelNodal();

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        virtual void Initialize();

        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual unsigned int iGetNumDof(void) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               FunctionCall func);

        virtual void
        GetRadialDeformation(doublereal& w,
                             doublereal& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;

        virtual void
        GetRadialDeformation(Gradient<0>& w,
                             Gradient<0>& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;

    private:
        enum DofMapIndex {
            DOF_PRESS = 0,
            DOF_DEF,
            DOF_LAST
        };
        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        virtual void Print(std::ostream& os) const;

        index_type iNumNodes, iNumModes;
        Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> C, D, E;
        Vector<doublereal, DYNAMIC_SIZE> w, dw_dt;
        const Modal* const pModalJoint;
        std::array<LocalDofMap, DOF_LAST> rgDofMap;
        ComplianceMatrixArray rgMatrices;
        grad::FunctionCall eCurrFunc;
    };

    class ComplianceModelNodalDouble: public ComplianceModel {
    public:
        enum DEhdBodyIdx {
            DEHD_BODY_FIXED = 0,
            DEHD_BODY_MOVING,
            DEHD_BODY_LAST
        };

        enum DEhdDeformationIdx {
            DEHD_DEF_TOTAL = DEHD_BODY_FIXED,
            DEHD_DEF_MOVING = DEHD_BODY_MOVING,
            DEHD_DEF_MOVING_INTERP
        };

        enum DEhdInterpolOption {
            INT_AXIAL_SMALL_DISP,
            INT_AXIAL_LARGE_DISP,
            INT_AXIAL_EXTRAPOLATE
        };
        
        typedef std::unique_ptr<ComplianceMatrix> ComplianceMatrixPtr;
        typedef std::array<ComplianceMatrixPtr, DEHD_BODY_LAST> ComplianceMatrixArray;
        typedef std::array<const Modal*, DEHD_BODY_LAST> ModalJointArray;

        explicit ComplianceModelNodalDouble(HydroMesh* pMesh,
                                            const ModalJointArray& rgModalJoints,
                                            doublereal dDefScale,
                                            doublereal dPressScale,
                                            ComplianceMatrixArray&& rgMatrices,
                                            const CylindricalBearing& oGeometry,
                                            DEhdInterpolOption eInterpolOption);
        virtual ~ComplianceModelNodalDouble();

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        virtual void Initialize();

        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual unsigned int iGetNumDof(void) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               FunctionCall func);

        virtual void
        GetRadialDeformation(doublereal& w,
                             doublereal& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;

        virtual void
        GetRadialDeformation(Gradient<0>& w,
                             Gradient<0>& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;
        virtual void
        GetRadialDeformation1(doublereal& w1,
                              const HydroUpdatedNode* pNode) const;
        virtual void
        GetRadialDeformation2(doublereal& w2,
                              const HydroUpdatedNode* pNode) const;
    private:
        enum DEhdFieldType {
            FT_TOTAL_PRESS,
            FT_DEF_MOVING
        };

        enum DEhdDofMapType {
            DOF_GLOBAL1 = 0,
            DOF_GLOBAL2,
            DOF_LOCAL,
            DOF_LAST
        };

        static constexpr index_type POLYORDER = 9;
        static constexpr index_type GRIDINTERP = 16;

        struct PolyData {
            template <typename MatExpr>
            PolyData(const MatExpr& pinvA,
                     const std::array<index_type, GRIDINTERP>& Cidx)
                :pinvA(pinvA),
                 Cidx(Cidx) {
            }

            Matrix<doublereal, POLYORDER, GRIDINTERP> pinvA;
            std::array<index_type, GRIDINTERP> Cidx;
        };

        template <typename G>
        void
        GetRadialDeformationTpl(G& w,
                                G& dw_dt,
                                doublereal dCoef,
                                FunctionCall func,
                                LocalDofMap* pDofMap,
                                const HydroUpdatedNode* pNode) const;

        void
        GetRadialDeformation(doublereal& wi,
                             doublereal& dwi_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             DEhdDeformationIdx eDefIndex,
                             index_type iCompIndex) const;

        void
        GetRadialDeformation(Gradient<0>& wi,
                             Gradient<0>& dwi_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             DEhdDeformationIdx eDefIndex,
                             index_type iCompIndex) const;

        template <DEhdFieldType eField, DEhdBodyIdx eMshSrc, DEhdBodyIdx eMshDst, typename T>
        void Interpolate(const index_type i,
                         const index_type j,
                         const Vector<T, 2>& dxm_g,
                         const doublereal dScale,
                         T& fij_f,
                         const doublereal dCoef,
                         const FunctionCall func,
                         LocalDofMap* const pDofMap) const;

        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        virtual void Print(std::ostream& os) const;

        inline void UpdateDefMovingInterp(index_type iCompIndex, doublereal wmi, doublereal dScale);
        void UpdateDefMovingInterp(index_type, const Gradient<0>&, doublereal){}

    private:
        std::array<index_type, DEHD_BODY_LAST> rgNumNodes;
        std::array<Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE>, DEHD_BODY_LAST> C, D, E;
        std::array<Vector<doublereal, DYNAMIC_SIZE>, DEHD_DEF_MOVING_INTERP + 1> w;
        std::array<Vector<doublereal, DYNAMIC_SIZE>, DEHD_DEF_MOVING + 1> dw_dt;
        std::array<std::vector<doublereal>, DEHD_BODY_LAST> xi, zi;
        std::array<std::vector<index_type>, DEHD_BODY_LAST> rgMatIdx;
        std::array<std::vector<ComplianceMatrix::GridIndex>, DEHD_BODY_LAST> rgActGridIdx;
        std::array<std::vector<PolyData>, DEHD_BODY_LAST> rgPolyData;
        const ModalJointArray rgModalJoints;
        doublereal dPressDofScale;
        std::array<LocalDofMap, DOF_LAST> rgDofMap;
        const doublereal dMeshRadius, dMinDistance_2;
        ComplianceMatrixArray rgMatrices;
        grad::FunctionCall eCurrFunc;
        const DEhdInterpolOption eInterpolOption;
        doublereal dAxialThreshold;
        static const std::array<integer, POLYORDER> px, pz;
        static const std::array<index_type, GRIDINTERP> xg, zg;
        static const index_type min_xg;
        static const index_type max_xg;
        static const index_type min_zg;
        static const index_type max_zg;
    };

    class ComplianceModelModal: public ComplianceModel {
    public:
        explicit ComplianceModelModal(HydroMesh* pMesh,
                                      doublereal dDefScale,
                                      doublereal dPressScale,
                                      const std::string& strFileName);
        virtual ~ComplianceModelModal();

        template <typename T>
        void AssRes(GradientAssVec<T>& WorkVec,
                    doublereal dCoef,
                    const GradientVectorHandler<T>& XCurr,
                    const GradientVectorHandler<T>& XPrimeCurr,
                    enum FunctionCall func);

        virtual void
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        virtual void
        AssJac(SparseSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               GradientAssVecBase::mode_t mode);

        template <typename T>
        void InitialAssRes(GradientAssVec<T>& WorkVec,
                           const GradientVectorHandler<T>& XCurr,
                           enum FunctionCall func);

        virtual void
        InitialAssRes(SubVectorHandler& WorkVec,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void
        InitialAssJac(SparseSubMatrixHandler& WorkMat,
                      const VectorHandler& XCurr,
                      GradientAssVecBase::mode_t mode);

        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const;

        virtual void Initialize();

        virtual void SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr);
        virtual unsigned int iGetNumDof(void) const;
        virtual unsigned int iGetInitialNumDof(void) const;
        virtual integer iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream&
        DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual std::ostream&
        DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;

        virtual void
        Update(const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr,
               doublereal dCoef,
               FunctionCall func);

        virtual void
        GetRadialDeformation(doublereal& w,
                             doublereal& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;

        virtual void
        GetRadialDeformation(Gradient<0>& w,
                             Gradient<0>& dw_dt,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap,
                             const HydroUpdatedNode* pNode) const;

    private:
        inline void
        GetModalDeformation(index_type iMode,
                            doublereal& qi,
                            doublereal dCoef,
                            FunctionCall func,
                            LocalDofMap* pDofMap) const;

        inline void
        GetModalDeformation(index_type iMode,
                            Gradient<0>& qi,
                            doublereal dCoef,
                            FunctionCall func,
                            LocalDofMap* pDofMap) const;

        inline void
        GetModalDeformationDer(index_type iMode,
                               doublereal& dqi_dt,
                               doublereal dCoef,
                               FunctionCall func,
                               LocalDofMap* pDofMap) const;

        inline void
        GetModalDeformationDer(index_type iMode,
                               Gradient<0>& dqi_dt,
                               doublereal dCoef,
                               FunctionCall func,
                               LocalDofMap* pDofMap) const;

        template <typename T>
        inline void
        GetRadialDeformationTpl(T& wi,
                                T& dwi_dt,
                                doublereal dCoef,
                                FunctionCall func,
                                LocalDofMap* pDofMap,
                                const HydroUpdatedNode* pNode) const;

        index_type iGetNumModes() const {
            return iNumModes;
        }

        template <typename T>
        void UnivAssRes(GradientAssVec<T>& WorkVec,
                        doublereal dCoef,
                        const GradientVectorHandler<T>& XCurr,
                        enum FunctionCall func);

        virtual void Print(std::ostream& os) const;

        Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> Phin, RPhiK;
        Vector<doublereal, DYNAMIC_SIZE> q, dq_dt;
        index_type iNumModes;
        LocalDofMap oDofModal, oDofPressure;
        const std::string strFileName;
        grad::FunctionCall eCurrFunc;
    };


    class HydroRootElement: virtual public Elem, public UserDefinedElem, public HydroRootBase
    {
    public:
        enum OutputFlags {
            OUTPUT_NOTHING                = 0x000000,
            OUTPUT_PRESSURE               = 0x000001, // 1
            OUTPUT_CONT_PRESSURE          = 0x000002, // 2
            OUTPUT_DENSITY                = 0x000004, // 3
            OUTPUT_CLEARANCE              = 0x000008, // 4
            OUTPUT_CLEARANCE_DER          = 0x000010, // 5
            OUTPUT_VELOCITY               = 0x000020, // 6
            OUTPUT_STRESS                 = 0x000040, // 7
            OUTPUT_CONT_STRESS            = 0x000080, // 8
            OUTPUT_TOTAL_DEFORMATION      = 0x000100, // 9
            OUTPUT_CENT_TEMPERATURE       = 0x000200, // 10
            OUTPUT_REACTION_FORCE         = 0x000400, // 11
            OUTPUT_FRICTION_LOSS          = 0x000800, // 12
            OUTPUT_VOLUME_FLUX_X          = 0x001000, // 13
            OUTPUT_VOLUME_FLUX_Z          = 0x002000, // 14
            OUTPUT_MASS_FLUX_X            = 0x004000, // 15
            OUTPUT_MASS_FLUX_Z            = 0x008000, // 16
            OUTPUT_HEAT_FLUX_X            = 0x010000, // 17
            OUTPUT_HEAT_FLUX_Z            = 0x020000, // 18
            OUTPUT_DEFORMATION1           = 0x040000, // 19
            OUTPUT_DEFORMATION2           = 0x080000, // 20
            OUTPUT_MESH                   = 0x100000, // 21
            OUTPUT_VOLUME_FLUX            = OUTPUT_VOLUME_FLUX_X | OUTPUT_VOLUME_FLUX_Z,
            OUTPUT_MASS_FLUX              = OUTPUT_MASS_FLUX_X | OUTPUT_MASS_FLUX_Z,
            OUTPUT_HEAT_FLUX              = OUTPUT_HEAT_FLUX_X | OUTPUT_HEAT_FLUX_Z,
            OUTPUT_FLUX_X                 = OUTPUT_VOLUME_FLUX_X | OUTPUT_MASS_FLUX_X | OUTPUT_HEAT_FLUX_X,
            OUTPUT_FLUX_Z                 = OUTPUT_VOLUME_FLUX_Z | OUTPUT_MASS_FLUX_Z | OUTPUT_HEAT_FLUX_Z,
            OUTPUT_FLUX                   = OUTPUT_VOLUME_FLUX | OUTPUT_MASS_FLUX | OUTPUT_HEAT_FLUX,
            OUTPUT_FIRST = OUTPUT_PRESSURE,
            OUTPUT_LAST = OUTPUT_MESH
        };

        enum ScaleType {
            SCALE_REYNOLDS_EQU = 0,
            SCALE_ELASTICITY_EQU,
            SCALE_ENERGY_EQ,
            SCALE_PRESSURE_DOF,
            SCALE_THETA_DOF,
            SCALE_TEMPERATURE_DOF,
            SCALE_LAST
        };

        enum InitialAssemblyDomain {
            INIT_ASS_NONE       = 0x0,
            INIT_ASS_HYDRAULIC  = 0x1,
            INIT_ASS_THERMAL    = 0x2,
            INIT_ASS_ELASTICITY = 0x4,
            INIT_ASS_ALL = INIT_ASS_HYDRAULIC | INIT_ASS_THERMAL | INIT_ASS_ELASTICITY
        };

        HydroRootElement(unsigned uLabel, const DofOwner *pDO,
                         DataManager* pDM, MBDynParser& HP);
        virtual ~HydroRootElement(void);
        virtual void Output(OutputHandler& OH) const;
        virtual void WorkSpaceDim(integer* piNumRows, integer* piNumCols) const;
        virtual unsigned int iGetNumDof(void) const;
        virtual DofOrder::Order GetDofType(unsigned int i) const;
        virtual DofOrder::Order GetEqType(unsigned int i) const;
        virtual std::ostream& DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const;
        virtual std::ostream& DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const;
        VariableSubMatrixHandler&
        AssJac(VariableSubMatrixHandler& WorkMat,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr);
        SubVectorHandler&
        AssRes(SubVectorHandler& WorkVec,
               doublereal dCoef,
               const VectorHandler& XCurr,
               const VectorHandler& XPrimeCurr);
        unsigned int iGetNumPrivData(void) const;
        virtual unsigned int iGetPrivDataIdx(const char *s) const;
        virtual doublereal dGetPrivData(unsigned int i) const;
        int GetNumConnectedNodes(void) const;
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
        virtual void AfterPredict(VectorHandler& X, VectorHandler& XP);
        virtual void Update(const VectorHandler& XCurr,
                            const VectorHandler& XPrimeCurr);
        virtual void AfterConvergence(const VectorHandler& X,
                                      const VectorHandler& XP);
        inline doublereal dGetTime() const;
        inline void AddNode(std::unique_ptr<Node2D>&& pNode);
        inline void AddElement(std::unique_ptr<HydroElement>&& pElement);
        inline void AddBoundaryCondition(std::unique_ptr<FluidStateBoundaryCond>&& pBoundCond);
        inline integer iGetNumNodes() const;
        inline Node2D* pGetNode(integer iNodeIndex) const;
        template <typename T>
        inline T* pGetNode(integer iNodeIndex) const;
        inline integer iGetNumElements() const;
        inline HydroElement* pGetElement(integer iElementIndex) const;
        inline integer iGetNumBoundaryConditions() const;
        inline FluidStateBoundaryCond* pGetBoundaryCondition(integer iBoundaryCondIndex) const;
        inline const HydroFluid* pGetFluid() const;
        inline const HydroMesh* pGetMesh() const;
        inline doublereal dGetStartupFactor() const;
        inline doublereal dGetScale(ScaleType eType) const;
        inline void AddFrictionLoss(enum FrictionLossType type, doublereal Pf);
        inline bool bUpdateFrictionLoss() const;
        inline void SetMaxTimeStep(doublereal dTimeStep);
        inline doublereal dGetMaxCFL() const;
        doublereal dGetMaxPressureGradient() const;
        inline void SetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintReal eType, doublereal dHint);
        inline integer GetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintInteger eType) const;
        inline doublereal GetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintReal eType) const;
        inline bool bInitialAssembly(InitialAssemblyDomain eDomain) const {
            return (uInitAssFlags & eDomain) != 0u;
        }
        unsigned uGetOutputFlags() const { return uOutputFlags; }
#if HYDRO_DEBUG > 0
        const DataManager* pGetDataManager() const { return pDM; }
#endif
    private:
        typedef std::vector<std::unique_ptr<Node2D>, GradientAllocator<std::unique_ptr<Node2D> > > NodeContainer;
        typedef std::map<integer, HydroDofOwner*, std::less<integer>, GradientAllocator<std::pair<integer, HydroDofOwner*> > > DofOwnerMap;
        typedef std::vector<HydroDofOwner*, GradientAllocator<HydroDofOwner*> > DofOwnerContainer;
        typedef std::vector<std::unique_ptr<HydroElement>, GradientAllocator<std::unique_ptr<HydroElement> > > ElementContainer;
        typedef std::vector<std::unique_ptr<FluidStateBoundaryCond>, GradientAllocator<std::unique_ptr<FluidStateBoundaryCond> > > BoundaryCondContainer;

        const DofOwnerMap& GetDofOwnerMap(grad::FunctionCall eFunc) const {
            HYDRO_ASSERT((eFunc & grad::STATE_MASK) == grad::REGULAR_FLAG
                         || (eFunc & grad::STATE_MASK) == grad::INITIAL_ASS_FLAG);

            return (eFunc & grad::REGULAR_FLAG)
                ? oDofOwnerReg
                : oDofOwnerInitAss;
        }

        DofOwnerMap& GetDofOwnerMap(grad::FunctionCall eFunc) {
            HYDRO_ASSERT((eFunc & grad::STATE_MASK) == grad::REGULAR_FLAG
                         || (eFunc & grad::STATE_MASK) == grad::INITIAL_ASS_FLAG);

            return (eFunc & grad::REGULAR_FLAG)
                ? oDofOwnerReg
                : oDofOwnerInitAss;
        }

        inline void UnivWorkSpaceDim(integer* piNumRows, integer* piNumCols, FunctionCall func) const;
        inline void AddDofOwner(HydroDofOwner* pDofOwner);
        void RebuildDofMap();
        inline void InitPrivData();
        inline const HydroDofOwner* pFindDofOwner(unsigned int i, grad::FunctionCall eFunc) const;
        void PrintNodeHeader(const Node2D* pNode, Node2D::NodeType eType, std::ostream& out) const;

        DataManager* const pDM;
        unsigned uOutputFlags;
        std::unique_ptr<HydroFluid> pFluid;
        std::unique_ptr<HydroMesh> pMesh;
        std::array<doublereal, SCALE_LAST> rgScale;
        std::unique_ptr<DriveCaller> pStartupFactor;
        NodeContainer rgNodes;
        DofOwnerMap oDofOwnerInitAss;
        DofOwnerMap oDofOwnerReg;
        DofOwnerContainer rgDofOwner;
        ElementContainer rgElements;
        BoundaryCondContainer rgBoundaryCond;
        LocalDofMap oNodeDofMap;
        doublereal dCFL;
        unsigned uInitAssFlags;

        mutable doublereal dMaxPressGradient;

        struct PrivDataVal {
            doublereal dCurr;
            doublereal dPrev;
        };

        static const index_type iNumNodeOutLoc = 4;
        static const Node2D::NodeType rgNodeOutLoc[iNumNodeOutLoc];
        static const int iNumFrictionLoss = 2;
        static const int iNumReactionForce = 12;
        static const int iNumPrivData = 12 + iNumFrictionLoss + iNumReactionForce;

        union PrivDataU {
            PrivDataVal a[iNumPrivData];
            struct PrivDataS {
                PrivDataVal MaxTimeStep;
                PrivDataVal rgPf[iNumFrictionLoss];
                PrivDataVal Maxp;
                PrivDataVal Maxpc;
                PrivDataVal Minh;
                PrivDataVal Minwtot;
                PrivDataVal Maxwtot;
                PrivDataVal Minrho;
                PrivDataVal Maxrho;
                PrivDataVal MinT;
                PrivDataVal MaxT;
                PrivDataVal MeanT;
                PrivDataVal ReactionForce[iNumReactionForce];
            } s;
        } PrivData;

        mutable bool bUpdatePrivData;

        static const struct PrivateData {
            char szName[8];
            doublereal dDefault;
        } rgPrivData[iNumPrivData];

#if CREATE_PROFILE == 1
        enum {
            PROF_RES = 0,
            PROF_JAC = 1
        };
        struct {
            doublereal dtRes;
            doublereal dtJac;
            doublereal dtUpdateGeom[2];
            doublereal dtUpdateNodes[2];
            doublereal dtElemAss[2];
            doublereal dtGeomAss[2];
        } profile;
#endif
    };

    const Node2D::NodeType HydroRootElement::rgNodeOutLoc[iNumNodeOutLoc] = {
        Node2D::HYDRAULIC_NODE,
        Node2D::THERMAL_NODE,
        Node2D::FLUX_NODE_X,
        Node2D::FLUX_NODE_Z
    };

    const struct HydroRootElement::PrivateData
    HydroRootElement::rgPrivData[iNumPrivData] = {
        {"max" "dt",  std::numeric_limits<doublereal>::max()},
        {"Pff",           0.},
        {"Pfc",           0.},
        {"max" "p",   -std::numeric_limits<doublereal>::max()},
        {"max" "pc",  -std::numeric_limits<doublereal>::max()},
        {"min" "h",    std::numeric_limits<doublereal>::max()},
        {"min" "wtot", std::numeric_limits<doublereal>::max()},
        {"max" "wtot", -std::numeric_limits<doublereal>::max()},
        {"min" "rho",  std::numeric_limits<doublereal>::max()},
        {"max" "rho", -std::numeric_limits<doublereal>::max()},
        {"min" "T",    std::numeric_limits<doublereal>::max()},
        {"max" "T",   -std::numeric_limits<doublereal>::max()},
        {"mean" "T",   0.},
        {"F1x", 0.},
        {"F1y", 0.},
        {"F1z", 0.},
        {"M1x", 0.},
        {"M1y", 0.},
        {"M1z", 0.},
        {"F2x", 0.},
        {"F2y", 0.},
        {"F2z", 0.},
        {"M2x", 0.},
        {"M2y", 0.},
        {"M2z", 0.}
    };

    HydroRootElement::HydroRootElement(
        unsigned uLabel, const DofOwner *pDO,
        DataManager* pDM, MBDynParser& HP)
        :       Elem(uLabel, flag(0)),
                UserDefinedElem(uLabel, pDO),
                pDM(pDM),
                uOutputFlags(OUTPUT_NOTHING),
                pFluid(nullptr),
                pMesh(nullptr),
                dCFL(1.),
                uInitAssFlags(INIT_ASS_ALL),
                dMaxPressGradient(-1.),
                bUpdatePrivData(false)
    {
        std::fill(rgScale.begin(), rgScale.end(), 1.);

        InitPrivData();
        // help
        if (HP.IsKeyWord("help"))
        {
            silent_cout(
                "\n"
                "Module:        HydrodynamicPlainBearingFD\n"
                "\n"
                "       This module implements a hydrodynamic plain bearing according to\n"
                "\n"
                "   Hans Juergen Butenschoen 1976\n"
                "   Das hydrodynamische zylindrische Gleitlager\n"
                "       endlicher Breite unter instationaerer Belastung\n"
                << std::endl);

            if (!HP.IsArg())
            {
                /*
                 * Exit quietly if nothing else is provided
                 */
                throw NoErr(MBDYN_EXCEPT_ARGS);
            }
        }

        std::unique_ptr<HydraulicFluid> pMbFluid{HP.GetHydraulicFluid()};

        const doublereal eta = pMbFluid->dGetViscosity();
        const doublereal pc = pMbFluid->dGetPres0();
        const doublereal T0 = pMbFluid->dGetTemp0();
        const doublereal rhoc = pMbFluid->dGetDensity();
        const doublereal drho_dp = pMbFluid->dGetDensityDPres();
        const doublereal drho_dT = pMbFluid->dGetDensityDTemp();
        const doublereal beta = drho_dp != 0. ? rhoc / drho_dp : std::numeric_limits<doublereal>::infinity();
        const doublereal gamma = -drho_dT / rhoc; // Definition according to Dirk Bartel 2009

        if (pc < 0) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                        << "): illegal value for reference pressure (fluid:"
                        << pMbFluid->GetLabel() << " " << pMbFluid->GetName()
                        << ") at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        HydroFluid::HydraulicType fluidType = (drho_dp == 0.)
            ? HydroFluid::INCOMPRESSIBLE
            : HydroFluid::COMPRESSIBLE;

        if (drho_dp != 0. && rhoc == 0.) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                        << "): reference density must not be zero for a compressible fluid at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        ThermalFluidModel oThermModel(T0, rhoc, eta, gamma);

        oThermModel.ParseInput(pDM, HP, this);

        doublereal etav = 0.;
        bool bGotEtaVapor = false;

        if (HP.IsKeyWord("viscosity" "vapor") ||  HP.IsKeyWord("viscosity" "vapour")) {
            if (HP.IsKeyWord("factor")) {
                etav = HP.GetReal() * eta;
            } else {
                etav = HP.GetReal();
            }
            bGotEtaVapor = true;
        }

        if (drho_dp != 0) {
            if (!bGotEtaVapor) {
                silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                            << "): keyword \"viscosity vapor\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (etav <= 0) {
                silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                            << "): viscosity vapor must be greater than zero at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (HP.IsKeyWord("bayada" "chupin")) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                        << "): fluid model \"bayad chupin\" is obsolte at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        } else {
            if (drho_dp == 0.) {
                pFluid.reset(new HydroIncompressibleFluid(pc, oThermModel));
            } else {
                pFluid.reset(new LinearCompressibleFluid(etav / eta, beta, pc, fluidType, oThermModel));
            }
        }

        if ( !HP.IsKeyWord("mesh")) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel() << "): keyword \"mesh\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("linear" "finite" "difference")) {
            pMesh.reset(new LinFDMesh(this));
        } else if (HP.IsKeyWord("quadratic" "finite" "element")) {
            pMesh.reset(new QuadFeIso9Mesh(this));
        } else {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel() << "): keyword \"linear finite difference\" or \"quadratic finite element\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        rgScale[SCALE_REYNOLDS_EQU] = 12. * eta / rhoc; // default scale for backward compatibility

        pMesh->ParseInput(pDM, HP);

        pStartupFactor.reset(HP.IsKeyWord("startup" "factor")
                             ? HP.GetDriveCaller()
                             : new OneDriveCaller());

        if (HP.IsKeyWord("initial" "assembly")) {
            static const struct DomainData {
                InitialAssemblyDomain eDomain;
                char szDomain[11];
            } rgDomains[] = {
                {INIT_ASS_NONE,       "none"},
                {INIT_ASS_ALL,        "all"},
                {INIT_ASS_HYDRAULIC,  "hydraulic"},
                {INIT_ASS_THERMAL,    "thermal"},
                {INIT_ASS_ELASTICITY, "elasticity"}
            };

            for (auto i = std::begin(rgDomains); i != std::end(rgDomains); ++i) {
                if (HP.IsKeyWord(i->szDomain)) {
                    if (HP.GetYesNoOrBool()) {
                        uInitAssFlags |= i->eDomain;
                    } else {
                        uInitAssFlags &= ~i->eDomain;
                    }
                }
            }
        }

        while (HP.IsArg()) {
            if (HP.IsKeyWord("pressure" "dof" "scale")) {
                rgScale[SCALE_PRESSURE_DOF] = HP.GetReal();
            } else if (HP.IsKeyWord("theta" "dof" "scale")) {
                rgScale[SCALE_THETA_DOF] = HP.GetReal();
            } else if (HP.IsKeyWord("temperature" "dof" "scale")) {
                rgScale[SCALE_TEMPERATURE_DOF] = HP.GetReal();
            } else if (HP.IsKeyWord("equation" "scale")) {
                rgScale[SCALE_ELASTICITY_EQU] =
                    rgScale[SCALE_REYNOLDS_EQU] =
                    rgScale[SCALE_ENERGY_EQ] = HP.GetReal();
            } else if (HP.IsKeyWord("reynolds" "equation" "scale")) {
                rgScale[SCALE_REYNOLDS_EQU] = HP.GetReal();
            } else if (HP.IsKeyWord("elasticity" "equation" "scale")) {
                rgScale[SCALE_ELASTICITY_EQU] = HP.GetReal();
            } else if (HP.IsKeyWord("energy" "balance" "scale")) {
                rgScale[SCALE_ENERGY_EQ] = HP.GetReal();
            } else if (HP.IsKeyWord("max" "cfl" "number")) {
                dCFL = HP.GetReal();
            } else if (HP.IsKeyWord("output" "pressure")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_PRESSURE;
                }
            } else if (HP.IsKeyWord("output" "density")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_DENSITY;
                }
            } else if (HP.IsKeyWord("output" "contact" "pressure")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_CONT_PRESSURE;
                }
            } else if (HP.IsKeyWord("output" "clearance")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_CLEARANCE;
                }
            } else if (HP.IsKeyWord("output" "clearance" "derivative")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_CLEARANCE_DER;
                }
            } else if (HP.IsKeyWord("output" "velocity")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_VELOCITY;
                }
            } else if (HP.IsKeyWord("output" "stress")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_STRESS;
                }
            } else if (HP.IsKeyWord("output" "contact" "stress")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_CONT_STRESS;
                }
            } else if (HP.IsKeyWord("output" "reaction" "force")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_REACTION_FORCE;
                }
            } else if (HP.IsKeyWord("output" "total" "deformation")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_TOTAL_DEFORMATION;
                }
            } else if (HP.IsKeyWord("output" "deformation" "shaft")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_DEFORMATION1;
                }
            } else if (HP.IsKeyWord("output" "deformation" "bearing")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_DEFORMATION2;
                }
            } else if (HP.IsKeyWord("output" "volume" "flux")) {
                bool bGotKeyWord = false;

                if (HP.IsKeyWord("x")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_VOLUME_FLUX_X;
                    }
                }

                if (HP.IsKeyWord("z")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_VOLUME_FLUX_Z;
                    }
                }

                if (!bGotKeyWord && HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_VOLUME_FLUX;
                }
            } else if (HP.IsKeyWord("output" "mass" "flux")) {
                bool bGotKeyWord = false;

                if (HP.IsKeyWord("x")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_MASS_FLUX_X;
                    }
                }

                if (HP.IsKeyWord("z")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_MASS_FLUX_Z;
                    }
                }

                if (!bGotKeyWord && HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_MASS_FLUX;
                }

            } else if (HP.IsKeyWord("output" "heat" "flux")) {
                bool bGotKeyWord = false;

                if (HP.IsKeyWord("x")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_HEAT_FLUX_X;
                    }
                }

                if (HP.IsKeyWord("z")) {
                    bGotKeyWord = true;

                    if (HP.GetYesNoOrBool()) {
                        uOutputFlags |= OUTPUT_HEAT_FLUX_Z;
                    }
                }

                if (!bGotKeyWord && HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_HEAT_FLUX;
                }
            } else if (HP.IsKeyWord("output" "temperature")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_CENT_TEMPERATURE;
                }
            } else if (HP.IsKeyWord("output" "mesh")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_MESH;
                }
            } else if (HP.IsKeyWord("output" "friction" "loss")) {
                if (HP.GetYesNoOrBool()) {
                    uOutputFlags |= OUTPUT_FRICTION_LOSS;
                    bUpdatePrivData = true;
                }
            } else if (HP.IsKeyWord("enable" "private" "data")) {
                if (HP.GetYesNoOrBool()) {
                    bUpdatePrivData = true;
                }
            } else {
                break;
            }
        }

        SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));

        if (!fToBeOutput()) {
            // The global output flags should overwrite all local flags
            uOutputFlags = OUTPUT_NOTHING;
        }

        rgNodes.resize(pMesh->iGetNumNodes()); // use resize because nodes can be generated in an arbitrary order
        rgElements.reserve(pMesh->iGetNumElements());
        rgDofOwner.reserve(rgNodes.capacity() + rgElements.capacity());
        rgBoundaryCond.reserve(pMesh->iGetNumBounaryConditions());
        pMesh->Generate();

        HYDRO_ASSERT(rgNodes.capacity() == rgNodes.size());
        HYDRO_ASSERT(rgNodes.size() == static_cast<size_t>(pMesh->iGetNumNodes()));

        for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
            HYDRO_ASSERT(*i != nullptr);

            AddDofOwner(dynamic_cast<HydroDofOwner*>(i->get()));
        }

        for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
            HYDRO_ASSERT(*i != nullptr);

            (*i)->Initialize();

            AddDofOwner(dynamic_cast<HydroDofOwner*>(i->get()));
        }

        pMesh->pGetGeometry()->Initialize();

        RebuildDofMap();

        HYDRO_TRACE("initial number of degrees of freedom: " << iGetInitialNumDof() << std::endl);
        HYDRO_TRACE("regular number of degrees of freedom: " << iGetNumDof() << std::endl);

        if (pedantic_output) {
            for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
                pedantic_cout("node(" << (*i)->iGetNodeNumber() + 1 << "):" << (*i)->GetPosition2D() << "\n");
            }
            for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                pedantic_cout("element(" << i - rgElements.begin() + 1 << "):" << **i << "\n");
            }
        }

        std::ostream& out = pDM->GetLogFile();

        out << "hydrodynamic plain bearing2: " << GetLabel() << ' ';

        out << eta << ' ';

        int iFlags = 0;

        for (unsigned f = OUTPUT_FIRST; f <= OUTPUT_LAST; f <<= 1) {
            ++iFlags;
        }

        out << iFlags << ' ';

        for (unsigned f = OUTPUT_FIRST; f <= OUTPUT_LAST; f <<= 1) {
            out << bool(uOutputFlags & f) << ' ';
        }

        std::array<index_type, iNumNodeOutLoc> rgNumNodes = {0};

        for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
            for (index_type j = 0; j < iNumNodeOutLoc; ++j) {
                if ((*i)->bIsNodeType(rgNodeOutLoc[j])) {
                    ++rgNumNodes[j];
                    break;
                }
            }
        }

        for (index_type j = 0; j < iNumNodeOutLoc; ++j) {
            out << rgNumNodes[j] << ' ';

            for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
                const Node2D* const pNode = i->get();
                PrintNodeHeader(pNode, rgNodeOutLoc[j], out);
            }
        }

        out << iGetNumElements() << ' ';

        for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
            const HydroElement* pElement = i->get();
            const integer iNumNodes = pElement->iGetNumNodes();

            out << i - rgElements.begin() + 1 << ' ' << iNumNodes << ' ';

            for (integer j = 0; j < iNumNodes; ++j) {
                out << pElement->pGetNode(j)->iGetNodeNumber() + 1 << ' ';
            }
        }

        pMesh->PrintLogFile(out);

        out << std::endl;

#if CREATE_PROFILE == 1
        memset(&profile, 0, sizeof(profile));
#endif
    }

    void HydroRootElement::AddDofOwner(HydroDofOwner* pDofOwner)
    {
        if (pDofOwner == nullptr) {
            return;
        }

        HYDRO_ASSERT(pDofOwner->iGetNumDof() != 0);

        rgDofOwner.push_back(pDofOwner);
    }

    void HydroRootElement::RebuildDofMap()
    {
        const index_type iNumDofMaps = 2;
        const std::array<grad::FunctionCall, iNumDofMaps> rgFuncs = {INITIAL_ASS_RES, REGULAR_RES};

        for (index_type j = 0; j < iNumDofMaps; ++j) {
            DofOwnerMap& oDofOwner = GetDofOwnerMap(rgFuncs[j]);
            index_type iOffsetIndex = 1; // start with one because we always need iGetFirstIndex() + 1:N

            oDofOwner.clear();

            for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
                HydroDofOwner* const pDofOwner = *i;

                const integer iNumDof = (rgFuncs[j] == INITIAL_ASS_RES)
                    ? pDofOwner->iGetInitialNumDof()
                    : pDofOwner->iGetNumDof();

                if (!iNumDof) {
                    continue;
                }

                pDofOwner->SetOffsetIndex(iOffsetIndex, rgFuncs[j]);

                const integer iFirstIndex = iOffsetIndex;
                const integer iLastIndex = iFirstIndex + iNumDof;

                for (integer j = iFirstIndex; j < iLastIndex; ++j) {
                    HYDRO_ASSERT(oDofOwner.find(j) == oDofOwner.end());
                    oDofOwner.insert(std::make_pair(j, pDofOwner));
                }

                iOffsetIndex += iNumDof;
            }
        }
    }

    HydroRootElement::~HydroRootElement(void)
    {
        // destroy private data

#if CREATE_PROFILE == 1
        std::cerr << "hydrodynamic plain bearing2(" << GetLabel() << ")" << std::endl;

        std::cerr << "dtRes=" << profile.dtRes << std::endl;

        std::cerr << "dtJac=" << profile.dtJac << std::endl;

        for (int i = PROF_RES; i <= PROF_JAC; ++i) {
            std::cerr << "dtUpdateGeom[" << i << "]=" << profile.dtUpdateGeom[i] << std::endl;
            std::cerr << "dtUpdateNodes[" << i << "]=" << profile.dtUpdateNodes[i] << std::endl;
            std::cerr << "dtElemAss[" << i << "]=" << profile.dtElemAss[i] << std::endl;
            std::cerr << "dtGeomAss[" << i << "]=" << profile.dtGeomAss[i] << std::endl;
        }

        std::cerr << std::endl;
#endif
    }

    void
    HydroRootElement::Output(OutputHandler& OH) const
    {
        if ( fToBeOutput()
             && OH.UseText(OutputHandler::LOADABLE)
             && uOutputFlags != OUTPUT_NOTHING) {

            std::ostream& os = OH.Loadable();

            os << std::setw(8) << GetLabel() << ' ';            // 0

            for (index_type j = 0; j < iNumNodeOutLoc; ++j) {
                for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
                    if ((*i)->bIsNodeType(rgNodeOutLoc[j])){
                        (*i)->Output(os, uOutputFlags);
                    }
                }
            }

            if (uOutputFlags & OUTPUT_REACTION_FORCE) {
                pMesh->pGetGeometry()->Output(os);
            }

            if (uOutputFlags & OUTPUT_FRICTION_LOSS) {
                for (int i = 0; i < iNumFrictionLoss; ++i) {
                    os << PrivData.s.rgPf[i].dPrev << ' ';
                }
            }

            if (uOutputFlags & OUTPUT_MESH) {
                pMesh->Output(os);
            }

            os << std::endl;
        }
    }



    void
    HydroRootElement::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
    {
        UnivWorkSpaceDim(piNumRows, piNumCols, REGULAR_JAC);
    }

    void HydroRootElement::UnivWorkSpaceDim(integer* piNumRows, integer* piNumCols, FunctionCall eFunc) const
    {
        integer iNumRowsTotal = 0;
        integer iNumItemsTotal = 0;

        integer iNumRowsElem = 0, iNumColsElem = 0;

        pMesh->pGetGeometry()->WorkSpaceDim(&iNumRowsElem, &iNumColsElem, eFunc);

        iNumItemsTotal += iNumRowsElem * iNumColsElem;
        iNumRowsTotal += iNumRowsElem;

        for (auto i = rgElements.cbegin(); i != rgElements.cend(); ++i) {
            (*i)->WorkSpaceDim(&iNumRowsElem, &iNumColsElem, eFunc);
            iNumItemsTotal += iNumRowsElem * iNumColsElem;
            iNumRowsTotal += iNumRowsElem;
        }

        *piNumRows = -iNumRowsTotal;
        *piNumCols = iNumItemsTotal / iNumRowsTotal + 1;

        HYDRO_TRACE("iNumRows=" << *piNumRows << std::endl);
        HYDRO_TRACE("iNumCols=" << *piNumCols << std::endl);
        HYDRO_TRACE("iMaxItems=" << *piNumRows * *piNumCols << std::endl);
    }

    unsigned int HydroRootElement::iGetNumDof(void) const
    {
        unsigned int iNumDof = 0;

        for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
            // Attention: don't use rgNodes here
            // because it could contain double entries
            // in case of periodic boundary conditions
            HYDRO_ASSERT((*i)->iGetNumDof() > 0);
            iNumDof += (*i)->iGetNumDof();
        }

        HYDRO_TRACE("iNumDof=" << iNumDof << std::endl);

        return iNumDof;
    }

    DofOrder::Order HydroRootElement::GetDofType(unsigned int i) const
    {
        // RebuildDofMap(REGULAR_RES);
        ++i; // we are using one based indices
        const HydroDofOwner* const pDO = pFindDofOwner(i, grad::REGULAR_RES);
        HYDRO_ASSERT(i >= unsigned(pDO->iGetOffsetIndex(grad::REGULAR_RES)));
        return pDO->GetDofType(i - pDO->iGetOffsetIndex(grad::REGULAR_RES));
    }

    DofOrder::Order HydroRootElement::GetEqType(unsigned int i) const
    {
        // RebuildDofMap(REGULAR_RES);
        ++i; // we are using one based indices
        const HydroDofOwner* const pDO = pFindDofOwner(i, grad::REGULAR_RES);
        HYDRO_ASSERT(i >= unsigned(pDO->iGetOffsetIndex(grad::REGULAR_RES)));
        return pDO->GetEqType(i - pDO->iGetOffsetIndex(grad::REGULAR_RES));
    }

    const HydroDofOwner* HydroRootElement::pFindDofOwner(unsigned int i, grad::FunctionCall eFunc) const
    {
        const DofOwnerMap& oDofOwnerMap = (eFunc & grad::REGULAR_FLAG)
            ? oDofOwnerReg
            : oDofOwnerInitAss;

        auto iNode = oDofOwnerMap.find(i);

        if (iNode == oDofOwnerMap.end()) {
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        HYDRO_ASSERT(iNode->first == integer(i));
        HYDRO_ASSERT(iNode->second->iGetOffsetIndex(eFunc) + iNode->second->iGetNumDof() > i);
        HYDRO_ASSERT(iNode->second->iGetOffsetIndex(eFunc) <= integer(i));

        return iNode->second;
    }

    void HydroRootElement::PrintNodeHeader(const Node2D* pNode, Node2D::NodeType eType, std::ostream& out) const
    {
        if (pNode->bIsNodeType(eType)) {
            const Vector<doublereal, 2>& x = pNode->GetPosition2D();
            const HydroDofOwner* const pDofOwner = dynamic_cast<const HydroDofOwner*>(pNode);
            const integer iOffset = pDofOwner && pDofOwner->iGetNumDof() ? pDofOwner->iGetOffsetIndex(grad::REGULAR_RES) : -1;
            out << pNode->iGetNodeNumber() + 1 << ' ' << x(1) << ' ' << x(2) << ' ' << iOffset << ' ';
        }
    }

    std::ostream& HydroRootElement::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
            (*i)->DescribeDof(out, prefix, bInitial);
        }

        return out;
    }

    std::ostream& HydroRootElement::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
            (*i)->DescribeEq(out, prefix, bInitial);
        }

        return out;
    }


    VariableSubMatrixHandler&
    HydroRootElement::AssJac(VariableSubMatrixHandler& WorkMatVar,
                             doublereal dCoef,
                             const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif
        oNodeDofMap.Reset(REGULAR_JAC);

        pMesh->Update(XCurr, XPrimeCurr, dCoef, REGULAR_JAC, &oNodeDofMap);

#if CREATE_PROFILE == 1
        profile.dtUpdateGeom[PROF_JAC] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
            (*i)->Update(XCurr, XPrimeCurr, dCoef, REGULAR_JAC, &oNodeDofMap);
        }

#if CREATE_PROFILE == 1
        profile.dtUpdateNodes[PROF_JAC] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        SparseSubMatrixHandler& WorkMat = WorkMatVar.SetSparse();
        WorkMat.Resize(1, 1);
        WorkMat.PutItem(1, 1, 1, 0); // FIXME: Avoid a segmentation fault if the matrix is empty

        enum GradientAssVecBase::mode_t mode = GradientAssVecBase::RESET;

        for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
            (*i)->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr, mode);
            mode = GradientAssVecBase::APPEND;
        }

#if CREATE_PROFILE == 1
        profile.dtElemAss[PROF_JAC] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        pMesh->pGetGeometry()->AssJac(WorkMat, dCoef, XCurr, XPrimeCurr, mode);

#if CREATE_PROFILE == 1
        profile.dtGeomAss[PROF_JAC] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        return WorkMatVar;
    }

    SubVectorHandler&
    HydroRootElement::AssRes(SubVectorHandler& WorkVec,
                             doublereal dCoef,
                             const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif
        if (bUpdatePrivData) {
            for (int i = 0; i < iNumPrivData; ++i) {
                // Reset the current state of friction loss and others
                // DriveCallers may still access PrivData.a[i].dPrev during AssRes
                PrivData.a[i].dCurr = rgPrivData[i].dDefault;
            }
        }

        pMesh->Update(XCurr, XPrimeCurr, dCoef, REGULAR_RES, 0);

#if CREATE_PROFILE == 1
        profile.dtUpdateGeom[PROF_RES] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        for (auto i = rgBoundaryCond.begin(); i != rgBoundaryCond.end(); ++i) {
            (*i)->Update();
        }

        dMaxPressGradient = -1.;

        for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
            (*i)->Update(XCurr, XPrimeCurr, dCoef, REGULAR_RES, 0);
        }

#if CREATE_PROFILE == 1
        profile.dtUpdateNodes[PROF_RES] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        enum GradientAssVecBase::mode_t mode = GradientAssVecBase::RESET;

        for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
            (*i)->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr, mode);
            mode = GradientAssVecBase::APPEND;
        }

#if CREATE_PROFILE == 1
        profile.dtElemAss[PROF_RES] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        pMesh->pGetGeometry()->AssRes(WorkVec, dCoef, XCurr, XPrimeCurr, mode);

#if CREATE_PROFILE == 1
        profile.dtGeomAss[PROF_RES] += mbdyn_clock_time() - start;
        start = mbdyn_clock_time();
#endif

        if (bUpdatePrivData) {
            index_type iNumT = 0;

            for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
                doublereal h;

                if ((*i)->bGetPrivateData(PD_CLEARANCE, h) && h < PrivData.s.Minh.dCurr) {
                    PrivData.s.Minh.dCurr = h;
                }

                doublereal p;

                if ((*i)->bGetPrivateData(PD_PRESSURE, p) && p > PrivData.s.Maxp.dCurr) {
                    PrivData.s.Maxp.dCurr = p;
                }

                doublereal pc;

                if ((*i)->bGetPrivateData(PD_CONT_PRESSURE, pc) &&  pc > PrivData.s.Maxpc.dCurr) {
                    PrivData.s.Maxpc.dCurr = pc;
                }

                doublereal rho;

                if ((*i)->bGetPrivateData(PD_DENSITY, rho)) {
                    if (rho < PrivData.s.Minrho.dCurr) {
                        PrivData.s.Minrho.dCurr = rho;
                    }

                    if (rho > PrivData.s.Maxrho.dCurr) {
                        PrivData.s.Maxrho.dCurr = rho;
                    }
                }

                doublereal T;

                if ((*i)->bGetPrivateData(PD_TEMPERATURE, T)) {
                    if (T < PrivData.s.MinT.dCurr) {
                        PrivData.s.MinT.dCurr = T;
                    }

                    if (T > PrivData.s.MaxT.dCurr) {
                        PrivData.s.MaxT.dCurr = T;
                    }

                    PrivData.s.MeanT.dCurr += T;
                    ++iNumT;
                }

                doublereal wtot;

                if ((*i)->bGetPrivateData(PD_TOTAL_DEFORMATION, wtot)) {
                    if (wtot < PrivData.s.Minwtot.dCurr) {
                        PrivData.s.Minwtot.dCurr = wtot;
                    }

                    if (wtot > PrivData.s.Maxwtot.dCurr) {
                        PrivData.s.Maxwtot.dCurr = wtot;
                    }
                }
            }

            for (int i = 0; i < iNumReactionForce; ++i) {
                pMesh->pGetGeometry()->bGetPrivateData(static_cast<PrivateDataType>(i + PD_F1x), PrivData.s.ReactionForce[i].dCurr);
            }

            if (iNumT) {
                PrivData.s.MeanT.dCurr /= iNumT;
            }

            for (int i = 0; i < iNumPrivData; ++i) {
                // DriveCallers will see up to date data from this point
                PrivData.a[i].dPrev = PrivData.a[i].dCurr;
            }
        }

        return WorkVec;
    }


    unsigned int
    HydroRootElement::iGetNumPrivData(void) const
    {
        return iNumPrivData;
    }

    unsigned int HydroRootElement::iGetPrivDataIdx(const char *s) const
    {
        for (int i = 0; i < iNumPrivData; ++i) {
            if (0 == strcmp(rgPrivData[i].szName, s)) {
                bUpdatePrivData = true;
                return i + 1;
            }
        }

        silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                    << "): private data \"" << s << "\" not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    doublereal HydroRootElement::dGetPrivData(unsigned int i) const
    {
        if (!bUpdatePrivData) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                        << "): private data " << i << " cannot be accessed "
                        "because \"enable private data, yes\" "
                        "was not specified!" << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        --i;

        if (i >= iNumPrivData) {
            silent_cerr("hydrodynamic plain bearing2(" << GetLabel()
                        << "): invalid private data index " << i + 1 << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return PrivData.a[i].dPrev;
    }

    int
    HydroRootElement::GetNumConnectedNodes(void) const
    {
        HYDRO_ASSERT(pMesh != nullptr);
        HYDRO_ASSERT(pMesh->pGetGeometry() != nullptr);

        return pMesh->pGetGeometry()->iGetNumConnectedNodes();
    }

    void
    HydroRootElement::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
    {
        HYDRO_ASSERT(pMesh != nullptr);
        HYDRO_ASSERT(pMesh->pGetGeometry() != nullptr);

        pMesh->pGetGeometry()->GetConnectedNodes(connectedNodes);
    }

    void
    HydroRootElement::SetValue(DataManager *pDM,
                               VectorHandler& X, VectorHandler& XP,
                               SimulationEntity::Hints *ph)
    {
        for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
            (*i)->SetValue(X, XP);
        }
    }

    std::ostream&
    HydroRootElement::Restart(std::ostream& out) const
    {
        out << "# hydrodynamic plain bearing2: restart not implemented\n";

        return out;
    }

    unsigned int
    HydroRootElement::iGetInitialNumDof(void) const
    {
        unsigned int iNumDof = 0;

        if (uInitAssFlags != INIT_ASS_NONE) {
            for (auto i = rgDofOwner.cbegin(); i != rgDofOwner.cend(); ++i) {
                iNumDof += (*i)->iGetInitialNumDof();
            }
        }

        return iNumDof;
    }

    void
    HydroRootElement::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
    {
        if (uInitAssFlags != INIT_ASS_NONE) {
            UnivWorkSpaceDim(piNumRows, piNumCols, INITIAL_ASS_JAC);
        } else {
            *piNumRows = 0;
            *piNumCols = 0;
        }
    }

    VariableSubMatrixHandler&
    HydroRootElement::InitialAssJac(
        VariableSubMatrixHandler& WorkMatVar,
        const VectorHandler& XCurr)
    {
        if (uInitAssFlags != INIT_ASS_NONE) {
            oNodeDofMap.Reset(INITIAL_ASS_JAC);

            const VectorHandler* const pXPrimeCurr = nullptr;

            pMesh->Update(XCurr, *pXPrimeCurr, 1., INITIAL_ASS_JAC, &oNodeDofMap);

            HYDRO_ASSERT(rgElements.size() > 0);

            for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
                (*i)->Update(XCurr, *pXPrimeCurr, 1., INITIAL_ASS_JAC, &oNodeDofMap);
            }

            SparseSubMatrixHandler& WorkMat = WorkMatVar.SetSparse();
            WorkMat.Resize(1, 1);
            WorkMat.PutItem(1, 1, 1, 0); // FIXME: Avoid a segmentation fault if the matrix is empty

            enum GradientAssVecBase::mode_t mode = GradientAssVecBase::RESET;

            for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                (*i)->InitialAssJac(WorkMat, XCurr, mode);
                mode = GradientAssVecBase::APPEND;
            }

            pMesh->pGetGeometry()->InitialAssJac(WorkMat, XCurr, mode);
        } else {
            WorkMatVar.SetNullMatrix();
        }

        return WorkMatVar;
    }

    SubVectorHandler&
    HydroRootElement::InitialAssRes(
        SubVectorHandler& WorkVec,
        const VectorHandler& XCurr)
    {
        if (uInitAssFlags != INIT_ASS_NONE) {
            HYDRO_ASSERT(rgElements.size() > 0);

            const VectorHandler* const pXPrimeCurr = nullptr;

            pMesh->Update(XCurr, *pXPrimeCurr, 1., INITIAL_ASS_RES, &oNodeDofMap);

            for (auto i = rgBoundaryCond.begin(); i != rgBoundaryCond.end(); ++i) {
                (*i)->Update();
            }

            for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
                (*i)->Update(XCurr, *pXPrimeCurr, 1., INITIAL_ASS_RES, 0);
            }

            enum GradientAssVecBase::mode_t mode = GradientAssVecBase::RESET;

            for (auto i = rgElements.begin(); i != rgElements.end(); ++i) {
                (*i)->InitialAssRes(WorkVec, XCurr, mode);
                mode = GradientAssVecBase::APPEND;
            }

            pMesh->pGetGeometry()->InitialAssRes(WorkVec, XCurr, mode);
        } else {
            WorkVec.ResizeReset(0);
        }

        return WorkVec;
    }

    void HydroRootElement::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
        for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
            (*i)->AfterPredict(X, XP);
        }

        for (auto i = rgElements.cbegin(); i != rgElements.cend(); ++i) {
            (*i)->AfterPredict(X, XP);
        }
    }

    void HydroRootElement::Update(const VectorHandler& XCurr,
                                  const VectorHandler& XPrimeCurr)
    {
        for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
            (*i)->DofUpdate(const_cast<VectorHandler&>(XCurr),
                            const_cast<VectorHandler&>(XPrimeCurr));
        }
    }

    void HydroRootElement::AfterConvergence(const VectorHandler& X,
                                            const VectorHandler& XP)
    {
        for (auto i = rgNodes.cbegin(); i != rgNodes.cend(); ++i) {
            (*i)->AfterConvergence(X, XP);
        }

        for (auto i = rgElements.cbegin(); i != rgElements.cend(); ++i) {
            (*i)->AfterConvergence(X, XP);
        }
    }

    inline doublereal HydroRootElement::dGetTime() const {
        return pDM->dGetTime();
    }

    inline void HydroRootElement::AddNode(std::unique_ptr<Node2D>&& pNode) {
        HYDRO_ASSERT(pNode != nullptr);
        HYDRO_ASSERT(pNode->iGetNodeNumber() >= 0);
        HYDRO_ASSERT(size_t(pNode->iGetNodeNumber()) < rgNodes.size());
        HYDRO_ASSERT(rgNodes[pNode->iGetNodeNumber()] == nullptr);

        rgNodes[pNode->iGetNodeNumber()] = std::move(pNode);
    }

    inline void HydroRootElement::AddElement(std::unique_ptr<HydroElement>&& pElement) {
        HYDRO_ASSERT(pElement != nullptr);
        HYDRO_ASSERT(rgElements.size() < rgElements.capacity());

        rgElements.push_back(std::move(pElement));
    }

    inline void HydroRootElement::AddBoundaryCondition(std::unique_ptr<FluidStateBoundaryCond>&& pBoundCond) {
        HYDRO_ASSERT(pBoundCond != nullptr);
        HYDRO_ASSERT(rgBoundaryCond.size() < rgBoundaryCond.capacity());

        rgBoundaryCond.push_back(std::move(pBoundCond));
    }

    inline integer HydroRootElement::iGetNumNodes() const {
        return rgNodes.size();
    }

    inline Node2D* HydroRootElement::pGetNode(integer iNodeIndex) const {
        HYDRO_ASSERT(iNodeIndex >= 0);
        HYDRO_ASSERT(size_t(iNodeIndex) < rgNodes.size());

        Node2D* pNode = rgNodes[iNodeIndex].get();

        HYDRO_ASSERT(pNode != nullptr);

        return pNode;
    }

    template <typename T>
    inline T* HydroRootElement::pGetNode(integer iNodeIndex) const {
        Node2D* pNode2D = pGetNode(iNodeIndex);

        HYDRO_ASSERT(pNode2D != nullptr);

        T* pNode = dynamic_cast<T*>(pNode2D);

        HYDRO_ASSERT(pNode != nullptr);

        if(pNode == nullptr) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return pNode;
    }

    inline integer HydroRootElement::iGetNumElements() const {
        return rgElements.size();
    }

    inline HydroElement* HydroRootElement::pGetElement(integer iElementIndex) const {
        HYDRO_ASSERT(iElementIndex >= 0);
        HYDRO_ASSERT(size_t(iElementIndex) < rgElements.size());

        return rgElements[iElementIndex].get();
    }

    inline integer HydroRootElement::iGetNumBoundaryConditions() const {
        return rgBoundaryCond.size();
    }

    inline FluidStateBoundaryCond* HydroRootElement::pGetBoundaryCondition(integer iBoundaryCondIndex) const {
        HYDRO_ASSERT(iBoundaryCondIndex >= 0);
        HYDRO_ASSERT(size_t(iBoundaryCondIndex) < rgBoundaryCond.size());

        return rgBoundaryCond[iBoundaryCondIndex].get();
    }

    inline const HydroFluid* HydroRootElement::pGetFluid() const {
        return pFluid.get();
    }

    inline const HydroMesh* HydroRootElement::pGetMesh() const {
        return pMesh.get();
    }

    inline doublereal HydroRootElement::dGetScale(ScaleType eType) const
    {
        HYDRO_ASSERT(eType >= 0);
        HYDRO_ASSERT(eType < SCALE_LAST);

        return rgScale[eType];
    }

    inline void HydroRootElement::AddFrictionLoss(enum FrictionLossType type, doublereal Pf) {
        HYDRO_ASSERT(type >= 0 && type < iNumFrictionLoss);
        PrivData.s.rgPf[type].dCurr += Pf;
    }

    inline bool HydroRootElement::bUpdateFrictionLoss() const {
        return bUpdatePrivData;
    }

    inline void HydroRootElement::SetMaxTimeStep(doublereal dTimeStep) {
        HYDRO_ASSERT(dTimeStep >= 0);

        if (dTimeStep < PrivData.s.MaxTimeStep.dCurr) {
            PrivData.s.MaxTimeStep.dCurr = dTimeStep;
        }
    }

    inline doublereal HydroRootElement::dGetMaxCFL() const {
        return dCFL;
    }

    void HydroRootElement::SetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintReal eType, doublereal dHint)
    {
        switch (eType) {
        case NonlinearSolver::LINESEARCH_LAMBDA_MAX: {
            NonlinearSolver* const pNLS = pDM->pGetNonlinearSolver();
            const doublereal dLambdaMin = pNLS->GetNonlinearSolverHint(eType);

            if (dHint < dLambdaMin) {
                pNLS->SetNonlinearSolverHint(eType, dHint);
            }
        } break;

        default:
            break;
        }
    }

    integer HydroRootElement::GetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintInteger eType) const
    {
        return pDM->pGetNonlinearSolver()->GetNonlinearSolverHint(eType);
    }

    doublereal HydroRootElement::GetNonlinearSolverHint(NonlinearSolver::NonlinearSolverHintReal eType) const
    {
        return pDM->pGetNonlinearSolver()->GetNonlinearSolverHint(eType);
    }

    doublereal HydroRootElement::dGetMaxPressureGradient() const {
        if (dMaxPressGradient < 0) {
            doublereal pmin = std::numeric_limits<doublereal>::max();
            doublereal pmax = -pmin;

            for (auto i = rgNodes.begin(); i != rgNodes.end(); ++i) {
                doublereal pi, paspi;

                if ((*i)->bGetPrivateData(PD_PRESSURE, pi)) {
                    if ((*i)->bGetPrivateData(PD_CONT_PRESSURE, paspi)) {
                        pi += paspi;
                    }

                    pmin = std::min(pi, pmin);
                    pmax = std::max(pi, pmax);
                }
            }

            dMaxPressGradient =  pmax - pmin;

            HYDRO_ASSERT(dMaxPressGradient >= 0);
        }

        return dMaxPressGradient;
    }

    inline void HydroRootElement::InitPrivData() {
        for (int i = 0; i < iNumPrivData; ++i) {
            PrivData.a[i].dPrev = PrivData.a[i].dCurr = rgPrivData[i].dDefault;
        }
    }

    inline doublereal HydroRootElement::dGetStartupFactor() const
    {
        doublereal dInitAss = pStartupFactor->dGet();

        if (dInitAss < 0.) {
            dInitAss = 0.;
        }

        if (dInitAss > 1.) {
            dInitAss = 1.;
        }

        return dInitAss;
    }

    Geometry2D::Geometry2D(const Vector<doublereal, 2>& x)
        :x(x)
    {

    }

    Geometry2D::~Geometry2D()
    {

    }

    std::unique_ptr<Geometry2D> Geometry2D::Read(HydroRootElement* pRoot, MBDynParser& HP)
    {
        if (!HP.IsKeyWord("position")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"position\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        Vector<doublereal, 2> xc;

        for (integer i = 1; i <= 2; ++i) {
            xc(i) = HP.GetReal();
        }

        if (HP.IsKeyWord("circle")) {
            if (!HP.IsKeyWord("radius")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"radius\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal r = HP.GetReal();

            return std::unique_ptr<Geometry2D>{new Circle2D{xc, r}};
        } else if (HP.IsKeyWord("rectangle")) {
            if (!HP.IsKeyWord("width")) {
                silent_cerr("hydrodynamic plain bearing2(" <<
                            pRoot->GetLabel()
                            << "): keyword \"width\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal w = HP.GetReal();

            if (!HP.IsKeyWord("height")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"height\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal h = HP.GetReal();

            return std::unique_ptr<Geometry2D>{new Rectangle2D{xc, w, h}};
        } else if (HP.IsKeyWord("complete" "surface")) {
            return std::unique_ptr<Geometry2D>{new CompleteSurface2D{xc}};
        } else if (HP.IsKeyWord("surface" "grid")) {
            Vector<doublereal, DYNAMIC_SIZE> x, z;
            Matrix<int8_t, DYNAMIC_SIZE, DYNAMIC_SIZE> status;

            if (!HP.IsKeyWord("x")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"x\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal tolx = HP.IsKeyWord("tolerance") ? HP.GetReal() : 0.;

            if (tolx < 0.) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): tolerance must be greater than or equal to zero at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iGridX = HP.GetInt();

            if (iGridX < 2) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): at least two grid points are required in x direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            x.Resize(iGridX);

            for (integer i = 1; i <= iGridX; ++i) {
                x(i) = HP.GetReal();

                if (i > 1 && x(i) <= x(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): x coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("z")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"z\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal tolz = HP.IsKeyWord("tolerance") ? HP.GetReal() : 0.;

            if (tolz < 0.) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): tolerance must be greater than or equal to zero at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iGridZ = HP.GetInt();

            if (iGridZ < 2) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): at least two grid points are required in z direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            z.Resize(iGridZ);

            for (integer i = 1; i <= iGridZ; ++i) {
                z(i) = HP.GetReal();

                if (i > 1 && z(i) <= z(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): z coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("status")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"status\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            status.Resize(iGridX - 1, iGridZ - 1);

            for (integer i = 1; i < iGridX; ++i) {
                for (integer j = 1; j < iGridZ; ++j) {
                    if (HP.IsKeyWord("active")) {
                        status(i, j) = true;
                    } else if (HP.IsKeyWord("inactive")) {
                        status(i, j) = false;
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pRoot->GetLabel()
                                    << "): keyword \"active\" or \"inactive\" expected at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }
            }

            return std::unique_ptr<Geometry2D>{new SurfaceGrid2D{xc, x, z, tolx, tolz, status}};
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"circle\" or \"rectangle\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    Circle2D::Circle2D(const Vector<doublereal, 2>& x, doublereal r)
        :Geometry2D(x), r(r)
    {

    }

    std::unique_ptr<Geometry2D> Circle2D::Clone(const Vector<doublereal, 2>& x) const
    {
        return std::unique_ptr<Geometry2D>{new Circle2D(x, r)};
    }

    bool Circle2D::bPointIsInside(const Vector<doublereal, 2>& p1) const
    {
        const Vector<doublereal, 2> v = p1 - x;

        return sqrt(Dot(v, v)) <= r;
    }

    Rectangle2D::Rectangle2D(const Vector<doublereal, 2>& x, doublereal w, doublereal h)
        :Geometry2D(x), w(w), h(h)
    {

    }

    std::unique_ptr<Geometry2D> Rectangle2D::Clone(const Vector<doublereal, 2>& x) const
    {
        return std::unique_ptr<Geometry2D>{new Rectangle2D(x, w, h)};
    }

    bool Rectangle2D::bPointIsInside(const Vector<doublereal, 2>& p1) const
    {
        const bool bInside = std::abs(p1(1) - x(1)) <= 0.5 * w
            && std::abs(p1(2) - x(2)) <= 0.5 * h;

        HYDRO_TRACE("point p1(" << p1 << ") is " << (bInside ? "inside" : "outside")
                    << " rectangle " << w << "x" << h << " at x(" << x << ")" << std::endl);

        return bInside;
    }

    CompleteSurface2D::CompleteSurface2D(const Vector<doublereal, 2>& x)
        :Geometry2D(x)
    {
    }

    std::unique_ptr<Geometry2D> CompleteSurface2D::Clone(const Vector<doublereal, 2>& x) const
    {
        return std::unique_ptr<Geometry2D>{new CompleteSurface2D{x}};
    }

    bool CompleteSurface2D::bPointIsInside(const Vector<doublereal, 2>& p1) const
    {
        // per definition everything is inside
        return true;
    }

    SurfaceGrid2D::SurfaceGrid2D(const Vector<doublereal, 2>& xc,
                                 const Vector<doublereal, DYNAMIC_SIZE>& x,
                                 const Vector<doublereal, DYNAMIC_SIZE>& z,
                                 doublereal tolx,
                                 doublereal tolz,
                                 const Matrix<int8_t, DYNAMIC_SIZE, DYNAMIC_SIZE>& status)
        :Geometry2D(xc),
         tolx(tolx),
         tolz(tolz),
         x(x),
         z(z),
         status(status)
    {
        HYDRO_ASSERT(x.iGetNumRows() >= 2);
        HYDRO_ASSERT(z.iGetNumRows() >= 2);
        HYDRO_ASSERT(x.iGetNumRows() - 1 == status.iGetNumRows());
        HYDRO_ASSERT(z.iGetNumRows() - 1 == status.iGetNumCols());
    }

    std::unique_ptr<Geometry2D> SurfaceGrid2D::Clone(const Vector<doublereal, 2>& xc) const
    {
        return std::unique_ptr<Geometry2D>{new SurfaceGrid2D{xc, x, z, tolx, tolz, status}};
    }

    bool SurfaceGrid2D::bPointIsInside(const Vector<doublereal, 2>& p1) const
    {
        index_type ix = x.iGetNumRows() - 1;

        for (index_type i = 1; i < x.iGetNumRows(); ++i) {
            if (p1(1) <= x(i + 1)) {
                ix = i;
                break;
            }
        }

        index_type iz = z.iGetNumRows() - 1;

        for (index_type j = 1; j < z.iGetNumRows(); ++j) {
            if (p1(2) <= z(j + 1)) {
                iz = j;
                break;
            }
        }

        if (status(ix, iz)) {
            return true;
        }

        if (ix > 1 && status(ix - 1, iz) && p1(1) - tolx <= x(ix - 1)) {
            return true;
        }

        if (ix < x.iGetNumRows() - 1 && status(ix + 1, iz) && p1(1) + tolx >= x(ix + 1)) {
            return true;
        }

        if (iz > 1 && status(ix, iz - 1) && p1(2) - tolz <= z(iz - 1)) {
            return true;
        }

        if (iz < z.iGetNumRows() - 1 && status(ix, iz + 1) && p1(2) + tolz >= z(iz + 1)) {
            return true;
        }

        return false;
    }

    LubricationGroove::LubricationGroove(integer iLabel, std::unique_ptr<Geometry2D>&& pGeometry)
        :iLabel(iLabel),
         pGeometry(std::move(pGeometry))
    {
        HYDRO_ASSERT(this->pGeometry != nullptr);
    }


    LubricationGroove::~LubricationGroove()
    {

    }

    std::unique_ptr<LubricationGrooveMaster> LubricationGroove::Read(integer iLabel, HydroRootElement* pRoot, BearingGeometry* pBearingGeometry, MBDynParser& HP)
    {
        enum Type eType = pBearingGeometry->ReadLubricationGrooveType(HP);

        std::unique_ptr<FluidStateBoundaryCond> pBoundaryCond(FluidStateBoundaryCond::Read(HP, pRoot));
        std::unique_ptr<Geometry2D> pGeometry(Geometry2D::Read(pRoot, HP));

        return std::unique_ptr<LubricationGrooveMaster>{new LubricationGrooveMaster(iLabel, std::move(pGeometry), pBoundaryCond.release(), eType)};
    }

    LubricationGrooveMaster::LubricationGrooveMaster(integer iLabel,
                                                     std::unique_ptr<Geometry2D>&& pGeometry,
                                                     FluidStateBoundaryCond* pBoundaryCond,
                                                     enum Type type)
        :LubricationGroove(iLabel, std::move(pGeometry)),
         eType(type),
         iNumNodes(0),
         pBoundaryCond(pBoundaryCond)
    {

    }

    LubricationGrooveMaster::~LubricationGrooveMaster()
    {

    }

    void LubricationGrooveMaster::AddNode(Node2D*)
    {
        iNumNodes++;
    }

    integer LubricationGrooveMaster::iGetNumNodes() const
    {
        return iNumNodes;
    }

    FluidStateBoundaryCond* LubricationGrooveMaster::pReleaseBoundaryCond()
    {
        HYDRO_ASSERT(pBoundaryCond.owner());

        return pBoundaryCond.release();
    }

    LubricationGroove::Type LubricationGrooveMaster::GetType() const
    {
        return eType;
    }

    FluidStateBoundaryCond* LubricationGrooveMaster::pGetBoundaryCond() const
    {
        return pBoundaryCond;
    }

    LubricationGrooveSlave::LubricationGrooveSlave(LubricationGrooveMaster* pMaster, const Vector<doublereal, 2>& x)
        :LubricationGroove(pMaster->iGetLabel(),
                           pMaster->pGetGeometry()->Clone(x)),
         pMaster(pMaster)
    {

    }

    LubricationGrooveSlave::~LubricationGrooveSlave()
    {

    }

    FluidStateBoundaryCond* LubricationGrooveSlave::pGetBoundaryCond() const
    {
        return pMaster->pGetBoundaryCond();
    }

    void LubricationGrooveSlave::AddNode(Node2D* pNode)
    {
        pMaster->AddNode(pNode);
    }

    integer LubricationGrooveSlave::iGetNumNodes() const
    {
        return pMaster->iGetNumNodes();
    }

    FluidStateBoundaryCond* LubricationGrooveSlave::pReleaseBoundaryCond()
    {
        return nullptr;
    }

    LubricationGroove::Type LubricationGrooveSlave::GetType() const
    {
        return pMaster->GetType();
    }

    Pocket::Pocket(std::unique_ptr<Geometry2D>&& pGeometry)
        :pGeometry(std::move(pGeometry))
    {

    }

    Pocket::~Pocket()
    {

    }

    std::unique_ptr<Pocket> Pocket::Read(HydroRootElement* pRoot, MBDynParser& HP, const CylindricalBearing* pParent)
    {
        std::unique_ptr<Geometry2D> pGeometry(Geometry2D::Read(pRoot, HP));

        if (!HP.IsKeyWord("pocket" "height")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"pocket height\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("const")) {
            const doublereal dy = HP.GetReal();

            return std::unique_ptr<Pocket>{new ConstHeightPocket(std::move(pGeometry), dy)};
        } else if (HP.IsKeyWord("linear")) {
            Vector<doublereal, 2> x, z;
            Matrix<doublereal, 2, 2> Deltay;

            if (!HP.IsKeyWord("x")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"x\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            for (integer i = 1; i <= 2; ++i) {
                x(i) = HP.GetReal();

                if (i > 1 && x(i) <= x(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): x coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("z")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"z\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            for (integer i = 1; i <= 2; ++i) {
                z(i) = HP.GetReal();

                if (i > 1 && z(i) <= z(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): z coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("delta" "y")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"delta y\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            for (integer i = 1; i <= 2; ++i) {
                for (integer j = 1; j <= 2; ++j) {
                    Deltay(i, j) = HP.GetReal();
                }
            }

            return std::unique_ptr<Pocket>{new RectangularPocket{std::move(pGeometry), x, z, Deltay}};
        } else if (HP.IsKeyWord("surface" "grid")) {
            Vector<doublereal, DYNAMIC_SIZE> x, z;
            Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE> Deltay;

            if (!HP.IsKeyWord("x")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"x\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iGridX = HP.GetInt();

            if (iGridX < 2) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): at least two grid points are required in x direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            x.Resize(iGridX);

            for (integer i = 1; i <= iGridX; ++i) {
                x(i) = HP.GetReal();

                if (i > 1 && x(i) <= x(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): x coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("z")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"z\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iGridZ = HP.GetInt();

            if (iGridZ < 2) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): at least two grid points are required in z direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            z.Resize(iGridZ);

            for (integer i = 1; i <= iGridZ; ++i) {
                z(i) = HP.GetReal();

                if (i > 1 && z(i) <= z(i - 1)) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pRoot->GetLabel()
                                << "): z coordinates must be in ascending order at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (!HP.IsKeyWord("delta" "y")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"delta y\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            Deltay.Resize(iGridX, iGridZ);

            for (integer i = 1; i <= iGridX; ++i) {
                for (integer j = 1; j <= iGridZ; ++j) {
                    Deltay(i, j) = HP.GetReal();
                }
            }

            return std::unique_ptr<Pocket>{new SurfaceGrid(std::move(pGeometry), x, z, Deltay)};
        } else if (HP.IsKeyWord("helical" "groove")) {
            if (!HP.IsKeyWord("profile")) {
                silent_cerr("hydrodynamic plain bearing2(" << pRoot->GetLabel()
                            << "): keyword \"profile\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            std::array<std::unique_ptr<DriveCaller>, 2> rgProfile;

            for (index_type i = 0; i < 2; ++i) {
                rgProfile[i].reset(HP.GetDriveCaller());

                if (!rgProfile[i]->bIsDifferentiable()) {
                    silent_cerr("hydrodynamic plain bearing2(" << pRoot->GetLabel()
                                << "): drive caller is not differentiable at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            doublereal beta = HP.IsKeyWord("pitch" "angle") ? HP.GetReal() : 0.;

            Matrix<doublereal, 2, 2> R0;
            R0(1, 1) = cos(beta); R0(1, 2) = -sin(beta);
            R0(2, 1) = sin(beta); R0(2, 2) = cos(beta);

            Vector<doublereal, 2> x0;

            if (HP.IsKeyWord("origin")) {
                for (index_type i = 1; i <= 2; ++i) {
                    x0(i) = HP.GetReal();
                }
            }

            const integer K = HP.IsKeyWord("number" "of" "segments") ? HP.GetInt() : 1;

            if (K < 1) {
                silent_cerr("hydrodynamic plain bearing2(" << pRoot->GetLabel()
                            << "): at least one segment is required at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const double P = 2. * M_PI * pParent->dGetMeshRadius() * sin(beta) / K;

            return std::unique_ptr<Pocket>{new HelicalGroove{std::move(pGeometry), std::move(rgProfile), R0, x0, P}};
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keywords \"const\" or \"linear\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    ConstHeightPocket::ConstHeightPocket(std::unique_ptr<Geometry2D>&& pGeometry, doublereal dy)
        :Pocket(std::move(pGeometry)), dy(dy)
    {

    }

    ConstHeightPocket::~ConstHeightPocket()
    {

    }

    void ConstHeightPocket::GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const
    {
        Deltay = dy;
    }

    void ConstHeightPocket::GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const
    {
        Deltay = dy;
    }

    void ConstHeightPocket::GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const
    {
        dDeltay_dx = 0.;
    }

    void ConstHeightPocket::GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const
    {
        dDeltay_dx = 0.;
    }

    void ConstHeightPocket::GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const
    {
        dDeltay_dz = 0.;
    }

    void ConstHeightPocket::GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const
    {
        dDeltay_dz = 0.;
    }

    std::unique_ptr<Pocket> ConstHeightPocket::Clone(const Vector<doublereal, 2>& x) const
    {
        return std::unique_ptr<Pocket>{new ConstHeightPocket(pGetGeometry()->Clone(x), dy)};
    }

    RectangularPocket::RectangularPocket(std::unique_ptr<Geometry2D>&& pGeometry,
                                         const Vector<doublereal, 2>& x,
                                         const Vector<doublereal, 2>& z,
                                         const Matrix<doublereal, 2, 2>& Deltay)
        :Pocket(std::move(pGeometry)),
         x(x),
         z(z),
         f(Deltay)
    {
        dfi1_dx = (f(2, 1) - f(1, 1)) / (x(2) - x(1));
        dfi2_dx = (f(2, 2) - f(1, 2)) / (x(2) - x(1));
    }

    RectangularPocket::~RectangularPocket()
    {

    }

    void RectangularPocket::GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void RectangularPocket::GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void RectangularPocket::GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const
    {
        GetHeightDerXTpl(x, dDeltay_dx);
    }

    void RectangularPocket::GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const
    {
        GetHeightDerXTpl(x, dDeltay_dx);
    }

    void RectangularPocket::GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const
    {
        GetHeightDerZTpl(x, dDeltay_dz);
    }

    void RectangularPocket::GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const
    {
        GetHeightDerZTpl(x, dDeltay_dz);
    }

    std::unique_ptr<Pocket> RectangularPocket::Clone(const Vector<doublereal, 2>& xgc) const
    {
        return std::unique_ptr<Pocket>{new RectangularPocket(pGetGeometry()->Clone(xgc), x, z, f)};
    }

    template <typename T> inline
    void RectangularPocket::GetHeightTpl(const Vector<T, 2>& xci, T& Deltay) const
    {
        const T& xi = xci(1);
        const T& zi = xci(2);

        const T dx = (xi - x(1)) / (x(2) - x(1));
        const T dz = (zi - z(1)) / (z(2) - z(1));
        const T fi1 = (f(2, 1) - f(1, 1)) * dx + f(1, 1);
        const T fi2 = (f(2, 2) - f(1, 2)) * dx + f(1, 2);
        Deltay = (fi2 - fi1) * dz + fi1;
    }

    template <typename T> inline
    void RectangularPocket::GetHeightDerXTpl(const Vector<T, 2>& xci, T& dDeltay_dx) const
    {
        const T& zi = xci(2);
        const T dz = (zi - z(1)) / (z(2) - z(1));

        dDeltay_dx = (dfi2_dx - dfi1_dx) * dz + dfi1_dx;
    }

    template <typename T> inline
    void RectangularPocket::GetHeightDerZTpl(const Vector<T, 2>& xci, T& dDeltay_dz) const
    {
        const T& xi = xci(1);

        const T dx = (xi - x(1)) / (x(2) - x(1));
        const T fi1 = (f(2, 1) - f(1, 1)) * dx + f(1, 1);
        const T fi2 = (f(2, 2) - f(1, 2)) * dx + f(1, 2);

        dDeltay_dz = (fi2 - fi1) / (z(2) - z(1));
    }

    SurfaceGrid::SurfaceGrid(std::unique_ptr<Geometry2D>&& pGeometry,
                             const Vector<doublereal, DYNAMIC_SIZE>& x,
                             const Vector<doublereal, DYNAMIC_SIZE>& z,
                             const Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE>& f)
        :Pocket(std::move(pGeometry)),
         x(x),
         z(z),
         f(f)
    {
        HYDRO_ASSERT(x.iGetNumRows() >= 2);
        HYDRO_ASSERT(z.iGetNumRows() >= 2);
        HYDRO_ASSERT(x.iGetNumRows() == f.iGetNumRows());
        HYDRO_ASSERT(z.iGetNumRows() == f.iGetNumCols());
    }

    SurfaceGrid::~SurfaceGrid()
    {
    }

    void SurfaceGrid::GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void SurfaceGrid::GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void SurfaceGrid::GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const
    {
        GetHeightDerXTpl(x, dDeltay_dx);
    }

    void SurfaceGrid::GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const
    {
        GetHeightDerXTpl(x, dDeltay_dx);
    }

    void SurfaceGrid::GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const
    {
        GetHeightDerZTpl(x, dDeltay_dz);
    }

    void SurfaceGrid::GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const
    {
        GetHeightDerZTpl(x, dDeltay_dz);
    }

    std::unique_ptr<Pocket> SurfaceGrid::Clone(const Vector<doublereal, 2>& xc) const
    {
        return std::unique_ptr<Pocket>{new SurfaceGrid{pGetGeometry()->Clone(xc), x, z, f}};
    }

    template <typename T> inline
    index_type SurfaceGrid::iFindGridX(const Vector<T, 2>& xci) const
    {
        index_type ix = x.iGetNumRows() - 1;

        for (index_type i = 1; i < x.iGetNumRows(); ++i) {
            if (xci(1) <= x(i + 1)) {
                ix = i;
                break;
            }
        }

        return ix;
    }

    template <typename T> inline
    index_type SurfaceGrid::iFindGridZ(const Vector<T, 2>& xci) const
    {
        index_type iz = z.iGetNumRows() - 1;

        for (index_type j = 1; j < z.iGetNumRows(); ++j) {
            if (xci(2) <= z(j + 1)) {
                iz = j;
                break;
            }
        }

        return iz;
    }

    template <typename T> inline
    void SurfaceGrid::GetHeightTpl(const Vector<T, 2>& xci, T& Deltay) const
    {
        const index_type ix = iFindGridX(xci);
        const index_type iz = iFindGridZ(xci);
        const T dx = (xci(1) - x(ix)) / (x(ix + 1) - x(ix));
        const T dz = (xci(2) - z(iz)) / (z(iz + 1) - z(iz));
        const T fi1 = (f(ix + 1, iz) - f(ix, iz)) * dx + f(ix, iz);
        const T fi2 = (f(ix + 1, iz + 1) - f(ix, iz + 1)) * dx + f(ix, iz + 1);

        Deltay = (fi2 - fi1) * dz + fi1;
    }

    template <typename T> inline
    void SurfaceGrid::GetHeightDerXTpl(const Vector<T, 2>& xci, T& dDeltay_dx) const
    {
        const index_type ix = iFindGridX(xci);
        const index_type iz = iFindGridZ(xci);
        const T dz = (xci(2) - z(iz)) / (z(iz + 1) - z(iz));
        const doublereal dfi1_dx = (f(ix + 1, iz) - f(ix, iz)) / (x(ix + 1) - x(ix));
        const doublereal dfi2_dx = (f(ix + 1, iz + 1) - f(ix, iz + 1)) / (x(ix + 1) - x(ix));

        dDeltay_dx = (dfi2_dx - dfi1_dx) * dz + dfi1_dx;
    }

    template <typename T> inline
    void SurfaceGrid::GetHeightDerZTpl(const Vector<T, 2>& xci, T& dDeltay_dz) const
    {
        const index_type ix = iFindGridX(xci);
        const index_type iz = iFindGridZ(xci);
        const T dx = (xci(1) - x(ix)) / (x(ix + 1) - x(ix));
        const T fi1 = (f(ix + 1, iz) - f(ix, iz)) * dx + f(ix, iz);
        const T fi2 = (f(ix + 1, iz + 1) - f(ix, iz + 1)) * dx + f(ix, iz + 1);

        dDeltay_dz = (fi2 - fi1) / (z(iz + 1) - z(iz));
    }

    HelicalGroove::HelicalGroove(std::unique_ptr<Geometry2D>&& pGeometry,
                                 std::array<std::unique_ptr<DriveCaller>, 2>&& rgProfile,
                                 const Matrix<doublereal, 2, 2>& R0,
                                 const Vector<doublereal, 2>& x0,
                                 doublereal P)
        :Pocket(std::move(pGeometry)),
         rgProfile(std::move(rgProfile)),
         R0(R0),
         x0(x0),
         P(P)
    {
    }

    HelicalGroove::~HelicalGroove()
    {
    }

    void HelicalGroove::GetHeight(const Vector<doublereal, 2>& x, doublereal& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void HelicalGroove::GetHeight(const Vector<Gradient<0>, 2>& x, Gradient<0>& Deltay) const
    {
        GetHeightTpl(x, Deltay);
    }

    void HelicalGroove::GetHeightDerX(const Vector<doublereal, 2>& x, doublereal& dDeltay_dx) const
    {
        GetHeightDerTpl(x, dDeltay_dx, 1);
    }

    void HelicalGroove::GetHeightDerX(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dx) const
    {
        GetHeightDerTpl(x, dDeltay_dx, 1);
    }

    void HelicalGroove::GetHeightDerZ(const Vector<doublereal, 2>& x, doublereal& dDeltay_dz) const
    {
        GetHeightDerTpl(x, dDeltay_dz, 2);
    }

    void HelicalGroove::GetHeightDerZ(const Vector<Gradient<0>, 2>& x, Gradient<0>& dDeltay_dz) const
    {
        GetHeightDerTpl(x, dDeltay_dz, 2);
    }

    std::unique_ptr<Pocket> HelicalGroove::Clone(const Vector<doublereal, 2>& x) const
    {
        std::array<std::unique_ptr<DriveCaller>, 2> rgProfileClone;

        for (index_type i = 0; i < 2; ++i) {
            rgProfileClone[i].reset(rgProfile[i]->pCopy());
        }

        return std::unique_ptr<Pocket>{new HelicalGroove{pGetGeometry()->Clone(x), std::move(rgProfileClone), R0, x0, P}};
    }

    template <typename T> inline
    void HelicalGroove::RelativePosition(const Vector<T, 2>& xci, Vector<T, 2>& x) const
    {
        x = Transpose(R0) * Vector<T, 2>(xci - x0);

        const doublereal z = dGetValue(x(2));

        int K = z / P + copysign(0.5, z);

        x(2) = T(x(2) - K * P);
    }

    template <typename T> inline
    void HelicalGroove::GetHeightTpl(const Vector<T, 2>& xci, T& Deltay) const
    {
        Vector<T, 2> x, f;

        RelativePosition(xci, x);

        for (index_type i = 0; i < 2; ++i) {
            rgProfile[i]->dGet(x(i + 1), f(i + 1));
        }

        Deltay = f(1) + f(2);
    }

    template <typename T> inline
    void HelicalGroove::GetHeightDerTpl(const Vector<T, 2>& xci, T& dDeltay_dx, index_type iDirection) const
    {
        Vector<T, 2> x, df_dx;

        RelativePosition(xci, x);

        for (index_type i = 0; i < 2; ++i) {
            rgProfile[i]->dGetP(x(i + 1), df_dx(i + 1));
        }

        dDeltay_dx = Dot(R0.GetRow(iDirection), df_dx);
    }

    template <typename T>
    KinematicsBoundaryCond<T>::KinematicsBoundaryCond()
        :bContact(false),
         h(0.),
         dh_dt(0.),
         pasp(0.),
         Pfc(0.)
    {

    }

    template <typename T>
    void KinematicsBoundaryCond<T>::Update(HydroUpdatedNode* pNode,
                                           doublereal dCoef,
                                           enum FunctionCall func,
                                           LocalDofMap* pDofMap) {

        const HydroMesh* const pMesh = pNode->pGetMesh();
        const BearingGeometry* const pGeometry = pMesh->pGetGeometry();
        const ContactModel* const pContact = pNode->pGetContactModel();
        FrictionModel* const pFriction = pNode->pGetFrictionModel();

        pGeometry->GetBoundaryConditions(pNode, h, dh_dt, U1, U2, U, dCoef, func, pDofMap);

        // Attention: stiction states for LuGre friction have to be updated
        // also if there is no contact at the current time step
        bContact = pContact != nullptr;

        if (pContact) {
            pContact->GetContactPressure(h, pasp);
        } else {
            pasp = 0.;
        }

        if (pFriction) {
            const Vector<T, 2> U = U1 - U2;

            pFriction->GetFrictionForce(h, U, pasp, tauc_0);

            Pfc = Dot(U, tauc_0);
        } else {
            for (index_type i = 1; i <= tauc_0.iGetNumRows(); ++i) {
                tauc_0(i) = 0.;
            }

            Pfc = 0.;
        }
    }

    template <typename T>
    void KinematicsBoundaryCond<T>::GetClearance(T& h, LocalDofMap* pDofMap) const
    {
        Copy(h, this->h, pDofMap);
    }

    template <typename T>
    void KinematicsBoundaryCond<T>::GetClearanceDerTime(T& dh_dt, LocalDofMap* pDofMap) const
    {
        Copy(dh_dt, this->dh_dt, pDofMap);
    }

    template <typename T>
    void KinematicsBoundaryCond<T>::GetVelocity(Vector<T, 2>& U1, Vector<T, 2>& U2, LocalDofMap* pDofMap) const
    {
        U1.Copy(this->U1, pDofMap);
        U2.Copy(this->U2, pDofMap);
    }

    template <typename T>
    void KinematicsBoundaryCond<T>::GetHydraulicVelocity(Vector<T, 2>& U, LocalDofMap* pDofMap) const
    {
        U.Copy(this->U, pDofMap);
    }

    template <typename T>
    bool KinematicsBoundaryCond<T>::GetContactPressure(T& pasp, LocalDofMap* pDofMap) const
    {
        Copy(pasp, this->pasp, pDofMap);

        return bContact;
    }

    template <typename T>
    bool KinematicsBoundaryCond<T>::GetContactStress(Vector<T, 2>& tauc_0, LocalDofMap* pDofMap) const
    {
        tauc_0.Copy(this->tauc_0, pDofMap);

        return bContact;
    }

    template <typename T>
    bool KinematicsBoundaryCond<T>::GetContactFrictionLossDens(T& Pfc, LocalDofMap* pDofMap) const
    {
        Copy(Pfc, this->Pfc, pDofMap);

        return bContact;
    }

    ThermWallBoundCond::ThermWallBoundCond()
    {
        std::fill(rgNodes.begin(), rgNodes.end(), nullptr);
    }

    ThermWallBoundCond::~ThermWallBoundCond()
    {
        // Nodes owned by DataManager
    }

    void ThermWallBoundCond::ParseInput(DataManager* pDM, MBDynParser& HP, const HydroRootElement* pParent)
    {
        static const char szWallName[2][8] = {"shaft", "bearing"};

        for (index_type i = 0; i < 2; ++i) {
            if (!HP.IsKeyWord(szWallName[i])) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pParent->GetLabel()
                            << "): keyword \""
                            << szWallName[i]
                            << "\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            rgNodes[i] = pDM->ReadNode<ThermalNode, Node::THERMAL>(HP);
        }
    }

    HydroDofOwner::HydroDofOwner()
    {
        std::fill(rgOffsetIndex.begin(), rgOffsetIndex.end(), grad::UNKNOWN_FUNC);
    }

    HydroDofOwner::~HydroDofOwner()
    {

    }

    integer HydroDofOwner::iGetOffsetIndex(grad::FunctionCall eFunc) const
    {
        return rgOffsetIndex[iFuncCallToIndex(eFunc)];
    }

    void HydroDofOwner::SetOffsetIndex(integer iOffset, grad::FunctionCall eFunc)
    {
        rgOffsetIndex[iFuncCallToIndex(eFunc)] = iOffset;
    }

    Node2D::Node2D(integer iNodeNo,
                   const Vector<doublereal, 2>& x,
                   HydroMesh* pParent,
                   integer iNodeFlags)
        :iNodeNo(iNodeNo),
         iNodeFlags(iNodeFlags),
         x(x),
         pMesh(pParent),
         pFluid(pParent->pGetParent()->pGetFluid())
    {
#if HYDRO_DEBUG > 0
        switch (iGetNodeFlags() & PHYSICS_MASK)
        {
        case HYDRAULIC_NODE:
        case THERMAL_NODE:
        case FLUX_NODE_X:
        case FLUX_NODE_Z:
            break;
        default:
            HYDRO_ASSERT(0);
        }

        switch (iGetNodeFlags() & LOCATION_MASK)
        {
        case CORNER_NODE:
        case CENTRAL_NODE:
            break;
        default:
            HYDRO_ASSERT(0);
        }

        switch (iGetNodeFlags() & ACTIVE_MASK)
        {
        case ACTIVE_NODE:
        case PASSIVE_NODE:
        case COUPLED_NODE:
        case COMPUTED_NODE:
            break;
        default:
            HYDRO_ASSERT(0);
        }

        switch (iGetNodeFlags() & MASTER_SLAVE_MASK)
        {
        case MASTER_NODE:
        case SLAVE_NODE:
            break;
        default:
            HYDRO_ASSERT(0);
        }

        switch (iGetNodeFlags() & COMPRESSIBLE_MASK)
        {
        case 0u:
        case COMPRESSIBLE_NODE:
        case INCOMPRESSIBLE_NODE:
            break;
        default:
            HYDRO_ASSERT(0);
        }
#endif
    }

    Node2D::~Node2D()
    {
    }

    integer Node2D::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return 0;
    }

    void Node2D::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
    }

    void Node2D::DofUpdate(VectorHandler& X, VectorHandler& XP)
    {
    }

    void Node2D::AfterConvergence(const VectorHandler& X,
                                  const VectorHandler& XP)
    {
    }

    ThermoHydrNode::ThermoHydrNode(integer iNodeNo,
                                   const Vector<doublereal, 2>& x,
                                   HydroMesh* pMesh,
                                   integer iNodeFlags)
        :Node2D(iNodeNo, x, pMesh, THERMAL_NODE | CORNER_NODE | iNodeFlags)
    {
    }

    ThermoHydrNode::~ThermoHydrNode()
    {
    }

    bool ThermoHydrNode::bGetPrivateData(HydroRootBase::PrivateDataType eType, doublereal& dPrivData) const
    {
        return false;
    }

    void ThermoHydrNode::Output(std::ostream& os, unsigned uOutputFlags) const
    {
        if (uOutputFlags & HydroRootElement::OUTPUT_CENT_TEMPERATURE) {
            doublereal T, dT_dt;

            GetTemperature(T);
            GetTemperatureDerTime(dT_dt);

            os << T << ' ' << dT_dt << ' ';
        }
    }

    ThermalActiveNode::ThermalActiveNode(integer iNodeNo,
                                         const Vector<doublereal, 2>& x,
                                         HydroMesh* pParent,
                                         doublereal T0,
                                         bool bDoInitAss)
        :ThermoHydrNode(iNodeNo, x, pParent, ACTIVE_NODE | MASTER_NODE),
         eCurrFunc(grad::INITIAL_ASS_FLAG),
         T(T0),
         dT_dt(0.),
         s(pParent->pGetParent()->dGetScale(HydroRootElement::SCALE_TEMPERATURE_DOF)),
         bDoInitAss(bDoInitAss)
    {
    }

    ThermalActiveNode::~ThermalActiveNode()
    {
    }

    void ThermalActiveNode::GetTemperature(doublereal& T, LocalDofMap*, doublereal) const
    {
        T = this->T;
    }

    void ThermalActiveNode::GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (eCurrFunc & grad::REGULAR_FLAG) {
            T.SetValuePreserve(this->T);
            T.DerivativeResizeReset(pDofMap,
                                    iGetFirstDofIndex(eCurrFunc),
                                    MapVectorBase::GLOBAL,
                                    -dCoef * s);
        } else {
            T.SetValue(this->T);
        }
    }

    void ThermalActiveNode::GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap*, doublereal) const
    {
        dT_dt = this->dT_dt;
    }

    void ThermalActiveNode::GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (bDoInitAss || (eCurrFunc & grad::REGULAR_FLAG)) {
            dT_dt.SetValuePreserve(this->dT_dt);
            dT_dt.DerivativeResizeReset(pDofMap,
                                        iGetFirstDofIndex(eCurrFunc),
                                        MapVectorBase::GLOBAL,
                                        -s);
        } else {
            dT_dt.SetValue(this->dT_dt);
        }
    }

    integer ThermalActiveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return iGetFirstDofIndex(eFunc);
    }

    integer ThermalActiveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(iGetOffsetIndex(eFunc) != UNKNOWN_OFFSET);
        HYDRO_ASSERT(unsigned(iGetOffsetIndex(eFunc)) <= pGetMesh()->pGetParent()->iGetNumDof());
        HYDRO_ASSERT((eFunc & grad::REGULAR_FLAG) || bDoInitAss);

        return pGetMesh()->pGetParent()->iGetFirstIndex() + iGetOffsetIndex(eFunc);
    }

    integer ThermalActiveNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        if (eFunc & grad::REGULAR_FLAG) {
            return iGetNumDof();
        } else {
            return iGetInitialNumDof();
        }
    }

    void
    ThermalActiveNode::Update(const VectorHandler& XCurr,
                              const VectorHandler& XPrimeCurr,
                              doublereal dCoef,
                              enum FunctionCall func,
                              LocalDofMap* pDofMap)
    {
        HYDRO_ASSERT(iGetOffsetIndex(func) > 0);

        const integer iIndex = iGetFirstDofIndex(func);

        HYDRO_ASSERT(iIndex > 0);
        HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

        if (func & grad::REGULAR_FLAG) {
            T = s * XCurr(iIndex);
            dT_dt = s * XPrimeCurr(iIndex);
        } else if (bDoInitAss) {
            dT_dt = s * XCurr(iIndex);
        }
    }

    void ThermalActiveNode::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);

        eCurrFunc = grad::REGULAR_FLAG;

        HYDRO_ASSERT(iGetOffsetIndex(eCurrFunc) > 0);

        const integer iIndex = iGetFirstDofIndex(eCurrFunc);

        HYDRO_ASSERT(iIndex > 0);
        HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

        XCurr.PutCoef(iIndex, T / s);
        XPrimeCurr.PutCoef(iIndex, dT_dt / s);
    }

    unsigned int ThermalActiveNode::iGetNumDof(void) const
    {
        return 1u;
    }

    unsigned int ThermalActiveNode::iGetInitialNumDof(void) const
    {
        return bDoInitAss ? 1u : 0u;
    }

    DofOrder::Order ThermalActiveNode::GetDofType(unsigned int i) const
    {
        return DofOrder::DIFFERENTIAL;
    }

    DofOrder::Order ThermalActiveNode::GetEqType(unsigned int i) const
    {
        return GetDofType(i);
    }

    std::ostream&
    ThermalActiveNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        integer iIndex = iGetFirstDofIndex(bInitial
                                           ? grad::INITIAL_ASS_FLAG
                                           : grad::REGULAR_FLAG);

        if (bInitial) {
            if (bDoInitAss) {
                out << prefix << iIndex << ": dT" << iGetNodeNumber() + 1 << "/dt" << std::endl;
            }
        } else {
            out << prefix << iIndex << ": T" << iGetNodeNumber() + 1 << std::endl;
        }

        return out;
    }

    std::ostream&
    ThermalActiveNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iIndex = iGetFirstDofIndex(bInitial
                                               ? grad::INITIAL_ASS_FLAG
                                               : grad::REGULAR_FLAG);

            out << prefix << iIndex << ":  energy balance node " << iGetNodeNumber() + 1 << std::endl;
        }

        return out;
    }

    ThermalCoupledNode::ThermalCoupledNode(integer iNodeNo,
                                           const Vector<doublereal, 2>& x,
                                           HydroMesh* pMesh,
                                           ThermalNode* pExtThermNode)
        :ThermoHydrNode(iNodeNo, x, pMesh, COUPLED_NODE | MASTER_NODE),
         pExtThermNode(pExtThermNode),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {
        HYDRO_ASSERT(pExtThermNode != nullptr);
    }

    ThermalCoupledNode::~ThermalCoupledNode()
    {
    }

    void ThermalCoupledNode::GetTemperature(doublereal& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pExtThermNode->GetX(T, dCoef, eCurrFunc, pDofMap);
    }

    void ThermalCoupledNode::GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pExtThermNode->GetX(T, dCoef, eCurrFunc, pDofMap);
    }

    void ThermalCoupledNode::GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pExtThermNode->GetXPrime(dT_dt, dCoef, eCurrFunc, pDofMap);
    }

    void ThermalCoupledNode::GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pExtThermNode->GetXPrime(dT_dt, dCoef, eCurrFunc, pDofMap);
    }

    integer ThermalCoupledNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(eFunc & grad::REGULAR_FLAG);

        return pExtThermNode->iGetFirstRowIndex() + 1;
    }

    integer ThermalCoupledNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(eFunc & grad::REGULAR_FLAG);

        return pExtThermNode->iGetFirstColIndex() + 1;
    }

    integer ThermalCoupledNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return eFunc & grad::REGULAR_FLAG ? 1 : 0;
    }

    void
    ThermalCoupledNode::Update(const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               doublereal dCoef,
                               enum FunctionCall func,
                               LocalDofMap* pDofMap)
    {
        if (func & grad::REGULAR_FLAG) {
            eCurrFunc = grad::REGULAR_FLAG;
        }
    }

    ThermalInletNode::ThermalInletNode(integer iNodeNo,
                                       const Vector<doublereal, 2>& x,
                                       HydroMesh* pParent,
                                       ThermalNode* pExtThermNode,
                                       bool bDoInitAss)
        :ThermalActiveNode(iNodeNo, x, pParent, pExtThermNode->dGetX(), bDoInitAss),
         oInletNode(iNodeNo, x, pParent, pExtThermNode)
    {

    }

    ThermalInletNode::~ThermalInletNode()
    {
    }

    void
    ThermalInletNode::Update(const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr,
                             doublereal dCoef,
                             enum FunctionCall func,
                             LocalDofMap* pDofMap)
    {
        ThermalActiveNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
        oInletNode.Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
    }

    integer ThermalInletNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return ThermalActiveNode::iGetNumColsWorkSpace(eFunc)
            + oInletNode.iGetNumColsWorkSpace(eFunc);
    }

    ThermalPassiveNode::ThermalPassiveNode(integer iNodeNo,
                                           const Vector<doublereal, 2>& x,
                                           HydroMesh* pParent,
                                           const FluidStateBoundaryCond* pBoundCond)
        :ThermoHydrNode(iNodeNo, x, pParent, PASSIVE_NODE | MASTER_NODE),
         pBoundCond(pBoundCond)
    {
        HYDRO_ASSERT(pBoundCond->bIncludeNode(GetNodePhysics()));
    }

    ThermalPassiveNode::~ThermalPassiveNode()
    {
    }

    void ThermalPassiveNode::GetTemperature(doublereal& T, LocalDofMap*, doublereal) const
    {
        T = pBoundCond->dGetTemperature();
    }

    void ThermalPassiveNode::GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        T.SetValue(pBoundCond->dGetTemperature());
    }

    void ThermalPassiveNode::GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap*, doublereal) const
    {
        dT_dt = pBoundCond->dGetTemperatureDerTime();
    }

    void ThermalPassiveNode::GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        dT_dt.SetValue(pBoundCond->dGetTemperatureDerTime());
    }

    integer ThermalPassiveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    integer ThermalPassiveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    void
    ThermalPassiveNode::Update(const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               doublereal dCoef,
                               enum FunctionCall func,
                               LocalDofMap* pDofMap)
    {
    }

    ThermalSlaveNode::ThermalSlaveNode(integer iNodeNo,
                                       const Vector<doublereal, 2>& x,
                                       ThermoHydrNode* pMasterNode)
        :ThermoHydrNode(iNodeNo, x, pMasterNode->pGetMesh(), (pMasterNode->iGetNodeFlags() & ~MASTER_NODE) | SLAVE_NODE),
         pMasterNode(pMasterNode)
    {
    }

    ThermalSlaveNode::~ThermalSlaveNode()
    {
    }

    integer ThermalSlaveNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetNumColsWorkSpace(eFunc);
    }

    void ThermalSlaveNode::GetTemperature(doublereal& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetTemperature(T, pDofMap, dCoef);
    }

    void ThermalSlaveNode::GetTemperature(Gradient<0>& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetTemperature(T, pDofMap, dCoef);
    }

    void ThermalSlaveNode::GetTemperatureDerTime(doublereal& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
    }

    void ThermalSlaveNode::GetTemperatureDerTime(Gradient<0>& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
    }

    integer ThermalSlaveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetFirstEquationIndex(eFunc);
    }

    integer ThermalSlaveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetFirstDofIndex(eFunc);
    }

    void
    ThermalSlaveNode::Update(const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr,
                             doublereal dCoef,
                             enum FunctionCall func,
                             LocalDofMap* pDofMap)
    {
        NO_OP;
    }

    FluxNode::FluxNode(integer iNodeNo,
                       HydroMesh* pMesh,
                       const std::array<const HydroNode*, iNumNodes>& rgNodes,
                       PressureSource ePressSrc,
                       NodeDataReq eNodeDataReq)
        :Node2D(iNodeNo,
                (rgNodes[0]->GetPosition2D() + rgNodes[1]->GetPosition2D()) * 0.5,
                pMesh,
                (iDirectionFromNodes(rgNodes) == 1 ? FLUX_NODE_X : FLUX_NODE_Z) |
                CENTRAL_NODE |
                COMPUTED_NODE |
                MASTER_NODE),
         iDirection(iDirectionFromNodes(rgNodes)),
         du(0.),
         rgNodes(rgNodes),
         ePressSource(ePressSrc),
         uNodeDataReq(eNodeDataReq)
    {
#if HYDRO_DEBUG > 0
        HYDRO_ASSERT(ePressSource >= 0);
        HYDRO_ASSERT(ePressSource < iNumPressSources);
        HYDRO_ASSERT(rgNodes[0]->pGetMesh() == pGetMesh());
        HYDRO_ASSERT(rgNodes[1]->pGetMesh() == pGetMesh());
#endif
        du = pMesh->pGetGeometry()->dGetNodeDistance2D(rgNodes[iNodeDown], rgNodes[iNodeUp], iDirection);
    }

    FluxNode::~FluxNode()
    {
    }

    index_type FluxNode::iDirectionFromNodes(const std::array<const HydroNode*, iNumNodes>& rgNodes)
    {
        index_type iDirection = -1;

        for (index_type i = 1; i <= 2; ++i) {
            if (rgNodes[0]->GetPosition2D()(i) != rgNodes[1]->GetPosition2D()(i)) {
                iDirection = i;
                break;
            }
        }

        HYDRO_ASSERT(iDirection != -1);

        return iDirection;
    }

    void FluxNode::GetEnergyBalance(doublereal& Qu, LocalDofMap*) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_THERMAL);

        Qu = oNode.Qu;
    }

    void FluxNode::GetEnergyBalance(Gradient<0>& Qu, LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_THERMAL);

        Qu.Copy(oNode_grad.Qu, pDofMap);
    }

    void FluxNode::GetDissipationFactors(doublereal& A0, doublereal& Ah,  doublereal& Ac, LocalDofMap*) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_THERMAL_WALL);

        A0 = oNode.A0;
        Ah = oNode.Ah;
        Ac = oNode.Ac;
    }

    void FluxNode::GetDissipationFactors(Gradient<0>& A0, Gradient<0>& Ah, Gradient<0>& Ac, LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_THERMAL_WALL);

        A0.Copy(oNode_grad.A0, pDofMap);
        Ah.Copy(oNode_grad.Ah, pDofMap);
        Ac.Copy(oNode_grad.Ac, pDofMap);
    }

    void FluxNode::GetVolumeFluxDens(doublereal& qu, LocalDofMap*, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        qu = rgFlux[ePressSrc].qu;
    }

    void FluxNode::GetVolumeFluxDens(Gradient<0>& qu, LocalDofMap* pDofMap, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        qu.Copy(rgFlux_grad[ePressSrc].qu, pDofMap);
    }

    void FluxNode::GetVelocityAvg(doublereal& wu, LocalDofMap*, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        wu = rgFlux[ePressSrc].wu;
    }

    void FluxNode::GetVelocityAvg(Gradient<0>& wu, LocalDofMap* pDofMap, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        wu.Copy(rgFlux_grad[ePressSrc].wu, pDofMap);
    }

    void FluxNode::GetMassFluxDens(doublereal& mdotu, LocalDofMap*, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        mdotu = rgFlux[ePressSrc].mdotu;
    }

    void FluxNode::GetMassFluxDens(Gradient<0>& mdotu, LocalDofMap* pDofMap, PressureSource ePressSrc) const
    {
        HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);
        HYDRO_ASSERT(ePressSrc <= ePressSource);

        mdotu.Copy(rgFlux_grad[ePressSrc].mdotu, pDofMap);
    }

    void FluxNode::RequestPressureSource(PressureSource ePressSrcReq)
    {
        HYDRO_ASSERT(ePressSrcReq >= 0);
        HYDRO_ASSERT(ePressSrcReq < iNumPressSources);

        if (ePressSrcReq > ePressSource) {
            ePressSource = ePressSrcReq;
        }
    }

    void FluxNode::RequestNodeData(NodeDataReq eFlag)
    {
        uNodeDataReq |= eFlag;
    }

    void
    FluxNode::Update(const VectorHandler& XCurr,
                     const VectorHandler& XPrimeCurr,
                     doublereal dCoef,
                     enum FunctionCall func,
                     LocalDofMap*)
    {
        if (uNodeDataReq & ND_HYDRAULIC) {
            switch (func) {
            case REGULAR_RES:
            case INITIAL_ASS_RES:
            case INITIAL_DER_RES:
                UpdateTpl(oNode,
                          rgFlux,
                          dCoef,
                          func,
                          &oDofMap);
                break;

            case REGULAR_JAC:
            case INITIAL_ASS_JAC:
            case INITIAL_DER_JAC:
                UpdateTpl(oNode_grad,
                          rgFlux_grad,
                          dCoef,
                          func,
                          &oDofMap);
                break;

            default:
                HYDRO_ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }
    }

    integer FluxNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);

        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }


    integer FluxNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);

        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    bool
    FluxNode::bGetPrivateData(HydroRootBase::PrivateDataType eType,
                              doublereal& dPrivData) const
    {
        return false;
    }

    void
    FluxNode::Output(std::ostream& os, unsigned uOutputFlags) const
    {
        if (((uOutputFlags & HydroRootElement::OUTPUT_VOLUME_FLUX_X) && iDirection == 1) ||
            ((uOutputFlags & HydroRootElement::OUTPUT_VOLUME_FLUX_Z) && iDirection == 2)) {
            HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);

            doublereal qu;

            GetVolumeFluxDens(qu);

            os << qu << ' ';
        }

        if (((uOutputFlags & HydroRootElement::OUTPUT_MASS_FLUX_X) && iDirection == 1) ||
            ((uOutputFlags & HydroRootElement::OUTPUT_MASS_FLUX_Z) && iDirection == 2)) {
            HYDRO_ASSERT(uNodeDataReq & ND_HYDRAULIC);

            doublereal mdotu;

            GetMassFluxDens(mdotu);

            os << mdotu << ' ';
        }

        if (((uOutputFlags & HydroRootElement::OUTPUT_HEAT_FLUX_X) && iDirection == 1) ||
            ((uOutputFlags & HydroRootElement::OUTPUT_HEAT_FLUX_Z) && iDirection == 2)) {
            HYDRO_ASSERT(uNodeDataReq & ND_THERMAL);

            doublereal Qu;

            GetEnergyBalance(Qu);

            os << Qu << ' ';
        }
    }

    template <typename G>
    void FluxNode::UpdateTpl(NodeData<G>& oNode,
                             std::array<FluxData<G>, iNumPressSources>& rgFlux,
                             doublereal dCoef,
                             FunctionCall func,
                             LocalDofMap* pDofMap) const
    {
        std::array<G, iNumNodes> pi, hi, etai;
        std::array<Vector<G, 2>, iNumNodes> Ui;
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        for (index_type i = 0; i < iNumNodes; ++i) {
            rgNodes[i]->GetClearance(hi[i], pDofMap);
            pGeometry->GetNonNegativeClearance(hi[i], hi[i]);
            rgNodes[i]->GetViscosity(etai[i], pDofMap, dCoef);
            rgNodes[i]->GetHydraulicVelocity(Ui[i], pDofMap);
        }

        const G h = 0.5 * (hi[iNodeDown] + hi[iNodeUp]);
        const G eta = 0.5 * (etai[iNodeDown] + etai[iNodeUp]);
        const G U = 0.5 * (Ui[iNodeDown](iDirection) + Ui[iNodeUp](iDirection));

        for (index_type j = 0; j <= ePressSource; ++j) {
            for (index_type i = 0; i < iNumNodes; ++i) {
                switch (j) {
                case PRESSURE_FROM_NODE:
                    rgNodes[i]->GetPressure(pi[i], pDofMap, dCoef);
                    break;

                case PRESSURE_FROM_MESH:
                    pGetMesh()->GetPressure(rgNodes[i], pi[i], pDofMap, dCoef);
                    break;
                };
            }

            const G dp_du = (pi[iNodeUp] - pi[iNodeDown]) / du;

            const G a0 = h * h / (12. * eta) * dp_du;

            rgFlux[j].wu = U - a0;
            rgFlux[j].qu = h * rgFlux[j].wu;

            const index_type iUpwindu = rgFlux[j].qu >= 0. ? iNodeDown : iNodeUp;
            G rho;

            rgNodes[iUpwindu]->GetDensity(rho, pDofMap, dCoef);

            rgFlux[j].mdotu = rho * rgFlux[j].qu;

            if ((uNodeDataReq & ND_THERMAL) && j == PRESSURE_FROM_NODE) {
                std::array<G, iNumNodes> Ti;
                std::array<Vector<G, 2>, iNumNodes> U1i, U2i;

                for (index_type i = 0; i < iNumNodes; ++i) {
                    rgNodes[i]->GetTemperature(Ti[i], pDofMap, dCoef);
                    rgNodes[i]->GetVelocity(U1i[i], U2i[i], pDofMap);
                }

                const G dU = 0.5 * (U1i[iNodeDown](iDirection) + U1i[iNodeUp](iDirection)
                                    - U2i[iNodeDown](iDirection) - U2i[iNodeUp](iDirection));

                doublereal beta = 0.;

                if (pi[iNodeUp] > pGetFluid()->dGetRefPressure() &&
                    pi[iNodeDown] > pGetFluid()->dGetRefPressure()) {
                    doublereal rhoc, drhoc_dT;

                    pGetFluid()->GetDensity(pGetFluid()->dGetRefPressure(),
                                            pGetFluid()->dGetRefTemperature(),
                                            rhoc,
                                            nullptr,
                                            &drhoc_dT);

                    beta = -drhoc_dT / rhoc;
                }

                G cp;

                pGetFluid()->GetSpecificHeat(pi[iUpwindu],
                                             Ti[iUpwindu],
                                             rho,
                                             cp,
                                             HydroFluid::SPEC_HEAT_TRUE);

                oNode.Qu = rgFlux[j].qu * (beta * Ti[iUpwindu] * dp_du
                                           - rho * cp * (Ti[iNodeUp] - Ti[iNodeDown]) / du)
                    + h * a0 * dp_du + eta * dU * dU / h;

                if (uNodeDataReq & ND_THERMAL_WALL) {
                    std::array<G, iNumNodes> Pfci;

                    for (index_type i = 0; i < iNumNodes; ++i) {
                        rgNodes[i]->GetContactFrictionLossDens(Pfci[i], pDofMap);
                    }

                    const G a1 = 0.5 * h / eta * dp_du;

                    const G Psi0 = dU / h - a1;
                    const G Psih = dU / h + a1;

                    oNode.A0 = eta * Psi0 * Psi0;
                    oNode.Ah = eta * Psih * Psih;
                    oNode.Ac = 0.5 * (Pfci[iNodeUp] + Pfci[iNodeDown]) / h;
                }
            }
        }
    }

    HydroNode::HydroNode(integer iNodeNo,
                         const Vector<doublereal, 2>& x,
                         HydroMesh* pParent,
                         integer iNodeFlags)
        :Node2D(iNodeNo, x, pParent, HYDRAULIC_NODE | CORNER_NODE | iNodeFlags),
         pThermalNode(nullptr)
    {

    }

    HydroNode::~HydroNode()
    {

    }

    HydroFluid::CavitationState HydroNode::GetCavitationState() const
    {
        return HydroFluid::FULL_FILM_REGION;
    }

    template <typename G>
    void HydroNode::GetTemperature(G& T, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (pThermalNode) {
            pThermalNode->GetTemperature(T, pDofMap, dCoef);
        } else {
            T = pGetFluid()->dGetRefTemperature();
        }
    }

    template <typename G>
    void HydroNode::GetTemperatureDerTime(G& dT_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (pThermalNode) {
            pThermalNode->GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
        } else {
            dT_dt = 0.;
        }
    }

    template <typename G>
    void HydroNode::GetViscosity(G& eta, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        G rho, T;

        GetDensity(rho, pDofMap, dCoef);
        GetTemperature(T, pDofMap, dCoef);

        pGetFluid()->GetViscosity(rho, T, eta);
    }

    doublereal HydroNode::dGetClearance(const FluidStateBoundaryCond* const pBoundCond, doublereal* const dh_dt) const
    {
        doublereal h;

        if (pBoundCond->bNeedClearance()) {
            GetClearance(h);

            if (dh_dt) {
                GetClearanceDerTime(*dh_dt);
            }

            pGetMesh()->pGetGeometry()->GetNonNegativeClearance(h, h, dh_dt, dh_dt);
        } else {
            h = 0.;

            if (dh_dt) {
                *dh_dt = 0.;
            }
        }

        return h;
    }


    bool
    HydroNode::bGetPrivateData(HydroRootBase::PrivateDataType eType,
                               doublereal& dPrivData) const
    {
        switch (eType) {
        case HydroRootBase::PD_CLEARANCE:
            GetClearance(dPrivData);
            return true;

        case HydroRootBase::PD_TOTAL_DEFORMATION: {
            doublereal dw_dt;
            GetRadialDeformation(dPrivData, dw_dt);
        } return true;

        case HydroRootBase::PD_PRESSURE:
            pGetMesh()->GetPressure(this, dPrivData);
            return true;

        case HydroRootBase::PD_CONT_PRESSURE:
            GetContactPressure(dPrivData);
            return true;

        case HydroRootBase::PD_DENSITY:
            GetDensity(dPrivData);
            return true;

        case HydroRootBase::PD_TEMPERATURE:
            GetTemperature(dPrivData);
            return true;

        default:
            return false;
        }
    }

    void
    HydroNode::Output(std::ostream& os, unsigned uOutputFlags) const
    {
        if (uOutputFlags & HydroRootElement::OUTPUT_PRESSURE) {
            doublereal p;

            pGetMesh()->GetPressure(this, p);

            os << p << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_CONT_PRESSURE) {
            doublereal pasp;
            GetContactPressure(pasp);
            os << pasp << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_DENSITY) {
            doublereal rho;

            GetDensity(rho);

            os << rho << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_CLEARANCE) {
            doublereal h;

            GetClearance(h);

            os << h << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_CLEARANCE_DER) {
            doublereal dh_dt;

            GetClearanceDerTime(dh_dt);

            os << dh_dt << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_VELOCITY) {
            Vector<doublereal, 2> U1, U2;

            GetVelocity(U1, U2);

            for (integer i = 1; i <= 2; ++i) {
                os << U1(i) << ' ';
            }

            for (integer i = 1; i <= 2; ++i) {
                os << U2(i) << ' ';
            }
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_STRESS) {
            doublereal tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h;

            GetStress(tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);

            os << tau_xy_0 << ' ' << tau_yz_0 << ' ' << tau_xy_h << ' ' << tau_yz_h << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_CONT_STRESS) {
            Vector<doublereal, 2> tauc_0;

            GetContactStress(tauc_0);

            for (index_type i = 1; i <= tauc_0.iGetNumRows(); ++i) {
                os << tauc_0(i) << ' ';
            }
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_TOTAL_DEFORMATION) {
            doublereal w, dw_dt;

            GetRadialDeformation(w, dw_dt);

            os << w << ' ' << dw_dt << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_DEFORMATION1) {
            doublereal w1;

            GetRadialDeformation1(w1);

            os << w1 << ' ';
        }

        if (uOutputFlags & HydroRootElement::OUTPUT_DEFORMATION2) {
            doublereal w2;

            GetRadialDeformation2(w2);

            os << w2 << ' ';
        }
    }

    HydroSlaveNode::HydroSlaveNode(integer iNodeNo,
                                   const Vector<doublereal, 2>& x,
                                   HydroMesh* pMesh,
                                   HydroNode* pMasterNode)
        :HydroNode(iNodeNo, x, pMesh, (pMasterNode->iGetNodeFlags() & ~MASTER_NODE) | SLAVE_NODE),
         pMasterNode(pMasterNode)
    {
        HYDRO_ASSERT(pMasterNode->bIsNodeType(MASTER_NODE));
    }

    HydroSlaveNode::~HydroSlaveNode()
    {

    }

    const Vector<doublereal, 3>&
    HydroSlaveNode::GetPosition3D() const
    {
        return pMasterNode->GetPosition3D();
    }

    const Matrix<doublereal, 3, 3>&
    HydroSlaveNode::GetTangentCoordSys() const
    {
        return pMasterNode->GetTangentCoordSys();
    }

    integer HydroSlaveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetFirstEquationIndex(eFunc);
    }

    integer HydroSlaveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetFirstDofIndex(eFunc);
    }

    integer HydroSlaveNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return pMasterNode->iGetNumColsWorkSpace(eFunc);
    }

    void
    HydroSlaveNode::Update(const VectorHandler& XCurr,
                           const VectorHandler& XPrimeCurr,
                           doublereal dCoef,
                           enum FunctionCall func,
                           LocalDofMap* pDofMap)
    {
        NO_OP;
    }

    void HydroSlaveNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const
    {
        pMasterNode->GetPressure(p);
    }

    void HydroSlaveNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetPressure(p, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetPressureDerTime(dp_dt, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetPressureDerTime(dp_dt, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetDensity(doublereal& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensity(rho, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensity(rho, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensityDerTime(drho_dt, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensityDerTime(drho_dt, pDofMap, dCoef);
    }

    void HydroSlaveNode::GetStress(doublereal& tau_xy_0, doublereal& tau_yz_0, doublereal& tau_xy_h, doublereal& tau_yz_h) const
    {
        pMasterNode->GetStress(tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);
    }

    void HydroSlaveNode::SetStress(doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h)
    {
        pMasterNode->SetStress(tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);
    }

    bool HydroSlaveNode::GetContactPressure(doublereal& pasp, LocalDofMap* pDofMap) const
    {
        return pMasterNode->GetContactPressure(pasp, pDofMap);
    }

    bool HydroSlaveNode::GetContactPressure(Gradient<0>& pasp, LocalDofMap* pDofMap) const
    {
        return pMasterNode->GetContactPressure(pasp, pDofMap);
    }

    void HydroSlaveNode::GetContactStress(Vector<doublereal, 2>& tauc_0, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetContactStress(tauc_0, pDofMap);
    }

    void HydroSlaveNode::GetContactStress(Vector<Gradient<0>, 2>& tauc_0, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetContactStress(tauc_0, pDofMap);
    }

    void HydroSlaveNode::GetContactFrictionLossDens(doublereal& Pfc, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetContactFrictionLossDens(Pfc, pDofMap);
    }

    void HydroSlaveNode::GetContactFrictionLossDens(Gradient<0>& Pfc, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetContactFrictionLossDens(Pfc, pDofMap);
    }

    void HydroSlaveNode::GetClearance(doublereal& h, LocalDofMap*) const
    {
        pMasterNode->GetClearance(h, 0);
    }

    void HydroSlaveNode::GetClearance(Gradient<0>& h, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetClearance(h, pDofMap);
    }

    void HydroSlaveNode::GetClearanceDerTime(doublereal& dh_dt, LocalDofMap*) const
    {
        pMasterNode->GetClearanceDerTime(dh_dt, 0);
    }

    void HydroSlaveNode::GetClearanceDerTime(Gradient<0>& dh_dt, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetClearanceDerTime(dh_dt, pDofMap);
    }

    void HydroSlaveNode::GetRadialDeformation(doublereal& w, doublereal& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap);
    }

    void HydroSlaveNode::GetRadialDeformation(Gradient<0>& w, Gradient<0>& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap);
    }

    void HydroSlaveNode::GetRadialDeformation1(doublereal& w1) const
    {
        pMasterNode->GetRadialDeformation1(w1);
    }

    void HydroSlaveNode::GetRadialDeformation2(doublereal& w2) const
    {
        pMasterNode->GetRadialDeformation2(w2);
    }

    void HydroSlaveNode::GetVelocity(Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetVelocity(U1, U2, pDofMap);
    }

    void HydroSlaveNode::GetVelocity(Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetVelocity(U1, U2, pDofMap);
    }

    void HydroSlaveNode::GetHydraulicVelocity(Vector<doublereal, 2>& U, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetHydraulicVelocity(U, pDofMap);
    }

    void HydroSlaveNode::GetHydraulicVelocity(Vector<Gradient<0>, 2>& U, LocalDofMap* pDofMap) const
    {
        pMasterNode->GetHydraulicVelocity(U, pDofMap);
    }

    const FluidStateBoundaryCond* HydroSlaveNode::pGetMovingPressBoundCond() const
    {
        return pMasterNode->pGetMovingPressBoundCond();
    }

    void HydroSlaveNode::SetMovingPressBoundCond(const FluidStateBoundaryCond* pBoundCond)
    {
        pMasterNode->SetMovingPressBoundCond(pBoundCond);
    }

    index_type HydroSlaveNode::iGetComplianceIndex() const
    {
        return pMasterNode->iGetComplianceIndex();
    }

    HydroMasterNode::HydroMasterNode(integer iNodeNo,
                                     const Vector<doublereal, 2>& x,
                                     HydroMesh* pMesh,
                                     integer iNodeFlags)
        :HydroNode(iNodeNo, x, pMesh, iNodeFlags | MASTER_NODE),
         pMovingBoundCond(nullptr)
    {

    }

    HydroMasterNode::~HydroMasterNode()
    {

    }

    void HydroMasterNode::Update(const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 doublereal dCoef,
                                 enum FunctionCall func,
                                 LocalDofMap* pDofMap)
    {
        pMovingBoundCond = nullptr;
    }

    const FluidStateBoundaryCond* HydroMasterNode::pGetMovingPressBoundCond() const
    {
        return pMovingBoundCond;
    }

    void HydroMasterNode::SetMovingPressBoundCond(const FluidStateBoundaryCond* pBoundCond)
    {
        pMovingBoundCond = pBoundCond;
    }

    HydroUpdatedNode::HydroUpdatedNode(integer iNodeNo,
                                       const Vector<doublereal, 2>& x,
                                       HydroMesh* pMesh,
                                       ContactModel* pContactModel,
                                       std::unique_ptr<FrictionModel>&& pFrictionModel,
                                       integer iNodeFlags)
        :HydroMasterNode(iNodeNo, x, pMesh, iNodeFlags | UPDATED_NODE),
         sum_tau_xy_0(0.),
         sum_tau_yz_0(0.),
         sum_tau_xy_h(0.),
         sum_tau_yz_h(0.),
         iNumStressEval(0),
         pContactModel(pContactModel),
         pFrictionModel(std::move(pFrictionModel)),
         pComplianceModel(pMesh->pGetComplianceModel()),
         iComplianceIndex(-1)
    {
        if (pComplianceModel) {
            iComplianceIndex = pComplianceModel->iGetNumNodes() + 1;
            pComplianceModel->SetNode(iComplianceIndex, this);
        }

        const BearingGeometry* pGeometry = pGetMesh()->pGetGeometry();
        pGeometry->GetPosition3D(GetPosition2D(), v);
        pGeometry->GetTangentCoordSys(GetPosition2D(), Rt);
    }

    HydroUpdatedNode::~HydroUpdatedNode()
    {

    }

    const Vector<doublereal, 3>&
    HydroUpdatedNode::GetPosition3D() const
    {
        return v;
    }

    const Matrix<doublereal, 3, 3>&
    HydroUpdatedNode::GetTangentCoordSys() const
    {
        return Rt;
    }

    void
    HydroUpdatedNode::Update(const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr,
                             doublereal dCoef,
                             enum FunctionCall func,
                             LocalDofMap* pDofMap)
    {
        HydroMasterNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);

        switch (func) {
        case REGULAR_RES:
        case INITIAL_ASS_RES:
        case INITIAL_DER_RES:
            sum_tau_xy_0 = 0.;
            sum_tau_yz_0 = 0.;
            sum_tau_xy_h = 0.;
            sum_tau_yz_h = 0.;
            iNumStressEval = 0;
            oBoundary.Update(this, dCoef, func, pDofMap);
            break;

        case REGULAR_JAC:
        case INITIAL_ASS_JAC:
        case INITIAL_DER_JAC:
            oBoundary_grad.Update(this, dCoef, func, pDofMap);
            break;

        default:
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    void HydroUpdatedNode::AfterPredict(VectorHandler& X,
                                        VectorHandler& XP)
    {
        if (pFrictionModel) {
            pFrictionModel->AfterPredict(X, XP);
        }
    }

    void HydroUpdatedNode::AfterConvergence(const VectorHandler& X,
                                            const VectorHandler& XP)
    {
        if (pFrictionModel) {
            pFrictionModel->AfterConvergence(X, XP);
        }
    }

    void HydroUpdatedNode::GetStress(doublereal& tau_xy_0, doublereal& tau_yz_0, doublereal& tau_xy_h, doublereal& tau_yz_h) const
    {
        // avoid division by zero
        const integer n = iNumStressEval > 0 ? iNumStressEval : 1;

        tau_xy_0 = sum_tau_xy_0 / n;
        tau_yz_0 = sum_tau_yz_0 / n;
        tau_xy_h = sum_tau_xy_h / n;
        tau_yz_h = sum_tau_yz_h / n;
    }

    void HydroUpdatedNode::SetStress(doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h)
    {
        sum_tau_xy_0 += tau_xy_0;
        sum_tau_yz_0 += tau_yz_0;
        sum_tau_xy_h += tau_xy_h;
        sum_tau_yz_h += tau_yz_h;

        ++iNumStressEval;
    }

    bool HydroUpdatedNode::GetContactPressure(doublereal& pasp, LocalDofMap* pDofMap) const
    {
        return oBoundary.GetContactPressure(pasp, pDofMap);
    }

    bool HydroUpdatedNode::GetContactPressure(Gradient<0>& pasp, LocalDofMap* pDofMap) const
    {
        return oBoundary_grad.GetContactPressure(pasp, pDofMap);
    }

    void HydroUpdatedNode::GetContactStress(Vector<doublereal, 2>& tauc_0, LocalDofMap* pDofMap) const
    {
        oBoundary.GetContactStress(tauc_0, pDofMap);
    }

    void HydroUpdatedNode::GetContactStress(Vector<Gradient<0>, 2>& tauc_0, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetContactStress(tauc_0, pDofMap);
    }

    void HydroUpdatedNode::GetContactFrictionLossDens(doublereal& Pfc, LocalDofMap* pDofMap) const
    {
        oBoundary.GetContactFrictionLossDens(Pfc, pDofMap);
    }

    void HydroUpdatedNode::GetContactFrictionLossDens(Gradient<0>& Pfc, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetContactFrictionLossDens(Pfc, pDofMap);
    }

    void HydroUpdatedNode::GetClearance(doublereal& h, LocalDofMap*) const
    {
        oBoundary.GetClearance(h, 0);
    }

    void HydroUpdatedNode::GetClearance(Gradient<0>& h, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetClearance(h, pDofMap);
    }

    void HydroUpdatedNode::GetClearanceDerTime(doublereal& dh_dt, LocalDofMap*) const
    {
        oBoundary.GetClearanceDerTime(dh_dt, 0);
    }

    void HydroUpdatedNode::GetClearanceDerTime(Gradient<0>& dh_dt, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetClearanceDerTime(dh_dt, pDofMap);
    }

    void HydroUpdatedNode::GetRadialDeformation(doublereal& w, doublereal& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const
    {
        if (pComplianceModel) {
            pComplianceModel->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap, this);
        } else {
            w = 0.;
            dw_dt = 0.;
        }
    }

    void HydroUpdatedNode::GetRadialDeformation(Gradient<0>& w, Gradient<0>& dw_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap) const
    {
        if (pComplianceModel) {
            pComplianceModel->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap, this);
        } else {
            w.SetValue(0.);
            dw_dt.SetValue(0.);
        }
    }

    void HydroUpdatedNode::GetRadialDeformation1(doublereal& w1) const
    {
        if (pComplianceModel) {
            pComplianceModel->GetRadialDeformation1(w1, this);
        } else {
            w1 = 0.;
        }
    }

    void HydroUpdatedNode::GetRadialDeformation2(doublereal& w2) const
    {
        if (pComplianceModel) {
            pComplianceModel->GetRadialDeformation2(w2, this);
        } else {
            w2 = 0.;
        }
    }

    void HydroUpdatedNode::GetVelocity(Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, LocalDofMap* pDofMap) const
    {
        oBoundary.GetVelocity(U1, U2, pDofMap);
    }

    void HydroUpdatedNode::GetVelocity(Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetVelocity(U1, U2, pDofMap);
    }

    void HydroUpdatedNode::GetHydraulicVelocity(Vector<doublereal, 2>& U, LocalDofMap* pDofMap) const
    {
        oBoundary.GetHydraulicVelocity(U, pDofMap);
    }

    void HydroUpdatedNode::GetHydraulicVelocity(Vector<Gradient<0>, 2>& U, LocalDofMap* pDofMap) const
    {
        oBoundary_grad.GetHydraulicVelocity(U, pDofMap);
    }

    index_type HydroUpdatedNode::iGetComplianceIndex() const
    {
        return iComplianceIndex;
    }

    HydroActiveNode::HydroActiveNode(integer iNodeNo,
                                     const Vector<doublereal, 2>& x,
                                     HydroMesh* pParent,
                                     ContactModel* pContactModel,
                                     std::unique_ptr<FrictionModel>&& pFrictionModel)
        :HydroIncompressibleNode(iNodeNo,
                                 x,
                                 pParent,
                                 pContactModel,
                                 std::move(pFrictionModel),
                                 ACTIVE_NODE),
         p(0.),
         dp_dt(0.),
         s(pParent->pGetParent()->dGetScale(HydroRootElement::SCALE_PRESSURE_DOF)),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {

    }

    HydroActiveNode::~HydroActiveNode()
    {

    }

    integer HydroActiveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return iGetFirstDofIndex(eFunc);
    }

    integer HydroActiveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(iGetOffsetIndex(eFunc) != HydroDofOwner::UNKNOWN_OFFSET);

        return pGetMesh()->pGetParent()->iGetFirstIndex() + iGetOffsetIndex(eFunc);
    }

    integer HydroActiveNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        integer iNumCols = HydroIncompressibleNode::iGetNumColsWorkSpace(eFunc);

        if (eFunc & grad::REGULAR_FLAG) {
            iNumCols += iGetNumDof();
        } else {
            iNumCols += iGetInitialNumDof();
        }

        return iNumCols;
    }

    void
    HydroActiveNode::Update(const VectorHandler& XCurr,
                            const VectorHandler& XPrimeCurr,
                            doublereal dCoef,
                            enum FunctionCall func,
                            LocalDofMap* pDofMap)
    {
        HYDRO_ASSERT(iGetOffsetIndex(func) > 0);

        const integer iIndex = iGetFirstDofIndex(func);

        HYDRO_ASSERT(iIndex > 0);
        HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

        p = s * XCurr(iIndex);

        switch (func) {
        case REGULAR_RES:
        case REGULAR_JAC:
            HYDRO_ASSERT(eCurrFunc == grad::REGULAR_FLAG);

            dp_dt = s * XPrimeCurr(iIndex);
            break;

        case INITIAL_ASS_RES:
        case INITIAL_ASS_JAC:
            HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);
            break;

        default:
            break;
        }

        HydroIncompressibleNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
    }

    void HydroActiveNode::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);

        eCurrFunc = grad::REGULAR_FLAG;

        const integer iIndex = iGetFirstDofIndex(eCurrFunc);

        HYDRO_ASSERT(iIndex > 0);
        HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

        XCurr(iIndex) = p / s;
        XPrimeCurr(iIndex) = dp_dt / s;
    }

    void HydroActiveNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const {
        p = this->p;
    }

    void HydroActiveNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const {
        p.SetValuePreserve(this->p);
        p.DerivativeResizeReset(pDofMap,
                                iGetFirstDofIndex(eCurrFunc),
                                MapVectorBase::GLOBAL,
                                -s);
    }

    void HydroActiveNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap*, doublereal) const {
        dp_dt = this->dp_dt;
    }

    void HydroActiveNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const {
        HYDRO_ASSERT(eCurrFunc == grad::REGULAR_FLAG || eCurrFunc == grad::INITIAL_ASS_FLAG);

        if (eCurrFunc == grad::REGULAR_FLAG) {
            dp_dt.SetValuePreserve(this->dp_dt);
            // We assume that db0Algebraic == db0Differential
            // In case of multistep and hope methods this is correct
            // only if algebraic and differential spectral radii are the same!
            dp_dt.DerivativeResizeReset(pDofMap,
                                        iGetFirstDofIndex(eCurrFunc),
                                        MapVectorBase::GLOBAL,
                                        -dCoef * s);
        } else {
            dp_dt.SetValue(this->dp_dt);
        }
    }

    unsigned int HydroActiveNode::iGetNumDof(void) const
    {
        return 1u;
    }

    unsigned int HydroActiveNode::iGetInitialNumDof(void) const
    {
        return 1u;
    }

    DofOrder::Order HydroActiveNode::GetDofType(unsigned int i) const
    {
        HYDRO_ASSERT(i == 0);

        return DofOrder::ALGEBRAIC;
    }

    DofOrder::Order HydroActiveNode::GetEqType(unsigned int i) const
    {
        HYDRO_ASSERT(i == 0);

        return DofOrder::ALGEBRAIC;
    }

    std::ostream&
    HydroActiveNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        integer iIndex = iGetFirstDofIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

        out << prefix << iIndex << ": hydrodynamic pressure of node number " << iGetNodeNumber() + 1 << std::endl;

        return out;
    }

    std::ostream&
    HydroActiveNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        integer iIndex = iGetFirstEquationIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

        out << prefix << iIndex << ": incompressible Reynolds equation for node number " << iGetNodeNumber() + 1 << std::endl;

        return out;
    }

    HydroPassiveNode::HydroPassiveNode(integer iNodeNo,
                                       const Vector<doublereal, 2>& x,
                                       HydroMesh* pParent,
                                       ContactModel* pContactModel,
                                       std::unique_ptr<FrictionModel>&& pFrictionModel,
                                       const FluidStateBoundaryCond* pBoundaryCond)
        :HydroIncompressibleNode(iNodeNo,
                                 x,
                                 pParent,
                                 pContactModel,
                                 std::move(pFrictionModel),
                                 PASSIVE_NODE),
         pBoundaryCond(pBoundaryCond)
    {
        HYDRO_ASSERT(pBoundaryCond->bIncludeNode(GetNodePhysics()));
    }

    HydroPassiveNode::~HydroPassiveNode()
    {

    }

    integer HydroPassiveNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    integer HydroPassiveNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    void HydroPassiveNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const
    {
        p = pGetBoundCond()->dGetPressure();
    }

    void HydroPassiveNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        p.SetValue(pGetBoundCond()->dGetPressure());
    }

    void HydroPassiveNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap*, doublereal) const
    {
        dp_dt = pGetBoundCond()->dGetPressureDerTime();
    }

    void HydroPassiveNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        dp_dt.SetValue(pGetBoundCond()->dGetPressureDerTime());
    }

    const FluidStateBoundaryCond*
    HydroPassiveNode::pGetBoundCond() const
    {
        const FluidStateBoundaryCond* pBound = pGetMovingPressBoundCond();

        if (!pBound) {
            pBound = pBoundaryCond;
        }

        return pBound;
    }

    HydroIncompressibleNode::HydroIncompressibleNode(integer iNodeNo,
                                                     const Vector<doublereal, 2>& x,
                                                     HydroMesh* pParent,
                                                     ContactModel* pContactModel,
                                                     std::unique_ptr<FrictionModel>&& pFrictionModel,
                                                     integer iNodeFlags)
        :HydroUpdatedNode(iNodeNo,
                          x,
                          pParent,
                          pContactModel,
                          std::move(pFrictionModel),
                          iNodeFlags | INCOMPRESSIBLE_NODE)
    {

    }

    HydroIncompressibleNode::~HydroIncompressibleNode()
    {

    }

    void HydroIncompressibleNode::GetDensity(doublereal& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        rho = oState.rho;
    }

    void HydroIncompressibleNode::GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        rho.Copy(oState_grad.rho, pDofMap);
    }

    void HydroIncompressibleNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        drho_dt = oState.drho_dt;
    }

    void HydroIncompressibleNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        drho_dt.Copy(oState_grad.drho_dt, pDofMap);
    }

    void
    HydroIncompressibleNode::Update(const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr,
                                    doublereal dCoef,
                                    enum FunctionCall func,
                                    LocalDofMap* pDofMap)
    {
        HydroUpdatedNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);

        switch (func) {
        case INITIAL_ASS_RES:
        case REGULAR_RES:
            UpdateState(oState, pDofMap, dCoef);
            break;

        case INITIAL_ASS_JAC:
        case REGULAR_JAC:
            UpdateState(oState_grad, pDofMap, dCoef);
            break;

        default:
            NO_OP;
        }
    }

    template <typename G>
    inline void
    HydroIncompressibleNode::UpdateState(FluidState<G>& oState, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        G p, T, dT_dt, drho_dp, drho_dT;

        GetPressure(p, pDofMap, dCoef);
        GetTemperature(T, pDofMap, dCoef);
        GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
        pGetFluid()->GetDensity(p, T, oState.rho, &drho_dp, &drho_dT);
        oState.drho_dt = drho_dT * dT_dt;

        HYDRO_ASSERT(drho_dp == 0.);
    }

    HydroCoupledNode::HydroCoupledNode(integer iNodeNo,
                                       const Vector<doublereal, 2>& x,
                                       HydroMesh* pParent,
                                       ContactModel* pContactModel,
                                       std::unique_ptr<FrictionModel>&& pFrictionModel,
                                       const PressureNode* pNode)
        :HydroIncompressibleNode(iNodeNo,
                                 x,
                                 pParent,
                                 pContactModel,
                                 std::move(pFrictionModel),
                                 COUPLED_NODE),
         pExtNode(pNode),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {

    }

    HydroCoupledNode::~HydroCoupledNode()
    {

    }

    integer HydroCoupledNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return pExtNode->iGetFirstRowIndex() + 1;
    }

    integer HydroCoupledNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        return pExtNode->iGetFirstColIndex() + 1;
    }

    integer HydroCoupledNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return 1 + HydroIncompressibleNode::iGetNumColsWorkSpace(eFunc);
    }

    void HydroCoupledNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const
    {
        p = pExtNode->dGetX();
    }

    void HydroCoupledNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (eCurrFunc & grad::REGULAR_FLAG) {
            p.SetValuePreserve(pExtNode->dGetX());
            p.DerivativeResizeReset(pDofMap,
                                    iGetFirstDofIndex(grad::UNKNOWN_FUNC),
                                    MapVectorBase::GLOBAL,
                                    -1.);   // ScalarNodes are inactive during initial assembly
        } else {
            p.SetValue(pExtNode->dGetX());
        }
    }

    void HydroCoupledNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap*, doublereal) const
    {
        dp_dt = pExtNode->dGetXPrime();
    }

    void HydroCoupledNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (eCurrFunc & grad::REGULAR_FLAG) {
            // We assume that db0Algebraic == db0Differential
            // In case of the multistep and hope methods this is true only if algebraic and differential spectral radii are the same!

            dp_dt.SetValuePreserve(pExtNode->dGetXPrime());
            dp_dt.DerivativeResizeReset(pDofMap,
                                        iGetFirstDofIndex(grad::UNKNOWN_FUNC),
                                        MapVectorBase::GLOBAL,
                                        -dCoef); // ScalarNodes are inactive during initial assembly
        } else {
            dp_dt.SetValue(pExtNode->dGetXPrime());
        }
    }

    void
    HydroCoupledNode::Update(const VectorHandler& XCurr,
                             const VectorHandler& XPrimeCurr,
                             doublereal dCoef,
                             enum FunctionCall func,
                             LocalDofMap* pDofMap)
    {
        if (func & grad::REGULAR_FLAG) {
            eCurrFunc = grad::REGULAR_FLAG;
        }

        HydroIncompressibleNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
    }

    HydroCompressibleNode::HydroCompressibleNode(integer iNodeNo,
                                                 const Vector<doublereal, 2>& x,
                                                 HydroMesh* pParent,
                                                 ContactModel* pContactModel,
                                                 std::unique_ptr<FrictionModel>&& pFrictionModel,
                                                 integer iNodeFlags)
        :HydroUpdatedNode(iNodeNo,
                          x,
                          pParent,
                          pContactModel,
                          std::move(pFrictionModel),
                          iNodeFlags | COMPRESSIBLE_NODE)
    {

    }

    HydroCompressibleNode::~HydroCompressibleNode()
    {

    }

    HydroActiveComprNode::HydroActiveComprNode(integer iNodeNo,
                                               const Vector<doublereal, 2>& x,
                                               HydroMesh* pParent,
                                               ContactModel* pContactModel,
                                               std::unique_ptr<FrictionModel>&& pFrictionModel)
        :HydroCompressibleNode(iNodeNo,
                               x,
                               pParent,
                               pContactModel,
                               std::move(pFrictionModel),
                               ACTIVE_NODE),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {
        std::array<HydroRootElement::ScaleType, iNumDofMax> rgScale = {
            HydroRootElement::SCALE_PRESSURE_DOF,
            HydroRootElement::SCALE_THETA_DOF
        };

        oCurrState.t = pGetMesh()->pGetParent()->dGetTime();
        oCurrState.eCavitationState = HydroFluid::FULL_FILM_REGION;

        for (index_type i = 0; i < iNumDofMax; ++i) {
            oCurrState.Theta[i] = pGetFluid()->GetTheta0(i);
            oCurrState.dTheta_dt[i] = 0.;
            s[i] = pParent->pGetParent()->dGetScale(rgScale[i]);
        }

        oRefState = oIncState  = oCurrState;

        UpdateState(oState);
    }

    HydroActiveComprNode::~HydroActiveComprNode()
    {

    }

    HydroFluid::CavitationState HydroActiveComprNode::GetCavitationState() const
    {
        return oCurrState.eCavitationState;
    }

    void HydroActiveComprNode::GetTheta(std::array<doublereal, iNumDofMax>& Theta, LocalDofMap*, doublereal) const
    {
        for (index_type i = 0; i < iNumDofMax; ++i) {
            Theta[i] = oCurrState.Theta[i];
        }
    }

    void HydroActiveComprNode::GetTheta(std::array<Gradient<0>, iNumDofMax>& Theta, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG || eCurrFunc == grad::REGULAR_FLAG);

        HYDRO_ASSERT(dCoef == 1. || eCurrFunc != grad::INITIAL_ASS_FLAG);

        const index_type iNumDofInit = iGetInitialNumDof();

        for (index_type i = 0; i < iNumDofMax; ++i) {
            Theta[i].SetValuePreserve(oCurrState.Theta[i]);

            doublereal dDeriv = 0;

            if (eCurrFunc == grad::REGULAR_FLAG) {
                dDeriv = -dCoef * s[i];
            } else if (i < iNumDofInit) {
                dDeriv = -s[i];
            }

            Theta[i].DerivativeResizeReset(pDofMap,
                                           iGetFirstDofIndex(eCurrFunc) + i,
                                           MapVectorBase::GLOBAL,
                                           dDeriv);
        }
    }

    void HydroActiveComprNode::GetThetaDerTime(std::array<doublereal, iNumDofMax>& dTheta_dt, LocalDofMap*, doublereal) const
    {
        for (index_type i = 0; i < iNumDofMax; ++i) {
            dTheta_dt[i] = oCurrState.dTheta_dt[i];
        }
    }

    void HydroActiveComprNode::GetThetaDerTime(std::array<Gradient<0>, iNumDofMax>& dTheta_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG || eCurrFunc == grad::REGULAR_FLAG);

        for (index_type i = 0; i < iNumDofMax; ++i) {
            if (eCurrFunc & grad::REGULAR_FLAG) {
                dTheta_dt[i].SetValuePreserve(oCurrState.dTheta_dt[i]);
                dTheta_dt[i].DerivativeResizeReset(pDofMap,
                                                   iGetFirstDofIndex(eCurrFunc) + i,
                                                   MapVectorBase::GLOBAL,
                                                   -s[i]);
            } else {
                dTheta_dt[i].SetValue(oCurrState.dTheta_dt[i]);
            }
        }
    }

    template <typename G>
    inline void HydroActiveComprNode::UpdateState(FluidState<G>& oState, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        GetTheta(oState.Theta, pDofMap, dCoef);
        GetThetaDerTime(oState.dTheta_dt, pDofMap, dCoef);
        GetTemperature(oState.T, pDofMap, dCoef);
        GetTemperatureDerTime(oState.dT_dt, pDofMap, dCoef);

        pGetFluid()->ThetaToPhysical(oState.Theta,
                                     oState.dTheta_dt,
                                     oState.T,
                                     oState.dT_dt,
                                     oState.p,
                                     oState.dp_dt,
                                     oState.rho,
                                     oState.drho_dt);

        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.Theta[0])));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.Theta[1])));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.dTheta_dt[0])));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.dTheta_dt[1])));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.T)));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(oState.dT_dt)));
    }

    integer HydroActiveComprNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return iGetFirstDofIndex(eFunc);
    }

    integer HydroActiveComprNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(iGetOffsetIndex(eFunc) != HydroDofOwner::UNKNOWN_OFFSET);

        return pGetMesh()->pGetParent()->iGetFirstIndex() + iGetOffsetIndex(eFunc);
    }

    integer HydroActiveComprNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        integer iNumCols = HydroCompressibleNode::iGetNumColsWorkSpace(eFunc);

        if (eFunc & grad::REGULAR_FLAG) {
            iNumCols += iGetNumDof();
        } else {
            iNumCols += iGetInitialNumDof();
        }

        return iNumCols;
    }

    void HydroActiveComprNode::UpdateTheta(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr)
    {
        oCurrState.t = pGetMesh()->pGetParent()->dGetTime();

        HYDRO_ASSERT(iGetOffsetIndex(eCurrFunc) > 0);

        for (index_type i = 0; i < iNumDofMax; ++i) {
            const integer iIndex = iGetFirstDofIndex(eCurrFunc) + i;

            HYDRO_ASSERT(iIndex > 0);
            HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

            oCurrState.Theta[i] = s[i] * XCurr(iIndex);
            oCurrState.dTheta_dt[i] = s[i] * XPrimeCurr(iIndex);
        }
    }

    void HydroActiveComprNode::UpdateCavitationState()
    {
        if (oCurrState.eCavitationState == HydroFluid::FULL_FILM_REGION) {
            if (oCurrState.Theta[0] < 0.) {
                oCurrState.eCavitationState = HydroFluid::CAVITATION_REGION;
            }
        } else {
            if (oCurrState.Theta[1] > 1.) {
                oCurrState.eCavitationState = HydroFluid::FULL_FILM_REGION;
            }
        }
    }

    void HydroActiveComprNode::ResolveCavitationState(VectorHandler& X, VectorHandler& XP)
    {
        if (oCurrState.eCavitationState == HydroFluid::CAVITATION_REGION) {
            oCurrState.Theta[0] = 0.;
            oCurrState.dTheta_dt[0] = 0.;
        } else {
            oCurrState.Theta[1] = 1.;
            oCurrState.dTheta_dt[1] = 0.;
        }

        for (index_type i = 0; i < iNumDofMax; ++i) {
            const integer iIndex = iGetFirstDofIndex(eCurrFunc) + i;

            HYDRO_ASSERT(iIndex > 0);
            HYDRO_ASSERT(iIndex <= X.iGetSize());

            X.PutCoef(iIndex, oCurrState.Theta[i] / s[i]);
            XP.PutCoef(iIndex, oCurrState.dTheta_dt[i] / s[i]);
        }
    }

    void
    HydroActiveComprNode::Update(const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 doublereal dCoef,
                                 enum FunctionCall func,
                                 LocalDofMap* pDofMap)
    {
        oDofMap.Reset(func);

        if (func & grad::INITIAL_ASS_FLAG) {
            oCurrState.Theta[0] = s[0] * XCurr(iGetFirstDofIndex(func));
            HYDRO_ASSERT(std::isfinite(oCurrState.Theta[0]));
        }

        HYDRO_ASSERT(eCurrFunc == (func & grad::STATE_MASK));

        switch (func) {
        case INITIAL_ASS_RES:
        case REGULAR_RES:
            UpdateState(oState, &oDofMap, dCoef);
            break;

        case INITIAL_ASS_JAC:
        case REGULAR_JAC:
            UpdateState(oState_grad, &oDofMap, dCoef);
            break;

        default:
            NO_OP;
        }

        HydroCompressibleNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
    }

    void HydroActiveComprNode::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
        UpdateTheta(X, XP);
        UpdateCavitationState();
        ResolveCavitationState(X, XP);
    }

    void HydroActiveComprNode::DofUpdate(VectorHandler& X, VectorHandler& XP)
    {
        const State oPrevState = oCurrState;

        UpdateTheta(X, XP);

        HydroRootElement* const pParent = pGetMesh()->pGetParent();
        const integer iIter = pParent->GetNonlinearSolverHint(NonlinearSolver::LINESEARCH_ITERATION_CURR);

        if (iIter == 0) {
            HYDRO_ASSERT(pParent->GetNonlinearSolverHint(NonlinearSolver::LINESEARCH_LAMBDA_CURR) == 1.);

            oIncState = oCurrState;
            oRefState = oPrevState;

            static const doublereal Thetax[iNumDofMax] = {0., 1.};

            if (oCurrState.eCavitationState != oPrevState.eCavitationState) {
                const index_type i = oPrevState.eCavitationState == HydroFluid::FULL_FILM_REGION ? 0 : 1;
                doublereal dLambdaLimit = (Thetax[i] - oPrevState.Theta[i]) / (oCurrState.Theta[i] - oPrevState.Theta[i]);

                if (dLambdaLimit > 0 && dLambdaLimit < 1) {
                    pedantic_cout("hydrodynamic plain bearing2("
                                  << pGetMesh()->pGetParent()->GetLabel()
                                  << "): fluid state Theta[" << i << "] of node " << iGetNodeNumber() + 1
                                  << " changed from " << oPrevState.Theta[i]
                                  << " to " << oCurrState.Theta[i]
                                  << " at t=" << oCurrState.t
                                  << " lambda=" << dLambdaLimit << std::endl);
                    dLambdaLimit *= (1. + sqrt(std::numeric_limits<doublereal>::epsilon())); // Make sure that state will not jump back
                    pParent->SetNonlinearSolverHint(NonlinearSolver::LINESEARCH_LAMBDA_MAX, dLambdaLimit);
                }
            }
        } else {
            const doublereal dLambda = pParent->GetNonlinearSolverHint(NonlinearSolver::LINESEARCH_LAMBDA_CURR);

            for (index_type i = 0; i < iNumDofMax; ++i) {
                oCurrState.Theta[i] = oRefState.Theta[i] + (oIncState.Theta[i] - oRefState.Theta[i]) * dLambda;
                oCurrState.dTheta_dt[i] = oRefState.dTheta_dt[i] + (oIncState.dTheta_dt[i] - oRefState.dTheta_dt[i]) * dLambda;
            }
        }

        UpdateCavitationState();
        ResolveCavitationState(X, XP);
    }

    void HydroActiveComprNode::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
    {
    }

    void HydroActiveComprNode::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);
        HYDRO_ASSERT(oCurrState.eCavitationState == HydroFluid::FULL_FILM_REGION);

        eCurrFunc = grad::REGULAR_FLAG;
        UpdateCavitationState();
        ResolveCavitationState(XCurr, XPrimeCurr);
    }

    void HydroActiveComprNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const
    {
        p = oState.p;
    }

    void HydroActiveComprNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        p.Copy(oState_grad.p, pDofMap);
    }

    void HydroActiveComprNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap*, doublereal) const
    {
        dp_dt = oState.dp_dt;
    }

    void HydroActiveComprNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        dp_dt.Copy(oState_grad.dp_dt, pDofMap);
    }

    void HydroActiveComprNode::GetDensity(doublereal& rho, LocalDofMap*, doublereal) const
    {
        rho = oState.rho;
    }

    void HydroActiveComprNode::GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        rho.Copy(oState_grad.rho, pDofMap);
    }

    void HydroActiveComprNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap*, doublereal) const
    {
        drho_dt = oState.drho_dt;
    }

    void HydroActiveComprNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        drho_dt.Copy(oState_grad.drho_dt, pDofMap);
    }

    unsigned int HydroActiveComprNode::iGetNumDof(void) const
    {
        return iNumDofMax;
    }

    unsigned int HydroActiveComprNode::iGetInitialNumDof(void) const
    {
        static_assert(iNumDofMax >= 1, "not enough degrees of freedome");

        return 1u;
    }

    DofOrder::Order HydroActiveComprNode::GetDofType(unsigned int i) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(i < iNumDofMax);

        switch (i) {
        case 0:
        case 1:
            // The thermal model requires dp/dt, so pressure should be declared as differential variable
            return DofOrder::DIFFERENTIAL;

        default:
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    DofOrder::Order HydroActiveComprNode::GetEqType(unsigned int i) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(i < iNumDofMax);

        switch (i) {
        case 0:
            return DofOrder::DIFFERENTIAL;

        case 1:
            return DofOrder::ALGEBRAIC;

        default:
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    std::ostream&
    HydroActiveComprNode::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        const integer iIndex = iGetFirstDofIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

        const index_type iNumDof = bInitial ? iGetInitialNumDof() : iGetNumDof();

        for (index_type i = 0; i < iNumDof; ++i) {
            out << prefix << iIndex + i << ": Theta" << iGetNodeNumber() + 1 << "[" << i << "]" << std::endl;
        }

        return out;
    }

    std::ostream&
    HydroActiveComprNode::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        const integer iIndex = iGetFirstDofIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

        static const char szEqName[][31] = {
            "compressible Reynolds equation",
            "equation of fluid state"
        };

        const index_type iNumDof = bInitial ? iGetInitialNumDof() : iGetNumDof();

        for (index_type i = 0; i < iNumDof; ++i) {
            out << prefix << iIndex + i << ": " << szEqName[i] << " " << iGetNodeNumber() + 1 << std::endl;
        }

        return out;
    }

    HydroPassiveComprNode::HydroPassiveComprNode(integer iNodeNo,
                                                 const Vector<doublereal, 2>& x,
                                                 HydroMesh* pParent,
                                                 ContactModel* pContactModel,
                                                 std::unique_ptr<FrictionModel>&& pFrictionModel,
                                                 const FluidStateBoundaryCond* pBoundaryCond)
        :HydroCompressibleNode(iNodeNo,
                               x,
                               pParent,
                               pContactModel,
                               std::move(pFrictionModel),
                               PASSIVE_NODE),
         pBoundaryCond(pBoundaryCond)
    {
        HYDRO_ASSERT(pBoundaryCond->bIncludeNode(GetNodePhysics()));
    }

    HydroPassiveComprNode::~HydroPassiveComprNode()
    {

    }

    integer HydroPassiveComprNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    integer HydroPassiveComprNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(0);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    HydroFluid::CavitationState HydroPassiveComprNode::GetCavitationState() const
    {
        // Should not be called anyway
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }

    void HydroPassiveComprNode::GetPressure(doublereal& p, LocalDofMap*, doublereal) const
    {
        const FluidStateBoundaryCond* const pBound = pGetFluidBoundCond();
        const doublereal h = dGetClearance(pBound);
        p = pBound->dGetPressure(h);
    }

    void HydroPassiveComprNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        doublereal ps;

        GetPressure(ps);
        p.SetValue(ps);
    }

    void HydroPassiveComprNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap*, doublereal) const
    {
        const FluidStateBoundaryCond* const pBound = pGetFluidBoundCond();
        doublereal h, dh_dt;

        h = dGetClearance(pBound, &dh_dt);
        dp_dt = pBound->dGetPressureDerTime(h, dh_dt);
    }

    void HydroPassiveComprNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        doublereal dps_dt;

        GetPressureDerTime(dps_dt);
        dp_dt.SetValue(dps_dt);
    }

    void HydroPassiveComprNode::GetDensity(doublereal& rho, LocalDofMap*, doublereal) const
    {
        const FluidStateBoundaryCond* const pBound = pGetFluidBoundCond();
        const doublereal h = dGetClearance(pBound);
        rho = pBound->dGetDensity(h);
    }

    void HydroPassiveComprNode::GetDensity(Gradient<0>& rho, LocalDofMap*, doublereal) const
    {
        doublereal rhos;

        GetDensity(rhos);
        rho.SetValue(rhos);
    }

    void HydroPassiveComprNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap*, doublereal) const
    {
        const FluidStateBoundaryCond* const pBound = pGetFluidBoundCond();
        doublereal dh_dt;
        const doublereal h = dGetClearance(pBound, &dh_dt);

        drho_dt = pBound->dGetDensityDerTime(h, dh_dt);
    }

    void HydroPassiveComprNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap*, doublereal) const
    {
        doublereal drhos_dt;

        GetDensityDerTime(drhos_dt);
        drho_dt.SetValue(drhos_dt);
    }

    const FluidStateBoundaryCond* HydroPassiveComprNode::pGetFluidBoundCond() const
    {
        const FluidStateBoundaryCond* pBound = pGetMovingPressBoundCond();

        if (pBound == nullptr) {
            pBound = pBoundaryCond;
        }

        return pBound;
    }

    HydroComprOutletNode::HydroComprOutletNode(integer iNodeNo,
                                               const Vector<doublereal, 2>& x,
                                               HydroMesh* pParent,
                                               ContactModel* pContactModel,
                                               std::unique_ptr<FrictionModel>&& pFrictionModel,
                                               const FluidStateBoundaryCond* pBoundaryCond,
                                               const HydroMasterNode* pMasterNode)
    :HydroPassiveComprNode(iNodeNo, x, pParent, pContactModel, std::move(pFrictionModel), pBoundaryCond),
     pMasterNode(pMasterNode)
    {
    }

    HydroComprOutletNode::~HydroComprOutletNode()
    {
    }

    void HydroComprOutletNode::GetDensity(doublereal& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensity(rho, pDofMap, dCoef);
    }

    void HydroComprOutletNode::GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensity(rho, pDofMap, dCoef);
    }

    void HydroComprOutletNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensityDerTime(drho_dt, pDofMap, dCoef);
    }

    void HydroComprOutletNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pMasterNode->GetDensityDerTime(drho_dt, pDofMap, dCoef);
    }

    HydroCoupledComprNode::HydroCoupledComprNode(integer iNodeNo,
                                                 const Vector<doublereal, 2>& x,
                                                 HydroMesh* pParent,
                                                 ContactModel* pContactModel,
                                                 std::unique_ptr<FrictionModel>&& pFrictionModel,
                                                 PressureNode* pExtNode)
        :HydroCompressibleNode(iNodeNo, x, pParent, pContactModel, std::move(pFrictionModel), COUPLED_NODE),
         pExtNode(pExtNode),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {
        if (pExtNode->dGetX() < pGetFluid()->dGetRefPressure()) {
            pExtNode->SetX(pGetFluid()->dGetRefPressure());
        }

        pext = pExtNode->dGetX();
        dpext_dt = 0.;
    }

    HydroCoupledComprNode::~HydroCoupledComprNode()
    {
    }

    integer HydroCoupledComprNode::iGetFirstEquationIndex(grad::FunctionCall eFunc) const
    {
        return pExtNode->iGetFirstRowIndex() + 1;
    }

    integer HydroCoupledComprNode::iGetFirstDofIndex(grad::FunctionCall eFunc) const
    {
        return pExtNode->iGetFirstColIndex() + 1;
    }

    integer HydroCoupledComprNode::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return eFunc & grad::REGULAR_FLAG ? 1 : 0;
    }

    void HydroCoupledComprNode::GetPressure(doublereal& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        p = oState.p;
    }

    void HydroCoupledComprNode::GetPressure(Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        p.Copy(oState_grad.p, pDofMap);
    }

    void HydroCoupledComprNode::GetPressureDerTime(doublereal& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        dp_dt = oState.dp_dt;
    }

    void HydroCoupledComprNode::GetPressureDerTime(Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        dp_dt.Copy(oState_grad.dp_dt, pDofMap);
    }

    void HydroCoupledComprNode::GetDensity(doublereal& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        rho = oState.rho;
    }

    void HydroCoupledComprNode::GetDensity(Gradient<0>& rho, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        rho.Copy(oState_grad.rho, pDofMap);
    }

    HydroFluid::CavitationState HydroCoupledComprNode::GetCavitationState() const
    {
        return oState.eCavitationState;
    }

    void HydroCoupledComprNode::GetDensityDerTime(doublereal& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        drho_dt = oState.drho_dt;
    }

    void HydroCoupledComprNode::GetDensityDerTime(Gradient<0>& drho_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        drho_dt.Copy(oState_grad.drho_dt, pDofMap);
    }

    void
    HydroCoupledComprNode::Update(const VectorHandler& XCurr,
                                  const VectorHandler& XPrimeCurr,
                                  doublereal dCoef,
                                  enum FunctionCall func,
                                  LocalDofMap* pDofMap)
    {
        oDofMap.Reset(func);

        if (func & grad::REGULAR_FLAG) {
            eCurrFunc = grad::REGULAR_FLAG;
            const integer iIndex = iGetFirstDofIndex(func);

            HYDRO_ASSERT(iIndex > 0);
            HYDRO_ASSERT(iIndex <= XCurr.iGetSize());

            pext = XCurr(iIndex);
            dpext_dt = XPrimeCurr(iIndex);
            eCurrFunc = grad::REGULAR_FLAG;
        }

        switch (func) {
        case REGULAR_RES:
        case INITIAL_ASS_RES:
            UpdateState(oState, &oDofMap, dCoef, func);
            break;

        case REGULAR_JAC:
        case INITIAL_ASS_JAC:
            UpdateState(oState_grad, &oDofMap, dCoef, func);
            break;

        default:
            break;
        }

        HydroCompressibleNode::Update(XCurr, XPrimeCurr, dCoef, func, pDofMap);
    }

    void HydroCoupledComprNode::GetExtPressure(doublereal& p, doublereal& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        p = pext;
        dp_dt = dpext_dt;
    }

    void HydroCoupledComprNode::GetExtPressure(Gradient<0>& p, Gradient<0>& dp_dt, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        if (eCurrFunc == grad::REGULAR_FLAG) {
            const index_type iIndex = iGetFirstDofIndex(eCurrFunc);

            p.SetValuePreserve(pext);
            p.DerivativeResizeReset(pDofMap, iIndex, MapVectorBase::GLOBAL, -1.);

            // We assume that db0Algebraic == db0Differential
            // In case of the multistep and hope methods this is true only if algebraic and differential spectral radii are the same!

            dp_dt.SetValuePreserve(dpext_dt);
            dp_dt.DerivativeResizeReset(pDofMap, iIndex, MapVectorBase::GLOBAL, -dCoef);
        } else {
            HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);
            HYDRO_ASSERT(dpext_dt == 0.);

            p.SetValue(pext); // ScalarNodes are inactive during initial assembly
            dp_dt.SetValue(dpext_dt);
        }
    }

    template <typename G>
    void HydroCoupledComprNode::UpdateState(FluidState<G>& oStateCurr, LocalDofMap* pDofMap, doublereal dCoef, FunctionCall func) const
    {
        G T, dT_dt, drho_dp, drho_dT;

        GetExtPressure(oStateCurr.p, oStateCurr.dp_dt, pDofMap, dCoef);
        GetTemperature(T, pDofMap, dCoef);
        GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
        pGetFluid()->GetDensity(oStateCurr.p, T, oStateCurr.rho, &drho_dp, &drho_dT);
        oStateCurr.eCavitationState = pGetFluid()->Cavitation(oStateCurr.p, &oStateCurr.dp_dt);

        oStateCurr.drho_dt = drho_dp * oStateCurr.dp_dt + drho_dT * dT_dt;
    }

    ComplianceModel::ComplianceModel(HydroMesh* pMesh,
                                     doublereal dDefScale,
                                     doublereal dPressScale)
        :HydroElement(pMesh, COMPLIANCE_ELEM),
         dDefScale(dDefScale),
         dPressScale(dPressScale),
         bDoInitAss(false)
    {
        HYDRO_ASSERT(pMesh->iGetNumNodes() > 0);

        rgNodes.reserve(pMesh->iGetNumNodes());
    }

    ComplianceModel::~ComplianceModel()
    {
    }

    int ComplianceModel::iGetNumNodes() const
    {
        return rgNodes.size();
    }

    void ComplianceModel::SetNode(int iNode, HydroNode* pNode)
    {
        HYDRO_ASSERT(size_t(iNode) == rgNodes.size() + 1);

        HydroUpdatedNode* pUpdatedNode = dynamic_cast<HydroUpdatedNode*>(pNode);

        HYDRO_ASSERT(pUpdatedNode != nullptr);
        HYDRO_ASSERT(rgNodes.capacity() >= rgNodes.size() + 1);

        rgNodes.push_back(pUpdatedNode);
    }

    HydroNode* ComplianceModel::pGetNode(int iNode) const
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(size_t(iNode) < rgNodes.size());

        return rgNodes[iNode];
    }

    integer ComplianceModel::iGetFirstIndex(grad::FunctionCall eFunc) const
    {
        HYDRO_ASSERT(unsigned(iGetOffsetIndex(eFunc)) <= pGetMesh()->pGetParent()->iGetNumDof());

        return pGetMesh()->pGetParent()->iGetFirstIndex() + iGetOffsetIndex(eFunc);
    }

    void ComplianceModel::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
    }

    void ComplianceModel::AfterConvergence(const VectorHandler& X,
                                           const VectorHandler& XP)
    {
    }

    void ComplianceModel::Initialize()
    {
        const HydroRootElement* pRootElem = pGetMesh()->pGetParent();
        const integer iTotNumElem = pRootElem->iGetNumElements();

        struct NodeLayout {
            index_type rgNodes[4];
        };

        static const NodeLayout rgQuad4[1] = {
            {0, 1, 2, 3},
        };

        static const NodeLayout rgQuad9[4] = {
            {0, 4, 8, 7},
            {4, 1, 5, 8},
            {8, 5, 2, 6},
            {7, 8, 6, 3}
        };

        integer iMaxNumElem = 0;

        for (integer i = 0; i < iTotNumElem; ++i) {
            const HydroElement* pElem = pRootElem->pGetElement(i);

            if (pElem->GetElementType() == FRICTION_ELEM) {
                switch (pElem->iGetNumNodes()) {
                case 4:
                    iMaxNumElem += sizeof(rgQuad4) / sizeof(rgQuad4[0]);
                    break;
                case 9:
                    iMaxNumElem += sizeof(rgQuad9) / sizeof(rgQuad9[0]);
                    break;
                }
            }
        }

        rgElements.reserve(iMaxNumElem);

        PressureElement oElem;

        for (integer i = 0; i < iTotNumElem; ++i) {
            const HydroElement* pElem = pRootElem->pGetElement(i);

            if (pElem->GetElementType() == FRICTION_ELEM) {
                const index_type iNumNodes = pElem->iGetNumNodes();

                decltype(std::begin(rgQuad4)) pBegQuad, pEndQuad;

                switch (iNumNodes){
                case 4:
                    pBegQuad = std::begin(rgQuad4);
                    pEndQuad = std::end(rgQuad4);
                    break;

                case 9:
                    pBegQuad = std::begin(rgQuad9);
                    pEndQuad = std::end(rgQuad9);
                    break;

                default:
                    throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
                }

                for (auto i = pBegQuad; i != pEndQuad; ++i) {
                    for (index_type j = 0; j < oElem.iGetNumNodes(); ++j) {
                        oElem.SetNode(j, pElem->pGetNode(i->rgNodes[j]));
                    }

                    rgElements.push_back(oElem);
                }
            }
        }
    }

    void
    ComplianceModel::GetRadialDeformation1(doublereal& w1, const HydroUpdatedNode* pNode) const
    {
        w1 = 0.;
    }

    void
    ComplianceModel::GetRadialDeformation2(doublereal& w2, const HydroUpdatedNode* pNode) const
    {
        w2 = 0.;
    }

    ComplianceModelNodal::ComplianceModelNodal(HydroMesh* pMesh,
                                               const Modal* pModalJoint,
                                               doublereal dDefScale,
                                               doublereal dPressScale,
                                               ComplianceMatrixArray&& rgMatArg)
        :ComplianceModel(pMesh, dDefScale, dPressScale),
         iNumNodes(-1), iNumModes(-1),
         pModalJoint(pModalJoint),
         rgMatrices(std::move(rgMatArg)),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {
    }

    ComplianceModelNodal::~ComplianceModelNodal()
    {

    }

    void
    ComplianceModelNodal::AssRes(SubVectorHandler& WorkVec,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    template <typename T>
    void ComplianceModelNodal::AssRes(GradientAssVec<T>& WorkVec,
                                      doublereal dCoef,
                                      const GradientVectorHandler<T>& XCurr,
                                      const GradientVectorHandler<T>& XPrimeCurr,
                                      enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename T>
    void ComplianceModelNodal::InitialAssRes(GradientAssVec<T>& WorkVec,
                                             const GradientVectorHandler<T>& XCurr,
                                             enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    void
    ComplianceModelNodal::AssJac(SparseSubMatrixHandler& WorkMat,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &rgDofMap[DOF_PRESS],
                                             mode);
    }

    void
    ComplianceModelNodal::InitialAssRes(SubVectorHandler& WorkVec,
                                        const VectorHandler& XCurr,
                                        GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<doublereal>::InitialAssRes(this,
                                                      WorkVec,
                                                      XCurr,
                                                      INITIAL_ASS_RES,
                                                      mode);
        }
    }

    void
    ComplianceModelNodal::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                        const VectorHandler& XCurr,
                                        GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                        WorkMat,
                                                        XCurr,
                                                        INITIAL_ASS_JAC,
                                                        &rgDofMap[DOF_PRESS],
                                                        mode);
        }
    }
    
    template <typename T>
    void ComplianceModelNodal::UnivAssRes(GradientAssVec<T>& WorkVec,
                                          doublereal dCoef,
                                          const GradientVectorHandler<T>& XCurr,
                                          enum FunctionCall func)
    {
        HYDRO_ASSERT(iNumNodes > 0);
        HYDRO_ASSERT(iNumModes >= 0);
        HYDRO_ASSERT(C.iGetNumRows() == iNumNodes || C.iGetNumRows() == 0);
        HYDRO_ASSERT(C.iGetNumCols() == iNumNodes || C.iGetNumCols() == 0);
        
        rgDofMap[DOF_DEF].Reset(func);
        
        Vector<T, DYNAMIC_SIZE> f1(iNumNodes), f2(E.iGetNumRows());
        T dw_dt;

        const doublereal dEquationScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_ELASTICITY_EQU);

        for (index_type i = 1; i <= f1.iGetNumRows(); ++i) {
            GetRadialDeformation(f1(i), dw_dt, dCoef, func, &rgDofMap[DOF_DEF], rgNodes[i - 1]);
            f1(i) *= dEquationScale;
            HYDRO_TRACE("f1(" << i << ")=" << f1(i) << std::endl);
        }

        if (pModalJoint) {
            HYDRO_ASSERT(pModalJoint->uGetNModes() == E.iGetNumRows());
            HYDRO_ASSERT(pModalJoint->uGetNModes() == D.iGetNumCols());
            HYDRO_ASSERT(iNumNodes == E.iGetNumCols());
            HYDRO_ASSERT(iNumNodes == D.iGetNumRows());
            
            for (index_type i = 1; i <= f2.iGetNumRows(); ++i) {
                SetDofMap(f2(i), &rgDofMap[DOF_DEF]);
            }
        }

        integer iEqIndex = iGetFirstIndex(func);

        T p, pasp;

        for (index_type j = 1; j <= iNumNodes; ++j) {
            rgDofMap[DOF_PRESS].Reset(); // reduces the size of the local vector stored in p

            pGetMesh()->GetPressure(rgNodes[j - 1], p, &rgDofMap[DOF_PRESS], dCoef);

            if (rgNodes[j - 1]->GetContactPressure(pasp, &rgDofMap[DOF_PRESS])) {
                p += pasp;
            }

            p *= dPressScale;

            HYDRO_TRACE("p=" << p << "\n");

            for (index_type i = 1; i <= C.iGetNumRows(); ++i) {
                HYDRO_TRACE("before: f1(" << i << ")=" <<  f1(i) << "\n");

                f1(i) -= dEquationScale * C(i, j) * p;

                HYDRO_TRACE("after: f1(" << i << ")=" <<  f1(i) << "\n");
            }
            
            if (pModalJoint) {
                for (index_type i = 1; i <= E.iGetNumRows(); ++i) {
                    f2(i) += E(i, j) * p;
                }
            }
        }

        if (pModalJoint) {
            HYDRO_ASSERT(pModalJoint->uGetNModes() == D.iGetNumCols());
            HYDRO_ASSERT(D.iGetNumRows() == f1.iGetNumRows());

            T aj;
            
            for (index_type j = 1; j <= D.iGetNumCols(); ++j) {
                rgDofMap[DOF_PRESS].Reset(); // reduces the size of the local vector stored in a(j)
                pModalJoint->GetACurr(j, aj, dCoef, func, &rgDofMap[DOF_PRESS]);
                
                for (index_type i = 1; i <= f1.iGetNumRows(); ++i) {
                    f1(i) -= dEquationScale * D(i, j) * aj;
                }
            }
        }

        for (index_type i = 1; i <= f1.iGetNumRows(); ++i) {
            CHECK_NUM_COLS_WORK_SPACE(this, func, f1(i), iEqIndex);
            
            WorkVec.AddItem(iEqIndex++, f1(i));

            HYDRO_TRACE("f1(" << i << ")=" << f1(i) << std::endl);
        }

        if (pModalJoint) {
            iEqIndex = pModalJoint->iGetFirstIndex() + pModalJoint->uGetNModes() + 1;
            
            for (index_type i = 1; i <= f2.iGetNumRows(); ++i) {
                CHECK_NUM_COLS_WORK_SPACE(this, func, f2(i), iEqIndex);
                
                WorkVec.AddItem(iEqIndex++, f2(i));

                HYDRO_TRACE("f2(" << i << ")=" << f2(i) << std::endl);
            }
        }
    }

    void ComplianceModelNodal::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = iNumNodes;
        *piNumCols = iNumNodes + (eFunc & grad::REGULAR_FLAG ? iGetNumDof() : iGetInitialNumDof());

        if (pModalJoint) {
            *piNumRows += E.iGetNumRows();
            *piNumCols += pModalJoint->uGetNModes();
        }

        *piNumCols += pGetMesh()->pGetGeometry()->iGetNumColsWorkSpace(eFunc);
    }

    void ComplianceModelNodal::Initialize()
    {
        ComplianceModel::Initialize();

        iNumNodes = iGetNumNodes();
        iNumModes = pModalJoint ? pModalJoint->uGetNModes() : 0;

        w.Resize(iNumNodes);
        dw_dt.Resize(iNumNodes);

        ComplianceMatrix::MatrixData oMatData{ComplianceMatrix::LOC_MESH_FIXED,
                pModalJoint,
                &C,
                &D,
                &E};

        for (index_type i = 0; i < 2; ++i) {
            if (rgMatrices[i].get()) {
                rgMatrices[i]->AddCompliance(oMatData,
                                             pGetMesh(),
                                             rgElements,
                                             rgNodes,
                                             dPressScale);
                HYDRO_ASSERT(iGetNumNodes() == static_cast<ssize_t>(oMatData.rgActGridIdx.size()));
            }
        }

#if HYDRO_DEBUG > 0
        if (pModalJoint) {
            HYDRO_ASSERT(D.iGetNumRows() == iNumNodes);
            HYDRO_ASSERT(D.iGetNumCols() == iNumModes);
            HYDRO_ASSERT(E.iGetNumRows() == iNumModes);
            HYDRO_ASSERT(E.iGetNumCols() == iNumNodes);
        }

        if (C.iGetNumRows() == 0 && C.iGetNumCols() == 0) {
            HYDRO_ASSERT(pModalJoint);
        } else {
            HYDRO_ASSERT(C.iGetNumRows() == iNumNodes);
            HYDRO_ASSERT(C.iGetNumCols() == iNumNodes);
        }
#endif
    }

    void ComplianceModelNodal::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);

        eCurrFunc = grad::REGULAR_FLAG;

        integer iDofIndex = iGetFirstIndex(eCurrFunc);

        for (index_type i = 1; i <= w.iGetNumRows(); ++i, ++iDofIndex) {
            HYDRO_ASSERT(i <= XCurr.iGetSize());

            XCurr(iDofIndex) = w(i) / dDefScale;
            XPrimeCurr(iDofIndex) = dw_dt(i) / dDefScale;
        }
    }

    unsigned int ComplianceModelNodal::iGetNumDof() const
    {
        return rgNodes.size();
    }

    unsigned int ComplianceModelNodal::iGetInitialNumDof() const
    {
        return bDoInitAss ? iGetNumDof() : 0u;
    }

    integer ComplianceModelNodal::iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const
    {
        return (eFunc & grad::REGULAR_FLAG ? 1 : bDoInitAss) * iNumNodes; // Number of columns per node in this case
    }

    DofOrder::Order ComplianceModelNodal::GetDofType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::DIFFERENTIAL;
    }

    DofOrder::Order ComplianceModelNodal::GetEqType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::ALGEBRAIC;
    }

    std::ostream&
    ComplianceModelNodal::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iDofIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

            for (index_type i = 1; i <= w.iGetNumRows(); ++i, ++iDofIndex) {
                out << prefix << iDofIndex << ": Elasticity dof #" << i << " of node number " << rgNodes[i - 1]->iGetNodeNumber() + 1 << std::endl;
            }
        }

        return out;
    }

    std::ostream&
    ComplianceModelNodal::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iDofIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

            for (index_type i = 1; i <= w.iGetNumRows(); ++i, ++iDofIndex) {
                out << prefix << iDofIndex << ": Elasticity definition equation #" << i << " for node number " << rgNodes[i - 1]->iGetNodeNumber() + 1 << std::endl;
            }
        }

        return out;
    }

    void ComplianceModelNodal::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr, doublereal dCoef, FunctionCall func)
    {
        if ((func & grad::REGULAR_FLAG) || bDoInitAss) {
            integer iDofIndex = iGetFirstIndex(func);

            for (index_type i = 1; i <= w.iGetNumRows(); ++i, ++iDofIndex) {
                HYDRO_ASSERT(iDofIndex >= 1);
                HYDRO_ASSERT(iDofIndex <= XCurr.iGetSize());

                w(i) = XCurr(iDofIndex) * dDefScale;

                if (func & grad::REGULAR_FLAG) {
                    HYDRO_ASSERT(iDofIndex <= XPrimeCurr.iGetSize());
                    dw_dt(i) = XPrimeCurr(iDofIndex) * dDefScale;
                }
            }
        }
    }

    void ComplianceModelNodal::GetRadialDeformation(doublereal& wi, doublereal& dwi_dt, doublereal, FunctionCall, LocalDofMap*, const HydroUpdatedNode* pNode) const
    {
        const index_type iCompIndex = pNode->iGetComplianceIndex();

        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= iGetNumNodes());
        HYDRO_ASSERT(rgNodes[iCompIndex - 1] == pNode);

        wi = w(iCompIndex);
        dwi_dt = dw_dt(iCompIndex);
    }

    void ComplianceModelNodal::GetRadialDeformation(Gradient<0>& wi, Gradient<0>& dwi_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap, const HydroUpdatedNode* pNode) const
    {
        const index_type iCompIndex = pNode->iGetComplianceIndex();

        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= iGetNumNodes());
        HYDRO_ASSERT(rgNodes[iCompIndex - 1] == pNode);

        if ((func & grad::REGULAR_FLAG) || bDoInitAss) {
            index_type iDofIndex = iGetFirstIndex(func) + iCompIndex - 1;

            wi.SetValuePreserve(w(iCompIndex));

            wi.DerivativeResizeReset(pDofMap,
                                     iDofIndex,
                                     MapVectorBase::GLOBAL,
                                     -dCoef * dDefScale);

            if (func & grad::REGULAR_FLAG) {
                dwi_dt.SetValuePreserve(dw_dt(iCompIndex));

                dwi_dt.DerivativeResizeReset(pDofMap,
                                             iDofIndex,
                                             MapVectorBase::GLOBAL,
                                             -dDefScale);
            } else {
                HYDRO_ASSERT(dCoef == 1);
                dwi_dt.SetValue(dw_dt(iCompIndex));
            }
        } else {
            wi.SetValue(w(iCompIndex));
            dwi_dt.SetValue(dw_dt(iCompIndex));
        }
    }

    void ComplianceModelNodal::Print(std::ostream& os) const
    {
        HydroElement::Print(os);

        os << C.iGetNumRows() << '\n';

        for (index_type i = 1; i <= C.iGetNumRows(); ++i) {
            for (index_type j = 1; j <= C.iGetNumCols(); ++j) {
                os << C(i, j) << '\t';
            }
            os << '\n';
        }
    }

    const std::array<integer, ComplianceModelNodalDouble::POLYORDER> ComplianceModelNodalDouble::px = {3, 2, 1, 0, 2, 1, 0, 1, 0};
    const std::array<integer, ComplianceModelNodalDouble::POLYORDER> ComplianceModelNodalDouble::pz = {0, 1, 2, 3, 0, 1, 2, 0, 1};
    const std::array<index_type, ComplianceModelNodalDouble::GRIDINTERP> ComplianceModelNodalDouble::xg = {-1, -1, -1, -1, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2};
    const std::array<index_type, ComplianceModelNodalDouble::GRIDINTERP> ComplianceModelNodalDouble::zg = {-1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2, -1, 0, 1, 2};
    const index_type ComplianceModelNodalDouble::min_xg = *std::min_element(xg.begin(), xg.end());
    const index_type ComplianceModelNodalDouble::max_xg = *std::max_element(xg.begin(), xg.end());
    const index_type ComplianceModelNodalDouble::min_zg = *std::min_element(zg.begin(), zg.end());
    const index_type ComplianceModelNodalDouble::max_zg = *std::max_element(zg.begin(), zg.end());

    ComplianceModelNodalDouble::ComplianceModelNodalDouble(HydroMesh* pMesh,
                                                           const ModalJointArray& rgModalJoints,
                                                           doublereal dDefScale,
                                                           doublereal dPressScale,
                                                           ComplianceMatrixArray&& rgMatrices,
                                                           const CylindricalBearing& oGeometry,
                                                           DEhdInterpolOption eInterpolOption)
    :ComplianceModel(pMesh, dDefScale, dPressScale),
     rgModalJoints(rgModalJoints),
     dPressDofScale(0.),
     dMeshRadius(oGeometry.dGetMeshRadius()),
     dMinDistance_2(std::pow(std::numeric_limits<doublereal>::epsilon(), 2./6.)),
     rgMatrices(std::move(rgMatrices)),
     eCurrFunc(grad::INITIAL_ASS_FLAG),
     eInterpolOption(eInterpolOption),
     dAxialThreshold(0.)
    {
        std::fill(rgNumNodes.begin(), rgNumNodes.end(), -1);
    }

    ComplianceModelNodalDouble::~ComplianceModelNodalDouble()
    {
    }

    template <typename T>
    void ComplianceModelNodalDouble::AssRes(GradientAssVec<T>& WorkVec,
                                            doublereal dCoef,
                                            const GradientVectorHandler<T>& XCurr,
                                            const GradientVectorHandler<T>& XPrimeCurr,
                                            enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    void
    ComplianceModelNodalDouble::AssRes(SubVectorHandler& WorkVec,
                                       doublereal dCoef,
                                       const VectorHandler& XCurr,
                                       const VectorHandler& XPrimeCurr,
                                       GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    ComplianceModelNodalDouble::AssJac(SparseSubMatrixHandler& WorkMat,
                                       doublereal dCoef,
                                       const VectorHandler& XCurr,
                                       const VectorHandler& XPrimeCurr,
                                       GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &rgDofMap[DOF_GLOBAL1],
                                             mode);
    }

    template <typename T>
    void ComplianceModelNodalDouble::InitialAssRes(GradientAssVec<T>& WorkVec,
                                                   const GradientVectorHandler<T>& XCurr,
                                                   enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    void
    ComplianceModelNodalDouble::InitialAssRes(SubVectorHandler& WorkVec,
                                              const VectorHandler& XCurr,
                                              GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<doublereal>::InitialAssRes(this,
                                                      WorkVec,
                                                      XCurr,
                                                      INITIAL_ASS_RES,
                                                      mode);
        }
    }

    void
    ComplianceModelNodalDouble::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                              const VectorHandler& XCurr,
                                              GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                        WorkMat,
                                                        XCurr,
                                                        INITIAL_ASS_JAC,
                                                        &rgDofMap[DOF_GLOBAL1],
                                                        mode);
        }
    }

    void ComplianceModelNodalDouble::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        if ((eFunc & grad::REGULAR_FLAG) || bDoInitAss) {
            *piNumRows = rgNumNodes[DEHD_BODY_FIXED] + rgNumNodes[DEHD_BODY_MOVING];
            *piNumCols = rgNumNodes[DEHD_BODY_FIXED] + (eFunc & grad::REGULAR_FLAG ? iGetNumDof() : iGetInitialNumDof());

            for (index_type i = 0; i < DEHD_BODY_LAST; ++i) {
                if (rgModalJoints[i]) {
                    *piNumRows += rgModalJoints[i]->uGetNModes();
                    *piNumCols += rgModalJoints[i]->uGetNModes();
                }
            }

            *piNumCols += pGetMesh()->pGetGeometry()->iGetNumColsWorkSpace(eFunc);
        } else {
            *piNumRows = *piNumCols = 0;
        }
    }

    void ComplianceModelNodalDouble::Initialize()
    {
        ComplianceModel::Initialize();

        dPressDofScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_PRESSURE_DOF);

        const std::array<ComplianceMatrix::MeshLocation, DEHD_BODY_LAST> rgMeshLoc = {ComplianceMatrix::LOC_MESH_FIXED,
                                                                                      ComplianceMatrix::LOC_MESH_MOVING};
        for (index_type i = 0; i < DEHD_BODY_LAST; ++i) {
            ComplianceMatrix::MatrixData oMatData{rgMeshLoc[i], rgModalJoints[i], &C[i], &D[i], &E[i]};

            if (!rgMatrices[i].get()) {
                HYDRO_ASSERT(0);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            rgMatrices[i]->AddCompliance(oMatData,
                                         pGetMesh(),
                                         rgElements,
                                         rgNodes,
                                         dPressScale);

            xi[i] = std::move(oMatData.rgGridX);
            zi[i] = std::move(oMatData.rgGridZ);
            rgMatIdx[i] = std::move(oMatData.rgMatIdx);
            rgActGridIdx[i] = std::move(oMatData.rgActGridIdx);
            rgNumNodes[i] = rgActGridIdx[i].size();
            
#ifdef HYDRO_DEBUG
            if (rgModalJoints[i]) {
                HYDRO_ASSERT(D[i].iGetNumRows() == rgNumNodes[i]);
                HYDRO_ASSERT(D[i].iGetNumCols() == rgModalJoints[i]->uGetNModes());
                HYDRO_ASSERT(E[i].iGetNumRows() == rgModalJoints[i]->uGetNModes());
                HYDRO_ASSERT(E[i].iGetNumCols() == rgNumNodes[i]);
            }

            HYDRO_ASSERT(C[i].iGetNumRows() == rgNumNodes[i] || C[i].iGetNumRows() == 0);
            HYDRO_ASSERT(C[i].iGetNumCols() == rgNumNodes[i] || C[i].iGetNumCols() == 0);
            HYDRO_ASSERT((C[i].iGetNumRows() > 0 && C[i].iGetNumCols() > 0) || rgModalJoints[i] != nullptr);
#endif
        }

        if (eInterpolOption == INT_AXIAL_LARGE_DISP) {
            HYDRO_ASSERT(zi[DEHD_BODY_MOVING].size() > 1);
            
            dAxialThreshold = (zi[DEHD_BODY_MOVING].back() - zi[DEHD_BODY_MOVING].front()) / (zi[DEHD_BODY_MOVING].size() - 1);
        }
        
        static const DEhdBodyIdx rgBodyIndexDef[] = {
            DEHD_BODY_FIXED,
            DEHD_BODY_MOVING,
            DEHD_BODY_FIXED
        };
        
        for (index_type i = DEHD_DEF_TOTAL; i <= DEHD_DEF_MOVING_INTERP; ++i) {
            w[i].Resize(rgNumNodes[rgBodyIndexDef[i]]);
        }

        for (index_type i = DEHD_DEF_TOTAL; i <= DEHD_DEF_MOVING; ++i) {
            dw_dt[i].Resize(rgNumNodes[rgBodyIndexDef[i]]);
        }

        HYDRO_ASSERT(!rgModalJoints[DEHD_BODY_MOVING] || D[DEHD_BODY_MOVING].iGetNumCols() > 0);
        HYDRO_ASSERT(!rgModalJoints[DEHD_BODY_FIXED] || D[DEHD_BODY_FIXED].iGetNumCols() > 0);
        HYDRO_ASSERT(!rgModalJoints[DEHD_BODY_MOVING] || E[DEHD_BODY_MOVING].iGetNumRows() > 0);
        HYDRO_ASSERT(!rgModalJoints[DEHD_BODY_FIXED] || E[DEHD_BODY_FIXED].iGetNumRows() > 0);

        for (index_type k = 0; k < DEHD_BODY_LAST; ++k) {
            HYDRO_ASSERT(xi[k].size() >= 2);
            HYDRO_ASSERT(zi[k].size() >= 2);
            const index_type nx = xi[k].size();
            const index_type nz = zi[k].size();

            HYDRO_ASSERT(rgPolyData[k].size() == 0);

            rgPolyData[k].reserve(nx * nz);

            for (index_type i = 0; i < nx; ++i) {
                index_type ox = i;

                if (ox + max_xg >= nx) {
                    ox = nx - max_xg - 1;
                }

                if (ox + min_xg < 0) {
                    ox = -min_xg;
                }

                if (!((ox + min_xg >= 0) && (ox + max_xg < nx))) {
                    silent_cerr("hydrodynamic plain bearing2(" << pGetMesh()->pGetParent()->GetLabel()
                                << "): not enough nodes provided in x-direction for interpolation" << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                for (index_type j = 0; j < nz; ++j) {
                    HYDRO_ASSERT(static_cast<ssize_t>(rgPolyData[k].size()) == i * nz + j);

                    index_type oz = j;

                    if (oz + max_zg >= nz) {
                        oz = nz - max_zg - 1;
                    }

                    if (oz + min_zg < 0) {
                        oz = -min_zg;
                    }

                    if (!((oz + min_zg >= 0) && (oz + max_zg < nz))) {
                        silent_cerr("hydrodynamic plain bearing2(" << pGetMesh()->pGetParent()->GetLabel()
                                    << "): not enough nodes provided in z-direction for interpolation" << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    Matrix<doublereal, GRIDINTERP, POLYORDER> A;
                    std::array<index_type, GRIDINTERP> Cidx;

                    for (index_type l = 0; l < GRIDINTERP; ++l) {
                        const index_type xgl = xg[l] + ox;
                        const index_type zgl = zg[l] + oz;

                        HYDRO_ASSERT(xgl >= 0);
                        HYDRO_ASSERT(xgl < nx);
                        HYDRO_ASSERT(zgl >= 0);
                        HYDRO_ASSERT(zgl < nz);

                        const doublereal dxl = xi[k][xgl] - xi[k][i];
                        const doublereal dzl = zi[k][zgl] - zi[k][j];

                        Cidx[l] = rgMatIdx[k][xgl * nz + zgl];

                        for (index_type m = 0; m < POLYORDER; ++m) {
                            A(l + 1, m + 1) = std::pow(dxl, px[m]) * std::pow(dzl, pz[m]);
                        }
                    }

                    Matrix<doublereal, POLYORDER, POLYORDER> A_TA = Transpose(A) * A;

                    integer INFO;
                    std::array<integer, POLYORDER> IPIV;

                    const integer M = A_TA.iGetNumCols();
                    const integer N = A_TA.iGetNumRows();

                    __FC_DECL__(dgetrf)(&M, &N, &A_TA(1, 1), &M, &IPIV[0], &INFO);

                    if (INFO != 0) {
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    const integer NRHS = A.iGetNumRows();

                    __FC_DECL__(dgetrs)("N", &N, &NRHS, &A_TA(1, 1), &M, &IPIV[0], &A(1, 1), &N, &INFO);

                    if (INFO != 0) {
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    rgPolyData[k].emplace_back(Transpose(A), Cidx);

                    HYDRO_ASSERT(static_cast<ssize_t>(rgPolyData[k].size()) == i * nz + j + 1);
                }
            }
        }

#if HYDRO_DEBUG > 1
        if (pedantic_output) {
            for (size_t i = 0; i < xi[DEHD_BODY_FIXED].size(); ++i) {
                HYDRO_TRACE("xi[DEHD_BODY_FIXED][" << i << "]=" << xi[DEHD_BODY_FIXED][i] << std::endl);
            }

            for (size_t i = 0; i < zi[DEHD_BODY_FIXED].size(); ++i) {
                HYDRO_TRACE("zi[DEHD_BODY_FIXED][" << i << "]=" << zi[DEHD_BODY_FIXED][i] << std::endl);
            }

            for (size_t i = 0; i < rgMatIdx[DEHD_BODY_FIXED].size(); ++i) {
                HYDRO_TRACE("rgMatIdx[DEHD_BODY_FIXED][" << i << "]=" << rgMatIdx[DEHD_BODY_FIXED][i] << std::endl);
            }

            HYDRO_TRACE("rgMatIdx[DEHD_BODY_FIXED]:\n");

            HYDRO_TRACE(std::setw(8) << "row/col" << ' ');

            for (size_t j = 0; j < zi[DEHD_BODY_FIXED].size(); ++j) {
                HYDRO_TRACE(std::setw(8) << j << ' ');
            }

            HYDRO_TRACE(std::endl);

            for (size_t i = 0; i < xi[DEHD_BODY_FIXED].size(); ++i) {
                HYDRO_TRACE(std::setw(8) << i + 1 << ' ');
                for (size_t j = 0; j < zi[DEHD_BODY_FIXED].size(); ++j) {
                    HYDRO_TRACE(std::setw(8) << rgMatIdx[DEHD_BODY_FIXED][i * zi[DEHD_BODY_FIXED].size() + j] << ' ');
                }
                HYDRO_TRACE(std::endl);
            }

            for (index_type i = 0; i < iNumNodes; ++i) {
                const Node2D* pNode = pRootElem->pGetNode(i);
                const HydroNode* const pHydrNode = dynamic_cast<const HydroNode*>(pNode);
                index_type iCompIndex = pHydrNode ? pHydrNode->iGetComplianceIndex() : -1;
                HYDRO_TRACE("node #" << pNode->iGetNodeNumber() << "(" << pNode->GetPosition2D()(1) << "," << pNode->GetPosition2D()(2) << "):" << iCompIndex << std::endl);
            }
        }
#endif
    }

    void ComplianceModelNodalDouble::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);

        eCurrFunc = grad::REGULAR_FLAG;

        integer iDofIndex = iGetFirstIndex(eCurrFunc);

        for (index_type k = DEHD_DEF_TOTAL; k <= DEHD_DEF_MOVING; ++k) {
            for (index_type i = 1; i <= w[k].iGetNumRows(); ++i, ++iDofIndex) {
                HYDRO_ASSERT(iDofIndex <= XCurr.iGetSize());

                XCurr(iDofIndex) = w[k](i) / dDefScale;
                XPrimeCurr(iDofIndex) = dw_dt[k](i) / dDefScale;
            }
        }
    }

    unsigned int ComplianceModelNodalDouble::iGetNumDof(void) const
    {
        return rgNumNodes[DEHD_BODY_FIXED] + rgNumNodes[DEHD_BODY_MOVING];
    }

    unsigned int ComplianceModelNodalDouble::iGetInitialNumDof(void) const
    {
        return bDoInitAss ? iGetNumDof() : 0u;
    }

    integer ComplianceModelNodalDouble::iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const
    {
        return (eFunc & grad::REGULAR_FLAG ? 1 : bDoInitAss) * iNumNodes; // Number of columns per node in this case
    }

    DofOrder::Order ComplianceModelNodalDouble::GetDofType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::DIFFERENTIAL;
    }

    DofOrder::Order ComplianceModelNodalDouble::GetEqType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::ALGEBRAIC;
    }

    std::ostream&
    ComplianceModelNodalDouble::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iDofIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

            static const char rgDesc[][39] = {" for total deformation",
                                              " for deformation opposite to mesh side"};
            for (index_type k = 0; k < DEHD_BODY_LAST; ++k) {
                for (index_type i = 1; i <= w[k].iGetNumRows(); ++i, ++iDofIndex) {
                    out << prefix << iDofIndex << ": elasticity dof #" << i << rgDesc[k] << std::endl;
                }
            }
        }

        return out;
    }

    std::ostream&
    ComplianceModelNodalDouble::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iEqIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);

            static const char rgDesc[][39] = {" for total deformation",
                                              " for deformation opposite to mesh side"};

            for (index_type k = 0; k < DEHD_BODY_LAST; ++k) {
                for (index_type i = 1; i <= w[k].iGetNumRows(); ++i, ++iEqIndex) {
                    out << prefix << iEqIndex << ": equation of elasticity #" << i << rgDesc[k] << std::endl;
                }
            }
        }

        return out;
    }

    void
    ComplianceModelNodalDouble::Update(const VectorHandler& XCurr,
                                       const VectorHandler& XPrimeCurr,
                                       doublereal dCoef,
                                       FunctionCall func)
    {
        switch (func) {
        case INITIAL_ASS_RES:
        case INITIAL_DER_RES:
        case REGULAR_RES:
            if ((func & grad::REGULAR_FLAG) || bDoInitAss) {
                integer iDofIndex = iGetFirstIndex(func);

                for (index_type k = DEHD_DEF_TOTAL; k <= DEHD_DEF_MOVING; ++k) {
                    for (index_type i = 1; i <= w[k].iGetNumRows(); ++i, ++iDofIndex) {
                        HYDRO_ASSERT(iDofIndex >= 1);
                        HYDRO_ASSERT(iDofIndex <= XCurr.iGetSize());

                        w[k](i) = XCurr(iDofIndex) * dDefScale;

                        if (func & grad::REGULAR_FLAG) {
                            HYDRO_ASSERT(iDofIndex <= XPrimeCurr.iGetSize());
                            dw_dt[k](i) = XPrimeCurr(iDofIndex) * dDefScale;
                        }
                    }
                }
            }
            break;

        default:
            break;
        }
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation(doublereal& wi,
                                                     doublereal& dwi_dt,
                                                     doublereal dCoef,
                                                     FunctionCall func,
                                                     LocalDofMap* pDofMap,
                                                     const HydroUpdatedNode* pNode) const
    {
        GetRadialDeformationTpl(wi, dwi_dt, dCoef, func, pDofMap, pNode);
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation(Gradient<0>& wi,
                                                     Gradient<0>& dwi_dt,
                                                     doublereal dCoef,
                                                     FunctionCall func,
                                                     LocalDofMap* pDofMap,
                                                     const HydroUpdatedNode* pNode) const
    {
        GetRadialDeformationTpl(wi, dwi_dt, dCoef, func, pDofMap, pNode);
    }

    template <typename G>
    void
    ComplianceModelNodalDouble::GetRadialDeformationTpl(G& wi,
                                                        G& dwi_dt,
                                                        doublereal dCoef,
                                                        FunctionCall func,
                                                        LocalDofMap* pDofMap,
                                                        const HydroUpdatedNode* pNode) const
    {
        const index_type iCompIndex = pNode->iGetComplianceIndex();

        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= iGetNumNodes());
        HYDRO_ASSERT(rgNodes[iCompIndex - 1] == pNode);

        GetRadialDeformation(wi, dwi_dt, dCoef, func, pDofMap, DEHD_DEF_TOTAL, iCompIndex);
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation(doublereal& wi,
                                                     doublereal& dwi_dt,
                                                     doublereal dCoef,
                                                     FunctionCall func,
                                                     LocalDofMap* pDofMap,
                                                     DEhdDeformationIdx eDefIndex,
                                                     index_type iCompIndex) const
    {
        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= w[eDefIndex].iGetNumRows());
        HYDRO_ASSERT(eDefIndex >= DEHD_DEF_TOTAL);
        HYDRO_ASSERT(eDefIndex <= DEHD_DEF_MOVING);

        wi = w[eDefIndex](iCompIndex);
        dwi_dt = dw_dt[eDefIndex](iCompIndex);
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation(Gradient<0>& wi,
                                                     Gradient<0>& dwi_dt,
                                                     doublereal dCoef,
                                                     FunctionCall func,
                                                     LocalDofMap* pDofMap,
                                                     DEhdDeformationIdx eDefIndex,
                                                     index_type iCompIndex) const
    {
        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= w[eDefIndex].iGetNumRows());
        HYDRO_ASSERT(eDefIndex >= DEHD_DEF_TOTAL);
        HYDRO_ASSERT(eDefIndex <= DEHD_DEF_MOVING);

        if ((func & grad::REGULAR_FLAG) || bDoInitAss) {
            index_type iDofIndex = iGetFirstIndex(func) + iCompIndex - 1;

            for (index_type k = 0; k < eDefIndex; ++k) {
                iDofIndex += w[k].iGetNumRows();
            }

            wi.SetValuePreserve(w[eDefIndex](iCompIndex));

            wi.DerivativeResizeReset(pDofMap,
                                     iDofIndex,
                                     MapVectorBase::GLOBAL,
                                     -dCoef * dDefScale);

            if (func & grad::REGULAR_FLAG) {
                dwi_dt.SetValuePreserve(dw_dt[eDefIndex](iCompIndex));

                dwi_dt.DerivativeResizeReset(pDofMap,
                                             iDofIndex,
                                             MapVectorBase::GLOBAL,
                                             -dDefScale);
            } else {
                HYDRO_ASSERT(dCoef == 1);
                dwi_dt.SetValue(dw_dt[eDefIndex](iCompIndex));
            }
        } else {
            wi.SetValue(w[eDefIndex](iCompIndex));
            dwi_dt.SetValue(dw_dt[eDefIndex](iCompIndex));
        }
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation1(doublereal& w1, const HydroUpdatedNode* pNode) const
    {
        const index_type iCompIndex = pNode->iGetComplianceIndex();

        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= w[DEHD_DEF_TOTAL].iGetNumRows());
        HYDRO_ASSERT(rgNodes[iCompIndex - 1] == pNode);

        switch (pGetMesh()->pGetGeometry()->GetType()) {
        case BearingGeometry::CYLINDRICAL_MESH_AT_SHAFT:
            w1 = w[DEHD_DEF_TOTAL](iCompIndex) - w[DEHD_DEF_MOVING_INTERP](iCompIndex);
            break;

        case BearingGeometry::CYLINDRICAL_MESH_AT_BEARING:
            w1 = w[DEHD_DEF_MOVING_INTERP](iCompIndex);
            break;

        default:
            HYDRO_ASSERT(false);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    void
    ComplianceModelNodalDouble::GetRadialDeformation2(doublereal& w2, const HydroUpdatedNode* pNode) const
    {
        const index_type iCompIndex = pNode->iGetComplianceIndex();

        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= w[DEHD_DEF_TOTAL].iGetNumRows());
        HYDRO_ASSERT(rgNodes[iCompIndex - 1] == pNode);

        switch (pGetMesh()->pGetGeometry()->GetType()) {
        case BearingGeometry::CYLINDRICAL_MESH_AT_SHAFT:
            w2 = w[DEHD_DEF_MOVING_INTERP](iCompIndex);
            break;

        case BearingGeometry::CYLINDRICAL_MESH_AT_BEARING:
            w2 = w[DEHD_DEF_TOTAL](iCompIndex) - w[DEHD_DEF_MOVING_INTERP](iCompIndex);
            break;

        default:
            HYDRO_ASSERT(false);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    template <typename IT, typename T>
    IT find_nearest(IT first, IT last, T x)
    {
        const IT i2 = std::lower_bound(first, last, x);

        if (i2 > first && i2 < last) {
            const IT i1 = i2 - 1;
            return fabs(x - *i1) < fabs(x - *i2) ? i1 : i2;
        } else {
            return i2 > first ? i2 - 1 : i2;
        }
    }

    template <ComplianceModelNodalDouble::DEhdFieldType eField,
              ComplianceModelNodalDouble::DEhdBodyIdx eMshSrc,
              ComplianceModelNodalDouble::DEhdBodyIdx eMshDst,
              typename T>
    void ComplianceModelNodalDouble::Interpolate(const index_type i,
                                                 const index_type j,
                                                 const Vector<T, 2>& dxm_g,
                                                 const doublereal dScale,
                                                 T& fij_f,
                                                 const doublereal dCoef,
                                                 const FunctionCall func,
                                                 LocalDofMap* const pDofMap) const {
        
        const Vector<T, 2> dxm(dxm_g, pDofMap);
        
        T alpha(1.);
        
        static_assert(eMshSrc != eMshDst, "source mesh must be different from destination mesh");
        static_assert((eMshSrc == DEHD_BODY_FIXED && eField == FT_TOTAL_PRESS) || (eMshSrc == DEHD_BODY_MOVING && eField == FT_DEF_MOVING), "interpolation not supported");

        const std::array<doublereal, DEHD_BODY_LAST> rgDirInterp = {-1., 1};

        HYDRO_ASSERT(fabs(rgDirInterp[eMshDst]) == 1.);

        T x_f = xi[eMshDst][i] + rgDirInterp[eMshDst] * dxm(1);

        if (x_f < xi[eMshSrc].front()) {
            x_f += 2. * dMeshRadius * M_PI;
        } else if (x_f > xi[eMshSrc].back()) {
            x_f -= 2. * dMeshRadius * M_PI;
        }

        HYDRO_ASSERT(x_f >= xi[eMshSrc].front());
        HYDRO_ASSERT(x_f <= xi[eMshSrc].back());

        const index_type i1_f = find_nearest(xi[eMshSrc].cbegin(), xi[eMshSrc].cend(), x_f) - xi[eMshSrc].cbegin();

        HYDRO_ASSERT(i1_f >= 0);
        HYDRO_ASSERT(i1_f < static_cast<ssize_t>(xi[eMshSrc].size()));

        index_type i2_f;

        if (x_f < xi[eMshSrc][i1_f]) {
            i2_f = (i1_f - 1 >= 0) ? i1_f - 1 : i1_f + 1;
        } else {
            i2_f = (i1_f + 1 < static_cast<ssize_t>(xi[eMshSrc].size())) ? i1_f + 1 : i1_f - 1;
        }

        HYDRO_ASSERT(i2_f >= 0);
        HYDRO_ASSERT(i2_f < static_cast<ssize_t>(xi[eMshSrc].size()));

        T z_f = zi[eMshDst][j] + rgDirInterp[eMshDst] * dxm(2);

        HYDRO_ASSERT(zi[eMshSrc].size());
        HYDRO_ASSERT(zi[eMshSrc].back() > zi[eMshSrc].front());

        if (eInterpolOption == INT_AXIAL_LARGE_DISP) {            
            if (z_f < zi[eMshSrc].front() - dAxialThreshold || z_f > zi[eMshSrc].back() + dAxialThreshold) {
                fij_f = 0.;
                return;
            } else if (z_f < zi[eMshSrc].front() + dAxialThreshold) {
                alpha = (z_f - (zi[eMshSrc].front() - dAxialThreshold)) / (2. * dAxialThreshold);
            } else if (z_f > zi[eMshSrc].back() - dAxialThreshold) {
                alpha = ((zi[eMshSrc].back() + dAxialThreshold) - z_f) / (2. * dAxialThreshold);               
            }
        } else if (eInterpolOption == INT_AXIAL_SMALL_DISP) {
            if (z_f < zi[eMshSrc].front()) {
                z_f = zi[eMshSrc].front();
            } else if (z_f > zi[eMshSrc].back()) {
                z_f = zi[eMshSrc].back();
            }
        }

        HYDRO_ASSERT(alpha >= 0.);
        HYDRO_ASSERT(alpha <= 1.);        

        const index_type j1_f = find_nearest(zi[eMshSrc].cbegin(), zi[eMshSrc].cend(), z_f) - zi[eMshSrc].cbegin();

        HYDRO_ASSERT(j1_f >= 0);
        HYDRO_ASSERT(j1_f < static_cast<ssize_t>(zi[eMshSrc].size()));

        index_type j2_f;

        if (z_f < zi[eMshSrc][j1_f]) {
            j2_f = (j1_f - 1 >= 0) ? j1_f - 1 : j1_f + 1;
        } else {
            j2_f = (j1_f + 1 < static_cast<ssize_t>(zi[eMshSrc].size())) ? j1_f + 1 : j1_f - 1;
        }

        HYDRO_ASSERT(j2_f >= 0);
        HYDRO_ASSERT(j2_f < static_cast<ssize_t>(zi[eMshSrc].size()));

        constexpr index_type GRIDINTERPX = 2;
        constexpr index_type GRIDINTERPZ = 2;

        const std::array<index_type, GRIDINTERPX> i_f = {i1_f, i2_f};
        const std::array<index_type, GRIDINTERPZ> j_f = {j1_f, j2_f};

        fij_f = 0.;

        T wij_f{0.};

        Vector<T, GRIDINTERP + 1> f;
        Vector<T, GRIDINTERP> Ck_f;

        for (index_type kx = 0; kx < GRIDINTERPX; ++kx) {
            const index_type ik_f = i_f[kx];

            for (index_type kz = 0; kz < GRIDINTERPZ; ++kz) {
                const index_type jk_f = j_f[kz];

                HYDRO_ASSERT(ik_f >= 0);
                HYDRO_ASSERT(ik_f < static_cast<ssize_t>(xi[eMshSrc].size()));
                HYDRO_ASSERT(jk_f >= 0);
                HYDRO_ASSERT(jk_f < static_cast<ssize_t>(zi[eMshSrc].size()));

                const index_type ijk_f = ik_f * zi[eMshSrc].size() + jk_f;

                HYDRO_ASSERT(ijk_f >= 0);
                HYDRO_ASSERT(ijk_f < static_cast<ssize_t>(rgMatIdx[eMshSrc].size()));

                const index_type dijk_f = rgMatIdx[eMshSrc][ijk_f];

                HYDRO_ASSERT(dijk_f >= 1);
                HYDRO_ASSERT(dijk_f <= w[eMshSrc].iGetNumRows());

                const PolyData& oPolyData = rgPolyData[eMshSrc][ijk_f];

                static_assert(eField == FT_TOTAL_PRESS || eField == FT_DEF_MOVING, "interpolation not supported");

                for (index_type l = 0; l <= GRIDINTERP; ++l) {
                    const index_type iCompIndex = (l < GRIDINTERP ? oPolyData.Cidx[l] : dijk_f);

                    switch (eField) {
                    case FT_TOTAL_PRESS: {
                        T p, pasp;

                        HYDRO_ASSERT(eMshSrc == DEHD_BODY_FIXED);
                        HYDRO_ASSERT(iCompIndex >= 1);
                        HYDRO_ASSERT(iCompIndex <= w[eMshSrc].iGetNumRows());
                        HYDRO_ASSERT(rgNodes[iCompIndex - 1]->iGetComplianceIndex() == iCompIndex);

                        pGetMesh()->GetPressure(rgNodes[iCompIndex - 1], p, pDofMap, dCoef);
                        rgNodes[iCompIndex - 1]->GetContactPressure(pasp, pDofMap);

                        f(l + 1) = (p + pasp) * dScale;
                    } break;
                    case FT_DEF_MOVING: {
                        T dw_dt;

                        HYDRO_ASSERT(eMshSrc == DEHD_BODY_MOVING);
                        HYDRO_ASSERT(iCompIndex >= 1);
                        HYDRO_ASSERT(iCompIndex <= w[eMshSrc].iGetNumRows());

                        GetRadialDeformation(f(l + 1),
                                             dw_dt,
                                             dCoef,
                                             func,
                                             pDofMap,
                                             DEHD_DEF_MOVING,
                                             iCompIndex);

                        f(l + 1) *= dScale;
                    } break;
                    }
                }

                for (index_type l = 1; l <= GRIDINTERP; ++l) {
                    Ck_f(l) = f(l) - f(GRIDINTERP + 1);
                }

                const Vector<T, POLYORDER> ak_f = oPolyData.pinvA * Ck_f;

                const T dxk_f = x_f - xi[eMshSrc][ik_f];
                const T dzk_f = z_f - zi[eMshSrc][jk_f];

                const T Ck_m = f(GRIDINTERP + 1)
                    + ak_f(1) * pow(dxk_f, 3)
                    + ak_f(2) * pow(dxk_f, 2) * dzk_f
                    + ak_f(3) * dxk_f * pow(dzk_f, 2)
                    + ak_f(4) * pow(dzk_f, 3)
                    + ak_f(5) * pow(dxk_f, 2)
                    + ak_f(6) * dxk_f * dzk_f
                    + ak_f(7) * pow(dzk_f, 2)
                    + ak_f(8) * dxk_f
                    + ak_f(9) * dzk_f;

                T dk_f_2 = dxk_f * dxk_f + dzk_f * dzk_f;

                if (dk_f_2 < dMinDistance_2) {
                    dk_f_2 = dMinDistance_2;
                }

                const T dk_f = sqrt(dk_f_2);
                const T wk_f = (dMeshRadius - dk_f) / (dMeshRadius * dk_f);
                const T wk_f3 = wk_f * wk_f * wk_f;

                fij_f += wk_f3 * Ck_m;
                wij_f += wk_f3;
            }
        }

        fij_f *= alpha / wij_f;
    }

    void ComplianceModelNodalDouble::UpdateDefMovingInterp(index_type iCompIndex, doublereal wmi, doublereal dScale)
    {
        HYDRO_ASSERT(iCompIndex >= 1);
        HYDRO_ASSERT(iCompIndex <= w[DEHD_DEF_MOVING_INTERP].iGetNumRows());

        w[DEHD_DEF_MOVING_INTERP](iCompIndex) = wmi / dScale;
    }

    template <typename T>
    void ComplianceModelNodalDouble::UnivAssRes(GradientAssVec<T>& WorkVec,
                                                doublereal dCoef,
                                                const GradientVectorHandler<T>& XCurr,
                                                enum FunctionCall func)
    {
#if HYDRO_DEBUG > 0
        for (index_type i = 0; i < DEHD_BODY_LAST; ++i) {
            HYDRO_ASSERT(rgNumNodes[i] > 0);
            HYDRO_ASSERT(rgNumNodes[i] == C[i].iGetNumRows() || C[i].iGetNumRows() == 0);
            HYDRO_ASSERT(rgNumNodes[i] == C[i].iGetNumCols() || C[i].iGetNumCols() == 0);
            HYDRO_ASSERT((rgNumNodes[i] == E[i].iGetNumCols() && rgModalJoints[i] != nullptr) || (E[i].iGetNumCols() == 0 && rgModalJoints[i] == nullptr));
            HYDRO_ASSERT((rgNumNodes[i] == D[i].iGetNumRows() && rgModalJoints[i] != nullptr) || (D[i].iGetNumRows() == 0 && rgModalJoints[i] == nullptr));
        }
#endif
        
        rgDofMap[DOF_GLOBAL2].Reset(func);

        const doublereal dEquationScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_ELASTICITY_EQU);
        integer iEqIndexHydro = iGetFirstIndex(func);

        Vector<T, 2> dxm;

        pGetMesh()->pGetGeometry()->GetMovingMeshOffset(dxm, &rgDofMap[DOF_GLOBAL1]);

        T dummy;

        Vector<T, DYNAMIC_SIZE> wm1(rgNumNodes[DEHD_BODY_MOVING]);

        for (index_type i = 1; i <= wm1.iGetNumRows(); ++i) {
            GetRadialDeformation(wm1(i), dummy, dCoef, func, &rgDofMap[DOF_GLOBAL1], DEHD_DEF_MOVING, i);
        }

        std::array<Vector<T, DYNAMIC_SIZE>, DEHD_BODY_LAST> a;

        for (index_type k = 0; k < DEHD_BODY_LAST; ++k) {
            if (rgModalJoints[k]) {
                HYDRO_ASSERT(rgModalJoints[k]->uGetNModes() == D[k].iGetNumCols());
                HYDRO_ASSERT((D[k].iGetNumRows() == C[k].iGetNumRows()) || C[k].iGetNumRows() == 0);
                HYDRO_ASSERT(D[k].iGetNumRows() == w[k].iGetNumRows());
                
                a[k].Resize(rgModalJoints[k]->uGetNModes());

                for (index_type j = 1; j <= D[k].iGetNumCols(); ++j) {
                    rgModalJoints[k]->GetACurr(j, a[k](j), dCoef, func, &rgDofMap[DOF_GLOBAL1]);
                }
            }
        }

        {
            Vector<T, DYNAMIC_SIZE> p(rgNumNodes[DEHD_BODY_FIXED]), pasp(rgNumNodes[DEHD_BODY_FIXED]);

            for (index_type j = 1; j <= p.iGetNumRows(); ++j) {
                // Do not add fluid and contact pressure for performance reasons (number of FLOPS caused by size of local slices)
                pGetMesh()->GetPressure(rgNodes[j - 1], p(j), &rgDofMap[DOF_GLOBAL1], dCoef);

                p(j) *= dPressScale;
            }

            for (index_type j = 1; j <= pasp.iGetNumRows(); ++j) {
                rgNodes[j - 1]->GetContactPressure(pasp(j), &rgDofMap[DOF_GLOBAL2]);

                pasp(j) *= dPressScale;
            }

            T f;
            
            for (index_type i = 1; i <= rgNumNodes[DEHD_BODY_FIXED]; ++i) {
                GetRadialDeformation(f, dummy, dCoef, func, &rgDofMap[DOF_GLOBAL1], rgNodes[i - 1]);

                f *= dEquationScale;

                for (index_type j = 1; j <= C[DEHD_BODY_FIXED].iGetNumCols(); ++j) {
                    const doublereal Cij = dEquationScale * C[DEHD_BODY_FIXED](i, j);
                    f -= Cij * p(j); // Do not use `Cij * (p(j) + pasp(j))` because it would be an O(N^2) operation
                    f -= Cij * pasp(j);
                }

                {
                    rgDofMap[DOF_LOCAL].Reset(); // Reduce size of local vectors slices during interpolation

                    const ComplianceMatrix::GridIndex& oGridIdx = rgActGridIdx[DEHD_BODY_FIXED][i - 1];

                    T wm2i;

                    Interpolate<FT_DEF_MOVING, DEHD_BODY_MOVING, DEHD_BODY_FIXED>(oGridIdx.ix,
                                                                                  oGridIdx.iz,
                                                                                  dxm,
                                                                                  dEquationScale,
                                                                                  wm2i,
                                                                                  dCoef,
                                                                                  func,
                                                                                  &rgDofMap[DOF_LOCAL]);

                    f -= wm2i;

                    UpdateDefMovingInterp(i, wm2i, dEquationScale);
                }

                if (rgModalJoints[DEHD_BODY_FIXED]) {
                    for (index_type j = 1; j <= D[DEHD_BODY_FIXED].iGetNumCols(); ++j) {
                        f -= dEquationScale * D[DEHD_BODY_FIXED](i, j) * a[DEHD_BODY_FIXED](j);
                    }
                }

                WorkVec.AddItem(iEqIndexHydro++, f);
            }

            if (rgModalJoints[DEHD_BODY_FIXED]) {
                index_type iEqIndexModal1 = rgModalJoints[DEHD_BODY_FIXED]->iGetFirstIndex() + rgModalJoints[DEHD_BODY_FIXED]->uGetNModes() + 1;

                for (index_type i = 1; i <= E[DEHD_BODY_FIXED].iGetNumRows(); ++i) {
                    f = 0.;

                    SetDofMap(f, &rgDofMap[DOF_GLOBAL1]);
                    
                    for (index_type j = 1; j <= E[DEHD_BODY_FIXED].iGetNumCols(); ++j) {
                        const doublereal Eij = E[DEHD_BODY_FIXED](i, j);

                        f += Eij * p(j);
                        f += Eij * pasp(j);
                    }

                    WorkVec.AddItem(iEqIndexModal1++, f);
                }
            }
        }

        {
            Vector<T, DYNAMIC_SIZE> fm1 = wm1 * dEquationScale;
            Vector<T, DYNAMIC_SIZE> fm2(E[DEHD_BODY_MOVING].iGetNumRows());

            for (index_type i = 1; i <= fm2.iGetNumRows(); ++i) {
                SetDofMap(fm2(i), &rgDofMap[DOF_GLOBAL1]);
            }

            for (index_type j = 1; j <= rgNumNodes[DEHD_BODY_MOVING]; ++j) {
                rgDofMap[DOF_LOCAL].Reset();

                const ComplianceMatrix::GridIndex& oGridIdx = rgActGridIdx[DEHD_BODY_MOVING][j - 1];
                T pmij;

                Interpolate<FT_TOTAL_PRESS, DEHD_BODY_FIXED, DEHD_BODY_MOVING>(oGridIdx.ix,
                                                                               oGridIdx.iz,
                                                                               dxm,
                                                                               dPressScale,
                                                                               pmij,
                                                                               dCoef,
                                                                               func,
                                                                               &rgDofMap[DOF_LOCAL]);

                for (index_type i = 1; i <= C[DEHD_BODY_MOVING].iGetNumRows(); ++i) {
                    const doublereal Cij = dEquationScale * C[DEHD_BODY_MOVING](i, j);
                    fm1(i) -= Cij * pmij;
                }

                if (rgModalJoints[DEHD_BODY_MOVING]) {
                    for (index_type i = 1; i <= E[DEHD_BODY_MOVING].iGetNumRows(); ++i) {
                        const doublereal Eij = E[DEHD_BODY_MOVING](i, j);
                        fm2(i) += Eij * pmij;
                    }
                }
            }

            if (rgModalJoints[DEHD_BODY_MOVING]) {
                for (index_type i = 1; i <= fm1.iGetNumRows(); ++i) {
                    for (index_type j = 1; j <= D[DEHD_BODY_MOVING].iGetNumCols(); ++j) {
                        fm1(i) -= dEquationScale * D[DEHD_BODY_MOVING](i, j) * a[DEHD_BODY_MOVING](j);
                    }
                }
            }

            for (index_type i = 1; i <= fm1.iGetNumRows(); ++i) {
                WorkVec.AddItem(iEqIndexHydro++, fm1(i));
            }

            if (rgModalJoints[DEHD_BODY_MOVING]) {
                index_type iEqIndexModal2 = rgModalJoints[DEHD_BODY_MOVING]->iGetFirstIndex()
                                          + rgModalJoints[DEHD_BODY_MOVING]->uGetNModes() + 1;

                for (index_type i = 1; i <= fm2.iGetNumRows(); ++i) {
                    WorkVec.AddItem(iEqIndexModal2++, fm2(i));
                }
            }
        }
    }

    void ComplianceModelNodalDouble::Print(std::ostream& os) const
    {
        HydroElement::Print(os);

        for (index_type k = 0; k < DEHD_BODY_LAST; ++k) {
            os << C[k].iGetNumRows() << '\n';

            for (index_type i = 1; i <= C[k].iGetNumRows(); ++i) {
                for (index_type j = 1; j <= C[k].iGetNumCols(); ++j) {
                    os << C[k](i, j) << '\t';
                }
                os << '\n';
            }
        }
    }

    ComplianceModelModal::ComplianceModelModal(HydroMesh* pMesh,
                                               doublereal dDefScale,
                                               doublereal dPressScale,
                                               const std::string& strFileName)
        :ComplianceModel(pMesh, dDefScale, dPressScale),
         iNumModes(0),
         strFileName(strFileName),
         eCurrFunc(grad::INITIAL_ASS_FLAG)
    {
    }

    ComplianceModelModal::~ComplianceModelModal()
    {
    }

    void ComplianceModelModal::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        const bool bDoAssembly = bDoInitAss || (eFunc & grad::REGULAR_FLAG);

        *piNumRows = iGetNumModes() * bDoAssembly;
        *piNumCols = (rgNodes.size() + iGetNumModes() + pGetMesh()->pGetGeometry()->iGetNumColsWorkSpace(eFunc)) * bDoAssembly;
    }

    void ComplianceModelModal::Initialize()
    {
        ComplianceModel::Initialize();

        ComplianceMatrix::MatrixData oMatData{ComplianceMatrix::LOC_MESH_FIXED,
                nullptr,
                &RPhiK,
                &Phin};

        ComplianceMatrixFileParser oParser(oMatData);

        oParser.Parse(strFileName, pGetMesh(), rgElements, rgNodes, dPressScale);

        iNumModes = RPhiK.iGetNumRows();

        q.Resize(iNumModes);
        dq_dt.Resize(iNumModes);

        HYDRO_ASSERT(iNumModes > 0);
        HYDRO_ASSERT(iNumModes < std::numeric_limits<index_type>::max());
        HYDRO_ASSERT(RPhiK.iGetNumRows() == iNumModes);
        HYDRO_ASSERT(Phin.iGetNumCols() == iNumModes);
        HYDRO_ASSERT(Phin.iGetNumRows() == iGetNumNodes());
        HYDRO_ASSERT(RPhiK.iGetNumCols() == iGetNumNodes());
    }

    void ComplianceModelModal::SetValue(VectorHandler& XCurr, VectorHandler& XPrimeCurr)
    {
        HYDRO_ASSERT(eCurrFunc == grad::INITIAL_ASS_FLAG);

        eCurrFunc = grad::REGULAR_FLAG;

        integer iDofIndex = iGetFirstIndex(eCurrFunc);
        const index_type iNumModes = iGetNumModes();

        for (index_type i = 1; i <= iNumModes; ++i, ++iDofIndex) {
            HYDRO_ASSERT(i <= XCurr.iGetSize());

            XCurr(iDofIndex) = q(i);
            XPrimeCurr(iDofIndex) = dq_dt(i);
        }
    }

    unsigned int ComplianceModelModal::iGetNumDof() const
    {
        return iGetNumModes();
    }

    unsigned int ComplianceModelModal::iGetInitialNumDof() const
    {
        return bDoInitAss ? iGetNumModes() : 0u;
    }

    integer ComplianceModelModal::iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type iNumNodes) const
    {
        return (eFunc & grad::REGULAR_FLAG) || bDoInitAss ? iGetNumModes() : 0;
    }

    DofOrder::Order ComplianceModelModal::GetDofType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::DIFFERENTIAL;
    }

    DofOrder::Order ComplianceModelModal::GetEqType(unsigned int i) const
    {
        HYDRO_ASSERT(i < iGetNumDof());

        return DofOrder::ALGEBRAIC;
    }

    std::ostream&
    ComplianceModelModal::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iDofIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);
            const index_type iNumModes = iGetNumModes();

            for (index_type i = 1; i <= iNumModes; ++i, ++iDofIndex) {
                out << prefix << iDofIndex << ": Modal elasticity dof for mode " << i << std::endl;
            }
        }

        return out;
    }

    std::ostream&
    ComplianceModelModal::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
    {
        if (bDoInitAss || !bInitial) {
            integer iDofIndex = iGetFirstIndex(bInitial ? grad::INITIAL_ASS_FLAG : grad::REGULAR_FLAG);
            const index_type iNumModes = iGetNumModes();

            for (index_type i = 1; i <= iNumModes; ++i, ++iDofIndex) {
                out << prefix << iDofIndex << ": Modal elasticity definition for mode " << i << std::endl;
            }
        }

        return out;
    }


    void ComplianceModelModal::Update(const VectorHandler& XCurr, const VectorHandler& XPrimeCurr, doublereal dCoef, FunctionCall func)
    {
        if (bDoInitAss || (func & grad::REGULAR_FLAG)) {
            integer iDofIndex = iGetFirstIndex(func);
            const index_type iNumModes = iGetNumModes();

            for (index_type i = 1; i <= iNumModes; ++i, ++iDofIndex) {
                HYDRO_ASSERT(iDofIndex >= 1);
                HYDRO_ASSERT(iDofIndex <= XCurr.iGetSize());

                q(i) = XCurr(iDofIndex);

                if (func & grad::REGULAR_FLAG) {
                    HYDRO_ASSERT(iDofIndex <= XPrimeCurr.iGetSize());
                    dq_dt(i) = XPrimeCurr(iDofIndex);
                } else {
                    HYDRO_ASSERT(dCoef == 1.);
                }
            }
        }
    }

    void
    ComplianceModelModal::AssRes(SubVectorHandler& WorkVec,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    template <typename T>
    void ComplianceModelModal::AssRes(GradientAssVec<T>& WorkVec,
                                      doublereal dCoef,
                                      const GradientVectorHandler<T>& XCurr,
                                      const GradientVectorHandler<T>& XPrimeCurr,
                                      enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename T>
    void ComplianceModelModal::InitialAssRes(GradientAssVec<T>& WorkVec,
                                             const GradientVectorHandler<T>& XCurr,
                                             enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    void
    ComplianceModelModal::AssJac(SparseSubMatrixHandler& WorkMat,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofModal,
                                             mode);
    }

    void
    ComplianceModelModal::InitialAssRes(SubVectorHandler& WorkVec,
                                        const VectorHandler& XCurr,
                                        GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<doublereal>::InitialAssRes(this,
                                                      WorkVec,
                                                      XCurr,
                                                      INITIAL_ASS_RES,
                                                      mode);
        }
    }

    void
    ComplianceModelModal::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                        const VectorHandler& XCurr,
                                        GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                        WorkMat,
                                                        XCurr,
                                                        INITIAL_ASS_JAC,
                                                        &oDofModal,
                                                        mode);
        }
    }

    template <typename T>
    void ComplianceModelModal::GetRadialDeformationTpl(T& wi, T& dwi_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap, const HydroUpdatedNode* pNode) const
    {
        const index_type iNodeIndex = pNode->iGetComplianceIndex();
        const index_type iNumModes = iGetNumModes();

        HYDRO_ASSERT(iNodeIndex >= 1);
        HYDRO_ASSERT(iNodeIndex <= iGetNumNodes());
        HYDRO_ASSERT(rgNodes[iNodeIndex - 1] == pNode);

        wi = 0.;
        dwi_dt = 0.;

        T qi;

        for (index_type iMode = 1; iMode <= iNumModes; ++iMode) {
            GetModalDeformation(iMode, qi, dCoef, func, pDofMap);

            wi += Phin(iNodeIndex, iMode) * qi;
        }

        if (func & grad::REGULAR_FLAG) {
            T dqi_dt;

            for (index_type iMode = 1; iMode <= iNumModes; ++iMode) {
                GetModalDeformationDer(iMode, dqi_dt, dCoef, func, pDofMap);

                dwi_dt += Phin(iNodeIndex, iMode) * dqi_dt;
            }
        } else {
            HYDRO_ASSERT(dCoef == 1.);
        }
    }

    void ComplianceModelModal::GetRadialDeformation(doublereal& wi, doublereal& dwi_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap, const HydroUpdatedNode* pNode) const
    {
        GetRadialDeformationTpl(wi, dwi_dt, dCoef, func, pDofMap, pNode);
    }

    void ComplianceModelModal::GetRadialDeformation(Gradient<0>& wi, Gradient<0>& dwi_dt, doublereal dCoef, FunctionCall func, LocalDofMap* pDofMap, const HydroUpdatedNode* pNode) const
    {
        GetRadialDeformationTpl(wi, dwi_dt, dCoef, func, pDofMap, pNode);
    }

    template <typename T>
    void ComplianceModelModal::UnivAssRes(GradientAssVec<T>& WorkVec,
                                          doublereal dCoef,
                                          const GradientVectorHandler<T>& XCurr,
                                          enum FunctionCall func)
    {
        const index_type iNumModes = iGetNumModes();
        const index_type iNumNodes = rgNodes.size();
        const doublereal dEquationScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_ELASTICITY_EQU);

        Vector<T, DYNAMIC_SIZE> f(iNumModes);

        for (index_type iMode = 1; iMode <= iNumModes; ++iMode) {
            GetModalDeformation(iMode, f(iMode), dCoef, func, &oDofModal);
        }

        T p, pasp;

        for (index_type jNode = 1; jNode <= iNumNodes; ++jNode) {
            oDofPressure.Reset();

            pGetMesh()->GetPressure(rgNodes[jNode - 1], p, &oDofPressure, dCoef);

            if (rgNodes[jNode - 1]->GetContactPressure(pasp, &oDofPressure)) {
                p += pasp;
            }

            p *= dPressScale;

            for (index_type iMode = 1; iMode <= iNumModes; ++iMode) {
                f(iMode) -= RPhiK(iMode, jNode) * p;
            }
        }

        integer iEqIndex = iGetFirstIndex(func);

        for (index_type iMode = 1; iMode <= iNumModes; ++iMode) {
            f(iMode) *= dEquationScale;

            CHECK_NUM_COLS_WORK_SPACE(this, func, f(iMode), iEqIndex);

            WorkVec.AddItem(iEqIndex++, f(iMode));
        }
    }

    void
    ComplianceModelModal::GetModalDeformation(index_type iMode,
                                              doublereal& qi,
                                              doublereal dCoef,
                                              FunctionCall func,
                                              LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(iMode >= 1);
        HYDRO_ASSERT(iMode <= iGetNumModes());

        qi = q(iMode);
    }

    void
    ComplianceModelModal::GetModalDeformation(index_type iMode,
                                              Gradient<0>& qi,
                                              doublereal dCoef,
                                              FunctionCall func,
                                              LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(iMode >= 1);
        HYDRO_ASSERT(iMode <= iGetNumModes());

        if (bDoInitAss || (func & grad::REGULAR_FLAG)) {
            const index_type iDofIndex = iGetFirstIndex(eCurrFunc) + iMode - 1;

            qi.SetValuePreserve(q(iMode));
            qi.DerivativeResizeReset(pDofMap,
                                     iDofIndex,
                                     MapVectorBase::GLOBAL,
                                     -dCoef);
        } else {
            qi.SetValue(q(iMode));
        }
    }

    void
    ComplianceModelModal::GetModalDeformationDer(index_type iMode,
                                                 doublereal& dqi_dt,
                                                 doublereal dCoef,
                                                 FunctionCall func,
                                                 LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(iMode >= 1);
        HYDRO_ASSERT(iMode <= iGetNumModes());

        dqi_dt = dq_dt(iMode);
    }

    void
    ComplianceModelModal::GetModalDeformationDer(index_type iMode,
                                                 Gradient<0>& dqi_dt,
                                                 doublereal dCoef,
                                                 FunctionCall func,
                                                 LocalDofMap* pDofMap) const
    {
        HYDRO_ASSERT(iMode >= 1);
        HYDRO_ASSERT(iMode <= iGetNumModes());

        if (func & grad::REGULAR_FLAG) {
            const index_type iDofIndex = iGetFirstIndex(eCurrFunc) + iMode - 1;

            dqi_dt.SetValuePreserve(dq_dt(iMode));
            dqi_dt.DerivativeResizeReset(pDofMap,
                                         iDofIndex,
                                         MapVectorBase::GLOBAL,
                                         -1.);
        } else {
            dqi_dt.SetValue(dq_dt(iMode));
        }
    }

    void ComplianceModelModal::Print(std::ostream& os) const
    {
        HydroElement::Print(os);

        os << RPhiK.iGetNumRows() << '\n';

        for (index_type i = 1; i <= RPhiK.iGetNumRows(); ++i) {
            for (index_type j = 1; j <= RPhiK.iGetNumCols(); ++j) {
                os << RPhiK(i, j) << '\t';
            }
            os << '\n';
        }

        os << '\n' << Phin.iGetNumRows() << '\n';

        for (index_type i = 1; i <= Phin.iGetNumRows(); ++i) {
            for (index_type j = 1; j <= Phin.iGetNumCols(); ++j) {
                os << Phin(i, j) << '\t';
            }
            os << '\n';
        }
    }

    void ComplianceMatrixFileParser::Parse(const std::string& strFileName,
                                           const HydroMesh* pMesh,
                                           const ElementContainer& rgElements,
                                           const NodesContainer& rgNodes,
                                           doublereal dPressScale)
    {
        std::ifstream oFile(strFileName.c_str(), std::ios::in);

        if (!oFile.good()) {
            silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                        << "): failed to open file \"" << strFileName << "\"\n");
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        std::vector<NodeRec> rgNodesHyd;

        rgNodesHyd.reserve(rgNodes.size());

        for (auto ppNode = rgNodes.begin(); ppNode != rgNodes.end(); ++ppNode) {
            const HydroUpdatedNode* const pNode = *ppNode;
            const Vector<doublereal, 2>& X = pNode->GetPosition2D();
            const index_type iCompIndex = pNode->iGetComplianceIndex();

            HYDRO_ASSERT(iCompIndex >= 1);
            HYDRO_ASSERT(iCompIndex <= static_cast<ssize_t>(rgNodes.size()));

            rgNodesHyd.emplace_back(iCompIndex,
                                    X(1),
                                    X(2));
        }

        std::sort(rgNodesHyd.begin(), rgNodesHyd.end());

        HYDRO_ASSERT(rgNodesHyd.size() == rgNodes.size());

        pedantic_cout("hydrodynamic plain bearing2("
                      << pMesh->pGetParent()->GetLabel()
                      << ") listing of updated nodes:\n");

        for (index_type i = 0; i < static_cast<ssize_t>(rgNodesHyd.size()); ++i) {
            pedantic_cout("\thydrodynamic node(" << i + 1 << "):" << rgNodesHyd[i].x << " " << rgNodesHyd[i].z << " " << rgNodesHyd[i].iComplianceIndex << "\n");
        }

        index_type iLineNo = 0;

        try {
            pedantic_cout("hydrodynamic plain bearing2("
                          << pMesh->pGetParent()->GetLabel()
                          << "): reading FEM file \"" << strFileName << "\"...\n");

            oFile.exceptions(std::ios::failbit | std::ios::eofbit | std::ios::badbit);

            std::vector<index_type> rgCompIndexHyd;
            std::string strTag;

            enum MatMask {
                MATMA_NODAL1 = 0x1,
                MATMA_NODAL2 = 0x2,
                MATMA_NODAL3 = 0x4,
                MATMA_MODAL1 = 0x8,
                MATMA_ALL = MATMA_NODAL1 | MATMA_NODAL2 | MATMA_NODAL3 | MATMA_MODAL1
            };

            unsigned uCurrMatMask = MATMA_ALL;
            
            enum TagType {
                FILE_FORMAT = 0,
                BEARING_DIAMETER,
                BEARING_WIDTH,
                GRID_X,
                GRID_Z,
                NODES,
                REF_PRESSURE,
                MATRIX_C1,
                MATRIX_C2,
                MATRIX_D2,
                MATRIX_E2,
                MATRIX_D3,
                MATRIX_E3,
                MATRIX_RPhiK,
                MATRIX_Phin,
                LAST_TAG
            };

            static const struct {
                TagType eTag;
                MatMask eMatMask;
                char szName[32];
            } rgTags[] = {
                {FILE_FORMAT,      MATMA_ALL,    "file format"},
                {BEARING_DIAMETER, MATMA_ALL,    "bearing diameter"},
                {BEARING_WIDTH,    MATMA_ALL,    "bearing width"},
                {GRID_X,           MATMA_ALL,    "circumferential grid"},
                {GRID_Z,           MATMA_ALL,    "axial grid"},
                {NODES,            MATMA_ALL,    "nodes"},
                {REF_PRESSURE,     MATMA_ALL,    "reference pressure"},
                {MATRIX_C1,        MATMA_NODAL1, "compliance matrix"},
                {MATRIX_C2,        MATMA_NODAL2, "compliance matrix substruct"},
                {MATRIX_D2,        MATMA_NODAL2, "substruct contrib matrix"},
                {MATRIX_E2,        MATMA_NODAL2, "substruct residual matrix"},
                {MATRIX_D3,        MATMA_NODAL3, "substruct total contrib matrix"},
                {MATRIX_E3,        MATMA_NODAL3, "substruct total residual matrix"},
                {MATRIX_RPhiK,     MATMA_MODAL1, "modal load"},
                {MATRIX_Phin,      MATMA_MODAL1, "mode shapes"}
            };

            enum FileFormatType {
                FILE_FORMAT_1_0 = 10,
                FILE_FORMAT_1_1 = 11,
                FILE_FORMAT_1_2 = 12
            };

            unsigned uFileFormat = FILE_FORMAT_1_0;

            constexpr index_type iNumTags = sizeof(rgTags) / sizeof(rgTags[0]);
            static_assert(iNumTags == LAST_TAG, "number of tags does not match");
            std::array<bool, iNumTags> rgTagsParsed = {false};
            bool bAllTagsParsed = false;
            index_type iNumModes = -1;
            const doublereal dMeshSideDiameter = 2. * pMesh->pGetGeometry()->dGetMeshRadius();

            ++iLineNo;

            while (!bAllTagsParsed) {
                int ch = oFile.get();
                switch (ch) {
                case ':':
                    for (index_type iTag = 0; iTag < iNumTags; ++iTag) {
                        if (strTag == rgTags[iTag].szName) {
                            pedantic_cout("\"" << strFileName << "\":"
                                          << iLineNo << ":" << strTag << "\n");

                            if ((uCurrMatMask & rgTags[iTag].eMatMask) == 0) {
                                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                            }
                            
                            uCurrMatMask = uCurrMatMask & rgTags[iTag].eMatMask;
                            
                            switch (rgTags[iTag].eTag) {
                            case FILE_FORMAT:
                                oFile >> uFileFormat;
                                pedantic_cout("\tFEM file format=" << uFileFormat << "\n");

                                if (uFileFormat != FILE_FORMAT_1_2) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): invalid file format: expected " << FILE_FORMAT_1_1
                                                << " but got " << uFileFormat << std::endl);
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                rgTagsParsed[FILE_FORMAT] = true;
                                break;

                            case BEARING_DIAMETER:
                                oFile >> oMatData.dBearingDiameter;
                                pedantic_cout("\tFEM bearing diameter=" << oMatData.dBearingDiameter << "\n");
                                rgTagsParsed[BEARING_DIAMETER] = true;
                                break;

                            case BEARING_WIDTH:
                                oFile >> oMatData.dBearingWidth;
                                pedantic_cout("\tFEM bearing width=" << oMatData.dBearingWidth << "\n");
                                rgTagsParsed[BEARING_WIDTH] = true;
                                break;

                            case GRID_X: {
                                index_type iSizeGridX;
                                index_type iGridX;
                                doublereal dGridX;

                                if (!rgTagsParsed[BEARING_DIAMETER]) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): missing record \"bearing diameter\"\n");
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oFile >> iSizeGridX;
                                pedantic_cout("\tFEM grid x: " << iSizeGridX << "\n");

                                if (iSizeGridX < 2) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oMatData.rgGridX.clear();
                                oMatData.rgGridX.reserve(iSizeGridX);

                                while (iSizeGridX--) {
                                    oFile >> iGridX >> dGridX;
                                    dGridX *= dMeshSideDiameter / oMatData.dBearingDiameter;

                                    pedantic_cout(iGridX << ":\t" << dGridX << "\n");

                                    if (oMatData.rgGridX.size() && dGridX <= oMatData.rgGridX.back()) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    oMatData.rgGridX.push_back(dGridX);
                                    ++iLineNo;
                                };

                                rgTagsParsed[GRID_X] = true;
                            } break;

                            case GRID_Z: {
                                index_type iSizeGridZ;
                                doublereal dGridZ;
                                index_type iGridZ;

                                if (!rgTagsParsed[BEARING_WIDTH]) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): missing record \"bearing width\"\n");
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oFile >> iSizeGridZ;
                                pedantic_cout("\tFEM grid z: " << iSizeGridZ << "\n");

                                if (iSizeGridZ < 2) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oMatData.rgGridZ.clear();
                                oMatData.rgGridZ.reserve(iSizeGridZ);

                                while (iSizeGridZ--) {
                                    oFile >> iGridZ >> dGridZ;

                                    pedantic_cout(iGridZ << ":\t" << dGridZ << "\n");

                                    if (oMatData.rgGridZ.size() && dGridZ <= oMatData.rgGridZ.back()) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    oMatData.rgGridZ.push_back(dGridZ);
                                    ++iLineNo;
                                }

                                rgTagsParsed[GRID_Z] = true;
                            } break;

                            case NODES: {
                                index_type iNumNodeFE, iNumDofFE;
                                InpNodeRec oNodeFE;

                                if (!rgTagsParsed[GRID_X]) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): missing record \"circumferential grid\"\n");
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                if (!rgTagsParsed[GRID_Z]) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): missing record \"axial grid\"\n");
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                if (!rgTagsParsed[FILE_FORMAT]) {
                                    silent_cerr("hydrodynamic plain bearing2("
                                                << pMesh->pGetParent()->GetLabel()
                                                << "): missing record \"file format\"\n");
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oFile >> iNumNodeFE >> iNumDofFE;

                                pedantic_cout("\tnumber of FEM nodes: " << iNumNodeFE << "\n");

                                if (oMatData.eMeshLocation == LOC_MESH_FIXED && iNumDofFE != static_cast<ssize_t>(rgNodesHyd.size())) {
                                    silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                                                << " number of FEM DOF in file \"" << strFileName << "\" = " << iNumDofFE
                                                << " does not match  number of hydrodynamic nodes = "
                                                << rgNodesHyd.size() << std::endl);
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                rgCompIndexHyd.clear();
                                rgCompIndexHyd.resize(iNumDofFE, -1);
                                oMatData.rgMatIdx.clear();
                                oMatData.rgMatIdx.resize(oMatData.rgGridX.size() * oMatData.rgGridZ.size(), -1);

                                oMatData.rgActGridIdx.clear();
                                oMatData.rgActGridIdx.resize(iNumDofFE, GridIndex(-1, -1));

                                const doublereal dTolNodePos = sqrt(std::numeric_limits<doublereal>::epsilon())
                                    * std::max(oMatData.dBearingWidth, M_PI * oMatData.dBearingDiameter);

                                for (index_type i = 0; i < iNumNodeFE; ++i) {
                                    oFile >> oNodeFE.iNode;
                                    oFile >> oNodeFE.iDof;
                                    oFile >> oNodeFE.iPosX;
                                    oFile >> oNodeFE.iPosZ;

                                    if (oNodeFE.iNode != i + 1) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    if (oNodeFE.iDof < 1 || oNodeFE.iDof > iNumDofFE) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    if (oNodeFE.iPosX < 1 || oNodeFE.iPosX > static_cast<ssize_t>(oMatData.rgGridX.size())) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    if (oNodeFE.iPosZ < 1 || oNodeFE.iPosZ > static_cast<ssize_t>(oMatData.rgGridZ.size())) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    ++iLineNo;

                                    if (oMatData.eMeshLocation == LOC_MESH_FIXED) {
                                        bool bNodeFound = false;

                                        for (index_type iNodeHyd = 0; iNodeHyd < static_cast<ssize_t>(rgNodesHyd.size()); ++iNodeHyd) {
                                            doublereal dx = fabs(oMatData.rgGridX[oNodeFE.iPosX - 1] - rgNodesHyd[iNodeHyd].x);
                                            doublereal dz = oMatData.rgGridZ[oNodeFE.iPosZ - 1] - rgNodesHyd[iNodeHyd].z;

                                            if (dx > 0.5 * dMeshSideDiameter * M_PI) {
                                                dx -= dMeshSideDiameter * M_PI;
                                            }

                                            if (fabs(dx) <= dTolNodePos && fabs(dz) <= dTolNodePos) {
                                                const index_type iCompIndex = rgNodesHyd[iNodeHyd].iComplianceIndex;

                                                HYDRO_ASSERT(iCompIndex >= 1);
                                                HYDRO_ASSERT(iCompIndex <= iNumDofFE);

                                                if (rgCompIndexHyd[oNodeFE.iDof - 1] == -1) {
                                                    rgCompIndexHyd[oNodeFE.iDof - 1] = iCompIndex;
                                                } else if (rgCompIndexHyd[oNodeFE.iDof - 1] != iCompIndex) {
                                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                                }

                                                const size_t j = (oNodeFE.iPosX - 1) * oMatData.rgGridZ.size() + oNodeFE.iPosZ - 1;

                                                if (oMatData.rgMatIdx[j] != -1) {
                                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                                }

                                                oMatData.rgMatIdx[j] = iCompIndex;

                                                GridIndex& oGridIdx = oMatData.rgActGridIdx[iCompIndex - 1];
                                                oGridIdx.ix = oNodeFE.iPosX - 1;
                                                oGridIdx.iz = oNodeFE.iPosZ - 1;

                                                bNodeFound = true;
                                                break;
                                            }
                                        }

                                        pedantic_cout("FEM node: " << oNodeFE.iNode << " "
                                                      << oMatData.rgGridX[oNodeFE.iPosX - 1] << " "
                                                      << oMatData.rgGridZ[oNodeFE.iPosZ - 1] << " " << rgCompIndexHyd[i] + 1 << "\n");

                                        if (!bNodeFound) {
                                            silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                                                        << "): FEM node " << i + 1 << " does not match any hydrodynamic node\n");
                                            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                        }
                                    } else {
                                        const index_type iCompIndex = oNodeFE.iDof;

                                        HYDRO_ASSERT(iCompIndex >= 1);
                                        HYDRO_ASSERT(iCompIndex <= iNumDofFE);

                                        if (rgCompIndexHyd[oNodeFE.iDof - 1] == -1) {
                                            rgCompIndexHyd[oNodeFE.iDof - 1] = iCompIndex;
                                        } else if (rgCompIndexHyd[oNodeFE.iDof - 1] != iCompIndex) {
                                            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                        }

                                        const size_t j = (oNodeFE.iPosX - 1) * oMatData.rgGridZ.size() + oNodeFE.iPosZ - 1;

                                        if (oMatData.rgMatIdx[j] != -1) {
                                            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                        }

                                        oMatData.rgMatIdx[j] = iCompIndex;

                                        GridIndex& oGridIdx = oMatData.rgActGridIdx[iCompIndex - 1];
                                        oGridIdx.ix = oNodeFE.iPosX - 1;
                                        oGridIdx.iz = oNodeFE.iPosZ - 1;
                                    }
                                }

                                auto i = std::find(rgCompIndexHyd.begin(), rgCompIndexHyd.end(), -1);

                                if (i != rgCompIndexHyd.end()) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                i = std::find(oMatData.rgMatIdx.begin(), oMatData.rgMatIdx.end(), -1);

                                if (i != oMatData.rgMatIdx.end()) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                auto j = std::find(oMatData.rgActGridIdx.begin(), oMatData.rgActGridIdx.end(), GridIndex(-1, -1));

                                if (j != oMatData.rgActGridIdx.end()) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }
                                
                                rgTagsParsed[NODES] = true;
                            } break;

                            case REF_PRESSURE:
                                oFile >> oMatData.dRefPressure;
                                pedantic_cout("\tFEM reference pressure=" << oMatData.dRefPressure << "\n");
                                rgTagsParsed[REF_PRESSURE] = true;
                                break;

                            case MATRIX_C1:
                            case MATRIX_C2:
                            case MATRIX_D2:
                            case MATRIX_E2:
                            case MATRIX_D3:
                            case MATRIX_E3:
                            case MATRIX_RPhiK:
                            case MATRIX_Phin: {
                                const TagType eCurrMatrix = rgTags[iTag].eTag;
                                MatrixKind eFileTypeCurr = MAT_INVALID;

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                case MATRIX_C2:
                                case MATRIX_D2:
                                case MATRIX_E2:
                                case MATRIX_D3:
                                case MATRIX_E3:
                                    eFileTypeCurr = MAT_FULL;
                                    break;

                                case MATRIX_RPhiK:
                                case MATRIX_Phin:
                                    eFileTypeCurr = MAT_MODAL;
                                    break;

                                default:
                                    HYDRO_ASSERT(0);
                                }

                                if (oMatData.rgMatrices.GetMatrixType() != eFileTypeCurr) {
                                    silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                                                << " invalid matrix type for selected compliance model in file \""
                                                << strFileName << "\" at line " << iLineNo << std::endl);
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                    if (oMatData.pModalJoint) {
                                        silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                                                    << " keyword \"modal element\" must not be used for matrix type \""
                                                    << rgTags[iTag].szName << "\" in file \""
                                                    << strFileName << "\" at line " << iLineNo << std::endl);
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                case MATRIX_C2:
                                case MATRIX_D2:
                                case MATRIX_E2:
                                case MATRIX_D3:
                                case MATRIX_E3:                                    
                                    if (!oMatData.pModalJoint) {
                                        silent_cerr("hydrodynamic plain bearing2(" << pMesh->pGetParent()->GetLabel()
                                                    << " keyword \"modal element\" is required for matrix type \""
                                                    << rgTags[iTag].szName << "\" in file \""
                                                    << strFileName << "\" at line " << iLineNo << std::endl);
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                default:
                                    break;
                                }

                                index_type iNumRows = -1, iNumCols = -1;

                                if (!rgTagsParsed[NODES] || !rgTagsParsed[REF_PRESSURE]) {
                                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                }

                                oFile >> iNumRows >> iNumCols;

                                pedantic_cout("\t" << rgTags[iTag].szName << " size=" << iNumRows << "x" << iNumCols << "\n");

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                case MATRIX_C2:
                                case MATRIX_D2:
                                case MATRIX_D3:
                                case MATRIX_Phin:
                                    if (iNumRows != static_cast<ssize_t>(rgCompIndexHyd.size())) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                case MATRIX_E2:
                                case MATRIX_E3:
                                    if (!oMatData.pModalJoint || iNumRows != oMatData.pModalJoint->uGetNModes()) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                case MATRIX_RPhiK:
                                    if (iNumModes < 0) {
                                        iNumModes = iNumRows;
                                    } else if (iNumModes != iNumRows) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                default:
                                    HYDRO_ASSERT(0);
                                }

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                case MATRIX_C2:
                                case MATRIX_E2:
                                case MATRIX_E3:
                                case MATRIX_RPhiK:
                                    if (iNumCols != static_cast<ssize_t>(rgCompIndexHyd.size())) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;
                                case MATRIX_D2:
                                case MATRIX_D3:
                                    if (!oMatData.pModalJoint || iNumCols != oMatData.pModalJoint->uGetNModes()) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;
                                case MATRIX_Phin:
                                    if (iNumModes < 0) {
                                        iNumModes = iNumCols;
                                    } else if (iNumModes != iNumCols) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }
                                    break;

                                default:
                                    HYDRO_ASSERT(0);
                                }

                                doublereal dVal, dScale = 0.;

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                case MATRIX_C2:
                                case MATRIX_E2:
                                case MATRIX_E3:
                                case MATRIX_RPhiK:
                                    dScale = 1. / (dPressScale * oMatData.dRefPressure);
                                    break;

                                case MATRIX_D2:
                                case MATRIX_D3:
                                case MATRIX_Phin:
                                    dScale = 1.; // we assume that Phi has already been scaled by L2 norm
                                    break;

                                default:
                                    HYDRO_ASSERT(0);
                                };

                                MatrixType* pCurrMat = nullptr;

                                switch (eCurrMatrix) {
                                case MATRIX_C1:
                                case MATRIX_C2:
                                    pCurrMat = oMatData.rgMatrices.C();
                                    pCurrMat->Resize(rgCompIndexHyd.size(), rgCompIndexHyd.size());

                                    if (oMatData.eMeshLocation != LOC_MESH_MOVING) {
                                        HYDRO_ASSERT(static_cast<size_t>(pCurrMat->iGetNumRows()) == rgNodes.size());
                                        HYDRO_ASSERT(static_cast<size_t>(pCurrMat->iGetNumCols()) == rgNodes.size());
                                        HYDRO_ASSERT(rgNodes.size() == rgCompIndexHyd.size());
                                    }
                                    break;

                                case MATRIX_D2:
                                case MATRIX_D3:
                                    HYDRO_ASSERT(oMatData.pModalJoint != nullptr);

                                    pCurrMat = oMatData.rgMatrices.D();
                                    pCurrMat->Resize(rgCompIndexHyd.size(), oMatData.pModalJoint->uGetNModes());
                                    
                                    if (oMatData.eMeshLocation != LOC_MESH_MOVING) {
                                        HYDRO_ASSERT(static_cast<size_t>(pCurrMat->iGetNumRows()) == rgNodes.size());
                                        HYDRO_ASSERT(pCurrMat->iGetNumCols() == oMatData.pModalJoint->uGetNModes());
                                        HYDRO_ASSERT(rgNodes.size() == rgCompIndexHyd.size());
                                    }
                                    break;

                                case MATRIX_E2:
                                case MATRIX_E3:
                                    HYDRO_ASSERT(oMatData.pModalJoint != nullptr);

                                    pCurrMat = oMatData.rgMatrices.E();
                                    pCurrMat->Resize(oMatData.pModalJoint->uGetNModes(), rgCompIndexHyd.size());
                                    
                                    if (oMatData.eMeshLocation != LOC_MESH_MOVING) {
                                        HYDRO_ASSERT(static_cast<size_t>(pCurrMat->iGetNumCols()) == rgNodes.size());
                                        HYDRO_ASSERT(pCurrMat->iGetNumRows() == oMatData.pModalJoint->uGetNModes());
                                        HYDRO_ASSERT(rgNodes.size() == rgCompIndexHyd.size());
                                    }
                                    break;

                                case MATRIX_RPhiK:
                                    if (oMatData.eMeshLocation != LOC_MESH_FIXED) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    pCurrMat = oMatData.rgMatrices.RPhiK();
                                    pCurrMat->Resize(iNumModes, rgNodes.size());
                                    break;

                                case MATRIX_Phin:
                                    if (oMatData.eMeshLocation != LOC_MESH_FIXED) {
                                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                                    }

                                    pCurrMat = oMatData.rgMatrices.Phin();
                                    pCurrMat->Resize(rgNodes.size(), iNumModes);
                                    break;

                                default:
                                    HYDRO_ASSERT(0);
                                }

                                for (index_type i = 0; i < iNumRows; ++i) {
                                    index_type iRowIndexMat = -1;

                                    switch (eCurrMatrix) {
                                    case MATRIX_C1:
                                    case MATRIX_C2:
                                    case MATRIX_D2:
                                    case MATRIX_D3:
                                    case MATRIX_Phin:
                                        iRowIndexMat = rgCompIndexHyd[i];
                                        break;

                                    case MATRIX_RPhiK:
                                    case MATRIX_E2:
                                    case MATRIX_E3:
                                        HYDRO_ASSERT(i < pCurrMat->iGetNumRows());
                                        iRowIndexMat = i + 1;
                                        break;

                                    default:
                                        HYDRO_ASSERT(0);
                                    }

                                    pedantic_cout('\t');

                                    for (index_type j = 0; j < iNumCols; ++j) {
                                        index_type iColIndexMat = -1;

                                        switch (eCurrMatrix) {
                                        case MATRIX_C1:
                                        case MATRIX_C2:
                                        case MATRIX_E2:
                                        case MATRIX_E3:
                                        case MATRIX_RPhiK:
                                            iColIndexMat = rgCompIndexHyd[j];
                                            break;

                                        case MATRIX_D2:
                                        case MATRIX_D3:
                                        case MATRIX_Phin:
                                            HYDRO_ASSERT(j < pCurrMat->iGetNumCols());
                                            iColIndexMat = j + 1;
                                            break;

                                        default:
                                            HYDRO_ASSERT(0);
                                        }

                                        oFile >> dVal;
                                        pedantic_cout(dVal << '\t');

                                        HYDRO_ASSERT(iRowIndexMat >= 1);
                                        HYDRO_ASSERT(iColIndexMat >= 1);
                                        HYDRO_ASSERT(iRowIndexMat <= pCurrMat->iGetNumRows());
                                        HYDRO_ASSERT(iColIndexMat <= pCurrMat->iGetNumCols());

                                        (*pCurrMat)(iRowIndexMat, iColIndexMat) += dVal * dScale;
                                    }

                                    pedantic_cout('\n');

                                    ++iLineNo;
                                }

                                for (index_type i = 0; i < iNumTags; ++i) {
                                    if ((rgTags[i].eMatMask & rgTags[iTag].eMatMask) == 0) {
                                        rgTagsParsed[i] = true;
                                    }
                                }
                                
                                rgTagsParsed[eCurrMatrix] = true;
                            } break;

                            default:
                                HYDRO_ASSERT(false);
                            };
                            break;
                        }
                    }

                    strTag.clear();

                    bAllTagsParsed = true;

                    for (index_type i = 0; i < iNumTags; ++i) {
                        if (!rgTagsParsed[i]) {
                            bAllTagsParsed = false;;
                        }
                    }
                    break;

                case '\n':
                    strTag.clear();
                    ++iLineNo;
                    break;

                default:
                    strTag += ch;
                }
            }

            HYDRO_ASSERT(bAllTagsParsed);
        } catch (const std::ios_base::failure& oErr) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pMesh->pGetParent()->GetLabel()
                        << "): IO error in file \""
                        << strFileName << "\" at line "
                        << iLineNo << "\n");

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        } catch (const std::exception& oErr) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pMesh->pGetParent()->GetLabel()
                        << "): syntax error in file \""
                        << strFileName << "\" at line "
                        << iLineNo << "\n");
            throw;
        }
    }

    ComplianceMatrix::ComplianceMatrix()
    {
    }

    ComplianceMatrix::~ComplianceMatrix()
    {
    }

    ElasticHalfSpace::ElasticHalfSpace(const Material& oMaterial1,
                                       const Material& oMaterial2)
        :Ered(oMaterial1.dGetReducedModulus(oMaterial2)) {
    }

    ElasticHalfSpace::~ElasticHalfSpace()
    {
    }

    void ElasticHalfSpace::AddCompliance(MatrixData& oMatData,
                                         const HydroMesh* const pMesh,
                                         const ElementContainer& rgElements,
                                         const NodesContainer& rgNodes,
                                         doublereal dPressScale) const
    {
        Matrix<doublereal, DYNAMIC_SIZE, DYNAMIC_SIZE>& C = *oMatData.rgMatrices.C();
        
        const doublereal alpha = 0.25 / (M_PI * Ered * dPressScale);

        typedef std::multimap<const HydroNode*, const PressureElement*> NodeToElemCont;

        NodeToElemCont mmNodeToElem;

        for (auto pElem = rgElements.begin(); pElem != rgElements.end(); ++pElem) {
            const integer iNumNodes = pElem->iGetNumNodes();

            for (integer i = 0; i < iNumNodes; ++i) {
                const HydroNode* pNode = pElem->pGetNode(i);
                mmNodeToElem.insert(std::make_pair(pNode, &*pElem));
            }
        }

        Vector<Vector<doublereal, 2>, 4> X;
        Vector<doublereal, 2> dX;
        const BearingGeometry* const pGeometry = pMesh->pGetGeometry();
        const index_type iNumNodes = rgNodes.size();

        std::set<doublereal> sx, sz;

        const HydroRootElement* const pRootElem = pMesh->pGetParent();
        
        for (index_type iNode = 0; iNode < pRootElem->iGetNumNodes(); ++iNode) {
            const HydroNode* const pNode = dynamic_cast<const HydroNode*>(pRootElem->pGetNode(iNode));

            if (!pNode) {
                continue;
            }

            auto x = pNode->GetPosition2D();
                
            sx.insert(x(1));
            sz.insert(x(2));
        }

        oMatData.rgGridX.clear();
        oMatData.rgGridZ.clear();
        oMatData.rgGridX.reserve(sx.size());
        oMatData.rgGridZ.reserve(sz.size());

        std::copy(sx.begin(), sx.end(), std::back_inserter(oMatData.rgGridX));
        std::copy(sz.begin(), sz.end(), std::back_inserter(oMatData.rgGridZ));
        
        oMatData.rgActGridIdx.clear();
        oMatData.rgActGridIdx.resize(iNumNodes, GridIndex(-1, -1));
        oMatData.rgMatIdx.clear();
        oMatData.rgMatIdx.resize(oMatData.rgGridX.size() * oMatData.rgGridZ.size(), -1);

        for (index_type iNode = 0; iNode < pRootElem->iGetNumNodes(); ++iNode) {
            const HydroNode* const pNode = dynamic_cast<const HydroNode*>(pRootElem->pGetNode(iNode));

            if (!pNode) {
                continue;
            }
            
            auto x = pNode->GetPosition2D();
            
            const auto ix = std::lower_bound(oMatData.rgGridX.begin(), oMatData.rgGridX.end(), x(1));

            if (ix == oMatData.rgGridX.end()) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const auto iPosX = ix - oMatData.rgGridX.begin();
            
            const auto iz = std::lower_bound(oMatData.rgGridZ.begin(), oMatData.rgGridZ.end(), x(2));

            if (iz == oMatData.rgGridZ.end()) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const auto iPosZ = iz - oMatData.rgGridZ.begin();
            
            const index_type iCompIndex = pNode->iGetComplianceIndex();

            const auto iMatIdx = iPosX * oMatData.rgGridZ.size() + iPosZ;

            if (oMatData.rgMatIdx[iMatIdx] != -1) {
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            oMatData.rgMatIdx[iMatIdx] = iCompIndex;

            auto& oGridIdx = oMatData.rgActGridIdx[iCompIndex - 1];
            
            oGridIdx.ix = iPosX;
            oGridIdx.iz = iPosZ;
        }

        auto iInvalid = std::find(oMatData.rgMatIdx.begin(), oMatData.rgMatIdx.end(), -1);

        if (iInvalid != oMatData.rgMatIdx.end()) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        
        if (C.iGetNumRows() == 0 && C.iGetNumCols() == 0) {
            C.Resize(iNumNodes, iNumNodes);
        } else if (C.iGetNumRows() != iNumNodes || C.iGetNumCols() != iNumNodes) {
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
        
        for (index_type iColIndex = 1; iColIndex <= iNumNodes; ++iColIndex) {
            const HydroUpdatedNode* const pNodeCol = rgNodes[iColIndex - 1];
            for (index_type iRowIndex = 1; iRowIndex <= iNumNodes; ++iRowIndex) {
                const HydroUpdatedNode* const pNodeRow = rgNodes[iRowIndex - 1];

                auto ppElemPress = mmNodeToElem.lower_bound(pNodeCol);
                const auto ppLastElemPress = mmNodeToElem.upper_bound(pNodeCol);

                for ( ; ppElemPress != ppLastElemPress; ++ppElemPress) {
                    const auto pElemPress = ppElemPress->second;

                    for (integer iNodeElemPress = 0; iNodeElemPress < pElemPress->iGetNumNodes(); ++iNodeElemPress) {
                        const HydroNode* const pNodePress = pElemPress->pGetNode(iNodeElemPress);
                        X(iNodeElemPress + 1) = pNodePress->GetPosition2D();
                    }

                    const doublereal dxElem = X(PressureElement::NODE_SE + 1)(1) - X(PressureElement::NODE_SW + 1)(1);
                    const doublereal dzElem = X(PressureElement::NODE_NW + 1)(2) - X(PressureElement::NODE_SW + 1)(2);
                    const doublereal dA = dxElem * dzElem;
                    const Vector<doublereal, 2> XElemCenter = (X(PressureElement::NODE_NE + 1) + X(PressureElement::NODE_SW + 1)) * 0.5;
                    pGeometry->GetClosestDistance2D(XElemCenter, pNodeRow->GetPosition2D(), dX);

                    C(iRowIndex, iColIndex) += alpha * dA / Norm(dX);
                }
            }
        }
    }

    ComplianceFromFile::ComplianceFromFile(const std::string& strFileName)
        :strFileName(strFileName)
    {
    }

    ComplianceFromFile::~ComplianceFromFile()
    {
    }

    void ComplianceFromFile::AddCompliance(MatrixData& oMatData,
                                           const HydroMesh* pMesh,
                                           const ElementContainer& rgElements,
                                           const NodesContainer& rgNodes,
                                           doublereal dPressScale) const
    {
        ComplianceMatrixFileParser oParser{oMatData};

        oParser.Parse(strFileName, pMesh, rgElements, rgNodes, dPressScale);
    }

    BearingGeometry::BearingGeometry(HydroRootElement* pParent)
        :pParent(pParent)
    {
#if CREATE_PROFILE == 1
        memset(&profile, 0, sizeof(profile));
#endif
    }

    BearingGeometry::~BearingGeometry()
    {
#if CREATE_PROFILE == 1
        for (int i = PROF_RES; i <= PROF_JAC; ++i) {
            std::cerr << "dtAddForce[" << i << "]=" << profile.dtAddForce[i] << std::endl;
            std::cerr << "dtOperatorPlus[" << i << "]=" << profile.dtOperatorPlus[i] << std::endl;
        }

#endif
    }

    void
    BearingGeometry::GetNonNegativeClearance(const doublereal& h,
                                             doublereal& hn,
                                             const doublereal* dh_dt,
                                             doublereal* dhn_dt) const
    {
        GetNonNegativeClearanceTpl(h, hn, dh_dt, dhn_dt);
    }

    void
    BearingGeometry::GetNonNegativeClearance(const Gradient<0>& h,
                                             Gradient<0>& hn,
                                             const Gradient<0>* dh_dt,
                                             Gradient<0>* dhn_dt) const
    {
        GetNonNegativeClearanceTpl(h, hn, dh_dt, dhn_dt);
    }

    template <typename T> inline void
    BearingGeometry::GetNonNegativeClearanceTpl(const T& h,
                                                T& hn,
                                                const T* dh_dt,
                                                T* dhn_dt) const
    {
#if HYDRO_DEBUG > 0
        if (dh_dt != 0) {
            HYDRO_ASSERT(dhn_dt != 0);
        } else {
            HYDRO_ASSERT(dhn_dt == 0);
        }
#endif
        const doublereal h0 = dGetMinClearance();

        const T H = h / h0;

        if (H >= 2) {
            hn = h;

            if (dh_dt != 0) {
                *dhn_dt = *dh_dt;
            }
        } else {
            const T alpha = exp(H - 2);

            hn = h0 * (1. + alpha);

            if (dh_dt != 0) {
                *dhn_dt = *dh_dt;
                *dhn_dt *= alpha;
            }
        }
    }

    void BearingGeometry::AddMovingLubrGroove(LubricationGroove* pGroove)
    {
        HYDRO_ASSERT(rgMovingGrooves.size() < rgMovingGrooves.capacity());

        rgMovingGrooves.push_back(pGroove);
    }

    void BearingGeometry::ReserveMovingLubrGrooves(size_t n)
    {
        rgMovingGrooves.reserve(n);
    }

    const LubricationGroove* BearingGeometry::pFindMovingLubricationGroove(const Vector<doublereal, 2>& x, Node2D::NodeType eNodeType) const
    {
        for (auto i = rgMovingGrooves.begin(); i != rgMovingGrooves.end(); ++i) {
            if ((*i)->pGetGeometry()->bPointIsInside(x) &&
                (*i)->pGetBoundaryCond()->bIncludeNode(eNodeType)) {
                return (*i);
            }
        }

        return nullptr;
    }

    doublereal BearingGeometry::dGetNodeDistance2D(const Node2D* pNode1, const Node2D* pNode2, index_type iDirection) const
    {
        const Vector<doublereal, 2>& x1 = pNode1->GetPosition2D();
        const Vector<doublereal, 2>& x2 = pNode2->GetPosition2D();

        HYDRO_ASSERT(iDirection >= 1);
        HYDRO_ASSERT(iDirection <= 2);

        const doublereal du = x2(iDirection) - x1(iDirection);
        // Use only the horizontal distance, do not account for vertical distance
        return du;
    }

    RigidBodyBearing::RigidBodyBearing(HydroRootElement* pParent)
        :BearingGeometry(pParent),
         pNode1(0),
         pNode2(0)
    {

    }

    RigidBodyBearing::~RigidBodyBearing()
    {

    }

    void RigidBodyBearing::ParseInput(DataManager* pDM, MBDynParser& HP)
    {
        if (!HP.IsKeyWord("shaft" "node")) {
            silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel() << "): keyword \"shaft\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        pNode1 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

        if (HP.IsKeyWord("offset")) {
            o1_R1 = HP.GetPosRel(ReferenceFrame(pNode1));
        } else {
            o1_R1 = Zero3;
        }

        if (HP.IsKeyWord("orientation")) {
            Rb1 = HP.GetRotRel(ReferenceFrame(pNode1));
        } else {
            Rb1 = Eye3;
        }

        if ( !HP.IsKeyWord("bearing" "node")) {
            silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel() << "): keyword \"bearing\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        pNode2 = pDM->ReadNode<const StructNode, Node::STRUCTURAL>(HP);

        if (HP.IsKeyWord("offset")) {
            o2_R2 = HP.GetPosRel(ReferenceFrame(pNode2));
        } else {
            o2_R2 = Zero3;
        }

        if (HP.IsKeyWord("orientation")) {
            Rb2 = HP.GetRotRel(ReferenceFrame(pNode2));
        } else {
            Rb2 = Eye3;
        }

#if HYDRO_DEBUG > 1
        HYDRO_TRACE("Rb1=\n");

        for (integer i = 1; i <= 3; ++i) {
            for (integer j = 1; j <= 3; ++j) {
                HYDRO_TRACE(Rb1(i, j) << " ");
            }
            HYDRO_TRACE(std::endl);
        }

        HYDRO_TRACE("o1=" << o1_R1 << std::endl);

        HYDRO_TRACE("Rb2=\n");

        for (integer i = 1; i <= 3; ++i) {
            for (integer j = 1; j <= 3; ++j) {
                HYDRO_TRACE(Rb2(i, j) << " ");
            }
            HYDRO_TRACE(std::endl);
        }

        HYDRO_TRACE("o2=" << o2_R2 << std::endl);
#endif
    }

    int RigidBodyBearing::iGetNumConnectedNodes(void) const
    {
        return 2;
    }

    void
    RigidBodyBearing::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
    {
        connectedNodes.resize(iGetNumConnectedNodes());
        connectedNodes[0] = pNode1;
        connectedNodes[1] = pNode2;
    }

    void RigidBodyBearing::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = 12;
        *piNumCols = iGetNumColsWorkSpace(eFunc);

        index_type iNumNodes = 0;

        for (index_type iNode = 0; iNode < pGetParent()->iGetNumNodes(); ++iNode) {
            const Node2D* pNode = pGetParent()->pGetNode(iNode);

            if (pNode->bIsNodeType(Node2D::HYDRAULIC_NODE) && pNode->bIsNodeType(Node2D::UPDATED_NODE)) {
                ++iNumNodes;
            }

            *piNumCols += pNode->iGetNumColsWorkSpace(eFunc);
        }

        const ComplianceModel* pComplianceModel = pGetParent()->pGetMesh()->pGetComplianceModel();

        if (pComplianceModel) {
            *piNumCols += pComplianceModel->iGetNumColsWorkSpace(eFunc, iNumNodes);
        }
    }

    integer RigidBodyBearing::iGetNumColsWorkSpace(grad::FunctionCall eFunc) const
    {
        return eFunc & grad::REGULAR_FLAG ? 12 : 24;
    }

    std::ostream& RigidBodyBearing::PrintLogFile(std::ostream& os) const
    {
        const long precision = os.precision();

        os.precision(16);

        os << pNode1->GetLabel() << ' ' << o1_R1 << ' ' << Rb1 << ' '
           << pNode2->GetLabel() << ' ' << o2_R2 << ' ' << Rb2 << ' ';

        os.precision(precision);

        return os;
    }

    void
    RigidBodyBearing::SaveReactionForce(const Vector<doublereal, 3>& F1,
                                        const Vector<doublereal, 3>& M1,
                                        const Vector<doublereal, 3>& F2,
                                        const Vector<doublereal, 3>& M2)
    {
        this->F1 = F1;
        this->M1 = M1;
        this->F2 = F2;
        this->M2 = M2;
    }

    void
    RigidBodyBearing::SaveReactionForce(const Vector<Gradient<0>, 3>& F1,
                                        const Vector<Gradient<0>, 3>& M1,
                                        const Vector<Gradient<0>, 3>& F2,
                                        const Vector<Gradient<0>, 3>& M2)
    {

    }

    std::ostream& RigidBodyBearing::Output(std::ostream& os) const
    {
        os << F1 << ' '
           << M1 << ' '
           << F2 << ' '
           << M2 << ' ';

        return os;
    }

    CylindricalBearing::CylindricalBearing(HydroRootElement* pParent)
        :RigidBodyBearing(pParent), hmin(0.)
    {

    }

    CylindricalBearing::~CylindricalBearing()
    {

    }

    void CylindricalBearing::ParseInput(DataManager* pDM, MBDynParser& HP) {
        if ( !HP.IsKeyWord("bearing" "width") )
        {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"bearing width\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        b = HP.GetReal();

        if (!HP.IsKeyWord("shaft" "diameter"))
        {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"shaft diameter\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        r = 0.5 * HP.GetReal();

        if ( !HP.IsKeyWord("bearing" "diameter") )
        {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"bearing diameter\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        R = 0.5 * HP.GetReal();

        if (HP.IsKeyWord("minimum" "clearance")) {
            if (HP.IsKeyWord("relative")) {
                hmin = HP.GetReal() * (R - r);
            } else if (HP.IsKeyWord("absolute")) {
                hmin = HP.GetReal();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"relative\" or \"absolute\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        } else {
            hmin = 1e-3 * (R - r);
        }

        if (HP.IsKeyWord("pockets")) {
            bool bHavePocket = false;

            if (HP.IsKeyWord("shaft")) {
                ReadPockets(HP, rgPocketsShaft);
                bHavePocket = true;
            }

            if (HP.IsKeyWord("bearing")) {
                ReadPockets(HP, rgPocketsBearing);
                bHavePocket = true;
            }

            if (!bHavePocket) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"bearing\" or keyword \"shaft\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        RigidBodyBearing::ParseInput(pDM, HP);
    }

    void CylindricalBearing::GetClosestDistance2D(const Vector<doublereal, 2>& x1, const Vector<doublereal, 2>& x2, Vector<doublereal, 2>& dx) const {
        dx = x2 - x1;

        const doublereal c = 2. * M_PI * dGetMeshRadius();

        if (std::fabs(dx(1) + c) < fabs(dx(1))) {
            dx(1) += c;
        } else if (std::fabs(dx(1) - c) < fabs(dx(1))) {
            dx(1) -= c;
        }
    }

    void CylindricalBearing::GetPosition3D(const Vector<doublereal, 2>& x1,
                                           Vector<doublereal, 3>& v1) const
    {
        GetPosition3DTpl(x1, v1);
    }

    void CylindricalBearing::GetPosition3D(const Vector<Gradient<0>, 2>& x1,
                                           Vector<Gradient<0>, 3>& v1) const
    {
        GetPosition3DTpl(x1, v1);
    }

    template <typename T>
    inline void
    CylindricalBearing::GetPosition3DTpl(const Vector<T, 2>& x1,
                                         Vector<T, 3>& v1) const
    {
        typedef Vector<doublereal, 2> CVec2;
        typedef Vector<T, 3> VVec3;
        typedef Matrix<doublereal, 3, 3> CMat3x3;

        const doublereal r = dGetMeshRadius();
        const T Phi1 = x1(1) / r;
        const T& z1 = x1(2);
        const Pocket* const pPocket = pFindMeshPocket(CVec2(x1));
        T dy;

        if (pPocket != 0) {
            pPocket->GetHeight(x1, dy);
        } else {
            dy = 0.;
        }

        const VVec3 v1_Rb((r + dy) * cos(Phi1),
                          (r + dy) * sin(Phi1),
                          z1);

        const CMat3x3& Rb = GetOrientationMeshNode();

        v1 = Rb * v1_Rb;
    }

    doublereal CylindricalBearing::dGetPocketHeightMesh(const Vector<doublereal, 2>& x) const
    {
        const Pocket* const pPocket = pFindMeshPocket(x);
        doublereal dy;

        if (pPocket) {
            pPocket->GetHeight(x, dy);
        } else {
            dy = 0.;
        }

        return dy;
    }

    void CylindricalBearing::GetTangentCoordSys(const Vector<doublereal, 2>& x,
                                                Matrix<doublereal, 3, 3>& Rt) const
    {
        GetTangentCoordSysTpl(x, Rt);
    }

    void
    CylindricalBearing::GetTangentCoordSys(const Vector<Gradient<0>, 2>& x,
                                           Matrix<Gradient<0>, 3, 3>& Rt) const
    {
        GetTangentCoordSysTpl(x, Rt);
    }

    template <typename T>
    void CylindricalBearing::GetTangentCoordSysTpl(const Vector<T, 2>& x1,
                                                   Matrix<T, 3, 3>& Rbt) const
    {
        typedef Matrix<doublereal, 3, 3> CMat3x3;
        typedef Vector<doublereal, 2> CVec2;
        typedef Matrix<T, 3, 3> Mat3x3;
        typedef Vector<T, 3> Vec3;

        const doublereal r = dGetMeshRadius();
        const T Phi1 = x1(1) / r;

        const Mat3x3 Rt{-sin(Phi1),  cos(Phi1), T{0.},
                -cos(Phi1), -sin(Phi1), T{0.},
                    T{0.},      T{0.}, T{1.}};

        const CMat3x3& Rb = GetOrientationMeshNode();

        const Pocket* const pPocket = pFindMeshPocket(CVec2(x1));

        if (pPocket) {
            T tan_beta, tan_gamma;

            pPocket->GetHeightDerX(x1, tan_beta);
            pPocket->GetHeightDerZ(x1, tan_gamma);

            const Mat3x3 dR{   T{1.},  -tan_beta,     T{0.},
                                           tan_beta,      T{1.}, tan_gamma,
                                                                     T{0.}, -tan_gamma,     T{1.}};

            Rbt = Rb * Mat3x3(dR * Rt);
        } else {
            Rbt = Rb * Rt;
        }
    }

    void CylindricalBearing::GetStructNodeOffset(const HydroNode* pHydroNode, Vector<doublereal, 3>& v) const
    {
        const Vector<doublereal, 2>& x = pHydroNode->GetPosition2D();
        const doublereal r = dGetMeshRadius();
        const doublereal Phi = x(1) / r;
        const doublereal z = x(2);
        const Pocket* const pPocket = pFindMeshPocket(x);

        doublereal dy;

        if (pPocket == nullptr) {
            dy = 0.;
        } else {
            pPocket->GetHeight(x, dy);
        }

        const Vector<doublereal, 3> v_Rb((r + dy) * cos(Phi),
                                         (r + dy) * sin(Phi),
                                         z);

        const Matrix<doublereal, 3, 3>& Rb = GetOrientationMeshNode();

        v = Rb * v_Rb;
    }

    void CylindricalBearing::Update(doublereal dCoef, enum FunctionCall func, LocalDofMap* pDofMap)
    {
        oDofMap.Reset(func);
    }

    std::ostream& CylindricalBearing::PrintLogFile(std::ostream& os) const
    {
        os << 2 * dGetShaftRadius() << ' '
           << 2 * dGetBearingRadius() << ' '
           << dGetBearingWidth() << ' ';

        RigidBodyBearing::PrintLogFile(os);

        return os;
    }

    doublereal CylindricalBearing::dGetMinClearance() const
    {
        return hmin;
    }

    doublereal CylindricalBearing::dGetReferenceClearance() const
    {
        return dGetBearingRadius() - dGetShaftRadius();
    }

    const Pocket* CylindricalBearing::pFindBearingPocket(const Vector<doublereal, 2>& x) const
    {
        return pFindPocket(x, rgPocketsBearing);
    }

    const Pocket* CylindricalBearing::pFindShaftPocket(const Vector<doublereal, 2>& x) const
    {
        return pFindPocket(x, rgPocketsShaft);
    }

    const Pocket* CylindricalBearing::pFindPocket(const Vector<doublereal, 2>& x, const PocketVector& rgPockets)
    {
        for (ConstPocketIterator i = rgPockets.begin(); i != rgPockets.end(); ++i) {
            if ((*i)->pGetGeometry()->bPointIsInside(x)) {
                return i->get();
            }
        }

        return 0;
    }

    void CylindricalBearing::ReadPockets(MBDynParser& HP, PocketVector& rgPockets)
    {
        const integer iNumPockets = HP.GetInt();

        rgPockets.reserve(2 * iNumPockets);

        const doublereal r = dGetMeshRadius();
        const doublereal c = 2 * M_PI * r;

        for (integer i = 0; i < iNumPockets; ++i) {
            std::unique_ptr<Pocket> pPocket(Pocket::Read(pGetParent(), HP, this));
            Vector<doublereal, 2> x = pPocket->pGetGeometry()->GetPosition();

            if (x(1) < 0. || x(1) > c) {
                // Note: This case makes no sense because of the periodic nature
                // of the cylindrical bearing we can always write
                //      x(1) = x(1) - c
                // or
                //      x(1) = x(1) + c
                // Note: x(2) may be outside the cylindrical bearing because
                // the pocket could be at the boundary of the bearing
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): pocket (" << i + 1
                            << ") is outside the bearing are at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (x(1) >= 0.5 * c) {
                x(1) -= c;
            } else {
                x(1) += c;
            }

            rgPockets.push_back(pPocket->Clone(x));
            rgPockets.push_back(std::move(pPocket));
        }
    }

    CylindricalMeshAtShaft::CylindricalMeshAtShaft(HydroRootElement* pParent)
        :CylindricalBearing(pParent),
         oBound(*this),
         oBound_grad(*this)
    {

    }

    CylindricalMeshAtShaft::~CylindricalMeshAtShaft()
    {

    }

    BearingGeometry::Type CylindricalMeshAtShaft::GetType() const
    {
        return CYLINDRICAL_MESH_AT_SHAFT;
    }

    void CylindricalMeshAtShaft::Initialize()
    {
        oBound.Initialize();
        oBound_grad.Initialize();

        integer iNumRows, iNumCols;

        // Will be always higher than WorkSpaceDim
        WorkSpaceDim(&iNumRows, &iNumCols, grad::INITIAL_ASS_FLAG);

        const unsigned iNumDof = pGetParent()->iGetNumDof() + iNumCols;

        for (integer i = 1; i <= 3; ++i) {
            oReaction_grad.F1_R1(i).Reserve(iNumDof);
            oReaction_grad.M1_R1(i).Reserve(iNumDof);
            oReaction_grad.F2_R1(i).Reserve(iNumDof);
            oReaction_grad.M2_R1(i).Reserve(iNumDof);

            oReaction_grad.F1_R1(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.M1_R1(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.F2_R1(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.M2_R1(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
        }
    }

    template <typename T>
    CylindricalMeshAtShaft::Boundary<T>::Boundary(const CylindricalMeshAtShaft& rParent)
        :rParent(rParent)
    {

    }

    template <typename T>
    void CylindricalMeshAtShaft::Boundary<T>::Initialize()
    {
        const CMat3x3& Rb2 = rParent.GetOrientationNode2();
        const CVec3& o2 = rParent.GetOffsetNode2();

        Rb2T_o2 = Transpose(Rb2) * o2;
    }

    template <typename T>
    void
    CylindricalMeshAtShaft::Boundary<T>::Update(
        doublereal dCoef,
        enum FunctionCall func,
        LocalDofMap* pDofMap)
    {
        const StructNode* const pNode1 = rParent.pGetNode1();
        const StructNode* const pNode2 = rParent.pGetNode2();

        pNode1->GetXCurr(X1, dCoef, func, pDofMap);
        pNode1->GetRCurr(R1, dCoef, func, pDofMap);
        pNode1->GetVCurr(X1P, dCoef, func, pDofMap);
        pNode1->GetWCurr(omega1, dCoef, func, pDofMap);

        pNode2->GetXCurr(X2, dCoef, func, pDofMap);
        pNode2->GetRCurr(R2, dCoef, func, pDofMap);
        pNode2->GetVCurr(X2P, dCoef, func, pDofMap);
        pNode2->GetWCurr(omega2, dCoef, func, pDofMap);

        const CMat3x3& Rb2 = rParent.GetOrientationNode2();

        Rb2T_R2T = Transpose(Rb2) * Transpose(R2);
    }

    template <typename T>
    void
    CylindricalMeshAtShaft::Boundary<T>::GetBoundaryConditions(HydroNode* pNode,
                                                               T& h,
                                                               T& dh_dt,
                                                               Vector<T, 2>& U1,
                                                               Vector<T, 2>& U2,
                                                               Vector<T, 2>& U,
                                                               doublereal dCoef,
                                                               enum FunctionCall func,
                                                               LocalDofMap* pDofMap) const
    {
        const CVec3& Rb1_v1 = pNode->GetPosition3D();
        const CVec3& o1 = rParent.GetOffsetNode1();
        const CVec3& o2 = rParent.GetOffsetNode2();
        const CMat3x3& Rb2 = rParent.GetOrientationNode2();

        const VVec3 a0 = R1 * VVec3(o1 + Rb1_v1);
        const VVec3 a2 = X1 - X2 + a0;
        const VVec3 a1 = X1P + Cross(omega1, a0)
            - X2P - Cross(omega2, a2);
        const VVec3 b = Rb2T_R2T * a2 - Rb2T_o2;
        const doublereal R = rParent.dGetBearingRadius();

        const T a5 = b(1) * b(1) + b(2) * b(2);
        const T a3 = sqrt(a5);

        const T cos_Phi2 = b(1) / a3;
        const T sin_Phi2 = b(2) / a3;
        const T& z2 = b(3);

        T Phi2 = atan2(sin_Phi2, cos_Phi2);

        if (Phi2 < 0.) {
            // Attention: pFindBearingPockets assumes that Phi1 is in range 0 ... 2 * pi
            Phi2 += 2 * M_PI;
        }

        const VVec2 x2(Phi2 * R, z2);
        const Pocket* const pPocket2 = rParent.pFindBearingPocket(CVec2(x2));

        T Deltay2, dDeltay2_dx2, dDeltay2_dz2;

        if (pPocket2 != 0) {
            pPocket2->GetHeight(x2, Deltay2);
            pPocket2->GetHeightDerX(x2, dDeltay2_dx2);
            pPocket2->GetHeightDerZ(x2, dDeltay2_dz2);
        } else {
            Deltay2 = 0.;
            dDeltay2_dx2 = 0.;
            dDeltay2_dz2 = 0.;
        }

        const LubricationGroove* const pMovingGroove = rParent.pFindMovingLubricationGroove(CVec2(x2), pNode->GetNodePhysics());

#if HYDRO_DEBUG > 1
        if (pMovingGroove) {
            HYDRO_TRACE("node(" << pNode->iGetNodeNumber() + 1
                        << ") affected by moving boundary condition at x2("
                        << CVec2(x2) << ")\n");
        }
#endif

        pNode->SetMovingPressBoundCond(pMovingGroove ? pMovingGroove->pGetBoundaryCond() : nullptr);

        T w, dw_dt;

        pNode->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap);

        const T h0 = R + Deltay2 - a3;
        h = h0 + w;

        const VVec3 db_dt = Rb2T_R2T * a1;
        const T dx2_dt = R * (b(1) * db_dt(2) - db_dt(1) * b(2)) / a5;
        const T& dz2_dt = db_dt(3);
        const T dDeltay2_dt = dDeltay2_dx2 * dx2_dt + dDeltay2_dz2 * dz2_dt;
        dh_dt = dDeltay2_dt - (b(1) * db_dt(1) + b(2) * db_dt(2)) / a3 + dw_dt;

        const VVec3 v2((R + Deltay2) * cos_Phi2,
                       (R + Deltay2) * sin_Phi2,
                       z2);

        const VVec3 vh(h0 * cos_Phi2,
                       h0 * sin_Phi2,
                       T(0.));

        const VVec3 P1Dot = X1P + Cross(omega1, a0);
        const VVec3 P2Dot = X2P + Cross(omega2, VVec3(R2 * VVec3(o2 + Rb2 * v2))) - Cross(omega1, VVec3(R2 * VVec3(Rb2 * vh)));

        const CMat3x3& Rbt1 = pNode->GetTangentCoordSys();

        const VVec3 P1Dot_R1 = Transpose(R1) * P1Dot;
        const VVec3 P2Dot_R1 = Transpose(R1) * P2Dot;

        U1(1) = Dot(Rbt1.GetCol(1), P1Dot_R1);
        U1(2) = Dot(Rbt1.GetCol(3), P1Dot_R1);
        U2(1) = Dot(Rbt1.GetCol(1), P2Dot_R1);
        U2(2) = Dot(Rbt1.GetCol(3), P2Dot_R1);

        U = (U2 - U1) * 0.5;
    }

    template <typename T>
    void CylindricalMeshAtShaft::Boundary<T>::GetMovingMeshOffset(Vector<T, 2>& x, LocalDofMap* pDofMap) const
    {
        const CVec3& o1 = rParent.GetOffsetNode1();
        const CVec3& o2 = rParent.GetOffsetNode2();
        const CMat3x3& Rb1 = rParent.GetOrientationNode1();
        const CMat3x3& Rb2 = rParent.GetOrientationNode2();
        const doublereal r = rParent.dGetMeshRadius();

        const VVec3 R2_Rb2_e1 = R2 * Rb2.GetCol(1);

        const T dx = r * atan2(Dot(R1 * Rb1.GetCol(2), R2_Rb2_e1),
                         Dot(R1 * Rb1.GetCol(1), R2_Rb2_e1));

        const T dz = Dot(Rb1.GetCol(3), Transpose(R1) * VVec3(X2 + R2 * o2 - X1) - o1);

        Copy(x(1), dx, pDofMap); // FIXME: Which class is responsible for copying Gradient's
        Copy(x(2), dz, pDofMap);
    }

    void CylindricalMeshAtShaft::GetMovingMeshOffset(Vector<doublereal, 2>& x, LocalDofMap* pDofMap) const
    {
        oBound.GetMovingMeshOffset(x, pDofMap);
    }

    void CylindricalMeshAtShaft::GetMovingMeshOffset(Vector<Gradient<0>, 2>& x, LocalDofMap* pDofMap) const
    {
        oBound_grad.GetMovingMeshOffset(x, pDofMap);
    }

    void CylindricalMeshAtShaft::GetBoundaryConditions(HydroNode* pNode, doublereal& h, doublereal& dh_dt, Vector<doublereal, 2>& U1, Vector<doublereal, 2>& U2, Vector<doublereal, 2>& U, doublereal dCoef, enum FunctionCall func, LocalDofMap* pDofMap) const
    {
        oBound.GetBoundaryConditions(pNode, h, dh_dt, U1, U2, U, dCoef, func, pDofMap);
    }

    void CylindricalMeshAtShaft::GetBoundaryConditions(HydroNode* pNode, Gradient<0>& h, Gradient<0>& dh_dt, Vector<Gradient<0>, 2>& U1, Vector<Gradient<0>, 2>& U2, Vector<Gradient<0>, 2>& U, doublereal dCoef, enum FunctionCall func, LocalDofMap* pDofMap) const
    {
        oBound_grad.GetBoundaryConditions(pNode, h, dh_dt, U1, U2, U, dCoef, func, pDofMap);
    }

    doublereal CylindricalMeshAtShaft::dGetMeshRadius() const
    {
        return dGetShaftRadius();
    }

    const Matrix<doublereal, 3, 3>&
    CylindricalMeshAtShaft::GetOrientationMeshNode() const
    {
        return GetOrientationNode1();
    }

    const Pocket*
    CylindricalMeshAtShaft::pFindMeshPocket(const Vector<doublereal, 2>& x) const
    {
        return pFindShaftPocket(x);
    }

    void
    CylindricalMeshAtShaft::AddReactionForce(const Vector<doublereal, 2>& x,
                                             const Vector<doublereal, 3>& v,
                                             const Matrix<doublereal, 3, 3>& Rt,
                                             const Vector<doublereal, 3>& dF_0_Rt,
                                             const Vector<doublereal, 3>& dF_h_Rt,
                                             const Vector<doublereal, 2>& dM_h_Rt,
                                             doublereal dCoef,
                                             enum FunctionCall func,
                                             LocalDofMap* pDofMap)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        AddReactionForce(x, v, Rt, dF_0_Rt, dF_h_Rt, dM_h_Rt, dCoef, func, pDofMap, oReaction);

#if CREATE_PROFILE == 1
        profile.dtAddForce[PROF_RES] += mbdyn_clock_time() - start;
#endif
    }

    void
    CylindricalMeshAtShaft::AddReactionForce(const Vector<doublereal, 2>& x,
                                             const Vector<doublereal, 3>& v,
                                             const Matrix<doublereal, 3, 3>& Rt,
                                             const Vector<Gradient<0>, 3>& dF_0_Rt,
                                             const Vector<Gradient<0>, 3>& dF_h_Rt,
                                             const Vector<Gradient<0>, 2>& dM_h_Rt,
                                             doublereal dCoef,
                                             enum FunctionCall func,
                                             LocalDofMap* pDofMap)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        AddReactionForce(x, v, Rt, dF_0_Rt, dF_h_Rt, dM_h_Rt, dCoef, func, pDofMap, oReaction_grad);

#if CREATE_PROFILE == 1
        profile.dtAddForce[PROF_JAC] += mbdyn_clock_time() - start;
#endif
    }

    template <typename T>
    void CylindricalMeshAtShaft::ReactionForce<T>::Reset()
    {
        F1_R1.Reset();
        M1_R1.Reset();
        F2_R1.Reset();
        M2_R1.Reset();
    }

    template <typename T>
    inline void
    CylindricalMeshAtShaft::AddReactionForce(const Vector<doublereal, 2>& x1,
                                             const Vector<doublereal, 3>& Rb1_v1,
                                             const Matrix<doublereal, 3, 3>& Rbt1,
                                             const Vector<T, 3>& dF2_0_Rt1,
                                             const Vector<T, 3>& dF1_h_Rt1,
                                             const Vector<T, 2>& dM1_h_Rt1,
                                             doublereal dCoef,
                                             enum FunctionCall func,
                                             LocalDofMap* pDofMap,
                                             ReactionForce<T>& oReact)
    {
        typedef Vector<T, 3> VVec3;

#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif
        const VVec3 dF1_R1 = Rbt1 * dF1_h_Rt1;
        const VVec3 dM1_h_R1 = Rbt1.GetCol(1) * dM1_h_Rt1(1) + Rbt1.GetCol(3) * dM1_h_Rt1(2);

        oReact.F1_R1 += dF1_R1;
        oReact.M1_R1 += Cross(Rb1_v1, dF1_R1) + dM1_h_R1;

        const VVec3 dF2_R1 = Rbt1 * dF2_0_Rt1;

        oReact.F2_R1 += dF2_R1;
        oReact.M2_R1 += Cross(Rb1_v1, dF2_R1) - dM1_h_R1;

#if GRADIENT_DEBUG >= 2
        for (integer i = 1; i <= 3; ++i)
        {
            std::cerr << "size(dF1_R1(" << i << "))=" << iGetGradientVectorSize(dF1_R1(i)) << std::endl;
            std::cerr << "size(dF2_R1(" << i << "))=" << iGetGradientVectorSize(dF2_R1(i)) << std::endl;
        }
#endif

#if CREATE_PROFILE == 1
        profile.dtOperatorPlus[func == REGULAR_JAC ? PROF_JAC : PROF_RES] += mbdyn_clock_time() - start;
#endif
    }

    const CylindricalMeshAtShaft::ReactionForce<doublereal>&
    CylindricalMeshAtShaft::GetReactionForce(const doublereal&) const
    {
        return oReaction;
    }

    const CylindricalMeshAtShaft::ReactionForce<Gradient<0> >&
    CylindricalMeshAtShaft::GetReactionForce(const Gradient<0>&) const
    {
        return oReaction_grad;
    }

    void
    CylindricalMeshAtShaft::Update(
        doublereal dCoef,
        enum FunctionCall func,
        LocalDofMap* pDofMap)
    {
        CylindricalBearing::Update(dCoef, func, pDofMap);

        switch (func) {
        case REGULAR_RES:
        case INITIAL_DER_RES:
        case INITIAL_ASS_RES:
            oBound.Update(dCoef, func, pDofMap);
            oReaction.Reset();
            break;

        case REGULAR_JAC:
        case INITIAL_DER_JAC:
        case INITIAL_ASS_JAC:
            oBound_grad.Update(dCoef, func, pDofMap);
            oReaction_grad.Reset();
            break;

        default:
            HYDRO_ASSERT(0);
        }
    }

    void
    CylindricalMeshAtShaft::AssRes(SubVectorHandler& WorkVec,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    CylindricalMeshAtShaft::AssJac(SparseSubMatrixHandler& WorkMat,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    CylindricalMeshAtShaft::InitialAssRes(SubVectorHandler& WorkVec,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    CylindricalMeshAtShaft::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }


    template <typename T>
    void CylindricalMeshAtShaft::AssRes(GradientAssVec<T>& WorkMat,
                                        doublereal dCoef,
                                        const GradientVectorHandler<T>& XCurr,
                                        const GradientVectorHandler<T>& XPrimeCurr,
                                        enum FunctionCall func)
    {
        UnivAssRes(WorkMat, dCoef, XCurr, func);
    }


    template <typename T>
    void CylindricalMeshAtShaft::InitialAssRes(GradientAssVec<T>& WorkMat,
                                               const GradientVectorHandler<T>& XCurr,
                                               enum FunctionCall func)
    {
        UnivAssRes(WorkMat, 1., XCurr, func);
    }

    template <typename T>
    void CylindricalMeshAtShaft::UnivAssRes(GradientAssVec<T>& WorkVec,
                                            doublereal dCoef,
                                            const GradientVectorHandler<T>& XCurr,
                                            enum FunctionCall func)
    {
#if GRADIENT_DEBUG >= 2
        std::cerr << "func=" << func << std::endl;
#endif

        const doublereal dInitAss = pGetParent()->dGetStartupFactor();

        typedef Vector<T, 3> VVec3;
        typedef Matrix<T, 3, 3> VMat3x3;
        typedef Vector<doublereal, 3> CVec3;
        typedef Matrix<doublereal, 3, 3> CMat3x3;

        const CVec3& o1 = GetOffsetNode1();
        const CVec3& o2 = GetOffsetNode2();
        const CMat3x3& Rb1 = GetOrientationNode1();
        const CMat3x3& Rb2 = GetOrientationNode2();

        VVec3 X1, X2;
        VMat3x3 R1, R2;

        const StructNode* const pNode1 = pGetNode1();

        pNode1->GetXCurr(X1, dCoef, func, &oDofMap);
        pNode1->GetRCurr(R1, dCoef, func, &oDofMap);

        const StructNode* const pNode2 = pGetNode2();

        pNode2->GetXCurr(X2, dCoef, func, &oDofMap);
        pNode2->GetRCurr(R2, dCoef, func, &oDofMap);

        const ReactionForce<T>& oReact = GetReactionForce(X1(1));

        const T lambda = -Dot(Rb1.GetCol(3), VVec3(Transpose(R1) * VVec3(X1 - X2 - R2 * o2) + o1))
            / Dot(Rb1.GetCol(3), VVec3(Transpose(R1) * VVec3(R2 * Rb2.GetCol(3))));
        const VVec3 F1 = (R1 * oReact.F1_R1) * dInitAss;
        const VVec3 M1 = (R1 * oReact.M1_R1) * dInitAss + Cross(VVec3(R1 * o1), F1);
        const VVec3 F2 = (R1 * oReact.F2_R1) * dInitAss;
        const VVec3 M2 = (R1 * oReact.M2_R1) * dInitAss + Cross(VVec3(R2 * VVec3(o2 + Rb2.GetCol(3) * lambda)), F2);

        const integer iFirstMomIndexNode1 = (func & INITIAL_ASS_FLAG) ? pNode1->iGetFirstPositionIndex() : pNode1->iGetFirstMomentumIndex();
        const integer iFirstMomIndexNode2 = (func & INITIAL_ASS_FLAG) ? pNode2->iGetFirstPositionIndex() : pNode2->iGetFirstMomentumIndex();

#if GRADIENT_DEBUG >= 2
        std::cerr << "F1=" << F1 << std::endl;
#endif

        SaveReactionForce(F1, M1, F2, M2);

        CHECK_NUM_COLS_WORK_SPACE(this, func, F1, iFirstMomIndexNode1 + 1);
        CHECK_NUM_COLS_WORK_SPACE(this, func, M1, iFirstMomIndexNode1 + 4);
        CHECK_NUM_COLS_WORK_SPACE(this, func, F2, iFirstMomIndexNode2 + 1);
        CHECK_NUM_COLS_WORK_SPACE(this, func, M2, iFirstMomIndexNode2 + 4);

        WorkVec.AddItem(iFirstMomIndexNode1 + 1, F1);
        WorkVec.AddItem(iFirstMomIndexNode1 + 4, M1);
        WorkVec.AddItem(iFirstMomIndexNode2 + 1, F2);
        WorkVec.AddItem(iFirstMomIndexNode2 + 4, M2);
    }

    enum LubricationGroove::Type CylindricalMeshAtShaft::ReadLubricationGrooveType(MBDynParser& HP) const
    {
        if (HP.IsKeyWord("at" "shaft")) {
            return LubricationGroove::FIXED;
        } else if (HP.IsKeyWord("at" "bearing")) {
            return LubricationGroove::MOVING;
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"at shaft\" or \"at bearing\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    bool
    CylindricalMeshAtShaft::bGetPrivateData(HydroRootBase::PrivateDataType eType,
                                            doublereal& dPrivData) const
    {
        const Matrix<doublereal, 3, 3>& Rb1 = GetOrientationNode1();

        switch (eType) {
        case HydroRootBase::PD_F1x:
        case HydroRootBase::PD_F1y:
        case HydroRootBase::PD_F1z: {
            const index_type iIndex = eType - HydroRootBase::PD_F1x + 1;
            dPrivData = Dot(Rb1.GetCol(iIndex), oReaction.F1_R1);
            return true;
        }
        case HydroRootBase::PD_M1x:
        case HydroRootBase::PD_M1y:
        case HydroRootBase::PD_M1z: {
            const index_type iIndex = eType - HydroRootBase::PD_M1x + 1;
            dPrivData = Dot(Rb1.GetCol(iIndex), oReaction.M1_R1);
            return true;
        }
        case HydroRootBase::PD_F2x:
        case HydroRootBase::PD_F2y:
        case HydroRootBase::PD_F2z: {
            const index_type iIndex = eType - HydroRootBase::PD_F2x + 1;
            dPrivData = Dot(Rb1.GetCol(iIndex), oReaction.F2_R1);
            return true;
        }
        case HydroRootBase::PD_M2x:
        case HydroRootBase::PD_M2y:
        case HydroRootBase::PD_M2z: {
            const index_type iIndex = eType - HydroRootBase::PD_M2x + 1;
            dPrivData = Dot(Rb1.GetCol(iIndex), oReaction.M2_R1);
            return true;
        }
        default:
            return false;
        }
    }

    CylindricalMeshAtBearing::CylindricalMeshAtBearing(HydroRootElement* pParent)
        :CylindricalBearing(pParent),
         oBound(*this),
         oBound_grad(*this)
    {

    }

    CylindricalMeshAtBearing::~CylindricalMeshAtBearing()
    {

    }

    BearingGeometry::Type CylindricalMeshAtBearing::GetType() const
    {
        return CYLINDRICAL_MESH_AT_BEARING;
    }

    void CylindricalMeshAtBearing::Initialize()
    {
        oBound.Initialize();
        oBound_grad.Initialize();

        integer iNumRows, iNumCols;

        // Will be always higher than WorkSpaceDim
        WorkSpaceDim(&iNumRows, &iNumCols, grad::INITIAL_ASS_FLAG);

        const unsigned iNumDof = pGetParent()->iGetNumDof() + iNumCols;

        for (integer i = 1; i <= 3; ++i) {
            oReaction_grad.F1_R2(i).Reserve(iNumDof);
            oReaction_grad.M1_R2(i).Reserve(iNumDof);
            oReaction_grad.F2_R2(i).Reserve(iNumDof);
            oReaction_grad.M2_R2(i).Reserve(iNumDof);

            oReaction_grad.F1_R2(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.M1_R2(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.F2_R2(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
            oReaction_grad.M2_R2(i).DerivativeResizeReset(&oDofMap, 0, MapVectorBase::LOCAL, 0.);
        }
    }

    void
    CylindricalMeshAtBearing::GetBoundaryConditions(HydroNode* pNode,
                                                    doublereal& h,
                                                    doublereal& dh_dt,
                                                    Vector<doublereal, 2>& U1,
                                                    Vector<doublereal, 2>& U2,
                                                    Vector<doublereal, 2>& U,
                                                    doublereal dCoef,
                                                    enum FunctionCall func,
                                                    LocalDofMap* pDofMap) const
    {
        oBound.GetBoundaryConditions(pNode, h, dh_dt, U1, U2, U, dCoef, func, pDofMap);
    }

    void
    CylindricalMeshAtBearing::GetBoundaryConditions(HydroNode* pNode,
                                                    Gradient<0>& h,
                                                    Gradient<0>& dh_dt,
                                                    Vector<Gradient<0>, 2>& U1,
                                                    Vector<Gradient<0>, 2>& U2,
                                                    Vector<Gradient<0>, 2>& U,
                                                    doublereal dCoef,
                                                    enum FunctionCall func,
                                                    LocalDofMap* pDofMap) const
    {
        oBound_grad.GetBoundaryConditions(pNode, h, dh_dt, U1, U2, U, dCoef, func, pDofMap);
    }

    void
    CylindricalMeshAtBearing::Update(
        doublereal dCoef,
        enum FunctionCall func,
        LocalDofMap* pDofMap)
    {
        CylindricalBearing::Update(dCoef, func, pDofMap);

        switch (func) {
        case REGULAR_RES:
        case INITIAL_DER_RES:
        case INITIAL_ASS_RES:
            oBound.Update(dCoef, func, pDofMap);
            oReaction.Reset();
            break;

        case REGULAR_JAC:
        case INITIAL_DER_JAC:
        case INITIAL_ASS_JAC:
            oBound_grad.Update(dCoef, func, pDofMap);
            oReaction_grad.Reset();
            break;

        default:
            HYDRO_ASSERT(0);
        }
    }

    void
    CylindricalMeshAtBearing::AssRes(SubVectorHandler& WorkVec,
                                     doublereal dCoef,
                                     const VectorHandler& XCurr,
                                     const VectorHandler& XPrimeCurr,
                                     GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    CylindricalMeshAtBearing::AssJac(SparseSubMatrixHandler& WorkMat,
                                     doublereal dCoef,
                                     const VectorHandler& XCurr,
                                     const VectorHandler& XPrimeCurr,
                                     GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    CylindricalMeshAtBearing::InitialAssRes(SubVectorHandler& WorkVec,
                                            const VectorHandler& XCurr,
                                            enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    CylindricalMeshAtBearing::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                            const VectorHandler& XCurr,
                                            GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }


    template <typename T>
    void CylindricalMeshAtBearing::AssRes(GradientAssVec<T>& WorkMat,
                                          doublereal dCoef,
                                          const GradientVectorHandler<T>& XCurr,
                                          const GradientVectorHandler<T>& XPrimeCurr,
                                          FunctionCall func)
    {
        UnivAssRes(WorkMat, dCoef, XCurr, func);
    }


    template <typename T>
    void CylindricalMeshAtBearing::InitialAssRes(GradientAssVec<T>& WorkMat,
                                                 const GradientVectorHandler<T>& XCurr,
                                                 FunctionCall func)
    {
        UnivAssRes(WorkMat, 1., XCurr, func);
    }

    template <typename T>
    void
    CylindricalMeshAtBearing::UnivAssRes(GradientAssVec<T>& WorkVec,
                                         doublereal dCoef,
                                         const GradientVectorHandler<T>& XCurr,
                                         FunctionCall func)
    {
        typedef Vector<T, 3> VVec3;
        typedef Matrix<T, 3, 3> VMat3x3;
        typedef Vector<doublereal, 3> CVec3;
        typedef Matrix<doublereal, 3, 3> CMat3x3;

        VVec3 X1, X2;
        VMat3x3 R1, R2;

        const StructNode* const pNode1 = pGetNode1();
        const StructNode* const pNode2 = pGetNode2();

        pNode1->GetXCurr(X1, dCoef, func, &oDofMap);
        pNode1->GetRCurr(R1, dCoef, func, &oDofMap);

        pNode2->GetXCurr(X2, dCoef, func, &oDofMap);
        pNode2->GetRCurr(R2, dCoef, func, &oDofMap);

        const CVec3& o1 = GetOffsetNode1();
        const CVec3& o2 = GetOffsetNode2();
        const CMat3x3& Rb1 = GetOrientationNode1();
        const CMat3x3& Rb2 = GetOrientationNode2();

        const ReactionForce<T>& oReact = GetReactionForce(X1(1));

        const doublereal dInitAss = pGetParent()->dGetStartupFactor();

        const T lambda = -Dot(Rb2.GetCol(3), VVec3(Transpose(R2) * VVec3(X1 + VVec3(R1 * o1) - X2) - o2))
            / Dot(Rb2.GetCol(3), VVec3(Transpose(R2) * VVec3(R1 * Rb1.GetCol(3))));
        const VVec3 F1 = (R2 * oReact.F1_R2) * dInitAss;
        const VVec3 M1 = (R2 * oReact.M1_R2) * dInitAss + Cross(VVec3(R1 * VVec3(o1 + Rb1.GetCol(3) * lambda)), F1);
        const VVec3 F2 = (R2 * oReact.F2_R2) * dInitAss;
        const VVec3 M2 = (R2 * oReact.M2_R2) * dInitAss + Cross(R2 * o2, F2);

        SaveReactionForce(F1, M1, F2, M2);

        const integer iFirstMomIndexNode1 = (func & INITIAL_ASS_FLAG) ? pNode1->iGetFirstPositionIndex() : pNode1->iGetFirstMomentumIndex();
        const integer iFirstMomIndexNode2 = (func & INITIAL_ASS_FLAG) ? pNode2->iGetFirstPositionIndex() : pNode2->iGetFirstMomentumIndex();

        CHECK_NUM_COLS_WORK_SPACE(this, func, F1, iFirstMomIndexNode1 + 1);
        CHECK_NUM_COLS_WORK_SPACE(this, func, M1, iFirstMomIndexNode1 + 4);
        CHECK_NUM_COLS_WORK_SPACE(this, func, F2, iFirstMomIndexNode2 + 1);
        CHECK_NUM_COLS_WORK_SPACE(this, func, M2, iFirstMomIndexNode2 + 4);

        WorkVec.AddItem(iFirstMomIndexNode1 + 1, F1);
        WorkVec.AddItem(iFirstMomIndexNode1 + 4, M1);
        WorkVec.AddItem(iFirstMomIndexNode2 + 1, F2);
        WorkVec.AddItem(iFirstMomIndexNode2 + 4, M2);
    }

    enum LubricationGroove::Type
    CylindricalMeshAtBearing::ReadLubricationGrooveType(MBDynParser& HP) const
    {
        if (HP.IsKeyWord("at" "shaft")) {
            return LubricationGroove::MOVING;
        } else if (HP.IsKeyWord("at" "bearing")) {
            return LubricationGroove::FIXED;
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"at shaft\" or \"at bearing\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    bool
    CylindricalMeshAtBearing::bGetPrivateData(HydroRootBase::PrivateDataType eType,
                                              doublereal& dPrivData) const
    {
        const Matrix<doublereal, 3, 3>& Rb2 = GetOrientationNode2();

        switch (eType) {
        case HydroRootBase::PD_F1x:
        case HydroRootBase::PD_F1y:
        case HydroRootBase::PD_F1z: {
            const index_type iIndex = eType - HydroRootBase::PD_F1x + 1;
            dPrivData = Dot(Rb2.GetCol(iIndex), oReaction.F1_R2);
            return true;
        }
        case HydroRootBase::PD_M1x:
        case HydroRootBase::PD_M1y:
        case HydroRootBase::PD_M1z: {
            const index_type iIndex = eType - HydroRootBase::PD_M1x + 1;
            dPrivData = Dot(Rb2.GetCol(iIndex), oReaction.M1_R2);
            return true;
        }
        case HydroRootBase::PD_F2x:
        case HydroRootBase::PD_F2y:
        case HydroRootBase::PD_F2z: {
            const index_type iIndex = eType - HydroRootBase::PD_F2x + 1;
            dPrivData = Dot(Rb2.GetCol(iIndex), oReaction.F2_R2);
            return true;
        }
        case HydroRootBase::PD_M2x:
        case HydroRootBase::PD_M2y:
        case HydroRootBase::PD_M2z: {
            const index_type iIndex = eType - HydroRootBase::PD_M2x + 1;
            dPrivData = Dot(Rb2.GetCol(iIndex), oReaction.M2_R2);
            return true;
        }
        default:
            return false;
        }
    }

    void
    CylindricalMeshAtBearing::AddReactionForce(const Vector<doublereal, 2>& x,
                                               const Vector<doublereal, 3>& v,
                                               const Matrix<doublereal, 3, 3>& Rt,
                                               const Vector<doublereal, 3>& dF_0_Rt,
                                               const Vector<doublereal, 3>& dF_h_Rt,
                                               const Vector<doublereal, 2>& dM_h_Rt,
                                               doublereal dCoef,
                                               enum FunctionCall func,
                                               LocalDofMap* pDofMap)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        AddReactionForce(x, v, Rt, dF_0_Rt, dF_h_Rt, dM_h_Rt, dCoef, func, pDofMap, oReaction);

#if CREATE_PROFILE == 1
        profile.dtAddForce[PROF_RES] += mbdyn_clock_time() - start;
#endif
    }

    void
    CylindricalMeshAtBearing::AddReactionForce(const Vector<doublereal, 2>& x,
                                               const Vector<doublereal, 3>& v,
                                               const Matrix<doublereal, 3, 3>& Rt,
                                               const Vector<Gradient<0>, 3>& dF_0_Rt,
                                               const Vector<Gradient<0>, 3>& dF_h_Rt,
                                               const Vector<Gradient<0>, 2>& dM_h_Rt,
                                               doublereal dCoef,
                                               enum FunctionCall func,
                                               LocalDofMap* pDofMap)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        AddReactionForce(x, v, Rt, dF_0_Rt, dF_h_Rt, dM_h_Rt, dCoef, func, pDofMap, oReaction_grad);

#if CREATE_PROFILE == 1
        profile.dtAddForce[PROF_JAC] += mbdyn_clock_time() - start;
#endif
    }

    template <typename T>
    void CylindricalMeshAtBearing::ReactionForce<T>::Reset()
    {
        F1_R2.Reset();
        M1_R2.Reset();
        F2_R2.Reset();
        M2_R2.Reset();
    }

    template <typename T>
    inline void
    CylindricalMeshAtBearing::AddReactionForce(const Vector<doublereal, 2>& x2,
                                               const Vector<doublereal, 3>& Rb2_v2,
                                               const Matrix<doublereal, 3, 3>& Rbt2,
                                               const Vector<T, 3>& dF2_0_Rt2,
                                               const Vector<T, 3>& dF1_h_Rt2,
                                               const Vector<T, 2>& dM1_h_Rt2,
                                               doublereal dCoef,
                                               enum FunctionCall func,
                                               LocalDofMap* pDofMap,
                                               ReactionForce<T>& oReact)
    {
        typedef Vector<T, 3> VVec3;

        const VVec3 dF1_R2 = Rbt2 * dF1_h_Rt2;
        const VVec3 dM1_h_R2 = Rbt2.GetCol(1) * dM1_h_Rt2(1) + Rbt2.GetCol(3) * dM1_h_Rt2(2);

        oReact.F1_R2 += dF1_R2;
        oReact.M1_R2 += Cross(Rb2_v2, dF1_R2) + dM1_h_R2;

        const VVec3 dF2_R2 = Rbt2 * dF2_0_Rt2;

        oReact.F2_R2 += dF2_R2;
        oReact.M2_R2 += Cross(Rb2_v2, dF2_R2) - dM1_h_R2;
    }

    const CylindricalMeshAtBearing::ReactionForce<doublereal>&
    CylindricalMeshAtBearing::GetReactionForce(const doublereal& dummy) const
    {
        return oReaction;
    }

    const CylindricalMeshAtBearing::ReactionForce<Gradient<0> >&
    CylindricalMeshAtBearing::GetReactionForce(const Gradient<0>& dummy) const
    {
        return oReaction_grad;
    }

    template <typename T>
    CylindricalMeshAtBearing::Boundary<T>::Boundary(const CylindricalMeshAtBearing& rParent)
        :rParent(rParent)
    {

    }

    template <typename T> void
    CylindricalMeshAtBearing::Boundary<T>::Initialize()
    {
        const CVec3& o1 = rParent.GetOffsetNode1();
        const CMat3x3& Rb1 = rParent.GetOrientationNode1();

        Rb1T_o1 = Transpose(Rb1) * o1;
    }

    template <typename T> void
    CylindricalMeshAtBearing::Boundary<T>::Update(doublereal dCoef,
                                                  enum FunctionCall func,
                                                  LocalDofMap* pDofMap)
    {
        const StructNode* const pNode1 = rParent.pGetNode1();
        const StructNode* const pNode2 = rParent.pGetNode2();

        pNode1->GetXCurr(X1, dCoef, func, pDofMap);
        pNode1->GetRCurr(R1, dCoef, func, pDofMap);
        pNode1->GetVCurr(X1P, dCoef, func, pDofMap);
        pNode1->GetWCurr(omega1, dCoef, func, pDofMap);

        pNode2->GetXCurr(X2, dCoef, func, pDofMap);
        pNode2->GetRCurr(R2, dCoef, func, pDofMap);
        pNode2->GetVCurr(X2P, dCoef, func, pDofMap);
        pNode2->GetWCurr(omega2, dCoef, func, pDofMap);

        const CMat3x3& Rb1 = rParent.GetOrientationNode1();

        Rb1T_R1T = Transpose(Rb1) * Transpose(R1);
    }

    template <typename T> void
    CylindricalMeshAtBearing::Boundary<T>::GetBoundaryConditions(HydroNode* pNode,
                                                                 T& h,
                                                                 T& dh_dt,
                                                                 Vector<T, 2>& U1,
                                                                 Vector<T, 2>& U2,
                                                                 Vector<T, 2>& U,
                                                                 doublereal dCoef,
                                                                 enum FunctionCall func,
                                                                 LocalDofMap* pDofMap) const
    {
        const CVec3& o1 = rParent.GetOffsetNode1();
        const CVec3& o2 = rParent.GetOffsetNode2();
        const CMat3x3& Rb1 = rParent.GetOrientationNode1();
        const CVec3& Rb2_v2 = pNode->GetPosition3D();

        const VVec3 a3 = R2 * CVec3(o2 + Rb2_v2);
        const VVec3 a1 = X2 + a3 - X1;
        const VVec3 b = Rb1T_R1T * a1 - Rb1T_o1;
        const doublereal r = rParent.dGetShaftRadius();
        const T a4 = b(1) * b(1) + b(2) * b(2);
        const T a0 = sqrt(a4);

        const T cos_Phi1 = b(1) / a0;
        const T sin_Phi1 = b(2) / a0;
        const T& z1 = b(3);

        T Phi1 = atan2(sin_Phi1, cos_Phi1);

        if (Phi1 < 0.) {
            // Attention: pFindShaftPockets assumes that Phi1 is in range 0 ... 2 * pi
            Phi1 += 2 * M_PI;
        }

        const VVec2 x1(r * Phi1, z1);
        const Pocket* const pPocket1 = rParent.pFindShaftPocket(CVec2(x1));

        T Deltay1, dDeltay1_dx1, dDeltay1_dz1;

        if (pPocket1 == nullptr) {
            Deltay1 = 0.;
            dDeltay1_dx1 = 0.;
            dDeltay1_dz1 = 0.;
        } else {
            pPocket1->GetHeight(x1, Deltay1);
            pPocket1->GetHeightDerX(x1, dDeltay1_dx1);
            pPocket1->GetHeightDerZ(x1, dDeltay1_dz1);
        }

        const LubricationGroove* const pMovingGroove = rParent.pFindMovingLubricationGroove(CVec2(x1), pNode->GetNodePhysics());

#if HYDRO_DEBUG > 1
        if (pMovingGroove) {
            HYDRO_TRACE("node(" << pNode->iGetNodeNumber() + 1
                        << ") affected by moving boundary condition at x1("
                        << CVec2(x1) << ")\n");
        }
#endif
        pNode->SetMovingPressBoundCond(pMovingGroove ? pMovingGroove->pGetBoundaryCond() : nullptr);

        T w, dw_dt;

        pNode->GetRadialDeformation(w, dw_dt, dCoef, func, pDofMap);

        const T h0 = a0 - r - Deltay1;
        h = h0 + w;

        const VVec3 v1((r + Deltay1) * cos_Phi1,
                       (r + Deltay1) * sin_Phi1,
                       z1);

        const VVec3 db_dt = Rb1T_R1T * VVec3(X2P + Cross(omega2, a3)
                                             - X1P - Cross(omega1, a1));

        const T dx1_dt = r * (b(1) * db_dt(2) - db_dt(1) * b(2)) / a4;
        const T& dz1_dt = db_dt(3);
        const T dDeltay1_dt = dDeltay1_dx1 * dx1_dt + dDeltay1_dz1 * dz1_dt;

        dh_dt = (b(1) * db_dt(1) + b(2) * db_dt(2)) / a0 - dDeltay1_dt + dw_dt;

        const VVec3 vh(h0 * cos_Phi1,
                       h0 * sin_Phi1,
                       T(0.));

        const VVec3 dP1_dt = X1P + Cross(omega1, VVec3(R1 * VVec3(o1 + Rb1 * v1))) + Cross(omega2, VVec3(R1 * VVec3(Rb1 * vh)));
        const VVec3 dP2_dt = X2P + Cross(omega2, a3);

        const CMat3x3 Rbt2 = pNode->GetTangentCoordSys();

        const VVec3 dP1_dt_R2 = Transpose(R2) * dP1_dt;
        const VVec3 dP2_dt_R2 = Transpose(R2) * dP2_dt;

        U1(1) = Dot(Rbt2.GetCol(1), dP1_dt_R2);
        U1(2) = Dot(Rbt2.GetCol(3), dP1_dt_R2);
        U2(1) = Dot(Rbt2.GetCol(1), dP2_dt_R2);
        U2(2) = Dot(Rbt2.GetCol(3), dP2_dt_R2);

        U = (U1 - U2) * 0.5;
    }

    template <typename T>
    void CylindricalMeshAtBearing::Boundary<T>::GetMovingMeshOffset(Vector<T, 2>& x, LocalDofMap* pDofMap) const
    {
        const CVec3& o1 = rParent.GetOffsetNode1();
        const CVec3& o2 = rParent.GetOffsetNode2();
        const CMat3x3& Rb1 = rParent.GetOrientationNode1();
        const CMat3x3& Rb2 = rParent.GetOrientationNode2();
        const doublereal r = rParent.dGetMeshRadius();

        const VVec3 R1_Rb1_e1 = R1 * Rb1.GetCol(1);

        const T dx = r * atan2(Dot(R2 * Rb2.GetCol(2), R1_Rb1_e1),
                         Dot(R2 * Rb2.GetCol(1), R1_Rb1_e1));

        const T dz = Dot(Rb2.GetCol(3), Transpose(R2) * VVec3(X1 + R1 * o1 - X2) - o2);

        Copy(x(1), dx, pDofMap); // FIXME: Which class is responsible for copying Gradient's
        Copy(x(2), dz, pDofMap);
    }

    void CylindricalMeshAtBearing::GetMovingMeshOffset(Vector<doublereal, 2>& x, LocalDofMap* pDofMap) const
    {
        oBound.GetMovingMeshOffset(x, pDofMap);
    }

    void CylindricalMeshAtBearing::GetMovingMeshOffset(Vector<Gradient<0>, 2>& x, LocalDofMap* pDofMap) const
    {
        oBound_grad.GetMovingMeshOffset(x, pDofMap);
    }

    doublereal CylindricalMeshAtBearing::dGetMeshRadius() const
    {
        return dGetBearingRadius();
    }

    const Matrix<doublereal, 3, 3>&
    CylindricalMeshAtBearing::GetOrientationMeshNode() const
    {
        return GetOrientationNode2();
    }

    const Pocket*
    CylindricalMeshAtBearing::pFindMeshPocket(const Vector<doublereal, 2>& x) const
    {
        return pFindBearingPocket(x);
    }

    void Material::ParseInput(integer iIndex, MBDynParser& HP, const HydroRootElement* pRoot)
    {
        std::ostringstream os;

        os << "E" << iIndex;

        if (!HP.IsKeyWord(os.str().c_str())) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"" << os.str() << "\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal E = HP.GetReal();

        os.str("");

        os << "nu" << iIndex;

        if (!HP.IsKeyWord(os.str().c_str())) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"" << os.str() << "\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal nu = HP.GetReal();

        *this = Material(E, nu);
    }

    ContactModel::ContactModel(HydroMesh* pMesh)
        :pMesh(pMesh)
    {

    }

    ContactModel::~ContactModel()
    {

    }

    void ContactModel::ParseInput(MBDynParser& HP)
    {
        if (HP.IsKeyWord("E")) {
            rgMaterials[0] = Material(HP.GetReal(), 0.);
            rgMaterials[1] = Material::Rigid();
        } else {
            for (index_type i = 0; i < iNumMaterials; ++i) {
                rgMaterials[i].ParseInput(i + 1, HP, pMesh->pGetParent());
            }
        }
    }

    GreenwoodTrippCM::GreenwoodTrippCM(HydroMesh* pMesh)
        :ContactModel(pMesh),
         k(0.),
         sigmaDelta(0.),
         H0(0.), Hoffset(0.), a0(0.), a1(0.)
    {

    }

    GreenwoodTrippCM::~GreenwoodTrippCM()
    {

    }

    void GreenwoodTrippCM::ParseInput(MBDynParser& HP)
    {
        ContactModel::ParseInput(HP);

        const doublereal E = GetMaterial(0).dGetReducedModulus(GetMaterial(1));

        enum
        {
            UNIT_m,
            UNIT_um
        } eUnits = UNIT_um;

        if (HP.IsKeyWord("base" "unit")) {
            if (HP.IsKeyWord("meters")) {
                eUnits = UNIT_m;
            } else if (HP.IsKeyWord("micro" "meters")) {
                eUnits = UNIT_um;
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetMesh()->pGetParent()->GetLabel()
                            << ") keyword \"meters\" or \"micro meters\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (!HP.IsKeyWord("M0")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"M0\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        doublereal M0 = HP.GetReal(); // assume unit (um)^2

        if (eUnits == UNIT_um) {
            M0 *= 1e-12;
        }

        if (!HP.IsKeyWord("M2")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"M2\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal M2 = HP.GetReal();     // dimension less

        if (!HP.IsKeyWord("M4")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"M4\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        doublereal M4 = HP.GetReal(); // assume unit 1/(um)^2

        if (eUnits == UNIT_um) {
            M4 *= 1e12;
        }

        if (HP.IsKeyWord("H0") || HP.IsKeyWord("linearize")) {
            if (HP.IsKeyWord("never")) {
                H0 = -std::numeric_limits<doublereal>::max();
            } else {
                H0 = HP.GetReal();
            }
        }

        if (HP.IsKeyWord("offset")) {
            Hoffset = HP.GetReal();
        }

        a0 = pow(4 - H0, 6.804);
        a1 = -6.804 * pow(4 - H0, 5.804);

        sigmaDelta = sqrt(M0);

        const doublereal etaasp = M4 / (M2 * 6. * M_PI * sqrt(3.));

        const doublereal betaasp = 3. * sqrt(M_PI) / (8. * sqrt(M4));

        k = 16. * sqrt(2) / 15. * M_PI * std::pow(etaasp * betaasp * sigmaDelta, 2) * E * sqrt(sigmaDelta / betaasp) * 4.4086e-5;
    }

    bool GreenwoodTrippCM::GetContactPressure(const doublereal h, doublereal& pasp) const
    {
        return ContactPressureTpl(h, pasp);
    }

    bool GreenwoodTrippCM::GetContactPressure(const Gradient<0>& h, Gradient<0>& pasp) const
    {
        return ContactPressureTpl(h, pasp);
    }

    template <typename T>
    bool GreenwoodTrippCM::ContactPressureTpl(const T& h, T& pasp) const
    {
        const T H = h / sigmaDelta - Hoffset; // assume sigmaDelta in SI units [m]

        if (H <= 4) {
            if (H >= H0) {
                pasp = k * pow(4. - H, 6.804);
            } else {
                pasp = k * (a0 + a1 * (H - H0));
            }

            return true;
        } else {
            pasp = 0.;
            return false;
        }
    }

    PenaltyCM::PenaltyCM(HydroMesh* pMesh, doublereal href)
        :ContactModel(pMesh),
         a(0.),
         b(0.),
         c(0.),
         h0(0.),
         h1(0.),
         href(href)
    {

    }

    PenaltyCM::~PenaltyCM()
    {

    }

    void PenaltyCM::ParseInput(MBDynParser& HP)
    {
        ContactModel::ParseInput(HP);

        doublereal dDefScale = 1;

        if (HP.IsKeyWord("reference" "gap" "height")) {
            dDefScale = HP.GetReal();
        } else {
            dDefScale = href;
        }

        const doublereal Ered = GetMaterial(0).dGetReducedModulus(GetMaterial(1));

        if (HP.IsKeyWord("compliance" "factor")) {
            b = Ered / (HP.GetReal() * dDefScale);
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"compliance factor\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (b < 0.) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): \"stiffness coefficient\" must not be negative "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("transition" "region")) {
            h0 = HP.GetReal() * dDefScale;

            if (h0 < 0.) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetMesh()->pGetParent()->GetLabel()
                            << "): \"transition region\" must not be negative "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (HP.IsKeyWord("offset")) {
            h1 = HP.GetReal() * dDefScale;
        }

        h0 += h1;
    }

    bool PenaltyCM::GetContactPressure(const doublereal h, doublereal& pasp) const
    {
        return ContactPressureTpl(h, pasp);
    }

    bool PenaltyCM::GetContactPressure(const Gradient<0>& h, Gradient<0>& pasp) const
    {
        return ContactPressureTpl(h, pasp);
    }

    template <typename T>
    bool PenaltyCM::ContactPressureTpl(const T& h, T& pasp) const
    {
        if (h <= h0) {
            if (h <= h1) {
                pasp = b * (h1 - h) + c;
            } else {
                const T dh = h0 - h;
                pasp = a * dh * dh;
            }

            return true;
        } else {
            pasp = 0.;
            return false;
        }
    }

    FrictionModel::FrictionModel(HydroMesh* pMesh)
        :pMesh(pMesh)
    {

    }

    FrictionModel::~FrictionModel()
    {

    }

    void FrictionModel::AfterPredict(const VectorHandler& X, const VectorHandler& XP)
    {

    }

    void FrictionModel::AfterConvergence(const VectorHandler& X, const VectorHandler& XP)
    {

    }

    CoulombFriction::CoulombFriction(HydroMesh* pMesh)
        :FrictionModel(pMesh),
         mu(0.),
         signumDeltaU(0.)
    {

    }

    CoulombFriction::~CoulombFriction()
    {

    }

    void CoulombFriction::ParseInput(MBDynParser& HP)
    {
        if (HP.IsKeyWord("mu")
            || HP.IsKeyWord("coulomb" "friction" "coefficient")) {
            mu = HP.GetReal();
        }

        if (HP.IsKeyWord("signum" "delta" "u")
            || HP.IsKeyWord("sliding" "velocity" "threshold")) {
            signumDeltaU = HP.GetReal();
        } else {
            silent_cerr("keyword \"sliding velocity threshold\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    void CoulombFriction::GetFrictionForce(const doublereal h, const Vector<doublereal, 2>& U, doublereal p, Vector<doublereal, 2>& tau)
    {
        GetFrictionForceTpl(h, U, p, tau);
    }

    void CoulombFriction::GetFrictionForce(const Gradient<0>& h, const Vector<Gradient<0>, 2>& U, const Gradient<0>& p, Vector<Gradient<0>, 2>& tau)
    {
        GetFrictionForceTpl(h, U, p, tau);
    }

    std::unique_ptr<FrictionModel> CoulombFriction::Clone() const
    {
        return std::unique_ptr<FrictionModel>{new CoulombFriction(*this)};
    }

    template <typename T>
    void CoulombFriction::GetFrictionForceTpl(const T& h, const Vector<T, 2>& U, const T& p, Vector<T, 2>& tau)
    {
        Vector<doublereal, 2> u(U); // Do not consider U in the Jacobian matrix for numerical reasons!

        const doublereal norm_u = sqrt(Dot(u, u));

        if (norm_u != 0.) {
            u /= norm_u;
        }

        if (signumDeltaU != 0.) {
            u *= tanh(2 * M_PI * norm_u / signumDeltaU);
        }

        tau = u * (mu * p);
    }

    LugreFriction::LugreFriction(HydroMesh* pMesh)
        :FrictionModel(pMesh),
         beta(1.),
         vs(0.),
         gamma(1.)
    {
        tCurr = tPrev = pGetMesh()->pGetParent()->dGetTime();
    }

    LugreFriction::~LugreFriction()
    {

    }

    void LugreFriction::ParseInput(MBDynParser& HP)
    {
        if (HP.IsKeyWord("method")) {
            if (HP.IsKeyWord("explicit" "euler")) {
                beta = 0.;
            } else if (HP.IsKeyWord("implicit" "euler")) {
                beta = 1.;
            } else if (HP.IsKeyWord("trapezoidal" "rule")) {
                beta = 0.5;
            } else if (HP.IsKeyWord("custom")) {
                beta = HP.GetReal();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetMesh()->pGetParent()->GetLabel()
                            << "): keyword \"explicit euler\", "
                            "\"implicit euler\", "
                            "\"trapezoidal rule\" or \"custom\" "
                            "expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (!(HP.IsKeyWord("coulomb" "friction" "coefficient")
              || HP.IsKeyWord("coulomb" "friction" "coefficient" "x"))) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"coulomb friction coefficient\""
                        " or \"coulomb friction coefficient x\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal mukx = HP.GetReal();

        doublereal muky;

        if (HP.IsKeyWord("coulomb" "friction" "coefficient" "y")) {
            muky = HP.GetReal();
        } else {
            muky = mukx;
        }

        doublereal musx, musy;

        if (HP.IsKeyWord("static" "friction" "coefficient")
            || HP.IsKeyWord("static" "friction" "coefficient" "x")) {
            musx = HP.GetReal();

            if (HP.IsKeyWord("static" "friction" "coefficient" "y")) {
                musy = HP.GetReal();
            } else {
                musy = musx;
            }
        } else {
            musx = mukx;
            musy = muky;
        }

        if (HP.IsKeyWord("sliding" "velocity" "coefficient")) {
            vs = HP.GetReal();
        } else {
            vs = 1.;
        }

        if (HP.IsKeyWord("sliding" "velocity" "exponent")) {
            gamma = HP.GetReal();
        } else {
            gamma = 1.;
        }

        if (!(HP.IsKeyWord("micro" "slip" "stiffness") || HP.IsKeyWord("micro" "slip" "stiffness" "x"))) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): keyword \"micro slip stiffness\" or \"micro slip stiffness x\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal sigma0x = HP.GetReal();

        doublereal sigma0y;

        if (HP.IsKeyWord("micro" "slip" "stiffness" "y")) {
            sigma0y = HP.GetReal();
        } else {
            sigma0y = sigma0x;
        }

        doublereal sigma1x, sigma1y;

        if (HP.IsKeyWord("micro" "slip" "damping") || HP.IsKeyWord("micro" "slip" "damping" "x")) {
            sigma1x = HP.GetReal();

            if (HP.IsKeyWord("micro" "slip" "damping" "y")) {
                sigma1y = HP.GetReal();
            } else {
                sigma1y = sigma1x;
            }
        } else {
            sigma1x = 0.;
            sigma1y = 0.;
        }

        Mk(1, 1) = mukx;
        Mk(2, 2) = muky;

        Mk2 = Mk * Mk;

        Ms(1, 1) = musx;
        Ms(2, 2) = musy;

        Ms2 = Ms * Ms;

        sigma0(1, 1) = sigma0x;
        sigma0(2, 2) = sigma0y;

        sigma1(1, 1) = sigma1x;
        sigma1(2, 2) = sigma1y;

        invMk2_sigma0 = Inv(Mk2) * sigma0;
    }

    void LugreFriction::GetFrictionForce(const doublereal h, const Vector<doublereal, 2>& U, doublereal p, Vector<doublereal, 2>& tau)
    {
        GetFrictionForceTpl(h, U, p, tau);
    }

    void LugreFriction::GetFrictionForce(const Gradient<0>& h, const Vector<Gradient<0>, 2>& U, const Gradient<0>& p, Vector<Gradient<0>, 2>& tau)
    {
        GetFrictionForceTpl(h, U, p, tau);
    }

    std::unique_ptr<FrictionModel> LugreFriction::Clone() const
    {
        return std::unique_ptr<FrictionModel>{new LugreFriction(*this)};
    }

    void LugreFriction::AfterConvergence(const VectorHandler& X,
                                         const VectorHandler& XP)
    {
        tPrev = tCurr;
        zPrev = zCurr;
        zPPrev = zPCurr;
    }

    template <typename T>
    void LugreFriction::GetFrictionForceTpl(const T& h, const Vector<T, 2>& U, const T& p, Vector<T, 2>& tau)
    {
        typedef Matrix<doublereal, 2, 2> CMat2x2;
        typedef Vector<doublereal, 2> CVec2;
        typedef Matrix<T, 2, 2> VMat2x2;
        typedef Vector<T, 2> VVec2;

        const VVec2 Ueff = U * doublereal(p > 0.);

        const T norm_Ueff = Dot(Ueff, Ueff);

        T kappa;

        if (norm_Ueff == 0.) {
            kappa = 0.;
        } else {
            const VVec2 Mk_U = Mk * Ueff;
            const VVec2 Ms_U = Ms * Ueff;
            const VVec2 Mk2_U = Mk2 * Ueff;
            const VVec2 Ms2_U = Ms2 * Ueff;
            const T norm_Mk2_U = sqrt(Dot(Mk2_U, Mk2_U));
            const T a0 = norm_Mk2_U / sqrt(Dot(Mk_U, Mk_U));
            const T a1 = sqrt(Dot(Ms2_U, Ms2_U)) / sqrt(Dot(Ms_U, Ms_U));
            const T g = a0 + (a1 - a0) * exp(-pow(sqrt(norm_Ueff) / vs, gamma));

            kappa = norm_Mk2_U / g;
        }

        tCurr = pGetMesh()->pGetParent()->dGetTime();

        const doublereal dt = tCurr - tPrev;
        const VMat2x2 A = invMk2_sigma0 * kappa;
        const VMat2x2 B = A * (beta * dt) + CMat2x2(1., 0., 0., 1);

        const Vector<T, 2> zP = Inv(B) * VVec2(Ueff - A * CVec2(zPrev + zPPrev * ((1 - beta) * dt)));
        const Vector<T, 2> z = zPrev + (zP * beta + zPPrev * (1 - beta)) * dt;

        SaveStictionState(z, zP);

        tau = (sigma0 * z + sigma1 * zP) * p;
    }

    void LugreFriction::SaveStictionState(const Vector<doublereal, 2>& z, const Vector<doublereal, 2>& zP)
    {
        zCurr = z;
        zPCurr = zP;

        HYDRO_TRACE( "dt=" << tCurr - tPrev << std::endl);
        HYDRO_TRACE( "z=" << zCurr << std::endl);
        HYDRO_TRACE( "zP=" << zPCurr << std::endl);
    }

    void LugreFriction::SaveStictionState(const Vector<Gradient<0>, 2>&, const Vector<Gradient<0>, 2>&)
    {
        // Do Nothing
    }

    HydroElement::HydroElement(HydroMesh* pMeshArg, ElementType eType)
        :pMesh(pMeshArg),
         pFluid(pMeshArg->pGetParent()->pGetFluid()),
         eType(eType)
    {

    }

    HydroElement::~HydroElement()
    {

    }

    integer HydroElement::iGetNumColsWorkSpace(grad::FunctionCall eFunc, index_type) const
    {
        const HydroMesh* const pMesh = pGetMesh();
        const BearingGeometry* const pGeometry = pMesh->pGetGeometry();
        const ComplianceModel* pComplianceModel = pMesh->pGetComplianceModel();

        integer iNumCols = pGeometry->iGetNumColsWorkSpace(eFunc);

        if (pComplianceModel) {
            iNumCols += pComplianceModel->iGetNumColsWorkSpace(eFunc, iGetNumNodes());
        }

        for (index_type i = 0; i < iGetNumNodes(); ++i) {
            const HydroNode* pNode = pGetNode(i);
            const ThermoHydrNode* pThermalNode = pNode->pGetThermalNode();

            iNumCols += pNode->iGetNumColsWorkSpace(eFunc);

            if (pThermalNode) {
                iNumCols += pThermalNode->iGetNumColsWorkSpace(eFunc);
            }
        }

        return iNumCols;
    }

    void HydroElement::Print(std::ostream& os) const
    {
        os << iGetNumNodes() << '(';

        for (int i = 0; i < iGetNumNodes(); ++i) {
            os << pGetNode(i)->iGetNodeNumber() + 1 << ' ';
        }

        os << ")\n";
    }

    void HydroElement::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
        // NO_OP
    }

    void HydroElement::AfterConvergence(const VectorHandler& X,
                                        const VectorHandler& XP)
    {
        // NO_OP
    }

    LinFD5Elem::LinFD5Elem(HydroMesh* pMesh, ElementType eType)
        :HydroElement(pMesh, eType), dx(0.), dz(0.)
    {
        std::fill(rgHydroNodes.begin(), rgHydroNodes.end(), nullptr);
        std::fill(rgFluxNodes.begin(), rgFluxNodes.end(), nullptr);
    }

    LinFD5Elem::~LinFD5Elem()
    {

    }

    int LinFD5Elem::iGetNumNodes() const
    {
        return iNumNodes;
    }

    void LinFD5Elem::SetNode(int iNode, HydroNode* pNode)
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgHydroNodes[iNode] == nullptr);
        HYDRO_ASSERT(pNode != nullptr);

        rgHydroNodes[iNode] = pNode;
    }

    HydroNode* LinFD5Elem::pGetNode(int iNode) const
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgHydroNodes[iNode] != nullptr);

        return rgHydroNodes[iNode];
    }

    void LinFD5Elem::SetFluxNode(int iNode, FluxNode* pFluxNode)
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumFluxNodes);
        HYDRO_ASSERT(rgFluxNodes[iNode] == nullptr);
        HYDRO_ASSERT(pFluxNode != nullptr);

        pFluxNode->RequestNodeData(FluxNode::ND_HYDRAULIC);

        rgFluxNodes[iNode] = pFluxNode;
    }

    void LinFD5Elem::Initialize()
    {
        for (int i = 0; i < iNumNodes; ++i) {
            HYDRO_ASSERT(rgHydroNodes[i] != nullptr);

            x[i] = rgHydroNodes[i]->GetPosition2D();
        }

        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        dx = 0.5 * (pGeometry->dGetNodeDistance2D(rgHydroNodes[iNodeWest], rgHydroNodes[iNodeCenter], 1)
                    + pGeometry->dGetNodeDistance2D(rgHydroNodes[iNodeCenter], rgHydroNodes[iNodeEast], 1));

        dz = 0.5 * (pGeometry->dGetNodeDistance2D(rgHydroNodes[iNodeSouth], rgHydroNodes[iNodeCenter], 2)
                    + pGeometry->dGetNodeDistance2D(rgHydroNodes[iNodeCenter], rgHydroNodes[iNodeNorth], 2));

        dA = dx * dz;

        HYDRO_ASSERT(rgFluxNodes[iNodeFlxWest]->pGetNode(FluxNode::iNodeDown) == rgHydroNodes[iNodeWest]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlxWest]->pGetNode(FluxNode::iNodeUp) == rgHydroNodes[iNodeCenter]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlxEast]->pGetNode(FluxNode::iNodeDown) == rgHydroNodes[iNodeCenter]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlxEast]->pGetNode(FluxNode::iNodeUp) == rgHydroNodes[iNodeEast]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlzSouth]->pGetNode(FluxNode::iNodeDown) == rgHydroNodes[iNodeSouth]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlzSouth]->pGetNode(FluxNode::iNodeUp) == rgHydroNodes[iNodeCenter]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlzNorth]->pGetNode(FluxNode::iNodeDown) == rgHydroNodes[iNodeCenter]);
        HYDRO_ASSERT(rgFluxNodes[iNodeFlzNorth]->pGetNode(FluxNode::iNodeUp) == rgHydroNodes[iNodeNorth]);

        HYDRO_ASSERT(x[iNodeEast](2) == x[iNodeWest](2));
        HYDRO_ASSERT(x[iNodeCenter](2) == x[iNodeWest](2));
        HYDRO_ASSERT(x[iNodeNorth](1) == x[iNodeSouth](1));
        HYDRO_ASSERT(x[iNodeCenter](1) == x[iNodeSouth](1));
        HYDRO_ASSERT(x[iNodeNorth](2) > x[iNodeSouth](2));
        HYDRO_ASSERT(x[iNodeNorth](2) > x[iNodeCenter](2));
        HYDRO_ASSERT(x[iNodeCenter](2) > x[iNodeSouth](2));
        HYDRO_ASSERT(x[iNodeEast](1) > x[iNodeWest](1));
        HYDRO_ASSERT(x[iNodeEast](1) > x[iNodeCenter](1));
        HYDRO_ASSERT(x[iNodeCenter](1) > x[iNodeWest](1));
    }

    LinFD4Elem::LinFD4Elem(HydroMesh* pMesh, ElementType eType)
        :HydroElement(pMesh, eType),
         dx(0.),
         dz(0.),
         dA(0.)
    {
        std::fill(rgHydroNodes.begin(), rgHydroNodes.end(), nullptr);
    }

    LinFD4Elem::~LinFD4Elem()
    {

    }

    int LinFD4Elem::iGetNumNodes() const
    {
        return iNumNodes;
    }

    void LinFD4Elem::SetNode(int iNode, HydroNode* pNode)
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgHydroNodes[iNode] == nullptr);

        rgHydroNodes[iNode] = pNode;
    }

    HydroNode* LinFD4Elem::pGetNode(int iNode) const
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgHydroNodes[iNode] != nullptr);

        return rgHydroNodes[iNode];
    }

    void LinFD4Elem::Initialize()
    {
#if HYDRO_DEBUG > 0
        for (int i = 0; i < iNumNodes; ++i) {
            HYDRO_ASSERT(rgHydroNodes[i] != nullptr);
        }

        const doublereal x1 = rgHydroNodes[iNode1NE]->GetPosition2D()(1);
        const doublereal x2 = rgHydroNodes[iNode2NW]->GetPosition2D()(1);
        const doublereal x3 = rgHydroNodes[iNode3SW]->GetPosition2D()(1);
        const doublereal x4 = rgHydroNodes[iNode4SE]->GetPosition2D()(1);
        const doublereal z1 = rgHydroNodes[iNode1NE]->GetPosition2D()(2);
        const doublereal z2 = rgHydroNodes[iNode2NW]->GetPosition2D()(2);
        const doublereal z3 = rgHydroNodes[iNode3SW]->GetPosition2D()(2);
        const doublereal z4 = rgHydroNodes[iNode4SE]->GetPosition2D()(2);

        // preconditions
        HYDRO_ASSERT(z1 == z2);
        HYDRO_ASSERT(z3 == z4);
        HYDRO_ASSERT(x1 == x4);
        HYDRO_ASSERT(x2 == x3);
#endif
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        dx = 0.5 * (pGeometry->dGetNodeDistance2D(rgHydroNodes[iNode2NW], rgHydroNodes[iNode1NE], 1)
                    + pGeometry->dGetNodeDistance2D(rgHydroNodes[iNode3SW], rgHydroNodes[iNode4SE], 1));

        dz = 0.5 * (pGeometry->dGetNodeDistance2D(rgHydroNodes[iNode4SE], rgHydroNodes[iNode1NE], 2)
                    + pGeometry->dGetNodeDistance2D(rgHydroNodes[iNode3SW], rgHydroNodes[iNode2NW], 2));

        dA = dx * dz;
    }

    LinFD5ReynoldsElem::LinFD5ReynoldsElem(HydroMesh* pMesh)
        :LinFD5Elem(pMesh, REYNOLDS_ELEM)
    {

#if CREATE_PROFILE == 1
        profile.pLastElem = this;
#endif
    }

    LinFD5ReynoldsElem::~LinFD5ReynoldsElem()
    {
#if CREATE_PROFILE == 1
        if (profile.pLastElem == this) {
            for (int i = PROF_RES; i <= PROF_JAC; ++i) {
                std::cerr << "LinFDReynoldsElem::dtAss[" << i << "]=" << profile.dtAss[i] << std::endl;
            }
        }
#endif
    }

    void LinFD5ReynoldsElem::AssRes(SubVectorHandler& WorkVec,
                                    doublereal dCoef,
                                    const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr,
                                    enum GradientAssVecBase::mode_t mode)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);

#if CREATE_PROFILE == 1
        profile.dtAss[PROF_RES] += mbdyn_clock_time() - start;
#endif
    }

    void LinFD5ReynoldsElem::AssJac(SparseSubMatrixHandler& WorkMat,
                                    doublereal dCoef,
                                    const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr,
                                    enum GradientAssVecBase::mode_t mode)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);

#if CREATE_PROFILE == 1
        profile.dtAss[PROF_JAC] += mbdyn_clock_time() - start;
#endif
    }

    void LinFD5ReynoldsElem::InitialAssRes(SubVectorHandler& WorkVec,
                                           const VectorHandler& XCurr,
                                           enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void LinFD5ReynoldsElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                           const VectorHandler& XCurr,
                                           enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }

    void LinFD5ReynoldsElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = 1;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());
    }

    template <typename T>
    void LinFD5ReynoldsElem::AssRes(GradientAssVec<T>& WorkVec,
                                    doublereal dCoef,
                                    const GradientVectorHandler<T>& XCurr,
                                    const GradientVectorHandler<T>& XPrimeCurr,
                                    enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename T>
    void LinFD5ReynoldsElem::InitialAssRes(GradientAssVec<T>& WorkVec,
                                           const GradientVectorHandler<T>& XCurr,
                                           enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    template <typename T>
    void LinFD5ReynoldsElem::UnivAssRes(GradientAssVec<T>& WorkVec,
                                        doublereal dCoef,
                                        const GradientVectorHandler<T>& XCurr,
                                        enum FunctionCall func)
    {
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();
        const FluidStateBoundaryCond* const pBoundaryCond = pGeometry->pGetMovingPressBoundCond(rgHydroNodes[iNodeCenter]);

        const integer iFirstIndex = rgHydroNodes[iNodeCenter]->iGetFirstEquationIndex(func);

        const doublereal dEquationScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_REYNOLDS_EQU) / dCoef;

        if (pBoundaryCond == nullptr) {
            T h, dh_dt, rho, drho_dt;

            std::array<T, iNumFluxNodes> mdot;

            for (index_type i = 0; i < iNumFluxNodes; ++i) {
                rgFluxNodes[i]->GetMassFluxDens(mdot[i], &oDofMap);
            }

            rgHydroNodes[iNodeCenter]->GetDensity(rho, &oDofMap, dCoef);
            rgHydroNodes[iNodeCenter]->GetClearance(h, &oDofMap);
            rgHydroNodes[iNodeCenter]->GetClearanceDerTime(dh_dt, &oDofMap);
            pGeometry->GetNonNegativeClearance(h, h, &dh_dt, &dh_dt);

            if (func & grad::REGULAR_FLAG) {
                rgHydroNodes[iNodeCenter]->GetDensityDerTime(drho_dt, &oDofMap, dCoef);
            } else {
                drho_dt = 0.;
            }

            const T Re = ((mdot[iNodeFlxEast] - mdot[iNodeFlxWest]) / dx
                          + (mdot[iNodeFlzNorth] - mdot[iNodeFlzSouth]) / dz
                          + (drho_dt * h + rho * dh_dt))
                * dEquationScale;

            CHECK_NUM_COLS_WORK_SPACE(this, func, Re, iFirstIndex);

            WorkVec.AddItem(iFirstIndex, Re);
        } else {
            const doublereal ppre = pBoundaryCond->dGetPressure();

            T pcenter(0);

            rgHydroNodes[iNodeCenter]->GetPressure(pcenter, &oDofMap, dCoef);

            const T f = (pcenter - ppre) * dEquationScale;

            HYDRO_TRACE("node=" << rgHydroNodes[iNodeCenter]->iGetNodeNumber() + 1
                        << ":" << rgHydroNodes[iNodeCenter]->iGetFirstEquationIndex(func)
                        << ":" << rgHydroNodes[iNodeCenter]->iGetFirstDofIndex(func) << std::endl);
            HYDRO_TRACE("ppre=" << ppre << std::endl);
            HYDRO_TRACE("p=" << pcenter << std::endl);
            HYDRO_TRACE("scale=" << dEquationScale << std::endl);
            HYDRO_TRACE("f(imposed pressure)=" << f << std::endl);

            CHECK_NUM_COLS_WORK_SPACE(this, func, f, iFirstIndex);

            WorkVec.AddItem(iFirstIndex, f);
        }
    }

    LinFD5CouplingElem::LinFD5CouplingElem(HydroMesh* pMesh)
        :LinFD5Elem(pMesh, COUPLING_ELEM)
    {

    }

    LinFD5CouplingElem::~LinFD5CouplingElem()
    {

    }

    void
    LinFD5CouplingElem::AssRes(SubVectorHandler& WorkVec,
                               doublereal dCoef,
                               const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    LinFD5CouplingElem::AssJac(SparseSubMatrixHandler& WorkMat,
                               doublereal dCoef,
                               const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    LinFD5CouplingElem::InitialAssRes(SubVectorHandler& WorkVec,
                                      const VectorHandler& XCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        NO_OP;
    }

    void
    LinFD5CouplingElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                      const VectorHandler& XCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        NO_OP;
    }

    void LinFD5CouplingElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        const bool bDoAssembly = eFunc & grad::REGULAR_FLAG;

        *piNumRows = 1 * bDoAssembly;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes()) * bDoAssembly;
    }

    template <typename T>
    void LinFD5CouplingElem::AssRes(GradientAssVec<T>& WorkVec,
                                    doublereal dCoef,
                                    const GradientVectorHandler<T>& XCurr,
                                    const GradientVectorHandler<T>& XPrimeCurr,
                                    enum FunctionCall func)
    {
        T h, dh_dt, rho, drho_dt;

        std::array<T, iNumFluxNodes> mdot;

        for (index_type i = 0; i < iNumFluxNodes; ++i) {
            rgFluxNodes[i]->GetMassFluxDens(mdot[i], &oDofMap);
        }

        rgHydroNodes[iNodeCenter]->GetDensity(rho, &oDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetClearance(h, &oDofMap);
        rgHydroNodes[iNodeCenter]->GetClearanceDerTime(dh_dt, &oDofMap);
        pGetMesh()->pGetGeometry()->GetNonNegativeClearance(h, h, &dh_dt, &dh_dt);
        rgHydroNodes[iNodeCenter]->GetDensityDerTime(drho_dt, &oDofMap, dCoef);

        const T dm_dt = (mdot[iNodeFlxEast] - mdot[iNodeFlxWest]) * dz
            + (mdot[iNodeFlzNorth] - mdot[iNodeFlzSouth]) * dx
            + (drho_dt * h + rho * dh_dt) * dA;

        const integer iFirstIndex = rgHydroNodes[iNodeCenter]->iGetFirstEquationIndex(func);

        CHECK_NUM_COLS_WORK_SPACE(this, func, dm_dt, iFirstIndex);

        WorkVec.AddItem(iFirstIndex, dm_dt);
    }

    LinFD4FrictionElem::LinFD4FrictionElem(HydroMesh* pMesh)
        :LinFD4Elem(pMesh, FRICTION_ELEM),
         dScaleEnergy(0.)
    {

#if CREATE_PROFILE == 1
        profile.pLastElem = this;
#endif
    }

    LinFD4FrictionElem::~LinFD4FrictionElem()
    {
#if CREATE_PROFILE == 1
        if (profile.pLastElem == this) {
            for (int i = PROF_RES; i <= PROF_JAC; ++i) {
                std::cerr << "LinFDFrictionElem::dtAss[" << i << "]=" << profile.dtAss[i] << std::endl;
            }
        }
#endif
    }

    void
    LinFD4FrictionElem::AssRes(SubVectorHandler& WorkVec,
                               doublereal dCoef,
                               const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               enum GradientAssVecBase::mode_t mode)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);

#if CREATE_PROFILE == 1
        profile.dtAss[PROF_RES] += mbdyn_clock_time() - start;
#endif
    }

    void
    LinFD4FrictionElem::AssJac(SparseSubMatrixHandler& WorkMat,
                               doublereal dCoef,
                               const VectorHandler& XCurr,
                               const VectorHandler& XPrimeCurr,
                               enum GradientAssVecBase::mode_t mode)
    {
#if CREATE_PROFILE == 1
        doublereal start = mbdyn_clock_time();
#endif

        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);

#if CREATE_PROFILE == 1
        profile.dtAss[PROF_JAC] += mbdyn_clock_time() - start;
#endif
    }

    void
    LinFD4FrictionElem::InitialAssRes(SubVectorHandler& WorkVec,
                                      const VectorHandler& XCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    LinFD4FrictionElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                      const VectorHandler& XCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }

    void LinFD4FrictionElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = *piNumCols = 0;
    }

    template <typename T> inline
    void LinFD4FrictionElem::AssRes(GradientAssVec<T>& WorkVec,
                                    doublereal dCoef,
                                    const GradientVectorHandler<T>& XCurr,
                                    const GradientVectorHandler<T>& XPrimeCurr,
                                    enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename T>
    void LinFD4FrictionElem::InitialAssRes(GradientAssVec<T>& WorkVec,
                                           const GradientVectorHandler<T>& XCurr,
                                           enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    template <typename T>
    void LinFD4FrictionElem::UnivAssRes(GradientAssVec<T>& WorkVec,
                                        doublereal dCoef,
                                        const GradientVectorHandler<T>& XCurr,
                                        enum FunctionCall func)
    {
        typedef Vector<T, 2> VVec2;
        typedef Vector<T, 3> VVec3;

        const bool bUpdateFriction = pGetMesh()->pGetParent()->bUpdateFrictionLoss();

        std::array<T, iNumNodes> pi, hi;
        std::array<VVec2, iNumNodes> U1i, U2i;
        T eta(0.), eta_tmp; // Note: viscosity could depend on density and temperature

        for (integer i = 0; i < iNumNodes; ++i) {
            pGetMesh()->GetPressure(rgHydroNodes[i], pi[i], &oDofMap, dCoef);
            rgHydroNodes[i]->GetClearance(hi[i], &oDofMap);
            rgHydroNodes[i]->GetVelocity(U1i[i], U2i[i], &oDofMap);
            rgHydroNodes[i]->GetViscosity(eta_tmp, &oDofMap, dCoef);
            eta += eta_tmp;
        }

        eta /= iNumNodes;

        // use fluid pressure only for pressure gradient
        const T dp_dx = ((pi[iNode1NE] - pi[iNode2NW]) + (pi[iNode4SE] - pi[iNode3SW])) * (0.5 / dx);
        const T dp_dz = ((pi[iNode1NE] - pi[iNode4SE]) + (pi[iNode2NW] - pi[iNode3SW])) * (0.5 / dz);

        VVec2 tauc_0, tauc_0_tmp;
        T dPfc{0.}, Pfc_tmp, pasp_tmp;
        index_type iNumNodesWithContact = 0;

        for (integer i = 0; i < iNumNodes; ++i) {
            if (rgHydroNodes[i]->GetContactPressure(pasp_tmp, &oDofMap)) {
                rgHydroNodes[i]->GetContactStress(tauc_0_tmp, &oDofMap);
                rgHydroNodes[i]->GetContactFrictionLossDens(Pfc_tmp, &oDofMap);
                pi[i] += pasp_tmp;
                tauc_0 += tauc_0_tmp;
                dPfc += Pfc_tmp;
                ++iNumNodesWithContact;
            }
        }

        tauc_0 /= iNumNodes;
        dPfc *= dA / iNumNodes;

        const VVec2 U = (U1i[iNode1NE] - U2i[iNode1NE] + U1i[iNode2NW] - U2i[iNode2NW]
                         + U1i[iNode3SW] - U2i[iNode3SW] + U1i[iNode4SE] - U2i[iNode4SE]) * 0.25;

        const T h = 0.25 * (hi[iNode1NE] + hi[iNode2NW] + hi[iNode3SW] + hi[iNode4SE]);
        const T ptot = 0.25 * (pi[iNode1NE] + pi[iNode2NW] + pi[iNode3SW] + pi[iNode4SE]);

        BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        T hx;

        pGeometry->GetNonNegativeClearance(h, hx);

        const T tau_xy_p_h = 0.5 * hx * dp_dx;
        const T tau_xy_U_h = eta * U(1) / hx;
        const T tau_yz_p_h = 0.5 * hx * dp_dz;
        const T tau_yz_U_h = eta * U(2) / hx;

        const T tau_xy_0 = (-tau_xy_p_h + tau_xy_U_h);
        const T tau_yz_0 = (-tau_yz_p_h + tau_yz_U_h);

        const T tau_xy_h = (tau_xy_p_h + tau_xy_U_h);
        const T tau_yz_h = (tau_yz_p_h + tau_yz_U_h);

        for (int i = 0; i < iNumNodes; ++i) {
            SetStress(rgHydroNodes[i], tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);
        }

        VVec3 dF_0_Rt(tau_xy_0 * dA,
                      -ptot * dA,
                      tau_yz_0 * dA);

        VVec3 dF_h_Rt(-tau_xy_h * dA,
                      ptot * dA,
                      -tau_yz_h * dA);

        const VVec2 dM_h_Rt((dA * dz / 24) * (pi[iNode3SW] - pi[iNode2NW] + pi[iNode4SE] - pi[iNode1NE]),
                            T(0.));

        if (bUpdateFriction) {
            AddFrictionLoss(HydroRootBase::FLUID_FRICTION, U2i, dF_0_Rt(1), dF_0_Rt(3));
            AddFrictionLoss(HydroRootBase::FLUID_FRICTION, U1i, dF_h_Rt(1), dF_h_Rt(3));
        }

        if (iNumNodesWithContact) {
            if (bUpdateFriction) {
                AddFrictionLoss(HydroRootBase::CONTACT_FRICTION, -dPfc);
            }

            static const int iReactionIdx[2] = {1, 3};

            for (int i = 1; i <= 2; ++i) {
                // Note: The contact pressure will be always
                //               normal to the surface of each part
                //           but the angle between the surface
                //               normal vectors might be slightly
                //               different from 180 degrees!

                // Note: In the same way the angle between
                //               the surface tangent vectors at
                //               two different parts might be different
                //               from zero. But it will be a small angle
                //               as long as the relative clearance is small
                //               in case of a cylindrical bearing.

                dF_0_Rt(iReactionIdx[i - 1]) += tauc_0(i) * dA;
                dF_h_Rt(iReactionIdx[i - 1]) -= tauc_0(i) * dA;
            }
        }

        pGeometry->AddReactionForce(xc,
                                    vc,
                                    Rtc,
                                    dF_0_Rt,
                                    dF_h_Rt,
                                    dM_h_Rt,
                                    dCoef,
                                    func,
                                    &oDofMap);
    }

    void LinFD4FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const std::array<Vector<doublereal, 2>, iNumNodes>& Ui, doublereal dTau_xy, doublereal dTau_yz) const
    {
        Vector<doublereal, 2> U;

        for (int i = 0; i < iNumNodes; ++i) {
            U += Ui[i];
        }

        U /= int(iNumNodes);

        const doublereal dPf = U(1) * dTau_xy + U(2) * dTau_yz;

        AddFrictionLoss(type, dPf);
    }

    void LinFD4FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const std::array<Vector<Gradient<0>, 2>, iNumNodes>& Ui, const Gradient<0>& dTau_xy, const Gradient<0>& dTau_yz) const
    {
        NO_OP;
    }

    void LinFD4FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, doublereal dPf) const
    {
        pGetMesh()->pGetParent()->AddFrictionLoss(type, dPf);
    }

    void LinFD4FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const Gradient<0>& dPf) const
    {
        NO_OP;
    }

    void LinFD4FrictionElem::SetStress(HydroNode* pNode, doublereal tau_xy_0, doublereal tau_yz_0, doublereal tau_xy_h, doublereal tau_yz_h)
    {
        pNode->SetStress(tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);
    }

    void LinFD4FrictionElem::SetStress(HydroNode*, const Gradient<0>&, const Gradient<0>&, const Gradient<0>&, const Gradient<0>&)
    {
        NO_OP;
    }

    void LinFD4FrictionElem::Initialize()
    {
        LinFD4Elem::Initialize();

        HYDRO_ASSERT(dA > 0.);
        HYDRO_ASSERT(dx > 0.);
        HYDRO_ASSERT(dz > 0.);

        const doublereal x1 = rgHydroNodes[iNode1NE]->GetPosition2D()(1);
        const doublereal x2 = rgHydroNodes[iNode2NW]->GetPosition2D()(1);
        const doublereal z1 = rgHydroNodes[iNode1NE]->GetPosition2D()(2);
        const doublereal z4 = rgHydroNodes[iNode4SE]->GetPosition2D()(2);

        xc(1) = 0.5 * (x1 + x2);
        xc(2) = 0.5 * (z1 + z4);

        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        pGeometry->GetPosition3D(xc, vc);
        pGeometry->GetTangentCoordSys(xc, Rtc);

        dScaleEnergy = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_ENERGY_EQ);
    }

    LinFD4MassFlowZ::LinFD4MassFlowZ(HydroMesh* pMesh)
        :LinFD4Elem(pMesh, COUPLING_ELEM)
    {
        std::fill(rgFluxNodes.begin(), rgFluxNodes.end(), nullptr);
    }

    LinFD4MassFlowZ::~LinFD4MassFlowZ()
    {
    }

    void
    LinFD4MassFlowZ::AssRes(SubVectorHandler& WorkVec,
                            doublereal dCoef,
                            const VectorHandler& XCurr,
                            const VectorHandler& XPrimeCurr,
                            enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    LinFD4MassFlowZ::AssJac(SparseSubMatrixHandler& WorkMat,
                            doublereal dCoef,
                            const VectorHandler& XCurr,
                            const VectorHandler& XPrimeCurr,
                            enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    LinFD4MassFlowZ::InitialAssRes(SubVectorHandler& WorkVec,
                                   const VectorHandler& XCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {

    }

    void
    LinFD4MassFlowZ::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                   const VectorHandler& XCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {

    }

    void LinFD4MassFlowZ::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        const bool bDoAssembly = eFunc & grad::REGULAR_FLAG;

        *piNumRows = 2 * bDoAssembly;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes()) * bDoAssembly;
    }

    template <typename T> inline
    void LinFD4MassFlowZ::AssRes(GradientAssVec<T>& WorkVec,
                                 doublereal dCoef,
                                 const GradientVectorHandler<T>& XCurr,
                                 const GradientVectorHandler<T>& XPrimeCurr,
                                 enum FunctionCall func)
    {
        static const index_type rgNodeIdxHydro[iNumFluxNodes] = {iNode1NE, iNode2NW};
        T mdotz;

        for (index_type i = 0; i < iNumFluxNodes; ++i) {
            rgFluxNodes[i]->GetMassFluxDens(mdotz, &oDofMap, FluxNode::PRESSURE_FROM_MESH); // correct mass flux even with non mass conserving cavitation model

            HYDRO_TRACE("mdotz(Equation " << rgHydroNodes[rgNodeIdxHydro[i]]->iGetFirstEquationIndex(func) <<  ")=" << mdotz << std::endl);
            HYDRO_ASSERT(rgHydroNodes[rgNodeIdxHydro[i]]->bIsNodeType(Node2D::COUPLED_NODE));

            CHECK_NUM_COLS_WORK_SPACE(this, func, mdotz, rgHydroNodes[rgNodeIdxHydro[i]]->iGetFirstEquationIndex(func));

            WorkVec.AddItem(rgHydroNodes[rgNodeIdxHydro[i]]->iGetFirstEquationIndex(func), -dx / iNumFluxNodes * mdotz);
        }
    }

    void LinFD4MassFlowZ::Initialize()
    {
        LinFD4Elem::Initialize();

        HYDRO_ASSERT(rgHydroNodes[iNode1NE] == rgFluxNodes[iFNodeEast]->pGetNode(dz > 0. ? FluxNode::iNodeUp : FluxNode::iNodeDown));
        HYDRO_ASSERT(rgHydroNodes[iNode2NW] == rgFluxNodes[iFNodeWest]->pGetNode(dz > 0. ? FluxNode::iNodeUp : FluxNode::iNodeDown));
        HYDRO_ASSERT(rgHydroNodes[iNode3SW] == rgFluxNodes[iFNodeWest]->pGetNode(dz > 0. ? FluxNode::iNodeDown : FluxNode::iNodeUp));
        HYDRO_ASSERT(rgHydroNodes[iNode4SE] == rgFluxNodes[iFNodeEast]->pGetNode(dz > 0. ? FluxNode::iNodeDown : FluxNode::iNodeUp));
        HYDRO_ASSERT(fabs(dx) > 0.);
        HYDRO_ASSERT(fabs(dz) > 0.);
        HYDRO_ASSERT(copysign(1., dx) == copysign(1., dz));
    }

    void LinFD4MassFlowZ::SetFluxNode(int iNode, FluxNode* pFluxNode)
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumFluxNodes);
        HYDRO_ASSERT(rgFluxNodes[iNode] == nullptr);
        HYDRO_ASSERT(pFluxNode != nullptr);

        pFluxNode->RequestPressureSource(FluxNode::PRESSURE_FROM_MESH);
        pFluxNode->RequestNodeData(FluxNode::ND_HYDRAULIC);

        rgFluxNodes[iNode] = pFluxNode;
    }

    LinFD5ComprReynoldsElem::LinFD5ComprReynoldsElem(HydroMesh* pMesh)
        :LinFD5Elem(pMesh, REYNOLDS_ELEM)
    {

    }

    LinFD5ComprReynoldsElem::~LinFD5ComprReynoldsElem()
    {

    }

    void
    LinFD5ComprReynoldsElem::AssRes(SubVectorHandler& WorkVec,
                                    doublereal dCoef,
                                    const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr,
                                    enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    LinFD5ComprReynoldsElem::AssJac(SparseSubMatrixHandler& WorkMat,
                                    doublereal dCoef,
                                    const VectorHandler& XCurr,
                                    const VectorHandler& XPrimeCurr,
                                    enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    LinFD5ComprReynoldsElem::InitialAssRes(SubVectorHandler& WorkVec,
                                           const VectorHandler& XCurr,
                                           enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    LinFD5ComprReynoldsElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                           const VectorHandler& XCurr,
                                           enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }

    void LinFD5ComprReynoldsElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = eFunc & grad::REGULAR_FLAG ? iNumDofMax : 1;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());
    }

    template <typename G>
    void LinFD5ComprReynoldsElem::AssRes(GradientAssVec<G>& WorkVec,
                                         doublereal dCoef,
                                         const GradientVectorHandler<G>& XCurr,
                                         const GradientVectorHandler<G>& XPrimeCurr,
                                         enum FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename G>
    void LinFD5ComprReynoldsElem::InitialAssRes(GradientAssVec<G>& WorkVec,
                                                const GradientVectorHandler<G>& XCurr,
                                                enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    template <typename G>
    void LinFD5ComprReynoldsElem::UnivAssRes(GradientAssVec<G>& WorkVec,
                                             doublereal dCoef,
                                             const GradientVectorHandler<G>& XCurr,
                                             enum FunctionCall func)
    {
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();
        const FluidStateBoundaryCond* const pBoundaryCond = pGeometry->pGetMovingPressBoundCond(rgHydroNodes[iNodeCenter]);

        const integer iFirstIndex = rgHydroNodes[iNodeCenter]->iGetFirstEquationIndex(func);

        const HydroRootElement* const pParent = pGetMesh()->pGetParent();
        const doublereal dEquationScale = pParent->dGetScale(HydroRootElement::SCALE_REYNOLDS_EQU);

        if (pBoundaryCond == nullptr) {
            G rho, drho_dt, h, dh_dt;
            std::array<G, iNumFluxNodes> mdot, w;

            for (index_type i = 0; i < iNumFluxNodes; ++i) {
                if (typeid(G) == typeid(doublereal)) {
                    rgFluxNodes[i]->GetVelocityAvg(w[i], &oDofMap);
                }
                rgFluxNodes[i]->GetMassFluxDens(mdot[i], &oDofMap);
            }

            rgHydroNodes[iNodeCenter]->GetDensity(rho, &oDofMap, dCoef);
            rgHydroNodes[iNodeCenter]->GetClearance(h, &oDofMap);
            rgHydroNodes[iNodeCenter]->GetClearanceDerTime(dh_dt, &oDofMap);
            pGetMesh()->pGetGeometry()->GetNonNegativeClearance(h, h, &dh_dt, &dh_dt);

            if (func & grad::REGULAR_FLAG) {
                rgHydroNodes[iNodeCenter]->GetDensityDerTime(drho_dt, &oDofMap, dCoef);
            } else {
                drho_dt = 0.;
            }

            const G Re = ((mdot[iNodeFlxEast] - mdot[iNodeFlxWest]) / dx
                          + (mdot[iNodeFlzNorth] - mdot[iNodeFlzSouth]) / dz
                          + (drho_dt * h + rho * dh_dt)) * dEquationScale;

            SetMaxTimeStep(w);

            CHECK_NUM_COLS_WORK_SPACE(this, func, Re, iFirstIndex);

            WorkVec.AddItem(iFirstIndex, Re);

            if (func & grad::REGULAR_FLAG) {
                G f, pc{pGetFluid()->dGetRefPressure()};

                if (rgHydroNodes[iNodeCenter]->GetCavitationState() == HydroFluid::FULL_FILM_REGION) {
                    G T, rhoc;

                    rgHydroNodes[iNodeCenter]->GetTemperature(T, &oDofMap, dCoef);

                    pGetFluid()->GetDensity(pc, T, rhoc); // Assume that pc is the cavitation pressure

                    f = (rho - rhoc) * (dEquationScale / (dCoef * pParent->dGetScale(HydroRootElement::SCALE_THETA_DOF)));
                } else {
                    G p;

                    rgHydroNodes[iNodeCenter]->GetPressure(p, &oDofMap, dCoef);

                    f = (p - pc) * (dEquationScale / (dCoef * pParent->dGetScale(HydroRootElement::SCALE_PRESSURE_DOF)));
                }

                CHECK_NUM_COLS_WORK_SPACE(this, func, f, iFirstIndex + 1);

                WorkVec.AddItem(iFirstIndex + 1, f);
            }
        } else {
            const doublereal h = rgHydroNodes[iNodeCenter]->dGetClearance(pBoundaryCond);
            const doublereal rho_pre = pBoundaryCond->dGetDensity(h);
            const doublereal p_pre = pBoundaryCond->dGetPressure(h);

            G rho, p;

            rgHydroNodes[iNodeCenter]->GetDensity(rho, &oDofMap, dCoef);
            rgHydroNodes[iNodeCenter]->GetPressure(p, &oDofMap, dCoef);

            CHECK_NUM_COLS_WORK_SPACE(this, func, p, iFirstIndex);

            WorkVec.AddItem(iFirstIndex, (p - p_pre) * dEquationScale);

            if (func & grad::REGULAR_FLAG) {
                CHECK_NUM_COLS_WORK_SPACE(this, func, rho, iFirstIndex + 1);

                WorkVec.AddItem(iFirstIndex + 1, (rho - rho_pre) * dEquationScale);
            }
        }
    }

    void LinFD5ComprReynoldsElem::SetMaxTimeStep(const std::array<doublereal, iNumFluxNodes>& w) const
    {
        doublereal wx = 0., wz = 0.;

        static const index_type rgNodeIdxFlx[2] = {iNodeFlxWest, iNodeFlxEast};
        static const index_type rgNodeIdxFlz[2] = {iNodeFlzSouth, iNodeFlzNorth};

        for (index_type i = 0; i < 2; ++i) {
            wx = std::max(wx, std::abs(w[rgNodeIdxFlx[i]]));
            wz = std::max(wz, std::abs(w[rgNodeIdxFlz[i]]));
        }

        HydroRootElement* const pParent = pGetMesh()->pGetParent();
        const doublereal CFL = pParent->dGetMaxCFL();
        const doublereal dtMax = CFL / (wx / dx + wz / dz);

        pParent->SetMaxTimeStep(dtMax);
    }

    void LinFD5ComprReynoldsElem::SetMaxTimeStep(const std::array<Gradient<0>, iNumFluxNodes>&) const
    {
        // Do nothing
    }

    LinFD5ThermalElem::LinFD5ThermalElem(HydroMesh* pMesh)
        :LinFD5Elem(pMesh, THERMAL_ELEM),
         dScale(0.)
    {
        std::fill(rgThermNodes.begin(), rgThermNodes.end(), nullptr);
    }

    LinFD5ThermalElem::~LinFD5ThermalElem()
    {
    }

    void LinFD5ThermalElem::SetNode(int iNode, HydroNode* pNode)
    {
        LinFD5Elem::SetNode(iNode, pNode);

        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iGetNumNodesTherm());
        HYDRO_ASSERT(rgThermNodes[iNode] == nullptr);
        HYDRO_ASSERT(iNode != iNodeCenter || dynamic_cast<ThermalActiveNode*>(pNode->pGetThermalNode()) != nullptr);

        rgThermNodes[iNode] = pNode->pGetThermalNode();
    }

    void LinFD5ThermalElem::SetFluxNode(int iNode, FluxNode* pFluxNode)
    {
        pFluxNode->RequestNodeData(FluxNode::ND_THERMAL);

        if (pGetMesh()->pGetThermWallBoundCond()) {
            pFluxNode->RequestNodeData(FluxNode::ND_THERMAL_WALL);
        }

        LinFD5Elem::SetFluxNode(iNode, pFluxNode);
    }

    int LinFD5ThermalElem::iGetNumNodesTherm() const
    {
        return rgThermNodes.size();
    }

    ThermoHydrNode* LinFD5ThermalElem::pGetNodeTherm(int iNode) const
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iGetNumNodesTherm());
        HYDRO_ASSERT(rgThermNodes[iNode] != nullptr);

        return rgThermNodes[iNode];
    }

    void LinFD5ThermalElem::Initialize()
    {
        LinFD5Elem::Initialize();

        HYDRO_ASSERT(dA > 0.);
        HYDRO_ASSERT(dx > 0.);
        HYDRO_ASSERT(dz > 0.);

#if HYDRO_DEBUG > 0
        for (index_type i = 0; i < iNumNodes; ++i) {
            HYDRO_ASSERT(rgHydroNodes[i]->GetPosition2D()(1) == rgThermNodes[i]->GetPosition2D()(1));
            HYDRO_ASSERT(rgHydroNodes[i]->GetPosition2D()(2) == rgThermNodes[i]->GetPosition2D()(2));
        }
#endif
        std::array<doublereal, 2> xlim, zlim;

        static const index_type idx[] = {iNodeWest, iNodeEast};
        static const index_type idz[] = {iNodeSouth, iNodeNorth};

        // Element size will be either to the center between the central node and a neighbour active node or to a neighbour passive or slave node.
        // This is needed for correct heat flux to thermal wall nodes
        for (index_type i = 0; i < 2; ++i) {
            if (rgThermNodes[idx[i]]->bIsNodeType(Node2D::ACTIVE_NODE) && rgThermNodes[idx[i]]->bIsNodeType(Node2D::MASTER_NODE)) {
                xlim[i] = 0.5 * (x[idx[i]](1) + x[iNodeCenter](1));
            } else {
                xlim[i] = x[idx[i]](1);
            }

            if (rgThermNodes[idz[i]]->bIsNodeType(Node2D::ACTIVE_NODE) && rgThermNodes[idz[i]]->bIsNodeType(Node2D::MASTER_NODE)) {
                zlim[i] = 0.5 * (x[idz[i]](2) + x[iNodeCenter](2));
            } else {
                zlim[i] = x[idz[i]](2);
            }
        }

        HYDRO_ASSERT(xlim[1] > xlim[0]);
        HYDRO_ASSERT(zlim[1] > zlim[0]);

        HYDRO_ASSERT((xlim[1] - xlim[0]) * (zlim[1] - zlim[0]) >= dA * (1. - sqrt(std::numeric_limits<doublereal>::epsilon())));

        dA = (xlim[1] - xlim[0]) * (zlim[1] - zlim[0]);

        dScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_ENERGY_EQ);

        HYDRO_ASSERT(dA > 0);
    }

    void LinFD5ThermalElem::SetMaxTimeStep(const std::array<doublereal, 2>& wxi, const std::array<doublereal, 2>& wzi) const
    {
        doublereal wx = 0., wz = 0.;

        for (index_type i = 0; i < 2; ++i) {
            wx = std::max(wx, std::abs(wxi[i]));
            wz = std::max(wz, std::abs(wzi[i]));
        }

        HydroRootElement* const pParent = pGetMesh()->pGetParent();
        const doublereal CFL = pParent->dGetMaxCFL();
        const doublereal dtMax = CFL / (wx / dx + wz / dz);

        pParent->SetMaxTimeStep(dtMax);
    }

    void LinFD5ThermalElem::UpwindWeight(const std::array<doublereal, 2>& q, std::array<doublereal, 2>& alpha)
    {
        // First order upwind differences
        alpha[0] = q[0] >= 0.;
        alpha[1] = q[1] < 0.;
    }

    template <typename G>
    void LinFD5ThermalElem::EnergyBalance(G& f, G& Qdot0, G& Qdoth, LocalDofMap* pDofMap, doublereal dCoef, grad::FunctionCall func) const
    {
        const index_type iFluxEval = 2;

        std::array<doublereal, iFluxEval> wx, wz, qx, qz, alphax, alphaz;

        static const index_type rgFluxIdx[iFluxEval] = {iNodeFlxWest, iNodeFlxEast};
        static const index_type rgFluxIdz[iFluxEval] = {iNodeFlzSouth, iNodeFlzNorth};

        for (index_type i = 0; i < iFluxEval; ++i) {
            rgFluxNodes[rgFluxIdx[i]]->GetVelocityAvg(wx[i], pDofMap);
            rgFluxNodes[rgFluxIdz[i]]->GetVelocityAvg(wz[i], pDofMap);
            rgFluxNodes[rgFluxIdx[i]]->GetVolumeFluxDens(qx[i], pDofMap);
            rgFluxNodes[rgFluxIdz[i]]->GetVolumeFluxDens(qz[i], pDofMap);
        }

        SetMaxTimeStep(wx, wz);

        UpwindWeight(qx, alphax);
        UpwindWeight(qz, alphaz);

        const ThermWallBoundCond* const pWall = pGetMesh()->pGetThermWallBoundCond();

        std::array<G, iFluxEval> Qx, Qz, A0x, A0z, Ahx, Ahz, Acx, Acz;

        for (index_type i = 0; i < iFluxEval; ++i) {
            rgFluxNodes[rgFluxIdx[i]]->GetEnergyBalance(Qx[i], pDofMap);
            rgFluxNodes[rgFluxIdz[i]]->GetEnergyBalance(Qz[i], pDofMap);

            if (pWall) {
                rgFluxNodes[rgFluxIdx[i]]->GetDissipationFactors(A0x[i], Ahx[i], Acx[i], pDofMap);
                rgFluxNodes[rgFluxIdz[i]]->GetDissipationFactors(A0z[i], Ahz[i], Acz[i], pDofMap);
            }
        }

        std::array<G, iNumNodes> T;

        for (index_type i = 0; i < iNumNodes; ++i) {
            rgThermNodes[i]->GetTemperature(T[i], pDofMap, dCoef);
        }

        G p, rho, h, dT_dt, dp_dt, cp, lambda, Pfc;

        rgHydroNodes[iNodeCenter]->GetPressure(p, pDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetDensity(rho, pDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetClearance(h, pDofMap);
        pGetMesh()->pGetGeometry()->GetNonNegativeClearance(h, h);

        rgThermNodes[iNodeCenter]->GetTemperatureDerTime(dT_dt, pDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetPressureDerTime(dp_dt, pDofMap, dCoef);
        pGetFluid()->GetSpecificHeat(p, T[iNodeCenter], rho, cp, HydroFluid::SPEC_HEAT_TRUE);
        pGetFluid()->GetThermalConductivity(T[iNodeCenter], rho, lambda);
        rgHydroNodes[iNodeCenter]->GetContactFrictionLossDens(Pfc, pDofMap);

        doublereal beta = 0.;

        if (p > pGetFluid()->dGetRefPressure()) {
            doublereal rhoc, drhoc_dT;

            pGetFluid()->GetDensity(pGetFluid()->dGetRefPressure(),
                                    pGetFluid()->dGetRefTemperature(),
                                    rhoc,
                                    nullptr,
                                    &drhoc_dT);

            beta = -drhoc_dT / rhoc;
        }

        const G dTE_dx = (T[iNodeEast] - T[iNodeCenter]) / (x[iNodeEast](1) - x[iNodeCenter](1));
        const G dTW_dx = (T[iNodeCenter] - T[iNodeWest]) / (x[iNodeCenter](1) - x[iNodeWest](1));
        const G dTN_dz = (T[iNodeNorth] - T[iNodeCenter]) / (x[iNodeNorth](2) - x[iNodeCenter](2));
        const G dTS_dz = (T[iNodeCenter] - T[iNodeSouth]) / (x[iNodeCenter](2) - x[iNodeSouth](2));

        f = h * (lambda * ((dTE_dx - dTW_dx) / (0.5 * (x[iNodeEast](1) - x[iNodeWest](1)))
                           + (dTN_dz - dTS_dz) / (0.5 * (x[iNodeNorth](2) - x[iNodeSouth](2))))
                 - rho * cp * dT_dt + beta * T[iNodeCenter] * dp_dt)
            + Pfc;

        G Phi0{0.}, Phih{0.};

        for (index_type i = 0; i < iFluxEval; ++i) {
            if (pWall && (func & grad::REGULAR_FLAG)) {
                // Should be most accurate in this way
                const G Ac = fabs(qx[i]) >= fabs(qz[i])
                    ? alphax[i] * Acx[i]
                    : alphaz[i] * Acz[i];

                Phi0 += alphax[i] * A0x[i]
                    + alphaz[i] * A0z[i] + Ac;

                Phih += alphax[i] * Ahx[i]
                    + alphaz[i] * Ahz[i] + Ac;
            }

            f += alphax[i] * Qx[i] + alphaz[i] * Qz[i];
        }

        if (pWall && (func & grad::REGULAR_FLAG)) {
            G T0, Th;

            pWall->GetWallTemperature(ThermWallBoundCond::TW_SHAFT, Th, dCoef, func, pDofMap);
            pWall->GetWallTemperature(ThermWallBoundCond::TW_BEARING, T0, dCoef, func, pDofMap);

            const G h2_lambda = h * h / lambda;

            const G a1 = h2_lambda * (1. / 8. * Phi0 - 1. / 24. * Phih) + 5 * T[iNodeCenter] - 3. / 2. * Th - 7. / 2. * T0;
            const G a2 = -0.5 * h2_lambda * Phi0;
            const G a3 = h2_lambda * (1. / 4. * Phih + 7. / 12. * Phi0) - 10. * T[iNodeCenter] + 5 * (Th + T0);
            const G a4 = -5. / 24 * h2_lambda * (Phih + Phi0) + 5 * T[iNodeCenter] - 5. / 2. * (Th + T0);

            Qdot0 = -a1 * lambda / h;
            Qdoth = lambda / h * (a1 + 2 * a2 + 3 * a3 + 4 * a4);

            f += Qdot0 + Qdoth;
        }
    }

    LinFD5ThermalElemImp::LinFD5ThermalElemImp(HydroMesh* pMesh, bool bDoInitAss)
        :LinFD5ThermalElem(pMesh),
         bDoInitAss(bDoInitAss)
    {

    }

    LinFD5ThermalElemImp::~LinFD5ThermalElemImp()
    {
    }

    void
    LinFD5ThermalElemImp::AssRes(SubVectorHandler& WorkVec,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    LinFD5ThermalElemImp::AssJac(SparseSubMatrixHandler& WorkMat,
                                 doublereal dCoef,
                                 const VectorHandler& XCurr,
                                 const VectorHandler& XPrimeCurr,
                                 enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    LinFD5ThermalElemImp::InitialAssRes(SubVectorHandler& WorkVec,
                                        const VectorHandler& XCurr,
                                        enum GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<doublereal>::InitialAssRes(this,
                                                      WorkVec,
                                                      XCurr,
                                                      INITIAL_ASS_RES,
                                                      mode);
        }
    }

    void
    LinFD5ThermalElemImp::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                        const VectorHandler& XCurr,
                                        enum GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                        WorkMat,
                                                        XCurr,
                                                        INITIAL_ASS_JAC,
                                                        &oDofMap,
                                                        mode);
        }
    }

    void LinFD5ThermalElemImp::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = 1;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());

        if ((eFunc & grad::REGULAR_FLAG) && pGetMesh()->pGetThermWallBoundCond()) {
            *piNumRows += 2;
            *piNumCols += 2;
        }
    }

    template <typename U>
    void LinFD5ThermalElemImp::AssRes(GradientAssVec<U>& WorkVec,
                                      doublereal dCoef,
                                      const GradientVectorHandler<U>& XCurr,
                                      const GradientVectorHandler<U>& XPrimeCurr,
                                      enum FunctionCall func)
    {
        UnivAssRes(WorkVec,
                   dCoef,
                   XCurr,
                   func);
    }

    template <typename T>
    void LinFD5ThermalElemImp::InitialAssRes(GradientAssVec<T>& WorkVec,
                                             const GradientVectorHandler<T>& XCurr,
                                             enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }


    template <typename G>
    void LinFD5ThermalElemImp::UnivAssRes(GradientAssVec<G>& WorkVec,
                                          doublereal dCoef,
                                          const GradientVectorHandler<G>& XCurr,
                                          enum FunctionCall func)
    {
        G f(0.), Qdot0(0.), Qdoth(0.);

        EnergyBalance(f, Qdot0, Qdoth, &oDofMap, dCoef, func);

        const ThermWallBoundCond* pWall = pGetMesh()->pGetThermWallBoundCond();

        if (pWall && (func & grad::REGULAR_FLAG)) {
            pWall->AddHeatFlux<G>(ThermWallBoundCond::TW_SHAFT, -dA * Qdoth, WorkVec);
            pWall->AddHeatFlux<G>(ThermWallBoundCond::TW_BEARING, -dA * Qdot0, WorkVec);
        }

        f *= dScale;

        CHECK_NUM_COLS_WORK_SPACE(this, func, f, rgThermNodes[iNodeCenter]->iGetFirstEquationIndex(func));

        WorkVec.AddItem(rgThermNodes[iNodeCenter]->iGetFirstEquationIndex(func), f);
    }

    LinFD5ThermalCouplingElem::LinFD5ThermalCouplingElem(HydroMesh* pMesh,
                                                         bool bDoInitAss)
        :LinFD5ThermalElem(pMesh),
         pInletNode(nullptr),
         bDoInitAss(bDoInitAss)
    {

    }

    LinFD5ThermalCouplingElem::~LinFD5ThermalCouplingElem()
    {
    }

    void LinFD5ThermalCouplingElem::Initialize()
    {
        LinFD5ThermalElem::Initialize();

        pInletNode = dynamic_cast<ThermalInletNode*>(rgThermNodes[iNodeCenter]);

        HYDRO_ASSERT(pInletNode != nullptr);
    }

    void
    LinFD5ThermalCouplingElem::AssRes(SubVectorHandler& WorkVec,
                                      doublereal dCoef,
                                      const VectorHandler& XCurr,
                                      const VectorHandler& XPrimeCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    LinFD5ThermalCouplingElem::AssJac(SparseSubMatrixHandler& WorkMat,
                                      doublereal dCoef,
                                      const VectorHandler& XCurr,
                                      const VectorHandler& XPrimeCurr,
                                      enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    LinFD5ThermalCouplingElem::InitialAssRes(SubVectorHandler& WorkVec,
                                             const VectorHandler& XCurr,
                                             enum GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<doublereal>::InitialAssRes(this,
                                                      WorkVec,
                                                      XCurr,
                                                      INITIAL_ASS_RES,
                                                      mode);
        }
    }

    void
    LinFD5ThermalCouplingElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                             const VectorHandler& XCurr,
                                             enum GradientAssVecBase::mode_t mode)
    {
        if (bDoInitAss) {
            GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                        WorkMat,
                                                        XCurr,
                                                        INITIAL_ASS_JAC,
                                                        &oDofMap,
                                                        mode);
        }
    }

    void LinFD5ThermalCouplingElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        if (eFunc & grad::REGULAR_FLAG) {
            *piNumRows = 3;
            *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());

            if (pGetMesh()->pGetThermWallBoundCond()) {
                *piNumRows += 2;
                *piNumCols += 2;
            }
        } else {
            *piNumRows = *piNumCols = 1;
        }
    }

    template <typename G>
    void LinFD5ThermalCouplingElem::AssRes(GradientAssVec<G>& WorkVec,
                                           doublereal dCoef,
                                           const GradientVectorHandler<G>& XCurr,
                                           const GradientVectorHandler<G>& XPrimeCurr,
                                           enum FunctionCall func)
    {
        std::array<G, iNumFluxNodes> mdot;

        for (index_type i = 0; i < iNumFluxNodes; ++i) {
            rgFluxNodes[i]->GetMassFluxDens(mdot[i], &oDofMap);
        }

        G h, dh_dt, rho, drho_dt, p, T, T_sup, cpm, cpm_sup;

        pInletNode->GetTemperature(T, &oDofMap, dCoef);
        pInletNode->pGetInletNode()->GetTemperature(T_sup, &oDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetPressure(p, &oDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetDensity(rho, &oDofMap, dCoef);
        rgHydroNodes[iNodeCenter]->GetClearance(h, &oDofMap);
        rgHydroNodes[iNodeCenter]->GetClearanceDerTime(dh_dt, &oDofMap);
        pGetMesh()->pGetGeometry()->GetNonNegativeClearance(h, h, &dh_dt, &dh_dt);
        rgHydroNodes[iNodeCenter]->GetDensityDerTime(drho_dt, &oDofMap, dCoef);

        pGetFluid()->GetSpecificHeat(p, T, rho, cpm, HydroFluid::SPEC_HEAT_AVERAGED);
        pGetFluid()->GetSpecificHeat(p, T_sup, rho, cpm_sup, HydroFluid::SPEC_HEAT_AVERAGED);

        const G mdot_sup = (mdot[iNodeFlxEast] - mdot[iNodeFlxWest]) * dz
            + (mdot[iNodeFlzNorth] - mdot[iNodeFlzSouth]) * dx
            + (drho_dt * h + rho * dh_dt) * (dx * dz);

        const G Qdot_sup = mdot_sup * (cpm_sup * T_sup - cpm * T);

        G f(0.), Qdot0(0.), Qdoth(0.);

        EnergyBalance(f, Qdot0, Qdoth, &oDofMap, dCoef, func);

        f += Qdot_sup / (dx * dz);
        f *= dScale;

        const ThermWallBoundCond* pWall = pGetMesh()->pGetThermWallBoundCond();

        if (pWall) {
            pWall->AddHeatFlux<G>(ThermWallBoundCond::TW_SHAFT, -dA * Qdoth, WorkVec);
            pWall->AddHeatFlux<G>(ThermWallBoundCond::TW_BEARING, -dA * Qdot0, WorkVec);
        }

        const integer iFirstIndexHydro = rgHydroNodes[iNodeCenter]->iGetFirstEquationIndex(func);
        const integer iFirstIndexTherm = rgThermNodes[iNodeCenter]->iGetFirstEquationIndex(func);
        const integer iFirstIndexExTherm = pInletNode->pGetInletNode()->iGetFirstEquationIndex(func);

        CHECK_NUM_COLS_WORK_SPACE(this, func, mdot_sup, iFirstIndexHydro);
        CHECK_NUM_COLS_WORK_SPACE(this, func, f, iFirstIndexTherm);
        CHECK_NUM_COLS_WORK_SPACE(this, func, Qdot_sup, iFirstIndexExTherm);

        WorkVec.AddItem(iFirstIndexHydro, mdot_sup);
        WorkVec.AddItem(iFirstIndexTherm, f);
        WorkVec.AddItem(iFirstIndexExTherm, -Qdot_sup);
    }

    template <typename G>
    void LinFD5ThermalCouplingElem::InitialAssRes(GradientAssVec<G>& WorkVec,
                                                  const GradientVectorHandler<G>& XCurr,
                                                  enum FunctionCall func)
    {
        G T, T_sup;

        pInletNode->GetTemperature(T, &oDofMap, 1.);
        pInletNode->pGetInletNode()->GetTemperature(T_sup, &oDofMap, 1.);

        const integer iFirstIndexTherm = pInletNode->iGetFirstEquationIndex(func);

        WorkVec.AddItem(iFirstIndexTherm, T - T_sup);
    }


    QuadFeIso9Elem::IntegrationRule::IntegrationRule(ElementType eElemType,
                                                     index_type iGaussMin,
                                                     index_type iGaussMax,
                                                     index_type iGaussStep)
        :iGaussMin(iGaussMin),
         iGaussMax(iGaussMax),
         iGaussStep(iGaussStep),
         eElemType(eElemType)
    {
        HYDRO_ASSERT(bInvariant());
    }

    void QuadFeIso9Elem::IntegrationRule::ParseInput(MBDynParser& HP, const HydroMesh* pMesh)
    {
        HYDRO_ASSERT(bInvariant());

        index_type iGaussLimMin = eElemType == REYNOLDS_ELEM ? 2 : 1;

        if (HP.IsKeyWord("min") ) {
            NO_OP;
        }

        iGaussMin = iGaussMax = HP.GetInt();

        if (iGaussMin < iGaussLimMin || iGaussMin > 10) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pMesh->pGetParent()->GetLabel()
                        << "): minimum number of gauss points must be in range " << iGaussLimMin << ":10 at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("max")) {
            iGaussMax = HP.GetInt();

            if (iGaussMax < iGaussMin) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pMesh->pGetParent()->GetLabel()
                            << "): maximum number of gauss points must not be less than "
                            "minimum number of gauss points at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (iGaussMax < iGaussLimMin || iGaussMax > 10) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pMesh->pGetParent()->GetLabel()
                            << "): maximum number of gauss points must be in range " << iGaussLimMin << ":10 at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (HP.IsKeyWord("step")) {
            iGaussStep = HP.GetInt();

            if (iGaussStep < 1 || iGaussStep > iGaussMax - iGaussMin) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pMesh->pGetParent()->GetLabel()
                            << "): gauss point step size must be in range 1:"
                            << iGaussMax - iGaussMin << " at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        HYDRO_ASSERT(bInvariant());
    }

#if HYDRO_DEBUG > 0
    bool QuadFeIso9Elem::IntegrationRule::bInvariant() const {
        HYDRO_ASSERT(iGaussMin >= 1);
        HYDRO_ASSERT(iGaussMax <= 10);
        HYDRO_ASSERT(iGaussStep >= 1);
        HYDRO_ASSERT(iGaussStep <= iGaussMax - iGaussMin + (iGaussMin == iGaussMax));
        HYDRO_ASSERT(iGaussMax >= iGaussMin);
        return true;
    }
#endif

    const QuadFeIso9Elem::NodeGroup QuadFeIso9Elem::rgNodeGroups[iNumNodeGroups] = {
        { 1,  1, {0, 4, 7, 8}},
        {-1,  1, {1, 4, 5, 8}},
        {-1, -1, {2, 5, 6, 8}},
        { 1, -1, {3, 6, 7, 8}}
    };

    QuadFeIso9Elem::QuadFeIso9Elem(HydroMesh* pMesh, const IntegrationRule& oIntegRule, ElementType eType)
        :HydroElement(pMesh, eType),
         iCurrIntegRule(-1)
    {
        rgGauss.reserve(oIntegRule.iGetGaussCount());

        for (index_type i = oIntegRule.iGetGaussFirst(); i <= oIntegRule.iGetGaussLast(); i = oIntegRule.iGetGaussNext(i)) {
            HYDRO_ASSERT(rgGauss.size() < rgGauss.capacity());
            rgGauss.emplace_back(i);
        }

        std::fill(rgNodes.begin(), rgNodes.end(), nullptr);

        if (pedantic_output) {
            for (index_type k = 0; k < iGetNumIntegrationRules(); ++k) {
                pedantic_cout("integration rule " << k << ":\n");

                for (index_type i = 1; i <= iGetNumGaussPoints(k); ++i) {
                    pedantic_cout("\tGauss point(" << i
                                  << "): r = " << std::setprecision(16) << dGetGaussPos(i, k)
                                  << ", alpha = " << std::setprecision(16) << dGetGaussWeight(i, k)
                                  << '\n');
                }
            }

            pedantic_cout(std::endl);
        }
    }

    QuadFeIso9Elem::~QuadFeIso9Elem()
    {
    }

    int QuadFeIso9Elem::iGetNumNodes() const
    {
        return iNumNodes;
    }

    void QuadFeIso9Elem::SetNode(int iNode, HydroNode* pNode)
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgNodes[iNode] == nullptr);

        rgNodes[iNode] = pNode;
    }

    HydroNode* QuadFeIso9Elem::pGetNode(int iNode) const
    {
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iNumNodes);
        HYDRO_ASSERT(rgNodes[iNode] != nullptr);

        return rgNodes[iNode];
    }

    void QuadFeIso9Elem::AfterPredict(VectorHandler& X, VectorHandler& XP)
    {
        iCurrIntegRule = -1;
    }

    void QuadFeIso9Elem::AfterConvergence(const VectorHandler& X,
                                          const VectorHandler& XP)
    {
        iCurrIntegRule = -1;
    }

    index_type QuadFeIso9Elem::iGetNumGaussPoints(index_type iIntegRule) const
    {
        HYDRO_ASSERT(iIntegRule >= 0);
        HYDRO_ASSERT(size_t(iIntegRule) < rgGauss.size());
        return rgGauss[iIntegRule].iGetNum();
    }

    index_type QuadFeIso9Elem::iGetNumIntegrationRules() const
    {
        return rgGauss.size();
    }

    index_type QuadFeIso9Elem::iGetGaussPointSize() const
    {
        integer iSize = iGetGaussPointSize1D();

        return iSize * iSize;
    }

    index_type QuadFeIso9Elem::iGetGaussPointSize1D() const
    {
        integer iSize = 0;

        for (auto i = rgGauss.begin(); i != rgGauss.end(); ++i) {
            iSize += i->iGetNum();
        }

        return iSize;
    }

    index_type QuadFeIso9Elem::iGetGaussPointIndex(index_type iGaussR, index_type iGaussS, index_type iIntegRule) const
    {
        HYDRO_ASSERT(iIntegRule >= 0);
        HYDRO_ASSERT(iIntegRule < iGetNumIntegrationRules());

        index_type iOffset = 0;

        for (index_type i = 0; i < iIntegRule; ++i) {
            const index_type iNumGaussCurr = iGetNumGaussPoints(i);
            iOffset +=  iNumGaussCurr * iNumGaussCurr;
        }

        HYDRO_ASSERT(iGaussR >= 1);
        HYDRO_ASSERT(iGaussR <= iGetNumGaussPoints(iIntegRule));
        HYDRO_ASSERT(iGaussS >= 1);
        HYDRO_ASSERT(iGaussS <= iGetNumGaussPoints(iIntegRule));

        const index_type iGaussPointIdx = iOffset + (iGaussR - 1) * iGetNumGaussPoints(iIntegRule) + iGaussS - 1;

        HYDRO_ASSERT(iGaussPointIdx >= 0);
        HYDRO_ASSERT(iGaussPointIdx <= iOffset + iGetNumGaussPoints(iIntegRule) * iGetNumGaussPoints(iIntegRule));

        return iGaussPointIdx;
    }

    index_type QuadFeIso9Elem::iGetGaussPointIndex1D(index_type iGaussR, index_type iIntegRule) const
    {
        HYDRO_ASSERT(iIntegRule >= 0);
        HYDRO_ASSERT(iIntegRule < iGetNumIntegrationRules());

        index_type iOffset = 0;

        for (index_type i = 0; i < iIntegRule; ++i) {
            const index_type iNumGaussCurr = iGetNumGaussPoints(i);
            iOffset +=  iNumGaussCurr;
        }

        HYDRO_ASSERT(iGaussR >= 1);
        HYDRO_ASSERT(iGaussR <= iGetNumGaussPoints(iIntegRule));

        const index_type iGaussPointIdx = iOffset + (iGaussR - 1);

        HYDRO_ASSERT(iGaussPointIdx >= 0);
        HYDRO_ASSERT(iGaussPointIdx <= iOffset + iGetNumGaussPoints(iIntegRule));

        return iGaussPointIdx;
    }

    doublereal QuadFeIso9Elem::dGetGaussWeight(index_type iGaussPoint, index_type iIntegRule) const
    {
        HYDRO_ASSERT(iIntegRule >= 0);
        HYDRO_ASSERT(iIntegRule < iGetNumIntegrationRules());
        HYDRO_ASSERT(iGaussPoint >= 1);
        HYDRO_ASSERT(iGaussPoint <= iGetNumGaussPoints(iIntegRule));

        return rgGauss[iIntegRule].dGetWght(iGaussPoint);
    }

    doublereal QuadFeIso9Elem::dGetGaussPos(index_type iGaussPoint, index_type iIntegRule) const
    {
        HYDRO_ASSERT(iIntegRule >= 0);
        HYDRO_ASSERT(iIntegRule < iGetNumIntegrationRules());
        HYDRO_ASSERT(iGaussPoint >= 1);
        HYDRO_ASSERT(iGaussPoint <= iGetNumGaussPoints(iIntegRule));

        return rgGauss[iIntegRule].dGetPnt(iGaussPoint);
    }

    index_type QuadFeIso9Elem::iSelectIntegrationRule(const Vector<doublereal, iNumNodes>& pe,
                                                      const Vector<doublereal, iNumNodes>*const paspe)
    {
        if (iCurrIntegRule == -1) {
            const index_type iMaxInteg = iGetNumIntegrationRules() - 1;

            if (iMaxInteg) {
                doublereal pmin = std::numeric_limits<doublereal>::max();
                doublereal pmax = -pmin;

                for (index_type i = 1; i <= pe.iGetNumRows(); ++i) {
                    doublereal pi = pe(i);

                    if (paspe) {
                        pi += (*paspe)(i);
                    }

                    pmax = std::max(pi, pmax);
                    pmin = std::min(pi, pmin);
                }

                const doublereal dp = (pmax - pmin) / std::max(1., pGetMesh()->dGetMaxPressureGradient());

                HYDRO_ASSERT(dp >= 0);
                HYDRO_ASSERT(dp <= 1);

                iCurrIntegRule = std::ceil(dp * iMaxInteg);

                if (iCurrIntegRule > iMaxInteg) {
                    iCurrIntegRule = iMaxInteg;
                }

                HYDRO_TRACE("new integration rule:" << iCurrIntegRule << "\n");
                HYDRO_TRACE("\tiGetNumIntegrationRules()=" << iGetNumIntegrationRules() << "\n");
                HYDRO_TRACE("\tdp=" << std::setw(16) << dp << "\n");
                HYDRO_TRACE("\tiCurrIntegRule=" << iCurrIntegRule << "\n\n");
            } else {
                iCurrIntegRule = 0;
            }
        }

        HYDRO_ASSERT(iCurrIntegRule >= 0);
        HYDRO_ASSERT(iCurrIntegRule < iGetNumIntegrationRules());

        return iCurrIntegRule;
    }

    index_type QuadFeIso9Elem::iSelectIntegrationRule(const Vector<Gradient<0>, iNumNodes>&,
                                                      const Vector<Gradient<0>, iNumNodes>*)
    {
        HYDRO_ASSERT(iCurrIntegRule >= 0);
        HYDRO_ASSERT(iCurrIntegRule < iGetNumIntegrationRules());

        return iCurrIntegRule;
    }


    void QuadFeIso9Elem::NodePositionMatrix(Matrix<doublereal, iNumNodes, 2>& xe) const
    {
        for (index_type i = 1; i <= iNumNodes; ++i) {
            const Vector<doublereal, 2>& xni = rgNodes[i - 1]->GetPosition2D();

            for (index_type j = 1; j <= 2; ++j) {
                xe(i, j) = xni(j);
            }
        }
    }

    void QuadFeIso9Elem::PressureInterpolMatrix(Vector<doublereal, iNumNodes>& N,
                                                const doublereal r,
                                                const doublereal s) const
    {
        N(1) = ((r*r+r)*s*s+(r*r+r)*s)/4.0e0;
        N(2) = ((r*r-r)*s*s+(r*r-r)*s)/4.0e0;
        N(3) = ((r*r-r)*s*s+(r-r*r)*s)/4.0e0;
        N(4) = ((r*r+r)*s*s+((-r*r)-r)*s)/4.0e0;
        N(5) = -((r*r-1)*s*s+(r*r-1)*s)/2.0e0;
        N(6) = -((r*r-r)*s*s-r*r+r)/2.0e0;
        N(7) = -((r*r-1)*s*s+(1-r*r)*s)/2.0e0;
        N(8) = -((r*r+r)*s*s-r*r-r)/2.0e0;
        N(9) = (r*r-1)*s*s-r*r+1;
    }

    void QuadFeIso9Elem::PressureInterpolMatrixDer(Matrix<doublereal, 2, iNumNodes>& dN_drv,
                                                   const doublereal r,
                                                   const doublereal s) const
    {
        dN_drv(1,1) = ((2*r+1)*s*s+(2*r+1)*s)/4.0e0;
        dN_drv(1,2) = ((2*r-1)*s*s+(2*r-1)*s)/4.0e0;
        dN_drv(1,3) = ((2*r-1)*s*s+(1-2*r)*s)/4.0e0;
        dN_drv(1,4) = ((2*r+1)*s*s+((-2*r)-1)*s)/4.0e0;
        dN_drv(1,5) = -(2*r*s*s+2*r*s)/2.0e0;
        dN_drv(1,6) = -((2*r-1)*s*s-2*r+1)/2.0e0;
        dN_drv(1,7) = -(2*r*s*s-2*r*s)/2.0e0;
        dN_drv(1,8) = -((2*r+1)*s*s-2*r-1)/2.0e0;
        dN_drv(1,9) = 2*r*s*s-2*r;
        dN_drv(2,1) = (2*(r*r+r)*s+r*r+r)/4.0e0;
        dN_drv(2,2) = (2*(r*r-r)*s+r*r-r)/4.0e0;
        dN_drv(2,3) = (2*(r*r-r)*s-r*r+r)/4.0e0;
        dN_drv(2,4) = (2*(r*r+r)*s-r*r-r)/4.0e0;
        dN_drv(2,5) = -(2*(r*r-1)*s+r*r-1)/2.0e0;
        dN_drv(2,6) = -(r*r-r)*s;
        dN_drv(2,7) = -(2*(r*r-1)*s-r*r+1)/2.0e0;
        dN_drv(2,8) = -(r*r+r)*s;
        dN_drv(2,9) = 2*(r*r-1)*s;
    }

    void QuadFeIso9Elem::PressureGradInterpolMatrix(Matrix<doublereal, 2, iNumNodes>& B,
                                                    doublereal& detJ,
                                                    const Matrix<doublereal, iNumNodes, 2>& xe,
                                                    const doublereal r,
                                                    const doublereal s) const
    {
        Matrix<doublereal, 2, iNumNodes> dN_drv;

        PressureInterpolMatrixDer(dN_drv, r, s);

        const Matrix<doublereal, 2, 2> J = dN_drv * xe;
        B = Inv(J) * dN_drv;
        detJ = Det(J);

        if (detJ <= 0.) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetMesh()->pGetParent()->GetLabel()
                        << "): element Jacobian matrix is singular\n");
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    const QuadFeIso9Elem::NodeGroup* QuadFeIso9Elem::pFindNodeGroup(doublereal r, doublereal s)
    {
        const doublereal signr = copysign(1, r);
        const doublereal signs = copysign(1, s);

        auto pNodeGroup = std::find_if(std::begin(rgNodeGroups),
                                       std::end(rgNodeGroups),
                                       [signr, signs] (const NodeGroup& g) { return g.r == signr && g.s == signs; });

        HYDRO_ASSERT(pNodeGroup != std::end(rgNodeGroups));

        return &*pNodeGroup;
    }

    QuadFeIso9ReynoldsElem::QuadFeIso9ReynoldsElem(HydroMesh* pMesh, const IntegrationRule& oIntegRule)
        :QuadFeIso9Elem(pMesh, oIntegRule, REYNOLDS_ELEM),
         rgGaussPntDat(iGetGaussPointSize()),
         dA(0.)
    {
    }

    QuadFeIso9ReynoldsElem::~QuadFeIso9ReynoldsElem()
    {
    }

    void
    QuadFeIso9ReynoldsElem::AssRes(SubVectorHandler& WorkVec,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    QuadFeIso9ReynoldsElem::AssJac(SparseSubMatrixHandler& WorkMat,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    QuadFeIso9ReynoldsElem::InitialAssRes(SubVectorHandler& WorkVec,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    QuadFeIso9ReynoldsElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }

    void QuadFeIso9ReynoldsElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = iNumNodes;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());
    }

    template <typename T>
    void QuadFeIso9ReynoldsElem::AssRes(GradientAssVec<T>& WorkMat,
                                        doublereal dCoef,
                                        const GradientVectorHandler<T>& XCurr,
                                        const GradientVectorHandler<T>& XPrimeCurr,
                                        enum FunctionCall func)
    {
        UnivAssRes(WorkMat, dCoef, XCurr, func);
    }

    template <typename T>
    void QuadFeIso9ReynoldsElem::InitialAssRes(GradientAssVec<T>& WorkVec,
                                               const GradientVectorHandler<T>& XCurr,
                                               enum FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    template <typename T>
    void QuadFeIso9ReynoldsElem::UnivAssRes(GradientAssVec<T>& WorkVec,
                                            doublereal dCoef,
                                            const GradientVectorHandler<T>& XCurr,
                                            enum FunctionCall func)
    {
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();
        Vector<T, iNumNodes> pe, he, hdote;
        Vector<Vector<T, 2>, iNumNodes> Ue;
        Vector<doublereal, iNumNodes> etae;
        std::array<index_type, iNumNodes> rgActiveNodes;

        struct MovingBoundCond {
            const FluidStateBoundaryCond* pMovingBound;
            index_type iNodeIdx;
        };

        std::array<MovingBoundCond, iNumNodes> rgMovingBound;
        index_type iNumActiveNodes = 0, iNumMovingBound = 0;

        for (index_type i = 1; i <= iNumNodes; ++i) {
            rgNodes[i - 1]->GetPressure(pe(i), &oDofMap, dCoef);
            rgNodes[i - 1]->GetHydraulicVelocity(Ue(i), &oDofMap);
            rgNodes[i - 1]->GetClearance(he(i), &oDofMap);
            rgNodes[i - 1]->GetClearanceDerTime(hdote(i), &oDofMap);
            rgNodes[i - 1]->GetViscosity(etae(i), &oDofMap);

            if (rgNodes[i - 1]->bIsNodeType(HydroNode::ACTIVE_NODE)) {
                const FluidStateBoundaryCond* const pMovingBoundCond = pGeometry->pGetMovingPressBoundCond(rgNodes[i - 1]);

                if (pMovingBoundCond == nullptr) {
                    rgActiveNodes[iNumActiveNodes++] = i - 1;
                } else {
                    rgMovingBound[iNumMovingBound].pMovingBound = pMovingBoundCond;
                    rgMovingBound[iNumMovingBound].iNodeIdx = i - 1;
                    ++iNumMovingBound;
                }
            }
        }

        const index_type iIntegRule = iSelectIntegrationRule(pe);
        const doublereal dEquationScale = pGetMesh()->pGetParent()->dGetScale(HydroRootElement::SCALE_REYNOLDS_EQU) / dCoef;

        if (iNumActiveNodes > 0) {
            const index_type iNumGauss = iGetNumGaussPoints(iIntegRule);
            Vector<T, iNumNodes> fe;

            for (index_type i = 1; i <= iNumGauss; ++i) {
                for (index_type j = 1; j <= iNumGauss; ++j) {
                    const index_type idxGauss = iGetGaussPointIndex(i, j, iIntegRule);
                    const doublereal alpha = rgGaussPntDat[idxGauss].detJ
                        * dGetGaussWeight(i, iIntegRule)
                        * dGetGaussWeight(j, iIntegRule)
                        * dEquationScale / dA; // do not scale the residual by surface area
                    const Vector<doublereal, iNumNodes>& N = rgGaussPntDat[idxGauss].N;
                    const Matrix<doublereal, 2, iNumNodes>& B = rgGaussPntDat[idxGauss].B;

                    T h = Dot(N, he);
                    T a0 = Dot(N, hdote);

                    pGeometry->GetNonNegativeClearance(h, h, &a0, &a0);

                    const T h3_alpha = h * h * h * alpha;

                    for (index_type k = 1; k <= 2; ++k) {
                        const T a1 = Dot(B.GetRow(k), pe) * h3_alpha;

                        for (index_type l = 0; l < iNumActiveNodes; ++l) {
                            const index_type iRowIdx = rgActiveNodes[l] + 1;
                            fe(iRowIdx) += B(k, iRowIdx) * a1;
                        }
                    }

                    const doublereal eta = Dot(N, etae);

                    for (index_type k = 1; k <= 2; ++k) {
                        T U(0);

                        for (index_type l = 1; l <= iNumNodes; ++l) {
                            U += N(l) * Ue(l)(k);
                        }
                        a0 += U * Dot(B.GetRow(k), he);
                    }

                    const T a2 = a0 * (12 * eta * alpha);

                    for (index_type l = 0; l < iNumActiveNodes; ++l) {
                        const index_type iRowIdx = rgActiveNodes[l] + 1;
                        fe(iRowIdx) += N(iRowIdx) * a2;
                    }
                }
            }

            for (index_type i = 0; i < iNumActiveNodes; ++i) {
                const index_type iNodeIdx = rgActiveNodes[i];

                CHECK_NUM_COLS_WORK_SPACE(this, func, fe(iNodeIdx + 1), rgNodes[iNodeIdx]->iGetFirstEquationIndex(func));

                WorkVec.AddItem(rgNodes[iNodeIdx]->iGetFirstEquationIndex(func), fe(iNodeIdx + 1));
            }
        }

        for (index_type i = 0; i < iNumMovingBound; ++i) {
            const index_type iNodeIdx = rgMovingBound[i].iNodeIdx;
            const doublereal ppre = rgMovingBound[i].pMovingBound->dGetPressure();
            const T f = (pe(iNodeIdx + 1) - ppre) * dEquationScale;

            CHECK_NUM_COLS_WORK_SPACE(this, func, f, rgNodes[iNodeIdx]->iGetFirstEquationIndex(func));

            WorkVec.AddItem(rgNodes[iNodeIdx]->iGetFirstEquationIndex(func), f);
        }
    }

    void QuadFeIso9ReynoldsElem::Initialize()
    {
        Matrix<doublereal, iNumNodes, 2> xe;

        NodePositionMatrix(xe);

        for (index_type k = 0; k < iGetNumIntegrationRules(); ++k) {
            const index_type iNumGauss = iGetNumGaussPoints(k);
            for (index_type i = 1; i <= iNumGauss; ++i) {
                for (index_type j = 1; j <= iNumGauss; ++j) {
                    const index_type idx = iGetGaussPointIndex(i, j, k);
                    const doublereal r = dGetGaussPos(i, k);
                    const doublereal s = dGetGaussPos(j, k);

                    PressureInterpolMatrix(rgGaussPntDat[idx].N, r, s);
                    PressureGradInterpolMatrix(rgGaussPntDat[idx].B, rgGaussPntDat[idx].detJ, xe, r, s);
                }
            }
        }

        dA = 0.;

        const index_type k = iGetNumIntegrationRules() - 1;
        const index_type iNumGaussPnt = iGetNumGaussPoints(k);

        for (index_type i = 1; i <= iNumGaussPnt; ++i) {
            for (index_type j = 1; j <= iNumGaussPnt; ++j) {
                const index_type idx = iGetGaussPointIndex(i, j, k);
                dA += rgGaussPntDat[idx].detJ * dGetGaussWeight(i, k) * dGetGaussWeight(j, k);
            }
        }

        pedantic_cout("element surface area: " << dA << '\n');
    }

    QuadFeIso9FrictionElem::QuadFeIso9FrictionElem(HydroMesh* pMesh,
                                                   const IntegrationRule& oIntegRule)
        :QuadFeIso9Elem(pMesh, oIntegRule, FRICTION_ELEM),
         rgGaussPntDat(iGetGaussPointSize())
    {
    }

    QuadFeIso9FrictionElem::~QuadFeIso9FrictionElem()
    {
    }

    void QuadFeIso9FrictionElem::Initialize()
    {
        Matrix<doublereal, iNumNodes, 2> xe;

        NodePositionMatrix(xe);

        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        for (index_type k = 0; k < iGetNumIntegrationRules(); ++k) {
            const index_type iNumGauss = iGetNumGaussPoints(k);

            for (index_type i = 1; i <= iNumGauss; ++i) {
                for (index_type j = 1; j <= iNumGauss; ++j) {
                    const index_type idx = iGetGaussPointIndex(i, j, k);
                    const doublereal r = dGetGaussPos(i, k);
                    const doublereal s = dGetGaussPos(j, k);

                    PressureInterpolMatrix(rgGaussPntDat[idx].N, r, s);
                    PressureGradInterpolMatrix(rgGaussPntDat[idx].B, rgGaussPntDat[idx].detJ, xe, r, s);

                    for (index_type k = 1; k <= 2; ++k) {
                        rgGaussPntDat[idx].xc(k) = Dot(rgGaussPntDat[idx].N, xe.GetCol(k));
                    }

                    pGeometry->GetPosition3D(rgGaussPntDat[idx].xc, rgGaussPntDat[idx].vc);
                    pGeometry->GetTangentCoordSys(rgGaussPntDat[idx].xc, rgGaussPntDat[idx].Rtc);
                }
            }
        }
    }

    void
    QuadFeIso9FrictionElem::AssRes(SubVectorHandler& WorkVec,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    QuadFeIso9FrictionElem::AssJac(SparseSubMatrixHandler& WorkMat,
                                   doublereal dCoef,
                                   const VectorHandler& XCurr,
                                   const VectorHandler& XPrimeCurr,
                                   enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    QuadFeIso9FrictionElem::InitialAssRes(SubVectorHandler& WorkVec,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::InitialAssRes(this,
                                                  WorkVec,
                                                  XCurr,
                                                  INITIAL_ASS_RES,
                                                  mode);
    }

    void
    QuadFeIso9FrictionElem::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                          const VectorHandler& XCurr,
                                          enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::InitialAssJac(this,
                                                    WorkMat,
                                                    XCurr,
                                                    INITIAL_ASS_JAC,
                                                    &oDofMap,
                                                    mode);
    }

    void QuadFeIso9FrictionElem::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = 0;
        *piNumCols = 0;
    }

    template <typename T>
    void QuadFeIso9FrictionElem::AssRes(GradientAssVec<T>& WorkVec,
                                        const doublereal dCoef,
                                        const GradientVectorHandler<T>& XCurr,
                                        const GradientVectorHandler<T>& XPrimeCurr,
                                        const FunctionCall func)
    {
        UnivAssRes(WorkVec, dCoef, XCurr, func);
    }

    template <typename T>
    void QuadFeIso9FrictionElem::InitialAssRes(GradientAssVec<T>& WorkVec,
                                               const GradientVectorHandler<T>& XCurr,
                                               const FunctionCall func)
    {
        UnivAssRes(WorkVec, 1., XCurr, func);
    }

    template <typename T>
    void QuadFeIso9FrictionElem::UnivAssRes(GradientAssVec<T>& WorkMat,
                                            doublereal dCoef,
                                            const GradientVectorHandler<T>& XCurr,
                                            enum FunctionCall func)
    {
        const bool bUpdateFriction = pGetMesh()->pGetParent()->bUpdateFrictionLoss();
        BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();

        Vector<T, iNumNodes> pe, paspe, he, Pfce;
        Vector<Vector<T, 2>, iNumNodes> U1e, U2e, tauc_0e;
        Vector<doublereal, iNumNodes> etae;
        const Vector<T, 2> dM_h_Rt; // FIXME: currently unused
        bool bContact = false;

        for (index_type i = 1; i <= iNumNodes; ++i) {
            rgNodes[i - 1]->GetPressure(pe(i), &oDofMap, dCoef); // use pressure from nodes instead of mesh

            if (rgNodes[i - 1]->GetContactPressure(paspe(i), &oDofMap)) {
                rgNodes[i - 1]->GetContactStress(tauc_0e(i), &oDofMap);
                if (bUpdateFriction) {
                    rgNodes[i - 1]->GetContactFrictionLossDens(Pfce(i), &oDofMap);
                }
                bContact = true;
            }

            rgNodes[i - 1]->GetVelocity(U1e(i), U2e(i), &oDofMap);
            rgNodes[i - 1]->GetClearance(he(i), &oDofMap);
            rgNodes[i - 1]->GetViscosity(etae(i), &oDofMap);
        }

        const index_type iIntegRule = iSelectIntegrationRule(pe, &paspe);
        const index_type iNumGaussPnt = iGetNumGaussPoints(iIntegRule);

        for (index_type i = 1; i <= iNumGaussPnt; ++i) {
            for (index_type j = 1; j <= iNumGaussPnt; ++j) {
                const index_type iGaussIdx = iGetGaussPointIndex(i, j, iIntegRule);
                const doublereal r = dGetGaussPos(i, iIntegRule);
                const doublereal s = dGetGaussPos(j, iIntegRule);
                const Vector<doublereal, iNumNodes>& N = rgGaussPntDat[iGaussIdx].N;
                const Matrix<doublereal, 2, iNumNodes>& B = rgGaussPntDat[iGaussIdx].B;
                const doublereal detJ = rgGaussPntDat[iGaussIdx].detJ;
                T p = Dot(N, pe), dp_dx, dp_dz;

                if (HydroFluid::CAVITATION_REGION == pGetFluid()->Cavitation(p)) {
                    dp_dx = 0.; // set pressure gradient to zero if pressure is negative according to Guembel boundary condition
                    dp_dz = 0.;
                } else {
                    dp_dx = Dot(B.GetRow(1), pe);
                    dp_dz = Dot(B.GetRow(2), pe);
                }

                T h = Dot(N, he);

                pGeometry->GetNonNegativeClearance(h, h);

                const doublereal eta = Dot(N, etae);

                Vector<T, 2> U1, U2, tauc_0;

                for (index_type k = 1; k <= 2; ++k) {
                    for (index_type l = 1; l <= iNumNodes; ++l) {
                        U1(k) += N(l) * U1e(l)(k);
                        U2(k) += N(l) * U2e(l)(k);
                    }
                }

                T pasp{0.};

                if (bContact) {
                    pasp = Dot(N, paspe);

                    if (pasp < 0.) { // Could happen in the transition region
                        pasp = 0.;
                    }

                    for (index_type k = 1; k <= 2; ++k) {
                        for (index_type l = 1; l <= iNumNodes; ++l) {
                            tauc_0(k) += N(l) * tauc_0e(l)(k);
                        }
                    }
                }

                const T ptot = p + pasp;

                const Vector<T, 2> dU = U1 - U2;

                const T tau_xy_p_h = 0.5 * h * dp_dx;
                const T tau_xy_U_h = eta * dU(1) / h;
                const T tau_yz_p_h = 0.5 * h * dp_dz;
                const T tau_yz_U_h = eta * dU(2) / h;

                const T tau_xy_0 = (-tau_xy_p_h + tau_xy_U_h);
                const T tau_yz_0 = (-tau_yz_p_h + tau_yz_U_h);

                const T tau_xy_h = (tau_xy_p_h + tau_xy_U_h);
                const T tau_yz_h = (tau_yz_p_h + tau_yz_U_h);

                SetStress(r, s, tau_xy_0, tau_yz_0, tau_xy_h, tau_yz_h);

                const doublereal dA = dGetGaussWeight(i, iIntegRule) * dGetGaussWeight(j, iIntegRule) * detJ;

                Vector<T, 3> dF_0_Rt(tau_xy_0 * dA,
                                     -ptot * dA,
                                     tau_yz_0 * dA);

                Vector<T, 3> dF_h_Rt(-tau_xy_h * dA,
                                     ptot * dA,
                                     -tau_yz_h * dA);

                if (bUpdateFriction) {
                    AddFrictionLoss(HydroRootBase::FLUID_FRICTION, U2, dF_0_Rt(1), dF_0_Rt(3));
                    AddFrictionLoss(HydroRootBase::FLUID_FRICTION, U1, dF_h_Rt(1), dF_h_Rt(3));
                }

                if (bContact) {
                    if (bUpdateFriction) {
                        AddFrictionLoss(HydroRootBase::CONTACT_FRICTION, -dA * Dot(N, Pfce));
                    }

                    static const int iReactionIdx[2] = {1, 3};

                    for (int i = 1; i <= 2; ++i) {
                        // Note: The contact pressure will be always
                        //               normal to the surface of each part
                        //           but the angle between the surface
                        //               normal vectors might be slightly
                        //               different from 180 degrees!

                        // Note: In the same way the angle between
                        //               the surface tangent vectors at
                        //               two different parts might be different
                        //               from zero. But it will be a small angle
                        //               as long as the relative clearance is small
                        //               in case of a cylindrical bearing.

                        dF_0_Rt(iReactionIdx[i - 1]) += tauc_0(i) * dA;
                        dF_h_Rt(iReactionIdx[i - 1]) -= tauc_0(i) * dA;
                    }
                }

                pGeometry->AddReactionForce(rgGaussPntDat[iGaussIdx].xc,
                                            rgGaussPntDat[iGaussIdx].vc,
                                            rgGaussPntDat[iGaussIdx].Rtc,
                                            dF_0_Rt,
                                            dF_h_Rt,
                                            dM_h_Rt,
                                            dCoef,
                                            func,
                                            &oDofMap);
            }
        }
    }

    void QuadFeIso9FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const Vector<doublereal, 2>& U, doublereal dTau_xy, doublereal dTau_yz) const
    {
        const doublereal dPf = U(1) * dTau_xy + U(2) * dTau_yz;

        AddFrictionLoss(type, dPf);
    }

    void QuadFeIso9FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const Vector<Gradient<0>, 2>& U, const Gradient<0>& dTau_xy, const Gradient<0>& dTau_yz) const
    {
        NO_OP;
    }

    void QuadFeIso9FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, doublereal dPf) const
    {
        pGetMesh()->pGetParent()->AddFrictionLoss(type, dPf);
    }

    void QuadFeIso9FrictionElem::AddFrictionLoss(HydroRootBase::FrictionLossType type, const Gradient<0>& dPf) const
    {
        NO_OP;
    }

    void QuadFeIso9FrictionElem::SetStress(doublereal r,
                                           doublereal s,
                                           doublereal tau_xy_0,
                                           doublereal tau_yz_0,
                                           doublereal tau_xy_h,
                                           doublereal tau_yz_h) const
    {
        const NodeGroup* const pGroup = pFindNodeGroup(r, s);

        for (auto i = std::begin(pGroup->nodes); i != std::end(pGroup->nodes); ++i) {
            rgNodes[*i]->SetStress(tau_xy_0,
                                   tau_yz_0,
                                   tau_xy_h,
                                   tau_yz_h);
        }
    }

    const index_type QuadFeIso9MassFlowZ::rgNodeIndexOutletBound[2][iNumNodesOutletBound] = {{3 - 1, 4 - 1, 7 - 1}, {1 - 1, 2 - 1, 5 - 1}};

    QuadFeIso9MassFlowZ::QuadFeIso9MassFlowZ(HydroMesh* pMesh, const IntegrationRule& oIntegRule, doublereal sref)
        :QuadFeIso9Elem(pMesh, oIntegRule, COUPLING_ELEM),
         rgGaussPntDat(iGetGaussPointSize1D()),
         sref(sref),
         pNodeIndexOutletBound{rgNodeIndexOutletBound[sref > 0.]}
    {
        HYDRO_ASSERT(fabs(sref) == 1.);
    }

    QuadFeIso9MassFlowZ::~QuadFeIso9MassFlowZ()
    {
    }

    void
    QuadFeIso9MassFlowZ::AssRes(SubVectorHandler& WorkVec,
                                doublereal dCoef,
                                const VectorHandler& XCurr,
                                const VectorHandler& XPrimeCurr,
                                enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<doublereal>::AssRes(this,
                                           WorkVec,
                                           dCoef,
                                           XCurr,
                                           XPrimeCurr,
                                           REGULAR_RES,
                                           mode);
    }

    void
    QuadFeIso9MassFlowZ::AssJac(SparseSubMatrixHandler& WorkMat,
                                doublereal dCoef,
                                const VectorHandler& XCurr,
                                const VectorHandler& XPrimeCurr,
                                enum GradientAssVecBase::mode_t mode)
    {
        GradientAssVec<Gradient<0> >::AssJac(this,
                                             WorkMat,
                                             dCoef,
                                             XCurr,
                                             XPrimeCurr,
                                             REGULAR_JAC,
                                             &oDofMap,
                                             mode);
    }

    void
    QuadFeIso9MassFlowZ::InitialAssRes(SubVectorHandler& WorkVec,
                                       const VectorHandler& XCurr,
                                       enum GradientAssVecBase::mode_t mode)
    {
    }

    void
    QuadFeIso9MassFlowZ::InitialAssJac(SparseSubMatrixHandler& WorkMat,
                                       const VectorHandler& XCurr,
                                       enum GradientAssVecBase::mode_t mode)
    {
    }

    void QuadFeIso9MassFlowZ::WorkSpaceDim(integer* piNumRows, integer* piNumCols, grad::FunctionCall eFunc) const
    {
        *piNumRows = iNumNodesOutletBound;
        *piNumCols = iGetNumColsWorkSpace(eFunc, iGetNumNodes());
    }

    template <typename T>
    void QuadFeIso9MassFlowZ::AssRes(GradientAssVec<T>& WorkVec,
                                     doublereal dCoef,
                                     const GradientVectorHandler<T>& XCurr,
                                     const GradientVectorHandler<T>& XPrimeCurr,
                                     enum FunctionCall func)
    {
        const BearingGeometry* const pGeometry = pGetMesh()->pGetGeometry();
        Vector<T, iNumNodes> pe, he;
        Vector<Vector<T, 2>, iNumNodes> Ue;
        Vector<doublereal, iNumNodes> etae, rhoe;

        for (index_type i = 1; i <= iNumNodes; ++i) {
            pGetMesh()->GetPressure(rgNodes[i - 1], pe(i), &oDofMap, dCoef);
            rgNodes[i - 1]->GetHydraulicVelocity(Ue(i), &oDofMap);
            rgNodes[i - 1]->GetClearance(he(i), &oDofMap);
            rgNodes[i - 1]->GetViscosity(etae(i), &oDofMap, dCoef);
            rgNodes[i - 1]->GetDensity(rhoe(i), &oDofMap, dCoef);
        }

        const index_type iIntegRule = iSelectIntegrationRule(pe);
        const index_type iNumGauss = iGetNumGaussPoints(iIntegRule);
        // Note: We have three nodes and only two Gauss points.
        // For that reason we are imposing the same average mass flux for all nodes.
        T mdotz{0.};

        for (index_type i = 1; i <= iNumGauss; ++i) {
            const index_type idxGauss = iGetGaussPointIndex1D(i, iIntegRule);
            const doublereal dx = rgGaussPntDat[idxGauss].detJr * dGetGaussWeight(i, iIntegRule);
            const Vector<doublereal, iNumNodes>& N = rgGaussPntDat[idxGauss].N;
            const Matrix<doublereal, 2, iNumNodes>& B = rgGaussPntDat[idxGauss].B;

            T h = Dot(N, he);

            pGeometry->GetNonNegativeClearance(h, h);

            T dp_dz = Dot(B.GetRow(2), pe);

            const doublereal eta = Dot(N, etae);
            const doublereal rho = Dot(N, rhoe);

            T Uz{0.};

            for (index_type l = 1; l <= iNumNodes; ++l) {
                Uz += N(l) * Ue(l)(2);
            }

            mdotz -= (rho * dx / iNumNodesOutletBound) * h * (Uz - h * h * dp_dz / (12. * eta));
        }

        HYDRO_TRACE("mdotz(" << sref << ")=" << dGetValue(mdotz) * iNumNodesOutletBound << std::endl);

        for (index_type i = 0; i < iNumNodesOutletBound; ++i) {
            HYDRO_ASSERT(rgNodes[pNodeIndexOutletBound[i]]->bIsNodeType(HydroNode::COUPLED_NODE));
            CHECK_NUM_COLS_WORK_SPACE(this, func, mdotz, rgNodes[pNodeIndexOutletBound[i]]->iGetFirstEquationIndex(func));
            WorkVec.AddItem(rgNodes[pNodeIndexOutletBound[i]]->iGetFirstEquationIndex(func), mdotz);
        }
    }

    void QuadFeIso9MassFlowZ::Initialize()
    {
        Matrix<doublereal, iNumNodes, 2> xe;
        Matrix<doublereal, 2, iNumNodes> dN_drv;

        NodePositionMatrix(xe);

        for (index_type k = 0; k < iGetNumIntegrationRules(); ++k) {
            const index_type iNumGauss = iGetNumGaussPoints(k);
            for (index_type i = 1; i <= iNumGauss; ++i) {
                const index_type idx = iGetGaussPointIndex1D(i, k);
                const doublereal r = dGetGaussPos(i, k);

                PressureInterpolMatrix(rgGaussPntDat[idx].N, r, sref);
                PressureInterpolMatrixDer(dN_drv, r, sref);
                PressureGradInterpolMatrix(rgGaussPntDat[idx].B, rgGaussPntDat[idx].detJ, xe, r, sref);

                doublereal ds = 0.;

                for (index_type j = 1; j <= 2; ++j) {
                    const doublereal dNxe_dr = Dot(dN_drv.GetRow(1), xe.GetCol(j));
                    ds += dNxe_dr * dNxe_dr;
                }

                rgGaussPntDat[idx].detJr = sqrt(ds) * sref;
            }
        }
    }

    ThermalFluidModel::ThermalFluidModel(doublereal T0,
                                         doublereal rho0,
                                         doublereal eta0,
                                         doublereal beta)
        :eType(ISOTHERMAL),
         T0(T0),
         cp0(0.),
         rho0(rho0),
         eta0(eta0),
         beta(beta),
         lambda0(0.),
         alphalambda(1.),
         Aeta2_Aeta3(0.),
         Aeta3(std::numeric_limits<doublereal>::max()),
         Ac1(0.),
         Ac2(0.),
         Ac3(0.),
         Ac4(0.),
         Ac5(0.),
         Alambda(0.) {
    }

    void ThermalFluidModel::ParseInput(DataManager* pDM, MBDynParser& HP, const HydroRootElement* pParent)
    {
        if (HP.IsKeyWord("isothermal")) {
            eType = ISOTHERMAL;
        } else if (HP.IsKeyWord("non" "isothermal")) {
            eType = NON_ISOTHERMAL;
        }

        if (HP.IsKeyWord("Aeta2")) {
            if (T0 < 0.) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pParent->GetLabel()
                            << "): illegal reference temperature at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal Aeta2 = HP.GetReal();

            if (!HP.IsKeyWord("Aeta3")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pParent->GetLabel()
                            << "): keyword \"Aeta3\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            Aeta3 = HP.GetReal();
            Aeta2_Aeta3 = Aeta2 / Aeta3;
        }

        if (HP.IsKeyWord("specific" "heat" "capacity") || HP.IsKeyWord("reference" "specific" "heat" "capacity")) {
            cp0 = HP.GetReal();
        }

        if (HP.IsKeyWord("Ac1")) {
            Ac1 = HP.GetReal();
        }

        if (HP.IsKeyWord("Ac2")) {
            Ac2 = HP.GetReal();
        }

        if (HP.IsKeyWord("Ac3")) {
            Ac3 = HP.GetReal();
        }

        if (HP.IsKeyWord("Ac4")) {
            Ac4 = HP.GetReal();
        }

        if (HP.IsKeyWord("Ac5")) {
            Ac5 = HP.GetReal();
        }

        if (HP.IsKeyWord("thermal" "conductivity")) {
            HP.IsKeyWord("liquid");

            lambda0 = HP.GetReal();

            if (HP.IsKeyWord("vapour")) {
                alphalambda = HP.GetReal() / lambda0;
            }

            if (lambda0 == 0.) {
                alphalambda = 0.;
            }
        }

        if (HP.IsKeyWord("Alambda")) {
            Alambda = HP.GetReal();
        }

        if (eType != ISOTHERMAL && !bValid()) {
            silent_cerr("hydrodynamic plain bearing2(" << pParent->GetLabel()
                        << "): \"density\", \"specific heat capacity\" and \"thermal conductivity\" "
                        "are required for non isothermal fluid models at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    template <typename U>
    U ThermalFluidModel::GetViscosityLiquid(const U& T) const
    {
        // Dirk Bartel 2009 equation (6-11)
        return eta0 * exp(Aeta2_Aeta3 * (T0 - T) / (Aeta3 + T - T0));
    }

    template <typename U>
    U ThermalFluidModel::GetSpecHeatPerVolume(const U& p, const U& T, HeatCapacityType eType) const
    {
        // alpha = 1 for rho * cp(T)
        // alpha = 0.5 for rho * integrate(cp(T), T, 0, T) / T

        static const doublereal alpha[] = {1., 0.5};

        static_assert(SPEC_HEAT_TRUE == 0, "index does not match");
        static_assert(SPEC_HEAT_AVERAGED == 1, "index does not match");

        // Dirk Bartel 2009 equation (6-8)
        return rho0 * cp0 * (1 + (Ac1 * p) / (1 + Ac2 * p)) * (1 + Ac3 * (1 + Ac4 * p + Ac5 * p * p) * (alpha[eType] * T - T0));
    }

    template <typename U>
    U ThermalFluidModel::GetDensityLiquid(const U& T, U* drho_dT) const
    {
        // Dirk Bartel 2009 equation (6-1)
        if (drho_dT) {
            *drho_dT = -rho0 * beta;
        }

        return rho0 * (1 - beta * (T - T0));
    }

    template <typename U>
    U ThermalFluidModel::GetThermalConductivityLiquid(const U& T) const
    {
        // Dirk Bartel 2009 equation (6-5)
        return lambda0 - Alambda * (T - T0);
    }

    template <typename U>
    U ThermalFluidModel::GetSpecificHeatLiquid(const U& p, const U& T, HeatCapacityType eType) const {
        return GetSpecHeatPerVolume(p, T, eType) / GetDensityLiquid(T);
    }

    template <typename U>
    U ThermalFluidModel::GetThermalConductivityMixture(const U& T, const U& rho) const {
        const U rholiq = GetDensityLiquid(T);
        const U lambdaliq = GetThermalConductivityLiquid(T);

        return lambdaliq * (alphalambda + (1. - alphalambda) * rho / rholiq);
    }

    bool ThermalFluidModel::bValid() const
    {
        return rho0 > 0. &&
            cp0 > 0. &&
            lambda0 >= 0. &&
            alphalambda >= 0. &&
            alphalambda <= 1. &&
            std::isfinite(Aeta2_Aeta3);
    }

    HydroFluid::HydroFluid(doublereal pc, const ThermalFluidModel& oThermModel)
        :pc(pc),
         oThermModel(oThermModel)
    {

    }

    HydroFluid::~HydroFluid()
    {

    }

    doublereal HydroFluid::dGetRefPressure() const
    {
        return pc;
    }

    doublereal HydroFluid::dGetRefDensity() const
    {
        return oThermModel.dGetRefDensity();
    }

    HydroIncompressibleFluid::HydroIncompressibleFluid(doublereal pc, const ThermalFluidModel& oThermModel)
        :HydroFluid(pc, oThermModel)
    {

    }

    HydroIncompressibleFluid::~HydroIncompressibleFluid()
    {

    }

    void HydroIncompressibleFluid::GetDensity(const doublereal& p, const doublereal& T, doublereal& rho, doublereal* drho_dp, doublereal* drho_dT) const
    {
        GetDensityTpl(p, T, rho, drho_dp, drho_dT);
    }

    void HydroIncompressibleFluid::GetDensity(const Gradient<0>& p, const Gradient<0>& T, Gradient<0>& rho, Gradient<0>* drho_dp, Gradient<0>* drho_dT) const
    {
        GetDensityTpl(p, T, rho, drho_dp, drho_dT);
    }

    void HydroIncompressibleFluid::GetPressure(const doublereal& rho, const doublereal& T, doublereal& p, doublereal* dp_drho, doublereal* dp_dT) const
    {
        GetPressureTpl(rho, T, p, dp_drho, dp_dT);
    }

    void HydroIncompressibleFluid::GetPressure(const Gradient<0>& rho, const Gradient<0>& T,  Gradient<0>& p, Gradient<0>* dp_drho, Gradient<0>* dp_dT) const
    {
        GetPressureTpl(rho, T, p, dp_drho, dp_dT);
    }

    void HydroIncompressibleFluid::GetViscosity(const doublereal& rho, const doublereal& T, doublereal& eta) const
    {
        GetViscosityTpl(rho, T, eta);
    }

    void HydroIncompressibleFluid::GetViscosity(const Gradient<0>& rho, const Gradient<0>& T, Gradient<0>& eta) const
    {
        GetViscosityTpl(rho, T, eta);
    }

    void
    HydroIncompressibleFluid::ThetaToPhysical(const std::array<doublereal, iNumDof>& Theta,
                                              const std::array<doublereal, iNumDof>& dTheta_dt,
                                              const doublereal& T,
                                              const doublereal& dT_dt,
                                              doublereal& p,
                                              doublereal& dp_dt,
                                              doublereal& rho,
                                              doublereal& drho_dt) const
    {
        ThetaToPhysicalTpl(Theta[0], dTheta_dt[0], T, dT_dt, p, dp_dt, rho, drho_dt);
    }

    void
    HydroIncompressibleFluid::ThetaToPhysical(const std::array<Gradient<0>, iNumDof>& Theta,
                                              const std::array<Gradient<0>, iNumDof>& dTheta_dt,
                                              const Gradient<0>& T,
                                              const Gradient<0>& dT_dt,
                                              Gradient<0>& p,
                                              Gradient<0>& dp_dt,
                                              Gradient<0>& rho,
                                              Gradient<0>& drho_dt) const
    {
        ThetaToPhysicalTpl(Theta[0], dTheta_dt[0], T, dT_dt, p, dp_dt, rho, drho_dt);
    }

    template <typename G>
    inline void
    HydroIncompressibleFluid::ThetaToPhysicalTpl(const G& Theta,
                                                 const G& dTheta_dt,
                                                 const G& T,
                                                 const G& dT_dt,
                                                 G& p,
                                                 G& dp_dt,
                                                 G& rho,
                                                 G& drho_dt) const
    {
        G drho_dT;
        rho = oThermModel.GetDensityLiquid(T, &drho_dT);
        drho_dt = drho_dT * dT_dt;
        p = Theta;
        dp_dt = dTheta_dt;
    }

    doublereal HydroIncompressibleFluid::GetTheta0(index_type iDofIndex) const
    {
        HYDRO_ASSERT(iDofIndex == 0);

        return pc;
    }

    HydroFluid::CavitationState HydroIncompressibleFluid::Cavitation(doublereal& p, doublereal* dp_dt) const
    {
        return CavitationTpl(p, dp_dt);
    }

    HydroFluid::CavitationState HydroIncompressibleFluid::Cavitation(Gradient<0>& p, Gradient<0>* dp_dt) const
    {
        return CavitationTpl(p, dp_dt);
    }

    template <typename G> inline void
    HydroIncompressibleFluid::GetDensityTpl(const G& p, const G& T, G& rho, G* drho_dp, G* drho_dT) const
    {
        rho = oThermModel.GetDensityLiquid(T, drho_dT);

        if (drho_dp) {
            *drho_dp = 0.;
        }
    }

    template <typename G> inline void
    HydroIncompressibleFluid::GetPressureTpl(const G& rho, const G& T, G& p, G* dp_drho, G* dp_dT) const
    {
        throw ErrNotUnique(MBDYN_EXCEPT_ARGS);
    }

    template <typename G> inline void
    HydroIncompressibleFluid::GetViscosityTpl(const G& rho, const G& T, G& eta) const
    {
        eta = oThermModel.GetViscosityLiquid(T);
    }

    template <typename T> inline HydroFluid::CavitationState
    HydroIncompressibleFluid::CavitationTpl(T& p, T* dp_dt) const
    {
        if (p < pc) {
            p = pc;

            if (dp_dt) {
                *dp_dt = 0.;
            }

            return CAVITATION_REGION;
        }

        return FULL_FILM_REGION;
    }

    HydroFluid::HydraulicType HydroIncompressibleFluid::GetHydraulicType() const
    {
        return INCOMPRESSIBLE;
    }

    LinearCompressibleFluid::LinearCompressibleFluid(doublereal etavap_etaliq, doublereal beta, const doublereal pc, HydraulicType type, const ThermalFluidModel& oThermModel)
        :HydroFluid(pc, oThermModel),
         etavap_etaliq(etavap_etaliq),
         beta(beta),
         type(type)
    {

    }

    LinearCompressibleFluid::~LinearCompressibleFluid()
    {

    }

    void LinearCompressibleFluid::GetDensity(const doublereal& p, const doublereal& T, doublereal& rho, doublereal* drho_dp, doublereal* drho_dT) const
    {
        GetDensityTpl(p, T, rho, drho_dp, drho_dT);
    }

    void LinearCompressibleFluid::GetDensity(const Gradient<0>& p, const Gradient<0>& T, Gradient<0>& rho, Gradient<0>* drho_dp, Gradient<0>* drho_dT) const
    {
        GetDensityTpl(p, T, rho, drho_dp, drho_dT);
    }

    void LinearCompressibleFluid::GetPressure(const doublereal& rho, const doublereal& T, doublereal& p, doublereal* dp_drho, doublereal* dp_dT) const
    {
        GetPressureTpl(rho, T, p, dp_drho, dp_dT);
    }

    void LinearCompressibleFluid::GetPressure(const Gradient<0>& rho, const Gradient<0>& T, Gradient<0>& p, Gradient<0>* dp_drho, Gradient<0>* dp_dT) const
    {
        GetPressureTpl(rho, T, p, dp_drho, dp_dT);
    }

    void LinearCompressibleFluid::GetViscosity(const doublereal& rho, const doublereal& T, doublereal& eta) const
    {
        GetViscosityTpl(rho, T, eta);
    }

    void LinearCompressibleFluid::GetViscosity(const Gradient<0>& rho, const Gradient<0>& T, Gradient<0>& eta) const
    {
        GetViscosityTpl(rho, T, eta);
    }

    void
    LinearCompressibleFluid::ThetaToPhysical(const std::array<doublereal, iNumDof>& Theta,
                                             const std::array<doublereal, iNumDof>& dTheta_dt,
                                             const doublereal& T,
                                             const doublereal& dT_dt,
                                             doublereal& p,
                                             doublereal& dp_dt,
                                             doublereal& rho,
                                             doublereal& drho_dt) const
    {
        ThetaToPhysicalTpl(Theta, dTheta_dt, T, dT_dt, p, dp_dt, rho, drho_dt);
    }

    void
    LinearCompressibleFluid::ThetaToPhysical(const std::array<Gradient<0>, iNumDof>& Theta,
                                             const std::array<Gradient<0>, iNumDof>& dTheta_dt,
                                             const Gradient<0>& T,
                                             const Gradient<0>& dT_dt,
                                             Gradient<0>& p,
                                             Gradient<0>& dp_dt,
                                             Gradient<0>& rho,
                                             Gradient<0>& drho_dt) const
    {
        ThetaToPhysicalTpl(Theta, dTheta_dt, T, dT_dt, p, dp_dt, rho, drho_dt);
    }

    template <typename G>
    inline void
    LinearCompressibleFluid::ThetaToPhysicalTpl(const std::array<G, iNumDof>& Theta,
                                                const std::array<G, iNumDof>& dTheta_dt,
                                                const G& T,
                                                const G& dT_dt,
                                                G& p,
                                                G& dp_dt,
                                                G& rho,
                                                G& drho_dt) const
    {
        static_assert(iNumDof >= 2, "number of degrees of freedome does not match");

        G drhoc_dT;
        const G rhoc = oThermModel.GetDensityLiquid(T, &drhoc_dT);

        p = pc + Theta[0];
        dp_dt = dTheta_dt[0];
        rho = rhoc * Theta[1];
        drho_dt = drhoc_dT * dT_dt * Theta[1] + rhoc * dTheta_dt[1];

        HYDRO_ASSERT(std::isfinite(grad::dGetValue(p)));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(dp_dt)));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(rho)));
        HYDRO_ASSERT(std::isfinite(grad::dGetValue(drho_dt)));
    }

    doublereal LinearCompressibleFluid::GetTheta0(index_type iDofIndex) const
    {
        switch(iDofIndex) {
        case 0:
            return 0.;

        case 1:
            return 1.;

        default:
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    HydroFluid::CavitationState LinearCompressibleFluid::Cavitation(doublereal& p, doublereal* dp_dt) const
    {
        return CavitationTpl(p, dp_dt);
    }

    HydroFluid::CavitationState LinearCompressibleFluid::Cavitation(Gradient<0>& p, Gradient<0>* dp_dt) const
    {
        return CavitationTpl(p, dp_dt);
    }

    template <typename G> inline void
    LinearCompressibleFluid::GetDensityTpl(const G& p, const G& T, G& rho, G* drho_dp, G* drho_dT) const
    {
        const G rhoc = oThermModel.GetDensityLiquid(T, drho_dT);

        if (p >= pc) {
            // incompressible fluid
            rho = rhoc;

            if (drho_dp) {
                *drho_dp = 0.;
            }
        } else {
            // Since pressure below the cavity pressure is not possible,
            // the pressure boundary condition will be interpreted as density boundary condition
            rho = p / pc * rhoc;

            if (drho_dp) {
                *drho_dp = rhoc / pc;
            }

            if (drho_dT) {
                *drho_dT *= p / pc;
            }
        }
    }

    template <typename G> inline void
    LinearCompressibleFluid::GetPressureTpl(const G& rho, const G& T, G& p, G* dp_drho, G* dp_dT) const
    {
        p = pc;

        if (dp_drho) {
            *dp_drho = 0.;
        }

        if (dp_dT) {
            *dp_dT = 0.;
        }
    }

    template <typename G> inline void
    LinearCompressibleFluid::GetViscosityTpl(const G& rho, const G& T, G& eta) const
    {
        const G rholiq = oThermModel.GetDensityLiquid(T);
        eta = oThermModel.GetViscosityLiquid(T);

        if (rho < rholiq) {
            eta *= ((1. - etavap_etaliq) * rho / rholiq + etavap_etaliq);
        }
    }

    template <typename T> inline HydroFluid::CavitationState
    LinearCompressibleFluid::CavitationTpl(T& p, T* dp_dt) const
    {
        if (p < pc) {
            p = pc;

            if (dp_dt) {
                *dp_dt = 0.;
            }

            return CAVITATION_REGION;
        }

        return FULL_FILM_REGION;
    }

    HydroFluid::HydraulicType LinearCompressibleFluid::GetHydraulicType() const
    {
        return type;
    }

    FluidStateBoundaryCond::FluidStateBoundaryCond(const HydroFluid* pFluid,
                                                   Type eType,
                                                   ExtrapMethod eExtrapMethod,
                                                   unsigned uNodeMask,
                                                   std::unique_ptr<DriveCaller>&& pTemp)
        :pFluid(pFluid),
         eType(eType),
         eExtrapMethod(eExtrapMethod),
         uNodeMask(uNodeMask),
         pTempDrv(std::move(pTemp))
    {

    }

    FluidStateBoundaryCond::~FluidStateBoundaryCond()
    {

    }

    std::unique_ptr<FluidStateBoundaryCond> FluidStateBoundaryCond::Read(MBDynParser& HP, const HydroRootElement* pRoot)
    {
        Type eBoundCondType = BC_PRESSURE; // default boundary condition for compatibility reasons
        unsigned uNodeMask = 0;

        if (HP.IsKeyWord("pressure")) {
            eBoundCondType = BC_PRESSURE;
            uNodeMask |= Node2D::HYDRAULIC_NODE;
        } else if (HP.IsKeyWord("density")) {
            eBoundCondType = BC_DENSITY;
            uNodeMask |= Node2D::HYDRAULIC_NODE;
        } else if (HP.IsKeyWord("filling" "ratio")) {
            eBoundCondType = BC_FILLING_RATIO;
            uNodeMask |= Node2D::HYDRAULIC_NODE;
        } else if (HP.IsKeyWord("temperature")) {
            uNodeMask |= Node2D::THERMAL_NODE;
        } else {
            uNodeMask |= Node2D::HYDRAULIC_NODE; // default boundary condition for backward compatibility
        }

        std::unique_ptr<DriveCaller> pPressDensFill{(uNodeMask & Node2D::HYDRAULIC_NODE)
                ? HP.GetDriveCaller()
                : new NullDriveCaller};

        std::unique_ptr<DriveCaller> pTemp;

        if ((uNodeMask & Node2D::THERMAL_NODE) || HP.IsKeyWord("temperature")) {
            uNodeMask |= Node2D::THERMAL_NODE;
            pTemp.reset(HP.GetDriveCaller());
        } else {
            pTemp.reset(new NullDriveCaller);
        }

        ExtrapMethod eExtrapMethod = EX_INLET;

        if (HP.IsKeyWord("type")) {
            if (HP.IsKeyWord("outlet")) {
                eExtrapMethod = EX_OUTLET;
            } else if (HP.IsKeyWord("inlet")) {
                eExtrapMethod = EX_INLET;
            } else {
                silent_cerr("hydrodynamic_plain_bearing2("
                            << pRoot->GetLabel()
                            << "): keyword \"outlet\" or \"inlet\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        switch (eBoundCondType) {
        case BC_PRESSURE:
        case BC_DENSITY:
            return std::unique_ptr<FluidStateBoundaryCond>{new FluidStateFunction(pRoot->pGetFluid(),
                                                                                  eBoundCondType,
                                                                                  eExtrapMethod,
                                                                                  uNodeMask,
                                                                                  std::move(pPressDensFill),
                                                                                  std::move(pTemp))};

        case BC_FILLING_RATIO:
            return std::unique_ptr<FluidStateBoundaryCond>{new FillingRatioFunction(pRoot->pGetFluid(),
                                                                                    eBoundCondType,
                                                                                    eExtrapMethod,
                                                                                    uNodeMask,
                                                                                    std::move(pPressDensFill),
                                                                                    std::move(pTemp))};

        default:
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    bool FluidStateBoundaryCond::bNeedClearance() const
    {
        return eType == BC_FILLING_RATIO;
    }

    inline doublereal FluidStateBoundaryCond::dGetTemperature() const
    {
        return T;
    }

    inline doublereal FluidStateBoundaryCond::dGetTemperatureDerTime() const
    {
        return dT_dt;
    }

    void FluidStateBoundaryCond::Update()
    {
        T = pTempDrv->dGet();
        dT_dt = pTempDrv->dGetP();
    }

    FluidStateFunction::FluidStateFunction(const HydroFluid* pFluid,
                                           Type eType,
                                           ExtrapMethod eExtrapMethod,
                                           unsigned uNodeMask,
                                           std::unique_ptr<DriveCaller>&& pPressDens,
                                           std::unique_ptr<DriveCaller>&& pTemp)
        :FluidStateBoundaryCond(pFluid, eType, eExtrapMethod, uNodeMask, std::move(pTemp)),
         pPressDensDrv(std::move(pPressDens)),
         p(0.),
         dp_dt(0.),
         rho(0.),
         drho_dt(0.)
    {


    }

    doublereal FluidStateFunction::dGetPressure(doublereal) const
    {
        return p;
    }

    doublereal FluidStateFunction::dGetPressureDerTime(doublereal, doublereal) const
    {
        return dp_dt;
    }

    doublereal FluidStateFunction::dGetDensity(doublereal) const
    {
        return rho;
    }

    doublereal FluidStateFunction::dGetDensityDerTime(doublereal, doublereal) const
    {
        return drho_dt;
    }

    FluidStateFunction::~FluidStateFunction()
    {

    }

    void FluidStateFunction::Update() {
        FluidStateBoundaryCond::Update();

        const doublereal T = dGetTemperature();
        const doublereal dT_dt = dGetTemperatureDerTime();

        switch (GetType()) {
        case BC_PRESSURE:
            p = pPressDensDrv->dGet();

            if (pPressDensDrv->bIsDifferentiable()) {
                dp_dt = pPressDensDrv->dGetP();
            }

            doublereal drho_dp, drho_dT;

            pGetFluid()->GetDensity(p, T, rho, &drho_dp, &drho_dT);

            pGetFluid()->Cavitation(p, &dp_dt);

            drho_dt = drho_dp * dp_dt + drho_dT * dT_dt;
            break;

        case BC_DENSITY:
            rho = pPressDensDrv->dGet();

            if (pPressDensDrv->bIsDifferentiable()) {
                drho_dt = pPressDensDrv->dGetP();
            }

            doublereal dp_drho, dp_dT;

            pGetFluid()->GetPressure(rho, T, p, &dp_drho, &dp_dT);

            dp_dt = dp_drho * drho_dt + dp_dT * dT_dt;
            break;

        default:
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }
    }

    FillingRatioFunction::FillingRatioFunction(const HydroFluid* pFluid,
                                               Type eType,
                                               ExtrapMethod eExtrapMethod,
                                               unsigned uNodeMask,
                                               std::unique_ptr<DriveCaller>&& pFill,
                                               std::unique_ptr<DriveCaller>&& pTemp)
        :FluidStateBoundaryCond(pFluid, eType, eExtrapMethod, uNodeMask, std::move(pTemp)),
         pFillRatioDrv(std::move(pFill)),
         h0(0.),
         dh0_dt(0.),
         rho(0.),
         drho_dt(0.)
    {
        HYDRO_ASSERT(eType == BC_FILLING_RATIO);
    }

    FillingRatioFunction::~FillingRatioFunction()
    {

    }

    doublereal FillingRatioFunction::dGetPressure(doublereal h) const
    {
        const doublereal rho = dGetDensity(h);
        const doublereal T = dGetTemperature();
        doublereal p;

        pGetFluid()->GetPressure(rho, T, p);

        return p;
    }

    doublereal FillingRatioFunction::dGetPressureDerTime(doublereal h, doublereal dh_dt) const
    {
        doublereal p, dp_drho, dp_dT;
        const doublereal T = dGetTemperature();
        const doublereal dT_dt = dGetTemperatureDerTime();

        pGetFluid()->GetPressure(rho, T, p, &dp_drho, &dp_dT);

        return dp_drho * drho_dt + dp_dT * dT_dt;
    }

    doublereal FillingRatioFunction::dGetDensity(doublereal h) const
    {
        doublereal alpha = h0 / h;

        if (alpha > 1.) {
            alpha = 1.;
        }

        return alpha * rho;
    }

    doublereal FillingRatioFunction::dGetDensityDerTime(doublereal h, doublereal dh_dt) const
    {
        doublereal alpha = h0 / h;
        doublereal dalpha_dt;

        if (alpha > 1.) {
            alpha = 1.;
            dalpha_dt = 0.;
        } else {
            dalpha_dt = dh0_dt / h - h0 / (h * h) * dh_dt;
        }

        return dalpha_dt * rho + alpha * drho_dt;
    }

    void FillingRatioFunction::Update()
    {
        FluidStateBoundaryCond::Update();

        h0 = pFillRatioDrv->dGet();

        if (pFillRatioDrv->bIsDifferentiable()) {
            dh0_dt = pFillRatioDrv->dGetP();
        }

        const doublereal T = dGetTemperature();
        const doublereal dT_dt = dGetTemperatureDerTime();
        const doublereal p = pGetFluid()->dGetRefPressure();
        doublereal drho_dT;

        pGetFluid()->GetDensity(p, T, rho, nullptr, &drho_dT);

        drho_dt = drho_dT * dT_dt;
    }

    PressureCouplingCond::PressureCouplingCond(integer iLabel, std::unique_ptr<Geometry2D>&& pGeometry)
        :iLabel(iLabel),
         pGeometry(std::move(pGeometry))
    {

    }

    PressureCouplingCond::~PressureCouplingCond()
    {

    }

    std::unique_ptr<PressureCouplingMaster> PressureCouplingCond::Read(integer iLabel,
                                                                       HydroRootElement* pRoot,
                                                                       DataManager* pDM,
                                                                       MBDynParser& HP)
    {
        std::unique_ptr<Geometry2D> pGeometry(Geometry2D::Read(pRoot, HP));

        if (!HP.IsKeyWord("pressure" "node")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"pressure node\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        PressureNode* pHydroNode = pDM->ReadNode<PressureNode, Node::HYDRAULIC>(HP);

        ThermalNode* pThermalNode = nullptr;

        const HydroFluid* pFluid = pRoot->pGetFluid();

        if (HP.IsKeyWord("thermal" "node")) {
            if (pFluid->GetHydraulicType() == HydroFluid::INCOMPRESSIBLE) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pRoot->GetLabel()
                            << "): thermal coupling not supported for incompressible fluids at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
            pThermalNode = pDM->ReadNode<ThermalNode, Node::THERMAL>(HP);
        } else if (pFluid->GetThermalType() != HydroFluid::ISOTHERMAL) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pRoot->GetLabel()
                        << "): keyword \"thermal node\" expected at line "
                        << HP.GetLineData() << std::endl);

            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return std::unique_ptr<PressureCouplingMaster>{new PressureCouplingMaster(iLabel,
                                                                                  pHydroNode,
                                                                                  pThermalNode,
                                                                                  std::move(pGeometry))};
    }

    PressureCouplingMaster::PressureCouplingMaster(integer iLabel,
                                                   PressureNode* pHydroNode,
                                                   ThermalNode* pThermalNode,
                                                   std::unique_ptr<Geometry2D>&& pGeometry)
        :PressureCouplingCond(iLabel, std::move(pGeometry)),
         pHydroNode(pHydroNode),
         pThermalNode(pThermalNode),
         iNumNodes(0)
    {

    }

    PressureCouplingMaster::~PressureCouplingMaster()
    {

    }

    PressureNode* PressureCouplingMaster::pGetNode() const {
        return pHydroNode;
    }

    ThermalNode* PressureCouplingMaster::pGetThermalNode() const {
        return pThermalNode;
    }

    void PressureCouplingMaster::AddNode(HydroNode* pNode) {
        ++iNumNodes;
    }

    integer PressureCouplingMaster::iGetNumNodes() const {
        return iNumNodes;
    }

    std::unique_ptr<PressureCouplingSlave> PressureCouplingMaster::Clone(const Vector<doublereal, 2>& x)
    {
        return std::unique_ptr<PressureCouplingSlave>{new PressureCouplingSlave(this, pGetGeometry()->Clone(x))};
    }

    PressureCouplingSlave::PressureCouplingSlave(PressureCouplingMaster* pMaster, std::unique_ptr<Geometry2D>&& pGeometry)
        :PressureCouplingCond(pMaster->iGetLabel(), std::move(pGeometry)),
         pMaster(pMaster)
    {

    }

    PressureNode* PressureCouplingSlave::pGetNode() const
    {
        return pMaster->pGetNode();
    }

    ThermalNode* PressureCouplingSlave::pGetThermalNode() const
    {
        return pMaster->pGetThermalNode();
    }

    void PressureCouplingSlave::AddNode(HydroNode* pNode)
    {
        pMaster->AddNode(pNode);
    }

    integer PressureCouplingSlave::iGetNumNodes() const
    {
        return pMaster->iGetNumNodes();
    }

    HydroMesh::HydroMesh(HydroRootElement* pParent)
        :pGeometry(nullptr),
         pCompliance(nullptr),
         bUseOutletAxial(false),
         bThermalModel(pParent->pGetFluid()->GetThermalType() != HydroFluid::ISOTHERMAL),
         pParent(pParent)
    {
    }

    HydroMesh::~HydroMesh()
    {
    }

    void HydroMesh::ParseGeometry(DataManager* pDM, MBDynParser& HP)
    {
        if (!HP.IsKeyWord("geometry")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"geometry\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!HP.IsKeyWord("cylindrical")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"cylindrical\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (!HP.IsKeyWord("mesh" "position")) {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"mesh position\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("at" "shaft")) {
            pGeometry.reset(new CylindricalMeshAtShaft(pGetParent()));
        } else if (HP.IsKeyWord("at" "bearing")) {
            pGeometry.reset(new CylindricalMeshAtBearing(pGetParent()));
        } else {
            silent_cerr("hydrodynamic plain bearing2("
                        << pGetParent()->GetLabel()
                        << "): keyword \"at shaft\" or \"at bearing\" expected at line "
                        << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        pGeometry->ParseInput(pDM, HP);
    }

    void HydroMesh::ParseBoundaryCond(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("pressure" "coupling" "conditions" "axial")) {
            for (index_type i = 0; i < iNumCouplingAxial; ++i) {
                rgOutletAxial[i].pExtHydroNode = pDM->ReadNode<PressureNode, Node::HYDRAULIC>(HP);
            }

            if (bThermalModel) {
                if (!HP.IsKeyWord("thermal" "coupling" "conditions" "axial")) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pGetParent()->GetLabel()
                                << "): keyword \"thermal coupling condition\" expected at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                for (index_type i = 0; i < iNumCouplingAxial; ++i) {
                    rgOutletAxial[i].pExtThermNode = pDM->ReadNode<ThermalNode, Node::THERMAL>(HP);
                }
            }

            bUseOutletAxial = true;
        } else {
            if (!HP.IsKeyWord("boundary" "conditions")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"boundary conditions\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            rgBoundaryCond.reserve(iNumCouplingAxial);

            for (integer i = 0; i < iNumCouplingAxial; ++i) {
                rgBoundaryCond.push_back(FluidStateBoundaryCond::Read(HP, pGetParent()));
            }
        }
    }

    void HydroMesh::ParsePressureCouplingCond(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("pressure" "coupling" "conditions" "radial")) {
            const integer iNumCouplings = HP.GetInt();

            for (integer i = 0; i < iNumCouplings; ++i) {
                std::unique_ptr<PressureCouplingMaster> pCoupling{PressureCouplingCond::Read(i + 1,
                                                                                             pGetParent(),
                                                                                             pDM,
                                                                                             HP)};

                Vector<doublereal, 2> x = pCoupling->pGetGeometry()->GetPosition();

                const doublereal c = 2 * M_PI * pGeometry->dGetMeshRadius();

                if (x(1) < 0. || x(1) > c) {
                    // x(2) may be outside the cylindrical bearing
                    silent_cerr("hydrodynamic plain bearing2("
                                << pGetParent()->GetLabel()
                                << "): pressure coupling condition(" << rgCouplingCond.size()
                                << ") is outside the bearing are at line "
                                << HP.GetLineData() << std::endl);

                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (x(1) >= 0.5 * c) {
                    x(1) -= c;
                } else {
                    x(1) += c;
                }

                rgCouplingCond.push_back(pCoupling->Clone(x));
                rgCouplingCond.push_back(std::move(pCoupling));
            }
        }
    }

    void HydroMesh::ParseContactModel(DataManager* pDM, MBDynParser& HP)
    {
        enum {
            FRIC_MOD_UNKNOWN,
            FRIC_MOD_NONE,
            FRIC_MOD_COULOMB,
            FRIC_MOD_LUGRE
        } eFrictionModel = FRIC_MOD_UNKNOWN;

        enum {
            CONT_MOD_NONE,
            CONT_MOD_PENALTY,
            CONT_MOD_GREENWOOD_TRIPP
        } eContactModel = CONT_MOD_NONE;

        if (HP.IsKeyWord("contact" "model")) {
            if (HP.IsKeyWord("greenwood" "williamson")
                || HP.IsKeyWord("greenwood" "williamson" "coulomb")) {
                pedantic_cerr("warning: contact model greenwood williamson is obsolete, use greenwood tripp instead at line "
                              << HP.GetLineData() << std::endl);

                eContactModel = CONT_MOD_GREENWOOD_TRIPP; // provided for compatibility
                eFrictionModel = FRIC_MOD_COULOMB;
            } else if (HP.IsKeyWord("greenwood" "williamson" "lugre")) {
                pedantic_cerr("warning: contact model greenwood williamson lugre is obsolete, use greenwood tripp instead at line "
                              << HP.GetLineData() << std::endl);

                eContactModel = CONT_MOD_GREENWOOD_TRIPP; // provided for compatibility
                eFrictionModel = FRIC_MOD_LUGRE;
            } else if (HP.IsKeyWord("penalty")) {
                eContactModel = CONT_MOD_PENALTY;
            } else if (HP.IsKeyWord("greenwood" "tripp")){
                eContactModel = CONT_MOD_GREENWOOD_TRIPP;
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"greenwood tripp\" or \"penalty\" expected at line "
                            << HP.GetLineData() << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        switch (eContactModel) {
        case CONT_MOD_PENALTY:
            pContact.reset(new PenaltyCM(this, pGetGeometry()->dGetReferenceClearance()));
            break;

        case CONT_MOD_GREENWOOD_TRIPP:
            pContact.reset(new GreenwoodTrippCM(this));
            break;

        case CONT_MOD_NONE:
            break;

        default:
            HYDRO_ASSERT(false);
        }

        if (eContactModel != CONT_MOD_NONE) {
            pContact->ParseInput(HP);

            if (eFrictionModel == FRIC_MOD_UNKNOWN) {
                if (HP.IsKeyWord("friction" "model")) {
                    if (HP.IsKeyWord("coulomb")) {
                        eFrictionModel = FRIC_MOD_COULOMB;
                    } else if (HP.IsKeyWord("lugre")) {
                        eFrictionModel = FRIC_MOD_LUGRE;
                    } else if (HP.IsKeyWord("none")) {
                        eFrictionModel = FRIC_MOD_NONE;
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): keyword \"coulomb\", \"lugre\" or \"none\" expected at line "
                                    << HP.GetLineData() << std::endl);

                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                } else {
                    eFrictionModel = FRIC_MOD_NONE;
                }
            }

            switch (eFrictionModel) {
            case FRIC_MOD_NONE:
            case FRIC_MOD_COULOMB:
                pFriction.reset(new CoulombFriction(this));
                break;

            case FRIC_MOD_LUGRE:
                pFriction.reset(new LugreFriction(this));
                break;

            default:
                HYDRO_ASSERT(false);
            }

            if (eFrictionModel != FRIC_MOD_NONE) {
                pFriction->ParseInput(HP);
            } // otherwise use coulomb friction with mu=0
        }
    }

    void HydroMesh::ParseLubricationGrooves(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("lubrication" "grooves")) {
            const integer iNumGrooves = HP.GetInt();
            integer iNumBoundaryCond = iNumGrooves + rgBoundaryCond.size();
            rgGrooves.reserve(2 * iNumGrooves); // periodic condition for cylindrical bearings
            rgBoundaryCond.reserve(iNumBoundaryCond);

            for (integer i = 0; i < iNumGrooves; ++i) {
                std::unique_ptr<LubricationGrooveMaster> pGroove{LubricationGroove::Read(i + 1,
                                                                                         pGetParent(),
                                                                                         pGeometry.get(),
                                                                                         HP)};
                rgBoundaryCond.emplace_back(pGroove->pReleaseBoundaryCond());

                // Take into account the periodic nature of the cylindrical bearing
                Vector<doublereal, 2> x = pGroove->pGetGeometry()->GetPosition();

                const doublereal c = 2 * M_PI * pGeometry->dGetMeshRadius();

                if (x(1) < 0. || x(1) > c) {
                    // x(2) may be outside the cylindrical bearing
                    silent_cerr("hydrodynamic plain bearing2("
                                << pGetParent()->GetLabel()
                                << "): lubrication groove (" << pGroove->iGetLabel()
                                << ") is outside the bearing are at line "
                                << HP.GetLineData() << std::endl);

                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (x(1) >= 0.5 * c) {
                    x(1) -= c;
                } else {
                    x(1) += c;
                }

                rgGrooves.emplace_back(new LubricationGrooveSlave(pGroove.get(), x));
                rgGrooves.push_back(std::move(pGroove));
            }
        }
    }

    void HydroMesh::ParseThermWallBoundCond(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("thermal" "wall" "boundary" "conditions")) {
            if (HP.IsKeyWord("coupled")) {
                pThermWallBoundCond.reset(new ThermWallBoundCond);
                pThermWallBoundCond->ParseInput(pDM, HP, pGetParent());
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"coupled\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }
    }

    std::string HydroMesh::ParseFileName(MBDynParser& HP)
    {
        const char* pszFileName = HP.GetFileName();

        if (!pszFileName || strlen(pszFileName) == 0) {
            silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                        << "): file name expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        return pszFileName;
    }

    void HydroMesh::ParseComplianceModel(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("compliance" "model")) {
            const index_type iNumMatricesMax = ComplianceModelNodalDouble::DEHD_BODY_LAST;
            std::array<Material, iNumMatricesMax> rgMaterials;
            std::array<std::unique_ptr<ComplianceMatrix>, iNumMatricesMax> rgMatrices;
            std::array<const Modal*, iNumMatricesMax> rgModalJoints = {nullptr};
            index_type iNumMatrices = iNumMatricesMax;
            index_type iNumMaterials = iNumMatrices;
            std::string strFileNameModal;
            auto eInterpolOption = ComplianceModelNodalDouble::INT_AXIAL_EXTRAPOLATE;
            ComplianceModel::Type eCompModType = ComplianceModel::COMP_MOD_UNKNOWN;

            if (pContact.get()) {
                for (index_type i = 0; i < iNumMatrices; ++i) {
                    rgMaterials[i] = pContact->GetMaterial(i);
                }
            }

            if (HP.IsKeyWord("elastic" "half" "space")) {
                eCompModType = ComplianceModel::COMP_MOD_NODAL;
                iNumMatrices = 1;

                if (!pContact.get()) {
                    for (index_type i = 0; i < 2; ++i) {
                        rgMaterials[i].ParseInput(i + 1, HP, pGetParent());
                    }
                }

                rgMatrices[0].reset(new ElasticHalfSpace(rgMaterials[0], rgMaterials[1]));
            } else if (HP.IsKeyWord("modal")) {
                eCompModType = ComplianceModel::COMP_MOD_MODAL;
                iNumMatrices = 1;

                if (!pContact.get()) {
                    iNumMaterials = 1;
                    rgMaterials[0].ParseInput(1, HP, pGetParent());
                }

                strFileNameModal = ParseFileName(HP);
            } else if (HP.IsKeyWord("double" "nodal")) {
                eCompModType = ComplianceModel::COMP_MOD_NODAL_DOUBLE;
                std::array<BearingGeometry::Type, ComplianceModelNodalDouble::DEHD_BODY_LAST> rgMeshPos;

                std::fill(std::begin(rgMeshPos), std::end(rgMeshPos), BearingGeometry::CYLINDRICAL_MESH_UNKNOWN);

                iNumMaterials = iNumMatrices = rgMeshPos.size();

                for (index_type i = 0; i < iNumMatrices; ++i) {
                    BearingGeometry::Type eMeshPos = BearingGeometry::CYLINDRICAL_MESH_UNKNOWN;
                    index_type iMeshIndex = -1;

                    if (HP.IsKeyWord("matrix" "at" "shaft")) {
                        eMeshPos = BearingGeometry::CYLINDRICAL_MESH_AT_SHAFT;
                        iMeshIndex = 1;
                    } else if (HP.IsKeyWord("matrix" "at" "bearing")) {
                        eMeshPos = BearingGeometry::CYLINDRICAL_MESH_AT_BEARING;
                        iMeshIndex = 2;
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): keywords \"matrix at shaft\" or \"matrix at bearing\" expected at "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    const index_type iIndex = eMeshPos == pGeometry->GetType()
                        ? ComplianceModelNodalDouble::DEHD_BODY_FIXED
                        : ComplianceModelNodalDouble::DEHD_BODY_MOVING;

                    if (rgMeshPos[iIndex] != BearingGeometry::CYLINDRICAL_MESH_UNKNOWN) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): the same mesh position has been used more than once at "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    rgMeshPos[iIndex] = eMeshPos;

                    if (HP.IsKeyWord("elastic" "half" "space")) {
                        if (!pContact.get()) {
                            rgMaterials[iIndex].ParseInput(iMeshIndex, HP, pGetParent());
                        }

                        rgMatrices[iIndex].reset(new ElasticHalfSpace(rgMaterials[iIndex], Material::Rigid()));
                    } else if (HP.IsKeyWord("from" "file")) {
                        const std::string strFileName = ParseFileName(HP);

                        if (!pContact.get()) {
                            // Needed for dPressScale only
                            rgMaterials[iIndex].ParseInput(iMeshIndex, HP, pGetParent());
                        }

                        rgMatrices[iIndex].reset(new ComplianceFromFile(strFileName));
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): keywords \"elastic half space\" or \"from file\" expected at "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    if (HP.IsKeyWord("modal" "element")) {
                        rgModalJoints[iIndex] = pDM->ReadElem<const Modal, const Joint, Elem::JOINT>(HP);
                    }
                }

                if (HP.IsKeyWord("axial" "displacement")) {
                    if (HP.IsKeyWord("small")) {
                        eInterpolOption = ComplianceModelNodalDouble::INT_AXIAL_SMALL_DISP;
                    } else if (HP.IsKeyWord("large")) {
                        eInterpolOption = ComplianceModelNodalDouble::INT_AXIAL_LARGE_DISP;                    
                    } else if (HP.IsKeyWord("extrapolate")) {
                        eInterpolOption = ComplianceModelNodalDouble::INT_AXIAL_EXTRAPOLATE;
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): keywords \"small\" or \"large\" expected at "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }
            } else {
                eCompModType = ComplianceModel::COMP_MOD_NODAL;

                for (index_type i = 0; i < iNumMatrices; ++i) {
                    if (!HP.IsKeyWord("matrix")) {
                        iNumMaterials = iNumMatrices = i;
                        break;
                    }

                    const integer iIndex = HP.GetInt();

                    if (iIndex <= 0 || iIndex > iNumMatrices) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): invalid index at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    if (HP.IsKeyWord("elastic" "half" "space")) {
                        if (!pContact.get()) {
                            rgMaterials[i].ParseInput(iIndex, HP, pGetParent());
                        }

                        rgMatrices[i].reset(new ElasticHalfSpace(rgMaterials[i], Material::Rigid()));
                    } else if (HP.IsKeyWord("from" "file")) {
                        const std::string strFileName = ParseFileName(HP);

                        if (!pContact.get()) {
                            // Needed for dPressScale only
                            rgMaterials[i].ParseInput(iIndex, HP, pGetParent());
                        }

                        rgMatrices[i].reset(new ComplianceFromFile(strFileName));
                    } else {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): keywords \"elastic half space\" or \"from file\" expected at "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }

                if (HP.IsKeyWord("modal" "element")) {
                    rgModalJoints[0] = pDM->ReadElem<const Modal, const Joint, Elem::JOINT>(HP);
                }
            }

            if (iNumMatrices < 1) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): at least one compliance matrix required at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            doublereal dPressScale = 0.;

            switch (iNumMaterials) {
            case 1:
                dPressScale = 1. / rgMaterials[0].dGetReducedModulus();
                break;
            case 2:
                dPressScale = 1. / rgMaterials[0].dGetReducedModulus(rgMaterials[1]);
                break;
            default:
                HYDRO_ASSERT(0);
            };

            doublereal dDefScale = pGetGeometry()->dGetReferenceClearance();

            if (HP.IsKeyWord("deformation" "dof" "scale")) {
                dDefScale = HP.GetReal();
            }

            switch (eCompModType) {
            case ComplianceModel::COMP_MOD_MODAL:
                pCompliance.reset(new ComplianceModelModal(this,
                                                           dDefScale,
                                                           dPressScale,
                                                           strFileNameModal));
                break;
            case ComplianceModel::COMP_MOD_NODAL: {
                pCompliance.reset(new ComplianceModelNodal(this,
                                                           rgModalJoints[0],
                                                           dDefScale,
                                                           dPressScale,
                                                           std::move(rgMatrices)));
            } break;
            case ComplianceModel::COMP_MOD_NODAL_DOUBLE:
                pCompliance.reset(new ComplianceModelNodalDouble(this,
                                                                 rgModalJoints,
                                                                 dDefScale,
                                                                 dPressScale,
                                                                 std::move(rgMatrices),
                                                                 *pGeometry,
                                                                 eInterpolOption));
                break;

            default:
                HYDRO_ASSERT(false);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }
    }

    void HydroMesh::GenerateBoundaryConditions()
    {
        HydroRootElement* const pParent = pGetParent();

        for (auto i = rgBoundaryCond.begin(); i != rgBoundaryCond.end(); ++i) {
            pParent->AddBoundaryCondition(std::move(*i));
        }

        rgBoundaryCond.clear();
    }

    void HydroMesh::GenerateMovingLubricationGrooves()
    {
        size_t iNumMovingLubrGrooves = 0;

        for (auto i = rgGrooves.begin(); i != rgGrooves.end(); ++i) {
            if ((*i)->GetType() == LubricationGroove::MOVING) {
                ++iNumMovingLubrGrooves;
            }
        }

        pGeometry->ReserveMovingLubrGrooves(iNumMovingLubrGrooves);

        for (auto i = rgGrooves.begin(); i != rgGrooves.end(); ++i) {
            if ((*i)->GetType() == LubricationGroove::MOVING) {
                pGeometry->AddMovingLubrGroove(i->get());
            }
        }
    }

    void HydroMesh::GenerateComplianceModel()
    {
        if (pCompliance) {
            pCompliance->EnableInitAss(pGetParent()->bInitialAssembly(HydroRootElement::INIT_ASS_ELASTICITY));

            pGetParent()->AddElement(std::unique_ptr<ComplianceModel>(pCompliance.release()));
        }
    }

    void HydroMesh::GetPressure(const HydroNode* pNode, doublereal& p, LocalDofMap*, doublereal) const
    {
        pNode->GetPressure(p);
        pGetParent()->pGetFluid()->Cavitation(p);
    }

    void HydroMesh::GetPressure(const HydroNode* pNode, Gradient<0>& p, LocalDofMap* pDofMap, doublereal dCoef) const
    {
        pNode->GetPressure(p, pDofMap, dCoef);
        pGetParent()->pGetFluid()->Cavitation(p);
    }

    integer HydroMesh::iGetNumBounaryConditions() const
    {
        return rgBoundaryCond.size();
    }

    CylindricalBearing* HydroMesh::pGetGeometry() const
    {
        return pGeometry.get();
    }

    ComplianceModel* HydroMesh::pGetComplianceModel() const
    {
        return pCompliance;
    }

    void
    HydroMesh::Update(const VectorHandler& XCurr,
                      const VectorHandler& XPrimeCurr, doublereal dCoef,
                      enum FunctionCall func,
                      LocalDofMap* pDofMap)
    {
        pGetGeometry()->Update(dCoef, func, pDofMap);

        if (pCompliance) {
            pCompliance->Update(XCurr, XPrimeCurr, dCoef, func);
        }
    }

    doublereal HydroMesh::dGetMaxPressureGradient() const
    {
        return pGetParent()->dGetMaxPressureGradient();
    }

    std::ostream& HydroMesh::PrintLogFile(std::ostream& os) const
    {
        const BearingGeometry* pGeometry = pGetGeometry();

        os << pGeometry->GetType() << ' ';

        pGeometry->PrintLogFile(os);

        return os;
    }

    LubricationGroove* HydroMesh::pFindGroove(const Vector<doublereal, 2>& x, Node2D::NodeType eNodeType, integer iNodeId)
    {
        LubricationGroove* pGroove = nullptr;

        for (auto it = rgGrooves.cbegin(); it != rgGrooves.cend(); ++it) {
            if ((*it)->GetType() == LubricationGroove::FIXED
                && (*it)->pGetGeometry()->bPointIsInside(x)
                && (*it)->pGetBoundaryCond()->bIncludeNode(eNodeType)) {

                if (pGroove == nullptr) {
                    pGroove = it->get();
                } else {
                    pedantic_cerr("hydrodynamic plain bearing2("
                                  << pGetParent()->GetLabel()
                                  << "): warning: node " << iNodeId + 1
                                  << " could be part of more than one groove" << std::endl);
                }
            }
        }

        return pGroove;
    }

    PressureCouplingCond* HydroMesh::pFindCouplingCond(const Vector<doublereal, 2>& x, integer iNodeId)
    {
        PressureCouplingCond* pCouplingCond = 0;

        for (auto it = rgCouplingCond.cbegin(); it != rgCouplingCond.cend(); ++it) {
            if ((*it)->pGetGeometry()->bPointIsInside(x)) {
                if (pCouplingCond == nullptr) {
                    pCouplingCond = it->get();
                } else {
                    pedantic_cerr("hydrodynamic plain bearing2("
                                  << pGetParent()->GetLabel()
                                  << "): warning node " << iNodeId + 1
                                  << " could be part of more than one pressure coupling conditions"
                                  << std::endl);
                }
            }
        }

        return pCouplingCond;
    }

    LinFDMesh::LinFDMesh(HydroRootElement* pParent)
        :HydroMesh(pParent),
         M(0),
         N(0),
         eElemType(CENT_DIFF_5)
    {
    }

    LinFDMesh::~LinFDMesh()
    {
    }

    void LinFDMesh::ParseInput(DataManager* pDM, MBDynParser& HP)
    {
        if (HP.IsKeyWord("element" "type")) {
            if (HP.IsKeyWord("central" "difference" "5")) {
                eElemType = CENT_DIFF_5;
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keywords \"central difference 5\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        ParseGeometry(pDM, HP);

        enum GridSpacing {
            UNIFORM,
            UNIFORM_INTERVALS
        } gridSpacing = UNIFORM;

        if (HP.IsKeyWord("grid" "spacing")) {
            if (HP.IsKeyWord("uniform")) {
                gridSpacing = UNIFORM;
            } else if (HP.IsKeyWord("uniform" "intervals")) {
                gridSpacing = UNIFORM_INTERVALS;
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keywords \"uniform\" or "
                            "\"uniform intervals\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        const doublereal b = pGeometry->dGetBearingWidth();
        const doublereal r = pGeometry->dGetMeshRadius();

        const integer iNumNodesZMin = 3;

        switch (gridSpacing) {
        case UNIFORM:
        {
            integer iNumNodesZ = -1;

            if (HP.IsKeyWord("mesh" "size" "axial")) {
                iNumNodesZ = floor(b / HP.GetReal() + 1.5);
            } else if (HP.IsKeyWord("number" "of" "nodes" "z")) {
                iNumNodesZ = HP.GetInt();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of nodes z\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!(iNumNodesZ >= iNumNodesZMin)) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): number of nodes in z-direction must be at least "
                            << iNumNodesZMin << " at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            M = iNumNodesZ - 1;

            integer iNumNodesPhi = -1;

            if (HP.IsKeyWord("mesh" "size" "circumferential")) {
                iNumNodesPhi = floor(2 * r * M_PI / HP.GetReal() + 2.5);
            } else if (HP.IsKeyWord("number" "of" "nodes" "Phi")) {
                iNumNodesPhi = HP.GetInt();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of nodes Phi\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!(iNumNodesPhi >= 4)) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): number of nodes in Phi-direction must be at least four at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            N = iNumNodesPhi - 1;

            const doublereal Deltax = 2 * M_PI * r / (N - 1);

            x.resize(N + 1);

            for (int j = 0; j <= N; ++j) {
                x[j] = (j - 1) * Deltax;
            }

            z.resize(M + 1);

            const doublereal Deltaz = b / M;

            for (int i = 0; i <= M; ++i) {
                z[i] = i * Deltaz - 0.5 * b;
            }
        }
        break;

        case UNIFORM_INTERVALS:
        {
            if (!HP.IsKeyWord("number" "of" "intervals" "z")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of intervals z\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const integer iNumIntZ = HP.GetInt();

            if (iNumIntZ < 1) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): at least one interval expected in z-direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            std::vector<doublereal> zi(iNumIntZ + 1);
            std::vector<integer> zni(iNumIntZ);
            M = 0;

            for (int i = 0; i <= iNumIntZ; ++i) {
                HYDRO_ASSERT(size_t(i) < zi.size());

                zi[i] = HP.GetReal();

                if (i > 0) {
                    if (zi[i] <= zi[i - 1]) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): z-coordinates must be in ascending order at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    if (zi[i] - zi[i - 1] > b) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): z-coordinate out of range b at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }

                if (i < iNumIntZ) {
                    HYDRO_ASSERT(size_t(i) < zni.size());

                    M += zni[i] = HP.GetInt();

                    if (zni[i] < 1) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): at least one subdivision expected at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }
            }

            if (M < iNumNodesZMin) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): at least " << iNumNodesZMin
                            << " nodes are required "
                            "in z-direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!HP.IsKeyWord("number" "of" "intervals" "x")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of intervals x\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const integer iNumIntX = HP.GetInt();

            if (iNumIntX < 1) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): at least one interval expected in x-direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            std::vector<doublereal> xi(iNumIntX + 1);
            std::vector<integer> xni(iNumIntX + 1);

            N = 1;

            for (int i = 0; i <= iNumIntX - 1; ++i) {
                HYDRO_ASSERT(size_t(i) < xi.size());

                xi[i] = HP.GetReal();

                if (xi[i] - xi[0] >= 2 * M_PI * r) {
                    silent_cerr("hydrodynamic plain bearing2(" <<
                                pGetParent()->GetLabel()
                                << "): x-coordinate exceeds range 0:2*pi*r at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }

                if (i > 0) {
                    if (xi[i] <= xi[i - 1]) {
                        silent_cerr("hydrodynamic plain bearing2("
                                    << pGetParent()->GetLabel()
                                    << "): x-coordinates must be in ascending order at line "
                                    << HP.GetLineData() << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }
                }

                HYDRO_ASSERT(size_t(i) < xni.size());

                N += xni[i] = HP.GetInt();

                if (xni[i] < 1) {
                    silent_cerr("hydrodynamic plain bearing2("
                                << pGetParent()->GetLabel()
                                << "): at least one subdivision expected at line "
                                << HP.GetLineData() << std::endl);
                    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                }
            }

            if (N < 3) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): at least three nodes are required in x-direction at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            x.resize(N + 1);
            z.resize(M + 1);

            int k = 0;
            z[k] = zi[k];

            for (int i = 1; i <= iNumIntZ; ++i) {
                HYDRO_ASSERT(size_t(i - 1) < zni.size());
                HYDRO_ASSERT(size_t(i) < zi.size());

                const int n = zni[i - 1];
                const doublereal Deltaz = (zi[i] - zi[i - 1]) / n;

                for (int j = 1; j <= n; ++j) {
                    HYDRO_ASSERT(size_t(k + 1) < z.size());

                    z[++k] = zi[i - 1] + j * Deltaz;
                }
            }

            k = 1;

            xi[iNumIntX] = xi[0] + 2 * M_PI * r;

            for (int i = 1; i <= iNumIntX; ++i) {
                HYDRO_ASSERT(size_t(i - 1) < xni.size());
                HYDRO_ASSERT(size_t(i) < xi.size());

                const int n = xni[i - 1];
                const doublereal Deltax = (xi[i] - xi[i - 1]) / n;

                for (int j = 1; j <= n; ++j) {
                    HYDRO_ASSERT(size_t(k + 1) < x.size());

                    x[k + 1] = x[k] + Deltax;
                    ++k;
                }
            }

            x[N] = 2 * x[N - 1] - x[N - 2];
            x[0] = x[1] - (x[N] - x[N - 1]);
        }
        break;

        default:
            HYDRO_ASSERT(0);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        HYDRO_ASSERT(x.size() == size_t(N + 1));
        HYDRO_ASSERT(z.size() == size_t(M + 1));

#if HYDRO_DEBUG > 1
        for (int i = 0; i <= N; ++i) {
            HYDRO_TRACE("x[" << std::setw(3) << i << "]=" << std::setw(10) << x[i]);

            if (i > 0) {
                HYDRO_TRACE(" x[" << std::setw(3) << i << "] - x[" << std::setw(3) << i - 1 << "]=" << std::setw(10) << x[i] - x[i - 1]);
            }
            HYDRO_TRACE(std::endl);
        }

        HYDRO_TRACE(std::endl);

        for (int i = 0; i <= M; ++i) {
            HYDRO_TRACE("z[" << std::setw(3) << i << "]=" << std::setw(10) << z[i]);

            if (i > 0) {
                HYDRO_TRACE(" z[" << std::setw(3) << i << "] - z[" << std::setw(3) << i - 1 << "]=" << std::setw(10) << z[i] - z[i - 1]);
            }

            HYDRO_TRACE(std::endl);
        }

        HYDRO_TRACE(std::endl);

        for (int i = 1; size_t(i) < x.size(); ++i) {
            HYDRO_ASSERT(x[i] > x[i - 1]);
            if (i > 1) {
                HYDRO_ASSERT(x[i] - x[1] <= 2 * M_PI * r * (1. + sqrt(std::numeric_limits<doublereal>::epsilon())));
            }
        }

        HYDRO_ASSERT(std::abs((x[N] - x[N - 1]) / (x[1] - x[0]) - 1.) < sqrt(std::numeric_limits<doublereal>::epsilon()));

        for (int i = 1; size_t(i) < z.size(); ++i) {
            HYDRO_ASSERT(z[i] > z[i - 1]);
        }
#endif
        ParseBoundaryCond(pDM, HP);
        ParsePressureCouplingCond(pDM, HP);
        ParseLubricationGrooves(pDM, HP);
        ParseThermWallBoundCond(pDM, HP);
        ParseContactModel(pDM, HP);
        ParseComplianceModel(pDM, HP);
    }

    integer LinFDMesh::iGetNumNodes() const
    {
        HYDRO_ASSERT(M > 0);
        HYDRO_ASSERT(N > 0);

        return (M + 1) * (N + 1) * (bThermalModel + 1)
            + (M - 1) * N // FluxNode qx
            + M * N; // FluxNode qz
    }

    integer LinFDMesh::iGetNumElements() const
    {
        HYDRO_ASSERT(M > 1);
        HYDRO_ASSERT(N > 1);

        return (M - 1) * (N - 1) * (bThermalModel + 1)
            + M * (N - 1)
            + 2 * (N - 1) * bUseOutletAxial
            + (pCompliance != nullptr);
    }

    void LinFDMesh::Generate()
    {
        HYDRO_ASSERT(M >= 2);
        HYDRO_ASSERT(N >= 3);

        HydroRootElement* const pParent = pGetParent();

        const HydroFluid* const pFluid = pParent->pGetFluid();

        const bool bIncompressible = pFluid->GetHydraulicType() == HydroFluid::INCOMPRESSIBLE;
        const bool bInitAssThermal = pGetParent()->bInitialAssembly(HydroRootElement::INIT_ASS_THERMAL);

        GenerateBoundaryConditions();
        GenerateMovingLubricationGrooves();

        // active pressure nodes
        for (integer i = 1; i <= M - 1; ++i) {
            for (integer j = 1; j <= N - 1; ++j) {
                const Vector<doublereal, 2> x = GetNodePosition(i, j);
                const integer iNodeIndex = iGetNodeIndexHydro(i, j);
                LubricationGroove* pGroove = pFindGroove(x, Node2D::HYDRAULIC_NODE, iNodeIndex);
                PressureCouplingCond* pCoupling = pFindCouplingCond(x, iNodeIndex);
                std::unique_ptr<HydroNode> pNode;
                std::unique_ptr<FrictionModel> pFrictionNode;

                if (pFriction) {
                    pFrictionNode = pFriction->Clone();
                }

                HYDRO_ASSERT(pGroove == nullptr || pGroove->pGetBoundaryCond()->bIncludeNode(Node2D::HYDRAULIC_NODE));

                if (pGroove == nullptr && pCoupling == nullptr) {
                    if (bIncompressible) {
                        pNode.reset(new HydroActiveNode(iNodeIndex,
                                                        x,
                                                        this,
                                                        pContact.get(),
                                                        std::move(pFrictionNode)));
                    } else {
                        pNode.reset(new HydroActiveComprNode(iNodeIndex,
                                                             x,
                                                             this,
                                                             pContact.get(),
                                                             std::move(pFrictionNode)));
                    }
                } else if (pCoupling != nullptr) {
                    if (pGroove != nullptr) {
                        pedantic_cerr("hydrodynamic plain bearing2("
                                      << pGetParent()->GetLabel()
                                      << "): warning node("
                                      << iNodeIndex
                                      << ") could be subject of lubrication groove ("
                                      << pGroove->iGetLabel()
                                      << ") or pressure coupling condition ("
                                      << pCoupling->iGetLabel()
                                      << ") but only pressure coupling will be effective!"
                                      << std::endl);
                    }

                    if (bIncompressible) {
                        pNode.reset(new HydroCoupledNode(iNodeIndex,
                                                         x,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         pCoupling->pGetNode()));
                    } else {
                        pNode.reset(new HydroCoupledComprNode(iNodeIndex,
                                                              x,
                                                              this,
                                                              pContact.get(),
                                                              std::move(pFrictionNode),
                                                              pCoupling->pGetNode()));
                    }

                    pCoupling->AddNode(pNode.get());
                } else {
                    HYDRO_ASSERT(pGroove != nullptr);
                    HYDRO_ASSERT(pCoupling == nullptr);

                    if (bIncompressible) {
                        pNode.reset(new HydroPassiveNode(iNodeIndex,
                                                         x,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         pGroove->pGetBoundaryCond()));
                    } else {
                        pNode.reset(new HydroPassiveComprNode(iNodeIndex,
                                                              x,
                                                              this,
                                                              pContact.get(),
                                                              std::move(pFrictionNode),
                                                              pGroove->pGetBoundaryCond()));
                    }

                    pGroove->AddNode(pNode.get());
                }

                HYDRO_ASSERT(pNode != nullptr);

                pParent->AddNode(std::move(pNode));
            }
        }

        // axial pressure boundary condition
        for (integer i = 0; i <= M; i += M) {
            for (integer j = 1; j <= N - 1; ++j) {
                const integer iNodeIndex = iGetNodeIndexHydro(i, j);
                const Vector<doublereal, 2> x = GetNodePosition(i, j);
                std::unique_ptr<HydroNode> pNode;
                std::unique_ptr<FrictionModel> pFrictionNode;

                if (pFriction) {
                    pFrictionNode = pFriction->Clone();
                }

                if (bUseOutletAxial) {
                    if (bIncompressible) {
                        pNode.reset(new HydroCoupledNode(iNodeIndex,
                                                         x,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         rgOutletAxial[i / M].pExtHydroNode));
                    } else {
                        pNode.reset(new HydroCoupledComprNode(iNodeIndex,
                                                              x,
                                                              this,
                                                              pContact.get(),
                                                              std::move(pFrictionNode),
                                                              rgOutletAxial[i / M].pExtHydroNode));
                    }
                } else {
                    LubricationGroove* pGroove = pFindGroove(x, Node2D::HYDRAULIC_NODE, iNodeIndex);
                    FluidStateBoundaryCond* pBoundaryCond = nullptr;

                    if (pGroove == nullptr) {
                        pBoundaryCond = pParent->pGetBoundaryCondition(i / M);
                    } else {
                        pBoundaryCond = pGroove->pGetBoundaryCond();
                    }

                    if (!pBoundaryCond->bIncludeNode(Node2D::HYDRAULIC_NODE)) {
                        silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                                    << "): pressure boundary condition not defined for node "
                                    << iNodeIndex + 1 << std::endl);
                        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                    }

                    if (bIncompressible) {
                        pNode.reset(new HydroPassiveNode(iNodeIndex,
                                                         x,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         pBoundaryCond));
                    } else if (pBoundaryCond->GetExtrapMethod() == FluidStateBoundaryCond::EX_OUTLET) {
                        pNode.reset(new HydroComprOutletNode(iNodeIndex,
                                                             x,
                                                             this,
                                                             pContact.get(),
                                                             std::move(pFrictionNode),
                                                             pBoundaryCond,
                                                             pParent->pGetNode<HydroMasterNode>(iGetNodeIndexHydro(i + 1 - 2 * i / M, j))));
                    } else {
                        pNode.reset(new HydroPassiveComprNode(iNodeIndex,
                                                              x,
                                                              this,
                                                              pContact.get(),
                                                              std::move(pFrictionNode),
                                                              pBoundaryCond));
                    }

                    if (pGroove != nullptr) {
                        pGroove->AddNode(pNode.get());
                    }
                }

                pParent->AddNode(std::move(pNode));
            }
        }

        // circumferential periodic pressure boundary condition
        for (integer i = 0; i <= M; ++i) {
            for (integer j = 0; j <= N; j += N) {
                const integer iMasterNode = iGetNodeIndexHydro(i, j == 0 ? N - 1 : 1);
                HYDRO_ASSERT(pParent->pGetNode(iMasterNode) != 0);

                pParent->AddNode(std::unique_ptr<HydroNode>{new HydroSlaveNode(iGetNodeIndexHydro(i, j),
                                                                               GetNodePosition(i, j),
                                                                               this,
                                                                               pParent->pGetNode<HydroNode>(iMasterNode))});
            }
        }

        std::array<const HydroNode*, 2> rgNodesFlux;

        FluxNode::NodeDataReq eFluxData = FluxNode::ND_NONE;

        if (bThermalModel && (pParent->uGetOutputFlags() & HydroRootElement::OUTPUT_HEAT_FLUX)) {
            eFluxData = FluxNode::ND_THERMAL;
        } else if (pParent->uGetOutputFlags() & (HydroRootElement::OUTPUT_MASS_FLUX | HydroRootElement::OUTPUT_VOLUME_FLUX)) {
            eFluxData = FluxNode::ND_HYDRAULIC;
        }

        // Flux in x-direction
        for (integer i = 0; i < M - 1; ++i) {
            for (integer j = 0; j < N; ++j) {
                const integer iNodeIndex = iGetNodeIndexFluxX(i, j);

                for (integer k = 0; k < 2; ++k) {
                    rgNodesFlux[k] = pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + 1, j + k));
                }

                std::unique_ptr<FluxNode> pNode{new FluxNode(iNodeIndex,
                                                             this,
                                                             rgNodesFlux,
                                                             FluxNode::PRESSURE_FROM_NODE,
                                                             eFluxData)};

                pParent->AddNode(std::move(pNode));
            }
        }

        // Flux in z-direction
        for (integer i = 0; i < M; ++i) {
            for (integer j = 0; j < N; ++j) {
                const integer iNodeIndex = iGetNodeIndexFluxZ(i, j);

                for (integer k = 0; k < 2; ++k) {
                    rgNodesFlux[k] = pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + k, j + 1));
                }

                std::unique_ptr<FluxNode> pNode{new FluxNode(iNodeIndex,
                                                             this,
                                                             rgNodesFlux,
                                                             FluxNode::PRESSURE_FROM_NODE,
                                                             eFluxData)};

                pParent->AddNode(std::move(pNode));
            }
        }

        if (bThermalModel) {
            const bool bInitAssThermal = pGetParent()->bInitialAssembly(HydroRootElement::INIT_ASS_THERMAL);

            const doublereal T0 = pGetParent()->pGetFluid()->dGetRefTemperature();

            for (integer i = 1; i <= M - 1; ++i) {
                for (integer j = 1; j <= N - 1; ++j) {
                    const Vector<doublereal, 2> x = GetNodePosition(i, j);
                    const integer iNodeIndex = iGetNodeIndexTherm(i, j);
                    LubricationGroove* pGroove = pFindGroove(x, Node2D::THERMAL_NODE, iNodeIndex);
                    PressureCouplingCond* pCoupling = pFindCouplingCond(x, iNodeIndex);
                    std::unique_ptr<ThermoHydrNode> pNode;

                    if (pGroove == nullptr && pCoupling == nullptr) {
                        pNode.reset(new ThermalActiveNode(iNodeIndex,
                                                          x,
                                                          this,
                                                          T0,
                                                          bInitAssThermal));
                    } else if (pCoupling != nullptr) {
                        HYDRO_ASSERT(pCoupling->pGetThermalNode() != nullptr);

                        pNode.reset(new ThermalInletNode(iNodeIndex,
                                                         x,
                                                         this,
                                                         pCoupling->pGetThermalNode(),
                                                         bInitAssThermal));
                    } else {
                        pNode.reset(new ThermalPassiveNode(iNodeIndex,
                                                           x,
                                                           this,
                                                           pGroove->pGetBoundaryCond()));

                        pGroove->AddNode(pNode.get());
                    }

                    HYDRO_ASSERT(pNode != nullptr);

                    pParent->AddNode(std::move(pNode));
                }
            }

            // axial temperature boundary condition
            for (integer i = 0; i <= M; i += M) {
                for (integer j = 1; j <= N - 1; ++j) {
                    const integer iNodeIndex = iGetNodeIndexTherm(i, j);
                    const Vector<doublereal, 2> x = GetNodePosition(i, j);
                    std::unique_ptr<ThermoHydrNode> pNode;

                    if (bUseOutletAxial) {
                        pNode.reset(new ThermalCoupledNode(iNodeIndex,
                                                           x,
                                                           this,
                                                           rgOutletAxial[i / M].pExtThermNode));
                    } else {
                        LubricationGroove* pGroove = pFindGroove(x, Node2D::THERMAL_NODE, iNodeIndex);
                        FluidStateBoundaryCond* pBoundaryCond = nullptr;

                        if (pGroove == nullptr) {
                            pBoundaryCond = pParent->pGetBoundaryCondition(i / M);
                        } else {
                            pBoundaryCond = pGroove->pGetBoundaryCond();
                        }

                        if (!pBoundaryCond->bIncludeNode(Node2D::THERMAL_NODE)) {
                            silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                                        << "): temperature boundary condition not defined for node "
                                        << iNodeIndex + 1 << std::endl);
                            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
                        }

                        if (pBoundaryCond->GetExtrapMethod() == FluidStateBoundaryCond::EX_OUTLET) {
                            pNode.reset(new ThermalSlaveNode(iNodeIndex,
                                                             x,
                                                             pParent->pGetNode<ThermoHydrNode>(iGetNodeIndexTherm(i + 1 - 2 * i / M, j))));

                        } else {
                            pNode.reset(new ThermalPassiveNode(iNodeIndex,
                                                               x,
                                                               this,
                                                               pBoundaryCond));
                        }

                        if (pGroove != nullptr) {
                            pGroove->AddNode(pNode.get());
                        }
                    }

                    pParent->AddNode(std::move(pNode));
                }
            }

            // circumferential periodic temperature boundary condition
            for (integer i = 0; i <= M; ++i) {
                for (integer j = 0; j <= N; j += N) {
                    const integer iMasterNode = iGetNodeIndexTherm(i, j == 0 ? N - 1 : 1);

                    HYDRO_ASSERT(pParent->pGetNode(iMasterNode) != nullptr);

                    pParent->AddNode(std::unique_ptr<ThermoHydrNode>{new ThermalSlaveNode(iGetNodeIndexTherm(i, j),
                                                                                          GetNodePosition(i, j),
                                                                                          pParent->pGetNode<ThermoHydrNode>(iMasterNode))});
                }
            }

            for (integer i = 0; i <= M; ++i) {
                for (integer j = 0; j <= N; ++j) {
                    const integer iNodeHydro = iGetNodeIndexHydro(i, j);
                    const integer iNodeTherm = iGetNodeIndexTherm(i, j);
                    pParent->pGetNode<HydroNode>(iNodeHydro)->SetThermalNode(pParent->pGetNode<ThermoHydrNode>(iNodeTherm));
                }
            }
        }

        for (integer i = 1; i <= M - 1; ++i) {
            for (integer j = 1; j <= N - 1; ++j) {
                HydroNode* const pCenterNode = pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j));
                std::unique_ptr<LinFD5Elem> pElement;

                if (typeid(*pCenterNode) == typeid(HydroActiveNode)) {
                    pElement.reset(new LinFD5ReynoldsElem(this));
                } else if (typeid(*pCenterNode) == typeid(HydroActiveComprNode)) {
                    pElement.reset(new LinFD5ComprReynoldsElem(this));
                } else if (typeid(*pCenterNode) == typeid(HydroCoupledNode) ||
                           typeid(*pCenterNode) == typeid(HydroCoupledComprNode)) {
                    if (bThermalModel) {
                        HYDRO_ASSERT(typeid(*pCenterNode->pGetThermalNode()) == typeid(ThermalInletNode));
                        pElement.reset(new LinFD5ThermalCouplingElem(this, bInitAssThermal));
                    } else {
                        pElement.reset(new LinFD5CouplingElem(this));
                    }
                } else {
                    HYDRO_ASSERT(typeid(*pCenterNode) == typeid(HydroPassiveNode)
                                 || typeid(*pCenterNode) == typeid(HydroPassiveComprNode));
                    continue;
                }

                pElement->SetNode(LinFD5Elem::iNodeCenter, pCenterNode);
                pElement->SetNode(LinFD5Elem::iNodeWest,   pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j - 1)));
                pElement->SetNode(LinFD5Elem::iNodeEast,   pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j + 1)));
                pElement->SetNode(LinFD5Elem::iNodeSouth,  pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i - 1, j)));
                pElement->SetNode(LinFD5Elem::iNodeNorth,  pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + 1, j)));

                pElement->SetFluxNode(LinFD5Elem::iNodeFlxWest, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxX(i - 1, j - 1)));
                pElement->SetFluxNode(LinFD5Elem::iNodeFlxEast, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxX(i - 1, j)));
                pElement->SetFluxNode(LinFD5Elem::iNodeFlzSouth, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(i - 1, j - 1)));
                pElement->SetFluxNode(LinFD5Elem::iNodeFlzNorth, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(i, j - 1)));

                HYDRO_TRACE("Element(" << pParent->iGetNumElements() << "):" << typeid(*pElement).name() << std::endl);

                pParent->AddElement(std::move(pElement));
            }
        }

        for (integer i = 0; i <= M - 1; ++i) {
            for (integer j = 1; j <= N - 1; ++j) {
                std::unique_ptr<HydroElement> pElement{new LinFD4FrictionElem(this)};

                pElement->SetNode(LinFD4Elem::iNode1NE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + 1, j + 1)));
                pElement->SetNode(LinFD4Elem::iNode2NW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + 1, j)));
                pElement->SetNode(LinFD4Elem::iNode3SW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j)));
                pElement->SetNode(LinFD4Elem::iNode4SE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j + 1)));

                pParent->AddElement(std::move(pElement));
            }
        }

        if (bUseOutletAxial) {
            // pressure coupling condition axial
            for (integer j = 1; j <= N - 1; ++j) {
                std::unique_ptr<LinFD4MassFlowZ> pElement{new LinFD4MassFlowZ(this)}; // z = b/2

                pElement->SetNode(LinFD4Elem::iNode1NE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(M, j + 1)));
                pElement->SetNode(LinFD4Elem::iNode2NW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(M, j)));
                pElement->SetNode(LinFD4Elem::iNode3SW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(M - 1, j)));
                pElement->SetNode(LinFD4Elem::iNode4SE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(M - 1, j + 1)));

                pElement->SetFluxNode(LinFD4MassFlowZ::iFNodeWest, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(M - 1, j - 1)));
                pElement->SetFluxNode(LinFD4MassFlowZ::iFNodeEast, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(M - 1, j)));

                pParent->AddElement(std::move(pElement));

                pElement.reset(new LinFD4MassFlowZ(this));   // z = -b/2

                pElement->SetNode(LinFD4Elem::iNode3SW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(1, j + 1)));
                pElement->SetNode(LinFD4Elem::iNode4SE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(1, j)));
                pElement->SetNode(LinFD4Elem::iNode1NE, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(0, j)));
                pElement->SetNode(LinFD4Elem::iNode2NW, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(0, j + 1)));

                pElement->SetFluxNode(LinFD4MassFlowZ::iFNodeEast, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(0, j - 1)));
                pElement->SetFluxNode(LinFD4MassFlowZ::iFNodeWest, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(0, j)));

                pParent->AddElement(std::move(pElement));
            }
        }

        for (auto i = rgGrooves.cbegin(); i != rgGrooves.cend(); ++i) {
            if ((*i)->GetType() == LubricationGroove::FIXED && (*i)->iGetNumNodes() == 0) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): No nodes affected by lubrication groove ("
                            << (*i)->iGetLabel() << ")" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        for (auto i = rgCouplingCond.cbegin(); i != rgCouplingCond.cend(); ++i) {
            if ((*i)->iGetNumNodes() == 0) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): No nodes affected by pressure coupling condition ("
                            << (*i)->iGetLabel() << ")" << std::endl);

                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (bThermalModel) {
            for (integer i = 1; i <= M - 1; ++i) {
                for (integer j = 1; j <= N - 1; ++j) {
                    ThermalActiveNode* pCenterNode = dynamic_cast<ThermalActiveNode*>(pParent->pGetNode(iGetNodeIndexTherm(i, j)));

                    if (pCenterNode != nullptr) {
                        std::unique_ptr<LinFD5ThermalElem> pElement{new LinFD5ThermalElemImp(this, bInitAssThermal)};
                        pElement->SetNode(LinFD5Elem::iNodeCenter, pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j)));
                        pElement->SetNode(LinFD5Elem::iNodeWest,   pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j - 1)));
                        pElement->SetNode(LinFD5Elem::iNodeEast,   pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i, j + 1)));
                        pElement->SetNode(LinFD5Elem::iNodeSouth,  pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i - 1, j)));
                        pElement->SetNode(LinFD5Elem::iNodeNorth,  pParent->pGetNode<HydroNode>(iGetNodeIndexHydro(i + 1, j)));

                        pElement->SetFluxNode(LinFD5Elem::iNodeFlxWest, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxX(i - 1, j - 1)));
                        pElement->SetFluxNode(LinFD5Elem::iNodeFlxEast, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxX(i - 1, j)));
                        pElement->SetFluxNode(LinFD5Elem::iNodeFlzSouth, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(i - 1, j - 1)));
                        pElement->SetFluxNode(LinFD5Elem::iNodeFlzNorth, pParent->pGetNode<FluxNode>(iGetNodeIndexFluxZ(i, j - 1)));

                        pParent->AddElement(std::move(pElement));
                    }
                }
            }
        }

        GenerateComplianceModel();
    }

    std::ostream& LinFDMesh::Output(std::ostream& os) const
    {
        return os;
    }

    Vector<doublereal, 2> LinFDMesh::GetNodePosition(integer i, integer j) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(size_t(i) < z.size());
        HYDRO_ASSERT(j >= 0);
        HYDRO_ASSERT(size_t(j) < x.size());

        return Vector<doublereal, 2>(x[j], z[i]);
    }

    integer LinFDMesh::iGetNodeIndexHydro(integer i, integer j) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(size_t(i) < z.size());
        HYDRO_ASSERT(j >= 0);
        HYDRO_ASSERT(size_t(j) < x.size());

        // Attention: thermal nodes must be updated before hydraulic nodes!
        return (M + 1) * (N + 1) * bThermalModel + i * (N + 1) + j;
    }

    integer LinFDMesh::iGetNodeIndexFluxX(integer i, integer j) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(i < M - 1);
        HYDRO_ASSERT(j >= 0);
        HYDRO_ASSERT(j < N);

        return (M + 1) * (N + 1) * (bThermalModel + 1)
            + i * N
            + j;
    }

    integer LinFDMesh::iGetNodeIndexFluxZ(integer i, integer j) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(i < M);
        HYDRO_ASSERT(j >= 0);
        HYDRO_ASSERT(j < N);

        return (M + 1) * (N + 1) * (bThermalModel + 1)
            + (M - 1) * N
            + i * N + j;
    }

    integer LinFDMesh::iGetNodeIndexTherm(integer i, integer j) const
    {
        HYDRO_ASSERT(i >= 0);
        HYDRO_ASSERT(size_t(i) < z.size());
        HYDRO_ASSERT(j >= 0);
        HYDRO_ASSERT(size_t(j) < x.size());

        // Attention: thermal nodes must be updated before hydraulic nodes!
        return i * (N + 1) + j;
    }

    QuadFeIso9Mesh::QuadFeIso9Mesh(HydroRootElement* pParent)
        :HydroMesh(pParent),
         oIntegRuleReynolds(QuadFeIso9Elem::REYNOLDS_ELEM),
         oIntegRuleFriction(QuadFeIso9Elem::FRICTION_ELEM),
         dSkewMesh(0.),
         eMeshGeometry(PLAIN_BEARING)
    {
    }

    QuadFeIso9Mesh::~QuadFeIso9Mesh()
    {
    }

    void QuadFeIso9Mesh::ParseInput(DataManager* pDM, MBDynParser& HP)
    {
        HydroRootElement* const pParent = pGetParent();
        const HydroFluid* const pFluid = pParent->pGetFluid();

        if (pFluid->GetHydraulicType() != HydroFluid::INCOMPRESSIBLE ||
            pFluid->GetThermalType() != HydroFluid::ISOTHERMAL)
        {
            silent_cerr("hydrodynamic plain bearing2(" <<
                        pParent->GetLabel()
                        << "): only incompressible isothermal fluids are supported "
                        "for isoparameteric 9 elements at line "
                        << HP.GetLineData() << std::endl);
            throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
        }

        if (HP.IsKeyWord("element" "type")) {
            if (!HP.IsKeyWord("isoparametric" "9")) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keywords \"isoparametric 9\""
                            " expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        if (HP.IsKeyWord("number" "of" "gauss" "points" "reynolds")) {
            oIntegRuleReynolds.ParseInput(HP, this);
        }

        if (HP.IsKeyWord("number" "of" "gauss" "points" "friction")) {
            oIntegRuleFriction.ParseInput(HP, this);
        }

        ParseGeometry(pDM, HP);

        if (HP.IsKeyWord("mesh" "geometry")) {
            if (HP.IsKeyWord("plain" "bearing")) {
                eMeshGeometry = PLAIN_BEARING;
            } else if (HP.IsKeyWord("helical" "groove" "pump")) {
                eMeshGeometry = HELICAL_GROOVE_PUMP;
            } else {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): keywords \"plain bearing\" or \"helical groove pump\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }
        }

        switch (eMeshGeometry) {
        case PLAIN_BEARING: {
            index_type iNumElemZ = -1;

            if (HP.IsKeyWord("number" "of" "nodes" "z")) {
                iNumElemZ = HP.GetInt() / 2;
            } else if (HP.IsKeyWord("number" "of" "elements" "z")) {
                iNumElemZ = HP.GetInt();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of elements z\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!(iNumElemZ >= 1)) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword number of elements must be at least one at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iNumNodesZ = 2 * iNumElemZ + 1;

            z.Resize(iNumNodesZ);

            const doublereal w = pGeometry->dGetBearingWidth();

            for (index_type i = 1; i <= iNumNodesZ; ++i) {
                z(i) = (static_cast<doublereal>(i - 1) / (iNumNodesZ - 1) - 0.5) * w;
            }

            index_type iNumElemPhi = -1;

            if (HP.IsKeyWord("number" "of" "nodes" "Phi")) {
                iNumElemPhi = (HP.GetInt() - 1) / 2;
            } else if (HP.IsKeyWord("number" "of" "elements" "Phi")) {
                iNumElemPhi = HP.GetInt();
            } else {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): keyword \"number of elements Phi\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!(iNumElemPhi >= 1)) {
                silent_cerr("hydrodynamic plain bearing2("
                            << pGetParent()->GetLabel()
                            << "): number of elements Phi must be at least one at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const index_type iNumNodesPhi = 2 * iNumElemPhi + 1;

            x.Resize(iNumNodesPhi);

            const doublereal c = 2 * M_PI * pGeometry->dGetMeshRadius();

            for (index_type i = 1; i <= iNumNodesPhi; ++i) {
                x(i) = static_cast<double>(i - 1) / (iNumNodesPhi - 1) * c;
            }
        } break;
        case HELICAL_GROOVE_PUMP: {
            if (!HP.IsKeyWord("pitch" "angle")) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): keyword \"pitch angle\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal beta = HP.GetReal();

            if (beta <= 0. || beta >= M_PI) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << ") beta must be greater than zero and less than pi at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const integer K = HP.IsKeyWord("number" "of" "segments") ? HP.GetInt() : 1;

            if (K < 1) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): at least one segment is required at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!HP.IsKeyWord("bar" "width")) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): keyword \"bar width\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal Ws = HP.GetReal();

            if (Ws <= 0.) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): bar width must be greater than zero at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!HP.IsKeyWord("chamfer" "width")) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): keyword \"chamfer width\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal Wc = HP.GetReal();

            if (Wc <= 0.) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): chamfer width must be greater than zero at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            if (!HP.IsKeyWord("mesh" "size")) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): keyword \"mesh size\" expected at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal dx = HP.GetReal();

            if (dx <= 0.) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): mesh size must be greater than zero at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const doublereal B = pGeometry->dGetBearingWidth();
            const doublereal D = 2. * pGeometry->dGetMeshRadius();

            const doublereal dx1 = Ws / sin(beta);
            const doublereal dx2 = Wc / sin(beta);
            const doublereal dx3 = D * M_PI / K - dx1 - 2 * dx2;

            dSkewMesh = B / tan(beta);

            if (dx3 <= 0.) {
                silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                            << "): groove width is below zero at line "
                            << HP.GetLineData() << std::endl);
                throw ErrGeneric(MBDYN_EXCEPT_ARGS);
            }

            const integer Nx1 = 2 * ceil(dx1 / dx) + 1;
            const integer Nx2 = 2 * ceil(dx2 / dx) + 1;
            const integer Nx3 = 2 * ceil(dx3 / dx) + 1;
            const integer Nx = K * (Nx1 + 2 * Nx2 + Nx3 - 4) + 1;
            const integer Nz = 2 * ceil(B / dx) + 1;
            const doublereal dz = B / (Nz - 1);

            x.Resize(Nx);

            index_type ix = 1;

            HYDRO_ASSERT(x(ix) == 0.);

            for (index_type i = 0; i < K; ++i) {
                for (index_type j = 1; j < Nx2; ++j) {
                    x(ix + 1) = x(ix) + dx2 / (Nx2 - 1);
                    ++ix;
                }

                for (index_type j = 1; j < Nx1; ++j) {
                    x(ix + 1) = x(ix) + dx1 / (Nx1 - 1);
                    ++ix;
                }

                for (index_type j = 1; j < Nx2; ++j) {
                    x(ix + 1) = x(ix) + dx2 / (Nx2 - 1);
                    ++ix;
                }

                for (index_type j = 1; j < Nx3; ++j) {
                    x(ix + 1) = x(ix) + dx3 / (Nx3 - 1);
                    ++ix;
                }
            }

            HYDRO_ASSERT(ix == x.iGetNumRows());

            x(ix) = D * M_PI;

            HYDRO_ASSERT(fabs((x(ix) - x(ix - 1)) / (dx3 / (Nx3 - 1)) - 1.) < sqrt(std::numeric_limits<doublereal>::epsilon()));

            z.Resize(Nz);

            for (index_type i = 1; i <= Nz; ++i) {
                z(i) = (i - 1) * dz - 0.5 * B;
            }

            HYDRO_ASSERT(fabs((z(Nz) - z(1)) / B - 1) < sqrt(std::numeric_limits<doublereal>::epsilon()));
        } break;
        default:
            HYDRO_ASSERT(false);
        }

        ParseBoundaryCond(pDM, HP);

        ParsePressureCouplingCond(pDM, HP);

        if (rgCouplingCond.size()) {
            silent_cerr("hydrodynamic plain bearing2(" << pGetParent()->GetLabel()
                        << "): pressure coupling conditions radial not implemented "
                        "for quadratic finite element mesh at line "
                        << HP.GetLineData() << std::endl);
            throw ErrNotImplementedYet(MBDYN_EXCEPT_ARGS);
        }

        ParseLubricationGrooves(pDM, HP);
        ParseContactModel(pDM, HP);
        ParseComplianceModel(pDM, HP);
    }

    integer QuadFeIso9Mesh::iGetNumNodes() const
    {
        return x.iGetNumRows() * z.iGetNumRows();
    }

    integer QuadFeIso9Mesh::iGetNumElements() const
    {
        index_type iNumElem = 2 * (x.iGetNumRows() - 1) * (z.iGetNumRows() - 1) / 4;

        if (bUseOutletAxial) {
            iNumElem += x.iGetNumRows() - 1;
        }

        if (pCompliance) {
            ++iNumElem;
        }

        return iNumElem;
    }

    void QuadFeIso9Mesh::Generate()
    {
        HYDRO_ASSERT(x.iGetNumRows() >= 3);
        HYDRO_ASSERT(z.iGetNumRows() >= 3);
        HYDRO_ASSERT(x.iGetNumRows() & 1);
        HYDRO_ASSERT(z.iGetNumRows() & 1);

        GenerateBoundaryConditions();
        GenerateMovingLubricationGrooves();

        std::unique_ptr<HydroNode> pNode;

        for (index_type i = 1; i <= x.iGetNumRows(); ++i) {
            for (index_type j = 1; j <= z.iGetNumRows(); ++j) {
                const Vector<doublereal, 2> xn = GetNodePosition(i, j);
                const index_type iNodeIndex = iGetNodeIndex(i, j);
                LubricationGroove* const pGroove = pFindGroove(xn, Node2D::HYDRAULIC_NODE, iNodeIndex);
                const FluidStateBoundaryCond* pBoundCond = pGroove ? pGroove->pGetBoundaryCond() : pFindBoundaryCond(i, j);
                std::unique_ptr<FrictionModel> pFrictionNode;

                if (pFriction) {
                    pFrictionNode = pFriction->Clone();
                }

                HYDRO_ASSERT(pGroove == nullptr || pGroove->pGetBoundaryCond()->bIncludeNode(Node2D::HYDRAULIC_NODE));

                if (i < x.iGetNumRows()) {
                    if (bUseOutletAxial && (j == 1 || j == z.iGetNumRows())) {
                        pNode.reset(new HydroCoupledNode(iNodeIndex,
                                                         xn,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         rgOutletAxial[j / z.iGetNumRows()].pExtHydroNode));
                    } else if (!pBoundCond) {
                        pNode.reset(new HydroActiveNode(iNodeIndex,
                                                        xn,
                                                        this,
                                                        pContact.get(),
                                                        std::move(pFrictionNode)));

                    } else {
                        pNode.reset(new HydroPassiveNode(iNodeIndex,
                                                         xn,
                                                         this,
                                                         pContact.get(),
                                                         std::move(pFrictionNode),
                                                         pBoundCond));
                        if (pGroove) {
                            pGroove->AddNode(pNode.get());
                        }
                    }
                } else {
                    const index_type iNodeMaster = iGetNodeIndex(1, j);

                    pNode.reset(new HydroSlaveNode(iNodeIndex,
                                                   xn,
                                                   this,
                                                   pGetParent()->pGetNode<HydroNode>(iNodeMaster)));
                }

                pGetParent()->AddNode(std::move(pNode));
            }
        }

        static const struct NodeLayoutIso9 {
            index_type iNode;
            index_type iOffsetX;
            index_type iOffsetZ;
        } rgNodeLayout[QuadFeIso9Elem::iNumNodes] = {
            {0, 2, 2},
            {1, 0, 2},
            {2, 0, 0},
            {3, 2, 0},
            {4, 1, 2},
            {5, 0, 1},
            {6, 1, 0},
            {7, 2, 1},
            {8, 1, 1}
        };

        enum ElemType {
            REYNOLDS_ELEM,
            FRICTION_ELEM,
            COUPLING_ELEM
        };

        std::unique_ptr<QuadFeIso9Elem> pElem;

        for (index_type l = REYNOLDS_ELEM; l <= COUPLING_ELEM; ++l) {
            for (index_type i = 1; i <= x.iGetNumRows() - 2; i += 2) {
                for (index_type j = 1; j <= z.iGetNumRows() - 2; j += 2) {
                    switch (l) {
                    case REYNOLDS_ELEM:
                        pElem.reset(new QuadFeIso9ReynoldsElem(this, oIntegRuleReynolds));
                        break;
                    case FRICTION_ELEM:
                        pElem.reset(new QuadFeIso9FrictionElem(this, oIntegRuleFriction));
                        break;
                    case COUPLING_ELEM:
                        if (bUseOutletAxial && (j == 1 || j == z.iGetNumRows() - 2)) {
                            pElem.reset(new QuadFeIso9MassFlowZ(this, oIntegRuleReynolds, j == 1 ? -1. : 1.));
                        } else {
                            continue;
                        }
                        break;
                    default:
                        HYDRO_ASSERT(false);
                    }

                    bool bAddElem = (l == FRICTION_ELEM || l == COUPLING_ELEM);

                    for (auto k = std::begin(rgNodeLayout); k != std::end(rgNodeLayout); ++k) {
                        index_type iNodeIdx = iGetNodeIndex(i + k->iOffsetX, j + k->iOffsetZ);
                        HydroNode* pNode = pGetParent()->pGetNode<HydroNode>(iNodeIdx);
                        bAddElem = bAddElem || pNode->bIsNodeType(HydroNode::ACTIVE_NODE);
                        pElem->SetNode(k->iNode, pNode);
                    }

                    if (bAddElem) {
                        pGetParent()->AddElement(std::move(pElem));
                        HYDRO_ASSERT(pElem == nullptr);
                    }
                }
            }
        }

        GenerateComplianceModel();
    }

    std::ostream& QuadFeIso9Mesh::Output(std::ostream& os) const
    {
        return os;
    }

    index_type QuadFeIso9Mesh::iGetNodeIndex(index_type i, index_type j) const
    {
        HYDRO_ASSERT(i >= 1);
        HYDRO_ASSERT(i <= x.iGetNumRows());
        HYDRO_ASSERT(j >= 1);
        HYDRO_ASSERT(j <= z.iGetNumRows());
        const index_type iNode = (i - 1) * z.iGetNumRows() + j - 1;
        HYDRO_ASSERT(iNode >= 0);
        HYDRO_ASSERT(iNode < iGetNumNodes());
        return iNode;
    }

    Vector<doublereal, 2> QuadFeIso9Mesh::GetNodePosition(index_type i, index_type j) const
    {
        return Vector<doublereal, 2>(x(i) + dSkewMesh * z(j) / pGeometry->dGetBearingWidth(), z(j));
    }

    const FluidStateBoundaryCond* QuadFeIso9Mesh::pFindBoundaryCond(index_type i, index_type j) const
    {
        if (!bUseOutletAxial && (j == 1 || j == z.iGetNumRows())) {
            return pGetParent()->pGetBoundaryCondition(j > 1);
        }

        return nullptr;
    }


#if HYDRO_DEBUG > 0
    const HydroRootElement* pGetRootElement(const HydroElement* pElement)
    {
        return pElement->pGetMesh()->pGetParent();
    }

    const HydroRootElement* pGetRootElement(const BearingGeometry* pGeometry)
    {
        return pGeometry->pGetParent();
    }

    template <typename ElementType, index_type N_SIZE>
    bool bCheckNumColsWorkSpace(const ElementType* pElem, grad::FunctionCall eFunc, const Gradient<N_SIZE>& g, index_type iRowIndex)
    {
        index_type iNumCols = 0;

        for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
            if (g.dGetDerivativeLocal(i)) {
                ++iNumCols;
            }
        }

        integer iMaxRows, iMaxCols;

        pElem->WorkSpaceDim(&iMaxRows, &iMaxCols, eFunc);

        if (iNumCols <= iMaxCols) {
            return true;
        } else {
            silent_cerr("error: expected "
                        << iMaxCols
                        << " but got " << iNumCols << " nonzeros in row "
                        << iRowIndex << std::endl);


            std::cerr << std::setw(12) << g.dGetValue() << std::endl;
            const LocalDofMap* pDofMap = g.pGetDofMap();

            if (pDofMap) {
                std::vector<index_type> rgGlobalDof;
                rgGlobalDof.reserve(g.iGetLocalSize());
                for (index_type i = g.iGetStartIndexLocal(); i < g.iGetEndIndexLocal(); ++i) {
                    rgGlobalDof.push_back(g.iGetGlobalDof(i));
                }

                std::sort(rgGlobalDof.begin(), rgGlobalDof.end());

                for (size_t i = 0; i < rgGlobalDof.size(); ++i) {
                    const index_type iGlobal = rgGlobalDof[i];

                    std::cerr << iGlobal
                              << ": "
                              << pGetRootElement(pElem)->pGetDataManager()->GetDofDescription(iGlobal)
                              << " = "
                              << g.dGetDerivativeGlobal(iGlobal)
                              << std::endl;
                }

                std::cerr << std::endl;
            }

            return false;
        }
    }
#endif

    class GrooveShapeDriveCaller : public DriveCaller
    {
    public:
        GrooveShapeDriveCaller(const DriveHandler* pDH, doublereal Ws, doublereal Wc, doublereal Hg);

        virtual ~GrooveShapeDriveCaller();
        bool bIsDifferentiable(void) const;
        virtual std::ostream& Restart(std::ostream& out) const;
        doublereal dGet(const doublereal& dVar) const;
        virtual doublereal dGetP(const doublereal& dVar) const;
        virtual DriveCaller* pCopy(void) const;

    private:
        template <typename T>
        void dGet(const T& x, T& y) const;

        const doublereal Ws, Wc, Hg;
    };

    struct GrooveShapeDCR : public DriveCallerRead {
        DriveCaller *
        Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred);
    };

    GrooveShapeDriveCaller::GrooveShapeDriveCaller(const DriveHandler* pDH, doublereal Ws, doublereal Wc, doublereal Hg)
        : DriveCaller(pDH),
          Ws(Ws), Wc(Wc), Hg(Hg)
    {
        NO_OP;
    }

    GrooveShapeDriveCaller::~GrooveShapeDriveCaller()
    {
        NO_OP;
    }

    template <typename T>
    void
    GrooveShapeDriveCaller::dGet(const T& x, T& y) const
    {
        // y = -0.5 * Hg * (1. + tanh(M_PI * (2. * (fabs(x) - 0.5 * Ws) / Wc - 1.)));

        if (fabs(x) <= 0.5 * Ws) {
            y = 0.;
        } else if (fabs(x) <= 0.5 * Ws + Wc) {
            y = -Hg * (fabs(x) - 0.5 * Ws) / Wc;
        } else {
            y = -Hg;
        }
    }

    doublereal GrooveShapeDriveCaller::dGet(const doublereal& x) const
    {
        doublereal y;

        dGet(x, y);

        return y;
    }

    doublereal GrooveShapeDriveCaller::dGetP(const doublereal& dVar) const
    {
        Gradient<1> x, y;

        x.SetValuePreserve(dVar);
        x.DerivativeResizeReset(0, 0, MapVectorBase::LOCAL, 1.);

        dGet(x, y);

        return y.dGetDerivativeLocal(0);
    }

    bool GrooveShapeDriveCaller::bIsDifferentiable(void) const
    {
        return true;
    }

    std::ostream&
    GrooveShapeDriveCaller::Restart(std::ostream& out) const
    {
        return out;
    }

    DriveCaller*
    GrooveShapeDriveCaller::pCopy(void) const
    {
        DriveCaller* pDC = 0;

        SAFENEWWITHCONSTRUCTOR(pDC,
                               GrooveShapeDriveCaller,
                               GrooveShapeDriveCaller(pGetDrvHdl(), Ws, Wc, Hg));

        return pDC;
    }

    DriveCaller *
    GrooveShapeDCR::Read(const DataManager* pDM, MBDynParser& HP, bool bDeferred)
    {
        NeedDM(pDM, HP, bDeferred, "helical groove shape");

        const DriveHandler* pDrvHdl = 0;

        if (pDM != 0) {
            pDrvHdl = pDM->pGetDrvHdl();
        }

        DriveCaller *pDC = 0;

        if (!HP.IsKeyWord("bar" "width")) {
            silent_cerr("helical groove shape: keyword \"bar width\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal Ws = HP.GetReal();

        if (!HP.IsKeyWord("chamfer" "width")) {
            silent_cerr("helical groove shape: keyword \"chamfer width\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal Wc = HP.GetReal();

        if (!HP.IsKeyWord("groove" "depth")) {
            silent_cerr("helical groove shape: keyword \"groove depth\" expected at line " << HP.GetLineData() << std::endl);
            throw ErrGeneric(MBDYN_EXCEPT_ARGS);
        }

        const doublereal Hg = HP.GetReal();

        SAFENEWWITHCONSTRUCTOR(pDC,
                               GrooveShapeDriveCaller,
                               GrooveShapeDriveCaller(pDrvHdl, Ws, Wc, Hg));

        return pDC;
    }

}
#endif

bool hydrodynamic_plain_bearing2_set(void)
{
#if defined(USE_AUTODIFF) && __cplusplus >= 201103L
    UserDefinedElemRead *pElemRead = new UDERead<HydroRootElement>;

    if (!SetUDE("hydrodynamic" "plain" "bearing2", pElemRead))
    {
        delete pElemRead;
        return false;
    }

    DriveCallerRead *pDriveRead = new GrooveShapeDCR;

    if (!SetDriveCallerData("helical" "groove" "shape", pDriveRead)) {
        delete pDriveRead;
        return false;
    }

    return true;
#else
    return false;
#endif
}

#ifndef STATIC_MODULES

extern "C"
{
    int module_init(const char *module_name, void *pdm, void *php)
    {
        if (!hydrodynamic_plain_bearing2_set())
        {
            silent_cerr("hydrodynamic plain bearing2: "
                        "module_init(" << module_name << ") "
                        "failed" << std::endl);

            return -1;
        }

        return 0;
    }
}

#endif // ! STATIC_MODULE
