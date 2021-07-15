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
/* module-rotor_disc
 * Author: Matteo Daniele
 *
 * Copyright (C) 2008-2021
 *
 * Matteo Daniele <matteo.daniele@polimi.it>
 * 
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 */
// Author: Matteo Daniele <matteo.daniele@polimi.it>
// tail rotor module for mbdyn

// model as follower force composed of a part depending on xapedal and a part depending on twist and inflow
// T = T0*(T_theta+T_lambda_twist)
#include "mbconfig.h" 		/* This goes first in every *.c,*.cc file */

#include <iostream>
#include <cmath>
#include <iostream>
#include <fstream>

#include "module-rotor_disc.h"


RotorDisc::RotorDisc( unsigned int uLabel, const DofOwner *pDO,
                                DataManager* pDM, MBDynParser& HP)
: Elem(uLabel, flag(0)), UserDefinedElem(uLabel, pDO)
{
    if (HP.IsKeyWord("help")) {
        silent_cout("\nModule: rotor disc\n"
        "Author: Matteo Daniele <matteo.daniele@polimi.it>\n"
        "Organization:	Dipartimento di Ingegneria Aerospaziale			\n"
        "		Politecnico di Milano					\n"
        "		http://www.aero.polimi.it				\n"
        "									\n"
        "	All rights reserved						\n"
        "Implementation of a disc rotor model with closed form inflow calculation\n"
        "Input:\n"
        "   hub structural node label,\n"
        "   rotor relative arm wrt structural node,\n"
        "   (the rotor force is always oriented as the Z of the node),\n"
        "   control input driver (collective input [rad]),\n"
        "   air density driver (rho [kg/m3]),\n"
        "   rotor angular velocity driver [rad/s],\n"
        "   rotor radius [m],\n"
        "   disc area [m2],\n"
        "   rotor solidity [-],\n"
        "   blade ClAlpha [1/rad],\n"
        "   blade twist [rad],\n"
        "   rotor distance wrt reference point [m]\n"
        "   (e.g.: for tail rotor, the distance between main and tail hub centers),\n"
        "   alpha stall min (stall) [rad],\n"
        "   alpha stall max (stall) [rad],\n"
        "   control input minvalue (saturation) [rad],\n"
        "   control input maxvalue (saturation) [rad],\n"
        "   reference power  [shp],\n"
        "   (e.g.: for tail rotor, the main rotor nominal power),\n"
        "   reference rotor angular velocity [rad/s],\n"
        "   (e.g.: for tail rotor, the main rotor angular velocity),\n"
        "\n"
        "Output:\n"
        "   1)  element label\n"
        "   2-4)forces\n"
        "   5-7)moments\n"
        "   8)  thrust value  [N]\n"
        "   9)  induced drag  [N]\n"
        "   10) induced power [W]\n"
        "   11) pitch control input [rad]\n"
        << std::endl);

        if (!HP.IsArg()){
            throw NoErr(MBDYN_EXCEPT_ARGS);
        }
    }

    // read rotor hub node
    pHubNode = dynamic_cast<const StructNode*>(pDM->ReadNode(HP, Node::STRUCTURAL));
    ReferenceFrame rf;
    if (pHubNode)
    {
        rf = ReferenceFrame(pHubNode);
    }
    
    // distanza dal nodo (vettore di 3 elementi)
    // if wrong reference node label
    if (pHubNode == 0)
    {
        silent_cerr("rotordisc(" << uLabel << "): "
				"invalid node type at line " << HP.GetLineData() << std::endl);
			    throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    // when keyword "position" is indicated, the next 3 elements are the 
    // relative position components of the point of application of the force
    // wrt the node
    /*force: label, type,
        node_label,
        position, [vec3 indicating relative position | reference, smt, vec3 of relative pos wrt smt]
        [vec3 of relative direction + driver force | vec3 of drive callers in the rf of the node]
    */
    std::cout<< "reading rotor disc element..."<<std::endl;
    if (HP.IsKeyWord("position")){
        // obtain relative arm wrt node
        HubNodeArm = HP.GetPosRel(rf);
    }
    std::cout<< "rotordisc("<< uLabel << "): thrust relative position wrt node "<< pHubNode->GetLabel() <<" : "<< HubNodeArm <<std::endl;

    if (HP.IsKeyWord("orientation")){
        // obtain relative orientation wrt node
        RThrustOrientation = HP.GetRotRel(rf);
    }
    std::cout<<"rotordisc("<< uLabel << "): thrust relative orientation "<<std::endl;
    std::cout<<"wrt node "<< pHubNode->GetLabel() <<" : " <<std::endl;
    std::cout<< RThrustOrientation.GetRow(1) <<std::endl;
    std::cout<< RThrustOrientation.GetRow(2) <<std::endl;
    std::cout<< RThrustOrientation.GetRow(3) <<std::endl;
    // direction of thrust in the reference frame of hub node
    // the magnitude of the drivers is not important since it only
    // indicates the direction of the force wrt hub reference frame
    // pDCForceDir = ReadDCVecRel(pDM, HP, rf);    
    // initialize the drive owner of the force direction
    // DOForceDir.Set(pDCForceDir);

    // drive di input (collective pitch)
    pXColl  = HP.GetDriveCaller(); // da adoperare assieme a pXColl->dGet()
    // drive di input air density
    pRho    = HP.GetDriveCaller();
    // drive di input tail rotor angular speed;
    pOmega  = HP.GetDriveCaller();
    // dati di progetto del rotore
    RotorRadius             = HP.GetReal();
    DiscArea                = HP.GetReal();
    RotorSolidity           = HP.GetReal();
    ClAlpha                 = HP.GetReal();
    BladeTwist              = HP.GetReal();
    AOAStallMin             = HP.GetReal();
    AOAStallMax             = HP.GetReal();
    thetaCollMin            = HP.GetReal();
    thetaCollMax            = HP.GetReal();
    hubs_distance           = HP.GetReal();
    mr_nominal_power_shp    = HP.GetReal();
    mr_nominal_omega        = HP.GetReal();
    
    // pay attention to possible negative values
    if (hubs_distance <= std::numeric_limits<doublereal>::epsilon())
    {
        hubs_distance *= -1.0;
    }
    // compute nominal main rotor torque in hover
    doublereal mr_nominal_power_w   = mr_nominal_power_shp*sHP2W;
    doublereal mr_nominal_torque_Nm = mr_nominal_power_w/mr_nominal_omega;
    Th                              = mr_nominal_torque_Nm/hubs_distance;
    // v1h depending on costant rotor parameters
    v1hPart             = sqrt(Th/(2.0*DiscArea));
    doublereal v1hInit  = v1hPart/sqrt(1.225);

    std::cout <<"Tail rotor initialized:" << std::endl;
    std::cout <<"Radius [m]: " << RotorRadius << std::endl;
    std::cout <<"Area [m2]: " << DiscArea << std::endl;
    std::cout <<"Sigma [-]: " << RotorSolidity << std::endl;
    std::cout <<"ClAlpha [1/rad]: " << ClAlpha << std::endl;
    std::cout <<"twist [rad]: " << BladeTwist << std::endl;
    std::cout <<"AOAStallMin [rad]: " << AOAStallMin << std::endl;
    std::cout <<"AOAStallMax [rad]: " << AOAStallMax << std::endl;
    std::cout <<"thetaCollMin [rad]: " << thetaCollMin << std::endl;
    std::cout <<"thetaCollMax [rad]: " << thetaCollMax << std::endl;
    std::cout <<"ref distance [m]: " << hubs_distance << std::endl;
    std::cout <<"ref nominal power [rad]: " << mr_nominal_power_shp << std::endl;
    std::cout <<"ref omega [rad/s]: " << mr_nominal_omega << std::endl;
    std::cout <<"###########################" << std::endl;
    std::cout <<"Induced velocity in Hover [m/s]: " << v1hInit << std::endl;
    std::cout <<"Thrust required for Hover [N]: " << Th << std::endl;
    std::cout <<"###########################" << std::endl;


    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    out << "rotordisc: " << uLabel
        << " " << pHubNode->GetLabel() // node label
        << " " << HubNodeArm
        << " " << RotorRadius
        << " " << DiscArea
        << " " << RotorSolidity                
        << " " << ClAlpha              
        << " " << BladeTwist             
        << " " << AOAStallMin       
        << " " << AOAStallMax
        << " " << thetaCollMin
        << " " << thetaCollMax
        << " " << hubs_distance        
        << " " << mr_nominal_power_shp 
        << " " << mr_nominal_omega                  
        << std::endl;
}

RotorDisc::~RotorDisc()
{
    NO_OP;
}

void RotorDisc::Output(OutputHandler& OH) const
{
    if (bToBeOutput())
    {
        std::ostream& out = OH.Loadable();
        out << std::setw(8) << GetLabel() // 1: label
            << " " << F            // 2-4: force
            << " " << M            // 5-7: moment
            << " " << Thrust              //   8: thrust value [N]
            << " " << DragInduced         //   9: induced drag [N]
            << " " << PowerInduced        //  10:induced power [W]
            << " " << thetaColl  //  11: pitch control input [rad]
            << " " << rho          //  12: air density [kg/m3]
            << " " << RotorOmega  //  13: rotor angular speed [rad/s]
            << std::endl; 
    }
}


// partial derivatives for Jacobian assembly
void RotorDisc::dMuCalc()
{
    // Vtip: rotor tip speed
    // V: airspeed vector (u, v, w)
    // partial derivatives wrt u,v,w,omega,rho
    dMu[0] = u/(Vtip*Vtot);    // = dMudU;
    dMu[1] = v/(Vtip*Vtot);    // = dMudV;
    dMu[2] = w/(Vtip*Vtot);    // = dMudW;
    dMu[3] = -(Vtot*RotorRadius)/(Vtip2);  // = dMudOmega;
    dMu[4] = 0.0;                          // = dMudRho;
}

void RotorDisc::dLambdaCalc()
{
    doublereal v1h4 = pow(v1h, 4.0);
    doublereal dLambdaRhoDen = 2.0*Vtip*V1*a_v1*rho;

    dLambda[0]  = (0.5*u/(Vtip*V1))*(1.0-Vtot2/(2.0*a_v1));  // = dLambdadU;
    dLambda[1]  = (0.5*v/(Vtip*V1))*(1.0-Vtot2/(2.0*a_v1));  // = dLambdadV;
    dLambda[2]  = (0.5*w/(Vtip*V1))*(1.0-Vtot2/(2.0*a_v1));  // = dLambdadW;
    dLambda[3]  = (V1*RotorRadius)/Vtip2;                                                      // = dLambdadOmega;
    dLambda[4]  = v1h4/dLambdaRhoDen;                                                                                     // = dLambdadRho;
}

void RotorDisc::dT0Calc()
{
    // advance ratio
    doublereal c0 = 0.5*rho*Vtip2*DiscArea*RotorSolidity*ClAlpha;
    doublereal c1Num = -(3.0*mu);
    doublereal c1Den = pow((1.0+1.5*mu2),2.0); 
    doublereal c1 = c1Num/c1Den;
    doublereal c2 = c0*c1;

    // for dT0dOmega and dT0dRho
    doublereal c3 = 0.5*DiscArea*RotorSolidity*ClAlpha;
    doublereal c4 = 1.0+1.5*mu2;
    doublereal c5 = c4*2.0*Vtip*RotorRadius;
    doublereal c6 = -3.0*Vtip2*mu;

    doublereal dT0dOmegaNum = rho*c3*(c5+c6*dMu[3]);
    doublereal dT0dOmegaDen = (pow(c4,2.0));

    dT0[0] = c2*dMu[0];             // = dT0dU;
    dT0[1] = c2*dMu[1];             // = dT0dV;
    dT0[2] = c2*dMu[2];             // = dT0dW;
    dT0[3] = dT0dOmegaNum/dT0dOmegaDen;   // = dT0dOmega;
    dT0[4] = c3*Vtip2/(c4);  // = dT0dRho;
}

void RotorDisc::dTThetaCalc()
{
    doublereal p0 = thetaColl*(6*mu2-4/3)*mu;

    dTTheta[0] = p0*dMu[0]; // = dTThetadU;
    dTTheta[1] = p0*dMu[1]; // = dTThetadV; 
    dTTheta[2] = p0*dMu[2]; // = dTThetadW;
    dTTheta[3] = p0*dMu[3]; // = dTThetadOmega;
    dTTheta[4] = 0.0;       // = dTThetadRho;
}

void RotorDisc::dTLambdaCalc()
{
    doublereal l0 = 1.0-0.5*mu2;
    doublereal l1 = 1.5*(2.0*mu2-1.0)*mu*BladeTwist;
    doublereal l2 = l1-lambda*mu;

    dTLambda[0] = l0*dLambda[0] + l2*dMu[0];// = dTLambdadU    
    dTLambda[1] = l0*dLambda[1] + l2*dMu[1];// = dTLambdadV    
    dTLambda[2] = l0*dLambda[2] + l2*dMu[2];// = dTLambdadW    
    dTLambda[3] = l0*dLambda[3] + l2*dMu[3];// = dTLambdadOmega
    dTLambda[4] = l0*dLambda[4];            // = dTLambdadRho                
}

void RotorDisc::T0Calc()
{
    doublereal t0 = 0.5*rho*Vtip2*DiscArea*RotorSolidity*ClAlpha;
                    
    doublereal t1 = 1.0 + 1.5*mu2;
    // update T0
    T0 = t0/t1;

}

void RotorDisc::TThetaCalc()
{
    // update TTheta
    TTheta = (2.0/3.0 - 2.0/3.0*mu2 + 1.5*mu4)*thetaColl;
}

void RotorDisc::TLambdaCalc()
{
    doublereal i0 = (1.0-0.5*mu2);
    doublereal i1 = 0.5*(1.0-1.5*mu2 + 1.5*mu4);
    // update TLambda
    TLambda = i0*lambda + i1*BladeTwist;
}

void RotorDisc::ThrustCalc()
{
    //T = T0*(TTheta+TLambda);
    Thrust = (T0)*(TTheta+TLambda);
    // update induced power
    PowerInduced = Thrust*V1;
    // update induced drag
    DragInduced = PowerInduced/Vtot;

    OutputThrust[0] = 0.0;
    OutputThrust[1] = 0.0;
    OutputThrust[2] = Thrust;
}

void RotorDisc::dTCalc()
{
    // parameters constant for each partial derivative
    doublereal c0 = (TTheta + TLambda);
    doublereal c1 = T0;

    // update thrust partial derivatives
    int n = *(&Thrust + 1) - Thrust;
    for (int i=0; i<n; i++)
    {
        dThrust[i] = c0*dT0[i] + c1*(dTTheta[i] + dTLambda[i]);
    }
}

void RotorDisc::updateStatesDeps()
{
    u = VTrHub[0];
    v = VTrHub[1];
    w = VTrHub[2];

    Vtot = sqrt(pow(u, 2.0) + pow(v, 2.0) + pow(w, 2.0));

    Vtot2 = pow(Vtot, 2.0);
    // Vtip
    Vtip = RotorOmega*RotorRadius;
    Vtip2 = pow(Vtip,2.0);
    // adv ratio
    mu = Vtot/Vtip;
    mu2 = pow(mu,2.0);
    mu4 = pow(mu,4.0);
    // v1h: induced velocity in hover = sqrt(Th/2A)*sqrt(1/rho)
    v1h = v1hPart*sqrt(1.0/rho);
    // CONSTANT MOMENTUM INDUCED VELOCITY
    a_v1 = sqrt(pow(0.5*Vtot2,2.0) + pow(v1h, 4.0));
    V1 = sqrt(-0.5*Vtot2 + a_v1);
    
    // inflow ratio
    lambda = -(V1/Vtip);

}

void RotorDisc::inputSaturation()
{
    // std::cout<< "ThetaCollBeforeSat = "<< thetaColl*180.0/M_PI << std::endl;
    // check if collective input saturation 
    // the minimum value between input and saturation upper limit
    thetaColl = std::min(thetaColl, thetaCollMax);
    // the maximum value between input and saturation lower limit
    thetaColl = std::max(thetaColl, thetaCollMin);
    // std::cout<< "ThetaCollAfterSat = "<< thetaColl*180.0/M_PI << std::endl;

    
}

void RotorDisc::assemblyJacobian()
{
    // compute partial derivatives of adv ratio
    dMuCalc();
    // compute partial derivatives of inflow ratio
    dLambdaCalc();
    // compute partial derivatives components of the jacobian 
    dT0Calc();
    dTThetaCalc();
    dTLambdaCalc();
    // compute the components of thrust (and thrust itself)
    computeRotorThrust(); // gives also the value of induced power and drag
    // and finally compute Jacobian of tail rotor thrust
    dTCalc();
}


void RotorDisc::computeRotorThrust()
{
    // compute the components of thrust (and thrust itself)
    T0Calc();
    TThetaCalc();
    TLambdaCalc();
    // compute thrust
    ThrustCalc();
}

////////////////////////////////////////////////////////////////////////////////////
void RotorDisc::WorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 6;
	*piNumCols = 3;
}

void RotorDisc::InitialWorkSpaceDim(integer* piNumRows, integer* piNumCols) const
{
	*piNumRows = 12;
	*piNumCols = 6;
}

//contributo al file di restart
std::ostream& RotorDisc::Restart(std::ostream& out) const
{
    return out << "# not implemented yet" << std::endl;
}

VariableSubMatrixHandler&
RotorDisc::AssJac(VariableSubMatrixHandler& Workmat,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
    // sicuri sicuri che qua non va niente?
    Workmat.SetNullMatrix();
    return Workmat;
    /*
    FullSubMatrixHandler& WM = Workmat.SetFull();
    // reset and workmatrix dimension
    integer iNumRows = 0;
    integer iNumCols = 0;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WM.ResizeReset(iNumRows, iNumCols);

    integer iFirstRotationIndex = pHubNode->iGetFirstPositionIndex() + 3;
    integer iFirstMomentumIndex = pHubNode->iGetFirstMomentumIndex();
    for (integer iCnt = 1; iCnt <= 3; iCnt++)
    {
        WM.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt); // force
        WM.PutRowIndex(3 + iCnt, iFirstMomentumIndex + 3 + iCnt); // couple
        WM.PutColIndex(iCnt, iFirstRotationIndex + iCnt); // rotation
    }

    

    // data
    const Mat3x3& R(pHubNode->GetRRef());
    Vec3 TmpDir(R*(f.Get()*dCoef));
    Vec3 TmpArm(R*HubNodeArm);
    /// |    F/\   |           |   F  |
	 // |          | Delta_g = |      |
	 // | (d/\F)/\ |           | d/\F |
	 ///
    WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
    WM.Add(4, 1, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));

    return Workmat;
    */
    
}

unsigned int RotorDisc::iGetNumPrivData(void) const
{
    // number of private data that can be extracted from the module
    // Thrust      
    // DragInduced 
    // PowerInduced
    // thetaColl
    // rho         
    // RotorOmega  
    return 6;
}

unsigned int RotorDisc::iGetPrivDataIdx(const char* s) const
{

    ASSERT(s != NULL);

    struct 
    {
        const char* s;
        int i;
    } sPrivData[] = {
        {"Thrust", THRUST},
        {"DragInd", DRAGINDUCED},
        {"PowerInd", POWERINDUCED},
        {"theta0", THETA0},
        {"rho", RHO},
        {"omega", OMEGA},
        {0}
    };

    for (int i = 0; sPrivData[i].s != 0; i++)
    {
        if (strcasecmp(s, sPrivData[i].s) == 0)
        {
            return sPrivData[i].i;
        }
    }

    return 0;
    

    //unsigned idx = 0;
    //switch(s[0])
    //{
    //    case 'T'://hrust
    //        idx += 1;
    //        break;
    //    case 'D'://ragInduced
    //        idx += 2;
    //        break;
    //    case 'P'://owerInduced
    //        idx += 3;
    //        break;
    //    case 't'://heta0
    //        idx += 4;
    //        break;
    //    case 'r'://ho
    //        idx += 5;
    //        break;
    //    case 'o'://mega
    //        idx += 6;
    //        break;
    //    default:
    //        return 0;
    //}
    //return idx;
}

doublereal RotorDisc::dGetPrivData(unsigned int i) const
{
    if (i <= 0 || i >= LASTPRIVDATA)
    {
        silent_cerr("rotordisc("<<GetLabel()<<"): "
        "private data "<< i << "not available" << std::endl);
        throw ErrGeneric(MBDYN_EXCEPT_ARGS);
    }
    else
    {
        switch (i)
        {
            case THRUST         : return Thrust;
            case DRAGINDUCED    : return DragInduced;
            case POWERINDUCED   : return PowerInduced;
            case THETA0         : return thetaColl;
            case RHO            : return rho;
            case OMEGA          : return RotorOmega;
        }
    }

    return 0.;
}

int RotorDisc::iGetNumConnectedNodes(void) const
{
    // hub node
    return 1;
}

void RotorDisc::GetConnectedNodes(std::vector<const Node *>& connectedNodes) const
{
    // hub node
    connectedNodes.resize(1);
    connectedNodes[0] = pHubNode;
}

void RotorDisc::SetValue(DataManager *pDM,
	VectorHandler& X, VectorHandler& XP,
	SimulationEntity::Hints *ph)
{
   	NO_OP;
}

unsigned RotorDisc::iGetInitialNumDof(void) const
{
    return 0;
}

SubVectorHandler&
RotorDisc::AssRes(SubVectorHandler& WorkVec,
                        doublereal dCoef,
                        const VectorHandler& XCurr,
                        const VectorHandler& XPrimeCurr)
{
    integer iNumRows;
    integer iNumCols;
    WorkSpaceDim(&iNumRows, &iNumCols);
    WorkVec.ResizeReset(iNumRows);
    // indices of node unkowns
    integer iFirstMomentumIndex = pHubNode->iGetFirstMomentumIndex();
    for (integer iCnt = 1; iCnt <= 6; iCnt++)
     {
         WorkVec.PutRowIndex(iCnt, iFirstMomentumIndex + iCnt);
     }

    // get velocity of the hub point in the reference frame of the node
    Vec3 VTrHubTemp = pHubNode->GetVCurr() + pHubNode->GetWCurr().Cross(HubNodeArm);
    // rotate the velocity in the reference frame of the rotor disc
    VTrHub = RThrustOrientation*VTrHubTemp;
    // the force is oriented along the ideal z of node reference frame
    //VTrHub = pHubNode->GetVCurr() + pHubNode->GetWCurr().Cross(HubNodeArm);
    // air density and rotor angular velocity
    RotorOmega  = pOmega->dGet();
    rho         = pRho->dGet();
    // rotor collective pitch input
    thetaColl   = pXColl->dGet();
    
    // check for saturation
    inputSaturation();
    // update state-dependent parameters (lambda, mu, ecc)
    updateStatesDeps();
    // compute rotor thrust (direction wrt to node hub given by constructor)
    computeRotorThrust();
    // OutputThrust = Zero3;
    const Mat3x3& R(pHubNode->GetRCurr());
    //Vec3 TmpDir = f.Get();
    //Vec3 F(R*TmpDir);
    //Vec3 M(R*HubNodeArm.Cross(TmpDir));
    // bring back rotor disc forces to node reference frame
    OutputThrust = RThrustOrientation.MulTV(OutputThrust);
    // std::cout<< "T0 = "<< T0 << std::endl;
    // std::cout<< "TTheta = "<< TTheta << std::endl;
    // std::cout<< "TLambda = "<< TLambda << std::endl;
    
    //OutputThrust = RThrustOrientation*OutputThrust;
    // from node reference frame to global reference frame
    F = R*OutputThrust;
    M = R*HubNodeArm.Cross(OutputThrust);

    WorkVec.Add(1, F);
    WorkVec.Add(4, M);

    return WorkVec;
}

VariableSubMatrixHandler& 
RotorDisc::InitialAssJac(VariableSubMatrixHandler& WorkMat,
                                const VectorHandler& XCurr)
{
    WorkMat.SetNullMatrix();
    return WorkMat;
    /*
    FullSubMatrixHandler& WM = WorkMat.SetFull();
    // workmatrix reset
    WM.ResizeReset(12, 6);

    integer iFirstPositionIndex = pHubNode->iGetFirstPositionIndex();
    integer iFirstVelocityIndex = iFirstPositionIndex + 6;
    for (integer iCnt = 1; iCnt <= 3; iCnt)
    {
        WM.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WM.PutRowIndex(3 + iCnt, iFirstPositionIndex + 3 + iCnt);
		WM.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
		WM.PutRowIndex(9 + iCnt, iFirstVelocityIndex + 3 + iCnt);
		WM.PutColIndex(iCnt, iFirstPositionIndex + 3 + iCnt);
		WM.PutColIndex(3 + iCnt, iFirstVelocityIndex + 3 + iCnt);
    }

    // data
	const Mat3x3& R(pHubNode->GetRRef());
	Vec3 TmpArm(R*HubNodeArm);
	Vec3 TmpDir = R*f.Get();
	const Vec3& Omega(pHubNode->GetWRef());

	/// |    F/\   |           |   F  |
	 // |          | Delta_g = |      |
	 // | (d/\F)/\ |           | d/\F |
	 ///

	WM.Add(1, 1, Mat3x3(MatCross, TmpDir));
	WM.Add(4, 1, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));
	WM.Add(7, 1, Mat3x3(MatCrossCross, Omega, TmpDir));
	WM.Add(7, 4, Mat3x3(MatCross, TmpDir));
	WM.Add(10, 1, Mat3x3(MatCrossCross, Omega, TmpArm.Cross(TmpDir)));
	WM.Add(10, 4, Mat3x3(MatCross, TmpArm.Cross(TmpDir)));

	return WorkMat;
    */
}

SubVectorHandler&
RotorDisc::InitialAssRes(SubVectorHandler& WorkVec, const VectorHandler& XCurr)
{
    WorkVec.Resize(0);
    return WorkVec;
    /*
    integer iNumRows;
	integer iNumCols;
	InitialWorkSpaceDim(&iNumRows, &iNumCols);
	WorkVec.ResizeReset(iNumRows);

	// Indici delle incognite del nodo
	integer iFirstPositionIndex = pHubNode->iGetFirstPositionIndex();
	integer iFirstVelocityIndex = iFirstPositionIndex + 6;
	for (integer iCnt = 1; iCnt <= 6; iCnt++) {
		WorkVec.PutRowIndex(iCnt, iFirstPositionIndex + iCnt);
		WorkVec.PutRowIndex(6 + iCnt, iFirstVelocityIndex + iCnt);
	}

	//Dati
	const Mat3x3& R(pHubNode->GetRCurr());
	Vec3 TmpDir(R*f.Get());
	Vec3 TmpArm(R*HubNodeArm);
	const Vec3& Omega(pHubNode->GetWCurr());

	WorkVec.Add(1, TmpDir);
	WorkVec.Add(4, TmpArm.Cross(TmpDir));
	WorkVec.Add(7, Omega.Cross(TmpDir));
	WorkVec.Add(10, (Omega.Cross(TmpArm)).Cross(TmpDir)
		+ TmpArm.Cross(Omega.Cross(TmpDir)));

	return WorkVec;
    */
}

bool RotorDisc_set(void)
{
    UserDefinedElemRead *rf = new UDERead<RotorDisc>;

    if (!SetUDE("rotordisc", rf))
    {
        delete rf;
        return false;
    }

    return true;
}

//#ifndef STATIC_MODULES
extern "C" int
module_init(const char *module_name, void *pdm, void *php)
{
    if (!RotorDisc_set())
    {
        silent_cerr("rotordisc: "
                    "module_init("<< module_name << ") "
                    "failed" << std::endl);
        return -1;
    }

    return 0;

}
//#endif // ! STATIC_MODULES