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
#include <fstream>

#include "module-rotor_disc.h"


RotorDisc::RotorDisc( unsigned int uLabel, const DofOwner *pDO,
                                DataManager* pDM, MBDynParser& HP)
//: AerodynamicElem(uLabel, pDO, flag(0))
//: Elem(uLabel, flag(0)), AerodynamicElem(uLabel, pDO, flag(0)), UserDefinedElem(uLabel, pDO)
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
        "   rotor relative position wrt structural node,\n"
        "   (the rotor force is always oriented as the Z of the reference frame considered),\n"
        "   control input driver (collective input [rad]),\n"
        "   rotor angular velocity driver [rad/s],\n"
        "   rotor radius [m],\n"
        "   rotor solidity [-],\n"
        "   blade Cl0 [1/rad],\n"
        "   blade ClAlpha [1/rad],\n"
        "   blade twist [rad],\n"
        "   alpha stall min (stall) [rad],\n"
        "   alpha stall max (stall) [rad],\n"
        "   control input minvalue (saturation) [rad],\n"
        "   control input maxvalue (saturation) [rad],\n"
        "   [main rotor data | MTOW]\n"
        "   |____ main rotor data,\n"
        "   |       hubs distance,\n"
        "   |       main rotor nominal power,\n"
        "   |       main rotor angular speed,\n"
        "   |____ MTOW\n"
        "\n"
        "Output:\n"
        "   1)      element label\n"
        "   2-4)    forces\n"
        "   5-7)    moments\n"
        "   8)      thrust value [N]\n"
        "   9)      induced drag [N]\n"
        "  10)      induced power [W]\n"
        "  11)      pitch control input [rad]\n"
        "  12)      air density [kg/m3]\n"
        "  13)      rotor angular speed [rad/s]\n"
        "  14)      alpha tip pat plane [rad]\n"
        "  15)      inflow ratio [-]\n"
        "  16)      inflow ratio - iterative [-]\n"
        "  17)      constant induced velocity on rotor disc [m/s]\n"
        "  18-20)   hub velocity [m/s]\n"
        "  21)      rotor induced velocity in hover [m/s]\n"
        "  22)      rotor induced power in hover [W]\n"
        "  23)      CT [-]\n"
        "  24)      CP [-]\n"
        "  25)      FOM [-]\n"
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
    DEBUGCOUT( "reading rotor disc element..."<<std::endl);
    if (HP.IsKeyWord("position")){
        // obtain relative arm wrt node
        HubNodeArm = HP.GetPosRel(rf);
    }
    DEBUGCOUT( "rotordisc("<< uLabel << "): thrust relative position wrt node "<< pHubNode->GetLabel() <<" : "<< HubNodeArm <<std::endl);

    if (HP.IsKeyWord("orientation")){
        // obtain relative orientation wrt node
        RThrustOrientation = HP.GetRotRel(rf);
    }
    DEBUGCOUT("rotordisc("<< uLabel << "): thrust relative orientation "<<std::endl);
    DEBUGCOUT("wrt node "<< pHubNode->GetLabel() <<" : " <<std::endl);
    DEBUGCOUT( RThrustOrientation.GetRow(1) <<std::endl);
    DEBUGCOUT( RThrustOrientation.GetRow(2) <<std::endl);
    DEBUGCOUT( RThrustOrientation.GetRow(3) <<std::endl);
    // direction of thrust in the reference frame of hub node
    // the magnitude of the drivers is not important since it only
    // indicates the direction of the force wrt hub reference frame
    // pDCForceDir = ReadDCVecRel(pDM, HP, rf);
    // initialize the drive owner of the force direction
    // DOForceDir.Set(pDCForceDir);
    // drive di input (collective pitch)
    if (HP.IsKeyWord("collective" "input")){
        pXColl  = HP.GetDriveCaller(); // da adoperare assieme a pXColl->dGet()
    }
    // drive di input air density
    // pRho    = HP.GetDriveCaller();

    // air density taken by exploiting class inheritance on aerodynamic element
    // const Vec3& XPosAbs(pHubNode->GetXCurr());
    DEBUGCOUT("rotordisc("<< uLabel << "): rho = "<< rho <<std::endl);
    rho = 1.225 ; // dGetAirDensity(pHubNode->GetXCurr());
    DEBUGCOUT("rotordisc("<< uLabel << "): rho = "<< rho <<std::endl);

    bool bGotOmega      = false;
    bool bGotChord      = false;
    bool bGotNBlades    = false;
    bool bGotSigma      = false;
    // drive di input rotor angular speed;
    if (HP.IsKeyWord("angular" "velocity")){
        pOmega  = HP.GetDriveCaller();
        bGotOmega = true;
    }
    // dati di progetto del rotore
    if (HP.IsKeyWord("radius")){
        RotorRadius             = HP.GetReal();
    }


    // TODO: AGGIUNGERE METODO CHE DECIDA COME SI VOGLIONO CALCOLARE I DATI
    if (HP.IsKeyWord("chord")){
        Chord = HP.GetReal();
        bGotChord = true;
    }
    if (HP.IsKeyWord("blade" "number")){
        NBlades = HP.GetInt();
        bGotNBlades = true;
    }
    // if i have both then useless to look for solidity in the data
    if (bGotChord & bGotNBlades){
        RotorSolidity = doublereal(NBlades)*Chord/(M_PI*RotorRadius);
    }
    else if (HP.IsKeyWord("sigma")){
        RotorSolidity = HP.GetReal();
    }

    if (HP.IsKeyWord("cl0")){
        Cl0 = HP.GetReal();
    }
    else { 
        std::cout << "rotordisc(" << uLabel << "): Cl0 not provided, assuming null Cl0"<< std::endl;
        Cl0 = 0.0;
    }

    if (HP.IsKeyWord("clalpha")){
        ClAlpha = HP.GetReal();
    }

    if (HP.IsKeyWord("twist")){
        BladeTwist = HP.GetReal();
    }
    else {
        std::cout << "rotordisc(" << uLabel << "): twist not provided, assuming null twist"<< std::endl;
        BladeTwist = 0.0;
    }

    if (HP.IsKeyWord("stall" "limits")){
        AOAStallMin             = HP.GetReal();
        AOAStallMax             = HP.GetReal();
    }
    else {
        std::cout << "rotordisc(" << uLabel << "): stall limits not provided, stall will not be included in the model"<< std::endl;
    }

    if (HP.IsKeyWord("control" "limits")){
        thetaCollMin            = HP.GetReal();
        thetaCollMax            = HP.GetReal();
    }
    else {
        std::cout << "rotordisc(" << uLabel << "): control limits not provided, saturation not included"<< std::endl;
    }

    // disc area
    DiscArea = M_PI*pow(RotorRadius, 2.0);

    // initialize the stall effects on the rotor
    computeCLInit();

    bool bGotTailRotor = false;
    bool bGotMainRotor = false;

    // reference values if we are dealing with a tail rotor
    if (HP.IsKeyWord("main" "rotor" "data")){
        if (HP.IsKeyWord("hubs" "distance")){
            hubs_distance = HP.GetReal();
            // pay attention to possible negative values
            hubs_distance = std::abs(hubs_distance);
        }
        if (HP.IsKeyWord("main" "rotor" "nominal" "power")){
            mr_nominal_power_shp    = HP.GetReal();
        }
        if (HP.IsKeyWord("main" "rotor" "nominal" "angular" "speed")){
            mr_nominal_omega        = HP.GetReal();
            if (bGotOmega == false){
                std::cout << "rotordisc(" << uLabel << "): tail rotor angular speed not provided: "<< std::endl;
                std::cout << "using as default value: tail_rotor_omega = 4*main_rotor_omega"<< std::endl;
                RotorOmega = 4.0*mr_nominal_omega;
            }
        }


        // compute nominal main rotor torque in hover
        doublereal mr_nominal_power_w   = mr_nominal_power_shp*sHP2W;
        doublereal mr_nominal_torque_Nm = mr_nominal_power_w/mr_nominal_omega;
        Th                              = mr_nominal_torque_Nm/hubs_distance;

        bGotTailRotor = true;

    }
    // second case: we are dealing with a disc rotor, we need to compute the thrust in hover
    // at sea level to find the induced velocity in hover
    else if (HP.IsKeyWord("MTOW")){

        MTOW = HP.GetReal();
        Th = MTOW*g;

        bGotMainRotor = true;
    }

    // v1h depending on costant rotor parameters
    v1hPart             = sqrt(Th/(2.0*DiscArea));
    doublereal v1hInit  = v1hPart/sqrt(1.225);
    // induced power in hover
    Wh                  = Th*v1hInit;

    DEBUGCOUT("Tail rotor initialized:" << std::endl);
    DEBUGCOUT("Radius [m]: " << RotorRadius << std::endl);
    DEBUGCOUT("Area [m2]: " << DiscArea << std::endl);
    DEBUGCOUT("Sigma [-]: " << RotorSolidity << std::endl);
    DEBUGCOUT("Cl0 [-]: " << Cl0 << std::endl);
    DEBUGCOUT("ClAlpha [1/rad]: " << ClAlpha << std::endl);
    DEBUGCOUT("twist [rad]: " << BladeTwist << std::endl);
    DEBUGCOUT("AOAStallMin [rad]: " << AOAStallMin << std::endl);
    DEBUGCOUT("AOAStallMax [rad]: " << AOAStallMax << std::endl);
    DEBUGCOUT("thetaCollMin [rad]: " << thetaCollMin << std::endl);
    DEBUGCOUT("thetaCollMax [rad]: " << thetaCollMax << std::endl);
    if (bGotTailRotor){
        std::cout <<"ref distance [m]: " << hubs_distance << std::endl;
        std::cout <<"ref nominal power [rad]: " << mr_nominal_power_shp << std::endl;
        std::cout <<"ref omega [rad/s]: " << mr_nominal_omega << std::endl;
    }
    else if (bGotMainRotor){
        std::cout <<"MTOW [kg]: " << MTOW << std::endl;
    }
    std::cout <<"rotordisc(" << uLabel << "): Initial Properties"<< std::endl;
    std::cout <<"###########################" << std::endl;
    std::cout <<"Induced velocity in Hover [m/s]: " << v1hInit << std::endl;
    std::cout <<"Thrust required for Hover [N]: " << Th << std::endl;
    std::cout <<"###########################" << std::endl;


    SetOutputFlag(pDM->fReadOutput(HP, Elem::LOADABLE));
    std::ostream& out = pDM->GetLogFile();
    if (bGotTailRotor){
        out << "rotordisc: " << uLabel
            << " " << pHubNode->GetLabel() // node label
            << " " << HubNodeArm
            << " " << RotorRadius
            << " " << DiscArea
            << " " << RotorSolidity
            << " " << Cl0
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
    else if (bGotMainRotor){
        out << "rotordisc: " << uLabel
            << " " << pHubNode->GetLabel() // node label
            << " " << HubNodeArm
            << " " << RotorRadius
            << " " << DiscArea
            << " " << RotorSolidity
            << " " << Cl0
            << " " << ClAlpha
            << " " << BladeTwist
            << " " << AOAStallMin
            << " " << AOAStallMax
            << " " << thetaCollMin
            << " " << thetaCollMax
            << " " << MTOW
            << std::endl;
    }
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
            << " " << F            	// 2-4: force
            << " " << M            	// 5-7: moment
            << " " << Thrust        //   8: thrust value [N]
            << " " << DragInduced   //   9: induced drag [N]
            << " " << PowerInduced  //  10: induced power [W]
            << " " << thetaColl  	//  11: pitch control input [rad]
            << " " << rho          	//  12: air density [kg/m3]
            << " " << RotorOmega  	//  13: rotor angular speed [rad/s]
            << " " << alphaTPP  	//  14: alpha tip pat plane [rad]
            << " " << lambda  		//  15: inflow ratio [-]
            << " " << lambdaNewman	//  16: inflow ratio - iterative [-]
            << " " << V1  		    //  17: constant induced velocity on rotor disc [m/s]
            << " " << VTrHub        //  18-20: hub velocity [m/s]
            << " " << v1h           //  21: rotor induced velocity in hover [m/s]
            << " " << Wh            //  22: rotor induced power in hover [W]
            << " " << CT            //  23: CT [-]
            << " " << CP            //  24: CP [-]
            << " " << FOM           //  25: FOM [-]
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
    doublereal c0 = 0.5*rho*Vtip2*DiscArea*RotorSolidity*Cl;
    doublereal c1Num = -(3.0*mu);
    doublereal c1Den = pow((1.0+1.5*mu2),2.0);
    doublereal c1 = c1Num/c1Den;
    doublereal c2 = c0*c1;

    // for dT0dOmega and dT0dRho
    doublereal c3 = 0.5*DiscArea*RotorSolidity*Cl;
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
    doublereal t0 = 0.5*rho*Vtip2*DiscArea*RotorSolidity*Cl;

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
    // WARNING: THRUST IS APPLIED AS A FOLLOWER FORCE APPLIED ON THE LOCAL Z AXIS,
    // THE SIGN HAS YO BE CHANGED TO TAKE INTO ACCOUNT THE CORRECT DIRECTION OF APPLICATION
    OutputThrust[2] = - Thrust;

    CT      = abs(Thrust)/(rho*DiscArea*Vtip2);
    CP      = pow(CT,1.5)/sqrt(2.0);
    CP == 0.0 ? FOM = 0.0 : FOM = pow(CT,1.5)/(CP*sqrt(2.0));

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
    mu  = Vtot/Vtip;
    mu2 = pow(mu,2.0);
    mu4 = pow(mu,4.0);
    // v1h: induced velocity in hover = sqrt(Th/2A)*sqrt(1/rho)
    v1h = v1hPart*sqrt(1.0/rho);
    // CONSTANT MOMENTUM INDUCED VELOCITY
    a_v1 = sqrt(pow(0.5*Vtot2,2.0) + pow(v1h, 4.0));
    V1 = sqrt(-0.5*Vtot2 + a_v1);
    // alpha tip path plane
    alphaTPP = atan2(w,u);
    // inflow ratio
    lambda = -(V1/Vtip);

}

void RotorDisc::computeLambdaNewman(){

    
    doublereal vsx,vsz; // components of the inflow along x,z
    doublereal dfdl, f, lambdaIt, lambdaItPrev; // function used to find lamba in newton iterations
    
    doublereal tollLNMax = 1.e-5; // tolerance
    int nMax = 15;      // max number of iterations
    
    doublereal ct2 = 0.5*CT;
    // initial guess is lambda just found out
    lambdaItPrev = 0.0;
    lambdaIt = lambda;
    vsx = mu*cos(alphaTPP);

    doublereal tollLn = 1.0;
    int it=0;
    while ((it<nMax) and (tollLn>=tollLNMax)){
        
        vsz = mu*sin(alphaTPP)+lambdaIt;
        f = lambdaIt-ct2*pow((pow(vsx,2.0)+pow(vsz,2.0)),-0.5);
        dfdl = 1+ct2*vsz*pow((pow(vsx,2.0)+pow(vsz,2.0)),-1.5);
        
        lambdaIt = lambdaItPrev - f/dfdl;

        tollLn = abs(lambdaIt-lambdaItPrev);
        lambdaItPrev = lambdaIt;
        it++;
        
    }

    lambdaNewman = lambdaIt;

    //std::cout << "LAMBDA\tLAMBDA_NEWMAN" << std::endl;
    //std::cout << lambda << "\t" << lambdaNewman << std::endl;

}

// compute cl with stall effects
void RotorDisc::computeCLInit(doublereal RDecayIn)
{
    Cl = ClAlpha*thetaColl+Cl0;

    // stall effect: circumference tangent to linear cl in clmax point
    clMax = ClAlpha*AOAStallMax+Cl0;
    clMin = ClAlpha*AOAStallMin+Cl0;

    xpMin = AOAStallMin;
    ypMin = clMin;
    xpMax = AOAStallMax;
    ypMax = clMax;

    // stall decay + "delay" modeled by circumference of radius R (assigned here)
    RDecay=RDecayIn*acos(-1.0)/180.0;
    // TODO: perfezionare stall decay+delay
    a1=atan(-1/ClAlpha);
    // circ center (max alpha)
    x0Min = xpMin-RDecay*cos(a1);
    y0Min = ypMin-RDecay*sin(a1);

    x0Max = xpMax+RDecay*cos(a1);
    y0Max = ypMax+RDecay*sin(a1);

    xbMin = xpMin-2*RDecay*cos(a1);
    ybMin = ypMin;

    xbMax = xpMax+2*RDecay*cos(a1);
    ybMax = ypMax;

    mb = 1/tan(a1);
    qbMin = ybMin-mb*xbMin;
    qbMax = ybMax-mb*xbMax;

    AOAAfterDecayMin = xbMin;
    AOAAfterDecayMax = xbMax;
    
}

void RotorDisc::computeCL(){

    // for simplicity, the input is thetacoll
    if (thetaColl < AOAAfterDecayMin){
        // after decay min phase
        Cl=mb*thetaColl+qbMin;
    }
    else if ((AOAAfterDecayMin <= thetaColl) and (thetaColl < AOAStallMin)){
        // decay min phase
        Cl=y0Min-sqrt(pow(RDecay,2.0)-pow((thetaColl-x0Min),2.0));
    }
    else if (thetaColl > AOAAfterDecayMax){
        // after decay max phase
        Cl=mb*thetaColl+qbMax;
    }
    else if ((AOAStallMax < thetaColl) and (thetaColl <= AOAAfterDecayMax)){
        // decay max phase
        Cl=y0Max+sqrt(pow(RDecay,2.0)-pow((thetaColl-x0Max),2.0));
    }
    else {
        // linear field
        Cl=ClAlpha*thetaColl+Cl0;
    }

    // Cl=ClAlpha*thetaColl+Cl0;

    
    //std::cout << " " << thetaColl << 
    //            //" " << ClAlpha <<
    //            //" " << Cl0 <<
    //            //" " << x0Max <<" " << x0Min <<
    //            //" " << y0Max <<" " << y0Min <<
    //            " " << AOAAfterDecayMin <<
    //            " " << AOAAfterDecayMax <<
    //            " " << AOAStallMin <<
    //            " " << AOAStallMax <<
    //            " " << Cl << std::endl;

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
    // DragInd
    // PowerInd
    // theta0
    // rho
    // omega
    // alphatpp
    // lambda
    // vind
    // vinh
    // pindh
    // ct
    // cp
    // fom
    return 14;
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
        {"alphatpp", ATPP},
        {"lambda", LAMBDA},
        {"vind", VIND},
        {"vinh", VINDH},
        {"pindh", PINDH},
        {"ct", CTPRIV},
        {"cp", CPPRIV},
        {"fom", FOMPRIV},
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
            case ATPP           : return alphaTPP;
            case LAMBDA         : return lambda;
            case VIND           : return V1;
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
    // air density taken by exploiting class inheritance on aerodynamic element
    rho = dGetAirDensity(pHubNode->GetXCurr());
    // rotor collective pitch input
    thetaColl   = pXColl->dGet();
    // check for saturation
    inputSaturation();
    // compute CL with stall effects
    computeCL();
    // update state-dependent parameters (lambda, mu, ecc)
    updateStatesDeps();
    // compute rotor thrust (direction wrt to node hub given by constructor)
    computeRotorThrust();
    // compute lambda newman (as reference)
    computeLambdaNewman();
    // OutputThrust = Zero3;
    const Mat3x3& R(pHubNode->GetRCurr());
    // bring back rotor disc forces to node reference frame
    OutputThrust = RThrustOrientation.MulTV(OutputThrust);
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
