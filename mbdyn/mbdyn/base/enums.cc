/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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

/* Dichiarazione ed inizializzazione delle stringhe relative agli enums */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <stdlib.h>

/* Tipi di elementi */
const char* psElemNames[] = {
   "Air Properties",
   "Rotor",

   "Automatic Structural",

   "Gravity",
   "Rigid Body",
   "Joint",
	"Joint Regularization",
   "Beam",
   "Plate",
   
   "Force",
   
   "Inertia",
   
   "Electric Bulk",
   "Electric",

   "Thermal",
   
   "Hydraulic",
   
   "Bulk",
   "Loadable",
   "Driven",
   "External",

   "Aeromodal",
   "Aerodynamic Element",

   "GENEL",

   "Socket Stream Output",

   NULL
};

/* Tipi di elementi */
const char* psReadElemsElems[] = {
   "air" "properties",
   "rotor",

   "automatic" "structural",

   "gravity",
   "body",
   "joint",
	"joint" "regularization",
   "beam",
   "plate",
   
   "force",
   
   "inertia",
   
   "electric" "bulk",
   "electric",

   "thermal",
   
   "hydraulic",
   
   "bulk",
   "loadable",
   "driven",
   "external",

   "aeromodal",
   "aerodynamic" "element",
   
   "genel",

   "socket" "stream" "output",

   NULL
};


/* Nomi degli elementi quando vengono dichiarati nel blocco dati di controllo;
 * hanno lo stesso numero progressivo degli elementi, quindi vi si accede con 
 * l'enum Elem::Type */
const char* psReadControlElems[] = {
   "air" "properties",
   "induced" "velocity" "elements",

   "automatic" "structural" "elements",

   "gravity",
   "rigid" "bodies",
   "joints",
	"joint" "regularizations",
   "beams",
   "plates",
   
   "forces",

   "inertia",
   
   "electric" "bulk" "elements",
   "electric" "elements",

   "thermal" "elements",
   
   "hydraulic" "elements",
   
   "bulk" "elements",
   "loadable" "elements",
   "driven" "elements",
   "external" "elements",

   "aeromodals",
   "aerodynamic" "elements",
   
   "genels",

   "output" "elements",

   NULL
};


/* Codici Adams per denominazione entita' */
const char *psAdamsElemCode[] = {
   "AIRP",
   "ROTR",

   "AUTO",

   "FORC",
   
   "GRAV",
   "PART",
   "LINK",
   "BEAM",
   "PLAT",
   
   "ELBK",
   "ELEC",

   "THER",
   
   "HYDR",
   
   "BULK",
   "LOAD",
   "DRVN",

   "RMBX",
   "EXTN",
   
   "AERO",
   "AERM",
   
   "GENL",

   "OUTP",

   NULL
};


/* Tipi di Joint */
const char* psJointNames[] = {
   "Distance",
   "Distance With Offset",
   "Clamp",
   "Spherical Hinge",
   "Pin",
   "Cardano Hinge",
   "Cardano Rotation",
   "Cardano Pin",
   "Plane Hinge",
   "Plane Rotation",
   "Plane Pin",
   "Axial Rotation",
   "Plane Displacement Hinge",
   "Plane Displacement Pin",
   "In Plane",
   "In Plane Contact",
   "Inline",
   "Rod",
   "Deformable Hinge",
   "Deformable Displacement Joint",
   "Deformable Joint",
   "Deformable Axial Joint",
   "Viscous Body",
   "Linear Velocity",
   "Angular Velocity",
   "Linear Acceleration",
   "Angular Acceleration",
   "Prismatic",
   "Drive Hinge",
   "Drive Displacement",
   "Drive Displacement Pin",
   "Imposed Displacement",
   "Imposed Displacement Pin",
   "Imposed Kinematics",
   "Beam Slider",
   "Brake",
   "Gimbal Rotation",
   "Point Surface Contact",
   "Total Joint",
   "Total Pin Joint",
   "Total Equation",
   "Total Reaction",
   
   "Modal",
   "Screw",
   
   NULL
};


/* Tipi di Rod */
const char* psRodNames[] = {
   "Elastic Rod",
   "Visco-Elastic Rod",
   "Visco-Elastic Rod With Offset",
   
   NULL
};


/* Tipi di DeformableHinge */
const char* psConstLawNames[] = {
   "Elastic",
   "Viscous",
   "Visco-Elastic",
   
   NULL
};


/* Tipi di Genel */
const char* psGenelNames[] = {
   "Swash Plate",
   "Rotor Trim",
   "Clamp",
   "Distance",
   "Spring",
   "Spring Support",
   "Cross Spring Support",
   "Spring Damper",
   "Spring Damper Support",
   "Cross Spring Damper Support",
   "Mass",
   "Scalar Filter",
   "State Space SISO",
   "State Space MIMO",
   
   NULL
};


/* Tipi di Force */
const char* psForceNames[] = {
   "Abstract Force",
   "Abstract Reaction Force",

   "Conservative Force",
   "Follower Force",
   "Conservative Couple",
   "Follower Couple",

   "Conservative Reaction Force",
   "Follower Reaction Force",
   "Conservative Reaction Couple",
   "Follower Reaction Couple",

   "External Structural Force",
   
   NULL
};


/* Tipi di elementi elettrici */
const char* psElectricNames[] = {
   "Accelerometer",
   "Discrete Control",
   "Motor",
   
   NULL
};


/* Tipi di elementi idraulici */
const char* psHydraulicNames[] = {
   "Pipe",
   "Minor Losses",
   "Control Valve",
   "Dynamic Control Valve",
   "Pressure Valve",
   "Flow Valve",
   "Orifices",
   "Accumulator",
   "Tank",
   "Full Pipe",
   "Prova",
   "Pipe (New)",
   "Pipe (Tubo)",
   
   NULL
};


/* Tipi di Beam */
const char* psBeamNames[] = {
   "Elastic Beam",
   "Visco-Elastic Beam",
   "Piezo-Electric Beam",
   
   NULL
};


/* Tipi di Aero */
const char* psAeroNames[] = {
   "Rotor",
   "Aeromodal",
   "Aerodynamic Body",
   "Aerodynamic Beam",
   "Aerodynamic External",
   "Aerodynamic External Modal",
   
   NULL
};


/* Tipi di Rotor */
const char* psRotorNames[] = {
   "No Induced Speed Rotor",
   "Uniform Induced Speed Rotor",
   "Glauert Induced Speed Rotor",
   "Mangler Induced Speed Rotor",
   "Dynamic Inflow Induced Speed Rotor",
   
   NULL
};


/* Tipi di Bulk */
const char* psBulkNames[] = {
   "SpringSupport",
   "Spring",
   
   NULL
};


/* Tipi di Drive */
const char* psDriveNames[] = {
   "File",
   
   NULL
};

const char* psReadControlDrivers[] = {
   "file" "drivers",
   
   NULL
};

/* Tipi di FileDrive */
const char* psFileDriveNames[] = {
   "FixedStepFile",
   "Socket",
   
   NULL
};

/* Tipi di Node */
const char* psNodeNames[] = {   
   "Abstract",     
   "Structural",   
   "Electric",   
   "Thermal",   
   "Parameter",   
   "Hydraulic",
   
   NULL
};


/* Nomi dei nodi quando vengono dichiarati nel blocco dati di controllo;
 * hanno lo stesso numero progressivo dei nodi, quindi vi si accede con 
 * l'enum Node::Type */
const char* psReadControlNodes[] = {
   "abstract" "nodes",
   "structural" "nodes",   
   "electric" "nodes",   
   "thermal" "nodes",   
   "parameter" "nodes",   
   "hydraulic" "nodes",
   
   NULL
};

const char* psReadNodesNodes[] = {
   "abstract",
   "structural",   
   "electric",   
   "thermal",   
   "parameter",   
   "hydraulic",
   
   NULL
};


/* Tipi di StructNode */
const char* psStructNodeNames[] = {
   "Dynamic",
   "Static",
   "Modal",
   "Dummy",
   
   NULL
};


/* Tipi di DofOwner */

const char* psDofOwnerNames[] = {
   "Structural Node",
   "Electric Node",
   "Thermal Node",
   "Abstract Node",
   "Hydraulic Node",
   
   "Joint Element",
   "Plate Element",
   "GENEL Element",
   "Induced Velocity Element",
   "Aerodynamic Element",
   "Aerodynamic Modal Element",
   "Electric Bulk Element",
   "Electric Element",
   "Thermal Element",
   "Hydraulic Element",
   "Loadable Element",
   
   NULL
};   

