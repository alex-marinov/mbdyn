#moduleFMU
This is a module to import and cosimulation of FMUs, developed using FMI standards, for MBDyn.

#Requirement
1. MBDyn simulation software. (can be downloaded from https://www.mbdyn.org/?Software_Download)
2. This module is created by using FMILibrary. FMILibrary is a JModellica software  that enables integration of Functional Mock-up Units (FMUs) import in applications. It can be downloaded from http://www.jmodelica.org/FMILibrary
3. In Makefile.inc from modules/module-FMU/ in mbdyn-VERSION, replace   
	$FMIL: Location where build directory of FMILibrary exists.  
	$FMII: Location where FMILibrary is installed.  

#Input Syntax 
The syntax for using the module is:

module load: "libmodule-FMU";

user defined:  <br />
\<label\>, FMU, “\<location to FMU\>”, <br />
type, \<cosimulation/import\>, tolerance, \<tolerance value\>, <br />
\<fmu input variable\>, \<drive caller\>, <br />
\<fmu input variable\>, \<drive caller\>, <br />
..... <br />
output, yes; <br />

*Example Usage:* <br />
user defined:  <br />
99, FMU, "/location/to/VanderPol.fmu", <br />
type, cosimulation, <br />
"u", node, 1, structural, string, "X[3]", direct, <br />
"u2", node, 1, structural, string, "X[2]", direct, <br />
output, yes; <br />

#Testing
To make sure module is working fine,
1. Look for FMIL Error
2. MBDyn should terminate successfully
3. The .out file will have the output states.

#Limitations
1. Currently the module does not support the string drive.
2. The module gives a segmentation fault if FMU with version 1.0 is used with simulation type not supported by the FMU.
3. Constant timestep for cosimulation. 

#Test Cases:
All the FMUs which are tested are from:
https://trac.fmi-standard.org/browser/branches/public/Test_FMUs
