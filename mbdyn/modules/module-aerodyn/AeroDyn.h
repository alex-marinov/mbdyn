/* 
 * Copyright (C) 2003
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
 * This header file is free software; you can redistribute it at will,
 * under the same license conditions of the AeroDyn package.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 */

#ifndef AERO_DYN_H
#define AERO_DYN_H

/* 
 * compile f90 files with

	ifc -r8

 * for double precision; defaults to single precision
 * link C++ executable with

	g++ -L /opt/intel/ia32/lib/ -lF90 -lCEPCF90 -lintrins

 */

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#if defined(USE_SINGLE_PRECISION)
typedef float F_REAL;
#elif defined(USE_DOUBLE_PRECISION)
typedef double F_REAL;
#else /* !USE_SINGLE_PRECISION && !USE_DOUBLE_PRECISION */
#error "define either USE_SINGLE_PRECISION or USE_DOUBLE_PRECISION"
#endif /* !USE_SINGLE_PRECISION && !USE_DOUBLE_PRECISION */
typedef long int F_LOGICAL;
typedef char F_CHAR;
typedef long int F_INTEGER;

/*
 * Info from:

			USER'S GUIDE
	to the Wind Turbine Aerodynamics Computer Software
			AeroDyn
			
	David J. Laino
	A. Craig Hansen
	Windward Engineering, LC
	Salt Lake City, UT 84117
	www.windwardengineering.com

	Phone:	801-278-7852
	Fax:	801-272-4132
	email:	dlaino@windwardengineering.com
		chansen@windwardengineering.com
	
	Software date and version
	AeroDyn 12.43, 26-Apr-2002
	
	Prepared for the
	National Renewable Energy Laboratory
	under Subcontract No. TCX-9-29209-01

*/
 
/*
 * AeroDyn initialization; must be called as early as possible.
 */
int __FC_DECL__(adinputgate)(void);

/*
 * Returns the force and moment for a given element.
 */
int __FC_DECL__(aerofrcintrface)(F_LOGICAL *FirstLoop, F_INTEGER *JElem,
		F_REAL *DFN, F_REAL *DFT, F_REAL *PMA);

/*
 * Rotor parameters - called once per time step.
 */
int __FC_DECL__(getrotorparams)(F_REAL *Omega, F_REAL *gamma, F_REAL *VHUB,
		F_REAL *tau);

/*
 * Blade parameters - called once for each blade at each time step.
 */
int __FC_DECL__(getbladeparams)(F_REAL *psi);

/*
 * Element parameters - called once for each element at each time step.
 */
int __FC_DECL__(getelemparams)(F_INTEGER *MulTabLoc, F_REAL *phi,
		F_REAL *radius,
		F_REAL *XGRND, F_REAL *YGRND, F_REAL *ZGRND);
/*
 * Compute VT, VN{W,E} based on VX, VY, VZ of the wind.
 */
int __FC_DECL__(getvnvt)(F_REAL *VX, F_REAL *VY, F_REAL *VZ,
		F_REAL *VT, F_REAL *VNW, F_REAL *VNE);

/*
 * Write an error message to the appropriate stream
 *
 * FIXME: the "msg" and "level" arrays should be reset by the caller
 * before writing the message, otherwise they're not '\0' terminated
 */
int __FC_DECL__(usrmes)(F_LOGICAL *Logical, F_CHAR msg[],
		F_INTEGER *code, F_CHAR level[]);

/*
 * self explanatory :)
 */
int __FC_DECL__(abort)(void);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AERO_DYN_H */

