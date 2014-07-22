/* $Header$ */
/* 
 * Copyright (C) 2003-2014
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
 * NOTE to gfortran users:
 *
 * - get AeroDyn 12.58
 * - apply patch aerodyn-12.58-mbdyn.patch
 * - run
	gfortran -O -c *.f90
	ar ru libAeroDyn.a *.o
 * - place libAeroDyn.a where the linker can find it,
 *   or tweak Makefile.inc as appropriate
 */

/*
 * NOTE to icc users:
 *
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
extern int
__FC_DECL__(ad_inputgate)(F_CHAR *input_file);
extern int
__FC_DECL__(adinputgate)(void);

// ADDED BY JENS VAN SCHELVE TO PROVIDE AERODYN ELEMENT DATA OUTPUT
extern int
__FC_DECL__(elemout)(void);

/*
 * Returns the force and moment for a given element.
 */
extern int
__FC_DECL__(aerofrcintrface)(F_LOGICAL *FirstLoop, F_INTEGER *JElem,
		F_REAL *DFN, F_REAL *DFT, F_REAL *PMA);

/*
 * Rotor parameters - called once per time step.
 */
extern int
__FC_DECL__(getrotorparams)(F_REAL *Omega, F_REAL *gamma, F_REAL *VHUB,
		F_REAL *tau);

/*
 * Blade parameters - called once for each blade at each time step.
 */
extern int
__FC_DECL__(getbladeparams)(F_REAL *psi);

/*
 * Element parameters - called once for each element at each time step.
 */
extern int
__FC_DECL__(getelemparams)(F_INTEGER *MulTabLoc, F_REAL *phi,
		F_REAL *radius,
		F_REAL *XGRND, F_REAL *YGRND, F_REAL *ZGRND);

/*
 * Compute VT, VN{W,E} based on VX, VY, VZ of the wind.
 */
extern int
__FC_DECL__(getvnvt)(F_REAL *VX, F_REAL *VY, F_REAL *VZ,
		F_REAL *VT, F_REAL *VNW, F_REAL *VNE);

/*
 * Write an error message to the appropriate stream
 *
 * FIXME: the "msg" and "level" arrays should be reset by the caller
 * before writing the message, otherwise they're not '\0' terminated
 */
extern int
__FC_DECL__(usrmes)(F_LOGICAL *Logical, F_CHAR msg[],
		F_INTEGER *code, F_CHAR level[]);

/*
 * self explanatory :)
 */
extern int
__FC_DECL__(ad_abort)(void);

/*
 * MBDyn stuff initialization
 */
extern int
__FC_DECL__(mbdyn_init)(F_CHAR *Version, F_INTEGER *nblades, F_REAL *rotradius);

/*
 * AeroDyn initialization wrapper
 */
extern int
__FC_DECL__(mbdyn_ad_inputgate)(F_CHAR *ifname, F_INTEGER *ifnamelen, F_CHAR *efname, F_INTEGER *efnamelen);

/*
 * MBDyn logical true
 */
extern int
__FC_DECL__(mbdyn_true)(F_LOGICAL *val);

/*
 * MBDyn logical false
 */
extern int
__FC_DECL__(mbdyn_false)(F_LOGICAL *val);

/*
 * By Fanzhong Meng 21 Feb. 2008
 * MBDyn-AeroDyn common module defined data
 */
extern int
__FC_DECL__(mbdyn_com_data)(F_INTEGER *c_blade, F_INTEGER *c_elem);

/* 
 * This subroutine is to pass the current simulation time 
 * from MBDyn to AeroDyn!
 * c_time: current time
 * By Fanzhong MENG 19 June 2008
 */
extern int
__FC_DECL__(mbdyn_sim_time)(doublereal *c_time);

/* 
 * This subroutine is to pass the current simulation time step 
 * from MBDyn to AeroDyn!
 * dt: time step
 * By Fanzhong MENG 19 June 2008
 */
extern int
__FC_DECL__(mbdyn_time_step)(F_REAL *dt);


/* 
 * This Subroutine is to get the Tip loss constants.
 * RLOCAL: store the position of the blade element.
 * Cur_elem: store the number of current blade element.
 * By Fanzhong Meng 14 August 2008
 */
extern int
__FC_DECL__(mbdyn_get_tl_const)(F_REAL *RLOCAL, F_INTEGER *Cur_elem);


/* 
 * This Subroutine is to get the Hub loss constants.
 * RLOCAL: store the position of the blade element.
 * Cur_elem: store the number of current blade element.
 * RHub:     Store the hub radius.
 * By Fanzhong Meng 14 August 2008
 */
extern int
__FC_DECL__(mbdyn_get_hl_const)(F_REAL *RLOCAL, F_INTEGER *Cur_elem, F_REAL *RHub);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* AERO_DYN_H */

