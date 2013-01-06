/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2013
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

/*****************************************************************************
 *                                                                           *
 *                               HARWLIB                                     *
 *                                                                           *
 *****************************************************************************/

/* Dichiarazione delle routines include nel file <harwlib.c>, 
 * che sono un subset della libreria di calcolo HARWELL. 
 * In particolare sono pate delle routines di gestione di matrici sparse,
 * con assemblaggio sparso, fattorizzazione LU, soluzione ecc.
 * Convertite in C/C++ da Pierangelo Masarati con f2c */

/* Per il significato dei vari termini, consultare il file <harwlib.c>, 
 * in cui c'e' un echo dei commenti originali in FORTRAN */


#ifndef HARWLIB_H
#define HARWLIB_H

#include <ac/f2c.h>

#ifdef USE_UNDERSCORES
#define __U_S__ _
#else
#define __U_S__
#endif

#ifdef __cplusplus
extern "C" {
#endif
   
extern int __FC_DECL__(ma28ad)(integer *n, integer *nz, doublereal *a,
			       integer *licn, integer *irn, integer *lirn,
			       integer *icn, doublereal *u, integer *ikeep, 
			       integer *iw, doublereal *w, integer *iflag);
extern int __FC_DECL__(ma28bd)(integer *n, integer *nz, doublereal *a,
			       integer *licn, integer *ivect, integer *jvect,
			       integer *icn, integer *ikeep, integer *iw, 
			       doublereal *w, integer *iflag);
extern int __FC_DECL__(ma28cd)(integer *n, doublereal *a, integer *licn,
			       integer *icn, integer *ikeep, doublereal *rhs, 
			       doublereal *w, integer *mtype);
extern int __FC_DECL__(ma28dd)(integer *n, doublereal *a, integer *licn,
			       integer *ivect, integer *jvect, integer *nz,
			       integer *icn, integer *lenr, integer *lenrl, 
			       integer *lenoff, integer *ip, integer *iq,
			       integer *iw1, integer *iw, doublereal *w1,
			       integer *iflag);
extern int __FC_DECL__(ma28id)(integer *n, integer *nz, doublereal *aorg,
			       integer *irnorg, integer *icnorg, 
			       integer *licn, doublereal *a, integer *icn, 
			       integer *ikeep, doublereal *rhs, 
			       doublereal *x, doublereal *r, doublereal *w,
			       integer *mtype, doublereal *prec, 
			       integer *iflag);
extern int __FC_DECL__(ma30ad)(integer *nn, integer *icn, doublereal *a,
			       integer *licn, integer *lenr, integer *lenrl,
			       integer *idisp, integer *ip, integer *iq, 
			       integer *irn, integer *lirn, integer *lenc,
			       integer *ifirst, integer *lastr,
			       integer *nextr, integer *lastc, integer *nextc,
			       integer *iptr, integer *ipc, doublereal *u,
			       integer *iflag);
extern int __FC_DECL__(ma30bd)(integer *n, integer *icn, doublereal *a,
			       integer *licn, integer *lenr, integer *lenrl,
			       integer *idisp, integer *ip, integer *iq, 
			       doublereal *w, integer *iw, integer *iflag);
extern int __FC_DECL__(ma30cd)(integer *n, integer *icn, doublereal *a,
			       integer *licn, integer *lenr, integer *lenrl,
			       integer *lenoff, integer *idisp, integer *ip,
			       integer *iq, doublereal *x, doublereal *w,
			       integer *mtype);
extern int __FC_DECL__(ma30dd)(doublereal *a, integer *icn, integer *iptr,
			       integer *n, integer *iactiv, integer *itop, 
			       logical *reals);
extern int __FC_DECL__(mc13d)(integer *n, integer *icn, integer *licn,
			      integer *ip, integer *lenr, integer *ior, 
			      integer *ib, integer *num, integer *iw);
extern int __FC_DECL__(mc13e)(integer *n, integer *icn, integer *licn,
			      integer *ip, integer *lenr, integer *arp,
			      integer *ib, integer *num, integer *lowl, 
			      integer *numb, integer *prev);
extern int __FC_DECL__(mc19ad)(integer *n, integer *na, doublereal *a,
			       integer *irn, integer *icn, real *r, real *c,
			       real *w);
extern int __FC_DECL__(mc20ad)(integer *nc, integer *maxa, doublereal *a,
			       integer *inum, integer *jptr, integer *jnum,
			       integer *jdisp);
extern int __FC_DECL__(mc20bd)(integer *nc, integer *maxa, doublereal *a,
			       integer *inum, integer *jptr);
extern int __FC_DECL__(mc21a)(integer *n, integer *icn, integer *licn,
			      integer *ip, integer *lenr, integer *iperm,
			      integer *numnz, integer *iw);
extern int __FC_DECL__(mc21b)(integer *n, integer *icn, integer *licn,
			      integer *ip, integer *lenr, integer *iperm,
			      integer *numnz, integer *pr, integer *arp,
			      integer *cv, integer *out);
extern int __FC_DECL__(mc22ad)(integer *n, integer *icn, doublereal *a,
			       integer *nz, integer *lenrow, integer *ip,
			       integer *iq, integer *iw, integer *iw1);
extern int __FC_DECL__(mc23ad)(integer *n, integer *icn, doublereal *a,
			       integer *licn, integer *lenr, integer *idisp,
			       integer *ip, integer *iq, integer *lenoff,
			       integer *iw, integer *iw1);
extern int __FC_DECL__(mc24ad)(integer *n, integer *icn, doublereal *a,
			       integer *licn, integer *lenr, integer *lenrl, 
			       doublereal *w);
/* comlen ma28ed_ 16 */
/* comlen ma28fd_ 52 */
/* comlen ma28gd_ 8 */
/* comlen ma28hd_ 80 */
/* comlen ma30ed_ 16 */
/* comlen ma30fd_ 20 */
/* comlen ma30id_ 28 */
/* comlen mc23bd_ 20 */
/* comlen ma30gd_ 16 */
/* comlen ma30hd_ 8 */
/* comlen mc19bd_ 8 */

   
/* Common Block Declarations */
   
struct ma28ed_1_ {
    integer lp, mp;
    logical lblock, grow;
};

#define ma28ed_1 (*(struct ma28ed_1_ *) &ma28ed_)

struct ma28fd_1_ {
    doublereal eps, rmin, resid;
    integer irncp, icncp, minirn, minicn, irank;
    logical abort1, abort2;
};

#define ma28fd_1 (*(struct ma28fd_1_ *) &ma28fd_)
      
struct ma28hd_1_ {
    doublereal tol, themax, big, dxmax, errmax, dres, cgce;
    integer ndrop, maxit, noiter, nsrch, istart;
    logical lbig;
};

#define ma28hd_1 (*(struct ma28hd_1_ *) &ma28hd_)   
   
extern struct ext_ma28ed_ {
    integer e_1[2];
    logical e_2[2];
    } ma28ed_;

extern struct ext_ma28fd_ {
    doublereal e_1;
    integer fill_2[9];
    logical e_3[2];
    } ma28fd_;

extern struct ext_ma28hd_ {
    doublereal e_1;
    doublereal fill_2[5];
    doublereal e_3;
    integer fill_4[1];
    integer e_5;
    integer fill_6[1];
    integer e_7[2];
    logical e_8;
    } ma28hd_;
   
extern struct ext_ma28gd_ {   
   integer idisp[2];
} ma28gd_;
   
extern union ext_ma30fd_ {
    struct {
	integer mirncp, micncp, mirank, mirn, micn;
    } _1;
    struct {
	integer irncp, icncp, irank, minirn, minicn;
    } _2;
} ma30fd_;

   
   
#ifdef __cplusplus
}
#endif
#endif
