/* wheel.f -- translated by f2c (version 19960717).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#ifdef __cplusplus
extern "C" {
#endif
#include "f2c.h"

/* Common Block Declarations */

struct {
    doublereal dr;
} pneu_;

#define pneu_1 pneu_

struct {
    doublereal alfa, vass1, vass2, ccy, ssy, sr, pippo;
} velang_;

#define velang_1 velang_

/* Table of constant values */

static integer c__1 = 1;
static integer c__3 = 3;
static doublereal c_b59 = 2.;
static real c_b76 = (float)0.;
static integer c__9 = 9;
static integer c__5 = 5;
static integer c_n1 = -1;
static doublereal c_b93 = 1.;
static integer c__6 = 6;

/* Subroutine */ int gfosub_(integer *id, doublereal *time, doublereal *par, 
	integer *npar, logical *dflag, logical *iflag, doublereal *result)
{
    /* Builtin functions */
    integer i_dnnt(doublereal *);

    /* Local variables */
    static doublereal beta, omep;
    static integer ircw, npft;
    static doublereal bleng, betap;
    static integer idvar[3];
    static doublereal array[11];
    static integer istat;
    extern /* Subroutine */ int wvefr_(doublereal *, doublereal *, doublereal 
	    *, doublereal *, doublereal *, integer *, integer *, doublereal *,
	     doublereal *, integer *, doublereal *, doublereal *, doublereal *
	    );
    static doublereal cc[9], dx, dy, dz, rm, tabfat[6], vw[6], vx, vy, vz, wx,
	     wy, wz;
    static logical errflg;
    extern /* Subroutine */ int gtaray_(integer *, doublereal *, integer *, 
	    integer *);
    static integer number;
    extern /* Subroutine */ int errmes_(logical *, char *, integer *, char *, 
	    ftnlen, ftnlen), sysfnc_(char *, integer *, integer *, doublereal 
	    *, logical *, ftnlen);
    static integer ipc;
    static doublereal ome;
    static integer idv[3], idw[3];
    static doublereal res[6];
    static integer idarray, idv2;

/*   -------------------------Main------------------------ */
/*   Inizializza variabili */
    /* Parameter adjustments */
    --result;
    --par;

    /* Function Body */
    idvar[0] = i_dnnt(&par[1]);
    idvar[1] = i_dnnt(&par[2]);
    idvar[2] = i_dnnt(&par[3]);
    idarray = i_dnnt(&par[4]);
    idv[0] = i_dnnt(&par[6]);
    idv[1] = i_dnnt(&par[5]);
    idv[2] = i_dnnt(&par[5]);
    idv2 = i_dnnt(&par[6]);
    idw[0] = i_dnnt(&par[7]);
    idw[1] = i_dnnt(&par[6]);
    idw[2] = i_dnnt(&par[6]);
/*   Importa coseni direttori asse mozzo */
    sysfnc_("VARVAL", idvar, &c__1, cc, &errflg, 6L);
    errmes_(&errflg, "Error getting cos(x) in GFOSUB.", id, "STOP", 31L, 4L);
    sysfnc_("VARVAL", &idvar[1], &c__1, &cc[1], &errflg, 6L);
    errmes_(&errflg, "Error getting cos(y) in GFOSUB.", id, "STOP", 31L, 4L);
    sysfnc_("VARVAL", &idvar[2], &c__1, &cc[2], &errflg, 6L);
    errmes_(&errflg, "Error getting cos(z) in GFOSUB.", id, "STOP", 31L, 4L);
/*  Importa posizione mozzo relativamente al punto fisso 0,0,0 sul terreno
*/
    sysfnc_("DX", idv, &c__3, &dx, &errflg, 2L);
    errmes_(&errflg, "Error getting DX in GFOSUB.", id, "STOP", 27L, 4L);
    cc[3] = dx;
    sysfnc_("DY", idv, &c__3, &dy, &errflg, 2L);
    errmes_(&errflg, "Error getting DY in GFOSUB.", id, "STOP", 27L, 4L);
    cc[4] = dy;
    sysfnc_("DZ", idv, &c__3, &dz, &errflg, 2L);
    errmes_(&errflg, "Error getting DZ in GFOSUB.", id, "STOP", 27L, 4L);
    cc[5] = dz;
/*   Inizializza costanti */
    gtaray_(&idarray, array, &number, &istat);
    cc[6] = array[0];
    cc[7] = array[1];
    cc[8] = array[2];
    tabfat[0] = array[3];
    tabfat[1] = array[4];
    tabfat[2] = array[5];
    tabfat[3] = array[6];
    tabfat[4] = array[7];
    tabfat[5] = array[8];
    ipc = i_dnnt(&array[9]);
    npft = i_dnnt(&array[10]);
/*   Importa velocita lineari e angolari del mozzo */
    sysfnc_("VX", &idv2, &c__1, &vx, &errflg, 2L);
    errmes_(&errflg, "Error getting VX in GFOSUB.", id, "STOP", 27L, 4L);
    vw[0] = vx;
    sysfnc_("VY", &idv2, &c__1, &vy, &errflg, 2L);
    errmes_(&errflg, "Error getting VY in GFOSUB.", id, "STOP", 27L, 4L);
    vw[1] = vy;
    sysfnc_("VZ", &idv2, &c__1, &vz, &errflg, 2L);
    errmes_(&errflg, "Error getting VZ in GFOSUB.", id, "STOP", 27L, 4L);
    vw[2] = vz;
    sysfnc_("WX", &idv2, &c__1, &wx, &errflg, 2L);
    errmes_(&errflg, "Error getting WX in GFOSUB.", id, "STOP", 27L, 4L);
    vw[3] = wx;
    sysfnc_("WY", &idv2, &c__1, &wy, &errflg, 2L);
    errmes_(&errflg, "Error getting WY in GFOSUB.", id, "STOP", 27L, 4L);
    vw[4] = wy;
    sysfnc_("WZ", &idv2, &c__1, &wz, &errflg, 2L);
    errmes_(&errflg, "Error getting WZ in GFOSUB.", id, "STOP", 27L, 4L);
    vw[5] = wz;
/*   Importa velocita' angolare ruota in asse mozzo */
    sysfnc_("WZ", idw, &c__3, &ome, &errflg, 2L);
    errmes_(&errflg, "Error getting OME in GFOSUB.", id, "STOP", 28L, 4L);
/*   Importa accelerazione angolare ruota in asse mozzo */
    sysfnc_("WDTZ", idw, &c__3, &omep, &errflg, 4L);
    errmes_(&errflg, "Error getting OMEP in GFOSUB.", id, "STOP", 29L, 4L);
/*   Chiama WVEFR */
    wvefr_(vw, &ome, cc, res, &rm, &ircw, &ipc, &omep, tabfat, &npft, &beta, &
	    betap, &bleng);
/*   Restituisce vettore risultati RES: Fx,Fy,Fz,Tx,Ty,Tz */
    result[1] = res[0];
    result[2] = res[1];
    result[3] = res[2];
    result[4] = res[3];
    result[5] = res[4];
    result[6] = res[5];
    return 0;
} /* gfosub_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-82                      SUBROUTINE WVEFR */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int wvefr_(doublereal *vw, doublereal *ome, doublereal *cc, 
	doublereal *res, doublereal *rm, integer *ircw, integer *ipc, 
	doublereal *omep, doublereal *tabfat, integer *npft, doublereal *beta,
	 doublereal *betap, doublereal *bleng)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), sqrt(doublereal), acos(
	    doublereal);

    /* Local variables */
    static doublereal espr, rmum, rmus;
    static integer i__;
    static doublereal arrif;
    extern /* Subroutine */ int vercp_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), setpl_(doublereal *, doublereal *), tyfrc_(
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *);
    static doublereal press;
    extern /* Subroutine */ int whtor_(doublereal *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, doublereal *);
    static doublereal rmust, ap, rc, gp[3], cracmb, pi, rk, pp[6], rr, wp[12],
	     coevol, smumax, cdp[12]	/* was [3][4] */, xch, zcg[18]	/* 
	    was [6][3] */, alfarif, pen, vel[3], vol, trv, wrv;

/*      common /forces/res1,res2,res3 */
/*   Inizializza variabili */
    /* Parameter adjustments */
    --tabfat;
    --res;
    --cc;
    --vw;

    /* Function Body */
    *rm = (float)0.;
    press = (float)4.65e5;
    rc = cc[9];
    rr = cc[7];
/*      COEVOL = .3 */
/*      COEVOL = 1. */
    coevol = (float).8;
/*  coevol:vol rif ruota? */
    wrv = (float)40.;
    rmust = (float).1;
    xch = (float)0.;
    cracmb = (float)1.;
    smumax = (float).075;
    rmum = (float).7;
    rmus = (float).1;
    espr = (float)1.25;
    trv = (float).005;
    d__1 = rr - rc;
    arrif = sqrt(pow_dd(&rr, &c_b59) - pow_dd(&d__1, &c_b59)) * (float)3.7 * 
	    rc;
    rk = 1.;
    pi = acos((float)-1.);
    alfarif = pi / (float)15.;
    wp[0] = press;
    wp[1] = pi * (float)2. * pi * rc * rc * (rr - rc);
    wp[1] *= coevol;
    wp[2] = wrv;
    wp[3] = (float)0.;
    wp[4] = arrif;
    wp[5] = rmust;
    wp[6] = xch;
    wp[7] = cracmb;
    wp[8] = smumax;
    wp[9] = rmum;
    wp[10] = rmus;
    gp[0] = espr;
    gp[1] = trv;
    gp[2] = (float)0.;
    pp[0] = (float)1e11;
    pp[1] = (float).3;
    pp[2] = (float).5;
    pp[3] = (float)1.;
    pp[4] = (float).8;
    pp[5] = (float)1.;
    for (i__ = 1; i__ <= 6; ++i__) {
	res[i__] = (float)0.;
/* L5: */
    }
    setpl_(&cc[4], cdp);
    whtor_(&cc[1], cdp, ircw, zcg, &vol, &pen, &arrif);
    if (*ircw == 0) {
	return 0;
    }
    vercp_(&cc[4], &vw[1], &zcg[15], vel);
    tyfrc_(vel, ome, zcg, &vol, &ap, cdp, pp, wp, gp, &cc[1], &res[1], rm, 
	    ircw, ipc, omep, &rk, &tabfat[1], npft, beta, betap, bleng, &
	    rmust, &alfarif);
/*      res1 = res(1) */
/*      res2 = res(2) */
/*      res3 = res(3) */
    return 0;
} /* wvefr_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-82                      SUBROUTINE WHTOR */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int whtor_(doublereal *cc, doublereal *cdp, integer *ircw, 
	doublereal *zcg, doublereal *vol, doublereal *dl, doublereal *arrif)
{
    /* System generated locals */
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), sin(doublereal);

    /* Local variables */
    static doublereal alfa, tvet[9];
    extern /* Subroutine */ int prvn_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static doublereal a[3];
    static integer i__, k;
    static doublereal r__[2], x[3], a1, a2, a3, d1, d2;
    extern /* Subroutine */ int rmove_(doublereal *, doublereal *, integer *);
    static doublereal ff, cr[3], drc, drl;

/* ------------------------------ */
    /* Parameter adjustments */
    zcg -= 7;
    cdp -= 4;
    --cc;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	for (k = 1; k <= 6; ++k) {
	    zcg[k + i__ * 6] = (float)0.;
/* L5: */
	}
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	cr[i__ - 1] = cc[i__ + 3] + cc[i__] * (float).5 * (cc[8] + cc[9]);
/* L10: */
    }
    rmove_(&cdp[10], &tvet[6], &c__3);
    prvn_(&tvet[6], &cc[1], tvet, &ff);
    prvn_(&tvet[6], tvet, &tvet[3], &ff);
    d1 = (float)0.;
    d2 = (float)0.;
    a1 = (float)0.;
    a2 = (float)0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	d1 -= tvet[i__ + 5] * cdp[i__ + 12];
	d2 -= cc[i__] * cr[i__ - 1];
	a1 += cc[i__] * tvet[i__ + 2];
	a2 += cc[i__] * tvet[i__ + 5];
/* L15: */
    }
    a3 = (float)0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	x[i__ - 1] = -d1 * tvet[i__ + 5] - (d2 + a2 * d1) * tvet[i__ + 2] / 
		a1;
	a3 += tvet[i__ - 1] * (cr[i__ - 1] - x[i__ - 1]);
/* L20: */
    }
    for (i__ = 1; i__ <= 3; ++i__) {
	zcg[i__ + 9] = x[i__ - 1] + a3 * tvet[i__ - 1];
	for (k = 2; k <= 3; ++k) {
	    zcg[i__ + 3 + k * 6] = zcg[i__ + 9];
/* L25: */
	}
    }
    pneu_1.dr = (float)0.;
    for (i__ = 1; i__ <= 3; ++i__) {
/* Computing 2nd power */
	d__1 = zcg[i__ + 9] - cr[i__ - 1];
	pneu_1.dr += d__1 * d__1;
/* L30: */
    }
    pneu_1.dr = sqrt(pneu_1.dr);
    if (pneu_1.dr < cc[7]) {
	goto L35;
    }
    *ircw = 0;
    *dl = (float)0.;
    return 0;
L35:
    *ircw = 1;
    *dl = cc[7] - pneu_1.dr;
    a[2] = *arrif * *dl * (float)2. / (d__1 = cc[9] - cc[8], abs(d__1));
    r__[0] = (d__1 = cc[8] - cc[9], abs(d__1)) * (float).5;
    r__[1] = cc[7];
    for (i__ = 1; i__ <= 2; ++i__) {
	drl = r__[i__ - 1] - *dl;
	drc = sqrt(r__[i__ - 1] * r__[i__ - 1] - drl * drl);
	alfa = atan2(drc, drl) * (float)2.;
	a[i__ - 1] = r__[i__ - 1] * (float).5 * r__[i__ - 1] * (alfa - sin(
		alfa));
/* L40: */
    }
    for (k = 1; k <= 3; ++k) {
	for (i__ = 1; i__ <= 3; ++i__) {
	    zcg[i__ + k * 6] = a[k - 1] * tvet[(k - 1) * 3 + i__ - 1];
/* L45: */
	}
    }
    *vol = a[2] * *dl * (float).5;
    return 0;
} /* whtor_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   29- 9-83                      SUBROUTINE SETPL */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int setpl_(doublereal *xyz, doublereal *cdp)
{
    /* Builtin functions */
    double cos(doublereal), sin(doublereal);

    /* Local variables */
    static doublereal alfa;
    extern /* Subroutine */ int rset_(doublereal *, integer *, real *), prvu_(
	    doublereal *, doublereal *, doublereal *);
    static integer i__;

/* ------------------------------------- */
    /* Parameter adjustments */
    cdp -= 4;
    --xyz;

    /* Function Body */
    alfa = (float)0.;
    for (i__ = 1; i__ <= 4; ++i__) {
	rset_(&cdp[i__ * 3 + 1], &c__3, &c_b76);
/* L5: */
    }
    cdp[4] = cos(alfa);
    cdp[6] = sin(alfa);
    cdp[10] = sin(alfa);
    cdp[12] = -cos(alfa);
    prvu_(&cdp[10], &cdp[4], &cdp[7]);
    cdp[13] = xyz[1];
    cdp[14] = xyz[2];
    cdp[15] = (float)0.;
    return 0;
} /* setpl_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-82                      SUBROUTINE TYFRC */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int tyfrc_(doublereal *vel, doublereal *ome, doublereal *zcg,
	 doublereal *vol, doublereal *ap, doublereal *cdp, doublereal *pp, 
	doublereal *wp, doublereal *gp, doublereal *cc, doublereal *res, 
	doublereal *rm, integer *ircw, integer *ipc, doublereal *omep, 
	doublereal *rk, doublereal *tabfat, integer *npfat, doublereal *beta, 
	doublereal *betap, doublereal *bleng, doublereal *rmust, doublereal *
	alfarif)
{
    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle();
    double sqrt(doublereal), sin(doublereal), tanh(doublereal), d_sign(
	    doublereal *, doublereal *), acos(doublereal), atan2(doublereal, 
	    doublereal);

    /* Local variables */
    static doublereal area, dfcr, fbet, fric, cosa, sina, omer, rant, 
	    anderval, vass[2];
    extern /* Subroutine */ int prvn_(doublereal *, doublereal *, doublereal *
	    , doublereal *);
    static doublereal runt, a, b[3], f[6];
    static integer i__, ipcad;
    static doublereal omega, ander;
    extern doublereal tange_(doublereal *, doublereal *, integer *);
    extern /* Subroutine */ int tmmnt_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), trsfv_(doublereal *, doublereal *, doublereal *,
	     integer *);
    static doublereal bb, cd, fd, bi, br, rd[3], va, pi, ay, th, fx, fy, rr[3]
	    , rs[3], raggio;
    extern /* Subroutine */ int forces_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *);
    static doublereal omezer, va1, va2, aaa, rif, sig, vep[3], sro;
    extern /* Subroutine */ int mux_(doublereal *, doublereal *, doublereal *,
	     integer *, doublereal *), muy_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *);
    static doublereal vep1, vep2;

    /* Fortran I/O blocks */
    static cilist io___83 = { 0, 6, 0, 0, 0 };
    static cilist io___121 = { 0, 6, 0, 0, 0 };


/*                       ,VREL(2) */
/* ------------------------------------- */
    /* Parameter adjustments */
    --tabfat;
    --res;
    --cc;
    --gp;
    --wp;
    --pp;
    cdp -= 4;
    zcg -= 7;
    --vel;

    /* Function Body */
    omezer = (float)0.;
    ipcad = 0;
    bi = (float)0.;
    velang_1.sr = (float)0.;
/* ---------------------------------------------------- */
/*     ELASTIC,DAMPING,FRICTION AND PLOWING FORCES */
/* ---------------------------------------------------- */
    if (*ipc == 0) {
	goto L5;
    }
    if (ipcad == 0) {
	omezer = *ome;
    }
    ipcad = 1;
L5:
    forces_(&cc[1], &vel[1], &zcg[7], vol, &cdp[4], &pp[1], &wp[1], &wp[8], &
	    gp[1], &res[1], &cc[4], ap, &fbet, ircw);
    for (i__ = 1; i__ <= 6; ++i__) {
	f[i__ - 1] = (float)0.;
/* L10: */
    }
    if (wp[5] <= (float)0. || wp[6] <= (float)0.) {
	s_wsle(&io___83);
	do_lio(&c__9, &c__1, " FORZE TANGENTI RUOTA NON CALCOLATE !!", 38L);
	do_lio(&c__5, &c__1, (char *)&wp[5], (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&wp[6], (ftnlen)sizeof(doublereal));
	e_wsle();
	return 0;
    }
/* ----------------------------------------- */
/*     VELOCITY IN EQUIVALENT PLANE FRAME */
/* ----------------------------------------- */
    trsfv_(&cdp[4], &vel[1], vep, &c_n1);
/* Computing 2nd power */
    d__1 = vep[0];
/* Computing 2nd power */
    d__2 = vep[1];
    va = sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
    d__1 = zcg[19];
/* Computing 2nd power */
    d__2 = zcg[20];
/* Computing 2nd power */
    d__3 = zcg[21];
    area = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    for (i__ = 1; i__ <= 3; ++i__) {
/* L15: */
	b[i__ - 1] = zcg[i__ + 21] - cc[i__ + 3];
    }
    prvn_(b, &cc[1], rr, &a);
/*      DO 16 I=1,2 */
/*        VREL(I) = OME*A*RR(I)+VEL(I) */
/*        VREL(I) = OME*A*RR(I)+A*VEL(I)/.252 */
/*  16  CONTINUE */
/* ----------------------- */
/*     FRICTION FORCES */
/* ----------------------- */
    if (a == (float)0.) {
	goto L20;
    }
    th = (float)1.;
    if (*ipc == 1 && omezer == (float)0.) {
/* ------------------------------------------------ */
/*     PROVA DI CADUTA SENZA PRE-SPIN DELLA RUOTA */
/* ------------------------------------------------ */
	raggio = (float)0.;
	for (i__ = 1; i__ <= 3; ++i__) {
/* Computing 2nd power */
	    d__1 = zcg[i__ + 21] - cc[i__ + 3];
	    raggio += d__1 * d__1;
/* L25: */
	}
	raggio = sqrt(raggio);
	if ((d__1 = *beta * *betap, abs(d__1)) > (float)1e-7) {
	    omer = *bleng / raggio * *betap * sin(*beta);
	    aaa = *ome / omer;
	} else {
	    aaa = 1.;
	}
	th = tanh((d__1 = (*ome - omer) / (float)70., abs(d__1)));
	velang_1.sr = (float)1. - aaa;
	sig = (float)-1.;
	if (*ome > omer) {
	    sig = (float)1.;
	}
    } else if (*ipc != 1 || omezer != (float)0.) {
/* -------------------------- */
/*     TUTTI GLI ALTRI CASI */
/* -------------------------- */
	omega = (float)0.;
	if (velang_1.sr == (float)0. || *ome == omezer) {
	    goto L30;
	}
	bi = tange_(&velang_1.sr, &tabfat[1], npfat);
	if (*ome != (float)0.) {
	    if (*ome - omezer != (float)0.) {
		omega = *ome - *omep * area * *ap * bi / (*rk * (*ome - 
			omezer));
	    } else {
		omega = (float)0.;
	    }
	} else {
	    omega = *ome;
	}
L30:
	omer = (float)0.;
	for (i__ = 1; i__ <= 3; ++i__) {
	    omer += vel[i__] * rr[i__ - 1];
/* L35: */
	}
/*        OMER = OMER/A */
	velang_1.sr = (omer + *ome * a) / (omer - omega * (float).252);
	sig = d_sign(&c_b93, &velang_1.sr);
	th = (float)1.;
    }
    sro = velang_1.sr;
    velang_1.sr = abs(velang_1.sr);
    if (velang_1.sr > (float)1.) {
	velang_1.sr = (float)1.;
    }
/*      CALL RONT(SR,TABFAT,NPFAT,RANT) */
/* alcolo cosangolo tra piano ruota e asse x del sistema di riferimento */
/* C   terreno */
    velang_1.ccy = cc[2] / sin(acos(cc[3]));
/*           print *,cc(1),cc(2),cc(3) */
/* alcolo sinangolo tra piano ruota e asse x del sistema di riferimento */
/* C   terreno */
    velang_1.ssy = cc[1] / sin(acos(cc[3]));
/* C   Attenzione vale solo se asse non Š normale alla pista. in */
/* C   cosbeta (denominatore) non */
/* C   ce il segno perchŠ il suo valore Š positivo indipententemente */
/* C   dalla sua posizione */
/* alcolo angolo tra piano ruota e asse x del sistema di riferimento */
/* C   terreno */
    ay = atan2(velang_1.ssy, velang_1.ccy);
/* C   Non ci sono problenmi se l'asse mozzo non Š orientato come x */
    vass[0] = -vel[1] * velang_1.ccy + vel[2] * velang_1.ssy;
    vass[1] = -vel[2] * velang_1.ccy + vel[1] * velang_1.ssy;
/*           vassabs = sqrt(vass(1)**2+vass(2)**2) */
/* alcolo angolo di deriva */
    va1 = abs(vass[0]);
    velang_1.vass1 = vass[0];
    va2 = abs(vass[1]);
    velang_1.vass2 = vass[1];
    rif = (float)1e-6;
    if (va1 >= rif) {
	velang_1.alfa = atan2(vass[1], vass[0]);
    } else if (va1 < rif && va2 >= rif) {
	velang_1.alfa = atan2(vass[0], vass[1]);
	velang_1.alfa = d_sign(&c_b93, &vass[1]) * (pi / (float)2. - 
		velang_1.alfa);
    } else if (va1 < rif && va2 < rif) {
	velang_1.alfa = (float)0.;
    }
    mux_(&velang_1.sr, &velang_1.alfa, &tabfat[1], npfat, &rant);
    fric = th * rant;
    fx = -sig * pp[4] * fric * area * *ap;
    f[0] = fx * (-velang_1.ccy);
    f[1] = fx * velang_1.ssy;
    f[2] = (float)0.;
/*      print *,'ccy,ssy',ccy,ssy */
    tmmnt_(&zcg[22], f, &cc[4], &res[1]);
L20:
/* -------------------- */
/*     DRIVING FORCE */
/* -------------------- */
    if (va > gp[2]) {
	if (a > (float)0.) {
	    goto L45;
	}
	dfcr = wp[6] * (float).7 * pp[6];
	cd = dfcr * *ap * area;
	for (i__ = 1; i__ <= 3; ++i__) {
/* L50: */
	    f[i__ - 1] = -cd * (vep[0] * cdp[i__ + 3] + vep[1] * cdp[i__ + 6])
		     / va;
	}
	goto L55;
L45:
	prvn_(&cdp[10], rr, rd, &a);
	prvn_(rd, &cdp[10], rs, &a);
	trsfv_(&cdp[4], rd, rr, &c_n1);
/*      SINA  = (VREL(1)*RR(1)+VREL(2)*RR(2))/SQRT(VREL(1)**2+VREL(2)
**2) */
	vep1 = vep[0];
	vep2 = vep[1];
	sina = (vep[0] * rr[0] + vep[1] * rr[1]) / va;
/* Computing 2nd power */
	d__1 = sina;
	cosa = sqrt((float)1. - d__1 * d__1);
/*      ANDER = ASIN(SINA) */
	ander = atan2(sina, cosa);
	bb = area * (float)1.32 / wp[5] + (float).002;
/* Computing 2nd power */
	d__1 = bb - (float)1.75;
	br = (bb + (float)1.75 - sqrt(d__1 * d__1 + (float).01)) * (float).5;
	br = ((float)2. - br) * br;
	fd = br / (bb + (float)1e-4);
	anderval = ander;
	ander = fd * ander;
	dfcr = wp[6] * pp[6];
/* ------------------------------------------------------ */
/*      CD    = TANH(vep(2)/GP(2))*DFCR*AP*AREA*ABS(TANH(.7*ANDER/DFCR
)) */
/* ------------------------------------------------------ */
	cd = dfcr * *ap * area * tanh(ander * (float).7 / dfcr);
	velang_1.pippo = tanh(velang_1.alfa / *alfarif);
/*      print *,'pvep1,vep2,va,sina,cosa', vep(1),vep(2),va,sina,cosa 
*/
	muy_(&velang_1.sr, &velang_1.alfa, alfarif, rmust, &runt);
	fy = -runt * *ap * area;
	f[0] = fy * (-velang_1.ssy);
	f[1] = fy * (-velang_1.ccy);
	f[2] = (float)0.;
/*      print *,'Fy',fy,alfa,vass(2) */
	a = cc[7] * (float)2. * wp[7] * tanh(abs(ander) * (float).7) * cosa;
	for (i__ = 1; i__ <= 3; ++i__) {
	    b[i__ - 1] = (zcg[i__ + 21] - rs[i__ - 1] * a * sig) * fbet;
/* L65: */
	}
L55:
	tmmnt_(b, f, &cc[4], &res[1]);
    } else {
	s_wsle(&io___121);
	do_lio(&c__9, &c__1, " Calcoli non effettuati", 23L);
	do_lio(&c__5, &c__1, (char *)&va, (ftnlen)sizeof(doublereal));
	do_lio(&c__5, &c__1, (char *)&gp[2], (ftnlen)sizeof(doublereal));
	e_wsle();
    }
/* ---------------------------------- */
/*     CALCOLO MOMENTO SULLA RUOTA */
/* ---------------------------------- */
    for (i__ = 1; i__ <= 3; ++i__) {
	*rm += res[i__ + 3] * cc[i__];
/* L70: */
    }
    return 0;
} /* tyfrc_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-82                      SUBROUTINE FORCES */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int forces_(doublereal *cc, doublereal *vel, doublereal *zcg,
	 doublereal *vol, doublereal *cdp, doublereal *pp, doublereal *op, 
	doublereal *red, doublereal *gp, doublereal *res, doublereal *rp, 
	doublereal *ap, doublereal *fbet, integer *ircw)
{
    /* Initialized data */

    static doublereal pgrq = .785398;

    /* System generated locals */
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal), atan2(doublereal, doublereal), cos(doublereal), 
	    pow_dd(doublereal *, doublereal *), tanh(doublereal);

    /* Local variables */
    extern /* Subroutine */ int rset_(doublereal *, integer *, real *);
    static doublereal f[6];
    static integer i__, j, k;
    static doublereal x, a1;
    extern /* Subroutine */ int tmmnt_(doublereal *, doublereal *, doublereal 
	    *, doublereal *), trsfv_(doublereal *, doublereal *, doublereal *,
	     integer *);
    static doublereal cb, dc, fc, an, pc, sb, va, vm, gm1, bet;
    static integer igo;
    static doublereal rap, vep[3], opr;

    /* Parameter adjustments */
    --res;
    --gp;
    --op;
    --pp;
    cdp -= 4;
    zcg -= 7;
    --vel;
    --cc;

    /* Function Body */
/* -------------------------------------- */
    sb = (float)0.;
    for (i__ = 1; i__ <= 3; ++i__) {
	sb += cc[i__] * cdp[i__ + 9];
/* L5: */
    }
    sb = abs(sb);
/* Computing 2nd power */
    d__1 = sb;
    cb = sqrt((float)1. - d__1 * d__1);
    bet = atan2(sb, cb);
    *fbet = (float)0.;
    if (bet < pgrq) {
	*fbet = cos(bet * (float)2.);
    }
    opr = op[1] * (*red + ((float)1. - *red) * *fbet);
/* ------------------------ */
/*     ACTUAL PRESSURE */
/* ------------------------ */
    gm1 = (float)1. / gp[1];
    d__1 = opr / pp[1];
    rap = pow_dd(&d__1, &gm1);
    rap *= op[2] / pp[2];
    x = (*vol + rap * pp[2] - op[2]) / (rap + (float)1.);
    if (x > (float)0. && x < *vol) {
	goto L30;
    }
    if (x > (float)0.) {
	goto L20;
    }
/* L10: */
    igo = 2;
    *ircw = 1;
    x = (float)0.;
    d__1 = op[2] / (op[2] - *vol);
    *ap = opr * pow_dd(&d__1, &gp[1]);
    goto L40;
L20:
    igo = 1;
    *ircw = 3;
    x = *vol;
    d__1 = pp[2] / (pp[2] - *vol);
    *ap = pp[1] * pow_dd(&d__1, &gp[1]);
    goto L40;
L30:
    igo = 1;
    *ircw = 2;
    d__1 = op[2] / (op[2] - *vol + x);
    *ap = opr * pow_dd(&d__1, &gp[1]);
L40:
/* --------------------------------- */
/*     ZERO-ING RESULTANT FORCE */
/* --------------------------------- */
    rset_(&res[1], &c__6, &c_b76);
    rset_(f, &c__6, &c_b76);
/* ---------------------- */
/*     ELASTIC FORCE */
/* ---------------------- */
    for (i__ = 1; i__ <= 3; ++i__) {
	f[i__ - 1] = zcg[i__ + 18] * *ap;
/* L45: */
    }
    tmmnt_(&zcg[22], f, rp, &res[1]);
/* ---------------------------------------------------- */
/*     VELOCITY IN EQUIVALENT PLANE FRAME */
/* ---------------------------------------------------- */
    trsfv_(&cdp[4], &vel[1], vep, &c_n1);
/* Computing 2nd power */
    d__1 = vep[0];
/* Computing 2nd power */
    d__2 = vep[1];
    va = sqrt(d__1 * d__1 + d__2 * d__2);
/* Computing 2nd power */
    d__1 = vep[0];
/* Computing 2nd power */
    d__2 = vep[1];
/* Computing 2nd power */
    d__3 = vep[2];
    vm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
/* --------------------------------- */
/*             DAMPING FORCE */
/* --------------------------------- */

    dc = tanh(vep[2] / op[3]);
    for (i__ = 1; i__ <= 3; ++i__) {
	f[i__ - 1] = -dc * zcg[i__ + 18] * *ap;
/* L50: */
    }
    tmmnt_(&zcg[22], f, rp, &res[1]);
/* ------------------------- */
/*     FRICTION FORCE */
/* ------------------------- */
    if (va < gp[2]) {
	goto L55;
    }
/* Computing 2nd power */
    d__1 = zcg[19];
/* Computing 2nd power */
    d__2 = zcg[20];
/* Computing 2nd power */
    d__3 = zcg[21];
    an = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    fc = op[4] * pp[4] * an * *ap / va;
    for (i__ = 1; i__ <= 3; ++i__) {
	f[i__ - 1] = -fc * (vep[0] * cdp[i__ + 3] + vep[1] * cdp[i__ + 6]);
/* L60: */
    }
    tmmnt_(&zcg[22], f, rp, &res[1]);
L55:
/* --------------------- */
/*     PLOWING FORCE */
/* --------------------- */
    if (igo == 2) {
	goto L65;
    }
    if (va < gp[2]) {
	goto L65;
    }
    if (*vol <= (float)0.) {
	goto L65;
    }
    pc = x * (float)2. * *ap * pp[5] * tanh(va / pp[3]) / (va * *vol);
    for (j = 1; j <= 2; ++j) {
	a1 = (float)0.;
	for (k = 1; k <= 3; ++k) {
	    a1 += zcg[k + j * 6] * vel[k];
/* L75: */
	}
	for (i__ = 1; i__ <= 3; ++i__) {
	    f[i__ - 1] = -pc * a1 * cdp[i__ + j * 3];
/* L80: */
	}
	tmmnt_(&zcg[j * 6 + 4], f, rp, &res[1]);
/* L70: */
    }
L65:
    return 0;
} /* forces_ */

/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
/*                 30- 1-82                         FUNCTION TANGE */
/* ------------------------------------------------------------------- */
/* ------------------------------------------------------------------- */
doublereal tange_(doublereal *x, doublereal *xv, integer *nxv)
{
    /* System generated locals */
    integer i__1;
    real ret_val;

    /* Local variables */
    static integer i__;

/* ----------------------------------- */
    /* Parameter adjustments */
    --xv;

    /* Function Body */
    if (*nxv > 1) {
	i__1 = *nxv;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (*x < xv[i__]) {
		goto L10;
	    }
/* L5: */
	}
	i__ = *nxv;
L10:
	if (i__ == 1) {
	    ++i__;
	}
	ret_val = (xv[i__ + *nxv] - xv[i__ + *nxv - 1]) / (xv[i__] - xv[i__ - 
		1]);
    } else {
	ret_val = (float)0.;
    }
    return ret_val;
} /* tange_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE VERCP */
/* ---------------------------------------------------------------------- */
/* --------------------------------------------------------------------- */
/* Subroutine */ int vercp_(doublereal *o, doublereal *v, doublereal *p, 
	doublereal *vc)
{
    extern /* Subroutine */ int prvu_(doublereal *, doublereal *, doublereal *
	    );
    static integer i__;
    static doublereal pmo[3];

/* -------------------------------------------- */
    /* Parameter adjustments */
    --vc;
    --p;
    --v;
    --o;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	pmo[i__ - 1] = p[i__] - o[i__];
/*        print *,i,Pmo(1),v(i),v(3+i) */
/* L5: */
    }
    prvu_(&v[4], pmo, &vc[1]);
    for (i__ = 1; i__ <= 3; ++i__) {
	vc[i__] += v[i__];
/* L10: */
    }
    return 0;
} /* vercp_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE PRVN */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int prvn_(doublereal *x, doublereal *y, doublereal *t, 
	doublereal *d__)
{
    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    extern /* Subroutine */ int prvu_(doublereal *, doublereal *, doublereal *
	    );
    static integer i__;
    static doublereal z__[3];

/* -------------------------------------- */
    /* Parameter adjustments */
    --t;
    --y;
    --x;

    /* Function Body */
    prvu_(&x[1], &y[1], z__);
    *d__ = z__[0] * z__[0] + z__[1] * z__[1] + z__[2] * z__[2];
    if (*d__ > (float)1e-12) {
	goto L5;
    }
    *d__ = (float)0.;
    return 0;
L5:
    *d__ = sqrt(*d__);
    for (i__ = 1; i__ <= 3; ++i__) {
	t[i__] = z__[i__ - 1] / *d__;
/* L10: */
    }
    return 0;
} /* prvn_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE RMOVE */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int rmove_(doublereal *rvi, doublereal *rvo, integer *n)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* ----------------------------------------- */
    /* Parameter adjustments */
    --rvo;
    --rvi;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rvo[i__] = rvi[i__];
/* L5: */
    }
    return 0;
} /* rmove_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE RONT */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int ront_(doublereal *x, doublereal *xyt, integer *np, 
	doublereal *rant)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__, im1;

/* -------------------------------------- */
    /* Parameter adjustments */
    --xyt;

    /* Function Body */
    if (*np > 1) {
	goto L10;
    }
    *rant = xyt[2];
    return 0;
L10:
    i__1 = *np;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (*x < xyt[i__]) {
	    goto L30;
	}
/* L20: */
    }
    i__ = *np;
L30:
    im1 = i__ - 1;
    *rant = (xyt[*np + i__] - xyt[*np + im1]) / (xyt[i__] - xyt[im1]) * (*x - 
	    xyt[im1]);
    *rant = xyt[*np + im1] + *rant;
    return 0;
} /* ront_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE RSET */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int rset_(doublereal *rv, integer *n, doublereal *rval)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    static integer i__;

/* --------------------------------------- */
    /* Parameter adjustments */
    --rv;

    /* Function Body */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	rv[i__] = *rval;
/* L10: */
    }
    return 0;
} /* rset_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE TMMNT */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int tmmnt_(doublereal *pa, doublereal *fa, doublereal *pb, 
	doublereal *fb)
{
    extern /* Subroutine */ int prvu_(doublereal *, doublereal *, doublereal *
	    );
    static integer i__;
    static doublereal t[3];

/* ------------------------------------- */
    /* Parameter adjustments */
    --fb;
    --pb;
    --fa;
    --pa;

    /* Function Body */
    for (i__ = 1; i__ <= 3; ++i__) {
	fb[i__] += fa[i__];
/* L10: */
	t[i__ - 1] = pa[i__] - pb[i__];
    }
    prvu_(t, &fa[1], t);
    for (i__ = 1; i__ <= 3; ++i__) {
/* L20: */
	fb[i__ + 3] += t[i__ - 1];
    }
    return 0;
} /* tmmnt_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE TRSFV */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int trsfv_(doublereal *tm, doublereal *rv, doublereal *rn, 
	integer *iver)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    extern /* Subroutine */ int mv3_(doublereal *, doublereal *, doublereal *,
	     integer *);
    static integer itr;

/* ------------------------------------------ */
    /* Parameter adjustments */
    --rn;
    --rv;
    tm -= 4;

    /* Function Body */
    itr = (i__1 = (*iver - 1) / 2, abs(i__1));
    mv3_(&tm[4], &rv[1], &rn[1], &itr);
    return 0;
} /* trsfv_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE PRVU */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int prvu_(doublereal *x, doublereal *y, doublereal *t)
{
    static doublereal z1, z2, z3;

/* ------------------------------------------ */
    /* Parameter adjustments */
    --t;
    --y;
    --x;

    /* Function Body */
    z1 = x[2] * y[3] - x[3] * y[2];
    z2 = x[3] * y[1] - x[1] * y[3];
    z3 = x[1] * y[2] - x[2] * y[1];
    t[1] = z1;
    t[2] = z2;
    t[3] = z3;
    return 0;
} /* prvu_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                   14- 1-83                      SUBROUTINE MV3 */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int mv3_(doublereal *rm, doublereal *vin, doublereal *vout, 
	integer *it)
{
    static doublereal v1, v2, v3;

/* ---------------------------------------------- */
    /* Parameter adjustments */
    --vout;
    --vin;
    --rm;

    /* Function Body */
    v1 = vin[1];
    v2 = vin[2];
    v3 = vin[3];
    if (*it > 0) {
	goto L10;
    }
    vout[1] = rm[1] * v1 + rm[4] * v2 + rm[7] * v3;
    vout[2] = rm[2] * v1 + rm[5] * v2 + rm[8] * v3;
    vout[3] = rm[3] * v1 + rm[6] * v2 + rm[9] * v3;
    return 0;
L10:
    vout[1] = rm[1] * v1 + rm[2] * v2 + rm[3] * v3;
    vout[2] = rm[4] * v1 + rm[5] * v2 + rm[6] * v3;
    vout[3] = rm[7] * v1 + rm[8] * v2 + rm[9] * v3;
    return 0;
} /* mv3_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                                              sub mux */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int mux_(doublereal *x, doublereal *alfa, doublereal *xyt, 
	integer *np, doublereal *rant)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double acos(doublereal);

    /* Local variables */
    static doublereal crid;
    static integer i__, im1;

/* -------------------------------------- */
    /* Parameter adjustments */
    --xyt;

    /* Function Body */
    crid = (float)1. - (float)2. / acos((float)-1.) * *alfa;
    if (*np > 1) {
	goto L10;
    }
    *rant = xyt[2];
    return 0;
L10:
    i__1 = *np;
    for (i__ = 1; i__ <= i__1; ++i__) {
	xyt[i__ + *np] *= crid;
/*      print *,'xyt',I,xyt(I+np) */
/* L15: */
    }
    i__1 = *np;
    for (i__ = 2; i__ <= i__1; ++i__) {
	if (*x < xyt[i__]) {
	    goto L30;
	}
/* L20: */
    }
    i__ = *np;
L30:
    im1 = i__ - 1;
    *rant = (xyt[*np + i__] - xyt[*np + im1]) / (xyt[i__] - xyt[im1]) * (*x - 
	    xyt[im1]);
    *rant = xyt[*np + im1] + *rant;
/*      print *,'rant',rant */
    return 0;
} /* mux_ */

/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/*                                              sub muy */
/* ---------------------------------------------------------------------- */
/* ---------------------------------------------------------------------- */
/* Subroutine */ int muy_(doublereal *x, doublereal *alfa, doublereal *
	alfarif, doublereal *rmust, doublereal *runt)
{
    /* Builtin functions */
    double acos(doublereal), tanh(doublereal);

    /* Local variables */
    static doublereal crod0, crod1;

/* -------------------------------------- */
/* rod Dipende da SR */
/*      crod = 1.-x */
    crod1 = *alfa / (acos((float)-1.) / (float)2.) * *rmust;
    crod0 = *rmust * tanh(*alfa / *alfarif);
    *runt = crod0 + (crod1 - crod0) / (float)1. * *x;
/*      runt = 0. */
/*      print *,'runt',runt */
    return 0;
} /* muy_ */

#ifdef __cplusplus
	}
#endif
