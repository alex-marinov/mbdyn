#ifndef WHEEL_H
#define WHEEL_H
#ifdef __cplusplus
extern "C" {
#endif
extern int __FC_DECL__(wvefr)(doublereal *vw, 
			      doublereal *ome, 
			      doublereal *cc, 
			      doublereal *res,
			      doublereal *rm,
			      integer *ircw,
			      integer *ipc,
			      doublereal *omep,
			      doublereal *tabfat,
			      integer *npft,
			      doublereal *beta,
			      doublereal *betap,
			      doublereal *bleng,
			      doublereal *pp,
			      doublereal *wp,
			      doublereal *gp);
#ifdef __cplusplus
}
#endif
#endif
