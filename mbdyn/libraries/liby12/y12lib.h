#ifndef Y12LIB_H
#define Y12LIB_H

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

extern int 
__FC_DECL__(y12mbf) (integer *n, integer *z__, doublereal *a, integer *snr, 
		   integer *nn, integer *rnr, integer *nn1, integer *ha,
		   integer *iha, doublereal *aflag, integer *iflag,
		   integer *ifail);
		   
extern int 
__FC_DECL__(y12mcf) (integer *n, integer *z__, doublereal *a, integer *snr,
		   integer *nn, integer *rnr, integer *nn1, doublereal *pivot,
		   doublereal *b, integer *ha, integer *iha, doublereal *aflag,
		   integer *iflag, integer *ifail);
		   
extern int 
__FC_DECL__(y12mdf) (integer *n, doublereal *a, integer *nn, doublereal *b,
		   doublereal *pivot, integer *snr, integer *ha, integer *iha,
		   integer *iflag, integer *ifail);

#ifdef __cplusplus 
}
#endif /* __cplusplus */

#endif /* Y12LIB_H */

