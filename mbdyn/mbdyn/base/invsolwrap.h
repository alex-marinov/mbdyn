/*
 * invsolwrap.h
 *
 *  Created on: 23/giu/2016
 *      Author: Simone Conso
 */

#ifndef INVSOLWRAP_H_
#define INVSOLWRAP_H_

extern void *
mb_sol_create(const char *sIn, const char *sOut);

extern void *
mb_sol_create_inv(const char *sIn, const char *sOut);

extern int
mb_sol_prepare(void *p);

extern int
mb_sol_start(void *p);

extern int
mb_sol_advance(void *p);

extern int
mb_sol_destroy(void *p);

extern int
mb_sol_setbufin(void *p, unsigned uLabel, int iSize, double *pdBuf);

extern int
mb_sol_setbufout(void *p, unsigned uLabel, int iSize, double *pdBuf);

#endif /* INVSOLWRAP_H_ */
