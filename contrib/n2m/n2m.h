/******************************************************************************
 * 
 * NASTRAN 2 MBDyn bulk data conversion tool
 * 
 * Copyright (c) 2003
 * 
 * Author:         Pierangelo Masarati <masarati@aero.polimi.it>
 * Affiliation:    Dipartimento di Ingegneria Aerospaziale
 *                 Politecnico di Milano
 *                 via La Masa 34, 20156 Milano (Italy)
 *                 http://www.aero.polimi.it
 * 
 * This program is part of the MBDyn package
 * http://www.mbdyn.org
 * https://mbdyn.aero.polimi.it/~masarati/mbdyn-index.html
 * 
 *****************************************************************************/

#ifndef N2M_H
#define N2M_H

enum {
	NASTRAN_CARD_TYPE = 0,
	NASTRAN_FIRST = 0,
	NASTRAN_SECOND,
	NASTRAN_THIRD,
	NASTRAN_FOURTH,
	NASTRAN_FIFTH,
	NASTRAN_SIXTH,
	NASTRAN_SEVENTH,
	NASTRAN_EIGHTH,
	NASTRAN_NINTH,
	NASTRAN_TENTH,
	NASTRAN_LAST = NASTRAN_TENTH,
	NASTRAN_MAX_TOKEN_NUMBER
}; 

enum {
	NASTRAN_CARD_GRID = 0,
	NASTRAN_CARD_CONM2,
	NASTRAN_CARD_CORD2R,
	NASTRAN_CARD_PBEAM,

	NASTRAN_CARD_LAST		/* always goes last */
};

enum {
	NASTRAN_FILE_IN = 0,
	NASTRAN_FILE_OUT_CTL,
	NASTRAN_FILE_OUT_REF,
	NASTRAN_FILE_OUT_NODE,
	NASTRAN_FILE_OUT_ELEM,
	NASTRAN_FILE_OUT_ERR,
	NASTRAN_FILE_LAST
};

#define NASTRAN_TOKEN_SIZE 8
enum {
	N2M_GOT_TOKEN = 0,
	N2M_USE_DEFAULT,
	N2M_OUT_OF_RANGE,
	N2M_INSUFFICIENT_TOKEN_SIZE
};

#define N_BUFFER_SIZE 1024
struct n2m_buffer {
	char buf[N_BUFFER_SIZE];
	int len;
};

struct n2m_card_data {
	int num;
};

/* global data */
extern int make_rigid_bodies;
extern int reference_offset;

/* helpers */
extern int get_buf(FILE *fin, struct n2m_buffer *b);
extern char *get_string(struct n2m_buffer *b, int number, struct n2m_buffer *buf, char *def);
extern int get_int(struct n2m_buffer *b, int number, int *def);
extern double get_double(struct n2m_buffer *b, int number, double *def);

/* CONM2 */
extern int do_conm2(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form);
extern int do_conm2_cont(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form, int def);

/* CORD2R */
int do_cord2r(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form, double *coords);
int do_cord2r_cont(struct n2m_buffer *b, FILE **f, struct n2m_buffer *form, double *coords);

/* GRID */
extern int do_grid(struct n2m_buffer *b, FILE **f);

/* PBEAM */
int do_pbeam_cont(struct n2m_buffer *b, FILE **f, int *cont);
int do_pbeam(struct n2m_buffer *b, FILE **f);

#endif /* N2M_H */

