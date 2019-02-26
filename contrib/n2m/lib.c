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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include <n2m.h>

static int
token_start[] = {
	0,	/* Type/Continuation */
	8,	/* Args ... */
	16,
	24,
	32,
	40,
	48,
	56,
	64,
	72	/* Continuation */
};

static int
trim(char *s)
{
	while (isspace(*s)) {
		memmove(s, s+1, strlen(s));
		if (*s == '\0') {
			return 0;
		}
	}
	while (*s != '\0' && !isspace(*s)) {
		s++;
	}
	*s = '\0';
	
	return 0;
}

static int
double_exp(char *s)
{
	if (*s == '+' || *s == '-') {
		s++;
	}
	
	while (isdigit(*s) || *s == '.') {
		s++;
	}
	
	switch (*s) {
	case '\0':
		return 0;
        case 'D':
	case 'd':
	case 'E':
	case 'e':
	case 'F':
	case 'f':
	case 'G':
	case 'g':
		return 0;
	case '-':
	case '+':
		memmove(s+1, s, strlen(s)+1);
		*s = 'e';
		return 0;
	}
	return -1;
}

int
get_buf(FILE *fin, struct n2m_buffer *b)
{
	char *res;

	res = fgets(b->buf, sizeof(b->buf), fin);
	if (res != NULL) {
		b->len = strlen(b->buf);
	} else {
		b->len = 0;
	}
	return (res != NULL);
}

static int
get_token(struct n2m_buffer *b, char *token, int size, int number)
{
	if (number < 0 || number > NASTRAN_LAST) {
		return N2M_OUT_OF_RANGE;
	}
	
	if (token_start[number] >= b->len) {
		token[0] = '\0';
		return N2M_USE_DEFAULT;
	}
	
	if (size <= NASTRAN_TOKEN_SIZE) {
		return N2M_INSUFFICIENT_TOKEN_SIZE;
	}
	
	strncpy(token, b->buf+NASTRAN_TOKEN_SIZE*number, NASTRAN_TOKEN_SIZE);
	token[NASTRAN_TOKEN_SIZE] = '\0';

	return N2M_GOT_TOKEN;
}

char *
get_string(struct n2m_buffer *b, int number, struct n2m_buffer *buf, char *def)
{
	char token[9];
	switch (get_token(b, token, sizeof(token), number)) {
	case N2M_GOT_TOKEN:
		trim(token);
		strncpy(buf->buf, token, sizeof(buf->buf));
		return buf->buf;
	case N2M_USE_DEFAULT:
		if (def == NULL) {
			fprintf(stderr,
				"### no default value is allowed for string token n. %d of line\n%s\nbailing out ...\n",
				number+1, b->buf);
			exit(EXIT_FAILURE);
		}
		strncpy(buf->buf, def, sizeof(buf->buf));
		return buf->buf;
	case N2M_OUT_OF_RANGE:
	case N2M_INSUFFICIENT_TOKEN_SIZE:
	default:
		fprintf(stderr,
			"### internal error while reading int value from token n. %d of line\n%s\nbailing out ...\n", 
			number+1, b->buf);
		exit(EXIT_FAILURE);
	}
}


int
get_int(struct n2m_buffer *b, int number, int *def)
{
	char token[9];
	
	switch (get_token(b, token, sizeof(token), number)) {
	case N2M_GOT_TOKEN:
		trim(token);
		return atoi(token);
	case N2M_USE_DEFAULT:
		if (def == NULL) {
			fprintf(stderr,
				"### no default value is allowed for int token n. %d of line\n%s\nbailing out ...\n",
				number+1, b->buf);
			exit(EXIT_FAILURE);
		}
		return *def;
	case N2M_OUT_OF_RANGE:
	case N2M_INSUFFICIENT_TOKEN_SIZE:
	default:
		fprintf(stderr, 
			"### internal error while reading int value from token n. %d of line\n%s\nbailing out ...\n",
			number+1, b->buf);
		exit(EXIT_FAILURE);
	}
}

double
get_double(struct n2m_buffer *b, int number, double *def)
{
	char token[10];
	
	switch (get_token(b, token, sizeof(token), number)) {
	case N2M_GOT_TOKEN:
		trim(token);
		double_exp(token);
		return atof(token);
	case N2M_USE_DEFAULT:
		if (def == NULL) {
			fprintf(stderr,
				"### no default value is allowed for double token n. %d of line\n%s\nbailing out ...\n",
				number+1, b->buf);
			exit(EXIT_FAILURE);
		}
		return *def;
	case N2M_OUT_OF_RANGE:
	case N2M_INSUFFICIENT_TOKEN_SIZE:
	default:
		fprintf(stderr,
			"### internal error while reading double value from token n. %d of line\n%s\nbailing out ...\n",
			number+1, b->buf);
		exit(EXIT_FAILURE);
	}
}

