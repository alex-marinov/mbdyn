#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */

int
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

int
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
	return 1;
}

