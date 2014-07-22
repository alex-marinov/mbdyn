/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2014
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

/* Input */

#ifndef INPUT_H
#define INPUT_H

#include <iostream>
#include <myassert.h>

/* Filtro per la classe istream che conta il numero di righe. */

/* InputStream - begin */

class InputStream {
	friend InputStream& operator >> (InputStream& in, int& i);
	friend InputStream& operator >> (InputStream& in, long int& i);
	friend InputStream& operator >> (InputStream& in, short int& i);
	friend InputStream& operator >> (InputStream& in, unsigned int& i);
	friend InputStream& operator >> (InputStream& in, unsigned long int& i);
	friend InputStream& operator >> (InputStream& in, unsigned short int& i);
	friend InputStream& operator >> (InputStream& in, char& i);
	friend InputStream& operator >> (InputStream& in, float& i);
	friend InputStream& operator >> (InputStream& in, double& i);

private:
	std::istream& iStrm;
	unsigned long uLineNumber;      
   
public:
	/* Costruttore - inizializza il filtro con un reference ad un istream */
	InputStream(std::istream& in);
   
	/* Distruttore banale */
	~InputStream(void);
   
	/* Legge un carattere; se e' un fine-riga, aggiorna il contatore */
	inline char get(void);
   
	/* Legge un carattere; se e' un fine-riga, aggiorna il contatore */
	inline std::istream& get(char& ch);
   
	/* Esegue il putback di un carattere */
	inline InputStream& putback(char ch);
   
	/* Restituisce il valore del contatore */
	inline unsigned long int GetLineNumber(void) const;
   
	/* eof */
	inline bool eof(void) const;
   
	/* Restituisce l'istream */
	inline const std::istream& GetStream(void) const;
	inline std::istream& GetStream(void);
};

/* Overload dell'operatore di lettura */
extern InputStream& operator >> (InputStream& in, int& i);
extern InputStream& operator >> (InputStream& in, long int& i);
extern InputStream& operator >> (InputStream& in, short int& i);
extern InputStream& operator >> (InputStream& in, unsigned int& i);
extern InputStream& operator >> (InputStream& in, unsigned long int& i);
extern InputStream& operator >> (InputStream& in, unsigned short int& i);
extern InputStream& operator >> (InputStream& in, char& i);
extern InputStream& operator >> (InputStream& in, float& i);
extern InputStream& operator >> (InputStream& in, double& i);

/* Legge un carattere; se e' un fine-riga, aggiorna il contatore */
inline char
InputStream::get(void)
{
	char ch = iStrm.get();
	if (ch == '\n') {
		uLineNumber++;
	}
	return ch;
}
   
/* Legge un carattere; se e' un fine-riga, aggiorna il contatore */
inline std::istream&
InputStream::get(char& ch) 
{
	std::istream& i = iStrm.get(ch);
	if (ch == '\n') {
		uLineNumber++;
	}
	return i;
}

/* Esegue il putback di un carattere */
inline InputStream&
InputStream::putback(char ch)
{
	iStrm.putback(ch);
	if (ch == '\n') {
		uLineNumber--;
	}
	return *this;
}

/* Restituisce il valore del contatore */
inline unsigned long int
InputStream::GetLineNumber(void) const 
{
	return uLineNumber;
}

/* eof */
inline bool
InputStream::eof(void) const
{
	return iStrm.eof();
}

/* Restituisce l'istream */
inline const std::istream&
InputStream::GetStream(void) const
{
	return iStrm;
}

inline std::istream&
InputStream::GetStream(void)
{
	return iStrm;
}

/* InputStream - end */

#endif /* INPUT_H */

