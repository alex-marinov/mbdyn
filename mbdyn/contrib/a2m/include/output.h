#ifndef OUTPUT_H
#define OUTPUT_H

// FUNZIONI NECESSARIE ALLA CORRETTA VISUALIZZAZIONE DELL'OUTPUT

#include <defs.h>
#include <iostream.h>

ostream& operator << (ostream& out, const Joint&);
ostream& operator << (ostream& out, const Friction&) ;
ostream& operator << (ostream& out, const Joint_Primitive&);
ostream& operator << (ostream& out, const Direction_Mode&);
ostream& operator << (ostream& out, const Boolean&);
ostream& operator << (ostream& out, const coord_type&);

template <class T>
ostream& outvec (ostream& out, T* data, int idx)
{
   for (int k=0; k<idx; k++)
     out << data[k] << " ";
   return out;
}

#endif
