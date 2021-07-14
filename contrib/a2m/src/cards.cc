/*

MBDyn (C) is a multibody analysis code. 
http://www.mbdyn.org

Copyright (C) 1996-2007

Pierangelo Masarati	<masarati@aero.polimi.it>
Paolo Mantegazza	<mantegazza@aero.polimi.it>

Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
via La Masa, 34 - 20156 Milano, Italy
http://www.aero.polimi.it

Changing this copyright notice is forbidden.

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


------------------------------------------------------------------------------

ADAMS2MBDyn (C) is a translator from ADAMS/View models in adm format
into raw MBDyn input files.

Copyright (C) 1999-2007
Leonardo Cassan		<lcassan@tiscalinet.it>

*/


#include <cards.h>

// GENERIC CARD CONSTRUCTOR

extern inline const char* Find_Token (int);

s_card::s_card(void)	      : label (0), _remark_ (NULL) {}

void s_card::Comment (char * p)
{
   _remark_=new char[strlen(p)];
   for (int i=0;i<strlen(p);i++)
     _remark_[i]=p[i];
   return;
}

void s_card::Store_Formula (char* tk, char *txt)
{
   if (txt[0]!='\0') {
     char* frm = new char[strlen(txt)+1];
     for (int i=0; i<strlen(txt); i++)
       frm[i]=txt[i];
     frm[strlen(txt)]='\0';
     Recipient[tk]=frm;
	/*.insert (formula_entry(prova,frm));*/
   }
   return;
}

void s_card::Display_Formula(ostream& out) const
{
   formula_map copia = Recipient;
   p_formula_entry p;
   if (Recipient.size()) out << "      [FORMULA DEFINED IN THE CARD:]" << endl;
   for (p=copia.begin();p!=copia.end();p++)
     out << "      " << (*p).first << " " << (*p).second << endl;
   return;
}
