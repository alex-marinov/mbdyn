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
