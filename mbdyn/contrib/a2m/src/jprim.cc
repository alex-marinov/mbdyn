//                                 JPRIM.CC                                   

#include <jprim.h>
#include <output.h>

s_jprim::s_jprim (void)       : Node1(0),Node2(0),Type(_INPLANE),
                                _Node1(N),_Node2(N),_Type(N)
                                {}

s_jprim::~s_jprim (void)      {}


Boolean s_jprim::Test()
{
   const int err_before=nerr;
   // VERIFICA DELL'ESISTENZA DI ENTRAMBI GLI ID DEI MARKER I E J
   if (_Node1==N) out_error (27,"ID I");
   if (_Node2==N) out_error (27,"ID J");
   // VERIFICA CHE SIA STATO SPECIFICATO IL TIPO DI PRIMITIVA
   if (_Type==N) out_error (28,"");
   if (err_before != nerr) return Y; else return N;
}

inline const char* const s_jprim::Gettype (void) const
{
   return "JPRIM";
}

ostream& s_jprim::Print (ostream& out) const
{
   out << endl;
   out << "JPRIM:" << label << "      TYPE:" << Type << endl;
   out << "     Node1 [" << _Node1 << "] = " << Node1 << endl
       << "     Node2 [" << _Node2 << "] = " << Node2 << endl;
   return out;
}

void s_jprim::Translate (ostream& out)
{
   return;
}


