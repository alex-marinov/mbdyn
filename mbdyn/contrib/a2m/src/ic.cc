//                                  IC.CC                                      

#include <ic.h>

s_ic::s_ic (void)  :   Aerror(1.0E-4),Alimit(30,DEGREE),Amaxit(25),
                       Apattern(new Boolean[10] ), Pattern(new Boolean[10]),
                       ncoord_pattern(0),ncoord_apattern(0),
                       Error(1.0E-10),Maxit(25),Tlimit(1.0E10),Verror(1E-4),
                       _Aerror(N),_Alimit(N),_Amaxit(N),_Error(N),_Maxit(N),
                       _Tlimit(N),_Verror(N),_Apattern(N),_Pattern(N)
{
   for (int i=0; i<10; i++) {
      Apattern[i]=Y;
      Pattern[i]=Y;
   }
}

s_ic::~s_ic(void)      {}

inline const char* const s_ic::Gettype (void) const
{
   return "IC";
}

Boolean s_ic::Test()
{
   /* Non necessita di controlli incrociati */
   return N;
}

ostream& s_ic::Print (ostream& out) const
{
   out << endl;
   out << "IC:" << endl;
   out << "     Aerror [" << _Aerror << "] = " << Aerror << endl
       << "     Alimit [" << _Alimit << "] = " << Alimit << endl
       << "     Amaxit [" << _Amaxit << "] = " << Amaxit << endl
       << "     Error  [" << _Error << "] = " << Error << endl
       << "     Maxit  [" << _Maxit << "] = " << Maxit << endl
       << "     Tlimit [" << _Tlimit << "] = " << Tlimit << endl
       << "     Verror [" << _Verror << "] = " << Verror << endl
       << "     Apattern [" << _Apattern << "] = ";
   if (_Apattern==Y) outvec (out,Apattern,ncoord_apattern);
   out << endl;
   out << "     Pattern [" << _Pattern << "] = ";
   if (_Pattern==Y) outvec(out,Pattern,ncoord_pattern);
   out << endl;
   out << endl;
   return out;
}

void s_ic::Translate (ostream& out)
{
   return;
}



