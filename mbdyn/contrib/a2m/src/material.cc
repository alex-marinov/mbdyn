//                              MATERIAL.CC                                   

#include <material.h>

s_material::s_material (void) : Young(0),Density(0),Name(),P_ratio(0),
                               _Young(N),_Density(N),_Name(N),_P_ratio(N) 
                               { Name=new char[80]; }

s_material::~s_material (void)  {}


Boolean s_material::Test()
{
   /* Rimanda in un secondo tempo */
   return N;
}

inline const char* const s_material::Gettype (void) const
{
   return "MATERIAL";
}

ostream& s_material::Print (ostream& out) const
{
   out << endl;
   out << "MATERIAL:" << label << "       NAME:" << Name << endl
     << "     YOUNG MODULUS [" << _Young << "] = " << Young << endl
     << "     DENSITY [" << _Density << "] = " << Density << endl
     << "     POISSONS RATIO [" << _P_ratio << "] = " << P_ratio << endl;
     return out;
}

void s_material::Translate(ostream& out)
{
   return;
}


