//                                EQUILIBRIUM                                 

// Definizioni per la card Equilibrium

#include <equilibrium.h>

s_equilibrium::s_equilibrium (void) : Alimit(0,DEGREE),Error(1.0E-4),
                                      Imbalance(1.0E-4),Maxit(25),
                                      Stability(1.0E-5),Tlimit(1.0E10),
                                      Pattern(new Boolean[10]),
                                      ncoord_pattern(0)
{
   for (int i=0; i<10; i++) Pattern[i]=Y;
}

s_equilibrium::~s_equilibrium (void)  {
}

Boolean s_equilibrium::Test()
{
   /* Non si rende necessario l'utilizzo di controlli incrociati */
   return N;
}

inline const char* const s_equilibrium::Gettype(void) const
{
   return "EQUILIBRIUM";
}

ostream& s_equilibrium::Print (ostream& out) const
{
   out << endl;
   out << "EQUILIBRIUM:" << endl;
   out << "     " << "Alimit [" << _Alimit << "] = " << Alimit << endl;
   out << "     " << "Error [" << _Error << "] = " << Error << endl
       << "     " << "Imbalance [" << _Imbalance << "] = "
       << Imbalance << endl 
       << "     " << "Max iterations [" << _Maxit << "] = "
       << Maxit << endl
       << "     " << "Stability [" << _Stability << "] = "
       << Stability << endl
       << "     " << "Tlimit [" << _Tlimit << "] = " 
       << Tlimit << endl;
   out << "     " << "Pattern [" << _Pattern << "] = ";
   if (_Pattern==Y) outvec(out,Pattern,ncoord_pattern);
   out << endl << endl;
   return out;
}

void s_equilibrium::Translate (ostream& out)
{
   return;
}
