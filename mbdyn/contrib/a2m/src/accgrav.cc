//                              ACCGRAV.CC                                    

#include <accgrav.h>

extern MBDyn_deck MBAccs;

s_accgrav::s_accgrav (void)   : Igrav(0),Jgrav(0),Kgrav(0),
                               _Igrav(N),_Jgrav(N),_Kgrav(N)
                               {}

Boolean s_accgrav::Test()
{
   /* non necessita di controlli incrociati */
   return N;
}

inline const char* const s_accgrav::Gettype (void) const
{
   return "ACCGRAV";
}

ostream& s_accgrav::Print (ostream& out) const
{
   out << endl;
   out << "ACCGRAV" << endl
     << "     I COMPONENT [" << _Igrav << "] = " << Igrav << endl
     << "     J COMPONENT [" << _Jgrav << "] = " << Jgrav << endl
     << "     K COMPONENT [" << _Kgrav << "] = " << Kgrav << endl;
   return out;
}

void s_accgrav::Translate(ostream& out)
{
   Vec3 G(Igrav,Jgrav,Kgrav);
   /* La gravità è espressa nel sistema globale.
    * Il template drive_caller di tipo Vec3 è il risultato
    * di una moltiplicazione del vettore Accellerazione per
    * la costante 1 */
   MBDyn_gravity* GRAVITY = new MBDyn_gravity (G);
   MBAccs.insert (MBDyn_entry(999,(MBDyn_card*) GRAVITY));
   return;
}

