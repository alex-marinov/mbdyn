//                                OUTPUT.CC                                   


#include <output.h>


ostream& operator << (ostream& out, const coord_type& R)
{
   switch(R) {
    case _XC: out << "X"; break;
    case _YC: out << "Y"; break;
    case _ZC: out << "Z"; break;
    case _PSIC : out << "PSI"; break;
    case _THETAC : out << "THETA"; break;
    case _PHIC : out << "PHI"; break;
    default: out << "NULL"; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Joint& R)
{
   switch (R) {
    case _CYLINDRICAL   : out << "CYLINDRICAL"; break;
    case _CONVEL	: out << "CONVEL"; break;
    case _REVOLUTE	: out << "REVOLUTE"; break;
    case _SPHERICAL	: out << "SPHERICAL"; break;
    case _SCREW         : out << "SCREW"; break;
    case _FIXED         : out << "FIXED"; break;
    case _HOOKE         : out << "HOOKE"; break;
    case _TRANSLATIONAL : out << "TRANSLATIONAL"; break;
    case _UNIVERSAL     : out << "UNIVERSAL"; break;
    case _RACKPIN       : out << "RACKPIN"; break;
    case _PLANAR        : out << "PLANAR"; break;
    default		: out << "UNCLASSIFIED."; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Boolean& R)
{
   if (R==Y) out << "TRUE"; else out << "FALSE";
   return out;
}

ostream& operator << (ostream& out, const Friction& R)
{
   switch(R) {
    case _OFF         : out << "OFF"; break;
    case _ON          : out << "ON"; break;
    case _PRELOAD_ONLY: out << "PRELOAD_ONLY"; break;
   }
   return out;
}

ostream& operator << (ostream& out, const Joint_Primitive& R)
{
   switch (R) {
    case _ATPOINT: out << "ATPOINT"; break;
    case _INLINE: out << "INLINE"; break;
    case _INPLANE: out << "INPLANE"; break;
    case _ORIENTATION: out << "ORIENTATION"; break;
    case _PARALLEL_AXES: out << "PARALLEL AXES"; break;
    case _PERPENDICULAR: out << "PERPENDICULAR"; break;
   }
   return out;
}

ostream& operator <<  (ostream& out, const Direction_Mode& R)
{
   switch (R) {
    case _ROTATION: out << "ROTATION"; break;
    case _TRANSLATION: out << "TRANSLATION"; break;
   }
   return out;
}


