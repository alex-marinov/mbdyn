#ifndef POINTMASS_H
#define POINTMASS_H

#include <common.h>
#include <output.h>

struct s_pointmass: public s_card {
   double Mass;
   Id Cm;
   Euler Reuler;
   Vec3 Qg,Zg,Xg;
   double Vx,Vy,Vz;
   coord_type* Exact;
   int ncoord;
   Boolean _Mass, _Cm, _Reuler, _Qg, _Zg, _Xg, _Vx, _Vy, _Vz,
           _Exact;

   enum Deftype {
      ANGLES,
      POINT,
      EYEMATRIX
   };
   
   Deftype Mode;

   
   s_pointmass (void);
   ~s_pointmass (void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case MASS: Set_Param (Mass,value, _Mass); break;
	 case CM: Set_Param (Cm, value, _Cm); break;
	 case VX: Set_Param (Vx, value, _Vx); break;
	 case VY: Set_Param (Vy, value, _Vy); break;
	 case VZ: Set_Param (Vz, value, _Vz); break;
	 case QG: Set_Param (Qg, value, _Qg); break;
	 case XG: Set_Param (Xg, value, _Xg); break;
	 case ZG: Set_Param (Zg, value, _Zg); break;
	 case EXACT: Set_Array (Exact, value,ncoord, _Exact); break;
	 case REULER: Set_Param (Reuler, value, _Reuler); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
