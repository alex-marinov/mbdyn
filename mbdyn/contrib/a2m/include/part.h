
#ifndef PART_H
#define PART_H

#include <common.h>
#include <output.h>

struct s_part: public s_card {
   double Mass;
   Id Cm;
   Id Im;
   Vec3 Ip,AddIp;
   char* Material;
   Euler Reuler;
   Vec3 Qg,Zg,Xg;
   double Vx,Vy,Vz;
   double Wx,Wy,Wz;
   Id Wm,Vm;
   Boolean Ground;
   coord_type* Exact;
   int ncoord;
   
   Boolean _Mass,_Cm,_Im,_Material,_Qg,_Zg,_Xg,_Vm,_Wm,
           _Ip,_AddIp,_Exact,_Reuler,
           _Vx,_Vy,_Vz,_Wx,_Wy,_Wz;

   enum Deftype {
      ANGLES,
      POINT,
      EYEMATRIX
   };
   
   Deftype Mode;

   
   s_part (void);
   ~s_part (void);
   ostream& Print (ostream& out) const;
   void Translate (ostream& out);
   inline const char* const Gettype (void) const;   
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch(param) {
	 case MASS: Set_Param (Mass,value,_Mass); break;
	 case CM: Set_Param (Cm, value,_Cm); break;
	 case IM: Set_Param (Im, value,_Im); break;
	 case MATERIAL: Set_Param(Material,value,_Material); break;
	 case WM: Set_Param(Wm, value,_Wm); break;
	 case VM: Set_Param(Vm, value,_Vm); break;
	 case IP: Set_Param(Ip,value,_Ip); break;
	 case QG: Set_Param(Qg,value,_Qg); break;
	 case XG: Set_Param(Xg,value,_Xg); break;
	 case ZG: Set_Param(Zg,value,_Zg); break;
	 case VX: Set_Param(Vx,value,_Vx); break;
	 case VY: Set_Param(Vy,value,_Vy); break;
	 case VZ: Set_Param(Vz,value,_Vz); break;
	 // Case ADDIP: Non essendoci il token viene gestito direttamente
	 case WX: Set_Param(Wx,value,_Wx); break;
	 case WZ: Set_Param(Wz,value,_Wz); break;
	 case WY: Set_Param(Wy,value,_Wy); break;
	 case REULER: Set_Param(Reuler,value,_Reuler); break;
	 case EXACT: Set_Array (Exact,value,ncoord,_Exact); break;
	 case GROUND: Set_Param(Ground,Y,Ground); break;
	 default: out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
	return;
     }
};

#endif
