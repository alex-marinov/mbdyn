#ifndef MATERIAL_H
#define MATERIAL_H

#include <common.h>

struct s_material : public s_card {
   double Young;
   double Density;
   double P_ratio;
   char *Name;
   
   Boolean _Young,_Density,_Name,_P_ratio;
   
   s_material(void);
   ~s_material (void);
   ostream& Print (ostream& out) const;
   void Translate(ostream& out);
   inline const char* const Gettype (void) const;
   Boolean Test();
   
   template <class T>
     void Set_Part (int param, T value)
     {
	switch (param) {
	 case YOUNGS_MODULUS : Set_Param (Young,value,_Young); break;
	 case DENSITY : Set_Param (Density, value, _Density); break;
	 case NAME : Set_Param (Name, value, _Name); break;
	 case POISSONS_RATIO : Set_Param (P_ratio, value, _P_ratio); break; 
	 default : out_error (3,Find_Token(param)); break;
	}
	if (ES) out_error (1,Find_Token(param));
     }
};

#endif
