#ifndef ___HYDRODYNAMIC_PLAIN_BEARING_H__INCLUDED___
#define ___HYDRODYNAMIC_PLAIN_BEARING_H__INCLUDED___

#include <ac/f2c.h>

#ifdef __cplusplus
extern "C" 
{
#endif
	
static const int NBDIRSMAX = 6;

    struct bearing_data
    {
        doublereal b, d, Psi, eta, eps_max, s, a[9];
    };

    void hydrodynamic_plain_bearing_init(bearing_data& bdat);
    
// computes the hydrodynamic bearing force at the bearing
	
    void hydrodynamic_plain_bearing_force(const bearing_data& bdat,
					const doublereal omega[2], 
					const doublereal e[2], 
					const doublereal e_dot[2], 
					doublereal k[3],	
					doublereal &eps,
					doublereal &eps_dot,
					doublereal& delta,
					doublereal& SoD,
					doublereal& SoV,
					doublereal& my,
                                          doublereal& beta);
					
// computes the hydrodynamic stiffness matrix and damping matrix at the bearing side

    void hydrodynamic_plain_bearing_force_dv(const bearing_data& bdat,
						const doublereal omega[2],
						const doublereal omegad[2][NBDIRSMAX],
						const doublereal e[2], 
						const doublereal ed[2][NBDIRSMAX],
						const doublereal e_dot[2], 
						const doublereal e_dotd[2][NBDIRSMAX],
						doublereal k[3],
						doublereal kd[3][NBDIRSMAX], 
						doublereal& eps, 
						doublereal& eps_dot, 
						doublereal& delta, 
						doublereal& SoD, 
						doublereal& SoV, 
						doublereal& my, 
						doublereal& beta,
						const int& nbdirs);

#ifdef __cplusplus
} // extern "C"
#endif
												
#endif												
