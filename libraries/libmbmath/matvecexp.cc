/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2017
 * 
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati	<masarati@aero.polimi.it>
 * Paolo Mantegazza	<mantegazza@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2017
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 * Teodoro Merlini  <merlini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */


#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "matvecexp.h"

// NOTE: do not use Zero3, Zero3x3 or mb_zero<>()
// because they might be not initialized yet
const VecExp ZeroExp(0., 0., 0., 0., 0., 0.);
const MatExp EyeExp(Mat3x3(1., 0., 0., 0., 1., 0., 0., 0., 1.),
	Mat3x3(0., 0., 0., 0., 0., 0., 0., 0., 0.));



ScalExp
pow(const ScalExp &d, const doublereal &e)
{
#ifdef MBDYN_SINGLE_PRECISION
    double p = d.GetVec();
    return ScalExp(doublereal(pow(p, double(e))),
		    doublereal(d.GetMom()*(e*pow(p, double(e-1.)))));
#else /* ! MBDYN_SINGLE_PRECISION */
    doublereal p = d.GetVec();
    return ScalExp(pow(p, e), d.GetMom()*(e*pow(p, e-1.)));
#endif /* ! MBDYN_SINGLE_PRECISION */
};

ScalExp
sqrt(const ScalExp &d)
{
    doublereal p = sqrt(d.GetVec());
    return ScalExp(p, d.GetMom()/(p*2.));
};

ScalExp
sin(const ScalExp &d)
{
    doublereal p = d.GetVec();
    return ScalExp(sin(p), d.GetMom()*cos(p));
};

ScalExp
cos(const ScalExp &d)
{
    doublereal p = d.GetVec();
    return ScalExp(cos(p), d.GetMom()*(-sin(p)));
};

ScalExp
exp(const ScalExp &d)
{
    doublereal p = exp(d.GetVec());
    return ScalExp(p, d.GetMom()*p);
};




std::ostream&
ScalExp::Write(std::ostream& out, const char* sFill) const
{
	return out 
		<< GetVec() << sFill
		<< GetMom();
}

//ScalExp
//operator - (const ScalExp& v)
//{
//	return ScalExp(-v.GetVec(), -v.GetMom());
//}
//

std::ostream&
operator << (std::ostream& out, const ScalExp& v)
{
	return out 
		<< v.GetVec() << " "
		<< v.GetMom();    
}

std::ostream&
Write(std::ostream& out, const ScalExp& v, const char* sFill)
{
	return v.Write(out, sFill);
}

std::ostream&
VecExp::Write(std::ostream& out, const char* sFill) const
{
	return out 
		<< vec.dGet(1) << sFill
		<< vec.dGet(2) << sFill
		<< vec.dGet(3) << sFill
		<< mom.dGet(1) << sFill
		<< mom.dGet(2) << sFill
		<< mom.dGet(3);
}

VecExp
operator - (const VecExp& v)
{
	return VecExp(-v.GetVec(), -v.GetMom());
}

std::ostream&
operator << (std::ostream& out, const VecExp& v)
{
	const Vec3& vec = v.GetVec();
	const Vec3& mom = v.GetMom();
	return out 
		<< vec.dGet(1) << " "
		<< vec.dGet(2) << " "
		<< vec.dGet(3) << " " 
		<< mom.dGet(1) << " "
		<< mom.dGet(2) << " "
		<< mom.dGet(3);    
}

std::ostream&
Write(std::ostream& out, const VecExp& v, const char* sFill)
{
	return v.Write(out, sFill);
}

std::ostream&
MatExp::Write(std::ostream& out, const char* sFill, const char* sFill2) const
{
	char* sF2 = (char*)sFill2; 
	if (sFill2 == NULL) {      
		sF2 = (char*)sFill;
	}   
	
	return out 
		<< vec.dGet(1,1) << sFill
		<< vec.dGet(1,2) << sFill
		<< vec.dGet(1,3) << sFill
		<< mom.dGet(1,1) << sFill
		<< mom.dGet(1,2) << sFill
		<< mom.dGet(1,3) << sF2
		<< vec.dGet(2,1) << sFill
		<< vec.dGet(2,2) << sFill
		<< vec.dGet(2,3) << sFill
		<< mom.dGet(2,1) << sFill
		<< mom.dGet(2,2) << sFill
		<< mom.dGet(2,3) << sF2
		<< vec.dGet(3,1) << sFill
		<< vec.dGet(3,2) << sFill
		<< vec.dGet(3,3) << sFill
		<< mom.dGet(3,1) << sFill
		<< mom.dGet(3,2) << sFill
		<< mom.dGet(3,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(1,1) << sFill
		<< vec.dGet(1,2) << sFill
		<< vec.dGet(1,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(2,1) << sFill
		<< vec.dGet(2,2) << sFill
		<< vec.dGet(2,3) << sF2
		<< 0. << sFill
		<< 0. << sFill
		<< 0. << sFill
		<< vec.dGet(3,1) << sFill
		<< vec.dGet(3,2) << sFill
		<< vec.dGet(3,3);
}

std::ostream&
operator << (std::ostream& out, const MatExp& m)
{
	const Mat3x3& vec = m.GetVec();
	const Mat3x3& mom = m.GetMom();
	return out 
		<< vec.dGet(1, 1) << " "
		<< vec.dGet(1, 2) << " "
		<< vec.dGet(1, 3) << " " 
		<< mom.dGet(1, 1) << " "
		<< mom.dGet(1, 2) << " "
		<< mom.dGet(1, 3) << std::endl
		<< vec.dGet(2, 1) << " "
		<< vec.dGet(2, 2) << " "
		<< vec.dGet(2, 3) << " " 
		<< mom.dGet(2, 1) << " "
		<< mom.dGet(2, 2) << " "
		<< mom.dGet(2, 3) << std::endl
		<< vec.dGet(3, 1) << " "
		<< vec.dGet(3, 2) << " "
		<< vec.dGet(3, 3) << " " 
		<< mom.dGet(3, 1) << " "
		<< mom.dGet(3, 2) << " "
		<< mom.dGet(3, 3) << std::endl
		<< 0. << " "
		<< 0. << " "
		<< 0. << " " 
		<< vec.dGet(1, 1) << " "
		<< vec.dGet(1, 2) << " "
		<< vec.dGet(1, 3) << std::endl
		<< 0. << " "
		<< 0. << " " 
		<< 0. << " " 
		<< vec.dGet(2, 1) << " "
		<< vec.dGet(2, 2) << " "
		<< vec.dGet(2, 3) << std::endl
		<< 0. << " "
		<< 0. << " "
		<< 0. << " " 
		<< vec.dGet(3, 1) << " "
		<< vec.dGet(3, 2) << " "
		<< vec.dGet(3, 3);     
}

std::ostream&
Write(std::ostream& out, const MatExp& m, const char* sFill, const char* sFill2)
{
	return m.Write(out, sFill, sFill2);
}

#if 0
VecExp
MultRV(const VecExp& v, const Mat3x3& R)
{
	return VecExp(R*v.GetVec(),R*v.GetMom());
}

MatExp
MultRM(const MatExp& m, const Mat3x3& R)
{
	return MatExp(R*m.GetVec(), R*m.GetMom());
}


MatExp
MultMRt(const MatExp& m, const Mat3x3& R)
{
	Mat3x3 Rt = R.Transpose();
	return MatExp(m.GetVec()*Rt, m.GetMom()*Rt);
}

MatExp
MultRMRt(const MatExp& m, const Mat3x3& R)
{   
	Mat3x3 Rt = R.Transpose();
	return MatExp(R*m.GetVec()*Rt, R*m.GetMom()*Rt);
}
#endif /* 0 */

