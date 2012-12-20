/* $Header$ */
/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
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
 * Authors:	Pierangelo Masarati <masarati@aero.polimi.it>
 * 		Tingnan Zhang <tingnan1986@gatech.edu>
 */

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cfloat>
#include <vector>

#include <thrust/device_vector.h>
#include <thrust/transform.h>
#include <thrust/inner_product.h>
#include <thrust/functional.h>

#include "cudatest.h"

class ViscoElasticCUDATest : public CUDATest {
private:
	thrust::host_vector<double> m_k_host, m_r_host;
	thrust::device_vector<double> m_k_device, m_r_device, m_f_device;

public:
	ViscoElasticCUDATest(const std::vector<double>& k, const std::vector<double>& r);
	virtual ~ViscoElasticCUDATest(void);
	virtual void GetForce(Vec3& F, const Vec3& X, const Vec3& V);
};

ViscoElasticCUDATest::ViscoElasticCUDATest(const std::vector<double>& k,
const std::vector<double>& r)
: m_k_host(k), m_r_host(r),
m_k_device(k), m_r_device(r),
m_f_device(k.size())
{
}

ViscoElasticCUDATest::~ViscoElasticCUDATest(void)
{
}

// Kelvin-Voight functor
struct kv_functor
{
	const double x, v;

	kv_functor(double x, double v) : x(x), v(v) {};

	__host__ __device__
	double operator()(const double& k, const double& r) const { 
		return -(k*x + r*v);
        };
};

void
ViscoElasticCUDATest::GetForce(Vec3& F, const Vec3& X, const Vec3& V)
{
	thrust::transform(m_k_device.begin(), m_k_device.end(),
		m_r_device.begin(), m_f_device.begin(),
		kv_functor(X(1), V(1)));
	F = Vec3(thrust::reduce(m_f_device.begin(), m_f_device.end(),
			0., thrust::plus<double>()),
		0., 0.);
}

extern "C" void *
mbdyn_CUDATest_init(unsigned n, double *pk, double *pr)
{
	std::vector<double> k(n), r(n);
	for (unsigned i = 0; i < n; i++) {
		k[i] = pk[i];
		r[i] = pr[i];
	}

	return new ViscoElasticCUDATest(k, r);
}

