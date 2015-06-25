/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2015
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

#ifndef FRICTION_H
#define FRICTION_H

#include <cfloat>
#include <limits>

doublereal 
norm(const doublereal &v)
{
	return fabs(v);
}

doublereal
dir(const doublereal &v)
{
	return copysign(1., v);
}

doublereal
norm(const Vec3& v)
{
	return v.Norm();
}

Vec3
dir(const Vec3& v)
{
	doublereal d = v.Norm();
	if (d < std::numeric_limits<doublereal>::epsilon()) {
		return Zero3;
	}
	return v/d;
}

Vec3
tanh(const Vec3& v)
{
	doublereal d = v.Norm();
	if (d < std::numeric_limits<doublereal>::epsilon()) {
		return Zero3;
	}
	return v*(tanh(d)/d);
}

Vec3
sin(const Vec3& v)
{
	doublereal d = v.Norm();
	if (d < std::numeric_limits<doublereal>::epsilon()) {
		return Zero3;
	}
	return v*(sin(d)/d);
}

template <class T>
class Friction {
public:
	enum UpdateType { FIRST, ANY };

protected:
	T m_F;

public:
	virtual ~Friction(void) { NO_OP; };

	virtual void Update(const T &force, 
			const T &position, 
			const T &velocity, 
			Friction::UpdateType = ANY) = 0;

	virtual const T &F(void) const {
		return m_F;
	};
};

template <class T>
class StepFriction : public Friction<T> {
private:
	doublereal m_mu0;
public:
	StepFriction(const doublereal& mu0) 
	: m_mu0(mu0) {
		NO_OP;
	};
	
	virtual ~StepFriction(void) { NO_OP; };

	virtual void Update(const T &force,
			const T & /* position */, 
			const T & velocity, 
			typename Friction<T>::UpdateType = ANY) {
		m_F = dir(velocity)*(-norm(force)*m_mu0);
	};
};

template <class T>
class TanhFriction : public Friction<T> {
private:
	doublereal m_mu0;
	doublereal m_vRef;
public:
	TanhFriction(const doublereal& mu0, const doublereal vRef) 
	: m_mu0(mu0), m_vRef(vRef) {
		NO_OP;
	};
	
	virtual ~TanhFriction(void) { NO_OP; };

	virtual void Update(const T &force,
			const T & /* position */, 
			const T &velocity, 
			typename Friction<T>::UpdateType = ANY) {
#if 0
		m_F = dir(velocity)*(-norm(force)*m_mu0*tanh(norm(velocity)/m_vRef));
#else
		m_F = tanh(velocity/m_vRef)*(-norm(force)*m_mu0);
#endif
	};
};

template <class T>
class DiscStateFriction : public Friction<T> {
public:
	enum State {
		Stick,
		Slip
	};

private:
	State		m_state;
	T		m_maxForce;
	doublereal	m_stiffness;
	doublereal	m_damping;
	T		m_s0;
	doublereal 	m_velTreshold;

public:
	DiscStateFriction(const doublereal &maxForce, 
			const doublereal &stiffness, 
			const doublereal &damping, 
			const doublereal &velTreshold, 
			DiscStateFriction::State initialState = Stick)
	: m_state(initialState),
	m_maxForce(maxForce),
	m_stiffness(stiffness),
	m_damping(damping),
	m_s0(0.),
	m_velTreshold(velTreshold) {
		NO_OP;
	};

	virtual ~DiscStateFriction(void) { NO_OP; };

private:
	const char *S() const {
		if (m_state == Stick) {
			return "stick";
		}
		return "slip";
	};

	virtual void Update(const T &force,
			const T & position, 
			const T &velocity, 
			typename Friction<T>::UpdateType update = ANY) {

		// std::cerr << ">>> state: " << S() << "; p=" << position << "; v=" << velocity << "; m_s0=" << m_s0 << "; F=" << m_F << std::endl;
		
		switch (m_state) {
		case Stick:
			m_F = (m_s0-position)*m_stiffness;
			if (position > m_s0) {
				if (m_F < -m_maxForce) {
					m_F = -m_maxForce;
					if (update == FIRST) {
						m_state = Slip;
					}
				}
			} else {
				if (m_F > m_maxForce) {
					m_F = m_maxForce;
					if (update == FIRST) {
						m_state = Slip;
					}
				}
			}

			m_F -= velocity*m_damping;

			break;

		case Slip: {
			doublereal vs = dir(velocity);
			doublereal vn = norm(velocity);
			
			if (vn < m_velTreshold) {
				if (update == FIRST) {
					m_s0 = position + m_maxForce*(vs/m_stiffness);
					m_state = Stick;
				}
				m_F = m_maxForce*(-vs*sin(M_PI_2*vn/m_velTreshold));

			} else {
				m_F = m_maxForce*(-vs);
			}

			break;
		}
		}

		// std::cerr << "<<< state: " << S() << "; p=" << position << "; v=" << velocity << "; m_s0=" << m_s0 << "; F=" << m_F << std::endl;
	};
		
};

#endif /* FRICTION_H */

