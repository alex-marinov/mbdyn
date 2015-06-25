# $Header$
# MBDyn (C) is a multibody analysis code. 
# http://www.mbdyn.org
# 
# Copyright (C) 1996-2015
# 
# Pierangelo Masarati	<masarati@aero.polimi.it>
# Paolo Mantegazza	<mantegazza@aero.polimi.it>
# 
# Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
# via La Masa, 34 - 20156 Milano, Italy
# http://www.aero.polimi.it
# 
# Changing this copyright notice is forbidden.
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation (version 2 of the License).
# 
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

import sys
import os
from numpy import *

try:
        import mbc_py
except ImportError:
        print "Import Error: mbc_py module"

class mbcNodal:
	def __init__(self, path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels):
		""" initialize the module """
		self.id = mbc_py.mbc_py_nodal_initialize(path, host, port, timeout, verbose, data_and_next, refnode, nodes, labels, rot, accels);
		if self.id < 0:
			print "mbc_py_nodal_initialize: error";
			raise Exception;

	def negotiate(self):
		""" set pointers """
		if (mbc_py.mbc_py_nodal_negotiate(self.id) < 0):
			print "mbc_py_nodal_negotiate: error";
			raise Exception;

		self.r_k_label = mbc_py.cvar.mbc_r_k_label;
		self.r_x = mbc_py.cvar.mbc_r_x;
		self.r_theta = mbc_py.cvar.mbc_r_theta;
		self.r_r = mbc_py.cvar.mbc_r_r;
		self.r_euler_123 = mbc_py.cvar.mbc_r_euler_123;
		self.r_xp = mbc_py.cvar.mbc_r_xp;
		self.r_omega = mbc_py.cvar.mbc_r_omega;
		self.r_xpp = mbc_py.cvar.mbc_r_xpp;
		self.r_omegap = mbc_py.cvar.mbc_r_omegap;
		self.r_d_label = mbc_py.cvar.mbc_r_d_label;
		self.r_f = mbc_py.cvar.mbc_r_f;
		self.r_m = mbc_py.cvar.mbc_r_m;

		self.r_k_label_size = mbc_py.cvar.mbc_r_k_label_size;
		self.r_x_size = mbc_py.cvar.mbc_r_x_size;
		self.r_theta_size = mbc_py.cvar.mbc_r_theta_size;
		self.r_r_size = mbc_py.cvar.mbc_r_r_size;
		self.r_euler_123_size = mbc_py.cvar.mbc_r_euler_123_size;
		self.r_xp_size = mbc_py.cvar.mbc_r_xp_size;
		self.r_omega_size = mbc_py.cvar.mbc_r_omega_size;
		self.r_xpp_size = mbc_py.cvar.mbc_r_xpp_size;
		self.r_omegap_size = mbc_py.cvar.mbc_r_omegap_size;
		self.r_d_label_size = mbc_py.cvar.mbc_r_d_label_size;
		self.r_f_size = mbc_py.cvar.mbc_r_f_size;
		self.r_m_size = mbc_py.cvar.mbc_r_m_size;

		self.n_k_labels = mbc_py.cvar.mbc_n_k_labels;
		self.n_x = mbc_py.cvar.mbc_n_x;
		self.n_theta = mbc_py.cvar.mbc_n_theta;
		self.n_r = mbc_py.cvar.mbc_n_r;
		self.n_euler_123 = mbc_py.cvar.mbc_n_euler_123;
		self.n_xp = mbc_py.cvar.mbc_n_xp;
		self.n_omega = mbc_py.cvar.mbc_n_omega;
		self.n_xpp = mbc_py.cvar.mbc_n_xpp;
		self.n_omegap = mbc_py.cvar.mbc_n_omegap;
		self.n_d_labels = mbc_py.cvar.mbc_n_d_labels;
		self.n_f = mbc_py.cvar.mbc_n_f;
		self.n_m = mbc_py.cvar.mbc_n_m;

		self.n_k_labels_size = mbc_py.cvar.mbc_n_k_labels_size;
		self.n_x_size = mbc_py.cvar.mbc_n_x_size;
		self.n_theta_size = mbc_py.cvar.mbc_n_theta_size;
		self.n_r_size = mbc_py.cvar.mbc_n_r_size;
		self.n_euler_123_size = mbc_py.cvar.mbc_n_euler_123_size;
		self.n_xp_size = mbc_py.cvar.mbc_n_xp_size;
		self.n_omega_size = mbc_py.cvar.mbc_n_omega_size;
		self.n_xpp_size = mbc_py.cvar.mbc_n_xpp_size;
		self.n_omegap_size = mbc_py.cvar.mbc_n_omegap_size;
		self.n_d_labels_size = mbc_py.cvar.mbc_n_d_labels_size;
		self.n_f_size = mbc_py.cvar.mbc_n_f_size;
		self.n_m_size = mbc_py.cvar.mbc_n_m_size;

	def send(self, last):
		""" send forces to peer """
		return mbc_py.mbc_py_nodal_send(self.id, last);

	def recv(self):
		""" receive kinematics from peer """
		return mbc_py.mbc_py_nodal_recv(self.id);

	def destroy(self):
		""" destroy handler """
		return mbc_py.mbc_py_nodal_destroy(self.id);

class mbcModal:
	def __init__(self, path, host, port, timeout, verbose, data_and_next, refnode, modes):
		""" initialize the module """
		self.id = mbc_py.mbc_py_modal_initialize(path, host, port, timeout, verbose, data_and_next, refnode, modes);
		if self.id < 0:
			print "mbc_py_modal_initialize: error";
			raise Exception;

	def negotiate(self):
		""" set pointers """
		if (mbc_py.mbc_py_modal_negotiate(self.id) < 0):
			print "mbc_py_modal_negotiate: error";
			raise Exception;

		self.r_k_label = mbc_py.cvar.mbc_r_k_label;
		self.r_x = mbc_py.cvar.mbc_r_x;
		self.r_theta = mbc_py.cvar.mbc_r_theta;
		self.r_r = mbc_py.cvar.mbc_r_r;
		self.r_euler_123 = mbc_py.cvar.mbc_r_euler_123;
		self.r_xp = mbc_py.cvar.mbc_r_xp;
		self.r_omega = mbc_py.cvar.mbc_r_omega;
		self.r_xpp = mbc_py.cvar.mbc_r_xpp;
		self.r_omegap = mbc_py.cvar.mbc_r_omegap;
		self.r_d_label = mbc_py.cvar.mbc_r_d_label;
		self.r_f = mbc_py.cvar.mbc_r_f;
		self.r_m = mbc_py.cvar.mbc_r_m;

		self.r_k_label_size = mbc_py.cvar.mbc_r_k_label_size;
		self.r_x_size = mbc_py.cvar.mbc_r_x_size;
		self.r_theta_size = mbc_py.cvar.mbc_r_theta_size;
		self.r_r_size = mbc_py.cvar.mbc_r_r_size;
		self.r_euler_123_size = mbc_py.cvar.mbc_r_euler_123_size;
		self.r_xp_size = mbc_py.cvar.mbc_r_xp_size;
		self.r_omega_size = mbc_py.cvar.mbc_r_omega_size;
		self.r_xpp_size = mbc_py.cvar.mbc_r_xpp_size;
		self.r_omegap_size = mbc_py.cvar.mbc_r_omegap_size;
		self.r_d_label_size = mbc_py.cvar.mbc_r_d_label_size;
		self.r_f_size = mbc_py.cvar.mbc_r_f_size;
		self.r_m_size = mbc_py.cvar.mbc_r_m_size;

		self.m_q = mbc_py.cvar.mbc_m_q;
		self.m_qp = mbc_py.cvar.mbc_m_qp;
		self.m_p = mbc_py.cvar.mbc_m_p;

		self.m_q_size = mbc_py.cvar.mbc_m_q_size;
		self.m_qp_size = mbc_py.cvar.mbc_m_qp_size;
		self.m_p_size = mbc_py.cvar.mbc_m_p_size;

	def send(self, last):
		""" send forces to peer """
		return mbc_py.mbc_py_modal_send(self.id, last);

	def recv(self):
		""" receive kinematics from peer """
		return mbc_py.mbc_py_modal_recv(self.id);

	def destroy(self):
		""" destroy handler """
		return mbc_py.mbc_py_modal_destroy(self.id);

