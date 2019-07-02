# $Header$
#
# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2017
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

#===============================================================================
# Title:	Macro for MBDyn modal joint data creation with Code Aster
# Author:	Pierangelo Masarati <pierangelo.masarati at polimi.it>
# Date:		July-August 2008, during a visit to EDF (LaMSID at Clamart, F)
#
# This file is a template for the generation of MBDyn's modal joint data
# using Code Aster <http://www.code-aster.org/>, the free FEM solver
# developed by EDF.
# It was inspired by the work done by Giulio Romanelli and Elisa Serioli
# for their graduation thesis.  Lots of material was copied and modified
# from the distributed documentation and examples.  Portions come from
# examples written by J-L. Flejou and O. Boiteau.
# Many thanks to Michael Abbas and Sylvain Mazet for their support;
# thanks also to Patrick Massin for making the visit and the cooperation
# possible.

#===============================================================================
# LOG:
# 2016-03-12: Support for Code_Aster 12.5.0 added by r.resch@a1.net
# 2008-08-30: fix TOUT = 'OUI' using 'ENV_SPHERE' (needs work)
# 2008-08-25: release with MBDyn 1.3.4-Beta
# 2008-08-05: optionally prints rigid-body inertia matrix
# 2008-08-05: optionally prints diagonal of mass matrix
# 2008-08-04: first attempt to add diagonal of mass matrix
# 2008-07-31: turned into macro
# 2008-07-29: file creation works
# 2008-07-29: normal modes + Craig-Bampton works
#
#===============================================================================
# TODO:
# - check the possibility to add inertia relief to static shapes
# - modal analysis where rigid modes are either rejected
#   or (selectively) accepted (e.g. aileron example)

#===============================================================================
#===============================================================================
#===============================================================================
# Macro for MBDyn modal element data generation
#===============================================================================
#===============================================================================
#===============================================================================

#===============================================================================
# import mesh manipulation code
# NOTE: probably this stuff is not yet available in distributed code...
from Utilitai.partition import *
from Cata.cata import *
from numpy import *
from math import *

#===============================================================================
# NOTE: this function is basically extracted from example sdll123a,
#	authored by O.BOITEAU.  I'm afraid writing anything like that
#	requires an uncommon knowledge of Code Aster's internals.

# extract the diagonal of the mass matrix for each node listed in exposed_id
def cms_diag_mass(matrrr, exposed_id):
	# construction des vecteurs jeveux
	nommatr = matrrr.nom;
	lenm = len(nommatr);
	nommatr = nommatr + ' '*(8 - lenm + 1);
	vectrav = nommatr + '          .REFA        ';
	nom = aster.getvectjev(vectrav);
	nomnume = nom[1];
	typm = nom[8];
	assert(typm[0:2] == 'MS');
	lenm = len(nomnume);
	nomnume = nomnume[0:9];

	nvar = nommatr + '          .VALM';
	nadia = nomnume + '     .SMOS.SMDI        ';
	nrtt = nomnume + '     .NUME.DELG        ';
	nrtt2 = nomnume + '     .NUME.DEEQ        ';

	var = aster.getcolljev(nvar);
	adia = aster.getvectjev(nadia);
	rtt = aster.getvectjev(nrtt);
	rtt2 = aster.getvectjev(nrtt2);

	"""
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "typm = %s" % typm;
	print "var:";
	print var;
	print "####################################################";
	print "adia:";
	print adia;
	print "####################################################";
	print "rtt:";
	print rtt;
	print "####################################################";
	print "rtt2:";
	print rtt2;
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "####################################################";
	"""

	valr = var[1];

	vc = len(rtt);

	nexposed = len(exposed_id);
	nnodes = 0;
	for ii in range(vc):
		if (rtt2[2*ii] > 0) and (rtt2[2*ii+1] == 1):
			nnodes = nnodes + 1;

	nodes_id = zeros([nnodes, 6], int);
	for ii in range(vc):
		if (rtt2[2*ii] > 0) and (rtt2[2*ii+1] > 0):
			nodes_id[rtt2[2*ii] - 1][rtt2[2*ii+1] - 1] = adia[ii];

	diagm = zeros([nexposed, 6], double);
	for rr in range(nexposed):
		assert(exposed_id[rr] <= nnodes);
		for cc in range(6):
			ii = nodes_id[exposed_id[rr]][cc];
			if (ii > 0):
				diagm[rr][cc] = valr[ii - 1];

	return diagm;
# end of cms_diag_mass

# compute the rigid-body inertia matrix
def cms_rigb_mass(matrrr, coord):
	# construction des vecteurs jeveux
	nommatr = matrrr.nom;
	lenm = len(nommatr);
	nommatr = nommatr + ' '*(8 - lenm + 1);
	vectrav = nommatr + '          .REFA        ';
	nom = aster.getvectjev(vectrav);
	nomnume = nom[1];
	typm = nom[8];
	assert(typm[0:2] == 'MS');
	lenm = len(nomnume);
	nomnume = nomnume[0:9];

	nvar = nommatr + '          .VALM';
	nadia = nomnume + '     .SMOS.SMDI        ';
	nnuml = nomnume + '     .SMOS.SMHC        ';
	nrtt = nomnume + '     .NUME.DELG        ';
	nrtt2 = nomnume + '     .NUME.DEEQ        ';

	var = aster.getcolljev(nvar);
	adia = aster.getvectjev(nadia);
	numl = aster.getvectjev(nnuml);
	rtt = aster.getvectjev(nrtt);
	rtt2 = aster.getvectjev(nrtt2);

	"""
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "typm = %s" % typm;
	print "var:";
	print var;
	print "####################################################";
	print "adia:";
	print adia;
	print "####################################################";
	print "numl:";
	print numl;
	print "####################################################";
	print "rtt:";
	print rtt;
	print "####################################################";
	print "rtt2:";
	print rtt2;
	print "####################################################";
	print "####################################################";
	print "####################################################";
	print "####################################################";
	"""

	# array containing mass matrix values
	valr = var[1];

	# number of ddl
	vc = len(rtt);

	# make room for rigid body mass matrix
	rigbm = zeros([6, 6], double);

	# make room for nodal rigid body motion matrices
	Zi = zeros([6, 6], double);
	Zj = zeros([6, 6], double);

	# initialize the diagonal
	for kk in range(6):
		Zi[kk][kk] = 1.;
		Zj[kk][kk] = 1.;

	ni_cur = -1;
	nj_cur = -1;

	for ii in range(vc):
		# skip Lagrange multipliers

		# node number (1 to nnodes)
		ni = rtt2[2*ii];
		if (ni <= 0):
			continue;

		# component number (1 to 6)
		ci = rtt2[2*ii + 1];
		if (ci <= 0):
			continue;

		# deal with self contribution to matrix
		if (ni_cur != ni):
			ni_cur = ni;
			xi = coord[ni - 1][0];
			yi = coord[ni - 1][1];
			zi = coord[ni - 1][2];

			Zi[0][4] = -zi;
			Zi[0][5] = yi;
			Zi[1][3] = zi;
			Zi[1][5] = -xi;
			Zi[2][3] = -yi;
			Zi[2][4] = xi;

		# index of diagonal entry of this ddl (base 1)
		ai = adia[ii];

		# contribution of this ddl
		dm = valr[ai - 1];

		# rigbm = rigbm + ZiT*M*Zi;
		for rr in range(6):
			for cc in range(6):
				rigbm[rr][cc] = rigbm[rr][cc] + dm*Zi[ci - 1][rr]*Zi[ci - 1][cc];

		# first ddl only
		if (ii == 0):
			continue;

		# index of diagonal entry of previous ddl (base 1)
		aim1 = adia[ii - 1];

		# array of indexes connected to ai
		# NOTE: use aim1 as array lower bound (included) because arrays are 0-based, while aim1 is 1-based)
		# NOTE: skip numl[ai] because we already considered it (it's the diagonal term) (upper bound is excluded)
		air = numl[aim1:ai - 1];
		lair = len(air);

		# check whether any of the previous ddls is connected to ai
		for jj in range(ii):
			# skip Lagrange multipliers

			# node number (1 to nnodes)
			nj = rtt2[2*jj];
			if (nj <= 0):
				continue;

			# component number (1 to 6)
			cj = rtt2[2*jj + 1];
			if (cj <= 0):
				continue;

			# index of diagonal entry of this ddl (base 1)
			aj = adia[jj];

			# ignore if this ddl (jj + 1) is not connected to ai
			if ((jj + 1) not in air):
				continue;

			# find the index of the connecting coefficient
			gotit = 0;
			for kk in range(lair):
				if (air[kk] == (jj + 1)):
					gotit = 1;
					break;

			assert(gotit == 1);

			# contribution of this ddl coupled with ii
			dm = valr[aim1 + kk];

			# deal with cross-contribution to matrix
			if (nj_cur != nj):
				nj_cur = nj;
				xj = coord[nj - 1][0];
				yj = coord[nj - 1][1];
				zj = coord[nj - 1][2];

				Zj[0][4] = -zj;
				Zj[0][5] = yj;
				Zj[1][3] = zj;
				Zj[1][5] = -xj;
				Zj[2][3] = -yj;
				Zj[2][4] = xj;

			# rigbm = rigbm + ZiT*M*Zj + ZjT*MT*Zi;
			for rr in range(6):
				for cc in range(6):
					rigbm[rr][cc] = rigbm[rr][cc] + \
						dm*(Zj[cj - 1][rr]*Zi[ci - 1][cc] + Zi[ci - 1][rr]*Zj[cj - 1][cc]);

	return rigbm;
# end of cms_rigb_mass

def cms_extract_mode_shape(_sd_tab, eps):
	# extract normal modes or static mode shapes
	# NOTE: each array is actually a matrix, containing the number of the mode
	#	in the first column, and the value in the second column; for example,
	#	for a three nodes model:
	#	sd_dx = [ 0, dx_mode0_node0 ]
	#		[ 0, dx_mode0_node1 ]
	#		[ 0, dx_mode0_node2 ]
	#		[ 1, dx_mode1_node0 ]
	#		[ 1, dx_mode1_node1 ]
	#		[ 1, dx_mode1_node2 ]
	#		[ ... ]
	#	couldn't do anything more synthetic...

        dof_names = ('DX','DY','DZ','DRX','DRY','DRZ');
        
        for j in range(len(dof_names)):
                sd=_sd_tab.EXTR_TABLE().Array('NUME_ORDRE', dof_names[j]);
                
                if j == 0:
                        modes=zeros([sd.shape[0], len(dof_names)], double);
                        
                for i in range(modes.shape[0]):
                        val = sd[i, 1];

                        # This might be helpful in case of modal shape values close to zero at the modal node
                        if (abs(val) < eps):
                                val = 0.;

                        # If the node does not have rotational degrees of freedom we get nan
                        if (isnan(val)):
                                val = 0.;
                                
                        modes[i, j] = val;
        
        return modes;

# write CMS data in MBDyn format
def cms_write_mbdyn(data, maillage, cms_interface, cms_exposed_fact, \
		nshapes, ndynamic, sol_dynamic, nstatic, sol_static, \
		macm, mack, mmass):

	from Accas import _F

	# desired precision
	precision = data[1];

	diag_mass = (data[2] == 'OUI');
	rigb_mass = (data[3] == 'OUI');

        eps = data[4];
        
	# "cook" format specifiers based on precision
	IFMT = "%" + str(precision) + "d";
	RFMT = "%" + str(precision + 8) + "." + str(precision) + "e";

	# open output file
	# NOTE: overwrites existing files
	# NOTE: must be absolute path, otherwise I don't know
	# where it shows up...
	# if (data[0][0] != '/'):
	# 	print "***";
	# 	print "*** WARNING: MBDyn modal element data file needs absolute path"
	# 	print "***          file=\"" + data[0] + "\"";
	# 	print "***";
	outf = file(data[0], 'w');

	if (cms_exposed_fact['GROUP_NO'] != None):
		cms_exposed = cms_exposed_fact['GROUP_NO'];

	else:
		assert(cms_exposed_fact['TOUT'] == 'OUI');
		# NOTE: this name must be reserved to MBDyn (!?!)
		cms_exposed = 'MBDYN_TN';
		# NOTE: 1.e+38 to make sure we catch all (?!?)
		_ma = maillage

		# NOTE: hack to find out whether a mesh is 2D or 3D
		# create handler for mesh
		mm = MAIL_PY();
		mm.FromAster(maillage);

		# Coordonnees des noeuds
		coord = mm.cn;

		is_3D = 0;
		z0 = coord[0][2];
		for nn in range(1, mm.dime_maillage[0] - 1):
			if (coord[nn][2] != z0):
				is_3D = 1;
				break;

		if is_3D:
			# FIXME: only works if model is truly 3D :-(
			# TODO: test whether the mesh is 2D or z_cst
			_ma = DEFI_GROUP(reuse = _ma,
					 MAILLAGE = _ma,
						CREA_GROUP_NO = ( _F(
							NOM = cms_exposed,
							OPTION = 'ENV_SPHERE',
							POINT = ( 0.0, 0.0, 0.0 ),
							RAYON = 1.e+38,
							PRECISION = ( 1.e+38 ) ) ) );
		else:
			_ma = DEFI_GROUP(reuse = _ma,
					 MAILLAGE = _ma,
						CREA_GROUP_NO = ( _F(
							NOM = cms_exposed,
							OPTION = 'ENV_SPHERE',
							POINT = ( 0.0, 0.0 ),
							RAYON = 1.e+38,
							PRECISION = ( 1.e+38 ) ) ) );
		maillage = _ma;

	# create handler for mesh
	mm = MAIL_PY();
	mm.FromAster(maillage);

	# Coordonnees des noeuds
	coord = mm.cn;
	# Noms des noeuds
	linomno = list(mm.correspondance_noeuds);

	# Groupe de Noeuds
	exposed_id = mm.gno[cms_exposed];
	nexposed = len(exposed_id);

	# header
	outf.write("** MBDyn MODAL DATA FILE - generated by Aster\n");
	outf.write("** NODE SET '" + cms_exposed + "'\n");

	# record 1
	outf.write("** RECORD GROUP 1, HEADER\n");
	outf.write("**   REVISION,  NODE,  NORMAL, ATTACHMENT, CONSTRAINT, REJECTED MODES.\n");
	outf.write((" REV0       " + IFMT + "  " + IFMT + "  " + IFMT + "  " + IFMT + "  " + IFMT + "\n") \
		% (nexposed, ndynamic, nstatic, 0, 0) );
	outf.write("**\n");

	# record 2
	outf.write("** RECORD GROUP 2, FINITE ELEMENT NODE LIST\n")
	l = 0
	while l < nexposed:
		if ((l > 0) & ((l % 6) == 0)):
			outf.write("\n");
		outf.write(" " + linomno[exposed_id[l]]);
		l = l + 1;
	outf.write("\n");
	outf.write("**\n");

	# record 3 (optional, default to zero, so could be omitted)
	outf.write("** RECORD GROUP 3, INITIAL MODAL DISPLACEMENTS\n");
	for m in range(nshapes):
		if ((m > 0) & ((m % 6) == 0)):
			outf.write("\n");
		outf.write(RFMT % 0.0);
	outf.write("\n");
	outf.write("**\n");

	# record 4 (optional, default to zero, so could be omitted)
	outf.write("** RECORD GROUP 4, INITIAL MODAL VELOCITIES\n");
	for m in range(nshapes):
		if ((m > 0) & ((m % 6) == 0)):
			outf.write("\n");
		outf.write(RFMT % 0.0);
	outf.write("\n");
	outf.write("**\n");

	# record 5
	outf.write("** RECORD GROUP 5, NODAL X COORDINATES\n");
	n = 0;
	while n < nexposed:
		outf.write((RFMT + "\n") % coord[exposed_id[n]][0]);
		n = n + 1;
	outf.write("**\n");

	# record 6
	outf.write("** RECORD GROUP 6, NODAL Y COORDINATES\n");
	n = 0;
	while n < nexposed:
		outf.write((RFMT + "\n") % coord[exposed_id[n]][1]);
		n = n + 1;
	outf.write("**\n");

	# record 7
	outf.write("** RECORD GROUP 7, NODAL Z COORDINATES\n");
	n = 0;
	while n < nexposed:
		outf.write((RFMT + "\n") % coord[exposed_id[n]][2]);
		n = n + 1;
	outf.write("**\n");

	# record 8
	outf.write("** RECORD GROUP 8, MODE SHAPES\n");

	for m in range(ndynamic):
                idx = 0;
                _sd_tab = POST_RELEVE_T( ACTION = _F(INTITULE = 'Normal modes',
						GROUP_NO = cms_exposed,
						RESULTAT = sol_dynamic,
						NOM_CHAM = 'DEPL',
					             NUME_ORDRE = (m + 1,),
						TOUT_CMP = 'OUI',
						OPERATION = 'EXTRACTION' ) );

                # extract normal modes
                modes_sd = cms_extract_mode_shape(_sd_tab, eps);
                
		outf.write("**    NORMAL MODE SHAPE #  %d\n" % (m + 1));
		n = 0;
		while n < nexposed:
                        for j in range(modes_sd.shape[1]):
                                outf.write(RFMT % modes_sd[idx, j]);
                                
                        outf.write("\n");
                        
			n = n + 1;
			idx = idx + 1;

                del _sd_tab;
                del modes_sd;


	if nstatic > 0:
		for m in range(nstatic):
                        _ss_tab = POST_RELEVE_T(ACTION = _F(INTITULE = 'Static Shapes',
							GROUP_NO = cms_exposed,
							RESULTAT = sol_static,
                                                            NUME_ORDRE = (m + 1,),
							NOM_CHAM = 'DEPL',
							TOUT_CMP = 'OUI',
							OPERATION = 'EXTRACTION' ) );

                        # extract static shapes
                        modes_ss = cms_extract_mode_shape(_ss_tab, eps);

		idx = 0;
			outf.write("**    NORMAL MODE SHAPE #  %d (STATIC SHAPE #  %d)\n" % (ndynamic + m + 1, m + 1));
			n = 0;
			while n < nexposed:
                                for j in range(modes_ss.shape[1]):
                                        outf.write(RFMT % modes_ss[idx, j]);
                                
                                outf.write("\n");
                        
				n = n + 1;
				idx = idx + 1;

                        del _ss_tab;
                        del modes_ss;

	outf.write("**\n");

	# record 9
	outf.write("** RECORD GROUP 9, MODAL MASS MATRIX\n");
	for r in range(nshapes):
		for c in range(nshapes):
			outf.write(RFMT % macm[r, c]);
		outf.write("\n");
	outf.write("**\n");

	# record 10
	outf.write("** RECORD GROUP 10, MODAL STIFFNESS MATRIX\n");
	for r in range(nshapes):
		for c in range(nshapes):
			outf.write(RFMT % mack[r, c]);
		outf.write("\n");
	outf.write("**\n");

	# NOTE: records 11 and 12 used to be mutually exclusive;
	#	anyway, this macro allows to set both
	if diag_mass:
		# record 11
		outf.write("** RECORD GROUP 11, DIAGONAL OF LUMPED MASS MATRIX\n");
		diagm = cms_diag_mass(mmass, exposed_id);

		for n in range(nexposed):
			for c in range(6):
				outf.write(RFMT % diagm[n][c]);
			outf.write("\n");
		outf.write("**\n");

	if rigb_mass:
		# record 12
		outf.write("** RECORD GROUP 12, RIGID BODY INERTIA MATRIX\n");
		rigbm = cms_rigb_mass(mmass, coord);

		m = rigbm[0][0];
		Xcm = zeros(3, double);
		J = zeros([3, 3], double);
		if (m > 0.):
			Xcm[0] = rigbm[2][4]/m;
			Xcm[1] = rigbm[0][5]/m;
			Xcm[2] = rigbm[1][3]/m;

			JJ = m*(Xcm[0]*Xcm[0] + Xcm[1]*Xcm[1] + Xcm[2]*Xcm[2]);

			for rr in range(3):
				for cc in range(3):
					J[rr][cc] = rigbm[3 + rr][3 + cc] + m*Xcm[rr]*Xcm[cc];
				J[rr][rr] = J[rr][rr] - JJ;

		outf.write((RFMT + "\n") % m);
		outf.write((RFMT + RFMT + RFMT + "\n") % (Xcm[0], Xcm[1], Xcm[2]));
		outf.write((RFMT + RFMT + RFMT + "\n") % (J[0][0], J[0][1], J[0][2]));
		outf.write((RFMT + RFMT + RFMT + "\n") % (J[1][0], J[1][1], J[1][2]));
		outf.write((RFMT + RFMT + RFMT + "\n") % (J[2][0], J[2][1], J[2][2]));
		outf.write("**\n");
	outf.close();

	return 0;
# end of cms_write_mbdyn

# write CMS data
def cms_write(gen_model, mmass, maillage, cms_interface, cms_exposed_fact, sol_dynamic, sol_static, type, data, options):
	rc = -1;
	if (type == 'MBDYN'):
		#===============================================================================
		# extract generalized problem info and generalized matrices
		ndynamic = options['NMAX_FREQ'];
		gen_k = gen_model.EXTR_MATR_GENE( 'RIGI_GENE' );
		gen_m = gen_model.EXTR_MATR_GENE( 'MASS_GENE' );
                nshapes = gen_k.shape[0];
		nstatic = nshapes - ndynamic;

		rc = cms_write_mbdyn(data, maillage, cms_interface, cms_exposed_fact, \
			nshapes, ndynamic, sol_dynamic, nstatic, sol_static, \
			gen_m, gen_k, mmass);

	# end of type == 'MBDYN'
	return rc;

# end of cms_write

# compute CMS data
def cms_ops(self, MAILLAGE, INTERFACE, EXPOSED, MODELE, CARA_ELEM, CHAM_MATER, CHAR_MECA, OPTIONS, OUT, **args):

	from Accas import _F

	# Define output and initialize error counter
	self.set_icmd(1)
	self.DeclareOut('GENMOD', self.sd)


	maillage = MAILLAGE;
	cms_interface = INTERFACE;
	cms_exposed_fact = EXPOSED;
	model = MODELE;
	caele = CARA_ELEM;
	chmat = CHAM_MATER;
	chbc = CHAR_MECA;
	options = OPTIONS;
	out = OUT;

	if options != None:
		cms_nmax_freq = options['NMAX_FREQ'];

	if out['TYPE'] == 'MBDYN':
		type = 'MBDYN';
		fichier = "fort.%d" % out['UNITE'];
		precision = out['PRECISION'];
		diag_mass = out['DIAG_MASS'];
		rigb_mass = out['RIGB_MASS'];
                epsilon = out['MODE_SHAPE_EPSILON'];
		data = (fichier, precision, diag_mass, rigb_mass, epsilon);
	else:
		ier = 1;
		return ier;

	# if an interface exists, use it
	do_craig_bampton = 0;
	if cms_interface != None:
		do_craig_bampton = 1;	# the real one
		# do_craig_bampton = 0;	# for debugging

	if do_craig_bampton:
		# clamp all interface nodes
		# NOTE: the model might need to have other parts clamped
		_cms_chbc = AFFE_CHAR_MECA(MODELE = model,
						DDL_IMPO = ( _F( GROUP_NO = cms_interface,
								DX = 0., DY = 0., DZ = 0.,
								DRX = 0., DRY = 0., DRZ = 0. ) ) );

	#===============================================================================
	# compute element contributions to matrices, number equations
	# and prepare stiffness and mass matrices
	# NOTE: specific to CMS, adds only boundary conditions required
	# 	for Craig-Bampton
	if chbc == None:
		if do_craig_bampton:
			_matlock = CALC_MATR_ELEM(MODELE = model,
							CARA_ELEM = caele,
							CHAM_MATER = chmat,
							CHARGE = _cms_chbc,
							OPTION = 'RIGI_MECA' );
		else:
			_matlock = CALC_MATR_ELEM(MODELE = model,
							CARA_ELEM = caele,
							CHAM_MATER = chmat,
							OPTION = 'RIGI_MECA' );

	else:
		if do_craig_bampton:
			_matlock = CALC_MATR_ELEM(MODELE = model,
							CARA_ELEM = caele,
							CHAM_MATER = chmat,
							CHARGE = ( chbc, _cms_chbc),
							OPTION = 'RIGI_MECA' );
		else:
			_matlock = CALC_MATR_ELEM(MODELE = model,
							CARA_ELEM = caele,
							CHAM_MATER = chmat,
							CHARGE = chbc,
							OPTION = 'RIGI_MECA' );

	_matlocm = CALC_MATR_ELEM(MODELE = model,
					CARA_ELEM = caele,
					CHAM_MATER = chmat,
					OPTION = 'MASS_MECA' );

	_matmdiag = CALC_MATR_ELEM(MODELE = model,
				   CARA_ELEM = caele,
				   CHAM_MATER = chmat,
				   OPTION = 'MASS_MECA_DIAG' );
        
	_num = NUME_DDL(MATR_RIGI = _matlock);
        
	_matassk = ASSE_MATRICE(MATR_ELEM = _matlock,
	                        NUME_DDL = _num);
        
	_matassm = ASSE_MATRICE(MATR_ELEM = _matlocm,
	                        NUME_DDL = _num);

        _matassmdiag = ASSE_MATRICE(MATR_ELEM = _matmdiag,
	                            NUME_DDL = _num);

        _sol_dyn=CALC_MODES(TYPE_RESU='DYNAMIQUE',
                            OPTION='PLUS_PETITE',
                            SOLVEUR_MODAL=_F(METHODE='TRI_DIAG',
                                             DIM_SOUS_ESPACE=5*cms_nmax_freq,),
                            MATR_RIGI=_matassk,
                            MATR_MASS=_matassm,
                            CALC_FREQ=_F(NMAX_FREQ=cms_nmax_freq,),
                            NORM_MODE=_F(NORME='MASS_GENE',),
                            VERI_MODE=_F(STOP_ERREUR='NON',),);


	#===============================================================================
	# put together the modal and static solutions
	_sol_stat = 0;
	if do_craig_bampton:
		_sol_stat = MODE_STATIQUE(MATR_RIGI = _matassk,
						MATR_MASS = _matassm,
						MODE_STAT = _F( GROUP_NO = cms_interface,
							TOUT_CMP = 'OUI' ),
						INFO = 2 );

		_interfa = DEFI_INTERF_DYNA(NUME_DDL = _num,
						INTERFACE = _F(	NOM = 'INTERF',
								TYPE = 'CRAIGB',
								GROUP_NO = cms_interface ) );

		# NOTE: RITZ no longer allows MODE_STAT in CodeAster 10.X
		_bm = DEFI_BASE_MODALE(RITZ=(_F(MODE_MECA = _sol_dyn ),
					     _F( MODE_INTF = _sol_stat)),
				       INTERF_DYNA = _interfa,
				       NUME_REF = _num);

	else:
		_bm = DEFI_BASE_MODALE(RITZ = _F( MODE_MECA = _sol_dyn),
				       NUME_REF = _num );

	#===============================================================================
	# create macro-element
	GENMOD = MACR_ELEM_DYNA(BASE_MODALE = _bm,
				MATR_RIGI = _matassk,
				MATR_MASS = _matassm );

	rc = cms_write(GENMOD, _matassmdiag, maillage, cms_interface, cms_exposed_fact, _sol_dyn, _sol_stat, type, data, options);

	return rc;
# end of cms_gen

# ===============================================================================
# MACRO CATALOGUE DEFINITION
# ===============================================================================
CMS = MACRO(	nom		= "CMS",
		op		= cms_ops,
		sd_prod		= macr_elem_dyna,
		fr		= "Create Component Mode Synthesis (CMS) data for an external program",
		MAILLAGE	= SIMP( statut = 'o', typ = maillage_sdaster,
					fr = "The maillage (mesh)" ),
		INTERFACE	= SIMP( statut = 'd', typ = 'TXM', defaut=None,
					fr = "Name of the interface nodes group" ),
		EXPOSED		= FACT(statut = 'o', regles = ( UN_PARMI( 'TOUT', 'GROUP_NO') ),
			TOUT		= SIMP( statut = 'f', typ = 'TXM', into = ( 'OUI' ),
						fr = "Expose all nodes" ),
			GROUP_NO	= SIMP( statut = 'f', typ = 'TXM',
						fr = "Name of the exposed nodes group" )
		),
		MODELE		= SIMP( statut = 'o', typ = modele_sdaster,
					fr = "The model" ),
		CARA_ELEM	= SIMP( statut = 'o', typ = cara_elem,
					fr = "The elements" ),
		CHAM_MATER	= SIMP( statut = 'o', typ = cham_mater,
					fr = "The materials" ),
		CHAR_MECA	= SIMP( statut = 'd', typ = char_meca, defaut=None,
					fr = "Extra boundary conditions" ),
		OPTIONS		= FACT( statut = 'd', defaut = None,
			NMAX_FREQ	= SIMP( statut = 'd', typ = 'I', defaut = 15,
						fr = "Maximum number of frequencies" ),
		),
		OUT		= FACT(
			TYPE		= SIMP( statut = 'o', typ = 'TXM',
						into = ( "MBDYN" ), fr = "Type of output" ),
			b_data		= BLOC(condition = "TYPE == 'MBDYN'",
						UNITE = SIMP(statut = 'd', typ = 'I',
							fr = "unit number for output file" ),
						PRECISION = SIMP(statut = 'd', typ = 'I',
							defaut = 8,
							fr = "Output precision (number of significant digits)" ),
                                                MODE_SHAPE_EPSILON = SIMP(statut = 'd', typ = 'R',
							defaut = 0.,
							fr = "eigenvector components below this value will be set to zero" ),
						DIAG_MASS = SIMP(statut = 'd', typ = 'TXM',
							into = ( 'OUI', 'NON' ), defaut = 'NON',
							fr = "Add diagonal of mass matrix to output; only makes sense when all nodes are exposed (not implemented yet)" ),
						RIGB_MASS = SIMP(statut = 'd', typ = 'TXM',
							into = ( 'OUI', 'NON' ), defaut = 'NON',
							fr = "Add rigid-body inertia to output" )
			)
		)
	)
# ===============================================================================
# END OF MACRO CATALOGUE DEFINITION
# ===============================================================================

