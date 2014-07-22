# $Header$
#
# MBDyn (C) is a multibody analysis code.
# http://www.mbdyn.org
#
# Copyright (C) 1996-2014
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
# Title:	Macro to generate dynamic loads from acceleration field
# Author:	Pierangelo Masarati <pierangelo.masarati at polimi.it>
# Date:		July-August 2008, during a visit to EDF (LaMSID at Clamart, F)
#
#===============================================================================
# LOG:
# 2007-08-05: first implementation
# 2007-08-06: works
#
#===============================================================================
# TODO:
# - allow other means to pass the nodal accelerations field

#===============================================================================
#===============================================================================
#===============================================================================
# Macro to generate dynamic loads from nodal accelerations field
#===============================================================================
#===============================================================================
#===============================================================================

#===============================================================================
# import mesh manipulation code
# NOTE: probably this stuff is not yet available in distributed code...
from Utilitai.partition import *
from Cata.cata import *

from Numeric import *
from LinearAlgebra import *

from Accas import _F

# read nodal accelerations field from file
def dynforces_read_acc_from_file(fichier, mm):
	# input file
	fin = file(fichier, "r");

	# prepare reverse indexing of node labels
	nodelabels = list(mm.correspondance_noeuds);
	nnodes = len(nodelabels);

	nodeidx = {};
	gotten = {};
	for nn in range(0, nnodes):
		nodeidx[nodelabels[nn].rsplit()[0]] = nn;

	acc = Numeric.zeros([nnodes, 6], Numeric.Float);

	cnt = 0;
	lcnt = 0;
	for line in fin.readlines():
		# increment line counter
		lcnt = lcnt + 1;

		# parse line
		flds = line.rsplit();
		if (len(flds) == 0 or flds[0][0] == '%'):
			continue;

		# check the right number of fields is present
		if (len(flds) != 7):
			# error?
			print "AFFE_DYNFORCES: line %d: WARNING: ignored" % lcnt;
			continue;

		# idx from maillage, based on label = flds[0]
		label = flds[0];

		# check if nodes are repeated
		if (gotten.get(label, None) is not None):
			# error?
			print "AFFE_DYNFORCES: line %d: WARNING: node %s repeated" % (lcnt, label);
			continue;
		gotten[label] = 1;

		# node index
		idx = nodeidx.get(label, None);
		if (idx is None):
			# error?
			print "AFFE_DYNFORCES: line %d: WARNING: node %s undefined (ignoring)" % (lcnt, label);
			continue;

		# read accelerations
		for cc in range(0, 6):
			# FIXME: check if number?
			acc[idx][cc] = eval(flds[cc + 1]);

		# increment node counter
		cnt = cnt + 1;

	# sanity check
	if (cnt != nnodes):
		print "AFFE_DYNFORCES: WARNING: expected %d nodes, got %d from file \"%s\"" % (nnodes, cnt, fichier);

	for rr in range(0, nnodes):
		print "node[%d]=%s %16.8e%16.8e%16.8e%16.8e%16.8e%16.8e" % (rr + 1, nodelabels[rr], acc[rr][0], acc[rr][1], acc[rr][2], acc[rr][3], acc[rr][4], acc[rr][5]);

	fin.close();

	return acc;
# end of dynforces_read_acc_from_file

# compute the dynamic loads by multiplying the mass matrix
# times the nodal accelerations field
def dynforces_acc2frc(matrm, acc):
	# acc is a matrix containing three components of acceleration 
	# for each node; nodes are sorted as in the internal structure

	# construction des vecteurs jeveux
	nommatr = matrm.nom;
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

	# array containing mass matrix values
	valr = var[1];

	# number of ddl
	vc = len(rtt);

	nnodalddl = 0;
	nnodes = 0;
	nkk_cur = -1;
	for kk in range(0, vc):
		nkk = rtt2[2*kk];
		if (nkk <= 0):
			continue;

		ckk = rtt2[2*kk + 1];
		if (ckk <= 0):
			continue;

		if (nkk_cur != nkk):
			nkk_cur = nkk;
			nnodes = nnodes + 1;

		nnodalddl = nnodalddl + 1;

	# make room for load vector
	dynforces = Numeric.zeros([nnodes, 6], Numeric.Float);

	# for each ddl
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

		# index of diagonal entry of this ddl (base 1)
		ai = adia[ii];

		# contribution of this ddl
		dm = valr[ai - 1];

		# deal with self-contribution to load vector
		dynforces[ni - 1][ci - 1] = dynforces[ni - 1][ci - 1] - dm*acc[ni - 1][ci - 1];

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

			# deal with cross-contribution to load vector
			dynforces[ni - 1][ci - 1] = dynforces[ni - 1][ci - 1] - dm*acc[nj - 1][cj - 1];
			dynforces[nj - 1][cj - 1] = dynforces[nj - 1][cj - 1] - dm*acc[ni - 1][ci - 1];

	return dynforces;
# end of dynforces_acc2frc

def dynforces_ops(self, MODELE, MATR_MASS, MAILLAGE, FICHIER):
	# Define output and initialize error counter
	self.set_icmd(1)
	self.DeclareOut('CHDF', self.sd)

	model = MODELE;
	matrm = MATR_MASS;
	maillage = MAILLAGE;
	fichier = FICHIER;

	# create handler for mesh
	mm = MAIL_PY();
	mm.FromAster(maillage);

	# Noms des noeuds
	nodelabels = list(mm.correspondance_noeuds);
	nnodes = len(nodelabels);

	# right now, only read from file
	acc = dynforces_read_acc_from_file(fichier, mm);

	# multiply mass matrix times accelerations
	dynforces = dynforces_acc2frc(matrm, acc);

	# create dynamic loads right-hand side
	CHDF = AFFE_CHAR_MECA(	MODELE = model,
				FORCE_NODALE = ( [ _F(	NOEUD = nodelabels[nn], 
							FX = dynforces[nn][0],
							FY = dynforces[nn][1],
							FZ = dynforces[nn][2],
							MX = dynforces[nn][3],
							MY = dynforces[nn][4],
							MZ = dynforces[nn][5] )
							for nn in range(0, nnodes) ] ) );
# end of dynforces_char_meca

# ===============================================================================
# MACRO CATALOGUE DEFINITION
# ===============================================================================
AFFE_DYNFORCES = MACRO(	nom		= "AFFE_DYNFORCES",
			op		= dynforces_ops,
			sd_prod		= char_meca,
			reentrant	= 'n',
			fr		= "Create inertia loads from nodal acceleration field",
			MODELE		= SIMP( statut = 'o', typ = modele_sdaster,
						fr = "The model" ),
			MATR_MASS	= SIMP( statut = 'o', typ = matr_asse_depl_r,
						fr = "The mass matrix" ),
			MAILLAGE	= SIMP( statut = 'o', typ = maillage_sdaster,
						fr = "The maillage (mesh)" ),
			FICHIER		= SIMP( statut = 'o', typ = 'TXM',
						fr = "Input file name" )
	)
# ===============================================================================
# END OF MACRO CATALOGUE DEFINITION
# ===============================================================================

