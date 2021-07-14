/* $Header$ */
/* A simple interface to the Chaco graph partioning library.  It does
 * some small conversions so that it is easy to substitute for Metis.
 *
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003 Mississippi State University
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <vector>
#include <iostream>

#include "chacowrap.h"

extern "C" {
/*
 * NOTE: the name "interface" has been reported (by Patrick Rix) 
 * to cause issues with gcc 4.10
 */
int       interface(
		int       nvtxs,		/* number of vertices in full graph */
		int      *start,		/* start of edge list for each vertex */
		int      *adjacency,		/* edge list data */
		int      *vwgts,		/* weights for all vertices */
		float    *ewgts,		/* weights for all edges */
		float    *x, float *y, float *z,/* coordinates for inertial method */
		char     *outassignname,	/* name of assignment output file */
		char     *outfilename,		/* output file name */
		short    *assignment,		/* set number of each vtx (length n) */
		int       architecture,		/* 0 => hypercube, d => d-dimensional mesh */
		int       ndims_tot,		/* total number of cube dimensions to divide */
		int       mesh_dims[3],		/* dimensions of mesh of processors */
		double   *goal,			/* desired set sizes for each set */
		int       global_method,	/* global partitioning algorithm */
		int       local_method,		/* local partitioning algorithm */
		int       rqi_flag,		/* should I use RQI/Symmlq eigensolver? */
		int       vmax,			/* how many vertices to coarsen down to? */
		int       ndims,		/* number of eigenvectors (2^d sets) */
		double    eigtol,		/* tolerance on eigenvectors */
		long      seed			/* for random graph mutations */
	);
}

extern "C" int FREE_GRAPH;

void
chaco_interface(
		const int	iTotVertices,
		int		*start,
		int		*adjacency,
		int		*vertex_weights,
		int		*comm_weights,		/* unsupported */
		int		*edge_weights,
		const int	num_processors,
		int		*pParAmgProcs
	)
{
	std::vector<short> set_assignment(iTotVertices);

	const int architecture(1),	/* The computers are connected as a
				       	 * simple one-dimensional mesh
					 * 0 => hypercube, d => d-dimensional mesh */
		global_method(2),	/* Use Spectral decomposition */
		local_method(1),	/* With KL local refinement */
		rqi_flag(0),		/* Use Lanczos for solving the
					 * eigenvalue problem for spectral
					 * decomposition */
		vmax(100),		/* When to stop coarsening.  Not used */
		seed(7654321);		/* A random number seed */
	int	mesh_dims[] = {		/* if architecture > 0, mesh size goes here */
		num_processors, 1, 1
	};
	double eigtol(1.e-3);		/* */

	/* Chaco vertices are base 1 */
	for (int i = 0; i < start[iTotVertices]; i++) {
		adjacency[i]++;
	}

	/* weights = 0 are not allowed */
	std::vector<int> ivertex_weights(iTotVertices);
	if (vertex_weights) {
		for (int i = 0; i < iTotVertices; i++) {
			/* trying to emulate communication weights? */
			if (comm_weights) {
				vertex_weights[i] += comm_weights[i];
			}

			if (vertex_weights[i] == 0) {
				ivertex_weights[i] = 1;
			} else {
				ivertex_weights[i] = vertex_weights[i];
			}
		}
	}

	/* Chaco uses floats for the communication weights */
	std::vector<float> fedge_weights(2*start[iTotVertices]);
	if (edge_weights) {
		for (int i = 0; i < iTotVertices; i++) {
			fedge_weights[i] = edge_weights[i];
			fedge_weights[start[iTotVertices] + i] = edge_weights[i];
		}
	}

	/* NOTE: don't free() adjacency !!! (default in Chaco 2.2) */
	if (FREE_GRAPH) {
		FREE_GRAPH = 0;
	}
	
	int rc = interface(
			iTotVertices,		/* number of vertices in full graph */
			start,			/* start of edge list for each vertex */
			adjacency,		/* edge list data */
			vertex_weights ? &ivertex_weights[0] : 0,	/* weights for all vertices */
			edge_weights ? &fedge_weights[0] : 0,		/* weights for all edges */
			NULL,			/* coordinates for inertial method */
			NULL,			/* ... */
			NULL,			/* ... */
			NULL,			/* name of assignment output file */
			NULL,			/* output file name */
			&set_assignment[0],	/* set number of each vtx (length n) */
			architecture,		/* 0 => hypercube, d => d-dimensional mesh */
			num_processors,		/* total number of cube dimensions to divide */
			mesh_dims,		/* dimensions of mesh of processors */
			NULL,			/* desired set sizes for each set */
			global_method,		/* global partitioning algorithm */
			local_method,		/* local partitioning algorithm */
			rqi_flag,		/* should I use RQI/Symmlq eigensolver? */
			vmax,			/* how many vertices to coarsen down to? */
			1,			/* number of eigenvectors (2^d sets) */
			eigtol,			/* tolerance on eigenvectors */
			seed			/* for random graph mutations */
	);

	/* Chaco vertices are base 1 */
	for (int i = 0; i < start[iTotVertices]; i++) {
		adjacency[i]--;
	}

	if (rc != 0) {
		silent_cerr("chaco_interface(): partition failed" << std::endl);
	}

	for (int i = 0; i < iTotVertices; i++) {
		pParAmgProcs[i] = set_assignment[i];
    	}
}

