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

#include <vector>
#include <iostream>

using namespace std;

extern "C" {
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

void chaco_interface(const int iTotVertices, int *start, int *adjacency,
                     int *vertex_weights, float *edge_weights,
                     const int num_processors, int *pParAmgProcs)
{
#ifdef USE_CHACO
  vector<short> set_assignment(iTotVertices);

  const int architecture(1),      /* The computers are connected as a
                                     simple one-dimensional mesh */
    global_method(2),             /* Use Spectral decomposition */
    local_method(1),              /* With KL local refinement */
    rqi_flag(0),                  /* Use Lanczos for solving the
                                     eigenvalue problem for spectral
                                     decomposition */
    vmax(100),                    /* When to stop coarsening.  Not used */
    seed(7654321);                /* A random number seed */
  int mesh_dims[]={1,1,1};

  double eigtol(1e-3);

  
  if(interface(iTotVertices,start,adjacency,vertex_weights,edge_weights,
               NULL,NULL,NULL,NULL,NULL,&set_assignment[0],architecture,
               num_processors,mesh_dims,NULL,global_method,local_method,rqi_flag,
               vmax,1,eigtol,seed)!=0)
    {
      cerr << "Partition failed" << endl;
    }

  for(int i=0;i<iTotVertices;++i)
    {
      pParAmgProcs[i]=set_assignment[i];
    }
#endif /* USE_CHACO */
}
