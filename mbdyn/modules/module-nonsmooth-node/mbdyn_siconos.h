/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2012
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
 * Author: Matteo Fancello <matteo.fancello@gmail.com>
 * Nonsmooth dynamics element;
 * uses SICONOS <http://siconos.gforge.inria.fr/>
 */

enum LCPsolver
		{	
			QP, 	//quadratic programm formulation
			CPG, 	//CPG (Conjugated Projected Gradient) solver for LCP based on quadratic minimization.
			PGS, 	//PGS is a basic Projected Gauss-Seidel solver for LCP.
			RPGS, 	//Regularized Projected Gauss-Seidel, is a solver for LCP, 
				// able to handle matrices with null diagonal terms
			PSOR, 	//Projected Succesive over relaxation solver for LCP. See cottle, Pang Stone Chap 5
			NSQP, 	//quadratic programm formulation for solving an non symmetric LCP
			LATIN, 	//(LArge Time INcrements) is a basic latin solver for LCP.
			LATIN_W, 	//(LArge Time INcrements) is a basic latin solver with relaxation for LCP
			LEXICO_LEMKE, 	//direct solver for LCP based on pivoting method principle for degenerate problem.
					//Choice of pivot variable is performed via lexicographic ordering 
			NEWTON_MIN, 	//nonsmooth Newton method based on the min formulation (or max formulation) of the LCP
			NEWTON_FB 	//uses a nonsmooth newton method based on the Fischer-Bursmeister convex function
			//NSGS_SBM 	//Gauss-Seidel solver based on a Sparse-Block storage for the matrix M of the LCP.
					//Can't be used here because Matrix M of the LCP must be formulated as SparseBlockStructuredMatrix. 
		};

struct solver_parameters
{
	// input parameters
	LCPsolver solver; 
	double solvertoll;
	int solveritermax;
	
	// output 
	int info;
	double resulting_error;		// only for: CPG, PGS, RPGS, NEWTON, LATIN, PSOR

};


extern void siconos_LCP_call(int size, double M[], double blcp[], double zlem[], double wlem[], solver_parameters * solparam);

