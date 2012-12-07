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

#include "SiconosNumerics.h"
#include <iostream>
#include "mbdyn_siconos.h"


void
mbdyn_siconos_LCP_call(int sizep, double W_NN[], double bLCP[], double Pkp1[], double wlem[], solver_parameters& solparam)
{

	int info = 0;

	LCPsolver lcpsol = solparam.solver;
	double tolerance = solparam.solvertol;
	int maxiternum = solparam.solveritermax;

	// LCP problem description:
	LinearComplementarityProblem OSNSProblem;
	OSNSProblem.size = sizep;
	OSNSProblem.q = bLCP;

	NumericsMatrix MM;
	MM.storageType = 0;
	MM.matrix0 = W_NN; 		// W_delassus
	MM.size0 = sizep;
	MM.size1 = sizep;
	OSNSProblem.M = &MM;

	SolverOptions numerics_solver_options = { 0 };

	int processed_iterations = 0.;
	solparam.resulting_error = 0.;

  	// Solving LCP problem
	// FIXME: initialize once, keep alive?
	switch (lcpsol) {
	case LEXICO_LEMKE:
		linearComplementarity_lexicolemke_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		lcp_lexicolemke(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		break;

	case RPGS:
		linearComplementarity_rpgs_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_rpgs(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;

	case QP:
		linearComplementarity_qp_setDefaultSolverOptions(&numerics_solver_options);
//		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_qp(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		break;

	case CPG:
		linearComplementarity_cpg_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_cpg(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;

	case PGS:
		linearComplementarity_pgs_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_pgs(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;

	case PSOR:
		linearComplementarity_psor_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_psor(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;

	case NSQP:
		linearComplementarity_nsqp_setDefaultSolverOptions(&numerics_solver_options);
//		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_nsqp(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		break;

	case LATIN:
		linearComplementarity_latin_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_latin(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		break;

	case LATIN_W:
		linearComplementarity_latin_w_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_latin_w(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		break;

	case NEWTON_MIN:
		linearComplementarity_newton_min_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_newton_min(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;

	case NEWTON_FB:
		linearComplementarity_newton_FB_setDefaultSolverOptions(&numerics_solver_options);
		numerics_solver_options.iparam[0] = maxiternum;
		numerics_solver_options.dparam[0] = tolerance;
		lcp_newton_FB(&OSNSProblem, Pkp1, wlem, &info, &numerics_solver_options);
		processed_iterations = numerics_solver_options.iparam[1];
		solparam.resulting_error = numerics_solver_options.dparam[1];
		break;
	}

	// analyzing the LCP solver final status
	switch (info) {
	case 0:	// convergence
		break;

	case 1: // iter=itermax
//		std::cout << std::endl
//			<< "loadable element nonsmooth node: max iterations reached in LCP solver "
//			<< std::endl;
//		std::cout << "processed_iterations "<< processed_iterations
//			<< " resulting_error " << solparam.resulting_error <<std::endl;
		break;

	default: // other problem
//		std::cout << std::endl
//			<<"loadable element nonsmooth node: problem in solution of LCP"
//			<< std::endl;
//		std::cout << "processed_iterations "<< processed_iterations
//			<< " resulting_error " << solparam.resulting_error <<std::endl;
		break;
	}

	solparam.info = info;

	deleteSolverOptions(&numerics_solver_options);
}


