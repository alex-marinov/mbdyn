
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

