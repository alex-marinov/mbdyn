/* 
 * HmFe (C) is a FEM analysis code. 
 *
 * Copyright (C) 1996-2001
 *
 * Marco Morandini  <morandini@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 */
/* December 2001 
 * Modified to add a Sparse matrix in row form and to implement methods
 * to be used in the parallel MBDyn Solver.
 *
 * Copyright (C) 1996-2001
 *
 * Giuseppe Quaranta  <quaranta@aero.polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *      
 */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2000
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
#ifndef Umfpack3SparseLUSolutionManager_hh
#define Umfpack3SparseLUSolutionManager_hh

#ifdef USE_UMFPACK3

#include <vector>
extern "C" {
#include "umfpack.h"
}

#include "solman.h"
#include "spmapmh.h"


class Umfpack3SparseLUSolutionManager:public SolutionManager {
 
 private:
     mutable SpMapMatrixHandler A;
     MyVectorHandler *xVH, *bVH;
     std::vector<double> x;
     std::vector<double> b;
     std::vector<double> Ax;
     std::vector<int> Ai;
     std::vector<int> Ap;

     double t,t1;
	
     int status;

     void * Symbolic;
     double Control[UMFPACK_CONTROL];
//	Control[UMFPACK_PRL] = 6 ;
     double Info[UMFPACK_INFO];
     void * Numeric;

     
 protected:      

 public:
    Umfpack3SparseLUSolutionManager(integer Dim):A(Dim), x(Dim), b(Dim), Symbolic(0), 
    						 Numeric(0) {
	xVH = new MyVectorHandler(Dim,&(x[0]));
	bVH = new MyVectorHandler(Dim,&(b[0]));
	 
	umfpack_defaults(Control);
    };    
    
     
   /* Distruttore: dealloca le matrici e distrugge gli oggetti propri */
   virtual ~Umfpack3SparseLUSolutionManager(void) {
   	if (Symbolic != 0) umfpack_free_symbolic(&Symbolic);
   	if (Numeric != 0) umfpack_free_numeric(&Numeric);
	delete xVH;
	delete bVH;
   };

   virtual void IsValid(void) const {NO_OP;};
   
   /* Inizializzatore generico */
   virtual void MatrInit(const doublereal& d = 0.) {A.Reset(d);};
   
   /* Risolve il sistema  Fattorizzazione + Bacward Substitution*/
   virtual void Solve(void){
   	A.MakeCompressedColumnForm(Ax,Ai,Ap);
	const double*const Axp=&(Ax[0]);
	const int*const Aip=&(Ai[0]);
	const int*const App=&(Ap[0]);
       t = umfpack_timer ( ) ;
       int pippo; pippo = b.size();
       status = umfpack_symbolic(pippo,App,Aip,&Symbolic,Control,Info);
       if (status != UMFPACK_OK) {
           umfpack_report_info(Control, Info) ;
	   umfpack_report_status(Control, status) ;
	   std::cerr << "umfpack_symbolic failed" << std::endl;
	   return;
       }
       // umfpack_report_symbolic ("Symbolic factorization of A",
       // 	Symbolic, Control) ;
       umfpack_report_info(Control, Info);
       t1 = umfpack_timer() - t ;
       status = umfpack_numeric(App,Aip,Axp,Symbolic,&Numeric,Control,Info);
       if (status != UMFPACK_OK) {
           umfpack_report_info(Control, Info) ;
	   umfpack_report_status(Control, status) ;
	   std::cerr << "umfpack_numeric failed" << std::endl;
	   //de-allocate memory
	   umfpack_free_symbolic(&Symbolic);
	   return;
       }
	
       // umfpack_report_numeric ("Numeric factorization of A",
       // 	Numeric, Control) ;
       umfpack_report_info(Control, Info);
       t1 = umfpack_timer() - t ;
       BackSub(t);
   };
   
   /* Bacward Substitution */
   void BackSub(doublereal t_iniz = 0){
	const double*const Axp=&(Ax[0]);
	const int*const Aip=&(Ai[0]);
	const int*const App=&(Ap[0]);
	const double*const bp=&(b[0]);
	double*const xp=&(x[0]);
        t = t_iniz;
	status = umfpack_solve("Ax=b",App,Aip,Axp,xp,bp,Numeric,Control,Info);
	if (status != UMFPACK_OK) {
		umfpack_report_info(Control, Info) ;
		umfpack_report_status(Control, status) ;
		std::cerr << "umfpack_solve failed" << std::endl;
		//de-allocate memory
		umfpack_free_symbolic(&Symbolic);
		umfpack_free_numeric(&Numeric);
		return;
	}
	umfpack_report_info(Control, Info);
	t1 = umfpack_timer() - t ;
   };
   
   /* Rende disponibile l'handler per la matrice */
   virtual SpMapMatrixHandler* pMatHdl(void) const{
       return &A;   
   
   };

   /* Rende disponibile l'handler per il termine noto */
   virtual MyVectorHandler* pResHdl(void) const{
       return bVH;   
   };

   /* Rende disponibile l'handler per la soluzione (e' lo stesso 
    * del termine noto, ma concettualmente sono separati) */
   virtual MyVectorHandler* pSolHdl(void) const{
       return xVH;      
   };
      
};

#endif /* USE_UMFPACK3 */

#endif //Umfpack3SparseLUSolutionManager_hh

