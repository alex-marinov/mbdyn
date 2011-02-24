/* stuttura dati lista forze */

#ifndef LOAD3_HH
#define LOAD3_HH

#include <matvec3.h>

class CLoad {
	
	private:
		
		int labelMb,labelFem;	 // Nodo a cui è associata la forza nei due modelli

		Vec3 Loa;	  	 // Vettore del carico (modello MB)
		
		Vec3 posMb,posFem;	 // Punti di applicazione della forza (MB e FEM)
		
		bool couple;		 // True if load is a couple
		
	public:
	
		CLoad(int id,const Vec3& L,bool c=false,const Vec3& P = Zero3): labelMb(id),labelFem(0),Loa(L),posMb(P),posFem(),couple(c){};
		
		~CLoad(void) {};
		
		unsigned long int GetLabMb(void) {
			return labelMb;
		};
		
		unsigned long int GetLabFem(void) {
			return labelFem;
		};
		
		Vec3 GetLoad() const {			 
			return Loa;
    		};
    		
    		Vec3 GetPosMb() const {			 
			return posMb;
    		};
    		
    		Vec3 GetPosFem() const {			 
			return posFem;
    		};
    		
    		bool IsCouple() {
    			return couple;
    		};

		void PutLabFem(const int id){
			labelFem = id;
		};
    		
    		void PutPosMb(const Vec3& P) {	 
			posMb = P;
		};
		
		void PutPosFem(const Vec3& P) {	 
			posFem = P;
		};
		
		void PutLoad(const Vec3& L) {
			Loa = L;
		};

};

#endif // LOAD3_HH
