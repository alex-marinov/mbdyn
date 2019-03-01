/* stuttura dati lista nodi */

#ifndef NODE3_HH
#define NODE3_HH

#include <matvec3.h>
#include <vector>

class Node {
		  	  
  private:
		
	unsigned long int label;
	int Nadj;
	double dL;

        Vec3 Pos;
        std::vector<Vec3> PosAdj;
		
  public:
	
	Node() {};
	
	Node(unsigned long int id,const Vec3& P,int N = 0,double ref_len = 0.)
        :label(id),Nadj(N),dL(ref_len),Pos(P) {
        	if (Nadj!=0) this->AddAdjNod();
        };
 	 		
	~Node(void) {};
	
	// Copy constructor
	Node(const Node& N){
		label = N.label;
		Pos = N.Pos;
		Nadj = N.Nadj;
		PosAdj = N.PosAdj;
	}
	
	unsigned long int GetLabel(void) {
		return this->label;
	};

	/* per fare in modo che una lista di nodi possa essere ordinata */
	bool operator < (const Node& N) const {
		return (this->label < N.label);
	};
	
	void PutPos(const Vec3& P) {	 
		Pos = P;
	};
				
    	Vec3 GetPos(void) const {			 
		return Pos;
    	};
    	
    	int GetNadj(void) {
    		return Nadj;
    	}; 
    	
    	void AddAdjNod(){
		double N = 10.;
    		switch (Nadj) {
    			case 3: {
    				Vec3 a1(1., 0., 0.),a2(0., 1., 0.),a3 = a1.Cross(a2);
    				PosAdj.push_back(Pos+a1*dL/N);
				PosAdj.push_back(Pos+a2*dL/N);
				PosAdj.push_back(Pos+a3*dL/N);
				break;
				}
			case 4: {
				// Porcata
				Vec3 a1(0., 1., 0.),a2(0., 0., 1.);
				PosAdj.push_back(Pos+a1*dL/N);
				PosAdj.push_back(Pos-a1*dL/N);
				PosAdj.push_back(Pos+a2*dL/N);
				PosAdj.push_back(Pos-a2*dL/N);
				break;
				}
			case 6: {
				Vec3 a1(1., 0., 0.),a2(0., 1., 0.),a3 = a1.Cross(a2);
    				PosAdj.push_back(Pos+a1*dL/N);
    				PosAdj.push_back(Pos-a1*dL/N);
				PosAdj.push_back(Pos+a2*dL/N);
				PosAdj.push_back(Pos-a2*dL/N);
				PosAdj.push_back(Pos+a3*dL/N);
				PosAdj.push_back(Pos-a3*dL/N);
				break;
				}
			default:
				Nadj = 0;
				break;
		}		    				
    	};
    	
    	std::vector<Vec3> GetAdjPos(){
    		return PosAdj;
    	};
    	
    	std::vector<Vec3> AdjAcc(const Vec3& AngVel,const Vec3& AngAcc){
    		Mat3x3 omegap(MatCross, AngAcc);
    		Mat3x3 omegaomega(MatCrossCross, AngVel, AngVel);
    		std::vector<Vec3> delta_acc(Nadj);
    		for (int i_adj = 0; i_adj < Nadj; i_adj++){
    			// Calcolo del vettore distanza
    			Vec3 f = PosAdj[i_adj]-Pos;
    			delta_acc[i_adj] = omegap*f+omegaomega*f; 
    		}
    		return delta_acc;
    	};

	std::vector<Vec3> AdjPos(const Vec3& Rot){
		Mat3x3 theta(MatCross, Rot);
		std::vector<Vec3> delta_pos(Nadj);
		for (int i_adj = 0; i_adj < Nadj; i_adj++){
			// Vettore distanza
			Vec3 f = PosAdj[i_adj]-Pos;
			delta_pos[i_adj] = theta*f;
		}
		return delta_pos;
	}
    	
};

#endif // NODE3_HH
