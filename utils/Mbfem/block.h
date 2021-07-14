#ifndef BLOCK_HH
#define BLOCK_HH

#include <iostream>
#include <stdio.h>
#include <vector>
#include <map>

#include "node3.h"
#include "force3.h"
#include <matvec3.h>
#include <fullmh.h>
#include <spmapmh.h>

#include "mls.h"

class Block{

	private:

		int BlockId;
		int Nmb,Nfem,Nfem_l;
		int Nadj;
		int Nstep;
		double dL;
		bool beam;
		std::vector<Node> MbNodes,FemNodes,FemNodes_l;
		FullMatrixHandler* FemMM;
		std::map<int,int> idnode;
		std::map<int,FullMatrixHandler*> MbAcc,MbVLo,MbSLo,MbDisp;
		std::map< int , std::vector<CLoad> > MbCLo;
		InterpMethod* pIntMap;
		SpMapMatrixHandler* H;
		std::map<int,FullMatrixHandler*> FemAcc,FemDisp;
		std::map<int,FullMatrixHandler*> FemIneCou;
			
	public:
	
		Block();

		Block(int,int,double,int,bool);
		
		~Block(void);
		
		int GetId();
		
		int GetNmb();
		
		int GetNfem();

		bool BeOnBe();

		bool isMM();
		
		void SetInterpMethod(int,int,bool);
		
		void SetNode(int,const Vec3&,bool,bool);
		
		void ShowNodes();
		
		void SetMbAcc(int,int, const Vec3&,const Vec3&,const Vec3&);

		void SetMbDisp(int,int, const Vec3&,const Vec3&);
		
		void ShowAcc();
		
		void SetMbLoad(int,int,const Vec3&, bool,int);

		void SetFemMM(const Vec3&,int);
		
		void CreateMap();

		void InterpAcc();
		
		void InterpCLoads();

		void InterpDisp();

		FullMatrixHandler* GetFemAcc(int);
		
		Vec3 GetFemAcc(int,int);
		
		Vec3 GetMbAcc(int,int);

		Vec3 GetFemRAcc(int,int);

		Vec3 GetMbRAcc(int,int);

		Vec3 GetMbDisp(int,int);

		Vec3 GetFemDisp(int,int);

		Vec3 GetFemRot(int,int);

		Vec3 GetMbRot(int,int);

		FullMatrixHandler* GetFemInCou(int);
		
		SpMapMatrixHandler* GetH();
		
		std::vector<Node> GetMbNodes();

		std::vector<Node> GetFemNodes();
		
		std::vector<CLoad> GetCLoads(int);

		void Summary();
		
};

#endif // BLOCK_HH
		
			
