#pragma once
#include "gurobi_c++.h"
#include "graph.h"

class mtz_v1
{
	private:
		int n;
		int K;
		double p;
		std::vector<std::vector<double>> inconv; //inconvenience
		std::vector<std::vector<GRBVar>> x; //routing
		std::vector<GRBVar> u; //mtz order
		std::vector<GRBVar> y1; //follower decision (home)
		std::vector<std::vector<GRBVar>> y2; //follower decision (stores)
		std::vector<GRBVar> z1; //leader discount (home)
		std::vector<std::vector<GRBVar>> z2; //leader discount (stores)
		std::vector<GRBVar> w1; //linearization (home)
		std::vector<std::vector<GRBVar>> w2; //linearization (stores)
		std::vector<GRBVar> a; //follower dual
		GRBEnv* env;
		GRBModel* model;

	public:
		mtz_v1(graph* gr, int S);
		~mtz_v1();
		void solve();
		void printSol();
};

