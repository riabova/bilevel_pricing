#include<iostream>
#include "mtz_v1.h"

mtz_v1::mtz_v1(graph* g, int S) {
	//g = gr;
	n = g->getN();
	p = g->getP();
	inconv = g->getInc();
	K = n - S;
	env = new GRBEnv();
	model = new GRBModel(*env);

	x = std::vector<std::vector<GRBVar>> (n, std::vector<GRBVar>(n)); //routing
	u = std::vector<GRBVar> (n); //mtz order
	y1 = std::vector<GRBVar> (K); //follower decision
	y2 = std::vector<std::vector<GRBVar>> (K, std::vector<GRBVar>(S));
	z1 = std::vector<GRBVar> (K); //leader discount
	z2 = std::vector<std::vector<GRBVar>> (K, std::vector<GRBVar>(S));
	w1 = std::vector<GRBVar> (K); //linearization
	w2 = std::vector<std::vector<GRBVar>> (K, std::vector<GRBVar>(S));
	a = std::vector<GRBVar> (K); //follower dual

	for (int i = 0; i < n; i++) {
		u[i] = model->addVar(0, n, 0, GRB_CONTINUOUS, "u_" + std::to_string(i));
		for (int j = 0; j < n; j++) {
			x[i][j] = model->addVar(0, 1, g->getDist(i, j), GRB_BINARY, "x_" + std::to_string(i) + "," + std::to_string(j));
		}
	}

	for (int k = 1; k < K; k++) {
		y1[k] = model->addVar(0, 1, -p, GRB_CONTINUOUS, "y1_" + std::to_string(k));
		z1[k] = model->addVar(0, 1, 0, GRB_CONTINUOUS, "z1_" + std::to_string(k));
		w1[k] = model->addVar(0, 1, 1, GRB_CONTINUOUS, "w1_" + std::to_string(k));
		a[k] = model->addVar(0, 1, 0, GRB_CONTINUOUS, "a_" + std::to_string(k));
		for (int s = 0; s < S; s++) {
			y2[k][s] = model->addVar(0, 1, -p, GRB_CONTINUOUS, "y2_" + std::to_string(k) + "," + std::to_string(s));
			z2[k][s] = model->addVar(0, 1, 0, GRB_CONTINUOUS, "z2_" + std::to_string(k) + "," + std::to_string(s));
			w2[k][s] = model->addVar(0, 1, 1, GRB_CONTINUOUS, "z2_" + std::to_string(k) + "," + std::to_string(s));
		}
	}

	model->update();

	GRBLinExpr constr1_1 = 0; //balance pt.1 (need to visit); homes
	GRBLinExpr constr1_2 = 0; //balance pt.1 (need to visit); stores
	GRBLinExpr constr2_1 = 0; //balance pt.2 (need to leave if visited)
	GRBLinExpr constr2_2 = 0; //balance pt.2 (need to leave if visited)
	GRBLinExpr constr3 = 0; //duality 1 (optimality)
	GRBLinExpr constr4 = 0; //duality 2 (primal feasibility)
	GRBLinExpr constr5 = 0; //duality 2 (dual feasibility)

	//balance pt.1
	for (int k = 1; k < K; k++) {
		for (int j = 0; j < n; j++) {
			constr1_1 += x[j][k];
		}
		model->addConstr(constr1_1 == y1[k]); //balance1_1
		constr1_1 = 0;
		for (int s = 0; s < S; s++) {
			constr3 += w2[k][s] - inconv[k][s+1] * y2[k][s];
			constr4 += y2[k][s];
			model->addConstr(a[k] >= z2[k][s] - inconv[k][s+1]); //duality 2 (dual feasibility)
			model->addConstr(w2[k][s] <= z2[k][s]); //linearization 1
			model->addConstr(w2[k][s] <= p * y2[k][s]); //linearization 2
			model->addConstr(w2[k][s] >= z2[k][s] - p * (1 - y2[k][s])); //linearization 3
			for (int j = 0; j < n; j++) {
				constr1_2 += x[j][s + K];
			}
			model->addConstr(constr1_2 >= y2[k][s]); //balance1_2
			constr1_2 = 0;
		}
		model->addConstr(a[k] == w1[k] - inconv[k][0] * y1[k] + constr3); //duality 1 (optimality)
		model->addConstr(y1[k] + constr4 == 1); //duality 2 (primal feasibility)
		model->addConstr(a[k] >= z1[k] - inconv[k][0]); //duality 2 (dual feasibility)
		model->addConstr(w1[k] <= z1[k]); //linearization 1
		model->addConstr(w1[k] <= p * y1[k]); //linearization 2
		model->addConstr(w1[k] >= z1[k] - p * (1 - y1[k])); //linearization 3
		constr3 = 0;
		constr4 = 0;
	}

	//balance pt.2
	for (int i = 0; i < n; i++) {
		model->addConstr(x[i][i] == 0); //eliminating diag
		for (int j = 0; j < n; j++) {
			constr2_1 += x[j][i];
			constr2_2 += x[i][j];
		}
		model->addConstr(constr2_1 == constr2_2);
		constr2_1 = 0;
		constr2_2 = 0;
	}

	//subtour
	for (int i = 1; i < n; i++) {
		for (int j = 1; j < n; j++) {
			model->addConstr(u[i] - u[j] + n * x[i][j] <= n - 1);
		}
	}

	model->update();

}

void mtz_v1::solve() {
	try {
		model->optimize();
	}
	catch (GRBException e) {
		std::cout << "Error number: " << e.getErrorCode() << std::endl;
		std::cout << e.getMessage() << std::endl;
	}
	catch (...) {
		std::cout << "Error during optimization" << std::endl;
	}
}

void mtz_v1::printSol() {
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			if (x[i][j].get(GRB_DoubleAttr_X) > 0.5) {
				std::cout << i << " " << j << std::endl;
			}
		}
	}
}

mtz_v1::~mtz_v1() {
	delete env;
	delete model;
}