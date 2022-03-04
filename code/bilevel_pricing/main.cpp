#include <iostream>
#include "gurobi_c++.h"
#include "instance.h"
#include <sstream>

int main(int argc, char** argv)
{
    //std::cout << argv[1];

    Instance problem = Instance(argv[1], 3, 2);

    //mtz_v1 problem(&G, 3);
    //problem.solve();
    //problem.printSol();
}