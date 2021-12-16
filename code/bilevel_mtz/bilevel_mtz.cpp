// bilevel_mtz.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "gurobi_c++.h"
#include "graph.h"
#include "mtz_v1.h"
#include <sstream>

int main(int argc, char** argv)
{
    std::cout << argv[1];
    //std::vector< std::vector<int> > coords = { {13, 6}, {17, 19}, {3, 10}, {10, 5}, {14, 4}, {4, 3}, {9, 13}, {11, 9}, {0, 12}, {15, 16} };
    graph g = graph(argv[1]);

    mtz_v1 problem(&g, 3);
    problem.solve();
    problem.printSol();
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
