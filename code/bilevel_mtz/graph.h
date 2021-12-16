#pragma once
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>

class graph{
    private:
        int n; //number of vertices
        double p; //price
        std::vector<double> x; //x coordinates
        std::vector<double> y; //y coordinates
        std::vector<std::vector<double>> dist;
        std::vector<std::vector<double>> inc; //inconvenience

    public:
        graph(const char* filename);
        //~graph();
        graph() {};
        double getDist(int i, int j);
        inline int getN() { return n; };
        inline double getP() { return p; };
        inline std::vector<std::vector<double>> getInc() { return inc; };
};