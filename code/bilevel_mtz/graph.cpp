#include "graph.h"
#include <math.h>

//using namespace std;

graph::graph(const char* filename) {
    //string filename = "D:\\Study\\Ph.D\\IE 670 Logistics\\hw1\\TSP instances\\TSP_instance_n_10_s_1.dat";
    std::fstream file(filename);
    if (!file) {
        std::cerr << "ERROR: could not open file '" << filename
            << "' for reading" << std::endl;
        throw (-1);
    }

    file >> n;

    p = 10;
    x = std::vector<double> (n);
    y = std::vector<double> (n);
    inc = { {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3}, {0, 1, 2, 3} };

    dist = std::vector<std::vector<double>> (n, std::vector<double> (n));
    for (int i = 0; i < n; i++) {
        file >> x[i] >> y[i];
    }

    for (int i = 0; i < n; i++)
        for (int j = i; j < n; j++)
            dist[i][j] = sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
}

double graph::getDist(int i, int j) {
    return sqrt(pow(x[i] - x[j], 2) + pow(y[i] - y[j], 2));
}

/*double graph::getDist(int i, int j) {
    if (i <= j)
        return dist[i][j];
    else if (j < i)
        return dist[j][i];
    else
        return 0;
}*/

/*graph::~graph() {

    for (int i = 0; i < n; i++) {
        delete(dist[i]);
    }
    delete x;
    delete y;
}*/
