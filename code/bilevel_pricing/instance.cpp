#include "instance.h"

Instance::Instance(const char* filename, int s, int p) {
    std::ifstream file(filename);
    if (!file) {
        std::cerr << "ERROR: could not open file '" << filename
            << "' for reading" << std::endl;
        throw (-1);
    }
    S = s;
    P = p;
    // read from file
    std::string line;
    std::getline(file, line);
    std::cout << line.substr(line.find("-k") + 2, line.length());
    R = std::stoi(line.substr(line.find("-k") + 2, line.length()));
    std::getline(file, line);
    std::getline(file, line);
    std::getline(file, line);
    N = std::stoi(line.substr(line.find(": ") + 2, line.length()));
    K = N - S - 1;
    // decide which points are stores
    std::vector<int> stores(S);
    std::vector<int> storesX(S);
    std::vector<int> storesY(S);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distr(1, N);
    for (int i = 0; i < S; i++) {
        stores[i] = distr(gen);
    }
    std::sort(stores.begin(), stores.end());
    // continue reading
    std::getline(file, line);
    std::getline(file, line);
    q = std::stoi(line.substr(line.find(": ") + 2, line.length()));
    std::getline(file, line);
    std::vector<double> x(N);
    std::vector<double> y(N);
    int n;
    int j = 0;
    for (int s = 0; s < S; s++) {
        for (int i = j; i < stores[s]; i++) {
            file >> n >> x[i] >> y[i];
        }
        j = stores[s];
        file >> n >> storesX[s] >> storesY[s];
    }
    int a = 9;
}

Item::Item() { k = 0; w = 0; u = 0; };
Item::~Item() { };
Instance::~Instance() { };
