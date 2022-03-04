#pragma once
#include <vector>;
#include <fstream>;
#include <sstream>;
#include <iostream>;
#include <string>;
#include <vector>;
#include <random>;

struct Item {
	int k;
	int w;
	float u;
	Item();
	~Item();
};

struct Instance {
	int S; //num stores
	int N; //num nodes
	int P;
	int K; //num custs
	int R; //num vehs
	int q; //capac
	std::vector<std::vector<double>> dists;
	std::vector<std::vector<int>> inconvs;
	std::vector<Item> items;
	Instance(const char* filename, int s, int p);
	~Instance();
};

