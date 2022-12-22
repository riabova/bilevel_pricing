import math
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

class Graph:
    points = []
    name = ""

    def read(self, file: str, S: int, seed: int):
        f = open(file, "r")
        lines = f.readlines()
        self.name = lines[0].split(": ")[1]
        self.n = int(lines[3].split(": ")[1])
        self.r = int(lines[0].split("-k")[1])
        self.q = int(lines[5].split(": ")[1])
        np.random.seed(seed)
        self.stores = sorted(np.random.randint(1, self.n, S))
        i = -1
        for s in self.stores:
            for pt in range(7 + i + 1, 7 + s):
                self.points.append([float(lines[pt].split()[1]), float(lines[pt].split()[2])])
            i = s
        for pt in range(7 + i + 1, 7 + self.n):
            self.points.append([float(lines[pt].split()[1]), float(lines[pt].split()[2])])
        for s in self.stores:
            self.points.append([float(lines[7 + s].split()[1]), float(lines[7 + s].split()[2])])
        f.close()
        self.dist = {(i, j):
        math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2))) for i in range(self.n) for j in range(i)}

    def read1(self, file: str, S: int, seed: int):
        f = open(file, "r")
        lines = f.readlines()
        self.name = lines[0].split(": ")[1][:-1]
        self.n = int(lines[3].split(": ")[1])
        self.r = int(lines[0].split("-k")[1])
        self.q = int(lines[5].split(": ")[1])
        self.S = S
        self.K = self.n
        self.demands = []
        np.random.seed(seed)
        for pt in range(7, 7 + self.n):
            self.points.append([float(lines[pt].split()[1]), float(lines[pt].split()[2])])
        for k in range(7 + self.n + 2, 7 + self.n + self.K + 1):
            self.demands.append(int(lines[k].split()[1]))
        f.close()
        kmeans = KMeans(n_clusters=S)
        kmeans.fit(self.points[1:])
        stores  = kmeans.cluster_centers_
        #plt.plot(self.points[0][0], self.points[0][1], "o")
        #plt.plot(np.transpose(self.points)[0], np.transpose(self.points)[1], "o")
        #plt.plot(stores[:, 0], stores[:, 1], "o")
        #plt.show()
        for s in stores:
            self.points.append(list(s))
        self.n = self.n + S
        self.dist = {(i, j):
        math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2))) for i in range(self.n) for j in range(i)}

    def readWithDists(self, ssf, ccf, scf, c_coords_file, s_coords_file, q, r):
        f1 = open(ssf, "r")
        f2 = open(ccf, "r")
        f3 = open(scf, "r")
        self.dist = {}
        self.K = int(f2.readline())
        self.q = q
        self.r = r
        self.S = int(f1.readline())
        self.n = self.K + self.S
        for i in range(self.S):
            line = f1.readline().split(" ")
            for j in range(i):
                self.dist[(self.K + i, self.K + j)] = float(line[j])
        for i in range(self.K):
            line = f2.readline().split(" ")
            for j in range(i):
                self.dist[(i, j)] = float(line[j])
        for i in range(self.S):
            line = f3.readline().split(" ")
            for j in range(self.K):
                self.dist[(self.K + i, j)] = float(line[j])
        f4 = open(c_coords_file, "r")
        f5 = open(s_coords_file, "r")
        for i in range(self.K):
            line = f4.readline().split(", ")
            self.points.append([float(line[0]), float(line[1])])
        for i in range(self.S):
            line = f5.readline().split(", ")
            self.points.append([float(line[0]), float(line[1])])

    def readSampleWithDists(self, ssf, ccf, scf, c_coords_file, s_coords_file, K, S, q, r):
        f1 = open(ssf, "r")
        f2 = open(ccf, "r")
        f3 = open(scf, "r")
        self.dist = {}
        self.K = K
        self.q = q
        self.r = r
        self.S = S
        self.n = self.K + self.S
        f1.readline()
        f2.readline()
        for i in range(S):
            line = f1.readline().split(" ")
            for j in range(i):
                self.dist[(self.K + i, self.K + j)] = float(line[j])
        for i in range(self.K):
            line = f2.readline().split(" ")
            for j in range(i):
                self.dist[(i, j)] = float(line[j])
        for i in range(S):
            line = f3.readline().split(" ")
            for j in range(self.K):
                self.dist[(self.K + i, j)] = float(line[j])
        f4 = open(c_coords_file, "r")
        f5 = open(s_coords_file, "r")
        for i in range(self.K):
            line = f4.readline().split(", ")
            self.points.append([float(line[0]), float(line[1])])
        for i in range(self.S):
            line = f5.readline().split(", ")
            self.points.append([float(line[0]), float(line[1])])

    def get_coords(self):
        return self.points

    def get_dist(self, i, j):
        return math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2)))
