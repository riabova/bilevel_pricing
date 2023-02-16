import math
import numpy as np
import RouteBuilder
from sklearn.cluster import KMeans
import polyline

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

    def readWithDists(self, ssf, ccf, scf, c_coords_file, s_coords_file, q, r, rf1=None, rf2=None, rf3=None):
        f1 = open(ssf, "r")
        f2 = open(ccf, "r")
        f3 = open(scf, "r")
        self.dist = {}
        self.K = int(f2.readline())
        self.q = q
        self.r = r
        self.S = int(f1.readline())
        self.n = self.K + self.S
        self.max_dist = 0
        for i in range(self.S):
            line = f1.readline().split(" ")
            for j in range(i):
                self.dist[(self.K + i, self.K + j)] = float(line[j])
        for i in range(self.K):
            line = f2.readline().split(" ")
            for j in range(i):
                self.dist[(i, j)] = float(line[j])
                self.max_dist = max(self.max_dist, self.dist[(i, j)])
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
        self.routs = {}
        if rf1:
            f6 = open(rf1, "r")
            for i in range(self.K):
                for j in range(i):
                    self.routs[(i, j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f6.readline().split(")[")[1][:-2].split("), ")]
        if rf2:
            f7 = open(rf2, "r")
            for i in range(self.S):
                for j in range(self.K):
                    self.routs[(self.K + i, j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f7.readline().split(")[")[1][:-2].split("), ")]
        if rf3:
            f8 = open(rf3, "r")
            for i in range(self.S):
                for j in range(i):
                    self.routs[(self.K + i, self.K + j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f8.readline().split(")[")[1][:-2].split("), ")]

    def readSampleWithDists(self, ssf, ccf, scf, c_coords_file, s_coords_file, K, S, q, r, rf1=None, rf2=None, rf3=None, tot_custs=None, tot_stores=None):
        f1 = open(ssf, "r")
        f2 = open(ccf, "r")
        f3 = open(scf, "r")
        self.dist = {}
        self.K = K
        self.q = q
        self.r = r
        self.S = S
        self.n = self.K + self.S
        self.max_dist = 0
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
                self.max_dist = max(self.max_dist, self.dist[(i, j)])
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
        self.routs = {}
        if rf1:
            f6 = open(rf1, "r")
            for i in range(tot_custs):
                for j in range(i):
                    if i <= self.K:
                        self.routs[(i, j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f6.readline().split(")[")[1][:-2].split("), ")]
                    else:
                        f6.readline()
        if rf2:
            f7 = open(rf2, "r")
            for i in range(tot_stores):
                for j in range(tot_custs):
                    if i <= self.S and j <= self.K:
                        self.routs[(self.K + i, j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f7.readline().split(")[")[1][:-2].split("), ")]
                    else:
                        f7.readline()
        if rf3:
            f8 = open(rf3, "r")
            for i in range(tot_stores):
                for j in range(i):
                    if i <= self.S:
                        self.routs[(self.K + i, self.K + j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f8.readline().split(")[")[1][:-2].split("), ")]
                    else:
                        f8.readline()

    def readSampleWOstores(self, ccf, c_coords_file, K, S, q, r, method="kmeans", rho=1, rf1=None):
        f1 = open(ccf, "r")
        self.dist = {}
        self.K = K
        self.q = q
        self.r = r
        self.S = S
        self.n = self.K + self.S
        self.routs = {}
        f1.readline()
        #read dists
        self.max_dist = 0
        for i in range(self.K):
            line = f1.readline().split(" ")
            for j in range(i):
                self.dist[(i, j)] = np.round(float(line[j]), decimals=5)
                self.max_dist = max(self.max_dist, self.dist[(i, j)])
        f2 = open(c_coords_file, "r")
        for i in range(self.K):
            line = f2.readline().split(", ")
            self.points.append([np.round(float(line[0]), decimals=5), np.round(float(line[1]), decimals=5)])
        #generate stores
        stores = self.makeStores(S, method=method, r=rho)
        for s in stores:
            self.points.append(list(s))
        #calcuate dists to stores
        rb = RouteBuilder.RouteBuilder()
        for s in range(S):
            for j in range(s):
                rt = rb.makeRequest(stores[s], stores[j])
                self.routs[self.K + s, self.K + j] = polyline.decode(rt['geometry'])
                self.dist[self.K + s, self.K + j] = rt["distance"]
            for k in range(self.K):
                rt = rb.makeRequest(self.points[self.K + s], self.points[k])
                self.routs[self.K + s, k] = polyline.decode(rt['geometry'])
                self.dist[self.K + s, k] = rt["distance"]
        if rf1:
            f3 = open(rf1, "r")
            for i in range(self.K):
                for j in range(i):
                    self.routs[(i, j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f3.readline().split(")[")[1][:-2].split("), ")]

    def makeStores(self, S, method="kmeans", r=1):
        if method == "kmeans":
            kmeans = KMeans(n_clusters=S)
            kmeans.fit(self.points[1:])
            return kmeans.cluster_centers_
        elif method == "radial":
            origin = self.points[0]
            coords = [[pt[0] - origin[0], pt[1] - origin[1]] for pt in self.points[1:]]
            pcoords = [[r, np.arctan(pt[1]/pt[0])] for pt in coords]
            kmeans = KMeans(n_clusters=S)
            kmeans.fit(pcoords)
            stores = kmeans.cluster_centers_
            return [[np.round(pt[0] * np.cos(pt[1]) + origin[0], decimals=5), np.round(pt[0] * np.sin(pt[1]) + origin[1], decimals=5)] for pt in stores]

    def get_coords(self):
        return self.points

    def get_dist(self, i, j):
        return math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2)))
