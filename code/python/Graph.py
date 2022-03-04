import math
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt

class Graph:
    points = []

    def read(self, file: str, S: int, seed: int):
        f = open(file, "r")
        lines = f.readlines()
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
        self.n = int(lines[3].split(": ")[1])
        self.r = int(lines[0].split("-k")[1])
        self.q = int(lines[5].split(": ")[1])
        np.random.seed(seed)
        for pt in range(7, 7 + self.n):
            self.points.append([float(lines[pt].split()[1]), float(lines[pt].split()[2])])
        f.close()
        kmeans = KMeans(n_clusters=S)
        kmeans.fit(self.points)
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


    def get_coords(self):
        return self.points

    def get_dist(self, i, j):
        return math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2)))
