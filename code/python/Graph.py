import math
import numpy as np
import RouteBuilder
from sklearn.cluster import KMeans
import polyline
import folium
from folium.features import CustomIcon
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

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

    def readSampleWithDistsRC(self, ssf, ccf, scf, c_coords_file, s_coords_file, K, S, q, r, rf1=None, rf2=None, rf3=None, tot_custs=None, tot_stores=None, seed=5):
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
        self.points_tmp = {}
        f1.readline()
        f2.readline()
        np.random.seed(seed)
        self.k_nums = np.random.choice(range(1, tot_custs), size=K - 1, replace=False)
        self.k_nums.sort()
        self.k_nums = np.insert(self.k_nums, 0, 0)
        for i in range(S):
            line = f1.readline().split(" ")
            for j in range(i):
                self.dist[(self.K + i, self.K + j)] = float(line[j])
        i_new = 0
        for i in range(tot_custs):
            line = f2.readline().split(" ")
            if i in self.k_nums:
                j_new = 0
                for j in range(i):
                    if j in self.k_nums:
                        self.dist[(i_new, j_new)] = float(line[j])
                        self.max_dist = max(self.max_dist, self.dist[(i_new, j_new)])
                        j_new += 1
                i_new += 1
        for i in range(S):
            line = f3.readline().split(" ")
            j_new = 0
            for j in range(tot_custs):
                if j in self.k_nums:
                    self.dist[(self.K + i, j_new)] = float(line[j])
                    j_new += 1
        f4 = open(c_coords_file, "r")
        f5 = open(s_coords_file, "r")
        i_new = 0
        for i in range(tot_custs):
            line = f4.readline().split(", ")
            if i in self.k_nums:
                self.points.append([float(line[0]), float(line[1])])
                self.points_tmp[i_new] = [float(line[0]), float(line[1])]
                i_new += 1
        for i in range(self.S):
            line = f5.readline().split(", ")
            self.points.append([float(line[0]), float(line[1])])
        self.routs = {}
        i_new = 0
        if rf1:
            f6 = open(rf1, "r")
            for i in range(tot_custs):
                if i in self.k_nums:
                    j_new = 0
                    for j in range(i):
                        if j in self.k_nums:
                            self.routs[(i_new, j_new)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f6.readline().split(")[")[1][:-2].split("), ")]
                            j_new += 1
                        else:
                            f6.readline()
                    i_new += 1
                else:
                    for j in range(i):
                        f6.readline()
        if rf2:
            f7 = open(rf2, "r")
            for i in range(tot_stores):
                j_new = 0
                for j in range(tot_custs):
                    if j in self.k_nums:
                        self.routs[(self.K + i, j_new)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f7.readline().split(")[")[1][:-2].split("), ")]
                        j_new += 1
                    else:
                        f7.readline()
                i_new += 1
        if rf3:
            f8 = open(rf3, "r")
            for i in range(tot_stores):
                for j in range(i):
                    if i <= self.S:
                        self.routs[(self.K + i, self.K + j)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f8.readline().split(")[")[1][:-2].split("), ")]
                    else:
                        f8.readline()
        self.plotMap()


    def readRandSampleWithDists(self, ssf, ccf, scf, c_coords_file, s_coords_file, K, S, q, r, rf1=None, rf2=None, rf3=None, tot_custs=None, tot_stores=None, seed=5):
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
        self.points_tmp = {}
        f1.readline()
        f2.readline()
        np.random.seed(seed)
        #s_tmp = list(range(tot_stores))
        #k_tmp = list(range(1, tot_custs))
        self.s_nums = np.random.choice(range(tot_stores), size=S, replace=False)
        np.random.seed(seed)
        self.k_nums = np.random.choice(range(1, tot_custs), size=K - 1, replace=False)
        '''for s in range(S):
            i = np.random.randint(0, len(s_tmp))
            self.s_nums.append(s_tmp[i])
            del s_tmp[i]
        np.random.seed(seed)
        for k in range(K - 1):
            np.random.seed(seed)
            i = np.random.randint(0, len(k_tmp))
            self.k_nums.append(k_tmp[i])
            del k_tmp[i]'''
        #self.s_nums.sort()
        self.k_nums.sort()
        self.k_nums = np.insert(self.k_nums, 0, 0)
        #s_nums = np.random.choice(tot_stores, S, replace=False)
        #k_nums = np.random.choice(tot_custs, K, replace=False)
        i_new = 0
        for i in range(tot_stores):
            line = f1.readline().split(" ")
            if i in self.s_nums:
                j_new = 0
                for j in range(i):
                    if j in self.s_nums:
                        self.dist[(self.K + i_new, self.K + j_new)] = float(line[j])
                        j_new += 1
                i_new += 1
        i_new = 0
        for i in range(tot_custs):
            line = f2.readline().split(" ")
            if i in self.k_nums:
                j_new = 0
                for j in range(i):
                    if j in self.k_nums:
                        self.dist[(i_new, j_new)] = float(line[j])
                        self.max_dist = max(self.max_dist, self.dist[(i_new, j_new)])
                        j_new += 1
                i_new += 1
        i_new = 0
        for i in range(tot_stores):
            if i in self.s_nums:
                line = f3.readline().split(" ")
                j_new = 0
                for j in range(tot_custs):
                    if j in self.k_nums:
                        self.dist[(self.K + i_new, j_new)] = float(line[j])
                        j_new += 1
                i_new += 1
        f4 = open(c_coords_file, "r")
        f5 = open(s_coords_file, "r")
        i_new = 0
        for i in range(tot_custs):
            line = f4.readline().split(", ")
            if i in self.k_nums:
                self.points.append([float(line[0]), float(line[1])])
                self.points_tmp[i_new] = [float(line[0]), float(line[1])]
                i_new += 1
        i_new = 0
        for i in range(tot_stores):
            line = f5.readline().split(", ")
            if i in self.s_nums:
                self.points.append([float(line[0]), float(line[1])])
                self.points_tmp[i_new + self.K] = [float(line[0]), float(line[1])]
                i_new += 1
        self.routs = {}
        i_new = 0
        if rf1:
            f6 = open(rf1, "r")
            for i in range(tot_custs):
                if i in self.k_nums:
                    j_new = 0
                    for j in range(i):
                        if j in self.k_nums:
                            self.routs[(i_new, j_new)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f6.readline().split(")[")[1][:-2].split("), ")]
                            j_new += 1
                        else:
                            f6.readline()
                    i_new += 1
                else:
                    for j in range(i):
                        f6.readline()
        i_new = 0
        if rf2:
            f7 = open(rf2, "r")
            for i in range(tot_stores):
                if i in self.s_nums:
                    j_new = 0
                    for j in range(tot_custs):
                        if j in self.k_nums:
                            self.routs[(self.K + i_new, j_new)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f7.readline().split(")[")[1][:-2].split("), ")]
                            j_new += 1
                        else:
                            f7.readline()
                    i_new += 1
                else:
                    for j in range(tot_custs):
                        f7.readline()
        i_new = 0
        if rf3:
            f8 = open(rf3, "r")
            for i in range(tot_stores):
                if i in self.s_nums:
                    j_new = 0
                    for j in range(i):
                        if j in self.s_nums:
                            self.routs[(self.K + i_new, self.K + j_new)] = [(float(x[1:].split(", ")[0]), float(x[:-1].split(", ")[1])) for x in f8.readline().split(")[")[1][:-2].split("), ")]
                            j_new += 1
                        else:
                            f8.readline()
                    i_new += 1
                else:
                    for j in range(i):
                        f8.readline()
        self.plotMap()

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
        self.plotMap()
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
        self.writeStoreFiles(S, method=method, r=rho)

    def makeStores(self, S, method="kmeans", r=1):
        if method == "kmeans":
            kmeans = KMeans(n_clusters=S)
            kmeans.fit(self.points[1:])
            return kmeans.cluster_centers_
        elif method == "radial":
            origin = self.points[0]
            coords = [[pt[0] - origin[0], pt[1] - origin[1]] for pt in self.points[1:]]
            pcoords = [[r, np.arctan2(pt[1], pt[0])] for pt in coords]
            #self.points1 = [[np.round(pt[0] * np.cos(pt[1]) + origin[0], decimals=5), np.round(pt[0] * np.sin(pt[1]) + origin[1], decimals=5)] for pt in pcoords]
            #self.points2 = [[pt[0] + origin[0], pt[1] + origin[1]] for pt in coords]
            kmeans = KMeans(n_clusters=S)
            kmeans.fit(pcoords)
            stores = kmeans.cluster_centers_
            return [[np.round(r * np.cos(pt[1]) + origin[0], decimals=5), np.round(r * np.sin(pt[1]) + origin[1], decimals=5)] for pt in stores]
        
    def writeStoreFiles(self, S, method="kmeans", r=1):
        with open("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\Buffalo\\s_coords_" + method + str(r) + ".txt", 'w') as f:
            for s in range(S):
                f.write(", ".join([str(x) for x in self.points[self.K + s]]) + "\n")
        with open("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\Buffalo\\ss_routs_" + method + str(r) + ".txt", 'w') as f1:
            for s in range(S):
                for j in range(s):
                    f1.write(str((self.K + s, self.K + j)) + str(self.routs[self.K + s, self.K + j]) + "\n")
        with open("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\Buffalo\\sc_routs_" + method + str(r) + ".txt", 'w') as f1:
            for s in range(S):
                for j in range(self.K):
                    f1.write(str((self.K + s, j)) + str(self.routs[self.K + s, j]) + "\n")
        with open("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\Buffalo\\ss_dists_" + method + str(r) + ".txt", 'w') as f1:
            f1.write(str(S) + "\n")
            for s in range(S):
                for j in range(s):
                    f1.write(str(self.dist[self.K + s, self.K + j]) + " ")
                f1.write("\n")
        with open("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\Buffalo\\sc_dists_" + method + str(r) + ".txt", 'w') as f1:
            for s in range(S):
                for j in range(self.K):
                    f1.write(str(self.dist[self.K + s, j]) + " ")
                f1.write("\n")
        
    def plotMap(self):
        wh = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\warehouse.png", icon_size=(50, 25))
        m = folium.Map(location=(42.93, -78.79), zoom_start=11, tiles=None)
        base_map = folium.FeatureGroup(name='Basemap', overlay=True, control=False)
        folium.TileLayer(tiles='OpenStreetMap').add_to(base_map)
        base_map.add_to(m)
        fg = folium.FeatureGroup(name="stores", overlay=False, show=False)
        folium.Marker(location=self.points[0], icon=wh).add_to(fg)
        for s in range(self.S):
            store = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\store.png", icon_size=(20, 20))
            folium.Marker(location=self.points[self.K + s], icon=store, popup=s).add_to(fg)
        fg.add_to(m)
        fg1 = folium.FeatureGroup(name="custs", overlay=False, show=False)
        for i in range(1, self.K):
            home = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\home.png", icon_size=(20, 20))
            folium.Marker(location=self.points[i - 1], icon=home, popup=i).add_to(fg1)
        fg1.add_to(m)
        m.save("D:\Study\Ph.D\Projects\Bilevel Optimization\\code\\python\\out\\case\\img\\map_n%d_s%d.html" % (self.n, self.S))

    def get_coords(self):
        return self.points

    def get_dist(self, i, j):
        return math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2)))
