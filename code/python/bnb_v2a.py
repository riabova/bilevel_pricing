import queue
import math
import time
import random
import numpy as np
#import bilevel_v4
import bilevel_v5
import folium
from folium.features import CustomIcon
import matplotlib.pyplot as plt
from matplotlib.offsetbox import OffsetImage, AnnotationBbox

class BNB:

    class Node:
        childRight = 0
        childLeft = 0
        depth = 0
        status = 0
        #status:
        #0 = unprocessed
        #1 = violating (need to branch)
        #2 = candidate
        #3 = pruned by bound
        #4 = infeasible
    
        def __init__(self, name: int, parent, obj: float, var):
            self.name = name
            self.parent = parent
            self.obj = obj
            self.var = var
            if parent != -1:
                self.depth = self.parent.depth + 1
    
        def __lt__(self, other):
            
            return (self.obj > other.obj)
    
    numNodes = 0
    nodeQueue = queue.PriorityQueue()
    bestNode = -1
    tol = 0.000001
    depth = 0
    UB = 0

    def __init__(self, graph, model, x, y, z, w, p, items, Lk, ui, L, distThrshd, I_coef):
        self.G = graph
        self.model = model
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.p = p
        self.n = graph.n
        self.items = items
        self.Lk = Lk
        self.ui = ui
        self.K = len(Lk) + 1
        self.S = len(ui[0])
        self.L = L
        self.distThrshd = distThrshd
        self.I_coef = I_coef
        self.model.Params.OutputFlag = 0
        self.model.setParam('TimeLimit', 36*60*60)
        start = time.time()
        self.model.optimize(bilevel_v5.callback_mult)
        self.time =  (time.time() - start)
        #self.model.Params.OutputFlag = 0
        root = self.Node(self.numNodes, -1, self.model.objVal, -1)
        root.parent = root
        self.nodeQueue.put(root)
        self.numNodes += 1

    def findConflictPair(self):
        for k in range(1, self.K):
            for i in range(self.K, self.K + self.S):#self.n):
                if self.G.dist[i, k] > self.distThrshd:
                        flag = 0
                        for l in self.Lk[k - 1]:
                            if self.y[k, l].x > 0.8:
                                y1 = self.y[k, l]
                                flag += 1
                                break
                        if flag > 0:
                            for l in self.Lk[k - 1]:
                                if self.y[i, l].x > 0.8:
                                    y2 = self.y[i, l]
                                    flag += 1
                                if flag > 1:
                                    return y1, y2
                for j in range(self.K, i):
                    if self.G.dist[i, j] > self.distThrshd:
                        flag = 0
                        for l in self.Lk[k - 1]:
                            if self.y[i, l].x > 0.8:
                                y1 = self.y[i, l]
                                flag += 1
                                break
                        if flag > 0:
                            for l in self.Lk[k - 1]:
                                if self.y[j, l].x > 0.8:
                                    y2 = self.y[j, l]
                                    flag += 1
                                if flag > 1:
                                    return y1, y2
        return -1

    def findConflictPair1a(self): #restricts any 2 stores together
        for k in range(1, self.K):
            for l in self.Lk[k - 1]:
                flag = 0
                sLst = []
                for s in range(self.S):
                    if self.y[s + self.K, l].x > 0.5:
                        flag += 1
                        sLst.append(s + self.K)
                        if flag > 1:
                            return sLst[0], sLst[1]
        return -1

    def solve(self):
        while not self.nodeQueue.empty():
            node = self.nodeQueue.get()
            if self.bestNode == -1 or self.bestNode.obj > math.floor(node.obj):
                self.processNode(node)
            else:  #pruned by bound (UB, obj of parent)
                node.status = 3

    def processNode(self, node):
        current = node
        while (current.name != 0):  #fixing variable of a branch
            var = current.var
            var.UB = 0
            current = current.parent

        if(node.name != 0):
            self.model.optimize(bilevel_v5.callback_mult)
            node.obj = self.model.objVal
        if (self.model.status != 3):
            if (self.bestNode == -1 or math.floor(node.obj) > self.bestNode.obj):
                confl = self.findConflictPair1a()
                print(confl)
                if confl == -1:  #candidate (no conflict)
                    node.status = 2
                    self.bestNode = node
                else:  #branching
                    node.status = 1
                    node.leaf = 0
                    node_l = self.Node(self.numNodes, node, node.obj, confl[0])
                    self.nodeQueue.put(node_l)
                    node.childLeft = node_l
                    self.numNodes += 1
                    node_r = self.Node(self.numNodes, node, node.obj, confl[1]) 
                    self.nodeQueue.put(node_r)
                    node.childRight = node_r
                    self.numNodes += 1
            
                    if(self.depth < node_l.depth):
                        self.depth = node_l.depth
            else:  #pruned by bound
                node.status = 3

        else:  #infeasible
            node.status = 4   

        current = node
        while (current.name != 0):  #releasing variables of a branch
            var = current.var
            var.UB = 1
            current = current.parent

    def store_sol_info(self):
        self.profit = -self.model.objVal
        #self.rCost = sum([self.model.getVarByName("C_%g").x % r for r in range(self.G.r)])
        #self.rev = self.model.getVarByName("rev").x
        self.rCost = self.model.getVarByName("rCost").x
        self.gap = self.model.MIPGap
        self.discStats = []
        '''for k in range(1, self.K):
            prices = []
            disc = 0
            for l in self.Lk[k - 1]:
                prices.append(self.items[l].price)
                for s in range(self.S):
                    if self.y[self.K + s, l].x >= 0.5:
                        disc += self.items[l].price - self.z[self.K + s, l].x
            self.discStats.append([disc, prices])'''

    def printSol(self):
        print("Leader's objective: %g" % -self.model.objVal)
        print("Customers decisions:")
        for k in range(1, self.K):
            print("Customer %d - Obj: %g   part1: %g   part2: %g" % (k, sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.w[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]) 
            - sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S)), sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.w[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]), 
            sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S))))
            for l in self.Lk[k - 1]:
                print("%d, %d: %g %g %g"  % (k, self.items[l].type, self.y[self.items[l].k, l].x, self.items[l].ul, self.items[l].price))
                print("%d, %d: %g -"  % (0, self.items[l].type, self.y[0, l].x))
                for s in range(self.S):
                    print("%d, %d: %g %g %g %g"  % (s + self.K, self.items[l].type, self.y[s + self.K, l].x,  self.items[l].ul, self.ui[k - 1][s], self.z[s + self.K, l].x))

        print("Routing:")
        for i, j, r in self.x.keys():
            if self.x[i, j, r].x > 0.5:
                print("%d, %d, %d" % (i, j, r))
        print("Routing cost: %g" % self.model.getVarByName("rCost").x)#sum([self.model.getVarByName("C_%g").x % r for r in range(self.G.r)]))
        #print("Revenue: %g" % self.model.getVarByName("rev").x)

    def plotRoute(self):
        home = OffsetImage(plt.imread("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\home.png"), zoom=0.55)
        wh = OffsetImage(plt.imread("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\warehouse.png"), zoom=0.55)
        store = OffsetImage(plt.imread("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\store.png"), zoom=0.55)
        fig, ax = plt.subplots(figsize=(10, 10))
        ax.axis('off')
        coords = np.transpose(self.G.points)
        plt.xlim([min(coords[0]) - 3, max(coords[0]) + 3])
        plt.ylim([min(coords[1]) - 3, max(coords[1]) + 3])
        ab = AnnotationBbox(wh, (self.G.points[0][0], self.G.points[0][1]), frameon=False)
        ax.add_artist(ab)
        for r in range(self.G.r):
            rgb = (random.random(), random.random(), random.random())
            for i in range(self.n):
                for j in range(i):
                    if self.x[i, j, r].x > 0.5:
                        plt.plot([self.G.points[i][0], self.G.points[j][0]], [self.G.points[i][1], self.G.points[j][1]], color=rgb, linewidth=4.0 )
        for k in range(1, self.K):
            ab = AnnotationBbox(home, (self.G.points[k][0], self.G.points[k][1]), frameon=False)
            ax.add_artist(ab)
        for s in range(self.S):
            ab = AnnotationBbox(store, (self.G.points[s + self.K][0], self.G.points[s + self.K][1]), frameon=False)
            ax.add_artist(ab)
        plt.show()
        #plt.savefig("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\results\img\\routs\\routs_n%d_s%d_r%d.png" % (self.n, self.S, self.G.r))

    def plotRouteMap(self):
        wh = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\warehouse.png", icon_size=(50, 25))
        m = folium.Map(location=(42.93, -78.79), zoom_start=11, tiles=None)
        base_map = folium.FeatureGroup(name='Basemap', overlay=True, control=False)
        folium.TileLayer(tiles='OpenStreetMap').add_to(base_map)
        base_map.add_to(m)
        fg = folium.FeatureGroup(name="stores", overlay=False, show=False)
        folium.Marker(location=self.G.points[0], icon=wh).add_to(fg)
        for s in range(self.S):
            store = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\store.png", icon_size=(20, 20))
            folium.Marker(location=self.G.points[self.K + s], icon=store, popup=s).add_to(fg)
        fg.add_to(m)
        fg1 = folium.FeatureGroup(name="custs", overlay=False, show=False)
        for i in range(1, self.K):
            home = CustomIcon("D:\Study\Ph.D\Projects\Bilevel Optimization\papers\img\loc_pricing\home.png", icon_size=(20, 20))
            folium.Marker(location=self.G.points[i], icon=home, popup=i).add_to(fg1)
        fg1.add_to(m)
        fg2 = folium.FeatureGroup(name="routs", overlay=False, show=False)
        colors = ["#990000", "#101D6B", "#028A0F", "#016064"]#["#000000", "#0000ff", "#028A0F", "#990000"]
        for r in range(self.G.r):
            rgb = colors[r]#(random.randint(0, 255), random.randint(0, 255), random.randint(0, 255))
            for i in range(self.n):
                for j in range(i):
                    if self.x[i, j, r].x > 0.5:
                        folium.PolyLine(self.G.routs[i, j], color=rgb, weight=3).add_to(fg2)
        fg2.add_to(m)
        m.save("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\results\img\\routs\\routs_n%d_s%d_r%d_I_%g.html" % (self.n, self.S, self.G.r, self.I_coef))