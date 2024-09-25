from gurobipy import *
import Graph
from itertools import combinations
import copy
import numpy as np
import time
import matplotlib.pyplot as plt
# import bilevel_v5_approx
import bilevel_v5_full
# import bilevel_v5_tsp


def callback_mult(model, where):
    if where == GRB.Callback.MIPSOL:
        G = model._G
        g_sol = model.cbGetSolution(model._g)
        v_sol = model.cbGetSolution(model._v)
        y_sol = model.cbGetSolution(model._y)
        x_sol = model.cbGetSolution(model._x)
        for k in range(model._R):
            callback(model, G, x_sol, k)

def callback(model, G, x_sol, k):
    for sc in range(model._nscen):
        adj = {}
        for i in range(G.n):
            adj[i] = []
        for i in range(1, G.n):
            for j in range(1, i):
                if x_sol[i, j, k, sc] > 0.5:
                    adj[i].append(j)
                    adj[j].append(i)
        cc = DFS(G, adj)
        if len(cc) > 1:
            for comp in cc:
                #print(comp)
                if len(comp) > 1:
                    #gh = list(combinations(comp, 2))
                    #print(quicksum(model._x[max(a), min(a), k] for a in combinations(comp, 2)) <= len(comp) - 1)
                    model.cbLazy(quicksum(model._x[max(a), min(a), k, sc] for a in combinations(comp, 2)) <= len(comp) - 1)

def callbackTSP(model, where):
    if where == GRB.Callback.MIPSOL:
        for sc in range(model._nscen):
            adj = {}
            G = model._G
            for i in range(G.n):
                adj[i] = []
            x_sol = model.cbGetSolution(model._x)
            for i, j in model._G.dist.keys():
                if x_sol[i, j, sc] > 0.5:
                    adj[i].append(j)
                    adj[j].append(i)
            cc = DFS(G, adj)
            if len(cc) > 1:
                for comp in cc:
                    model.cbLazy(quicksum(model._x[i, j, sc] + model._x[j, i, sc] for i, j in combinations(comp, 2)) <= len(comp) - 1)

def callbackTSP2(model, where):
    if where == GRB.Callback.MIPSOL:
        for sc in range(model._nscen):
            adj = {}
            G = model._G
            x_sol = model.cbGetSolution(model._x)
            for i in range(G.n):
                adj[i] = []
            for i in range(1, G.n):
                for j in range(1, i):
                    if x_sol[i, j, sc] > 0.5:
                        adj[i].append(j)
                        adj[j].append(i)
            cc = DFS(G, adj)
            if len(cc) > 1:
                for comp in cc:
                    #print(comp)
                    if len(comp) > 1:
                        #gh = list(combinations(comp, 2))
                        #print(quicksum(model._x[max(a), min(a), k] for a in combinations(comp, 2)) <= len(comp) - 1)
                        model.cbLazy(quicksum(model._x[max(a), min(a), sc] for a in combinations(comp, 2)) <= len(comp) - 1)

def DFS(G, adj):
    color = {vert : 'white' for vert in range(G.n)}
    cc = []
    for vert in range(G.n):
        if color[vert] == 'white':
            comp = []
            cc.append(process(G, comp, adj, vert, color))
    return cc

def process(G, comp, adj, vert, color):
    color[vert] = 'black'
    comp.append(vert)
    for v in adj[vert]:
        if color[v] == 'white':
            comp = process(G, comp, adj, v, color)
    return comp



class UnsrtLearner:
    def __init__(self, G: Graph.Graph, items: list(), Lk: list(), ui, S: int, q: list(), c: float, I_range: list()):

        self.G = G
        self.S = S
        self.K = self.G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
        self.R = G.r
        self.items = items
        self.Lk = Lk
        self.q = q
        self.c = c
        self.uiT = ui #true inconv values
        self.I_range = I_range #unsrt range for inconv
        self.ui = []
        for ir_k in I_range:
            #self.ui.append([(ir_k[s][0] + ir_k[s][1])/2 for s in range(S)])
            self.ui.append([ir_k[s][1] for s in range(self.S)])
        #self.ui = [[(x[s][1] - x[s][0])/2 for s in range(S)] for x in I_range] #point estimate (avg) for inconv

    def setupMaster(self):
        #init main problem
        self.model = Model("bl5")
        self.model.setParam("OutputFlag", 0)
        self.model.setParam("TimeLimit", 3 * 60)
        self.x = {} #leader's routing
        self.C = {} #routs oredring by cost
        self.y  ={} #follower's decision (continuous for duality)
        self.z = {} #leader's discount
        self.v = {} #leader's product-location/vehicle mapping
        self.w = {} #linerization
        self.a = {} #follower's dual 1
        self.b = {} #follower's dual 2
        self.p = {} #follower visits i
        self.g = {}

        self.rCost = self.model.addVars(len(self.scenarios), vtype=GRB.CONTINUOUS, name="rCost")
        for i in range(1, self.G.n):
            for r in range(self.R):
                for sc in range(len(self.scenarios)):
                    self.x[i, 0, r, sc] = self.model.addVar(vtype=GRB.BINARY, ub=2, name="x_%d_%d_%d_%d" % (i, 0, r, sc))
                for j in range(1, i):
                    for sc in range(len(self.scenarios)):
                        self.x[i, j, r, sc] = self.model.addVar(vtype=GRB.BINARY, name="x_%d_%d_%d_%d" % (i, j, r, sc))
        for r in range(self.R):
            for sc in range(len(self.scenarios)):
                self.C[r, sc] = self.model.addVar(obj=1/len(self.scenarios), vtype=GRB.CONTINUOUS, name="C_%d_%d" % (r, sc))
            for i in range(self.G.n):
                for sc in range(len(self.scenarios)):
                    self.g[r, i, sc] = self.model.addVar(vtype = GRB.BINARY, name="g_%d_%d_%d" % (r, i, sc))
    
        for l in range(len(self.items)):
            for sc in range(len(self.scenarios)):
                self.y[self.items[l].k, l, sc] = self.model.addVar(obj=-self.items[l].price/len(self.scenarios), vtype=GRB.BINARY, name='y_%d_%d_%d' % (self.items[l].k, l, sc))
                self.a[l, sc] = self.model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%d_%d' % (l, sc))
            for r in range(self.R):
                for sc in range(len(self.scenarios)):
                    self.v[self.items[l].k, l, r, sc] = self.model.addVar(vtype=GRB.BINARY, name="v_%d_%d_%d_%d" % (self.items[l].k, l, r, sc))
            for s in range(self.S):
                self.z[s + self.K, l] = self.model.addVar(vtype=GRB.CONTINUOUS, lb=self.items[l].lb, ub=self.items[l].price, name='z_%d_%d' % (s + self.K, l)) #VARIABLE PRICE, NOT DISCOUNT
                for sc in range(len(self.scenarios)):
                    self.y[s + self.K, l, sc] = self.model.addVar(vtype=GRB.BINARY, name='y_%d_%d_%d' % (s +self. K, l, sc))
                    self.w[s + self.K, l, sc] = self.model.addVar(obj=-1/len(self.scenarios), vtype=GRB.CONTINUOUS, name='w_%d_%d_%d' % (s + self.K, l, sc))
                    self.b[s + self.K, l, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='b_%d_%d_%d' % (s + self.K, l, sc))
                for r in range(self.R):
                    for sc in range(len(self.scenarios)):
                        self.v[s + self.K, l, r, sc] = self.model.addVar(vtype=GRB.BINARY, name="v_%d_%d_%d_%d" % (s + self.K, l, r, sc))
        for k in range(1, self.K):
            for s in range(self.S):
                for sc in range(len(self.scenarios)):
                    self.p[k, s + self.K, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='p_%d_%d_%d' % (k, s + self.K, sc))
    
        self.model.update()
    
    
        self.model.addConstrs(self.rCost[sc] == quicksum(self.C[r, sc] for r in range(self.R)) for sc in range(len(self.scenarios)))
    
        #balance
        self.model.addConstrs(self.C[r, sc] <= self.C[r - 1, sc] for r in range(1, self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.C[r, sc] >= quicksum(self.c * self.G.dist[i, j] * self.x[i, j, r, sc] for i in range(self.G.n) for j in range(i)) for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[j, 0, r, sc] for j in range(1, self.G.n)) == 2 * self.g[r, 0, sc] for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[i, j, r, sc] for j in range(i)) + quicksum(self.x[j, i, r, sc] for j in range(i + 1, self.G.n)) == 2 * self.g[r, i, sc] for r in range(self.R) for i in range(self.G.n) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[r, 0, sc] >= self.v[self.items[l].k, l, r, sc] for l in range(len(self.items)) for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[r, 0, sc] >= self.v[s + self.K, l, r, sc] for s in range(self.S) for l in range(len(self.items)) for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[r, self.items[l].k, sc] >= self.v[self.items[l].k, l, r, sc] for l in range(len(self.items)) for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[r, k, sc] <= quicksum(self.v[k, l, r, sc] for l in self.Lk[k - 1]) for k in range(1, self.K) for r in range(self.R) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[r, s + self.K, sc] >= self.v[s + self.K, l, r, sc] for s in range(self.S) for l in range(len(self.items)) for r in range(self.R) for sc in range(len(self.scenarios)))
    
        #capacity
        self.model.addConstrs(quicksum(self.v[self.items[l].k, l, r, sc] for r in range(self.R)) == self.y[self.items[l].k, l, sc] for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.v[s + self.K, l, r, sc] for r in range(self.R)) == self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        # self.model.addConstrs(quicksum(self.items[l].w*(self.v[self.items[l].k, l, r, sc] + quicksum(self.v[s + self.K, l, r, sc] for s in range(self.S))) for l in range(len(self.items))) <= self.q[r] for r in range(self.R) for sc in range(len(self.scenarios)))
            
        #duality
        self.uConstrs1 = self.model.addConstrs(quicksum(self.a[l, sc] for l in self.Lk[k - 1]) == quicksum(self.items[l].ul - self.items[l].price * self.y[self.items[l].k, l, sc] - quicksum(self.w[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(self.scenarios[sc][k, s + self.K] * self.p[k, s + self.K, sc] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.y[self.items[l].k, l, sc] + quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) == 1 for l in range(len(self.items)) for sc in range(len(self.scenarios)))     
        self.model.addConstrs(self.y[s + self.K, l, sc] <= self.p[k, s + self.K, sc] for k in range(1, self.K) for l in self.Lk[k - 1] for s in range(self.S) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] + self.b[s + self.K, l, sc] >= self.items[l].ul - self.z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.uConstrs2 = self.model.addConstrs(quicksum(self.b[s + self.K, l, sc] for l in self.Lk[k - 1]) <= self.scenarios[sc][k, s + self.K] for s in range(self.S) for k in range(1, self.K) for sc in range(len(self.scenarios)))
    
        #linearization
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.items[l].price * self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.z[s + self.K, l] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] >= self.z[s + self.K, l] - 1.1 * self.items[l].price*(1 - self.y[s +self. K, l, sc]) for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
    
    
        self.model._x = self.x
        self.model._g = self.g
        self.model._v = self.v
        self.model._y = self.y
        self.model._G = self.G
        self.model._R = self.R
        self.model._nscen = len(self.scenarios)
        self.model.Params.lazyConstraints = 1

    def setupTSP2(self):
        #init main problem
        self.model = Model("bl5")
        self.model.setParam("OutputFlag", 0)
        self.model.setParam("TimeLimit", 60)
        self.x = {} #leader's routing
        self.C = {} #routs oredring by cost
        self.y  ={} #follower's decision (continuous for duality)
        self.z = {} #leader's discount
        self.v = {} #leader's product-location/vehicle mapping
        self.w = {} #linerization
        self.a = {} #follower's dual 1
        self.b = {} #follower's dual 2
        self.p = {} #follower visits i
        self.g = {}

        self.rCost = self.model.addVars(len(self.scenarios), obj=1/len(self.scenarios), vtype=GRB.CONTINUOUS, name="rCost")
        for i in range(1, self.G.n):
            for sc in range(len(self.scenarios)):
                self.x[i, 0, sc] = self.model.addVar(vtype=GRB.BINARY, ub=2, name="x_%d_%d_%d" % (i, 0, sc))
            for j in range(1, i):
                for sc in range(len(self.scenarios)):
                    self.x[i, j, sc] = self.model.addVar(vtype=GRB.BINARY, name="x_%d_%d_%d" % (i, j, sc))
        for i in range(self.G.n):
            for sc in range(len(self.scenarios)):
                self.g[i, sc] = self.model.addVar(vtype = GRB.BINARY, name="g_%d_%d" % (i, sc))
    
        for l in range(len(self.items)):
            for sc in range(len(self.scenarios)):
                self.y[self.items[l].k, l, sc] = self.model.addVar(obj=-self.items[l].price/len(self.scenarios), vtype=GRB.CONTINUOUS, name='y_%d_%d_%d' % (self.items[l].k, l, sc))
                self.a[l, sc] = self.model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%d_%d' % (l, sc))
            for s in range(self.S):
                self.z[s + self.K, l] = self.model.addVar(vtype=GRB.CONTINUOUS, lb=self.items[l].lb, ub=self.items[l].price, name='z_%d_%d' % (s + self.K, l)) #VARIABLE PRICE, NOT DISCOUNT
                for sc in range(len(self.scenarios)):
                    self.y[s + self.K, l, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='y_%d_%d_%d' % (s +self. K, l, sc))
                    self.w[s + self.K, l, sc] = self.model.addVar(obj=-1/len(self.scenarios), vtype=GRB.CONTINUOUS, name='w_%d_%d_%d' % (s + self.K, l, sc))
                    self.b[s + self.K, l, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='b_%d_%d_%d' % (s + self.K, l, sc))
        for k in range(1, self.K):
            for s in range(self.S):
                for sc in range(len(self.scenarios)):
                    self.p[k, s + self.K, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='p_%d_%d_%d' % (k, s + self.K, sc))
    
        self.model.update()
    
    
        self.model.addConstrs(self.rCost[sc] == quicksum(self.c * self.G.dist[i, j] * self.x[i, j, sc] for i in range(self.G.n) for j in range(i)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[j, 0, sc] for j in range(1, self.G.n)) == 2 * self.g[0, sc] for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[i, j, sc] for j in range(i)) + quicksum(self.x[j, i, sc] for j in range(i + 1, self.G.n)) == 2 * self.g[i, sc] for i in range(self.G.n) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[self.items[l].k, sc] >= self.y[self.items[l].k, l, sc] for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[s + self.K, sc] >= self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
    
        #duality
        self.uConstrs1 = self.model.addConstrs(quicksum(self.a[l, sc] for l in self.Lk[k - 1]) == quicksum(self.items[l].ul - self.items[l].price * self.y[self.items[l].k, l, sc] - quicksum(self.w[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(self.scenarios[sc][k, s + self.K] * self.p[k, s + self.K, sc] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.y[self.items[l].k, l, sc] + quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) == 1 for l in range(len(self.items)) for sc in range(len(self.scenarios)))     
        self.model.addConstrs(self.y[s + self.K, l, sc] <= self.p[k, s + self.K, sc] for k in range(1, self.K) for l in self.Lk[k - 1] for s in range(self.S) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] + self.b[s + self.K, l, sc] >= self.items[l].ul - self.z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.uConstrs2 = self.model.addConstrs(quicksum(self.b[s + self.K, l, sc] for l in self.Lk[k - 1]) <= self.scenarios[sc][k, s + self.K] for s in range(self.S) for k in range(1, self.K) for sc in range(len(self.scenarios)))
    
        #linearization
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.items[l].price * self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.z[s + self.K, l] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] >= self.z[s + self.K, l] - 1.1 * self.items[l].price*(1 - self.y[s +self. K, l, sc]) for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
    
    
        self.model._x = self.x
        self.model._g = self.g
        self.model._v = self.v
        self.model._y = self.y
        self.model._G = self.G
        self.model._R = self.R
        self.model._nscen = len(self.scenarios)
        self.model.Params.lazyConstraints = 1

    def setupTSP(self):
        #init main problem
        self.model = Model("bl_tsp")
        self.model.setParam("OutputFlag", 0)
        self.model.setParam("TimeLimit", 60)
        self.x = {} #leader's routing
        self.y  ={} #follower's decision (continuous for duality)
        self.z = {} #leader's discount
        self.v = {} #leader's product-location/vehicle mapping
        self.w = {} #linerization
        self.a = {} #follower's dual 1
        self.b = {} #follower's dual 2
        self.p = {} #follower visits i
        self.g = {}

        self.rCost = self.model.addVars(len(self.scenarios), vtype=GRB.CONTINUOUS, name="rCost")
        for i, j in self.G.dist.keys():
            for sc in range(len(self.scenarios)):
                self.x[i, j, sc] = self.model.addVar(obj=self.c * self.G.dist[i, j]/len(self.scenarios), vtype=GRB.BINARY, name='x_%d_%d_%d' % (i, j, sc))
                self.x[j, i, sc] = self.model.addVar(obj=self.c * self.G.dist[i, j]/len(self.scenarios), vtype=GRB.BINARY, name='x_%d_%d_%d' % (j, i, sc))
        for i in range(self.G.n):
            for sc in range(len(self.scenarios)):
                self.g[i, sc] = self.model.addVar(vtype = GRB.BINARY, name="g_%d_%d" % (i, sc))
    
        for l in range(len(self.items)):
            for sc in range(len(self.scenarios)):
                self.y[self.items[l].k, l, sc] = self.model.addVar(obj=-self.items[l].price/len(self.scenarios), vtype=GRB.CONTINUOUS, name='y_%d_%d_%d' % (self.items[l].k, l, sc))
                self.a[l, sc] = self.model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%d_%d' % (l, sc))
            for s in range(self.S):
                self.z[s + self.K, l] = self.model.addVar(vtype=GRB.CONTINUOUS, lb=self.items[l].lb, ub=self.items[l].price, name='z_%d_%d' % (s + self.K, l)) #VARIABLE PRICE, NOT DISCOUNT
                for sc in range(len(self.scenarios)):
                    self.y[s + self.K, l, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='y_%d_%d_%d' % (s +self. K, l, sc))
                    self.w[s + self.K, l, sc] = self.model.addVar(obj=-1/len(self.scenarios), vtype=GRB.CONTINUOUS, name='w_%d_%d_%d' % (s + self.K, l, sc))
                    self.b[s + self.K, l, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='b_%d_%d_%d' % (s + self.K, l, sc))
        for k in range(1, self.K):
            for s in range(self.S):
                for sc in range(len(self.scenarios)):
                    self.p[k, s + self.K, sc] = self.model.addVar(vtype=GRB.CONTINUOUS, name='p_%d_%d_%d' % (k, s + self.K, sc))
    
        self.model.update()
        

        self.model.addConstrs(self.rCost[sc] == quicksum(self.c * self.G.dist[i, j] * (self.x[i, j, sc] + self.x[j, i, sc]) for i, j in self.G.dist.keys()) for sc in range(len(self.scenarios)))
        #balance
        self.model.addConstrs(quicksum(self.x[0, j, sc] for j in range(1, self.G.n)) == 1 for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[j, 0, sc] for j in range(1, self.G.n)) == 1 for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[i, j, sc] for j in range(i + 1, self.G.n)) + quicksum(self.x[i, j, sc] for j in range(0, i)) == self.g[i, sc] for i in range(1, self.G.n) for sc in range(len(self.scenarios)))
        self.model.addConstrs(quicksum(self.x[j, i, sc] for j in range(i + 1, self.G.n)) + quicksum(self.x[j, i, sc] for j in range(0, i)) == self.g[i, sc] for i in range(1, self.G.n) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[self.items[l].k, sc] >= self.y[self.items[l].k, l, sc] for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.g[s + self.K, sc] >= self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))   
        #duality
        self.uConstrs1 = self.model.addConstrs(quicksum(self.a[l, sc] for l in self.Lk[k - 1]) == quicksum(self.items[l].ul - self.items[l].price * self.y[self.items[l].k, l, sc] - quicksum(self.w[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(self.scenarios[sc][k, s + self.K] * self.p[k, s + self.K, sc] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.y[self.items[l].k, l, sc] + quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) == 1 for l in range(len(self.items)) for sc in range(len(self.scenarios)))     
        self.model.addConstrs(self.y[s + self.K, l, sc] <= self.p[k, s + self.K, sc] for k in range(1, self.K) for l in self.Lk[k - 1] for s in range(self.S) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.a[l, sc] + self.b[s + self.K, l, sc] >= self.items[l].ul - self.z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.uConstrs2 = self.model.addConstrs(quicksum(self.b[s + self.K, l, sc] for l in self.Lk[k - 1]) <= self.scenarios[sc][k, s + self.K] for s in range(self.S) for k in range(1, self.K) for sc in range(len(self.scenarios)))
    
        #linearization
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.items[l].price * self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] <= self.z[s + self.K, l] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.model.addConstrs(self.w[s + self.K, l, sc] >= self.z[s + self.K, l] - 1.1*self.items[l].price*(1 - self.y[s +self. K, l, sc]) for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))

        self.model._x = self.x
        self.model._g = self.g
        self.model._v = self.v
        self.model._y = self.y
        self.model._G = self.G
        self.model._R = self.R
        self.model._nscen = len(self.scenarios)
        self.model.Params.lazyConstraints = 1

    def primalTest(self, z, I):
    
        self.model._x = self.x
        self.model._g = self.g
        self.model._v = self.v
        self.model._y = self.y
        self.model._G = self.G
        self.model._R = self.R
        self.model._nscen = len(self.scenarios)
        self.model.Params.lazyConstraints = 1

    def generateScenarios(self, n):
        np.random.seed(self.iter + 1)
        self.scenarios = [{} for i in range(n)]
        for i in range(n):
            for k in range(1, self.K):
                for s in range(self.S):
                    self.scenarios[i][k, s + self.K] = np.random.uniform(self.I_range[k - 1][s][0], self.I_range[k - 1][s][1])

    def setupStochastic(self):

        #init main problem
        self.s_model = Model("s_bl5")
        self.s_model.setParam("OutputFlag", 0)
        self.s_model.setParam("NonConvex", 2)
        self.s_model.setParam("TimeLimit", 3 * 60)
        self.y  ={} #follower's decision (continuous for duality)
        self.z = {} #leader's discount
        self.w = {} #linerization
        self.a = {} #follower's dual 1
        self.b = {} #follower's dual 2
        self.p = {} #follower visits i
        self.g = {}

        self.rCost = self.s_model.addVar(obj=1, vtype=GRB.CONTINUOUS, name="rCost")
        self.avgDists = self.s_model.addVars(len(self.scenarios), vtype=GRB.CONTINUOUS, name="avgDist") #for all scenarios
        self.rCostPt2_raw = self.s_model.addVars(len(self.scenarios), vtype=GRB.CONTINUOUS, name="rCostPt2_raw")
        self.rCostPt2_sqrt = self.s_model.addVars(len(self.scenarios), vtype=GRB.CONTINUOUS, name="rCostPt2_sqrt")
        for i in range(1, self.G.n):
            for sc in range(len(self.scenarios)):
                self.g[i, sc] = self.s_model.addVar(vtype = GRB.BINARY, name="g_%d_%d" % (i, sc))
    
        for l in range(len(self.items)):
            for sc in range(len(self.scenarios)):
                self.y[self.items[l].k, l, sc] = self.s_model.addVar(obj=-self.items[l].price/len(self.scenarios), vtype=GRB.BINARY, name='y_%d_%d_%d' % (self.items[l].k, l, sc))
                self.a[l, sc] = self.s_model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%d_%d' % (l, sc))
            for s in range(self.S):
                self.z[s + self.K, l] = self.s_model.addVar(obj = 0, vtype=GRB.CONTINUOUS, lb=self.items[l].lb, ub=self.items[l].price, name='z_%d_%d' % (s + self.K, l)) #VARIABLE PRICE, NOT DISCOUNT
                for sc in range(len(self.scenarios)):
                    self.y[s + self.K, l, sc] = self.s_model.addVar(vtype=GRB.BINARY, name='y_%d_%d_%d' % (s +self. K, l, sc))
                    self.w[s + self.K, l, sc] = self.s_model.addVar(obj=-1/len(self.scenarios), vtype=GRB.CONTINUOUS, name='w_%d_%d_%d' % (s + self.K, l, sc))
                    self.b[s + self.K, l, sc] = self.s_model.addVar(vtype=GRB.CONTINUOUS, name='b_%d_%d_%d' % (s + self.K, l, sc))
        for k in range(1, self.K):
            for s in range(self.S):
                for sc in range(len(self.scenarios)):
                    self.p[k, s + self.K, sc] = self.s_model.addVar(vtype=GRB.CONTINUOUS, name='p_%d_%d_%d' % (k, s + self.K, sc))
    
        self.s_model.update()
    
    
        self.s_model.addConstr(self.rCost == 1/len(self.scenarios) * sum((2 * self.avgDists[sc] * sum(item.w for item in self.items)/self.G.q + 0.57 * self.rCostPt2_sqrt[sc]) for sc in range(len(self.scenarios))))
        self.s_model.addConstrs(self.avgDists[sc] * (sum(self.g[k, sc] for k in range(1, self.K)) + sum(self.g[i, sc] for i in range(self.K, self.K + self.S))) == (sum(self.G.dist[k, 0] * self.g[k, sc] for k in range(1, self.K)) + sum(self.G.dist[i, 0] * self.g[i, sc] for i in range(self.K, self.K + self.S))) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.rCostPt2_raw[sc] == self.G.A * (sum(self.g[k, sc] for k in range(1, self.K)) + sum(self.g[i, sc] for i in range(self.K, self.K + self.S))) for sc in range(len(self.scenarios)))
        for sc in range(len(self.scenarios)):
            self.s_model.addGenConstrPow(self.rCostPt2_raw[sc], self.rCostPt2_sqrt[sc], 0.5)

        #balance
        self.s_model.addConstrs(self.g[self.items[l].k, sc] >= self.y[self.items[l].k, l, sc] for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.g[s + self.K, sc] >= self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        
        #duality
        self.uConstrs1s = self.s_model.addConstrs(quicksum(self.a[l, sc] for l in self.Lk[k - 1]) == quicksum(self.items[l].ul * self.y[self.items[l].k, l, sc] - self.items[l].price * self.y[self.items[l].k, l, sc] + quicksum(self.items[l].ul * self.y[s + self.K, l, sc] - self.w[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(self.scenarios[sc][k, s + self.K] * self.p[k, s + self.K, sc] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.y[self.items[l].k, l, sc] + quicksum(self.y[s + self.K, l, sc] for s in range(self.S)) == 1 for l in range(len(self.items)) for sc in range(len(self.scenarios)))     
        self.s_model.addConstrs(self.y[s + self.K, l, sc] <= self.p[k, s + self.K, sc] for k in range(1, self.K) for l in self.Lk[k - 1] for s in range(self.S) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.a[l, sc] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.a[l, sc] + self.b[s + self.K, l, sc] >= self.items[l].ul - self.z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.uConstrs2s = self.s_model.addConstrs(quicksum(self.b[s + self.K, l, sc] for l in self.Lk[k - 1]) <= self.scenarios[sc][k, s + self.K] for s in range(self.S) for k in range(1, self.K) for sc in range(len(self.scenarios)))
    
        #linearization
        self.s_model.addConstrs(self.w[s + self.K, l, sc] <= self.items[l].price * self.y[s + self.K, l, sc] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.w[s + self.K, l, sc] <= self.z[s + self.K, l] for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))
        self.s_model.addConstrs(self.w[s + self.K, l, sc] >= self.z[s + self.K, l] - 1.1*self.items[l].price*(1 - self.y[s +self. K, l, sc]) for s in range(self.S) for l in range(len(self.items)) for sc in range(len(self.scenarios)))

    def primalTest(self, z, I):

        self.m1 = Model("pTest")
        self.m1.setParam("OutputFlag", 0)
        self.y1 = {} #follower's decision (continuous for duality)
        self.w1 = {} #linerization
        self.p1 = {} #follower visits i
        self.a1 = {}
        self.b1 = {}

        for l in range(len(self.items)):
            self.y1[self.items[l].k, l] = self.m1.addVar(vtype=GRB.CONTINUOUS, name='y_%g_%g' % (self.items[l].k, l))
            self.a1[l] = self.m1.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%g' % l)
            for s in range(self.S):
                self.y1[s + self.K, l] = self.m1.addVar(vtype=GRB.CONTINUOUS, name='y_%g_%g' % (s + self.K, l))
                self.b1[s + self.K, l] = self.m1.addVar(vtype=GRB.CONTINUOUS, name='b_%g_%g' % (s + self.K, l))
        for k in range(1, self.K):
            for s in range(self.S):
                self.p1[k, s + self.K] = self.m1.addVar(vtype=GRB.CONTINUOUS, name='p_%g_%g' % (k, s + self.K))

        self.m1.update()

        # self.m1.setObjective(sum(quicksum(self.items[l].ul * self.y1[self.items[l].k, l] - self.items[l].price * self.y1[self.items[l].k, l] + quicksum(self.items[l].ul*self.y1[s + self.K, l] - z[s + self.K, l] * self.y1[s + self.K, l] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(I[k - 1][s] * self.p1[k, s + self.K] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y1[s + self.K, l] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K)))
        # self.m1.ModelSense = -1

        self.m1.addConstrs(quicksum(self.a1[l] for l in self.Lk[k - 1]) == quicksum((self.items[l].ul - self.items[l].price) * self.y1[self.items[l].k, l] + quicksum((self.items[l].ul - z[s + self.K, l]) * self.y1[s + self.K, l] for s in range(self.S)) for l in self.Lk[k - 1]) - quicksum(self.uiT[k - 1][s] * self.p1[k, s + self.K] for s in range(self.S)) + quicksum(self.items[l].inc * quicksum(self.y1[s + self.K, l] for s in range(self.S)) for l in self.Lk[k - 1]) for k in range(1, self.K))
        self.m1.addConstrs(self.y1[self.items[l].k, l] + quicksum(self.y1[s + self.K, l] for s in range(self.S)) == 1 for l in range(len(self.items)))     
        self.m1.addConstrs(self.y1[s + self.K, l] <= self.p1[k, s + self.K] for k in range(1, self.K) for l in self.Lk[k - 1] for s in range(self.S))
        self.m1.addConstrs(self.a1[l] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)))
        self.m1.addConstrs(self.a1[l] + self.b1[s + self.K, l]  >= self.items[l].ul - z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)))
        self.m1.addConstrs(quicksum(self.b1[s + self.K, l] for l in self.Lk[k - 1]) <= self.uiT[k - 1][s] for s in range(self.S) for k in range(1, self.K))
        
        self.m1.optimize()

        # print("-----Y's------")
        # for l in range(len(self.items)):
        #     print(self.y1[self.items[l].k, l].x)
        #     for s in range(self.S):
        #         print(self.y1[s + self.K, l].x)
        # print("------P's------")
        # for k in range(1, self.K):
        #     for s in range(self.S):
        #         print(self.p1[k, s + self.K].x)

    def dualTest(self, z, I):

        self.m2 = Model("dtest")

        self.a1 = {} #follower's dual 1
        self.b1 = {} #follower's dual 2

        for l in range(len(self.items)):
            self.a1[l] = self.m2.addVar(obj=1, lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%g' % l)
            for s in range(self.S):
                self.b1[s + self.K, l] = self.m2.addVar(vtype=GRB.CONTINUOUS, name='b_%g_%g' % (s + self.K, l))

        self.m2.update()

        self.m2.addConstrs(self.a1[l] >= self.items[l].ul - self.items[l].price for l in range(len(self.items)))
        self.m2.addConstrs(self.a1[l] + self.b1[s + self.K, l] >= self.items[l].ul - z[s + self.K, l] + self.items[l].inc for s in range(self.S) for l in range(len(self.items)))
        self.m2.addConstrs(quicksum(self.b1[s + self.K, l] for l in self.Lk[k - 1]) <= self.uiT[k - 1][s] for s in range(self.S) for k in range(1, self.K))

        self.m2.optimize()

    def setupSub(self):

        self.sub = Model("Sub")
        self.sub.setParam("OutputFlag", 0)
        self.I_var1 = {}
        self.I_var2 = {}
        for k in range(1, self.K):
            for s in range(self.S):
                self.I_var1[k, s + self.K] = self.sub.addVar(obj=1, lb=self.I_range[k - 1][s][0], ub=self.I_range[k - 1][s][1], name="I1_%d,%d" % (k, s))
                self.I_var2[k, s + self.K] = self.sub.addVar(obj=-1, lb=self.I_range[k - 1][s][0], ub=self.I_range[k - 1][s][1], name="I2_%d,%d" % (k, s))
        self.sub.update()
        for k in range(1, self.K):
            for s in range(self.S):
                self.sub.addConstr(self.I_var1[k, s + self.K] <= self.I_var2[k, s + self.K])

    def refineUncrtS(self, niter=3, nscen=100):
        self.iter = 0
        self.generateScenarios(nscen)
        # print(all(element == self.scenarios[0] for element in self.scenarios))
        # self.setupStochastic()
        self.setupMaster()
        # self.setupTSP2()
        # self.model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\vrp_k10.lp")
        self.setupSub()
        xmax = 1.07 * np.max(self.I_range)
        change = []
        uncrt = []
        fig,ax = plt.subplots(figsize=(8, 5))
        colors = ["#01949A", "#004369", "#DB1F48"]
        for s in range(self.S):
                plt.plot((0, 0), (0, 0), '|-', markersize=0.01, color=colors[s], label="store %d" % (s + 1))
        for k in range(1, self.K):
            for s in range(self.S):
                plt.plot((self.I_range[k - 1][s][0], self.I_range[k - 1][s][1]), (k + 0.1 * s, k + 0.1 * s), '|-', color=colors[s])
                for sc in range(len(self.scenarios)):
                    plt.vlines(x=self.scenarios[sc][k, s + self.K], ymin=k + 0.1 * s - 0.07, ymax=k + 0.1 * s + 0.07, color="black", linewidth=0.5)
                plt.plot(self.uiT[k - 1][s], k + 0.1 * s, "o", color=colors[s], markersize=3)
        plt.yticks(range(1, self.K))
        plt.xlim(0, xmax)
        plt.ylabel("Customer")
        plt.xlabel("Uncertainty range")
        plt.title("Round 0")
        plt.legend()
        plt.savefig("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\papers\\img\\uncrt\\k10_s3_sc5\\range_iter0.png", dpi=500)
        # plt.savefig("img/stoch/k%d_s%d_sc%d_i%d/range_iter0.png" %(self.K, self.S, len(self.scenarios), niter), dpi=500)
        for iter in range(niter):
            self.iter +=1
            plt.clf()
            start = time.time()
            # self.model.optimize(callbackTSP2)
            self.model.optimize(callback_mult)
            # if iter == 0:
            #     self.printOneScen()
            #     self.model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\vrp_k10.lp")
            print("Stoch obj: %g" % self.model.objVal)
            print("Solved stochastic: %g sec" % (time.time() - start))
            print([self.rCost[sc].x for sc in range(len(self.scenarios))])
            # print("TSP route:")
            # for i in range(self.G.n):
            #     for j in range(i):
            #         if self.x[i, j, 0].x > 0.5:
            #             print("%d, %d" % (i, j))
            #     for j in range(i + 1, self.G.n):
            #         if self.x[i, j, 0].x > 0.5:
            #             print("%d, %d" % (i, j))
            # print("Stoch obj: %g" % self.s_model.objVal)
            # print(self.rCost.x, self.avgDists[0].x, self.rCostPt2_raw[0].x, self.rCostPt2_sqrt[0].x)
            # for i in range(1, self.G.n):
            #     print(self.g[i, 0].x, self.y[i, i - 1, 0].x)
            print("Solved stochastic: gap %g" % self.model.MIPGap)
            # self.printOneScen()
            #print([[self.items[l].price, [self.z[s + self.K, l].x for s in range(self.S)]] for l in range(len(self.items))])
            z = {}
            for l in range(len(self.items)):
                for s in range(self.S):
                    z[s + self.K, l] = self.z[s + self.K, l].x
            # y_det = bilevel_v5_approx.getModel(self.G, self.items, self.Lk, self.uiT, self.S, z)
            y_det = bilevel_v5_full.getModel(self.G, self.items, self.Lk, self.uiT, self.S, self.G.r, self.q, self.c, z)
            # y_det2 = bilevel_v5_tsp.getModel(self.G, self.items, self.Lk, self.uiT, self.S, self.G.r, self.q, self.c, z)
            allSame = self.allRespSame()
            print("All scenario responces same ", allSame)
            # self.printOneScen(0)
            # self.printOneScen(69)
            # detDiff = False
            # for key in y_det.keys():
            #     if abs(y_det[key] - y_det2[key]) > 0.1:
            #         detDiff = True
            #         break
            # print("VRP and TSP same sol: ", (not detDiff))
            if allSame:
                sol_diff = 0
                sol_diff2 = 0
                for key in y_det.keys():
                    for sc in range(len(self.scenarios)):
                        sol_diff += abs(self.y[key[0], key[1], sc].x - y_det[key])
                        # sol_diff2 += abs(self.y[key[0], key[1], sc].x - y_det2[key])
                print("True choice coincides (VRP) ", (sol_diff < 1))
                # print("True choice coincides (TSP) ", (sol_diff2 < 1))
                if (sol_diff < 1):
                    break
            for l in range(len(self.items)):
                if y_det[self.items[l].k, l] > 0.5:
                    j = self.items[l].k
                    for s in range(self.S):
                        self.sub.addConstr(self.items[l].price <= self.z[s + self.K, l].x + self.I_var1[self.items[l].k, s + self.K] - self.items[l].inc)
                        self.sub.addConstr(self.items[l].price <= self.z[s + self.K, l].x + self.I_var2[self.items[l].k, s + self.K] - self.items[l].inc)
                else:
                    for s in range(self.S):
                        if y_det[s + self.K, l] > 0.5:
                            j = s + self.K
                            #better than home
                            self.sub.addConstr(self.items[l].price >= self.z[j, l].x + self.I_var1[self.items[l].k, j] - self.items[l].inc)
                            self.sub.addConstr(self.items[l].price >= self.z[j, l].x + self.I_var2[self.items[l].k, j] - self.items[l].inc)
                            #better than other stores
                            for s1 in range(s):
                                self.sub.addConstr(self.z[j, l].x + self.I_var1[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var1[self.items[l].k, s1 + self.K])
                                self.sub.addConstr(self.z[j, l].x + self.I_var1[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var2[self.items[l].k, s1 + self.K])
                                #self.sub.addConstr(self.z[j, l].x + self.I_var2[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var1[self.items[l].k, s1 + self.K])
                                self.sub.addConstr(self.z[j, l].x + self.I_var2[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var2[self.items[l].k, s1 + self.K])
                            for s1 in range(s + 1, self.S):
                                self.sub.addConstr(self.z[j, l].x + self.I_var1[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var1[self.items[l].k, s1 + self.K])
                                self.sub.addConstr(self.z[j, l].x + self.I_var1[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var2[self.items[l].k, s1 + self.K])
                                #self.sub.addConstr(self.z[j, l].x + self.I_var2[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var1[self.items[l].k, s1 + self.K])
                                self.sub.addConstr(self.z[j, l].x + self.I_var2[self.items[l].k, j] <= self.z[s1 + self.K, l].x + self.I_var2[self.items[l].k, s1 + self.K])
            self.sub.optimize()
            #self.sub.write("sub_%d.lp" % iter)
            change.append(0)
            uncrt.append(sum([self.I_range[k - 1][s][1] - self.I_range[k - 1][s][0] for k in range(1, self.K) for s in range(self.S)]))
            for k in range(1, self.K):
                for s in range(self.S):
                    if self.I_range[k - 1][s][0] != self.I_var1[k, s + self.K].x:
                        print("LB change: cust %d, store %d - %g, %g; true %g" % (k, s + self.K, self.I_range[k - 1][s][0], self.I_var1[k, s + self.K].x, self.uiT[k - 1][s]))
                        change[-1] += self.I_var1[k, s + self.K].x - self.I_range[k - 1][s][0]
                        self.I_range[k - 1][s][0] = self.I_var1[k, s + self.K].x
                        #self.I_var1.lb = self.I_var1[k, s + self.K].x
                    if self.I_range[k - 1][s][1] != self.I_var2[k, s + self.K].x:
                        print("UB change: cust %d, store %d - %g, %g; true %g" % (k, s + self.K, self.I_range[k - 1][s][1], self.I_var2[k, s + self.K].x, self.uiT[k - 1][s]))
                        change[-1] += self.I_range[k - 1][s][1] - self.I_var2[k, s + self.K].x
                        self.I_range[k - 1][s][1] = self.I_var2[k, s + self.K].x
            self.old_sc = self.scenarios.copy()
            self.generateScenarios(nscen)
            print("Done generating scenarios")
            for s in range(self.S):
                plt.plot((0, 0), (0, 0), '|-', markersize=0.01, color=colors[s], label="store %d" % (s + 1))
            for k in range(1, self.K):
                for s in range(self.S):
                    plt.plot((self.I_range[k - 1][s][0], self.I_range[k - 1][s][1]), (k + 0.1 * s, k + 0.1 * s), '|-', color=colors[s])
                    for sc in range(len(self.scenarios)):
                        plt.vlines(x=self.scenarios[sc][k, s + self.K], ymin=k + 0.1 * s - 0.07, ymax=k + 0.1 * s + 0.07, color="black", linewidth=0.5)
                    plt.plot(self.uiT[k - 1][s], k + 0.1 * s, "o", color=colors[s], markersize=3)
            plt.yticks(range(1, self.K))
            plt.xlim(0, xmax)
            plt.ylabel("Customer")
            plt.xlabel("Uncertainty range")
            plt.title("Round %d" % (iter + 1))
            plt.legend()
            plt.savefig("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\papers\\img\\uncrt\\k10_s3_sc5\\range_iter%d.png" % (iter + 1), dpi=500)
            # plt.savefig("img/stoch/k%d_s%d_sc%d_i%d/range_iter%d.png" %(self.K, self.S, len(self.scenarios), niter, iter + 1), dpi=500)
            print("----------------------------------")
            for k in range(1, self.K):
                for s in range(self.S):
                    for sc in range(len(self.scenarios)):
                        self.model.chgCoeff(self.uConstrs1[k, sc], self.p[k, s + self.K, sc], self.scenarios[sc][k, s + self.K])
                        self.uConstrs2[s, k, sc].rhs = self.scenarios[sc][k, s + self.K]
            print("done updating constrs")
        plt.clf()
        plt.plot(change, color=colors[0])
        plt.title("Learned per round")
        plt.savefig("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\papers\\img\\uncrt\\k10_s3_sc5\\learned.png")
        # plt.savefig("img/stoch/k%d_s%d_sc%d_i%d/learned.png" %(self.K, self.S, len(self.scenarios), niter), dpi=500)
        plt.clf()
        plt.plot(uncrt, color=colors[2])
        plt.title("Total uncertainty")
        plt.savefig("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\papers\\img\\uncrt\\k10_s3_sc5\\uncrt.png")
        # plt.savefig("img/stoch/k%d_s%d_sc%d_i%d/uncrt.png" %(self.K, self.S, len(self.scenarios), niter), dpi=500)
        print(self.I_range)
        print(self.uiT)
        self.printOneScen()
    
    def allRespSame(self):
        res = True
        for key in self.y.keys():
            val = self.y[key[0], key[1], 0].x
            if self.y[key].x - val > 0.01:
                res = False
                print("Scenario %d different from scenario 0" % key[2])
                break
        return res

    
    def printTrueChoice(self):
        for k in range(1, self.K):
            print("Customer %d obj: %g" % (k, sum(self.items[l].ul * self.y1[self.items[l].k, l].x - self.items[l].price * self.y1[self.items[l].k, l].x 
                  + sum(self.items[l].ul * self.y1[s + self.K, l].x - self.z[s + self.K, l].x * self.y1[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]) 
                  - sum(self.uiT[k - 1][s] * self.p1[k, s + self.K].x for s in range(self.S)) + sum(self.items[l].inc * sum(self.y1[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1])))
            for l in self.Lk[k - 1]:
                print("%d, %d: %g %g %g"  % (k, self.items[l].type, self.y1[self.items[l].k, l].x, self.items[l].ul, self.items[l].price))
                for s in range(self.S):
                    print("%d, %d, %d: %g %g %g %g %g"  % (s + self.K, self.items[l].type, l, self.y1[s + self.K, l].x, self.uiT[k - 1][s], self.z[s + self.K, l].x, self.items[l].ul - self.z[s + self.K, l].x - self.uiT[k - 1][s] + self.items[l].inc, self.items[l].inc))

    def printOneScen(self, sc=0):
        for k in range(1, self.K):
            print("Customer %d obj: %g, %g" % (k, sum((self.items[l].ul - self.items[l].price) * self.y[self.items[l].k, l, sc].x 
                  + sum((self.items[l].ul - self.z[s + self.K, l].x) * self.y[s + self.K, l, sc].x for s in range(self.S)) for l in self.Lk[k - 1]) 
                  - sum(self.old_sc[sc][k, s + self.K] * self.p[k, s + self.K, sc].x for s in range(self.S)) 
                  + sum(self.items[l].inc * sum(self.y[s + self.K, l, sc].x for s in range(self.S)) for l in self.Lk[k - 1]), sum(self.a[l, sc].x for l in self.Lk[k - 1])))
            for l in self.Lk[k - 1]:
                print("%d, %d: %g %g %g"  % (k, self.items[l].type, self.y[self.items[l].k, l, sc].x, self.items[l].ul, self.items[l].price))
                for s in range(self.S):
                    print("%d, %d, %d: %g %g %g %g %g"  % (s + self.K, self.items[l].type, l, self.y[s + self.K, l, sc].x, self.old_sc[sc][k, s + self.K], self.z[s + self.K, l].x, self.items[l].ul - self.z[s + self.K, l].x - self.old_sc[sc][k, s + self.K] + self.items[l].inc, self.items[l].inc))


    def printSol(self):
        print("Leader's objective: %g" % -self.model.objVal)
        print("Customers decisions:")
        for k in range(1, self.K):
            print("Customer %d - Obj: %g   part1: %g   part2: %g" % (k, sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.items[l].price * self.y[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]) 
            - sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S))
            + sum(self.items[l].inc * sum(self.y[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]), sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.items[l].price * self.y[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]), 
            sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S))))
            for l in self.Lk[k - 1]:
                print("%d, %d: %g %g %g"  % (k, self.items[l].type, self.y[self.items[l].k, l].x, self.items[l].ul, self.items[l].price))
                for s in range(self.S):
                    print("%d, %d, %d: %g %g %g %g %g"  % (s + self.K, self.items[l].type, l, self.y[s + self.K, l].x,  self.items[l].ul, self.ui[k - 1][s], self.z[s + self.K, l].x, self.items[l].price))

        print("Routing:")
        print([self.model.getVarByName("C_%g" % r).x for r in range(self.G.r)])
        for r in range(self.G.r):
            items = []
            for l in range(len(self.items)):
                if self.model.getVarByName("v_%g_%g_%g" % (self.items[l].k, l, r)).x > 0.5:
                    items.append([l, self.items[l].k])
                for s in range(self.S):
                    if self.model.getVarByName("v_%g_%g_%g" % (s + self.K, l, r)).x > 0.5:
                        items.append([l, self.items[l].k, s + self.K])
            print(items)
        for i, j, r in self.x.keys():
            if self.x[i, j, r].x > 0.5:
                print("%d, %d, %d" % (i, j, r))
        print("Routing cost: %g" % self.model.getVarByName("rCost").x)#sum([self.model.getVarByName("C_%g").x % r for r in range(self.G.r)]))
        #print("Revenue: %g" % self.model.getVarByName("rev").x)
    