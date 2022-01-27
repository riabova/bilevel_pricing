from gurobipy import *
from numpy.core import numeric
import Graph
import random
import numpy as np
from itertools import combinations
import matplotlib.pyplot as plt

def callback_mult(model, where):
    if where == GRB.Callback.MIPSOL:
        G = model._G
        x_sol = model.cbGetSolution(model._x)
        for k in range(model._R):
            callback(model, G, k, x_sol)

def callback(model, G, k, x_sol):
    adj = {}
    for i in range(G.n):
        adj[i] = []
    for i in range(1, G.n):
        for j in range(1, i):
            if x_sol[i, j, k] > 0.5:
                adj[i].append(j)
                adj[j].append(i)
        for j in range(i+1, G.n):
            if x_sol[i, j, k] > 0.5:
                adj[i].append(j)
                adj[j].append(i)
    cc = DFS(G, adj)
    if len(cc) > 1:
        for comp in cc:
            model.cbLazy(quicksum(model._x[i, j, k] + model._x[j, i, k] for i, j in combinations(comp, 2)) <= len(comp) - 1)


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

def getModel(G: Graph.Graph, S: int, L: int, R: int, ul: list(), ui: list, h: list(), q: list()):
    model = Model("bl2")
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    x = {} #leader's routing
    y  ={} #followe's decision (continuous for duality)
    z = {} #leader's discount
    v = {} #leader's product-location/vehicle mapping
    w = {} #linerization
    a = {} #follower dual
    M1 = 1000000
    M2 = 1000

    for i in range(G.n):
        for r in range(R):
            for j in range(i):
                x[i, j, r] = model.addVar(obj=G.dist[i, j], vtype=GRB.BINARY, name="x_%g_%g_%g" % (i, j, r))
            for j in range(i+1, G.n):
                x[i, j, r] = model.addVar(obj=G.dist[j, i], vtype=GRB.BINARY, name="x_%g_%g_%g" % (i, j, r))
    for i in range(1, G.n):
        for r in range(R):
            for l in range(L):
                v[i, l, r] = model.addVar(vtype=GRB.BINARY, name="v_%g_%g_%g" % (i, l, r))

    for k in range(1, K):
        for l in range(L):
            y[k, k, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g_%g' % (k, k, l))
            y[k, 0, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g_%g' % (k, 0, l))
            z[k, k, l] = model.addVar(vtype=GRB.CONTINUOUS, lb=30, name='z_%g_%g_%g' % (k, k, l))
            w[k, k, l] = model.addVar(obj=-1, vtype=GRB.CONTINUOUS, name='w_%g_%g_%g' % (k, k, l))
            a[k, l] = model.addVar(vtype=GRB.CONTINUOUS, name='a_%g_%g' % (k, l))
            for s in range(S):
                y[k, s + K, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g_%g' % (k, s + K, l))
                z[k, s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, lb=30, name='z_%g_%g_%g' % (k, s + K, l))
                w[k, s + K, l] = model.addVar(obj=-1, vtype=GRB.CONTINUOUS, name='w_%g_%g_%g' % (k, s + K, l))

    model.update()

    '''for r in range(R):
        x[0, 7, r].obj = 1
        x[7, 0, r].obj = 1'''

    #balance
    for r in range(R):
        model.addConstr(quicksum(x[0, j, r] for j in range(1, G.n)) <= 1)
        model.addConstr(quicksum(x[j, 0, r] for j in range(1, G.n)) <= 1)
        model.addConstr(quicksum(x[0, j, r] for j in range(1, G.n)) - quicksum(x[j, 0, r] for j in range(1, G.n)) == 0)
        for i in range(1, G.n):
            model.addConstr(quicksum(x[i, j, r] for j in range(i)) + quicksum(x[i, j, r] for j in range(i + 1, G.n)) >= quicksum(v[i, l, r] for l in range(L))/M1)
            model.addConstr(quicksum(x[i, j, r] for j in range(i)) + quicksum(x[i, j, r] for j in range(i + 1, G.n)) 
            - (quicksum(x[j, i, r] for j in range(i)) + quicksum(x[j, i, r] for j in range(i + 1, G.n))) == 0)

    #capacity
    model.addConstrs(quicksum(v[k, l, r] for r in range(R)) == y[k, k, l] for k in range(1, K) for l in range(L))
    model.addConstrs(quicksum(v[s + K, l, r] for r in range(R)) == quicksum(y[k, s + K, l] for k in range(1, K))  for s in range(S) for l in range(L))
    model.addConstrs(quicksum(h[l]*quicksum(v[i, l, r] for i in range(1, G.n)) for l in range(L)) <= q[r] for r in range(R))

    #duality
    model.addConstrs(a[k, l] == ul[k, l]*y[k, k, l] - w[k, k, l] + quicksum((ul[k, l] - ui[k, s])*y[k, s + K, l] - w[k, s + K, l] for s in range(S)) for l in range(L) for k in range(1, K))
    model.addConstrs(y[k, 0, l] + y[k, k, l] + quicksum(y[k, s + K, l] for s in range(S)) == 1 for l in range(L) for k in range(1, K))
    model.addConstrs(a[k, l] >= 0 for l in range(L) for k in range(1, K)) #utility of non-buying
    model.addConstrs(a[k, l] >= ul[k, l] - z[k, k, l] for l in range(L) for k in range(1, K))
    model.addConstrs(a[k, l] >= ul[k, l] - ui[k, s] - z[k, s + K, l] for s in range(S) for l in range(L) for k in range(1, K))

    #linearization
    model.addConstrs(w[k, k, l] <= z[k, k, l] for l in range(L) for k in range(1, K))
    model.addConstrs(w[k, s + K, l] <= z[k, s + K, l] for l in range(L) for s in range(S) for k in range(1, K))
    model.addConstrs(w[k, k, l] <= M2*y[k, k, l] for l in range(L) for k in range(1, K))
    model.addConstrs(w[k, s + K, l] <= M2*y[k, s + K, l] for s in range(S) for l in range(L) for k in range(1, K))
    model.addConstrs(w[k, k, l] >= z[k, k, l] - M2*(1 - y[k, k, l]) for l in range(L) for k in range(1, K))
    model.addConstrs(w[k, s + K, l] >= z[k, s + K, l] - M2*(1 - y[k, s + K, l]) for s in range(S) for l in range(L) for k in range(1, K))

    '''for k in range(1, K):
        for l in range(L):
            model.addConstr(y[k, 0, l] == 1)'''
    #model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\test1.lp")

    model._x = x
    model._G = G
    model._R = R
    model.Params.lazyConstraints = 1

    #model.optimize(callback_mult)

    '''for i, j, r in x.keys():
        if x[i, j, r].x > 0.5:
            print("%d, %d, %d" % (i, j, r))

    print("//////////////////////////")
    for k, i, l in y.keys():
        if i > 0:
            print("%d, %d, %d: %d %d %d"  % (k, i, l, np.round(z[k, i, l].x), np.round(w[k, i, l].x), np.round(y[k, i, l].x)))
        else:
            print("%d, %d, %d: %d %d %d" % (k, i, l, -1, -1, y[k, i, l].x))
        #print("---------------------")'''

    return model, y, z, K, L

'''G = Graph.Graph()
G.read("..\..\data\TSP_instance_n_10_s_1.dat")
S = 3
R = 2
q = [100]*R
L = 4
h = [3, 40, 5, 6]
ui = {}
ul = {}
f = open("..\\..\\data\\ui_n_10_s_1.txt", 'r')
f.readline()
k = 0
for line in f:
    k += 1
    utls = line.split()
    for i in range(S):
        ui[k, i] = int(utls[i])
f.close()
f = open("..\\..\\data\\ul_n_10_s_1.txt", 'r')
f.readline()
k = 0
for line in f:
    k += 1
    utls = line.split()
    for l in range(L):
        ul[k, l] = int(utls[l])
f.close()
mtz = Bilevel_v2(G, S, L, R, ul, ui, h, q)'''