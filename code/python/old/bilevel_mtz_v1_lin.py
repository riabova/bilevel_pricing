from gurobipy import *
from numpy.core import numeric
import Graph
import random
import numpy as np
import matplotlib.pyplot as plt

def Bilevel_v1(G: Graph.Graph, S: int, p: float, r: dict()):
    model = Model("bl1")
    #model.params.NonConvex = 2
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    u = {} #mtz subtour variables
    x = {} #leader routing
    y  ={} #follower decision (continuous for duality)
    z = {} #leader discount
    w = {} #linerization
    a = {} #follower dual

    for i, j in G.dist.keys():
        x[i, j] = model.addVar(vtype=GRB.BINARY, name='x_%g_%g' % (i, j))
        x[j, i] = model.addVar(vtype=GRB.BINARY, name='x_%g_%g' % (j, i))
    for i in range(G.n):
        u[i] = model.addVar(vtype=GRB.CONTINUOUS, name="u_%g" % i)

    for k in range(1, K):
        y[k, k] = model.addVar(vtype=GRB.CONTINUOUS, name='y_%g_%g' % (k, k))
        z[k, k] = model.addVar(vtype=GRB.CONTINUOUS, name='z_%g_%g' % (k, k))
        w[k, k] = model.addVar(vtype=GRB.CONTINUOUS, name='w_%g_%g' % (k, k))
        a[k] = model.addVar(vtype=GRB.CONTINUOUS, name='a_%g' % k)
        for s in range(S):
            y[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='y_%g_%g' % (k, s + K))
            z[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='z_%g_%g' % (k, s + k))
            w[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='w_%g_%g' % (k, s + k))

    model.update()

    obj = quicksum(p*y[k, i] - w[k, i] for k, i in y.keys()) - quicksum(G.dist[i, j]*x[i, j] + G.dist[i, j]*x[j, i] for i, j in G.dist.keys())
    model.setObjective(obj, GRB.MAXIMIZE)

    for k in range(1, K): #balance1
        model.addConstr(quicksum(x[j, k] for j in range(0, k)) + quicksum(x[j, k] for j in range(k + 1, G.n)) == y[k, k])
        model.addConstrs(quicksum(x[j, s + K] for j in range(0, s + K)) + quicksum(x[j, s + K] for j in range(s + K + 1, G.n)) >= y[k, s + K] for s in range(S))
    
    #balance2
    model.addConstrs(quicksum(x[i, j] for j in range(0, i)) + quicksum(x[i, j] for j in range(i + 1, G.n)) == quicksum(x[j, i] for j in range(0, i)) + quicksum(x[j, i] for j in range(i + 1, G.n)) for i in range(G.n))
    
    #subtour
    for i, j in x.keys():
        if i > 0 and j > 0:
            model.addConstr(u[i] - u[j] + G.n*x[i, j] <= G.n - 1)

    #zero discount for home (questionable)
    model.addConstrs(z[k, k] == 0 for k in range(1, K))

    #duality
    model.addConstrs(a[k] == w[k, k] - r[k][0]*y[k, k] + quicksum(w[k, s + K] - r[k][s + 1]*y[k, s + K] for s in range(S)) for k in range(1, K))
    model.addConstrs(y[k, k] + quicksum(y[k, s + K] for s in range(S)) == 1 for k in range(1, K))
    model.addConstrs(a[k] >= z[k, k] - r[k][0] for k in range(1, K))
    model.addConstrs(a[k] >= z[k, s + K] - r[k][s + 1] for s in range(S) for k in range(1, K))

    #linearization
    model.addConstrs(w[k, k] <= z[k, k] for k in range(1, K))
    model.addConstrs(w[k, s + K] <= z[k, s + K] for s in range(S) for k in range(1, K))
    model.addConstrs(w[k, k] <= p*y[k, k] for k in range(1, K))
    model.addConstrs(w[k, s + K] <= p*y[k, s + K] for s in range(S) for k in range(1, K))
    model.addConstrs(w[k, k] >= z[k, k] - p*(1 - y[k, k]) for k in range(1, K))
    model.addConstrs(w[k, s + K] >= z[k, s + K] - p*(1 - y[k, s + K]) for s in range(S) for k in range(1, K))

    model.optimize()

    for i, j in x.keys():
        if x[i, j].x > 0.5:
            print("%d, %d" % (i, j))

    print("//////////////////////////")
    for k, i in y.keys():
        print("%d, %d: %d" % (k, i, y[k, i].x))
        print("%d, %d: %d" % (k, i, z[k, i].x))
        print("---------------------")

G = Graph.Graph()
G.read("..\..\data\TSP_instance_n_10_s_1.dat")
a = np.transpose(G.points)
plt.scatter(a[0], a[1])
for i in range(G.n):
    plt.annotate(i, (a[0][i], a[1][i]))
plt.show()
S = 3
p = 20
#r = {i: [0] + [random.randrange(0, 10) for x in range(S)] for i in range(1, G.n - S)}
#print(r)
r = {i: [0, 1, 2, 3] for i in range(1, G.n - S)}
mtz = Bilevel_v1(G, S, p, r)