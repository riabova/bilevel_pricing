from gurobipy import *
import Graph
from itertools import combinations
import matplotlib.pyplot as plt

def getModel(G: Graph.Graph, items: list(), Lk: list(), ui, S: int, z_bd={}):
    model = Model("bl5")
    model.setParam("NonConvex", 2)
    model.setParam("OutputFlag", 0)
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    y  ={} #follower's decision (continuous for duality)
    z = {} #leader's discount
    w = {} #linerization
    a = {} #follower's dual 1
    b = {} #follower's dual 2
    #d = {} #follower's dual 3
    #gm = {} #follower's dual 4
    #th = {} #follower's dual 5
    p = {} #follower visits i
    g = {}
    
    rCost = model.addVar(obj=1, vtype=GRB.CONTINUOUS, name="rCost")
    avgDist = model.addVar(vtype=GRB.CONTINUOUS, name="avgDist")
    rCostPt2_raw = model.addVar(vtype=GRB.CONTINUOUS, name="rCostPt2_raw")
    rCostPt2_sqrt = model.addVar(vtype=GRB.CONTINUOUS, name="rCostPt2_sqrt")
    for i in range(G.n):
        g[i] = model.addVar(vtype = GRB.BINARY, name="g_%g" % (i))

    for l in range(len(items)):
        y[items[l].k, l] = model.addVar(obj=-items[l].price, vtype=GRB.BINARY, name='y_%g_%g' % (items[l].k, l))
        a[l] = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%g' % l)
        for s in range(S):
            y[s + K, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g' % (s + K, l))
            z[s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, lb=items[l].lb, name='z_%g_%g' % (s + K, l))
            w[s + K, l] = model.addVar(obj=-1, vtype=GRB.CONTINUOUS, name='w_%g_%g' % (s + K, l))
            b[s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, name='b_%g_%g' % (s + K, l))

    for k in range(1, K):
        for s in range(S):
            p[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='p_%g_%g' % (k, s + K))

    model.update()

    if len(z_bd) > 0:
        for key in z.keys():
            z[key].lb = z_bd[key]
            z[key].ub = z_bd[key]

    #model.addConstrs(y[s + K, l] == 0 for s in range(S) for l in range(len(items)))

    model.addConstr(rCost == 2 * avgDist * sum(item.w for item in items)/G.q + 0.57 * rCostPt2_sqrt)
    model.addConstr(avgDist * (sum(g[k] for k in range(1, G.K)) + sum(g[i] for i in range(G.K, G.K + G.S))) == (sum(G.dist[k, 0] * g[k] for k in range(1, G.K)) + sum(G.dist[i, 0] * g[i] for i in range(G.K, G.K + G.S))))
    model.addConstr(rCostPt2_raw == G.A * (sum(g[k] for k in range(1, G.K)) + sum(g[i] for i in range(G.K, G.K + G.S))))
    model.addGenConstrPow(rCostPt2_raw, rCostPt2_sqrt, 0.5)

    #balance
    model.addConstrs(g[items[l].k] >= y[items[l].k, l] for l in range(len(items)))
    model.addConstrs(g[s + G.K] >= y[s + G.K, l] for s in range(G.S) for l in range(len(items)))

    #non-surcharging
    model.addConstrs(z[s + K, l] <= items[l].price for l in range(len(items)) for s in range(S))

    #duality
    model.addConstrs(quicksum(a[l] for l in Lk[k - 1]) == quicksum((items[l].ul - items[l].price) * y[items[l].k, l] + quicksum(items[l].ul*y[s + K, l] - w[s + K, l] for s in range(S)) for l in Lk[k - 1]) - quicksum(ui[k - 1][s]*p[k, s + K] for s in range(S)) + quicksum(items[l].inc * quicksum(y[s + K, l] for s in range(S)) for l in Lk[k - 1]) for k in range(1, K))
    model.addConstrs(y[items[l].k, l] + quicksum(y[s + K, l] for s in range(S)) == 1 for l in range(len(items)))     
    model.addConstrs(y[s + K, l] <= p[k, s + K] for k in range(1, K) for l in Lk[k - 1] for s in range(S))
    #model.addConstrs(p[k, s + K] <= quicksum(y[s + K, l] for l in Lk[k - 1]) for s in range(S) for k in range(1, K))
    model.addConstrs(a[l] >= items[l].ul - items[l].price for l in range(len(items)))
    model.addConstrs(a[l] + b[s + K, l]  >= items[l].ul - z[s + K, l] + items[l].inc for s in range(S) for l in range(len(items)))
    model.addConstrs(quicksum(b[s + K, l] for l in Lk[k - 1]) <= ui[k - 1][s] for s in range(S) for k in range(1, K))

    #linearization
    model.addConstrs(w[s + K, l] <= z[s + K, l] for l in range(len(items)) for s in range(S))
    model.addConstrs(w[s + K, l] <= items[l].price*y[s + K, l] for s in range(S) for l in range(len(items)))
    model.addConstrs(w[s + K, l] >= z[s + K, l] - 1.1*items[l].price*(1 - y[s + K, l]) for s in range(S) for l in range(len(items)))

    '''for k in range(1, K):
        for l in range(L):
            model.addConstr(y[k, 0, l] == 1)'''
    # model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\appr_det.lp")

    model.optimize()

    y_vals = {}

    for k in range(1, G.K):
        # print("Customer %d obj: %g" % (k, sum(a[l].x for l in Lk[k - 1])))
        for l in Lk[k - 1]:
            # print("%d, %d: %g %g %g"  % (k, items[l].type, y[items[l].k, l].x, items[l].ul, items[l].price))
            y_vals[k, l] = y[k, l].x
            for s in range(G.S):
                # print("%d, %d, %d: %g %g %g %g %g"  % (s + G.K, items[l].type, l, y[s + K, l].x, ui[k - 1][s], z[s + K, l].x, items[l].price, items[l].inc))
                y_vals[s + G.K, l] = y[s + G.K, l].x

    return y_vals

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
prob = getModel(G, S, L, R, ul, ui, h, q)'''