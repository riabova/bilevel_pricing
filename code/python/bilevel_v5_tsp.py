from gurobipy import *
import Graph
from itertools import combinations
import matplotlib.pyplot as plt

def callback(model, where):
    if where == GRB.Callback.MIPSOL:
        adj = {}
        G = model._G
        for i in range(G.n):
            adj[i] = []
        x_sol = model.cbGetSolution(model._x)
        for i, j in x_sol.keys():
            if x_sol[i, j] > 0.5:
                adj[i].append(j)
                adj[j].append(i)
        cc = DFS(G, adj)
        if len(cc) > 1:
            for comp in cc:
                model.cbLazy(quicksum(model._x[i, j] + model._x[j, i] for i, j in combinations(comp, 2)) <= len(comp) - 1)


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

def getModel(G: Graph.Graph, items: list(), Lk: list(), ui, S: int, R: int, q: list(), c: float, z_bd={}):
    model = Model("bl5")
    model.setParam("OutputFlag", 0)
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    x = {} #leader's routing
    C = {} #routs oredring by cost
    y  ={} #follower's decision (continuous for duality)
    z = {} #leader's discount
    v = {} #leader's product-location/vehicle mapping
    w = {} #linerization
    a = {} #follower's dual 1
    b = {} #follower's dual 2
    #d = {} #follower's dual 3
    #gm = {} #follower's dual 4
    #th = {} #follower's dual 5
    p = {} #follower visits i
    g = {}
    
    #rev = model.addVar(vtype=GRB.CONTINUOUS, name="rev")
    rCost = model.addVar(vtype=GRB.CONTINUOUS, name="rCost")
    for i, j in G.dist.keys():
            x[i, j] = model.addVar(obj=c * G.dist[i, j], vtype=GRB.BINARY, name='x_%d_%d' % (i, j))
            x[j, i] = model.addVar(obj=c * G.dist[i, j], vtype=GRB.BINARY, name='x_%d_%d' % (j, i))
    for i in range(G.n):
            g[i] = model.addVar(vtype = GRB.BINARY, name="g_%d" % (i))

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

    model.addConstr(rCost == quicksum(c * G.dist[i, j] * (x[i, j] + x[j, i]) for i, j in G.dist.keys()))

    #balance
    model.addConstr(quicksum(x[0, j] for j in range(1, G.n)) == 1)
    model.addConstr(quicksum(x[j, 0] for j in range(1, G.n)) == 1)
    model.addConstrs(quicksum(x[i, j] for j in range(i + 1, G.n)) + quicksum(x[i, j] for j in range(0, i)) == g[i] for i in range(1, G.n))
    model.addConstrs(quicksum(x[j, i] for j in range(i + 1, G.n)) + quicksum(x[j, i] for j in range(0, i)) == g[i] for i in range(1, G.n))
    model.addConstrs(g[items[l].k] >= y[items[l].k, l] for l in range(len(items)))
    model.addConstrs(g[s + K] >= y[s + K, l] for s in range(S) for l in range(len(items)))

    #non-surcharging
    model.addConstrs(z[s + K, l] <= items[l].price for l in range(len(items)) for s in range(S))

    #duality
    model.addConstrs(quicksum(a[l] for l in Lk[k - 1]) == quicksum(items[l].ul*y[items[l].k, l] - items[l].price * y[items[l].k, l] + quicksum(items[l].ul*y[s + K, l] - w[s + K, l] for s in range(S)) for l in Lk[k - 1]) - quicksum(ui[k - 1][s]*p[k, s + K] for s in range(S)) + quicksum(items[l].inc * quicksum(y[s + K, l] for s in range(S)) for l in Lk[k - 1]) for k in range(1, K))
    model.addConstrs(y[items[l].k, l] + quicksum(y[s + K, l] for s in range(S)) == 1 for l in range(len(items)))     
    model.addConstrs(y[s + K, l] <= p[k, s + K] for k in range(1, K) for l in Lk[k - 1] for s in range(S))
    model.addConstrs(a[l]  >= items[l].ul - items[l].price for l in range(len(items)))
    model.addConstrs(a[l] + b[s + K, l]  >= items[l].ul - z[s + K, l] + items[l].inc for s in range(S) for l in range(len(items)))
    model.addConstrs(quicksum(b[s + K, l] for l in Lk[k - 1]) <= ui[k - 1][s] for s in range(S) for k in range(1, K))

    #linearization
    model.addConstrs(w[s + K, l] <= z[s + K, l] for l in range(len(items)) for s in range(S))
    model.addConstrs(w[s + K, l] <= items[l].price * y[s + K, l] for s in range(S) for l in range(len(items)))
    model.addConstrs(w[s + K, l] >= z[s + K, l] - 1.1*items[l].price*(1 - y[s + K, l]) for s in range(S) for l in range(len(items)))

    '''for k in range(1, K):
        for l in range(L):
            model.addConstr(y[k, 0, l] == 1)'''
    #model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\test_v4.lp")

    model._x = x
    model._g = g
    model._v = v
    model._y = y
    model._G = G
    model._R = R
    model.Params.lazyConstraints = 1

    model.optimize()

    y_vals = {}

    print("Obj: ", model.objVal)

    for k in range(1, G.K):
        print("Customer %d obj: %g" % (k, sum(a[l].x for l in Lk[k - 1])))
        for l in Lk[k - 1]:
            print("%d, %d: %g %g %g"  % (k, items[l].type, y[items[l].k, l].x, items[l].ul, items[l].price))
            y_vals[k, l] = y[k, l].x
            for s in range(G.S):
                print("%d, %d, %d: %g %g %g %g %g"  % (s + G.K, items[l].type, l, y[s + K, l].x, ui[k - 1][s], z[s + K, l].x, items[l].price, items[l].inc))
                y_vals[s + G.K, l] = y[s + G.K, l].x

    print("Routing: ")
    for i, j in x.keys():
        if x[i, j].x > 0.5:
            print("%d, %d" % (i, j))

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