from gurobipy import *
import Graph
from itertools import combinations
import matplotlib.pyplot as plt

def callback_mult(model, where):
    if where == GRB.Callback.MIPSOL:
        G = model._G
        g_sol = model.cbGetSolution(model._g)
        v_sol = model.cbGetSolution(model._v)
        y_sol = model.cbGetSolution(model._y)
        x_sol = model.cbGetSolution(model._x)
        for k in range(model._R):
            '''print("Routing:")
            for i in range(G.n):
                for j in range(i):
                    if x_sol[i, j, k] > 0.5:
                        print("%d, %d, %d: %d" % (i, j, k, x_sol[i, j, k]))'''
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
    cc = DFS(G, adj)
    if len(cc) > 1:
        for comp in cc:
            #print(comp)
            if len(comp) > 1:
                #gh = list(combinations(comp, 2))
                #print(quicksum(model._x[max(a), min(a), k] for a in combinations(comp, 2)) <= len(comp) - 1)
                model.cbLazy(quicksum(model._x[max(a), min(a), k] for a in combinations(comp, 2)) <= len(comp) - 1)


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

def getModel(G: Graph.Graph, items: list(), Lk: list(), ui, S: int, R: int, q: list(), c: float):
    model = Model("bl5")
    model.setParam("NonConvex", 2)
    M = 100 * max([x.price for x in items])
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    x = {} #leader's routing
    C = {} #routs oredring by cost
    y  ={} #follower's decision (continuous for duality)
    zp = {} #discount
    et = {} #leader's discount (binary)
    v = {} #leader's product-location/vehicle mapping
    w = {} #linerization
    a = {} #follower's dual 1
    b = {} #follower's dual 2
    # gm = {} #follower's dual 3
    # eps = {} #follower's dual 4
    p = {} #follower visits i
    g = {}
    d = {}
    
    #rev = model.addVar(vtype=GRB.CONTINUOUS, name="rev")
    rCost = model.addVar(vtype=GRB.CONTINUOUS, name="rCost")
    for i in range(1, G.n):
        for r in range(R):
            x[i, 0, r] = model.addVar(vtype=GRB.INTEGER, ub=2, name="x_%g_%g_%g" % (i, 0, r))
            for j in range(1, i):
                x[i, j, r] = model.addVar(vtype=GRB.BINARY, name="x_%g_%g_%g" % (i, j, r))
    for r in range(R):
        C[r] = model.addVar(vtype=GRB.CONTINUOUS, name="C_%g" % r)
        for i in range(G.n):
            g[r, i] = model.addVar(vtype = GRB.BINARY, name="g_%g_%g" % (r, i))

    for l in range(len(items)):
        #u_opt = max(items[l].ul - items[l].price, max([items[l].ul - items[l].price - ui[items[l].k - 1][s] + items[l].inc for s in range(S)]))
        y[items[l].k, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g' % (items[l].k, l))
        a[l] = model.addVar(lb=-GRB.INFINITY, vtype=GRB.CONTINUOUS, name='a_%g' % l)
        for r in range(R):
            v[items[l].k, l, r] = model.addVar(vtype=GRB.BINARY, name="v_%g_%g_%g" % (items[l].k, l, r))
        for s in range(S):
            #z[s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, name='z_%d_%d' % (s + K, l))#u_opt - (items[l].ul - items[l].price - ui[items[l].k - 1][s] + items[l].inc)#min(items[l].price + items[l].inc, items[l].price) #CONST
            y[s + K, l] = model.addVar(vtype=GRB.BINARY, name='y_%g_%g' % (s + K, l))
            b[s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, name='b_%g_%g' % (s + K, l))
            # gm[s + K, l] = model.addVar(vtype=GRB.CONTINUOUS, name='gm_%g_%g' % (s + K, l))
            for r in range(R):
                v[s + K, l, r] = model.addVar(vtype=GRB.BINARY, name="v_%g_%g_%g" % (s + K, l, r))
    for k in range(1, K):
        for s in range(S):
            zp[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='z_%d_%d' % (k, s + K))
            p[k, s + K] = model.addVar(vtype=GRB.BINARY, name='p_%g_%g' % (k, s + K))
            et[k, s + K] = model.addVar(vtype=GRB.BINARY, name='et_%g_%g' % (k, s + K))
            d[k, s + K] = model.addVar(vtype=GRB.BINARY, name='d_%g_%g' % (k, s + K))
            # eps[k, s + K] = model.addVar(vtype=GRB.CONTINUOUS, name='eps_%d_%d' % (k, s + K))

    model.update()

    model.addConstr(rCost == quicksum(C[r] for r in range(R)))

    #balance
    model.addConstrs(C[r] <= C[r - 1] for r in range(1, R))
    model.addConstrs(C[r] >= quicksum(c * G.dist[i, j]*x[i, j, r] for i in range(G.n) for j in range(i)) for r in range(R))
    model.addConstrs(quicksum(x[j, 0, r] for j in range(1, G.n)) == 2*g[r, 0] for r in range(R))
    model.addConstrs(quicksum(x[i, j, r] for j in range(i)) + quicksum(x[j, i, r] for j in range(i + 1, G.n)) == 2*g[r, i] for r in range(R) for i in range(G.n))
    model.addConstrs(g[r, 0] >= v[items[l].k, l, r] for l in range(len(items)) for r in range(R))
    model.addConstrs(g[r, 0] >= v[s + K, l, r] for s in range(S) for l in range(len(items)) for r in range(R))
    model.addConstrs(g[r, items[l].k] >= v[items[l].k, l, r] for l in range(len(items)) for r in range(R))
    model.addConstrs(g[r, k] <= quicksum(v[k, l, r] for l in Lk[k - 1]) for k in range(1, K) for r in range(R))
    model.addConstrs(g[r, s + K] >= v[s + K, l, r] for s in range(S) for l in range(len(items)) for r in range(R))

    #capacity
    model.addConstrs(quicksum(v[items[l].k, l, r] for r in range(R)) == y[items[l].k, l] for l in range(len(items)))
    model.addConstrs(quicksum(v[s + K, l, r] for r in range(R)) == y[s + K, l] for s in range(S) for l in range(len(items)))
    model.addConstrs(quicksum(items[l].w * (v[items[l].k, l, r] + quicksum(v[s + K, l, r] for s in range(S))) for l in range(len(items))) <= q[r] for r in range(R))

    #price lb
    # model.addConstrs(et[k, s + K] <= d[k, s + K] for k in range(1, K) for s in range(S))
    model.addConstrs(zp[k, s + K] >= ui[k - 1][s] - quicksum(items[l].inc * y[s + K, l] for l in Lk[k - 1]) for s in range(S) for k in range(1, K))
    model.addConstrs(zp[k, s + K] >= - ui[k - 1][s] + quicksum(items[l].inc * y[s + K, l] for l in Lk[k - 1]) for s in range(S) for k in range(1, K))
    # for k in range(1, K):
    #     for s in range(S):
            # model.addConstr((ui[k - 1][s] - quicksum(items[l].inc * y[s + K, l] for l in Lk[k - 1])) * et[k, s + K] >= 0)
            # model.addConstr(quicksum(items[l].price + items[l].inc * y[s + K, l] - items[l].lb for l in Lk[k - 1]) - ui[k - 1][s] <= M * d[k, s + K])
            # model.addConstr(quicksum(items[l].lb - (items[l].price + items[l].inc * y[s + K, l]) for l in Lk[k - 1]) + ui[k - 1][s] <= M * (1 - d[k, s + K]))
            # model.addConstr(float(sum(items[l].price for l in Lk[k - 1]) - ui[k - 1][s] - sum(items[l].lb for l in Lk[k - 1])) <= M * d[k, s + K])
            # model.addConstr(-float(sum(items[l].price for l in Lk[k - 1]) - ui[k - 1][s] - sum(items[l].lb for l in Lk[k - 1])) <= M * (1 - d[k, s + K]))

    #duality
    # model.addConstr(y[16, 0] == 1)
    # model.addConstr(y[16, 1] == 1)
    model.addConstrs(quicksum(a[l] for l in Lk[k - 1]) == 
                     quicksum((items[l].ul - items[l].price) * y[k, l] + quicksum((items[l].ul - items[l].price) * y[s + K, l] for s in range(S)) for l in Lk[k - 1]) 
                     - quicksum(ui[k - 1][s] * p[k, s + K] for s in range(S)) 
                     + quicksum(items[l].inc * quicksum(y[s + K, l] for s in range(S)) for l in Lk[k - 1])
                     + quicksum(zp[k, s + K] for s in range(S)) for k in range(1, K))
                     #+ quicksum(ui[k - 1][s] * et[k, s + K] * p[k, s + K] for s in range(S)) 
                     #- quicksum(items[l].inc * quicksum(y[s + K, l] * et[k, s + K] for s in range(S)) for l in Lk[k - 1]) for k in range(1, K)) 
    model.addConstrs(y[items[l].k, l] + quicksum(y[s + K, l] for s in range(S)) == 1 for l in range(len(items)))     
    model.addConstrs(y[s + K, l] <= p[k, s + K] for k in range(1, K) for l in Lk[k - 1] for s in range(S))
    # model.addConstrs(y[s + K, l] <= et[items[l].k, s + K] for l in range(len(items)) for s in range(S))
    model.addConstrs(a[l] >= items[l].ul - items[l].price for l in range(len(items)))
    model.addConstrs(a[l] + b[s + K, l] >= items[l].ul - items[l].price + items[l].inc for s in range(S) for l in range(len(items)))# - items[l].inc * et[items[l].k, s + K] for s in range(S) for l in range(len(items)))
    model.addConstrs(quicksum(b[s + K, l] for l in Lk[k - 1]) <= ui[k - 1][s] - zp[k, s + K] for s in range(S) for k in range(1, K))
    # model.addConstrs(eps[k, s + K] <= ui[k - 1][s] for s in range(S) for k in range(1, K))

    #obj
    # model.setObjective(- (quicksum(items[l].price for l in range(len(items))) 
    #                    - quicksum(ui[k - 1][s] * et[k, s + K] * p[k, s + K] for s in range(S) for k in range(1, K)) 
    #                    + quicksum(items[l].inc * et[items[l].k, s + K] * y[s + K, l] for l in range(len(items)) for s in range(S))
    #                    - quicksum(C[r] for r in range(R))))
    model.setObjective(- (quicksum(items[l].price for l in range(len(items))) 
                       - quicksum(zp[k, s + K] * p[k, s + K] for s in range(S) for k in range(1, K))
                       - quicksum(C[r] for r in range(R))))

    '''for k in range(1, K):
        for l in range(L):
            model.addConstr(y[k, 0, l] == 1)'''
    model.write("D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\code\\python\\models\\test_v4.lp")

    model._x = x
    model._g = g
    model._v = v
    model._y = y
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

    return model, x, y, zp, w, p, et