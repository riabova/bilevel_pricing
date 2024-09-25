from gurobipy import *
import Graph
from itertools import combinations
import matplotlib.pyplot as plt


def getModel(G: Graph.Graph, items: list(), Lk: list(), ui, S: int, R: int, q: list(), c: float):
    model = Model("bl5")
    model.setParam("NonConvex", 2)
    M = 100 * max([x.price for x in items])
    K = G.n - S #number of customers (assume stores are the last s locations, 0 is depot)
    x = {} #leader's routing
    C = {} #routs oredring by cost
    y  ={} #follower's decision (continuous for duality)
    z = {} #variable price (consts)
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

    for k in range(1, K):
        for s in range(S):
            d[k, s + K] = model.addVar(vtype=GRB.BINARY, name='d_%g_%g' % (k, s + K))

    model.update()


    #price lb
    # model.addConstrs(et[k, s + K] <= d[k, s + K] for k in range(1, K) for s in range(S))
    for k in range(1, K):
        for s in range(S):
                # model.addConstr(items[l].price - ui[k - 1][s] + quicksum(min(0, items[l].inc) * quicksum(y[s + K, l] for s in range(S)) for l in Lk[k - 1]) - items[l].lb <= M * d[k, s + K])
                # model.addConstr(items[l].lb - (items[l].price - ui[k - 1][s] + quicksum(min(0, items[l].inc) * quicksum(y[s + K, l] for s in range(S)) for l in Lk[k - 1])) <= M * (1 - d[k, s + K]))
            model.addConstr(float(sum(items[l].price for l in Lk[k - 1]) - ui[k - 1][s] - sum(items[l].lb for l in Lk[k - 1])) <= M * d[k, s + K])
            model.addConstr(-float(sum(items[l].price for l in Lk[k - 1]) - ui[k - 1][s] - sum(items[l].lb for l in Lk[k - 1])) <= M * (1 - d[k, s + K]))


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

    model.optimize()

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

    return model, x, y, z, w, p, et