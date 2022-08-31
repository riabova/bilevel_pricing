import bnb_v2
import Graph
import bilevel_v5
import bilevel_v4
import gen_utils
import time
import csv
import numpy as np
import os, sys
import matplotlib.pyplot as plt

class Item:
    def __init__(self, k, l, wt, util, bp, p):
        self.k = k
        self.type = l
        self.w = wt
        self.ul = util
        self.lb = bp
        self.price = p
        #self.ui = inc

def gen_utils(K: int, P: int, G: Graph.Graph, seed: int):
    np.random.seed(seed)
    d = max(G.dist.values())
    base_price = np.random.uniform(d/4, d, P)
    price = [np.random.uniform(base_price[i], 1.5 * base_price[i]) for i in range(P)]
    utils = np.transpose([np.random.beta(1, 5, K)*base_price[p] + 1.2*base_price[p] for p in range(P)])
    return base_price, price, utils

G = Graph.Graph()
S=3

#bnbTree.plotRoute()
def get_sol_info(G, p, s, l):
    S=s #change!
    q = [G.q]*G.r
    L = l
    h = np.random.uniform(0.01*G.q, 0.05*G.q, L)
    items = []
    Lk = []
    bprices, prices, utils = gen_utils(G.n - S - 1, L, G, seed=1)
    inconvs = [[p*np.mean(prices)**2/G.dist[max(k, s), min(k, s)] for s in range(G.n - S, G.n)] for k in range(1, G.n - S)]
    k = 0
    for u_p in utils:
        k += 1
        lk = []
        for l in range(L):
            if u_p[l] > 0:
                x = Item(k, l, h[l], u_p[l], bprices[l], prices[l])
                items.append(x)
                lk.append(len(items) - 1)
        Lk.append(lk)
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd)
    #bnbTree.solve()
    bnbTree.printSol()
    bnbTree.store_sol_info()
    return [bnbTree.profit, bnbTree.rCost, bnbTree.time, bnbTree.gap]

if __name__ == "__main__":
    #test over instances
    #input_path = sys.argv[1]
    s = 3
    l = 2
    script_dir = os.path.dirname(os.path.realpath(__file__))
    pct = 2
    G = Graph.Graph()
    G.read1("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\tests\E-n22-k4.vrp", S=S, seed=1)
    sol_info = get_sol_info(G, pct, s, l)
    '''rel_path = "\output\stats" + G.name + "_p%d.csv" % (100*pct)
    with open(script_dir + rel_path, "w", encoding="utf-16") as f1:
        f = csv.writer(f1, lineterminator="\n")
        f.writerow(["Instance", "Profit", "Routing Cost", "Runtime", "Gap"])
        sol_info = get_sol_info(G, pct, s, l)
        f.writerow([G.name] + sol_info)'''