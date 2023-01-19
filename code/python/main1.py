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

def get_sol_info(G, I_coef, l):
    S = G.S #change!
    q = [G.q]*G.r
    L = l
    h = np.random.uniform(0.01*G.q, 0.05*G.q, L)
    items = []
    Lk = []
    bprices, prices, utils = gen_utils(G.n - S - 1, L, G, seed=1)
    inconvs = [[I_coef * np.mean(prices)**2/G.dist[max(k, s), min(k, s)] for s in range(G.n - S, G.n)] for k in range(1, G.n - S)]
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
    #return items, inconvs, Lk
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd)
    #bnbTree.solve()
    bnbTree.printSol()
    bnbTree.store_sol_info()
    return [bnbTree.profit, bnbTree.rCost, bnbTree.time, bnbTree.gap]

def get_sol_info1a(G, I_coef, L, maxl, seed=7):
    S = G.S #change!
    q = [G.q]*G.r #vehicle capacities
    w = G.q * G.r/((G.K - 1) * maxl)
    h = np.random.uniform(0.1*w, 1.5 * w, L) #product weights
    items = []
    Lk = []
    d = max(G.dist.values())
    np.random.seed(seed)
    price_lb = np.random.uniform(d/4, d, L)
    prices = [np.random.uniform(price_lb[i], 1.5 * price_lb[i]) for i in range(L)]
    inconvs = []
    rng = np.random.default_rng()
    for k in range(1, G.K):
        i_rate = np.random.random()
        inconvs.append([i_rate * G.dist[max(k, s), min(k, s)] for s in range(G.n - S, G.n)])
        num_prod = np.random.randint(1, maxl + 1)#!!!!!!!!!!!!!!!!!!!REPLACE BY 1!!!!!!!!!!!!!!!!!!!!!!!!!
        prods = rng.choice(L, size=num_prod)
        lk = []
        for prod in prods:
            u_p = prices[prod] + 100*np.random.exponential(1) #check units
            items.append(Item(k, prod, h[prod], u_p, price_lb[prod], prices[prod]))
            lk.append(len(items) - 1)
        Lk.append(lk)
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd, I_coef)
    #bnbTree.solve()
    #bnbTree.printSol()
    bnbTree.store_sol_info()
    bnbTree.plotRouteMap()
    return [bnbTree.profit, bnbTree.rCost, bnbTree.time, bnbTree.gap, bnbTree.numNodes]

if __name__ == "__main__":
    #test over instances
    #input_path = sys.argv[1]
    s = 3
    l = 10 #all products
    maxl = 3 #max products in cart
    q = 100
    r = 4
    script_dir = os.path.dirname(os.path.realpath(__file__))
    I_coef = 0.1
    G = Graph.Graph()
    #G.read1("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\tests\A-n10-k1.dat", S=S, seed=1)
    f1 = "/home/gamma03/Projects/bilevel_pricing/data/Buffalo/ss_dists.txt"
    f2 = "/home/gamma03/Projects/bilevel_pricing/data/Buffalo/cc_dists.txt"
    f3 = "/home/gamma03/Projects/bilevel_pricing/data/Buffalo/sc_dists.txt"
    f4 = "/home/gamma03/Projects/bilevel_pricing/data/Buffalo/cust_coords.txt"
    f5 = "/home/gamma03/Projects/bilevel_pricing/data/Buffalo/ups_coords.txt"
    #G.readWithDists(f1, f2, f3, f4, f5, q, r)
    #G.readSampleWithDists(f1, f2, f3, f4, f5, 31, 4, q, r)
    #sol_info = get_sol_info1a(G, I_coef, l, maxl)
    #print(sol_info)
    I_coef = 0.3
    G.readSampleWithDists(f1, f2, f3, f4, f5, 11, 3, q, r)
    sol_info = get_sol_info1a(G, I_coef, l, maxl)
    print(sol_info)
    ''''insts = [(51, 10)]#, (63, 12)]#(11, 3), (21, 4), (31, 6), (41, 8), 
    rel_path = "\output\statsBuffalo" + "_p30_10_63.csv"# % (100*I_coef)
    with open(script_dir + rel_path, "w", encoding="utf-16") as file:
        f = csv.writer(file, lineterminator="\n")
        f.writerow(["Instance", "Profit", "Routing Cost", "Runtime", "Gap", "Nodes"])
        for I_coef in [0.1, 0.3]:
            G.readSampleWithDists(f1, f2, f3, f4, f5, 63, 10, q, r)
            #G.readWithDists(f1, f2, f3, f4, f5, q, r)
            sol_info = get_sol_info1a(G, I_coef, l, maxl)
            #sol_info = get_sol_info(G, I_coef, s, l)
            f.writerow(["k_%d" % 51] + sol_info)'''