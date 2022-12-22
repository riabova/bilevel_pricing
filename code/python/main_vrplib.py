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
    h = np.random.uniform(0.1*w, w, L) #product weights
    items = []
    Lk = []
    d = max(G.dist.values())
    np.random.seed(seed)
    price_lb = np.random.uniform(d/4, d, L)
    prices = [np.random.uniform(price_lb[i], 1.5 * price_lb[i]) for i in range(L)]
    inconvs = [[I_coef * G.dist[max(k, s), min(k, s)] for s in range(G.n - S, G.n)] for k in range(1, G.n - S)]
    rng = np.random.default_rng(seed=seed)
    for k in range(1, G.K):
        num_prod = np.random.randint(1, maxl + 1)
        prods = rng.choice(L, size=num_prod)
        lk = []
        wts = np.random.multinomial(G.demands[k - 1], [1/num_prod]*num_prod)
        for prod in range(num_prod):
            u_p = prices[prods[prod]] + 100*np.random.exponential(1) #check units
            items.append(Item(k, prods[prod], wts[prod], u_p, price_lb[prods[prod]], prices[prods[prod]]))
            lk.append(len(items) - 1)
        Lk.append(lk)
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd, I_coef)
    bnbTree.solve()
    #bnbTree.printSol()
    bnbTree.store_sol_info()
    #bnbTree.plotRouteMap()
    return [bnbTree.profit, bnbTree.rCost, bnbTree.time, bnbTree.gap, bnbTree.numNodes, bnbTree.discStats]

if __name__ == "__main__":
    #test over instances
    #input_path = sys.argv[1]
    s = 3
    l = 10 #all products
    maxl = 3 #max products in cart
    q = 100
    r = 4
    script_dir = os.path.dirname(os.path.realpath(__file__))
    I_coef = 1000
    G = Graph.Graph()
    #G.read1("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\CVRP_A\A-n32-k5.vrp", S=s, seed=1)
    '''f1 = "D:\Study\Ph.D\Projects\Bilevel Optimization\data\Buffalo\ss_dists.txt"
    f2 = "D:\Study\Ph.D\Projects\Bilevel Optimization\data\Buffalo\cc_dists.txt"
    f3 = "D:\Study\Ph.D\Projects\Bilevel Optimization\data\Buffalo\sc_dists.txt"
    f4 = "D:\Study\Ph.D\Projects\Bilevel Optimization\data\Buffalo\cust_coords.txt"
    f5 = "D:\Study\Ph.D\Projects\Bilevel Optimization\data\Buffalo\\ups_coords.txt"
    G.readSampleWithDists(f1, f2, f3, f4, f5, 11, 3, q, r)'''
    #sol_info = get_sol_info1a(G, I_coef, l, maxl)
    #print(sol_info)
    rel_path = "\output\statsSetA\\inf_inc.csv"# % (100*I_coef)
    rel_path1 = "\output\statsSetA\\inf_incDiscp.csv"# % (100*I_coef)
    with open(script_dir + rel_path, "w", encoding="utf-16") as file:
        f = csv.writer(file, lineterminator="\n")
        f.writerow(["Instance", "Profit", "Routing Cost", "Runtime", "Gap", "Nodes"])
        with open(script_dir + rel_path1, "w", encoding="utf-16") as disc_file:
            f1 = csv.writer(disc_file, lineterminator="\n")
            for inFile in os.listdir("D:\Study\Ph.D\Projects\Bilevel Optimization\\data\CVRP_A - Copy"):
                G.read1(os.path.join("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\CVRP_A - Copy", inFile), S=s, seed=1)
                sol_info = get_sol_info1a(G, I_coef, l, maxl)
                f.writerow([inFile] + sol_info[:-1])
                f1.writerow(sol_info[-1])