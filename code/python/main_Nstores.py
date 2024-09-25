import bnb_v2
import Graph
import bilevel_v5
import time
import csv
import numpy as np
import os, sys
import matplotlib.pyplot as plt

class Item:
    def __init__(self, k, l, wt, util, bp, p, inc):
        self.k = k
        self.type = l
        self.w = wt
        self.ul = util
        self.lb = bp
        self.price = p
        self.inc = inc

def gen_utils(K: int, P: int, G: Graph.Graph, seed: int):
    np.random.seed(seed)
    d = max(G.dist.values())
    base_price = np.random.uniform(d/4, d, P)
    price = [np.random.uniform(base_price[i], 1.5 * base_price[i]) for i in range(P)]
    utils = np.transpose([np.random.beta(1, 5, K)*base_price[p] + 1.2*base_price[p] for p in range(P)])
    return base_price, price, utils

def get_sol_info1a(G, I_coef, L, maxl, expt="", seed=7):
    S = G.S #change!
    q = [G.q]*G.r #vehicle capacities
    w = G.q * G.r/((G.K - 1) * maxl)
    np.random.seed(seed)
    h = np.random.uniform(0.1 * w, 1.5 * w, L) #product weights
    items = []
    Lk = []
    c = 0.1
    d = c * G.max_dist
    #np.random.seed(seed)
    price_lb = np.random.uniform(d/10, 2 * d, L)
    prices = [np.random.uniform(price_lb[i], 1.5 * price_lb[i]) for i in range(L)]
    inconvs = []
    rng = np.random.default_rng(seed=seed)
    prodInc = [2 * np.random.normal(- h[l]) * np.random.binomial(1, 0.3) + 0.01 * np.random.normal(prices[l]) * np.random.binomial(1, 0.2) for l in range(L)]
    prodInc = [prodInc[l] if -prodInc[l]/prices[l] > 0.03 or prodInc[l] > 0 else 0 for l in range(L)]
    irates = 3 * np.random.exponential(0.05, size=G.K-1)
    num_prods = np.random.randint(1, maxl + 1, size=G.K - 1)
    prods = [rng.choice(L, size=num_prods[k], replace=False) for k in range(G.K - 1)]
    devs = [100 * np.random.exponential(1, size=num_prods[k]) for k in range(G.K - 1)]
    u_ps = [[prices[prods[k][i]] + devs[k][i] for i in range(len(prods[k]))] for k in range(G.K - 1)]
    lk = []
    lk_count = 0
    Lk = []
    for k in range(G.K - 1):
        np.random.seed(seed)
        dropout = np.random.binomial(1, 1, S)
        np.random.seed(seed)
        inconvs.append([irates[k] * G.dist[max(k + 1 , s), min(k + 1 , s)] * dropout[s - G.n + S] for s in range(G.n - S, G.n)])#
        for i in range(len(prods[k])):
            items.append(Item(k, prods[k][i], h[prods[k][i]], u_ps[k][i], price_lb[prods[k][i]], prices[prods[k][i]], prodInc[prods[k][i]]))
            lk.append(lk_count)
            lk_count += 1
        Lk.append(lk)
        lk = []
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q, c)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd, I_coef)
    bnbTree.solve()
    bnbTree.printSol()
    bnbTree.store_sol_info()
    bnbTree.plotRouteMap(expt)
    spent = 0
    pc_stores_l = 0
    pc_stores_k = 0
    for l in range(len(items)):
        for s in range(G.S):
            if modelInf[4][s + G.K, l].x > 0.5:
                spent += (items[l].price - modelInf[4][s + G.K, l].x)
                print((items[l].price - modelInf[4][s + G.K, l].x))
                pc_stores_l += 1
    pc_stores_l = pc_stores_l/len(items)
    lbreak = False
    for k in range(1, G.K):
        for l in Lk[k - 1]:
            for s in range(G.S):
                if modelInf[4][s + G.K, l].x > 0.5:
                    pc_stores_k += 1
                    lbreak = True
                    break
            if lbreak:
                break
    pc_stores_k = pc_stores_k/(G.K - 1)
    pc_k_full = 0
    for k in range(1, G.K):
        sl = 0
        for l in Lk[k - 1]:
            for s in range(G.S):
                if modelInf[4][s + G.K, l].x > 0.5:
                    sl += 1
                    break
        if sl == len(Lk[k - 1]):
            pc_k_full += 1
    pc_k_full = pc_k_full/(G.K - 1)
    return [bnbTree.profit, bnbTree.rCost, bnbTree.time, bnbTree.gap, bnbTree.numNodes, spent, pc_stores_l, pc_stores_k, pc_k_full]

if __name__ == "__main__":
    #test over instances
    '''in1= sys.argv[1]
    in2= sys.argv[2]
    in3= sys.argv[3]
    in4= sys.argv[4]
    in5= sys.argv[5]
    in6= sys.argv[6]
    in7= sys.argv[7]
    in8= sys.argv[8]
    s = sys.argv[9]'''
    #s = 3
    l = 10 #all products
    maxl = 3 #max products in cart
    q = 100
    r = 4
    script_dir = os.path.dirname(os.path.realpath(__file__))
    tot_custs = 63
    tot_stores = 18
    G = Graph.Graph()
    #G.read1("D:\Study\Ph.D\Projects\Bilevel Optimization\data\\tests\A-n10-k1.dat", S=S, seed=1)
    f1 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\ss_dists.txt"
    f2 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\cc_dists.txt"
    f3 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\sc_dists.txt"
    f4 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\cust_coords.txt"
    f5 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\ups_coords.txt"
    f6 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\cc_routs.txt"
    f7 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\sc_routs.txt"
    f8 = "D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\ss_routs.txt"
    #G.readWithDists(in1, in2, in3, in4, in5, q, r, rf1=in6, rf2=in7, rf3=in8)
    #G.readWithDists(f1, f2, f3, f4, f5, q, r, rf1=f6, rf2=f7, rf3=f8)
    G.readRandSampleWithDists(f1, f2, f3, f4, f5, 11, 3, q, r, rf1=f6, rf2=f7, rf3=f8, tot_custs=tot_custs, tot_stores=tot_stores)
    #sol_info = get_sol_info1a(G, I_coef, l, maxl)
    #print(sol_info)
    I_coef = 11
    #G.readSampleWOstores(f2, f4, 63, 3, q, r, method="kmeans", rho=0.02, rf1="D:\Study\Ph.D\Projects\Bilevel Optimization\\data\\Buffalo\\cc_routs.txt")
    sol_info = get_sol_info1a(G, I_coef, l, maxl, expt="nstores")
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
