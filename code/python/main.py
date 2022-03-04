import bnb_v2
import Graph
import bilevel_v5
import bilevel_v4
import gen_utils
import time
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from dataclasses import dataclass

@dataclass
class Item:
    def __init__(self, k, l, wt, util, bp):
        self.k = k
        self.type = l
        self.w = wt
        self.ul = util
        self.lb = bp
        #self.ui = inc

def gen_utils(K: int, P: int, seed: int, write=False, path="D:\\Study\\Ph.D\\Projects\\Bilevel Optimization\\data\\utils\\"):
    np.random.seed(seed)
    base_price = np.random.uniform(50, 200, P)
    #utils = np.transpose([np.random.normal(mu, 30, K) for mu in base_price])
    #plt.plot(np.random.beta(2, 5, K))
    #plt.show()
    utils = np.transpose([np.random.beta(1, 5, K)*base_price[p] + base_price[p] for p in range(P)])
    if write:
        with open(path + 'up_k%d_p%d_s%d.txt' % (K, P, seed), 'w') as f:
            f.write("%d %d \n" % (K, P))
            for line in utils:
                f.write(" ".join([str(x) for x in line]) + "\n")
    return base_price, utils

G = Graph.Graph()
S=3
G.read1("..\..\data\CVRP_A\A-n10-k1.txt", S=S, seed=1)
#G.read("..\..\data\TSP_instance_n_10_s_1.dat", S=S, seed=1)
q = [G.q]*G.r
L = 4
h = np.random.uniform(0.01*G.q, 0.05*G.q, L)
items = []
Lk = []
'''f = open("..\\..\\data\\ui_n_10_s_1_test1.txt", 'r')
f.readline()
for line in f:
    inc = line.split()
    inconvs.append([int(x) for x in inc])
f.close()
f = open("..\\..\\data\\ul_n_10_s_1_test1.txt", 'r')
f.readline()
k = 0
for line in f:
    k += 1
    utls = [int(x) for x in line.split()]
    lk = []
    for l in range(L):
        if utls[l] > 0:
            x = Item(k, l, h[l], utls[l])
            items.append(x)
            lk.append(len(items) - 1)
    Lk.append(lk)
f.close()'''
kmS = 100*np.mean(list(G.dist.values())) #bad!
inconvs = [[kmS/G.dist[max(k, s), min(k, s)] for s in range(G.n - S, G.n)] for k in range(1, G.n - S)]
prices, utils = gen_utils(G.n - S - 1, L, seed=1, write=True)
k = 0
for u_p in utils:
    k += 1
    lk = []
    for l in range(L):
        if u_p[l] > 0:
            x = Item(k, l, h[l], u_p[l], prices[l])
            items.append(x)
            lk.append(len(items) - 1)
    Lk.append(lk)
modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q)
dThrshd = 2
start = time.time()
bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd)
#bnbTree.solve()
print("Solved in %g" % (time.time() - start))
bnbTree.printSol()