import bnb_v2
import Graph
import bilevel_v4
from dataclasses import dataclass

@dataclass
class Item:
    def __init__(self, k, l, wt, util):
        self.k = k
        self.type = l
        self.w = wt
        self.ul = util
        #self.ui = inc

G = Graph.Graph()
G.read("..\..\data\TSP_instance_n_10_s_1.dat")
S = 3
R = 5
q = [91]*R
L = 2
h = [45, 60]#, 45, 46]
items = []
inconvs = []
Lk = []
f = open("..\\..\\data\\ui_n_10_s_1_test1.txt", 'r')
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
f.close()
modelInf = bilevel_v4.getModel(G, items, Lk, inconvs, S, R, q)
dThrshd = 2
bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd)
bnbTree.solve()
bnbTree.printSol()