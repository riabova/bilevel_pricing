import bnb
import Graph
import bilevel_v3
import bilevel_v4a

G = Graph.Graph()
G.read("..\..\data\TSP_instance_n_10_s_1.dat")
S = 3
R = 5
q = [91]*R
L = 2
h = [45, 60]#, 45, 46]
ui = {}
ul = {}
f = open("..\\..\\data\\ui_n_10_s_1_test1.txt", 'r')
f.readline()
k = 0
for line in f:
    k += 1
    utls = line.split()
    for i in range(S):
        ui[k, i] = int(utls[i])
f.close()
f = open("..\\..\\data\\ul_n_10_s_1_v4_test1.txt", 'r')
f.readline()
k = 0
for line in f:
    k += 1
    utls = line.split()
    for l in range(L):
        ul[k, l] = int(utls[l])
f.close()
modelInf = bilevel_v3.getModel(G, S, L, R, ul, ui, h, q)
dThrshd = 2
bnbTree = bnb.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], modelInf[6], modelInf[7], L, dThrshd, ui, ul)
#bnbTree.solve()
bnbTree.printSol()