import uncrt
import Graph
import bilevel_v5
# import bilevel_v5_approx
# import bilevel_v5_full
# import bilevel_v5_tsp2
import bnb_v2
import numpy as np
import os
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

def get_sol_info(G, L, maxl, expt="", nscen=10, niter=10, seed=7):
    S = G.S #change!
    q = [G.q]*G.r #vehicle capacities
    w = G.q * G.r/((G.K - 1) * maxl)
    np.random.seed(seed)
    h = np.random.uniform(0.5*w, 2.5 * w, L) #product weights
    items = []
    Lk = []
    c = 0.1
    d = 3400 #c * G.max_dist
    price_lb = np.random.uniform(d/10, 2 * d, L)
    prices = [np.random.uniform(price_lb[i], 1.5 * price_lb[i]) for i in range(L)]
    inconvs = []
    I_range = []
    rng = np.random.default_rng(seed=seed)
    prodInc = [2 * np.random.normal(- h[l]) * np.random.binomial(1, 0.3) + 0.01 * np.random.normal(prices[l]) * np.random.binomial(1, 0.2) for l in range(L)]
    prodInc = [prodInc[l] if -prodInc[l]/prices[l] > 0.03 or prodInc[l] > 0 else 0 for l in range(L)]
    dropout = [[] for k in range(G.K - 1)]
    for s in range(S):
        for k in range(G.K - 1):
            dropout[k].append(max(0.3, np.random.binomial(1, 0.98)))
    np.random.seed(seed)
    for k in range(1, G.K):
        i_rate = 5 * np.random.exponential(0.05)
        inconvs.append([i_rate * c * G.dist[max(k, s + G.K), min(k, s + G.K)] * dropout[k - 1][s] for s in range(S)])#
        tmp_rng = []
        for s in range(S):
            e_l = inconvs[-1][s] * np.random.uniform() #need better hypothesis
            e_r = 2 * inconvs[-1][s] * np.random.uniform()
            tmp_rng.append([max(0, inconvs[-1][s] - e_l), inconvs[-1][s] + e_r])
        I_range.append(tmp_rng)
        num_prod = 1#np.random.randint(1, maxl + 1)
        prods = rng.choice(L, size=num_prod, replace=False)
        lk = []
        for prod in prods:
            u_p = prices[prod] + 100*np.random.exponential(1) #check units
            items.append(Item(k, prod, h[prod], u_p, price_lb[prod], prices[prod], prodInc[prod]))
            lk.append(len(items) - 1)
        Lk.append(lk)

    # z_bd = [[584.5735346854182, 586.5918910522746, 583.3545037974426], [2109.534097171978, 2168.203631703444, 2023.5065103190948], [6472.389345240785, 6629.863107874747, 6253.324858774016], [8481.578115777067, 8481.578115777067, 8481.578115777067], [500.8490099752481, 500.8490099752481, 500.8490099752481], [5248.201474348687, 5749.445803739658, 5428.0438574712225], [3123.4050518093454, 3123.4050518093454, 3123.4050518093454], [2827.575598517667, 3123.4050518093454, 3123.4050518093454], [8481.578115777067, 8481.578115777067, 8481.578115777067], [5886.957308266643, 5811.935826275022, 6147.471100491255], [8481.578115777067, 8481.578115777067, 8481.578115777067], [4484.159308893743, 4484.159308893743, 4484.159308893743], [8096.977428214164, 8160.285613647437, 8097.223420419893], [8160.285613647435, 8160.285613647435, 7740.833726156585], [3202.1075248481184, 3038.701004199791, 3202.1075248481184], [765.948801012147, 976.0063990151905, 765.948801012147], [3123.4050518093454, 2916.0221019729124, 3123.4050518093454], [617.9857883334305, 617.9857883334305, 568.1295963991627]]
    # modelInf = bilevel_v5_approx.getModel(G, items, Lk, inconvs, S)
    # modelInf = bilevel_v5_full.getModel(G, items, Lk, inconvs, S, G.r, q, c)
    # modelInf2 = bilevel_v5_tsp2.getModel(G, items, Lk, inconvs, S, G.r, q, c)
    modelInf = bilevel_v5.getModel(G, items, Lk, inconvs, S, G.r, q, c)
    dThrshd = 2 #change!
    bnbTree = bnb_v2.BNB(G, modelInf[0], modelInf[1], modelInf[2], modelInf[3], modelInf[4], modelInf[5], items, Lk, inconvs, L, dThrshd, 0)
    bnbTree.solve()
    bnbTree.printSol()

    # learner = uncrt.UnsrtLearner(G, items, Lk, inconvs, S, q, c, I_range)
    # learner.refineUncrtS(niter=niter, nscen=nscen)

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
    K = int(sys.argv[9])
    S = int(sys.argv[10])
    nsc = int(sys.argv[11])
    nit = int(sys.argv[12])
    expt = sys.argv[13]'''
    s = 3
    l = 10 #all products
    maxl = 3 #max products in cart
    q = 100
    r = 1
    script_dir = os.path.dirname(os.path.realpath(__file__))
    tot_custs = 41#63
    tot_stores = 18
    G = Graph.Graph()
    #G.read1("..\..\\data\\tests\A-n10-k1.dat", S=S, seed=1)
    f1 = "..\..\\data\\Buffalo\\ss_dists_new.txt"
    f2 = "..\..\\data\\Buffalo\\cc_dists_rand1.txt"
    f3 = "..\..\\data\\Buffalo\\sc_dists_rand1.txt"
    f4 = "..\..\\data\\Buffalo\\cust_coords_rand1.txt"
    f5 = "..\..\\data\\Buffalo\\store_coords.txt"
    f6 = "..\..\\data\\Buffalo\\cc_routs_rand1.txt"
    f7 = "..\..\\data\\Buffalo\\sc_routs_rand1.txt"
    f8 = "..\..\\data\\Buffalo\\ss_routs_new.txt"
    expt = ""
    G.readSampleWithDistsRC(f1, f2, f3, f4, f5, 11, 3, q, r, rf1=f6, rf2=f7, rf3=f8, tot_custs=tot_custs, tot_stores=5)
    # G.readSampleWithDistsRC(in1, in2, in3, in4, in5, K + 1, S, q, r, rf1=f6, rf2=f7, rf3=f8, tot_custs=tot_custs, tot_stores=5)
    # G.calcArea()
    #G.readRandSampleWithDists(in1, in2, in3, in4, in5, K + 1, S, q, r, rf1=in6, rf2=in7, rf3=in8, tot_custs=tot_custs, tot_stores=tot_stores)
    sol_info = get_sol_info(G, l, maxl, expt="nscen", nscen=5, niter=200)
    #print(sol_info)

