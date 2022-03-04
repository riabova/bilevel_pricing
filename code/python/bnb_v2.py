import queue
import math
import time
from dataclasses import dataclass
import bilevel_v4
import bilevel_v5

class BNB:

    @dataclass
    class Node:
        childRight = 0
        childLeft = 0
        depth = 0
        status = 0
        #status:
        #0 = unprocessed
        #1 = violating (need to branch)
        #2 = candidate
        #3 = pruned by bound
        #4 = infeasible
    
        def __init__(self, name: int, parent, obj: float, var):
            self.name = name
            self.parent = parent
            self.obj = obj
            self.var = var
            if parent != -1:
                self.depth = self.parent.depth + 1
    
        def __lt__(self, other):
            
            return (self.obj > other.obj)
    
    numNodes = 0
    nodeQueue = queue.PriorityQueue()
    bestNode = -1
    tol = 0.000001
    depth = 0
    UB = 0

    def __init__(self, graph, model, x, y, z, w, p, items, Lk, ui, L, distThrshd):
        self.G = graph
        self.model = model
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.p = p
        self.n = graph.n
        self.items = items
        self.Lk = Lk
        self.ui = ui
        self.K = len(Lk) + 1
        self.S = len(ui[0])
        self.L = L
        self.distThrshd = distThrshd
        self.model.Params.OutputFlag = 0
        #start = time.time()
        self.model.optimize(bilevel_v5.callback_mult)
        #print("Solved in %g" % (time.time() - start))
        #self.model.Params.OutputFlag = 0
        root = self.Node(self.numNodes, -1, self.model.objVal, -1)
        root.parent = root
        self.nodeQueue.put(root)
        self.numNodes += 1

    def findConflictPair(self):
        for k in range(1, self.K):
                for i in range(self.K, self.K + self.S):#self.n):
                    if self.G.dist[i, k] > self.distThrshd:
                            flag = 0
                            for l in self.Lk[k - 1]:
                                if self.y[k, l].x > 0.8:
                                    y1 = self.y[k, l]
                                    flag += 1
                                    break
                            if flag > 0:
                                for l in self.Lk[k - 1]:
                                    if self.y[i, l].x > 0.8:
                                        y2 = self.y[i, l]
                                        flag += 1
                                    if flag > 1:
                                        return y1, y2
                    for j in range(self.K, i):
                        if self.G.dist[i, j] > self.distThrshd:
                            flag = 0
                            for l in self.Lk[k - 1]:
                                if self.y[i, l].x > 0.8:
                                    y1 = self.y[i, l]
                                    flag += 1
                                    break
                            if flag > 0:
                                for l in self.Lk[k - 1]:
                                    if self.y[j, l].x > 0.8:
                                        y2 = self.y[j, l]
                                        flag += 1
                                    if flag > 1:
                                        return y1, y2
        return -1

    def solve(self):
        while not self.nodeQueue.empty():
            node = self.nodeQueue.get()
            if self.bestNode == -1 or self.bestNode.obj > math.floor(node.obj):
                self.processNode(node)
            else:  #pruned by bound (UB, obj of parent)
                node.status = 3

    def processNode(self, node):
        current = node
        while (current.name != 0):  #fixing variable of a branch
            var = current.var
            var.UB = 0
            current = current.parent

        if(node.name != 0):
            self.model.optimize(bilevel_v5.callback_mult)
            node.obj = self.model.objVal
        if (self.model.status != 3):
            if (self.bestNode == -1 or math.floor(node.obj) > self.bestNode.obj):
                confl = self.findConflictPair()
                print(confl)
                if confl == -1:  #candidate (no conflict)
                    node.status = 2
                    self.bestNode = node
                else:  #branching
                    node.status = 1
                    node.leaf = 0
                    node_l = self.Node(self.numNodes, node, node.obj, confl[0])
                    self.nodeQueue.put(node_l)
                    node.childLeft = node_l
                    self.numNodes += 1
                    node_r = self.Node(self.numNodes, node, node.obj, confl[1]) 
                    self.nodeQueue.put(node_r)
                    node.childRight = node_r
                    self.numNodes += 1
            
                    if(self.depth < node_l.depth):
                        self.depth = node_l.depth
            else:  #pruned by bound
                node.status = 3

        else:  #infeasible
            node.status = 4   

        current = node
        while (current.name != 0):  #releasing variables of a branch
            var = current.var
            var.UB = 1
            current = current.parent

    def printSol(self):
        print("Leader's objective: %g" % -self.model.objVal)
        print("Customers decisions:")
        for k in range(1, self.K):
            print("Customer %d - Obj: %g   part1: %g   part2: %g" % (k, sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.w[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]) 
            - sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S)), sum(self.items[l].ul*self.y[self.items[l].k, l].x - self.w[self.items[l].k, l].x 
            + sum(self.items[l].ul*self.y[s + self.K, l].x - self.w[s + self.K, l].x for s in range(self.S)) for l in self.Lk[k - 1]), 
            sum(self.ui[k - 1][s]*self.p[k, s + self.K].x for s in range(self.S))))
            for l in self.Lk[k - 1]:
                print("%d, %d: %g %g"  % (k, self.items[l].type, self.y[self.items[l].k, l].x, self.items[l].ul - self.z[self.items[l].k, l].x))
                print("%d, %d: %g -"  % (0, self.items[l].type, self.y[0, l].x))
                for s in range(self.S):
                    print("%d, %d: %g %g"  % (s + self.K, self.items[l].type, self.y[s + self.K, l].x,  self.items[l].ul - self.z[s + self.K, l].x))

        print("Routing:")
        for i, j, r in self.x.keys():
            if self.x[i, j, r].x > 0.5:
                print("%d, %d, %d" % (i, j, r))
        print("Routing cost: %g" % self.model.getVarByName("C_0").x)
