import gurobipy as grb
import queue
import math
import numpy as np
from dataclasses import dataclass
import bilevel_v2

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

    def __init__(self, graph, model, x, y, z, w, p, K, S, L, distThrshd, ui, ul):
        self.G = graph
        self.model = model
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.p = p
        self.n = graph.n
        self.K = K
        self.S = S
        self.L = L
        self.ui = ui
        self.ul = ul
        self.distThrshd = distThrshd
        self.model.Params.OutputFlag = 0
        self.model.optimize(bilevel_v2.callback_mult)
        root = self.Node(self.numNodes, -1, self.model.objVal, -1)
        root.parent = root
        self.nodeQueue.put(root)
        self.numNodes += 1

    def findConflictPair(self):
        for k in range(1, self.K):
                for i in range(self.K, self.K + self.S):#self.n):
                    if self.G.dist[i, k] > self.distThrshd:
                            flag = 0
                            for l in range(self.L):
                                if self.y[k, k, l].x > 0.8:
                                    y1 = self.y[k, k, l]
                                    flag += 1
                                    break
                            if flag > 0:
                                for l in range(self.L):
                                    if self.y[k, i, l].x > 0.8:
                                        y2 = self.y[k, i, l]
                                        flag += 1
                                    if flag > 1:
                                        return y1, y2
                    for j in range(self.K, i):
                        if self.G.dist[i, j] > self.distThrshd:
                            flag = 0
                            for l in range(self.L):
                                if self.y[k, i, l].x > 0.8:
                                    y1 = self.y[k, i, l]
                                    flag += 1
                                    break
                            if flag > 0:
                                for l in range(self.L):
                                    if self.y[k, j, l].x > 0.8:
                                        y2 = self.y[k, j, l]
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
            self.model.optimize(bilevel_v2.callback_mult)
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
            print("Obj: %g   part1: %g   part2: %g" % (sum(self.ul[k, l]*self.y[k, k, l].x - self.w[k, k, l].x 
            + sum(self.ul[k, l]*self.y[k, s + self.K, l].x - self.w[k, s + self.K, l].x for s in range(self.S)) for l in range(self.L)) 
            - sum(self.ui[k, s]*self.p[k, s + self.K].x for s in range(self.S)), sum(self.ul[k, l]*self.y[k, k, l].x - self.w[k, k, l].x 
            + sum(self.ul[k, l]*self.y[k, s + self.K, l].x - self.w[k, s + self.K, l].x for s in range(self.S)) for l in range(self.L)), 
            sum(self.ui[k, s]*self.p[k, s + self.K].x for s in range(self.S))))
            for l in range(self.L):
                print("%d, %d, %d: %g %g"  % (k, k, l, self.y[k, k, l].x, self.ul[k, l] - self.z[k, k, l].x))
                print("%d, %d, %d: %g -"  % (k, 0, l, self.y[k, 0, l].x))
                for s in range(self.S):
                    print("%d, %d, %d: %g %g"  % (k, s + self.K, l, self.y[k, s + self.K, l].x,  self.ul[k, l] - self.z[k, s + self.K, l].x))

        print("Routing:")
        for i, j, r in self.x.keys():
            if self.x[i, j, r].x > 0.5:
                print("%d, %d, %d" % (i, j, r))
