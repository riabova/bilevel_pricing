import math
from queue import PriorityQueue

class Graph:

    n = 0

    def read(self, file):
        f = open(file, 'r')
        self.n = int(f.readline())
        self.points = []

        for line in f:
            coords = line.split()
            self.points.append((float(coords[0]), float(coords[1])))

        self.dist = {(i, j):
        math.sqrt(sum((self.points[i][k] - self.points[j][k])**2 for k in range(2))) for i in range(self.n) for j in range(i)}

        f.close()

    def MST(self, start=0):
        q = PriorityQueue()
        q.put(([0, start]))
        cost = [0]
        color = ['white']
        parents = [None]
        for i in range(1, self.n-start):
            cost.append(float('inf'))
            color.append('white')
            parents.append(None)
        while not q.empty():
            v = q.get()[1]
            color[v-start] = 'black'
            for u in list(range(start, v)):
                if color[u-start] == 'white':
                    if cost[u-start] > self.dist[v, u]:
                        cost[u-start] = self.dist[v, u]
                        parents[u-start] = v
                        q.put(([cost[u-start], u]))
            for u in list(range(v+1, self.n)):
                if color[u-start] == 'white':
                    if cost[u-start] > self.dist[u, v]:
                        cost[u-start] = self.dist[u, v]
                        parents[u-start] = v
                        q.put(([cost[u-start], u]))
        return parents
