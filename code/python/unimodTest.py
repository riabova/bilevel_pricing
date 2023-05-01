from gurobipy import *

def test(S, ui, items, Lk, disc):
    model = Model("test")
    
    y  ={} #follower's decision (continuous for duality)
    p = {} #follower visits i

    #I = model.addVar(name="I")
    for l in Lk:
        y[S, l] = model.addVar(ub=1, vtype=GRB.CONTINUOUS, name='y_%g_%g' % (0, l))
        for s in range(S):
            y[s, l] = model.addVar(ub=1, vtype=GRB.CONTINUOUS, name='y_%g_%g' % (s, l))
    for s in range(S):
        p[s] = model.addVar(ub=1, vtype=GRB.CONTINUOUS, name='p_%g' % (s))
    model.update()

    model.addConstrs(y[S, l] + quicksum(y[s, l] for s in range(S)) == 1 for l in Lk)     
    model.addConstrs(y[s, l] <= p[s] for l in Lk for s in range(S))
    model.addConstrs(p[s] <= quicksum(y[s, l] for l in Lk) for s in range(S))

    model.setObjective(quicksum(items[l].ul*y[S, l] - items[l].price * y[S, l] + quicksum(items[l].ul*y[s, l]
                         - (items[l].price - disc[s]) * y[s, l] for s in range(S)) for l in Lk) - quicksum(ui[s]*p[s] for s in range(S))
                         + quicksum(items[l].inc * quicksum(p[s] for s in range(S)) for l in Lk), GRB.MAXIMIZE)
    
    model.update()
    model.optimize()

    for l in Lk:
        print(y[S, l].x)
        for s in range(S):
            print(y[s, l])
    for s in range(S):
        print(p[s])