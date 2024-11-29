from pyscipopt import Model
import numpy as np
import time


def read_adjacency_matrix(file_path):
    with open(file_path, 'r') as file:
        matrix = [list(map(int, line.split())) for line in file]
    return np.array(matrix)

def sparsest_cut_scipopt(graph_adj, demands):

    n = graph_adj.shape[0]
    model = Model("SparsestCut")
    model.hideOutput()
    x = {i: model.addVar(vtype="B", name=f"x_{i}") for i in range(n)}
    sparsity = model.addVar(vtype="C", name="sparsity", lb=0.0)

    
    cut_cost = model.addVar(vtype="C", name="cut_cost", lb=0.0)
    total_demand = model.addVar(vtype="C", name="total_demand", lb=1e-20)  

    model.addCons(
        cut_cost == sum(
            graph_adj[i, j] * (x[i] - x[i] * x[j])
            for i in range(n) for j in range(n) if graph_adj[i, j] >= 0
        ),
        name="CutCostConstraint",
    )
    
    model.addCons(
        total_demand == sum(
            demands[i, j] * (x[i] - x[i] * x[j])
            for i in range(n) for j in range(n)
        ),  
        name="TotalDemandConstraint",
    )
    model.addCons(
        sum(x[i] for i in range(n)) <= n-1,
        name="PartitionConstraint",
    )
    model.addCons(
        sum(x[i] for i in range(n)) >= 1,
        name="PartitionConstraint",
    )
    model.addCons(
        sparsity * total_demand >= cut_cost,
        name="SparsityConstraint",
    )
    model.setObjective(sparsity, "minimize")
    model.optimize()

    if model.getStatus() != "optimal":
        raise Exception("The problem does not have an optimal solution.")

    set_A = {i for i in range(n) if model.getVal(x[i]) < 0.5}
    set_B = {i for i in range(n) if model.getVal(x[i]) >= 0.5}
    best_sparsity = model.getVal(sparsity)

    return (set_A, set_B), best_sparsity


if __name__ == "__main__":
    file_path = "graph.txt"
    graph_adj = read_adjacency_matrix(file_path)
    demands = read_adjacency_matrix("demands.txt")
    start_time = time.time()
    cut, sparsity = sparsest_cut_scipopt(graph_adj, demands)
    end_time = time.time()
    # print("Partition:", cut)
    # print("Sparsity:", sparsity)
    print(sparsity)
    print((end_time-start_time))
