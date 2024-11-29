import numpy as np
from scipy.sparse.csgraph import shortest_path
from scipy.stats import median_abs_deviation
import math
from scipy.optimize import linprog
from scipy.sparse import lil_matrix
import cvxpy as cp
import networkx as nx
import time
Q = 10
NUM_TRIALS = 4


def read_adjacency_matrix(file_path):
    with open(file_path, 'r') as file:
        matrix = [list(map(int, line.split())) for line in file]
    return np.array(matrix)

def generate_frechet_embedding(graph_nx, demand_pairs):
    n = len(list(graph_nx.nodes()))
    distances = dict(nx.all_pairs_shortest_path_length(graph_nx))
    demand_set = []
    for (a, b, d) in demand_pairs:
        demand_set.append(a)
        demand_set.append(b)
    new_set = set(demand_set)
    ld = len(demand_set)
    ld_2n = 1 if ld == 0 else 2**(ld - 1).bit_length()
    num_attempts = 10*ld_2n
    
    while(len(new_set) < min(n,ld_2n) and num_attempts > 0):
        
        rand_num = np.random.randint(0, n)
        if rand_num not in demand_set:
            new_set.add(rand_num)
        num_attempts -= 1

    # if len(new_set) < ld_2n:
    #     raise Exception("Could not generate enough random points for Frechet embedding.")
    
    new_set = list(new_set)
    new_set.sort()

    tau = ld_2n.bit_length()-1
    L = Q * math.log(tau/2)

    frechet_sets = []
    for t in range(tau):
        for l in range(int(L) + 1):
            s = set()
            for i in range(2**(tau - t)):
                s.add(new_set[np.random.randint(0, len(new_set))])
            frechet_sets.append(s)

    
    
    embeddings = [[min([distances[i][j] for j in s]) for s in frechet_sets] for i in range(n)]
    return embeddings

def get_sets(graph_xe, embeddings):
    n = len(graph_xe)
    # print("len:  ", len(embeddings), len(embeddings[0]))
    embeddings_ind = [(embeddings[j], j) for j in range(len(embeddings))]
    F = []
    for i in range(len(embeddings[0])):
        emb_i = sorted(embeddings_ind, key=lambda x:x[0][i])

        subset = []
        for j in range(len(embeddings_ind)-1):
            subset.append(emb_i[j][1])
            F.append(set(subset.copy()))

    return F

def get_sparsity(graph_adj, demands, subset):
    c = 0
    d = 0
    if(len(subset)==len(graph_adj)):
        raise Exception("Complete set partition")
    
    for i in range(len(graph_adj)):
        for j in range(len(graph_adj)):
            c += graph_adj[i, j] if i in subset and j not in subset and graph_adj[i, j] >=0 else 0
            d += demands[i, j] if i in subset and j not in subset else 0

    if(d==0):
        return float('inf')
    
    return c/d
def scipy_solver(graph, demands):
    edges = list(graph.edges)
    nodes = list(graph.nodes)
    edge_to_index = {edge: i for i, edge in enumerate(edges)}
    node_to_index = {node: i for i, node in enumerate(nodes)}
    num_edges = len(edges)
    num_nodes = len(nodes)
    num_demands = len(demands)

    costs = np.array([graph[u][v]['weight'] for u, v in edges])
    num_vars = num_edges + num_demands + num_demands * num_nodes

    c = np.zeros(num_vars)
    c[:num_edges] = costs

    A_eq = lil_matrix((num_demands+1, num_vars))
    b_eq = np.zeros(num_demands+1)

    A_ub = lil_matrix((num_demands*(num_edges+1), num_vars))
    b_ub = np.zeros(num_demands*(num_edges+1))

    for i, (s, _, _) in enumerate(demands):
        d_idx = num_edges + num_demands + node_to_index[s]*num_demands + i
        # eq_row = np.zeros(num_vars)
        # eq_row[d_idx] = 1
        # A_eq.append(eq_row)
        # b_eq.append(0)
        A_eq[i, d_idx] = 1
        b_eq[i] = 0
    
    for i, (_, t, _) in enumerate(demands):
        d_idx = num_edges + num_demands + node_to_index[t]*num_demands + i
        y_idx = num_edges + i
        # up_row = np.zeros(num_vars)
        # up_row[d_idx] = -1
        # up_row[y_idx] = 1
        # A_ub.append(up_row)
        # b_ub.append(0)
        A_ub[i, d_idx] = -1
        A_ub[i, y_idx] = 1
        b_ub[i] = 0

    row_index = len(demands)
    for (u,v), e_idx in edge_to_index.items():
        for i in range(num_demands):
            d_u_idx = num_edges + num_demands + node_to_index[u]*num_demands + i
            d_v_idx = num_edges + num_demands + node_to_index[v]*num_demands + i
            # ub_row = np.zeros(num_vars)
            # ub_row[d_v_idx] = 1
            # ub_row[d_u_idx] = -1
            # ub_row[e_idx] = -1
            # A_ub.append(ub_row)
            # b_ub.append(0)
            A_ub[row_index, d_v_idx] = 1
            A_ub[row_index, d_u_idx] = -1
            A_ub[row_index, e_idx] = -1
            b_ub[row_index] = 0
            row_index += 1

    
    demand_weights = np.array([d for _, _, d in demands])
    # eq_row = np.zeros(num_vars)
    # eq_row[num_edges:num_edges+num_demands] = demand_weights
    # A_eq.append(eq_row)
    # b_eq.append(1)
    A_eq[-1, num_edges:num_edges+num_demands] = demand_weights
    b_eq[-1] = 1

    # A_eq = np.array(A_eq)
    # b_eq = np.array(b_eq)
    # A_ub = np.array(A_ub)
    # b_ub = np.array(b_ub)

    # print("Solving LP")
    result = linprog(c, A_eq=A_eq.tocsr(), b_eq=b_eq, A_ub=A_ub.tocsr(), b_ub=b_ub, method='highs')

    if not result.success:
        raise Exception("LP solver failed")
    
    x_values = result.x[:num_edges]
    y_values = result.x[num_edges:num_edges+num_demands]
    return {
        "x": {e: x_values[edge_to_index[e]] for e in edges},
        "y": {i: y_values[i] for i in range(num_demands)}
    }

def sparsest_cut_polynomial(nx_graph, demands_tuples):
    edges = list(nx_graph.edges())
    nodes = list(nx_graph.nodes())
    edges_to_index = {edges[i]: i for i in range(len(edges))}
    nodes_to_index = {nodes[i]: i for i in range(len(nodes))}
    n_edges = len(edges)
    n_nodes = len(nodes)
    n_demands = len(demands_tuples)
    x = cp.Variable(n_edges, nonneg = True)
    y = cp.Variable(n_demands, nonneg = True)
    dist = cp.Variable((n_nodes, n_demands), nonneg = True)
    costs = np.array([nx_graph[u][v].get('weight', 1) for u, v in edges])
    objective = cp.Minimize(cp.sum(cp.multiply(costs, x)))
    constraints = []

    for i, (s, t, d_i) in enumerate(demands_tuples):
        constraints.append(dist[nodes_to_index[s], i] == 0)
        constraints.append(dist[nodes_to_index[t], i] >= y[i])

    for (u,v), edge_idx in edges_to_index.items():
        for i in range(n_demands):
            u_idx = nodes_to_index[u]
            v_idx = nodes_to_index[v]
            constraints.append(dist[u_idx, i] <= dist[v_idx, i] + x[edge_idx])

    demand_weights = np.array([d_i for _, _, d_i in demands_tuples])
    constraints.append(cp.sum(cp.multiply(demand_weights, y)) == 1)

    prob = cp.Problem(objective, constraints)
    # print("Solving LP")
    prob.solve()

    x_values = {e: x.value[edges_to_index[e]] for e in edges}
    y_values = {f"y({i})": y.value[i] for i in range(n_demands)}

    return {"x": x_values, "y": y_values}

def solve_lp(graph_adj, demands):
    n = graph_adj.shape[0]
    c = np.array([graph_adj[i, j] for i in range(n) for j in range(n) if graph_adj[i, j] >= 0])
    d = np.array([demands[i, j] for i in range(n) for j in range(n) if demands[i, j] >= 0])
    c_map = [(i, j) for i in range(n) for j in range(n) if graph_adj[i, j] >= 0]
    d_map = [(i, j) for i in range(n) for j in range(n) if demands[i, j] >= 0]

    

    X = np.zeros(len(c) + len(d))

def dictionary_edges_to_adj(edges, n):
    adj = np.zeros((n, n))
    for (i, j) in edges:
        adj[i, j] = 1
        adj[j, i] = 1
    return adj

def adj_to_nx(graph_adj):
    n = graph_adj.shape[0]
    nx_graph = nx.Graph()
    for i in range(n):
        for j in range(i+1,n):
            if graph_adj[i, j] >= 0:
                nx_graph.add_edge(i, j, weight=graph_adj[i, j])
    return nx_graph

def adj_to_tuple(demands):
    n = demands.shape[0]
    demands_tuples = []
    for i in range(n):
        for j in range(i+1, n):
            if demands[i, j] > 0:
                demands_tuples.append((i, j, demands[i, j]))
    return demands_tuples
            
def xes_to_nx(xes, n):
    nx_graph = nx.Graph()
    for e in xes:
        i, j = e
        nx_graph.add_edge(i, j, weight=xes[e])
    return nx_graph

def fretched_sparsest_cut(graph_adj, demands, num_trials=10):
    n = graph_adj.shape[0]
    nx_graph = adj_to_nx(graph_adj)
    demands_tuples = adj_to_tuple(demands)

    best_cut = None
    best_ratio = float('inf')

    for _ in range(num_trials):
        # print("Starting trial")
        result = scipy_solver(nx_graph, demands_tuples)
        # print("Solved LP")
        x = result["x"]
        y = result["y"]
        xes = {e: x[e] for e in x}
        nx_sol = xes_to_nx(xes, n)
        fretchet = generate_frechet_embedding(nx_sol, demands_tuples)
        # print("Generated Frechet embedding")
        F = get_sets(graph_adj, fretchet)
        # print("Generated sets")
        for f in F:
            ratio = get_sparsity(graph_adj, demands, f)
            if ratio < best_ratio:
                best_ratio = ratio
                best_cut = f

    # print(fretchet)
    # print(F)
    return best_cut, best_ratio


if __name__ == "__main__":
    graph_file = "graph.txt"
    graph_adj = read_adjacency_matrix(graph_file)

    demands = read_adjacency_matrix("demands.txt")

    # print(get_sparsity(graph_adj, demands, {1, 2, 3, 4, 5, 7, 8}))

    start_time = time.time()
    cut, sparsity_ratio = fretched_sparsest_cut(graph_adj, demands, num_trials=NUM_TRIALS)
    end_time = time.time()
    print(sparsity_ratio)
    print(end_time-start_time)
    # print("Partition:", cut, set(range(len(graph_adj)))-cut)
    # print("Sparsity Ratio:", sparsity_ratio)




    

    

