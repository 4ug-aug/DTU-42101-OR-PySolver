"""
Script til at lave Fundamental Insight følsomhedsanalyse på et lineært programmeringsproblem.
"""
import numpy as np
import sympy as sp
from collections import defaultdict 

def find_circle(edges):
    """
    Helper function to find circles in a graph.
    """

    # reformat nodes: (i, j) -> Ri, Rj
    edges = [(f"R{i}", f"L{j}") for i, j in edges]

    # Create a dictionary to store the graph
    adj_list = defaultdict(list)

    # Add edges to the graph
    for i, j in edges:
        adj_list[i].append(j)
        adj_list[j].append(i)

    def dfs(node, adj_list, visited, path, cycles):
        """
        Helper function to perform a DFS on a graph.
        """
        # Add the node to the path
        path.append(node)
        # Add the node to the visited set
        visited.add(node)

        # For all neighbors of the node
        for neighbor in adj_list[node]:
            # If the neighbor is not visited
            if neighbor not in visited:
                # Perform a DFS on the neighbor
                dfs(neighbor, adj_list, visited, path, cycles)
            # If the neighbor is in the path
            elif neighbor in path:
                # Add the path to the cycles
                cycles.append(path[path.index(neighbor):])
                # print(f"Found cycle: {path[path.index(neighbor):]}")
        # Remove the node from the path
        path.remove(node)

    # Now we want to find all cycles where we end up at the same node
    # We do this by doing a DFS on the graph, and if we end up at the same node
    # we have found a cycle
    cycles = []
    visited = set()
    for node in adj_list:
        if node not in visited:
            dfs(node, adj_list, visited, [], cycles)

    # return cycles longer than 2
    return [cycle for cycle in cycles if len(cycle) > 2]

def sensitivity_analysis(S_star, b, var_index, y_star=None):
    """
    Perform sensitivity analysis on a linear programming problem in tableau form,
    using the Fundamental Insight.

    Parameters:
    S_star (numpy array): The tableau of the optimized problem.
    S_star_b (numpy array): The solution of the optimized problem.

    Returns:
    A dictionary containing the range for delta.
    """

    var_index -= 1

    assert var_index < S_star.shape[1], "var_index must be less than the number of variables in the problem."
    # assert S_star.shape[1] == 3, "S_star must be a 3x3 matrix."

    d = sp.Symbol('d')

    mat = sp.Matrix([0, 0, 0])
    # put d in the correct index
    mat[var_index] = d
    m1 = S_star @ mat
    m2 = np.array(S_star @ b).reshape(-1, 1)

    # Create the left-hand side of the equation
    lhs = m1 + m2
    # Create the inequality expression
    # Create equation: LHS >= 0
    sols = []
    for i in range(lhs.shape[0]):
        sols.append(sp.solve(lhs[i] >= 0, d))
        
    for sol in sols:
        print(f"solution: {sol}")

    if y_star is not None:
        print(f"y_star: {y_star}")

def construct_matrix(A, b, c, basis):
    """
    Function to construct and solve a linear programming problem in tableau form.
    """

    # correct for 1-indexing
    basis -= 1
    
    # Add slack variables
    c_padded = np.hstack((c, np.zeros(A.shape[0])))

    AI = np.hstack((A, np.eye(A.shape[0])))
    B = AI[:, basis]

    cB = c_padded[basis]

    try:
        S_star = np.linalg.inv(B)
    except Exception as e:
        print("Matrix singular, using adding small difference to diag")
        small_diff = 1e-3 + np.ones(B.shape[0])
        S_star = np.linalg.inv(B+small_diff)

    print("b: ")
    print(b)
    print()

    print("S^*: ")
    print(S_star)
    print()

    # conduct sensitivity analysis by constructing equations from
    # S_star @ b + S_star*[d, 0, 0] >= 0



    y_star = cB @ S_star


    # Construct the tableau on form:
    # [
    # [y_star*AI - c | y_star | y_star*b]
    # [S_star*AI | S_star | S_star*b]
    # ]

    y_star_A = y_star @ A
    y_star_b = y_star @ b

    S_star_A = S_star @ A
    S_star_b = S_star @ b

    # Construct the tableau
    tableau = np.vstack((
        np.hstack((
            y_star_A - c,
            y_star,
            y_star_b
        )),
        np.hstack((
            S_star_A,
            S_star,
            S_star_b.reshape(-1, 1)
        ))
    ))

    # round all values to 2 decimals
    tableau = np.round(tableau, 2)

    print(tableau)

    return S_star, y_star

def total_unimodularity(A, m1):
    """
    Function to take a matrix and sets m1 and m2 of row indices to partition.
    Checks if the partitions are of total unimodularity.
    """

    # correct for 1-indexing
    m1 -= 1

    # Partition the rows of A
    m2 = np.setdiff1d(np.arange(A.shape[0]), m1)
    A1 = A[m1, :]
    A2 = A[m2, :]

    # For all columns in A1, check if the sum of values minus values in A2 is 0
    for i in range(A1.shape[1]):
        col_sum = np.sum(A1[:, i]) - np.sum(A2[:, i])
        if col_sum != 0:
            print(f"Column {i} is not 0")
            return False
        
    print(f"Partition: {m1}, {m2} is total unimodular")
    return True

def transport_problem(C, edges, flow = None):
    """
    Function to solve a Transport problem.
    The function takes the cost matrix and the edges of the graph.
    It computes the dual variables and the reduced costs.

    Args:
        C (numpy array): The cost matrix.
        edges (list): The edges of the graph.

    Returns:
        reduced_costs (numpy array): The reduced costs.
    """

    # Adjust edges for 1-indexing
    edges = [(i-1, j-1) for i, j in edges]

    # vj = ui + cij
    # ui = vj - cij

    u = np.nan*np.ones(C.shape[0])
    v = np.nan*np.ones(C.shape[1])

    u[0] = 0

    # Compute the dual variables by following the edges
    # We compute the dual variables by following the edges
    # We check for nan in either u or v, and if both are nan we cannot compute the dual variables yet
    while np.isnan(u).any() or np.isnan(v).any():
        for i, j in edges:
            if np.isnan(u[i]) and not np.isnan(v[j]):
                u[i] = v[j] - C[i, j]
                print(f"u{i+1} = v{j+1} - c{i+1}{j+1} = {v[j]} - {C[i, j]} = {u[i]}")
            elif np.isnan(v[j]) and not np.isnan(u[i]):
                v[j] = u[i] + C[i, j]
                print(f"v{j+1} = u{i+1} + c{i+1}{j+1} = {u[i]} + {C[i, j]} = {v[j]}")
    
    print("Dual variables computed")
    print(f"u: {u}")
    print(f"v: {v}")
    print()

    # Compute the reduced costs for each arc not in the list of edges
    edges_not_in = [(i, j) for i in range(C.shape[0]) for j in range(C.shape[1]) if (i, j) not in edges]

    reduced_costs = {}
    for i, j in edges_not_in:
        reduced_costs[(i, j)] = C[i, j] + u[i] - v[j]
        print(f"Reduced cost for arc {i+1}, {j+1}: {reduced_costs[(i, j)]}")

    # Find circles in the graph if we add the most negative reduced cost arc
    if min(reduced_costs.values()) < 0:
        print("There is a negative reduced cost arc, we need to find a circle")
        print()

        # Find the most negative reduced cost arc
        i, j = min(reduced_costs, key=reduced_costs.get)
        min_reduced_cost = reduced_costs[(i, j)]
        print(f"Most negative reduced cost arc: ({i+1}, {j+1})")
        print()

        new_edges = edges + [(i, j)]
        print(f"New edges: {[(i+1, j+1) for i, j in new_edges]}")
        print()

        # Find the circles in the graph
        circle = find_circle(new_edges)
        print(f"Circles: {circle}")
        print("First circle: ", circle[0]+[circle[0][0]])
        # Find reverse arcs in circle, where we go from L to R
        reverse_flow = []
        for i in range(1, len(circle[0])):
            if circle[0][i][0] == "L" and circle[0][i-1][0] == "R":
                reverse_flow.append((circle[0][i], circle[0][i-1]))

        print(f"Reverse flow: {reverse_flow}")
        print()

        # Reformat reverse_flow to old format
        reverse_flow = [(int(i[1]), int(j[1])) for i, j in reverse_flow]

        # If we got flow as input, we can compute the flow on the circle
        if flow is not None:
            # Create flow matrix by dict of edges and flow
            least_flow = []
            flow_dict = dict(zip(edges, flow))
            for i, j in reverse_flow:
                print(f"Flow on arc ({i}, {j}): {flow_dict[j, i]}")
                least_flow.append(flow_dict[j, i])
            print(f"Least flow: {min(least_flow)}")
            least_flow = min(least_flow)
            print()
            
            old_cost = np.sum([C[i, j]*flow_dict[i, j] for i, j in edges])
            print(f"Old cost: {old_cost}")
            new_cost = old_cost - abs(least_flow*min_reduced_cost)
            print(f"New cost: {new_cost}")


if __name__ == "__main__":

    # ------------ Fundamental Insight ------------

    # Coefficients
    A = np.array([
        [5, 4, 0],
        [-1, 4, 3]
    ])

    # Constraints
    b = np.array([3, 2])

    # Objective function
    c = np.array([4, 8, 9])

    # # Basis variables - rettet til at være 1-indexed
    basis = np.array([1, 5])

    S_star, y_star = construct_matrix(A, b, c, basis)

    # ------------ Sensitivity Analysis ------------

    # # Perform sensitivity analysis on the first variable
    # # Rettet til at være 1-indexed
    var_index = 2

    # sensitivity_analysis(S_star, b, var_index, y_star)

    # ------------ Total Unimodularity ------------

    A_mod = np.array([
        [1, 0, 0, 0, 0],
        [1, 1, 1, 0, 0],
        [0, 1, 0, 1, 0],
        [0, 0, 0, 1, -1],
        [0, 0, -1, 0, 1],
    ])

    m1 = np.array([1, 3])

    # total_unimodularity(A_mod, m1)

    # ------------ Transport Problem ------------

    C_tp = np.array([
        [8, 3, 9],
        [6, 4, 7],
        [7, 3, 9]
    ])

    # Edges are 1-indexed
    edges = [(1, 1), (1, 2), (2, 1), (2, 3), (3, 3)]

    flow = [4, 2, 1, 4, 6]

    # transport_problem(C_tp, edges, flow)

    test_A = np.array([
        [1, 0, 1, 0],
        [1, 1, 0, 0],
        [0, -1, -1, 0],
        [0, 1, 1, 1]

    ])

    # print(np.linalg.det(test_A))