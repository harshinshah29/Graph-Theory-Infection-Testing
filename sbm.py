import numpy as np
import networkx as nx

def SBM(N,M,q0,q1):
    '''This function is designed to generate the Stochastic Block Model.
    input params (consistent with the project description):
    sizes (list of int): the sizes of the communities
    q0 (float): probability of intra-subgroup connection
    q1 (float): probability of inter-subgroup connection

    output:
    G (N*N): adjacency matrix of the generated graph
    '''
    #################################################
    community_sizes = [N // M] * M
    for i in range(N % M):  
        community_sizes[i] += 1

    G = np.zeros((N, N), dtype=int)
    community_bounds = np.cumsum(community_sizes)
    start_indices = np.insert(community_bounds[:-1], 0, 0)

    # Intra-community connections
    for idx, size in zip(start_indices, community_sizes):
        end = idx + size
        block = np.random.rand(size, size) < q0
        np.fill_diagonal(block, 0)
        G[idx:end, idx:end] = block

    # Inter-community connections
    for i in range(N):
        start_idx = np.searchsorted(community_bounds, i, side='right')
        for j in range(i + 1, N):
            end_idx = np.searchsorted(community_bounds, j, side='right')
            if start_idx != end_idx:
                if np.random.rand() < q1:
                    G[i, j] = G[j, i] = 1


   
    #################################################

    return G
