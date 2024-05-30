import numpy as np
import random

def infect_step(G,p1,individuals,N):
    '''The function serves as the infection model for each day.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p1: the probability each individual infects neighbours.
    '''

    ###################################################
    new_infections = np.zeros(N)
    for i in range(N):
        if individuals[i] == 1:  # Check if individual i is infected
            for j in range(N):
                if G[i, j] == 1 and individuals[j] == 0:  # Check if j is a neighbor and not infected
                    if random.random() < p1:  # Infect with probability p1
                        new_infections[j] = 1
    individuals_updated = np.maximum(individuals, new_infections)  # Update the infection status
    
                        
    ###################################################
    return individuals_updated




def infect(G,p0,p1,time_steps):
    '''The function serves as the infection model for each day.
    input params (consistent with the project description):
    G (ndarray N*N): the adjacency matrix.
    p0: the infection probability for initial status.
    p1: the probability each individual infects neighbours.
    time_steps: log N
    '''
    N = G.shape[0]
    individuals = np.random.rand(N) < p0
    ###################################################
    for _ in range(time_steps):
        individuals = infect_step(G, p1, individuals, N)  # Simulate infection spread

    ###################################################

    return individuals