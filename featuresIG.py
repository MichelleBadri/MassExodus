from __future__ import division
import random
import numpy as np
import networkx as nx
import igraph as ig
import numpy.linalg


import warnings
warnings.simplefilter("ignore", np.ComplexWarning) # Ignore Casting warning.

""" Computes features of the graph and implements infection models. """

# random.seed(10301949)

# Fill a ndarray vector with random integers between 0 and 1.
randArray = np.random.rand(1,100000)
# Convert 1D array to list for pop 
randArray = randArray[0].tolist()

def compute_features(G,checks=False):
    """ Computes features of the graph. """

    n = G.order()
    ### read in graph as IGraph network
    edge = list(G.edges())
    g = ig.Graph()
    for i in xrange(1,n+2):
        g.add_vertex(i)
    g.add_edges(edge)

    ## LARGEST EIGENVALUE
    # nmeval = max(nx.adjacency_spectrum(G)) / n  ### Original Calculation for largest Eigenval. Too slow 
    # IGRAPH REWORK
    e = g.evcent(directed=False, scale= False,return_eigenvalue=True)
    nmeval= round((max(e[0])),5) 
    
    ## NETWORKX CONNECTED COMPONENTS 
    # comps = nx.number_connected_components(G)  ### Original Calculation too slow 
    # IGRAPH REWORK
    comps = len(g.components()) - 1 

    ## MAXIMAL INDEPENDENT SET 
    mis = len(nx.maximal_independent_set(G)) / n ### Networkx version is fastest

    ## DENSITY 
    # density = nx.density(G) ### Original Calculation too slow
    # IGRAPH REWORK
    density = g.density() 

    ## CLUSTERING COEFFICIENT
    # cc = nx.average_clustering(G) ### Original Calculation is too slow
    # IGRAPH REWORK
    cc= g.transitivity_avglocal_undirected()

    # tris = sum(nx.triangles(G).values()) / n #### TOO SLOW
    # fracdeg1 = sum([len(G.neighbors(u)) == 1 for u in G]) / n
    # fracdeg0 = sum([len(G.neighbors(u)) == 0 for u in G]) / n

    # Gcc_list = list(nx.connected_component_subgraphs(G))
    # Gcc = Gcc_list[0]
    # ngcc = Gcc.order()
    # mgcc = Gcc.size()


    # return (nmeval,comps,mis,density,cc,tris,fracdeg1,fracdeg0,ngcc,mgcc)
    return (nmeval,comps,mis,density,cc)


#==============================================================================
#                              INFECTION MODELS
#==============================================================================
def compute_infected_set_sir(G,u,BETA):
    """ Computes infected set after infection of u in G using the SIR model. """
    assert 0 <= BETA <= 1
    
    if randArray.pop() < BETA: # Random number from numpy array
        Recovered = set()        
        Infected = set([u])       
        while len(Infected) > 0:
            infnode = Infected.pop()        
            Recovered.add(infnode)

            ToAdd = set()
            for infnodenb in G.neighbors(infnode):
                if infnodenb in Infected or infnodenb in Recovered: continue
                if randArray.pop() < BETA: # success.   # Random number from numpy array
                    ToAdd.add(infnodenb)
                
            Infected.update(ToAdd)
    else:
        Recovered = set()

    return Recovered

