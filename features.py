from __future__ import division
import numpy,random
import networkx as nx

import warnings
warnings.simplefilter("ignore", numpy.ComplexWarning) # Ignore Casting warning.

""" Computes features of the graph and implements infection models. """

random.seed(10301949)


def compute_features(G,checks=False):
    """ Computes features of the graph. """

    n = G.order()
    
    nmeval = max(nx.adjacency_spectrum(G)) / n
    comps = nx.number_connected_components(G)
    mis = len(nx.maximal_independent_set(G)) / n
    density = nx.density(G)
    cc = nx.average_clustering(G)
    tris = sum(nx.triangles(G).values()) / n
    fracdeg1 = sum([len(G.neighbors(u)) == 1 for u in G]) / n
    fracdeg0 = sum([len(G.neighbors(u)) == 0 for u in G]) / n

    Gcc_list = list(nx.connected_component_subgraphs(G))
    Gcc = Gcc_list[0]
    ngcc = Gcc.order()
    mgcc = Gcc.size()
    
    return (nmeval,comps,mis,density,cc,tris,fracdeg1,fracdeg0,ngcc,mgcc)


#==============================================================================
#                              INFECTION MODELS
#==============================================================================
def compute_infected_set_sir(G,u,BETA):
    """ Computes infected set after infection of u in G using the SIR model. """
    assert 0 <= BETA <= 1
    
    if random.random() < BETA:
        Recovered = set()        
        Infected = set([u])       
        while len(Infected) > 0:
            infnode = Infected.pop()        
            Recovered.add(infnode)

            ToAdd = set()
            for infnodenb in G.neighbors(infnode):
                if infnodenb in Infected or infnodenb in Recovered: continue
                if random.random() < BETA: # success.
                    ToAdd.add(infnodenb)
                
            Infected.update(ToAdd)
    else:
        Recovered = set()

    return Recovered