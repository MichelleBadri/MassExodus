#!/usr/bin/env python

from __future__ import division
import networkx as nx
import random,sys,time
import features as F
import argparse



""" Runs the SIR bursty node loss model using the duplication growth process. """ 

random.seed(10301949)


def doit(n,qmod,qcon,beta,T,outfile):
    G = nx.Graph()

    # Nodes that have gone through a DMC iteration and have not been burned.
    G.add_edge(1,2)
    InPlayPool = set([1,2])

    # Nodes that are initially isolated or have been burned.
    IsolatedPool = set()
    for i in xrange(3,n+1):
        G.add_node(i)
        IsolatedPool.add(i)

    # Header.
    print "#Iter\tNodes\tInPlay\tIso\tEdges\tNmEval\tComps\tMIS\tDensity\tCC\tTris\tFracDeg1\tFracDeg0\tNGcc\tMGcc"

    for iter in xrange(1,T+1):

        assert G.order() == len(InPlayPool) + len(IsolatedPool) == n # Always n nodes.

        # If this statement passes over, some nodes will get removed from burning
        # and then in the next iteration it will go through.
        if len(IsolatedPool) != 0:

            # 0. If all nodes are isolated, start over with the dumbell.
            if len(InPlayPool) == 0:
                u = IsolatedPool.pop()
                v = IsolatedPool.pop()
                G.add_edge(u,v)
                InPlayPool.add(u)
                InPlayPool.add(v)

            # 1. Select random node to add.
            v = IsolatedPool.pop()

            # 2. Select node to copy from.
            u = InPlayPool.pop()

            InPlayPool.add(v)
            InPlayPool.add(u) # Hack because sets can't return+retain random element.

            # 3. Run DMC iteration.
            Delete = set()
            for neighbor in G.neighbors(u):
                assert neighbor != u
                if random.random() < qmod: # modify the edge.
                    if random.random() < 0.5: # delete u->neighbor.
                        Delete.add(neighbor)
                        G.add_edge(v,neighbor)
                    #else: # delete v->neighbor -- already done.
                else: # don't modify the edge.
                    G.add_edge(v,neighbor)
                    assert v != neighbor

            for neighbor in Delete: G.remove_edge(u,neighbor)

            if random.random() < qcon:
                assert u != v
                assert not G.has_edge(u,v)
                G.add_edge(u,v)

        # 4. Burn and remove infected nodes.
        b = random.choice(G.nodes()) # random infected node.
        Infected = F.compute_infected_set_sir(G,b,beta)
        for b in Infected:
            if b in IsolatedPool: # Burnt nodes was already isolated.
                assert len(Infected) == 1 # no one else should be infected.
                continue
            else:
                InPlayPool.remove(b)
                IsolatedPool.add(b)

                # Remove and add-back as isolated.
                G.remove_node(b)
                G.add_node(b)

        # In case a burnt node / dmc iter causes an in-play node to become isolated.
        for u in nx.isolates(G):
            if u in InPlayPool:
                assert u not in IsolatedPool
                IsolatedPool.add(u)
                InPlayPool.remove(u)

        # 6. Compute distances after burning.
        (nmeval,comps,mis,density,cc,tris,fracdeg1,fracdeg0,ngcc,mgcc) = F.compute_features(G)
        m = G.size()
        
        # 8. Output results.
        print "%i\t%i\t%i\t%i\t%i\t%.5f\t%i\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%i\t%i" %(iter,G.order(),len(InPlayPool),len(IsolatedPool),m,nmeval,comps,mis,density,cc,tris,fracdeg1,fracdeg0,ngcc,mgcc)

    # 9. Last iteration: print component sizes for final graph.
    comps = ""
    for Gc in nx.connected_component_subgraphs(G):
        comps += "%i\t" %(Gc.order())
    print "# Components\t%s" %(comps.strip())

    # Print the final network as well.
    out = open(outfile,"w")
    #TODO: write header optionally
#    out.write("#nodes=%i\n#edges=%i\n#qmod=%.1f\n#qcon=%.1f\n#beta=%.2f\n" %(n,G.size(),qmod,qcon,beta))
#    out.write("#comps=%i\n#isolates=%i\n" %(nx.number_connected_components(G),len(nx.isolates(G))))
    for u,v in G.edges_iter(): out.write("%s\t%s\n" %(u,v))
    for u in nx.isolates(G): out.write("%s\t%s\n" %(u,u))       
    out.close()

def main():
    ## Parse options
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-n', '--nnodes', type=int, help='number of nodes', default=100)
    parser.add_argument('-s', '--qsteal', type=float, help='prob of stealing an edge, btw 0-1', default=0.3)
    parser.add_argument('-c', '--qcon', type=float, help='prob of connecting to a new node, btw 0-1', default=0.8)
    parser.add_argument('-b', '--beta', type=float, help='harshness parameter, btw 0-1', default=0.03)
    parser.add_argument('-t', '--niter', type=int, help='number of iterations', default=1000)
    parser.add_argument('-o', '--outfile', help='path to file to save graph')
    args = parser.parse_args()

    n    = args.nnodes
    qmod = args.qsteal # same as qsteal.
    qcon = args.qcon
    beta = args.beta
    T    = args.niter
    if args.outfile is None:
        outfile = "massexodus-%i-%.1f-%.1f-%.2f.graph" %(n,qmod,qcon,beta)
    else:
        outfile = args.outfile

    print outfile

    start = time.time()

#    n = int(sys.argv[1])
#    qmod = float(sys.argv[2]) # same as qsteal.
#    qcon = float(sys.argv[3])
#    beta = float(sys.argv[4])

    doit(n,qmod,qcon,beta,T,outfile)

    tot_time = (time.time()-start)/60
    print "# Time to run: %.3f (mins)" %(tot_time)

if __name__ == "__main__":
    main()
