"""Contains methods to assign classical and diameter-modified Strahler orders to
vessels. Can group segments to elements.
Difficulties that arise:
    * selfloops attached to endpoints --> endpoints not identified as such.
    * selfloops consisting of multiple points are not easily identified as such, 
      if they contain vertices of order > 2.
These difficulties are solved on the basis of a flow computation. True 
endpoints are set to high pressure, the root is set to low pressure. Edges with
(almost) no-flow are assigned order zero.
"""


from __future__ import division, print_function, with_statement
import numpy as np
import pylab
import vgm




def upstream_edges(G, edge):
    """Returns the indices of those edges that connect the common vertex with
    vectices of higher pressure.
    INPUT: G: VascularGraph.
           edge: Index of the starting edge.
    OUTPUT: Edge-indices of the upstream edges.
    """
    upstreamEdges = []
    referencePressure = G.es[edge]['pressure']
    vTuple = G.es[edge].tuple
    vertex = sorted(zip(G.vs[vTuple]['pressure'], vTuple))[1][1]
    for e in G.adjacent(vertex):
        if e != edge:
            if G.es[e]['pressure'] > referencePressure:
                upstreamEdges.append(e)
    return upstreamEdges           



def downstream_edges(G, edge):
    """Returns the indices of those edges that connect the common vertex with
    vectices of lower pressure.
    INPUT: G: VascularGraph.
           edge: Index of the starting edge.
    OUTPUT: Edge-indices of the downstream edges.
    """
    downstreamEdges = []
    referencePressure = G.es[edge]['pressure']
    vTuple = G.es[edge].tuple
    vertex = sorted(zip(G.vs[vTuple]['pressure'], vTuple))[0][1]
    for e in G.adjacent(vertex):        
        if e != edge:
            if G.es[e]['pressure'] <= referencePressure:
                downstreamEdges.append(e)
    return downstreamEdges           


def classical_strahler_order(G):
    """Assigns order numbers to the edges of a tree-structured VascularGraph.
    The order is determined using the Strahler method, as described in the work
    of Strahler 'Quantitative analysis of watershed geomorphology' (Trans Am 
    Geophys Union, 1957).
    INPUT: G: VascularGraph of tree-structure.
    OUTPUT: None, G is modified in-place.
    """
    # Remove order property, if it already exists:
    if 'order' in G.es.attribute_names():
        print('WARNING: removing pre-existing order!')
    G.es['order'] = [None for e in G.es]

    # Find root and upstream endpoints on the basis of topology (the endpoint 
    # with the biggest diameter edge is the root vertex, all other endpoints are
    # upstream endpoints):
    endpoints = G.get_endpoints()
    if 'attachmentVertex' in G.attributes():
        rootVertex = G['attachmentVertex']
    else:
        maxDiameter = 0.0
        for ep in endpoints:
            if G.es(G.adjacent(ep))['diameter'][0] > maxDiameter:
                maxDiameter = G.es(G.adjacent(ep))['diameter'][0]
                rootVertex = ep
        print(rootVertex)
    if rootVertex in endpoints:            
        endpoints.remove(rootVertex)        

    # Solve for pressure to create an upstream / downstream ordering:
    G.vs[endpoints]['pBC'] = [2.0 for e in endpoints]
    G.vs[rootVertex]['pBC'] = 1.0
    if 'conductance' not in G.es.attribute_names():
        G.es['conductance'] = [1.0 for e in G.es]
    LS = vgm.LinearSystem(G)
    LS.solve_direct()    

    # Order edges by pressure, establishing a direction of flow: 
    G.es['pressure'] = [(G.vs[e.source]['pressure'] + 
                        G.vs[e.target]['pressure']) / 2.0 for e in G.es]
    peList = sorted(zip(G.es['pressure'], xrange(G.ecount())), reverse=True)
    
    # Move downstream, assigning orders as we go along:
    for pressure, edge in peList:
        # If no upstream edges exist, the edge is either an inflow edge or 
        # its flow value must be close to zero (e.g. in the case of an endpoint
        # or a loop). Assign order zero.
        upstreamEdges = upstream_edges(G, edge)
        if np.allclose(G.es[edge]['flow'], 0.0) or len(upstreamEdges) == 0:
            G.es[edge]['order'] = 0
        # If there is only one upstream edge, its order is retained.
        # If the upstream orders differ, the maximum order is chosen.
        # If the upstream orders are equal, the new order is the upstream
        # order + 1:
        else:
            orders = G.es[upstreamEdges]['order']
            if min(orders) == max(orders) and len(orders) > 1:
                G.es[edge]['order'] = orders[0] + 1
            else:
                G.es[edge]['order'] = max(orders)



def modified_strahler_order(G, fraction=1.0, limit=True, 
                            maxIterations=10, fstep=0.01):
    """Assigns order numbers to the edges of a tree-structured VascularGraph.
    The order is determined using the diameter-modified Strahler method, as
    described in the work of Kassab et al. 'Morphometry of pig coronary arterial
    trees' (Am J Physiol, 1993).
    This function assigns orders based primarily on the tree topology, e.g.
    ending vessels are assigned order zero, even if their diameter is well
    outside the mean+std of the zero-order bin. Also, an order decrease in
    downstream direction is not allowed, even if the diameter of the vessel
    would suggest this.
    INPUT: G: VascularGraph of tree-structure.
           fraction: The fraction of of edge-orders that need to remain constant
                     from one iteration to the next in order to stop.
           limit: Should the order increase going downstream be limited to one?
                  (Boolean.)
           maxIterations: The maximum number of iterations before the order
                          assignment is aborted and the fraction of required 
                          unchanged orders is lowered.
           fstep: The amount by which fraction is to be lowered each time the
                  iteration count reaches maxIterations.               
    OUTPUT: Returns a bash-like exit code, i.e. 0 for a successful run, and 1 
            if the provided fraction had to be reduced to succeed. 
            The edges of G are given a new property 'order'.
    """

    # Begin by assigning the classical Strahler order:
    classical_strahler_order(G)

    
    # Iteratively assign modified Strahler orders until convergence:
    diameters = np.array(G.es['diameter']) 
    edgePressures = [(G.vs[e.source]['pressure'] + 
                      G.vs[e.target]['pressure']) / 2.0 for e in G.es]
    peList = sorted(zip(edgePressures, xrange(G.ecount())), reverse=True)    

    ec = 0
    counter = 0
    while True:
        if counter == maxIterations:
            fraction -= fstep
            ec = 1
            print(fraction)
            counter = 0
        else:
            counter += 1
        oldOrder = G.es['order']
        maxOrder = max(oldOrder)
        mean = []
        std = []
        for order in range(maxOrder+1):
            mean.append(np.mean(G.es(order_eq=order)['diameter']))
            std.append(np.std(G.es(order_eq=order)['diameter']))
        bounds = {}
        for order in range(maxOrder):
            bounds[order] = (mean[order] + std[order] + 
                             mean[order+1] - std[order+1]) / 2.0
        bounds[maxOrder] = 1e100       

        
        # Move downstream, assigning orders as we go along:
        for pressure, edge in peList:
            # If no upstream edges exist, the edge is either an inflow edge or 
            # its flow value must be close to zero (e.g. in the case of an endpoint
            # or a loop). Assign order zero.
            upstreamEdges = upstream_edges(G, edge)
            if np.allclose(G.es[edge]['flow'], 0.0) or len(upstreamEdges) == 0:
                G.es[edge]['order'] = 0
            # If the diameter surpasses the bound of the current maximum
            # upstream order, the new order is increased accordingly.
            # Else, the new order is the maximum upstream order:
            else:
                order = max(G.es[upstreamEdges]['order'])
                while diameters[edge] > bounds[order]:
                    order = order + 1
                    if limit:
                        break
                G.es[edge]['order'] = order                   
        
        # End iteration if a given fraction of orders remains constant:
        if sum([1 for z in zip(G.es['order'], oldOrder) if z[0] == z[1]]) >= \
           fraction * len(oldOrder):
            return ec



def modified_strahler_order_v2(G, startingOrder=1, enforceMonotony=True):
    """Assigns order numbers to the edges of a tree-structured VascularGraph.
    The order is determined using the diameter-modified Strahler method, as
    described in the work of Kassab et al. 'Morphometry of pig coronary arterial
    trees' (Am J Physiol, 1993).
    This function reassigns orders based primarily on the mean and std of the 
    order bins. The enforce-monotony correction, which takes the topological 
    ordering into account is added as a final reassingment step. Without this 
    correction, the diameter-ranges of the different orders (> startingOrder) 
    are non-overlapping.
    INPUT: G: VascularGraph of tree-structure.
           startingOrder: The lowest order that should be included in the
                          sorting process.
           enforceMonotony: Whether or not to enforce that orders will never
                            decrease going downstream.
    OUTPUT: None, the edges of G are given a new property 'order'
    """

    # Begin by assigning the classical Strahler order:
    classical_strahler_order(G)

    # Iteratively assign modified Strahler orders until convergence:
    diameters = np.array(G.es['diameter'])
    preserve = G.es(order_lt=startingOrder).indices

    while True:
        oldOrder = G.es['order']
        maxOrder = max(oldOrder)

        mean = []
        std = []
        for order in range(maxOrder+1):
            mean.append(np.mean(G.es(order_eq=order)['diameter']))
            std.append(np.std(G.es(order_eq=order)['diameter']))

        bounds = {}
        for order in range(1, maxOrder):
            bounds[order] = ([(mean[order-1] + std[order-1] + 
                               mean[order] - std[order]) / 2.0, 
                              (mean[order] + std[order] + 
                               mean[order+1] - std[order+1]) / 2.0])
        bounds[0] = [0.0, bounds[1][0]]
        bounds[maxOrder] = [bounds[maxOrder-1][1], 1e100]
        
        for order in range(startingOrder, maxOrder+1):
            eIndices = np.nonzero((diameters >= bounds[order][0]) * 
                                  (diameters < bounds[order][1]))[0].tolist()
            eIndices = [e for e in eIndices if e not in preserve]
            G.es[eIndices]['order'] = [order for e in eIndices]

        if G.es['order'] == oldOrder:
            break
    
    # Ensure that orders do not decrease in downstream direction:
    if enforceMonotony:
        while True:
            oldOrder = G.es['order']
            for edge in G.es:
                edge['order'] = max(pylab.flatten([edge['order'],
                                    G.es[upstream_edges(G, edge.index)]['order']]))
            if G.es['order'] == oldOrder:
                break




class Permutation: 
    """Implements all possible permutations of a given array.
    """
    def __init__(self, alist): 
        """Initializes the Permutation generator.
        INPUT: alist: List or array to be permuted.
        OUTPUT: None
        """
        self._data = alist[:] 
        self._current = [] 
    def __iter__(self): 
        return self.next() 
    def next(self): 
        for elem in self._data: 
            if elem not  in self._current: 
                self._current.append(elem) 
                if len(self._current) == len(self._data): 
                    yield self._current[:] 
                else: 
                    for v in self.next(): 
                        yield v 
                self._current.pop() 




def pair_edges(G, vertex):
    """Pairs edges incident to a VascularGraph vertex depending on Strahler 
    order, flow direction and diameter. Attempts to pair each upstream edge
    with a downstream edge of the same order and closely matching diameter.
    INPUT: G: VascularGraph with assigned Strahler order and pressure defined
              at the vertices.
           vertex: The index of the vertex at which edges are to be paired.
    OUTPUT: List of edge pairs. Unpaired edges are not returned.       
    """
    # Create empty pair list, which will be returned as such, if no pairings 
    # are made:
    pairs = []

    # Determine upstream and downstream edges incident to vertex:
    us = []
    ds = []
    pressure = G.vs[vertex]['pressure']
    
    # Abort unphysiological cases:
    if len(G.neighbors(vertex)) > 5:
        return pairs
        
    for nb, aj in zip(G.neighbors(vertex), G.adjacent(vertex)):
        if G.vs[nb]['pressure'] > pressure:
            us.append(aj)
        else:
            ds.append(aj)

    # Loop over all orders present at the junction. Connections are only 
    # possible between edges of the same order:
    orders = G.es[us]['order']
    orders.extend(G.es[ds]['order'])
    orders = np.unique(orders).tolist()
    for order in orders:
        uso = [e for e in us if G.es[e]['order'] == order]
        dso = [e for e in ds if G.es[e]['order'] == order]

        # Continue if either up- or downstream edges of the current order do 
        # not exist:
        if min(len(uso), len(dso)) == 0:
            continue

        counter = -1
        while len(uso) < len(dso):
            uso.append(counter)
            counter = counter - 1
        while len(dso) < len(uso):
            dso.append(counter)
            counter = counter - 1
 
        # Make a list of all upstream edge permutations:
        usoperms = list(Permutation(uso))
 
        # Determine the minimal sum of diameter differences to choose the best 
        # pairing of upstream and downstream edges:
        mindiff = 1e100
        for i, usoperm in enumerate(usoperms):
            tmpdiff = 0
            for u, d in zip(usoperm, dso):
                if min([u, d]) < 0:
                    continue
                tmpdiff = tmpdiff + abs(G.es[u]['diameter'] - G.es[d]['diameter'])
            if tmpdiff < mindiff:
                mindiff = tmpdiff
                bestusoperm = usoperm[:]
               
        for z in zip(bestusoperm, dso):       
            if min(z) < 0:
                continue
            else:
                pairs.append(z)
               
    return pairs       



def assign_elements(G):
    """Traverses the vertices of a VascularGraph whose edges have been given
    Strahler orders. Edges incident to each vertex are paired according to 
    order and diameter values, as well as flow direction. Consecutive edges 
    (segments) of the same order with suitable (closely matching) diameters 
    form a common element.
    The result is a new VascularGraph that reflects this ordering scheme and 
    can thus be used for statistical analysis.
    INPUT: G: VascularGraph with Strahler-ordered edges.
    OUTPUT: Gel: VascularGraph of connected elements.
    """
    pairs = []
    for v in xrange(G.vcount()):
        pairs.extend(pair_edges(G, v))
    
    Ge = vgm.VascularGraph(G.ecount())
    Ge.add_edges(pairs)
    
    co = Ge.components()

    Gel = vgm.VascularGraph(len(co))
    Gel.vs['edges'] = [c for c in co]
    Gel.vs['order'] = [G.es[c[0]]['order'] for c in co]
    Gel.vs['vertices'] = [G.get_edge_vertices(e) for e in Gel.vs['edges']]

    edges = []
    connectingVertices = []
    for v1 in Gel.vs:
        for v2 in Gel.vs(xrange(v1.index + 1, Gel.vcount())):
            for vertex in v1['vertices']:
                if vertex in v2['vertices']:
                    edges.append((v1.index, v2.index))
                    connectingVertices.append(vertex)
                    break
    Gel.add_edges(edges)
    Gel.es['cVertex'] = connectingVertices

    return Gel



