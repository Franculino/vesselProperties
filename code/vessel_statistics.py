from __future__ import division, print_function, with_statement
#import xlwt
import numpy as np
import scipy as sp
from scipy.spatial import kdtree
import os
import glob
import vgm
import string
import strahler
import convex_hull
#import matplotlib.pyplot as plt
#import pdb


def align_trees_with_z_axis(database, outdir):
    """Traverses the vessel database and aligns the vascular trees with the
    z-axis. If 'outdir' and 'database' are identical, the original trees will be
    overwritten.
    INPUT: database: The path to the vessel database.
           outdir: The output directory of the corrected vessels.
    OUTPUT: None (files written to disk).
    """
    # This step is already implemented in the database generation. Remove it
    # and implement it here, if needed.
    pass


def write_statistics(database, outfile, vType, minVcount=1):
    """Traverses the vessel database, assigns modified strahler orders to each
    vascular tree and groups segments of the same order into elements. Writes 
    the statistics based on each tree's elements to an excel file.
    INPUT: database: The path to the vessel database.
           outfile: The name of the output file containing the vessel 
                    statistics.
           vType: The type of trees. This can be either 'a' (arterial) or 'v'
                  (venous).
           minVcount: The minimum number of vertices in a tree to be considered
                      in the statistics. 
    OUTPUT: None (file written to disk).                
    """
    
    lengthDB = [[] for i in xrange(25)]
    vcountDB = []
    ## Prepare excel file for writing:
    #font0 = xlwt.Font()
    #font0.name = 'Times New Roman'
    #font0.colour_index = 2
    #font0.bold = True
    #style0 = xlwt.XFStyle()
    #style0.font = font0

    #font1 = xlwt.Font()
    #font1.name = 'Times New Roman'
    #font1.colour_index = 0
    #font1.bold = True
    #style1 = xlwt.XFStyle()
    #style1.font = font1

    #wb = xlwt.Workbook()

    # Loop over trees in database and write statistics to file:
    trees = [(int(os.path.basename(t)[:-4]), t) for t in 
             glob.glob(os.path.join(database, '*.pkl'))]
    trees = [t[1] for t in sorted(trees)]
    tenPercSteps = np.round(np.linspace(1, 10, 10) * len(trees) / 10.)
    percentages = np.linspace(10, 100, 10)
    print('Writing database...')
    treeCount = 0
    for i, tree in enumerate(trees):
        print(tree)
        # Read and process tree data:
        G = vgm.read_pkl(tree)
        vcountDB.append(G.vcount())
        if (i+1) in tenPercSteps:
            print('%i percent completed' % int(percentages[np.setmember1d(tenPercSteps, [i+1])][0]))
        if G.vcount() < minVcount:
            continue
#-------------------------------------------------------------
        # Delete order 2 vertices and add length and conductance 
        # (initially None for the new, joined edges). Rename the attachment
        # vertex (vertex indexing changes while order 2 vertices are removed):
        avr = G.vs[G['attachmentVertex']]['r']
        G.delete_order_two_vertices()    
        G.es['length'] = [sum(e['lengths']) for e in G.es]
        vgm.add_conductance(G, vType)
        # Delete selfloops:
        G.delete_selfloops()
        KDT = kdtree.KDTree(G.vs['r'])
        G['attachmentVertex'] = int(KDT.query(avr)[1])
#-------------------------------------------------------------

        strahler.modified_strahler_order(G, 0.99)
        Gel = strahler.assign_elements(G)
        maxOrder = max(Gel.vs['order'])
        cl = 0 # current line

        # Write tree ID, sample ID, vcount, ecount:
        ws = wb.add_sheet('Tree' + str(treeCount))
        treeCount = treeCount + 1
        ws.write(cl, 0, 'Tree ID', style0)
        ws.write(cl, 1, int(string.strip(os.path.basename(tree), '.pkl')))
        cl = cl + 1
        ws.write(cl, 0, 'Sample ID', style0)
        ws.write(cl, 1, G['sampleName'])
        cl = cl + 1        
        ws.write(cl, 0, 'Vertices', style0)
        ws.write(cl, 1, G.vcount())
        cl = cl + 1
        ws.write(cl, 0, 'Edges', style0)
        ws.write(cl, 1, G.ecount())
        cl = cl + 1
        ws.write(cl, 0, 'Dist to border', style0)
        ws.write(cl, 1, G['distanceToBorder'])
        cl = cl + 1
        ws.write(cl, 0, 'Root offset', style0)
        ws.write(cl, 1, G['avZOffset'])        
        cl = cl + 3
        
        # Write 'order' header:
        ws.write(cl, 0, 'Order', style0)
        for order in range(maxOrder + 1):
            ws.write(cl, order+1, order, style0)
        cl = cl + 1
        
        # Write element frequency 'N':
        ws.write(cl, 0, 'N', style1)
        for order in range(maxOrder + 1):
            ws.write(cl, order+1, len(Gel.vs(order_eq=order)))
        cl = cl + 1    
        
        # Write 'segments per element':
        # Note that this is the total number of elements of an order divided by
        # the total number of segments in that order (i.e. not on an element-by-
        # element basis).
        ws.write(cl, 0, 'SpE', style1)
        for order in range(maxOrder + 1):
            numElements = len(Gel.vs(order_eq=order))
            numSegments = sum([len(v['edges']) for v in Gel.vs(order_eq=order)])
            ws.write(cl, order+1, numSegments / numElements)
        cl = cl + 1

        # Write element 'length' and 'diameter':
        for eProperty in ['length', 'diameter']:
            ws.write(cl, 0, eProperty + ' (mean) [micron]', style1)
            ws.write(cl+1, 0, eProperty + ' (std) [micron]', style1)
            for order in range(maxOrder + 1):
                data = []
                for element in Gel.vs(order_eq=order):
                    if eProperty == 'length':
                        # Length of element is sum of segment-lengths:
                        data.append(sum(G.es[element['edges']][eProperty]))
                    else:
                        # Diameter of element is mean of segment-diameters:
                        data.append(np.mean(G.es[element['edges']][eProperty]))
                ws.write(cl, order+1, np.mean(data))
                ws.write(cl+1, order+1, np.std(data))
                if eProperty == 'length':
                    lengthDB[order].append(np.mean(data))
            cl = cl + 2
            
         
        # Add some additional whitespace:
        cl = cl + 2
         
        # Compute connectivity matrix and branching angles: 
        cm =  dict(zip(range(maxOrder + 1), 
                       [dict(zip(range(maxOrder + 1), 
                                 [[] for o in range(maxOrder + 1)]))
                         for o in range(maxOrder + 1)]))
        ba =  dict(zip(range(maxOrder + 1), 
                       [dict(zip(range(maxOrder + 1), 
                                 [[] for o in range(maxOrder + 1)]))
                         for o in range(maxOrder + 1)]))                         
        # Loop over the elements of G, which are the vertices of Gel:                                         
        for elm in Gel.vs:
            order = elm['order']
            orderCounter = dict(zip(range(maxOrder + 1), 
                                    [0 for o in range(maxOrder + 1)]))
            branchingAngle = dict(zip(range(maxOrder + 1), 
                                      [[] for o in range(maxOrder + 1)]))
                                      
            # Loop over neighboring elements:                        
            for nElm in Gel.vs(Gel.neighbors(elm.index)):
                # Find vertex common to both elements. 
                # Note that two elements can theoretically touch at more than 
                # one location, here only the first common vertex is considered.
                vertex = int(np.array(elm['vertices'])[np.setmember1d(elm['vertices'], nElm['vertices'])][0])
                # Find the edge that is part of the current element, adjacent,
                # and upstream to vertex. This ensures a 'natural', flow-based
                # choice of nElm-elm edge angle if the nElm edge connects to 
                # than one elm-edges.
                # Note that edges at equal mean pressures are not considered!
                edges = np.array(elm['edges'])[np.setmember1d(elm['edges'], G.adjacent(vertex))]
                pMax = G.vs[vertex]['pressure']
                edge = None
                for e in edges:
                    e = int(e) # Convert from numpy integer
                    meanP = np.mean(G.vs[G.es[e].tuple]['pressure'])
                    if meanP > pMax:
                        edge = e
                        pMax = meanP
                # Skip this element-junction if the current element does not 
                # contain upstream edges:        
                if edge is None:
                    continue       
                # Find the edge that is part of the neighboring element, 
                # adjacent, and downstream to vertex. This ensures a 'natural', 
                # flow-based choice of nElm-elm edge angle if the elm edge 
                # connects to more than one nElm-edges.
                # Note that edges at equal mean pressures are not considered!
                edges = np.array(nElm['edges'])[np.setmember1d(nElm['edges'], G.adjacent(vertex))]                                
                pMin = G.vs[vertex]['pressure']
                nEdge = None
                for e in edges:
                    e = int(e) # Convert from numpy integer
                    meanP = np.mean(G.vs[G.es[e].tuple]['pressure'])
                    if meanP < pMin:
                        nEdge = e
                        pMin = meanP
                # Skip this element-junction if the neighboring element does not 
                # contain downstream edges:        
                if nEdge is None:
                    continue                         

                # Increase the count of the connectivity matrix and compute the
                # branching angle:
                orderCounter[nElm['order']] += 1
                # A spherical region around the vertex is defined as the 
                # branching region. The radius of the sphere is equal to the 
                # maximum of the radii of the adjacent edges at the location of
                # the vertex.
                # The angle of an elements edge is then computed from the point
                # where the edge penetrates the sphere up to a distance of 
                # about twice the vessel diameter away from the point of 
                # penetration (unless of corse, the edge is shorter).
                # Radius of branching sphere:
                maxDiameter = 0
                for e in G.adjacent(vertex):
                    if G.es[e].source != vertex:
                        d = G.es[e]['diameters'][-1]
                    else:
                        d = G.es[e]['diameters'][0]
                    if d > maxDiameter:
                        maxDiameter = d
                radius = maxDiameter / 2.0
                # Vectors of the two branching edges:
                vectors = []
                for e in (edge, nEdge):
                    if G.es[e].source != vertex:
                        points = G.es[e]['points'][::-1]
                    else:
                        points = G.es[e]['points']
                    
                    startPoint = round(len(points) * 
                                       radius / G.es[e]['length'])
                    if startPoint > len(points)-1 or startPoint < 1:
                        startPoint = 1

                    # In order to determine the direction of a vessel leaving a
                    # bifurcation, we consider the points it comprises starting
                    # just outside the bifurcation (i.e. at maximum radius of  
                    # all vessels incident to the bifurcation) and ending two
                    # diameter lengths away from the bifurcation (angles 
                    # deviate least from the mean of angles determined from 1,
                    # 2, 3, 4 diameter lengths).
                    endPoint = min(len(points),
                                   round(len(points) * 
                                         (radius + 2.0 * G.es[e]['diameter']) /
                                         G.es[e]['length']))
                    endPoint = max(endPoint, startPoint+1)                     
                    points = points[startPoint-1:endPoint]
                    
                    vectors.append(vgm.average_path_direction(points))                           
                
                # The magnidudes of the vectors are one, hence the angle can
                # be computed without dividing by them. Note that an angle of 0
                # degrees means that both parent and child vessel are going in 
                # the same direction, whereas 180 degrees signifies opposite 
                # directions: 
                angle = np.pi - np.arccos(np.dot(vectors[0], vectors[1]))
                branchingAngle[nElm['order']].append(angle)          
     
            # Update the connectivity matrix and branching angle matrix with the
            # results of the current element:
            for nOrder in orderCounter.keys():
                cm[order][nOrder].append(orderCounter[nOrder])    
                ba[order][nOrder].extend(branchingAngle[nOrder])

        # Write connectivity matrix:        
        ws.write(cl, 0, 'Connectivity matrix', style1)
        cl = cl + 1
        ws.write(cl, 0, 'Order', style0)
        neo = [] # Number of elements of order o
        for j, order in enumerate(range(maxOrder + 1)):
            ws.write(cl, order+1, order, style0)
            ws.write(cl+j+1, 0, order, style0)
            neo.append(len(Gel.vs(order_eq=order)))
        cl = cl + 1                  
        for order in range(maxOrder + 1):
            for nOrder in range(order, maxOrder + 1):
            # Use the following line instead of the previous, if order zero
            # elements never connect to other order zero elements:
            # for nOrder in range(max(1, order), maxOrder + 1):
                if cm[order][nOrder] != []:
                    ws.write(order+cl, nOrder+1, sum(cm[order][nOrder]) / neo[nOrder])         
        # Equal spacing for all worksheets, irrespective of maximum order:
        cl = cl + 20
        #cl = cl + maxOrder + 3
        
        # Write branching angles:
        ws.write(cl, 0, 'Branching angles [deg]', style1)
        cl = cl + 1
        ws.write(cl, 0, 'Order', style0)
        neo = [] # Number of elements of order o
        for j, order in enumerate(range(maxOrder + 1)):
            ws.write(cl, order+1, order, style0)
            ws.write(cl+j+1, 0, order, style0)
            neo.append(len(Gel.vs(order_eq=order)))
        cl = cl + 1                  
        for order in range(maxOrder + 1):
            for nOrder in range(order, maxOrder + 1):
            # Use the following line instead of the previous, if order zero
            # elements never connect to other order zero elements:
            # for nOrder in range(max(1, order), maxOrder + 1):
                if ba[order][nOrder] != []:
                    ws.write(order+cl, nOrder+1, np.rad2deg(np.mean(ba[order][nOrder])))
        # Equal spacing for all worksheets, irrespective of maximum order:
        cl = cl + 20
        #cl = cl + maxOrder + 3
        
        # Compute irrigation / drainage volume:
        G.vs['z'] = [r[2] for r in G.vs['r']]
        zMin = np.min(G.vs['z'])
        zMax = np.max(G.vs['z'])
        cylinderLength = 200
        intervals = vgm.intervals(zMin, zMax, (zMax - zMin) / cylinderLength)
        pl = np.concatenate(G.es['points'], axis=0)
        z = pl[:,2]        
        
        if intervals == []:
            intervals = [[zMin, zMax]]
        volume = 0.0    
        for i, interval in enumerate(intervals):           
            points = pl[np.nonzero(map(all, zip(z>interval[0],
                                                z<interval[1])))[0], :] # use points
            points = np.array(zip(points[:,0], points[:,1]))
            #points = np.array(G.vs(z_ge=interval[0], z_le=interval[1])['r']) # use vertices
            
            if len(points) < 3:
                continue
            # Add slab volume: 
            com = np.mean(points, axis=0)
            hullpts = convex_hull.qhull(points)
            radius = np.mean([np.linalg.norm(com - hp) for hp in hullpts])
            volume += np.pi * radius**2 * cylinderLength           
                 
        ws.write(cl, 0, 'Irrigation / drainage volume [microL]', style1)
        cl = cl + 1
        ws.write(cl, 0, volume / 1e9) # Conversion micron**3 --> microL
        cl = cl + 3
        
    wb.save(outfile)
    
    #wb2 = xlwt.Workbook()
    #ws2 = wb2.add_sheet('Element Length')
    cl = 0
    maxOrder = 0
    maxN = []
    for order, ldb in enumerate(lengthDB):
        if len(ldb) > 0:
            maxOrder = order
        maxN.append(len(ldb))    
    for order in xrange(maxOrder + 1):
        ws2.write(0, order, order, style0)
        for i, l in enumerate(lengthDB[order]): 
            ws2.write(i+1, order, l, style1)
    wb2.save('elmLength_' + vType + '.xls')
