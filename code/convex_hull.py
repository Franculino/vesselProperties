import numpy as np


link = lambda a,b: np.concatenate((a,b[1:]))
edge = lambda a,b: np.concatenate(([a],[b]))


def dome(sample,base): 
    h, t = base
    dists = np.dot(sample-h, np.dot(((0,-1),(1,0)),(t-h)))
    outer = np.repeat(sample, dists>0, 0)

    if len(outer):
        pivot = sample[np.argmax(dists)]
        return link(dome(outer, edge(h, pivot)),
                    dome(outer, edge(pivot, t)))
    else:
        return base


def qhull(sample):
    if len(sample) > 2:
    	axis = sample[:,0]
    	base = np.take(sample, [np.argmin(axis), np.argmax(axis)], 0)
    	return link(dome(sample, base),
                    dome(sample, base[::-1]))
    else:
	    return sample


if __name__ == "__main__":
    #sample = 10*array([(x,y) for x in arange(10) for y in arange(10)])
    sample = 100*np.random.random((32,2))
    hull = qhull(sample)
	
    print "%!"
    print "100 500 translate 2 2 scale 0 0 moveto"
    print "/tick {moveto 0 2 rlineto 0 -4 rlineto 0 2 rlineto"
    print "              2 0 rlineto -4 0 rlineto 2 0 rlineto} def"
    for (x,y) in sample:
	    print x, y, "tick"
    print "stroke"
    print hull[0,0], hull[0,1], "moveto"
    for (x,y) in hull[1:]:
	    print x, y, "lineto"
    print "closepath stroke showpage"
