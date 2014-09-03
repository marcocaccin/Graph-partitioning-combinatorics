# little program to draw relationship between number of possible connected equipartitionings of a graph and some of his properties.
# Question: given a graph of size "total_nodes", let's try to divide it into "parts" parts of equal size "size", with the requirement that each part must be connected. How many possibilities are there?
# step 2) Can we correlate the possibilities with a number that is a property of a graph?


import networkx as nx
from numpy.random import random
import scipy as sp
import itertools
from sklearn.preprocessing import normalize
import sys, os
MACHINE_EPSILON = sp.finfo(sp.double).eps

def remove_repeated_vectors_from_matrix(a):
    # see http://stackoverflow.com/questions/16970982/find-unique-rows-in-numpy-array
    b = np.ascontiguousarray(a).view(np.dtype((np.void, a.dtype.itemsize * a.shape[1])))
    _, idx = np.unique(b, return_index=True)
    a = a[idx]
    return

def matching_elements(mylist):
    f = sp.array(mylist).flatten() 
    return sp.unique(f).size != f.size

def my_algebraic_connectivity(graph, normalise=False):
    if normalise:
        eigvals, eigvecs = sp.sparse.linalg.eigsh(nx.normalized_laplacian_matrix(graph).asfptype(), 2, which='SA')
        a = eigvals[1]
    else:
        eigvals, eigvecs = sp.sparse.linalg.eigsh(nx.laplacian_matrix(graph).asfptype(), 2, which='SA')
        a = eigvals[1]
    if a < MACHINE_EPSILON: a = 0.0
    return a


iterations = 300
normalise = True
if normalise:
    label = 'n'
else:
    label = ''
try:
    parts = int(sys.argv[1])
    size = int(sys.argv[2])
except ValueError:
    print "parts and size not understood."
    sys.exit()
total_nodes = parts * size

if os.path.exists("results_%dx%d_%s.npy" % (parts,size, label)):
    results = list(sp.load("results_%dx%d_%s.npy" % (parts,size, label)))
else:
    results = []

for iter in range(iterations):
    edge_creation_probability = 0.8 * random() + 0.1 # number between 0.1 and 0.9
    print "%0.3f" % edge_creation_probability
    # generate an initial random graph
    g = nx.fast_gnp_random_graph(total_nodes, edge_creation_probability)

    allpaths = []
    for start in g.nodes():
        for stop in g.nodes():
            # find all paths long at most "size"
            paths = nx.all_simple_paths(g, start, stop, cutoff = size)
            # purge all paths smaller than "size" and closed loops
            paths = [sorted(p) for p in paths if (len(p) == size and len(p) == len(sp.unique(p)))]
            try:
                # purge all paths that contain the same atoms
                remove_repeated_vectors_from_matrix(sp.array(paths))
            except:
                pass # print "only loops found, carrying on..."
            # add all the new connected subsets of atoms
            # allpaths += [p for p in paths if p not in allpaths]
            for p in paths:
                if p not in allpaths: allpaths.append(p)
    choices = []
    n_choices = 0
    for it in itertools.combinations(allpaths, parts):
        if not matching_elements(it):
            n_choices += 1
            # choices.append(it)
            # print it
    try:
        algebraic_connectivity = my_algebraic_connectivity(g, normalise=normalise)
    except:
        print "CAUGHT"
        algebraic_connectivity = nx.algebraic_connectivity(g, normalized=normalise)
    print n_choices, algebraic_connectivity, nx.average_clustering(g)
    results.append([n_choices, edge_creation_probability, algebraic_connectivity, nx.average_clustering(g)])
    sp.save("results_%dx%d_%s.npy" % (parts,size, label), sp.array(results))



