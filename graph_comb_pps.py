import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def maximum_combinations(n, p):
    k = n / p
    r = 1
    for i in range(k):
        r *= sp.misc.comb(n-i*p, p)
    return r / sp.misc.factorial(k)

def fitting_curve(x, sigma):
    return sp.exp(-(x - 1.0)**2 / (2 * sigma**2))


used = [(3,3), (5,2), (2,4), (4,2)]

plt.clf()
for i, (parts, size) in enumerate(used):
    results = sp.load('results_%dx%d_n.npy' %(parts,size))
    results = results.T
    results[0] = results[0] / maximum_combinations(parts*size, size)
    fitted_ncomb, cov = sp.optimize.curve_fit(fitting_curve, results[2], results[0])
    print parts, size, cov
    plt.plot(results[2], results[0], '+', label='%dx%d' %(parts,size), color = cm.spectral(1.0*i/len(used)))
    plt.plot(sp.linspace(0,1,200), fitting_curve(sp.linspace(0,1,200), fitted_ncomb), '-', label='%dx%d' %(parts,size), color = cm.spectral(1.0*i/len(used)))

plt.xlim(0,1.2)
plt.ylim(0,1.2)
plt.xlabel('Algebraic connectivity (normalised Laplacian)')
plt.ylabel('Combinations normalised wrt theoretical maximum')    
plt.legend(loc=2)
plt.show()

