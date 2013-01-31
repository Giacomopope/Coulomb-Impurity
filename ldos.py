import numpy as np
import math
import matplotlib
from matplotlib import pylab
from pylab import *


def ldoscalc():
    dostens = np.load("dostens.npy")
    cdtens = np.load("cdtens.npy")
    r = np.load("rvec.npy")
    m0 = len(cdtens)
    r0 = len(cdtens[0])
    i0 = len(dostens[0,0])
    n0 = len(dostens[0])
    cdmat = np.zeros((r0,n0))
    dosmat = np.zeros((n0,i0))
    ldosmat = np.zeros((r0,i0))
    for m in range (0, m0):
        print "Calculating ldos for:", (m+1), "/", m0
        cdmat = cdtens[m,:,:]
        dosmat = dostens[m,:,:]
        ldosmat += np.dot(cdmat, dosmat)
    ldosmat = np.transpose(ldosmat)
    ldosmat = ldosmat[:,:] / r
    np.savetxt("ldosmat.txt",ldosmat)
    return ldosmat

def ldosplot(ldosmat):
    r = np.load("rvec.npy")
    E = np.load("Evec.npy")
    print "Producing colour plot..."
    if True:
        pcolor(r,E,ldosmat, vmin=0.0, vmax=0.10)
        colorbar()
#        xvals = arange (0.3, 10.0, 0.1)
#        yvals = 5.0/xvals
#        plot (xvals, yvals, 'w--')
#        plot (xvals, -yvals, 'w--')
        ylim(-10.0, 10.0)
        xlim(0.0, 12.5)
        show()
        figure()

    for rs in [1.0, 2.0, 3.0, 4.0]:
        si = int (rs / 25.0 * 500)
        ri = r[si]
        print rs, si, ri
        
        plot (E, ldosmat[:, si], label='r = %g' % ri) 
        ra = 0.5
        rb = 1.5
        ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
        xvals = [E[t] for t in ivals]
        yvals = [ldosmat[t, si] for t in ivals]
        plot (xvals, yvals, label='slice %f' %si)
        grad = np.polyfit(xvals, yvals, 1)
        print grad, "for r =", ri 
    xlim(-2.5, 2.5)
    legend()
    show()


    return 0

if __name__ == '__main__':
    ldosmat =  ldoscalc()
    ldosplot(ldosmat)
    
