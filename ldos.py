import numpy as np
import math
import scipy
from scipy import special
import matplotlib
from matplotlib import pylab
from pylab import *


def ldoscalc():
    dostens = np.load("dostens.npy")
    cdtens = np.load("cdtens.npy")
    r = np.load("rvec.npy")
    nm = len(cdtens)
    nr = len(cdtens[0])
    ni = len(dostens[0,0])
    nn = len(dostens[0])
    cdmat = np.zeros((nr,nn))
    dosmat = np.zeros((nn,ni))
    ldosmat = np.zeros((nr,ni))
    for m in range (0, nm):
        print "Calculating ldos for:", (m+1), "/", nm
        cdmat = cdtens[m,:,:]
        dosmat = dostens[m,:,:]
        ldosmat += np.dot(cdmat, dosmat)
    ldosmat = np.transpose(ldosmat)
    ldosmat = ldosmat[:,:] 
    np.save("ldosmat",ldosmat)
    return ldosmat

def ldosplot(ldosmat):
    mlist = np.load("mlist.npy")
    cdtens = np.load("cdtens.npy")
    r = np.load("rvec.npy")
    dr = np.load("drvec.npy")
    nm = len(mlist)
    nldosmat = 2 * np.pi * dr * r * ldosmat
    E = np.load("Evec.npy")
    ldostot = np.zeros((len(E)))
    nu = np.zeros((len(E)))
    nu_m = np.zeros((len(E)))
    print "Producing colour plot..."
    if False:
        pcolor(r,E,ldosmat, vmin=0.0, vmax=0.15)
        colorbar()
        ylim(-10.0, 10.0)
        xlim(0.0, 12.5)
        show()
        figure()

    for r1 in range (0, len(r)):
        ldostot[:] += nldosmat[:, r1]
    figure()
    plot(E, ldostot, label='Global Density of States')
    legend()
    rq = 0.1
    rp = 0.5
    ivals = [t[0] for t in enumerate(E) if t[1] > rq and t[1] < rp]
    xvals = [E[t] for t in ivals]
    yvals = [ldostot[t] for t in ivals]
    grad1 = np.polyfit(xvals, yvals, 1)
    print "gradient of ldost", grad1[0]
    show()

    for rs in [1.0, 2.0, 3.0, 4.0]:
        si = int (rs / 25.0 * 300)
        ri = r[si]
        plot (E, ldosmat[:, si], label='r = %g' % ri) 
        ra = 0.5
        rb = 1.5
        ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
        xvals = [E[t] for t in ivals]
        yvals = [ldosmat[t, si] for t in ivals]
#        plot (xvals, yvals, label='slice %f' %si)
        grad = np.polyfit(xvals, yvals, 1)
        print grad[0], "for r =", ri 
    xlim(-2.5, 2.5)
    legend()
    #show()

    if True:
#        figure()
        print "Calculating analytical LDOS"
        for rs in [1.0, 2.0, 3.0, 4.0]:
            si = int (rs / 25.0 * 300)
            ri = r[si]
            kr = abs(E) * ri
            nu = zeros ((shape(E)))
            for m in range (0,nm):
                jm = special.jn(abs(mlist[m]),kr)
                jm1 = special.jn(abs(mlist[m]+1),kr)
                nu_m = abs(E) / 2.0 / np.pi * (jm**2 + jm1**2)
                nu += nu_m
            plot(E, nu, label='Theory Slice r= %g' %ri)
            ra = 0.1
            rb = 0.5
            ivals = [t[0] for t in enumerate(E) if t[1] < rb and t[1] > ra]
            xvals = [E[t] for t in ivals]
            yvals = [nu[t] for t in ivals]
            grad = np.polyfit(xvals, yvals, 1)
            print grad[0], "for theory r =", ri
    xlim(-2.5, 2.5)
    legend()
    show()

    return 0

if __name__ == '__main__':
    ldosmat =  ldoscalc()
    ldosplot(ldosmat)
    
