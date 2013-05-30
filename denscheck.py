import math
import numpy as np
import scipy
from scipy import special
from scipy import integrate
import matplotlib
from matplotlib import pylab
from pylab import *


def Pi_mono(kf,q):
    #if abs(kf) < 1e-2: return q / 16.0; 
    pimin = np.pi * q / 8.0
    kfabs = abs(kf)
    if abs(q) <= 2.0 * kfabs:
            piplus = kfabs - pimin
    else:
            piplus  = kf - kfabs / 2.0 * np.sqrt( 1 - ((2.0*kfabs / q)**2))
            piplus += -q / 4.0 * np.arcsin( 2.0 * kfabs / q)
    return (pimin + piplus) / 2.0 / np.pi

if False: 
   qvals   = np.linspace(0.0, 10.0, 1000)
   figure()
   for kf in [1.0, 2.0]: 
       pivals = np.vectorize(lambda q: Pi_mono(kf, q)) (qvals)
       plot (qvals, pivals, label='kF = %g' % kf)
   plot (qvals, 1.0 /16.0 * qvals, 'k--', label='kf = 0 limit')
   legend()
   show() 


def polaris(r, Elims,Ustr,r0):
    nr = len(r)
    func = np.zeros(nr)

    for a in range(0,nr):
        def f(q):
            Pi_1 = Pi_mono( Elims[1], q )
            Pi_0 = Pi_mono( Elims[0], q )        
            dPi = Pi_1 - Pi_0
            U = (2 * np.pi) * np.exp( - q * r0 ) / q * Ustr
            J = special.jn( 0, q * r[a] )
            return  q * J * U * dPi
        func[a], eps = integrate.quad(f,0,Inf)
    func = func / (2*np.pi)
    #plot(r, func, label='polarisation operator')
    return func


def drho(r,Elims,Ustr,r0):
    N = len(r)
    check = np.zeros((len(r)))
    Emax = Elims[1]
    Emin = Elims[0]
    rho0 = np.load('charge-density-U=0-B=0-m=1-N=%d-E=%g-%g.npy' %(N, Emin, Emax))
    U0 = -1 /np.sqrt(r**2 + r0**2)
    U = Ustr * U0
    rho = np.load('charge-density-U=%g-B=0-m=1-N=%d-E=%g-%g.npy' %(Ustr,N,Emin,Emax))
    difft = - (Emax - Emin) * U / (2.0 * np.pi)
    diffs = rho - rho0
#    figure()

    if False:
        cdtens0 = np.load("cdtens-U=0-B=0-ms=15-N=%d.npy" %N)
        cdtensU = np.load("cdtens-U=%g-B=0-ms=15-N=%d.npy" %(Ustr, N))

    plot(r, diffs, label='delta rho - sim')
    check = -1.0 / 16  * r0 / (r0**2 + r**2)**(1.5)
   # plot(r,check, label='check')
    title('U =- %g' %Ustr)
    ylim(0,0.1)
    legend()
#    figure()
#    loglog(r, diffs, label='delta rho - sim')
#    loglog(r, difft, label='delta rho - theory')
#    title('U =- %g/r Log-Log' %a)
#    legend()
    return diffs



if __name__ == '__main__':
    r = np.load('rvec.npy')
    r0 = r[10]
    Elims = np.load('Elims.npy')
    print 'Elims =', Elims
    Ustr =[0.002] #[-0.025, 0.025]  #[0.001, 0.01, 0.025, 0.05, 0.075, 0.10, 0.2]
    ratios = np.zeros((len(r),len(Ustr)))
    for a in range (0,len(Ustr)):
        U = Ustr[a]
        figure()
        diffp = polaris(r, Elims, U, r0)
        diffs = drho(r, Elims, U, r0)
        ratios[:,a] = diffp / diffs
    figure()
    np.save('ratios-N=%d-m=1-E=%g-%g' %(len(r), Elims[0], Elims[1]), ratios)
    for b in range (0,len(Ustr)):
        plot(r, ratios[:,b], label='U = %g' %Ustr[b])
        title('Ratio of Polarisation / Simulation')
        legend()
    show()
