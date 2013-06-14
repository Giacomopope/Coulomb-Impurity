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
            piplus  = kfabs - kfabs / 2.0 * np.sqrt( 1 - ((2.0*kfabs / q)**2))
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
    plot(r, func, label='polarisation operator')
    return func

#
# New version of linear response code which uses 
# the Fourier transform of the potential as an input
#
def polaris_generic(r, E0, Uq):
    nr = len(r)
    func = np.zeros(nr)
    for a in range(0,nr):
        def f(q):
            Pi_0 = Pi_mono( E0, q )        
            J = special.jn( 0, q * r[a] )
            return  q * J * Uq(q) * Pi_0
        func[a], eps = integrate.quad(f,0,Inf, limit=100)
    func = func / (2*np.pi)
    #plot(r, func, label='polarisation operator')
    return func

if False:
   #2 window Coulomb potential with r0
   r0 = 1.0
   rvals = np.arange(0.1, 50.0, 0.2)
   def Uq_r0 (q):
       return 2.0 * np.pi / q * np.exp(-q*r0)
   E_2 = -1.0
   E_1 = -0.5
   rhovals2 = polaris_generic(rvals, E_1, Uq_r0) - polaris_generic(rvals, E_2, Uq_r0)  
 
   rhovals0 = polaris(rvals, (E_2, E_1), 1.0, r0)
   plot (rvals, rhovals2, label="New")
   plot (rvals, rhovals0, label="Old")
   legend()
   show()

if False:
   #1 window Coulomb potential with r0
   r0 = 1.0
   E_F = 0.1
   rvals = np.arange(0.1, 50.0, 0.2)
   def Uq_r0 (q):
       return 2.0 * np.pi / q * np.exp(-q*r0)
   rhovals2 = polaris_generic(rvals, E_F, Uq_r0)  
 
   plot (rvals, rhovals2, label="New")
   rho_undoped = 1.0 / 16.0 / np.sqrt(r0**2 + rvals**2)**3
   plot (rvals, rho_undoped, label="Undoped")
   rho_doped = abs(E_F) / 2.0 / np.pi * 1.0 / np.sqrt(r0**2 + rvals**2)
   plot (rvals, rho_doped,   label="Doped")
   legend()
   show()



def drho(r,Elims,Ustr,r0):
    N = len(r)    
    check = np.zeros((len(r)))
    Emax = Elims[1]
    Emin = Elims[0]
    drho = np.zeros((N))
    nm = 15
    r0 = 1.0
#    U0 = -1 / np.sqrt(r**2+r0**2)
#    U0 = -1.0 * np.exp(-0.2 * r)
#    U = Ustr * U0
    U =np.load("potvec-cn-U=%g-grid=%g-r0=%g.npy" %(Ustr,len(r),r0))
    mdrho = np.load("mdrho-cn-U=%g-B=0-m=%d-grid=%d-E=%g-%g.npy" 
                   %(Ustr, nm, N, Emin, Emax))
    for m in range (0,nm):
         drho[:] += mdrho[:,m] 
    difft = - (abs(Emax) - abs(Emin)) * U / (2.0 * np.pi)
    plot(r, drho, label='delta rho - sim')
    plot(r, difft, 'r--', label='delta rho - TF')
    title('Change in charge density; U = %g, Energy window: %g to %g'  
          %(-Ustr, Emin, Emax))
#    ylim(0,0.1)
    legend()
    show()
    figure()
   # loglog(r, 1e-2*r**(-3.0), 'k--', label='r^-3')
   # loglog(r, 1e-2*r**(-2.0), 'm--', label='r^-2')
    loglog(r, 1e-1*r**(-1.0), 'c--', label='r^-1')
    if Ustr >= 0:
        loglog(r, -drho, label='-delta rho - sim')### minus sign!
        loglog(r, -difft, label='-delta rho - TF')### minus sign!
    else:
        loglog(r, drho, label='delta rho - sim')                             
        loglog(r, difft, label='delta rho - TF')
    title('U = %g/r Log-Log' %(-Ustr))
#    legend()
#    show()
    return drho


if __name__ == '__main__':
    r = np.load('rvec.npy')
    r0 = 1.0
    alpha = 0.2
    Elims = np.load('Elims.npy')
    print 'Elims =', Elims
    #Ustr = [-1.0, -0.9, -0.8, -0.7, -0.6, -0.5, -0.4, -0.3, -0.2, -0.1],
    Ustr = [1.0, 0.1, -0.1, -1.0]#, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    diffs = np.zeros((len(r), len(Ustr)))
    diffp = np.zeros((len(r), len(Ustr)))
    ratios = np.zeros((len(r),len(Ustr)))
    def Uq_Coulomb(q):
        return 2.0 * np.pi / q * np.exp( - q * r0)
    def Uq_exp(q):
        return 2.0 * np.pi * alpha / (q**2 + alpha**2)**(1.5)
    for a in range (0,len(Ustr)):
        U =  Ustr[a]
        print "Checking: ", -U
        figure()
        #diffp[:,a] = polaris(r, Elims, U, r0)
        if True:
            diffp[:,a] = U*(polaris_generic(r, Elims[1], Uq_Coulomb)
                            - polaris_generic(r, Elims[0], Uq_Coulomb))
            plot(r, diffp[:,a], label='polarisation operator')
        diffs[:,a] = drho(r, Elims, U, r0)
        ratios[:,a] = (diffs[:,a]/ U) - (diffp[:,0] / Ustr[0])
        if U >= 0:
            loglog(r, -diffp[:,a], label='-Polarisation Operator')### minus sign!
        else:
            loglog(r, diffp[:,a], label='Polarisation Operator')
        legend()
        show()
    figure()
    for a in range (0,len(Ustr)):
        U =  Ustr[a]
        U1 = np.load("potvec-cn-U=%g-grid=%g-r0=%g.npy" %(U,len(r),r0))
        plot(r, (diffs[:,a]/U), label='U = %g / r Sim' %(-U))
    plot(r, (diffp[:,0]/Ustr[0]), label='Linear response U = %g / r' %(-U))
    title("Change in Charge Density Over Inducing Potential")
    legend()
    show()
    np.save("cou-linear-array", diffp)
   # np.save('ratios-N=%d-m=1-E=%g-%g' %(len(r), Elims[0], Elims[1]), ratios)
    figure()
   # plot([r[0],r[-1]], [1,1], 'k--')
    for b in range (0,len(Ustr)):
        plot(r, ratios[:,b], label='U = %g' %Ustr[b])
        title('Difference Simulation - Polarisation')
        legend()
    show()
    fluff = diffp[:,a] / U
    np.save('generic-linear-response-cn-r0=1.0', fluff)
