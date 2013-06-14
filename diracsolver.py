import numpy as np
import scipy
from scipy import linalg
import cmath
import math
from math import *
from pylab import *
from scipy import special
from ldos import *
import time


def diracham(r,pot,mlist,B0):
    N = len(r)
    b = len(mlist)
    H = np.zeros((2*N,2*N), dtype=complex)
    P = np.zeros((2*N,2*N), dtype=complex)
    M = np.zeros((2*N,2*N), dtype=complex)
    U = np.zeros((2*N,2*N), dtype=complex)
    Ematp = np.zeros((2*N,b))
    Ematn = np.zeros((2*N,b))
    cdtens = np.zeros((b,N,2*N))
    cdtensp = np.zeros((b,N,2*N))
    cdtensn = np.zeros((b,N,2*N))
    psi_up = np.zeros((N))
    psi_down = np.zeros((N))
    modpsi = np.zeros((N))
    totmodpsi = np.zeros((N))
    dr = np.zeros((N))
    dr[1:] = r[1:] - r[0:-1]
    dr[0] = r[1] - r[0]
    np.save("drvec", dr)
    info = np.load("EMinfo.npy")
    j = 1j
    timestart1 = time.time() 

    for B in [B0, -B0]:
        for i_m, m in enumerate(mlist): 
            M[:, :] = 0.0
            print "Calculating Momentum Channel:", m
            for y in range (0,N):
                if y == 0:
                    a = r[y+1] - r[y]
                    P[0,1] = -1.0j / a
                    P[1,0] = 1.0j /a
                else:
                    a = r[y] - r[y-1]
                    P[2*y+1,2*y]= 1.0j /a   ####FIX THIS
                    P[2*y,2*y+1]= -1.0j / a
                    P[2*y,2*y-1]= 1.0j /a
                    P[2*y-1,2*y]= -1.0j / a
                # Calculating the m/r term. We split this 
                # into two steps, to allow for coupling to the 
                # left and right value of psi_1 --- AVS 
                if (y != N - 1): 
                    r_1 = r[y]
                    r_2 = r[y + 1]
                    y_next = y + 1
                    #y_next = y
                    #r_2 = r[y]
                else:
                    r_1 = r[-2]
                    r_2 = r[-1]
                    y_next = y
                r_mid = r_2 #(r_1 + 3.0*r_2)/4.0
                m_eff = m + 0.5 + B * r_mid**2 / 2.0
                M[2*y, 2*y_next + 1] += -0.0j * m_eff / r_mid
                #print "**", y, M[2 * y, 2*y_next + 1], 2*y, 2*y_next + 1
                M[2*y_next + 1, 2*y] = M[2*y, 2*y_next + 1].conj()
                #print M[2*y, 2*y_next + 1], M[2 * y_next + 1, 2 * y]
      
                if (y < N - 1): 
                    r_1 = r[y]
                    r_2 = r[y + 1]
                else:
                    r_1 = r[-2]
                    r_2 = r[-1]
                r_mid = r_1
                #r_mid = (3.0 * r_1 + r_2) / 4.0

                if y < 6 :
                    delta_m = np.sqrt(r_1) / (np.sqrt(r_1) + np.sqrt(r_2))
                else:
                    delta_m = 0.5
                m_eff = m + delta_m + B * r_mid**2 / 2.0     
                #print m_eff, r_mid
                M[2*y, 2*y + 1] += -1.0j * m_eff / r_mid
                M[2*y + 1, 2*y] = M[2*y, 2*y + 1].conj()
                #print M[2 * y, 2*y + 1], 2 * y, 2 * y + 1
                #print M[2*y, 2*y + 1], M[2 * y + 1, 2 * y]
                U[2*y, 2*y] = pot[y]
                U[2*y +1, 2*y +1] = pot[y]
                H = P + M + U
            print "diagonalising... "
            w, vr =  scipy.linalg.eigh(H)
            if B == B0:
                Ematp[:,i_m] = w[:]
            elif B == (-1.0 * B0) and B!= 0.0:
                Ematn[:,i_m] = w[:]
            if False:
                ea = -2
                eb = 2
                ivals = [t[0] for t in enumerate(w) if t[1] > ea and t[1] < eb]
                evals = [w[t] for t in ivals]
                if m == 0:
                    ens = list(evals)
                else:
                    ens.extend(evals)
            print "Resolving wavefunctions"
            iplot = []
            plot_wf = True
            for Eplot in [-0.24, -0.04]:
                Ei = list(enumerate (w))
                Ei.sort (lambda x, y: cmp(abs(x[1]- Eplot), abs(y[1]- Eplot)))
                iplot.append(Ei[0][0])
            for i in range (0,2*N):
                u = vr[:,i]
                u_up = u[::2]
                u_down = u[1::2]
                u_up_real = u_up.real
                u_up_imag = u_up.imag
                u_down_real = u_down.real
                u_down_imag = u_down.imag
                modpsi =(abs(u_up)**2+abs(u_down)**2)/(2*np.pi*r*dr)
                if B == B0:
                    cdtensp[i_m,:,i] = modpsi[:]
                elif B == (-1.0 * B0) and B!=0.0:
                    cdtensn[i_m,:,i] = modpsi[:]
                totmodpsi += modpsi
                if i in iplot and mlist[m]==0:
                #if  i >= (N - 2) and i <= (N + 1):
                    if plot_wf: #mlist[m] < 0:
                        print "Plotting state %d" %i
                        figure()
                        sqr = np.sqrt(r)
                        plot(r,u_up_real/sqr, label='psi up real')
                        plot(r,u_up_imag/sqr, label='psi up imag')
                        plot(r,u_down_real/sqr, label='psi down real')
                        plot(r,u_down_imag/sqr, label='psi down imag')
                        legend()
                        title('Momentum Channel %d' % m)
#                    figure()
#                    plot(r,totmodpsi,label='charge density m %f' %mlist[m])
#                    legend()
#            show()
            if plot_wf: show()
        timeend1 = time.time()
        print "Time taken:", timeend1 - timestart1
        if B == 0:
            break
    cdtens = cdtensp + cdtensn
    np.save("cdtens-cn-U=%g-B=%g-ms=%d-N=%d-dm6" %(info[0],B0,b,N),cdtens)
    np.save("totmodpsi-cn-U=%g-B=%g-ms=%d-N=%d-dm6" %(info[0],B0,b,N), totmodpsi)
    return Ematp, Ematn, cdtens

def DOS(Ematp, Ematn, mlist ,r):
    info = np.load("EMinfo.npy")
    N0 = len(Ematp)
    c = len(mlist)
    N = 4 * N0
    E = np.zeros((N))
    N0min = 0
    N0max = N0-1
    dos = np.zeros((N))
    dostens = np.zeros((c,N0,N))
    doschan = np.zeros((N))
    rmax = r[N0/2 -1.0]
    gam = np.pi * 0.8 / rmax
    Emax = 24.0
    Emin = -Emax
    wf = np.load("cdtens-cn-U=%g-B=%g-ms=%d-N=%d-dm6.npy" %(info[0],B0,len(mlist),len(r)))
    timestart2 = time.time()
    A = 2.0/np.pi*gam**3 
    for i in range (0,N):
        E[i] = Emin + float(i) * (Emax - Emin)/N
    print 'Filling Global Density of States'
    for m in range (0,c):
        mlab = mlist[m]
        for n in range (N0min, N0max):
            if False: # Emat[n,m] < 0.001 and Emat[n,m] > -0.001:          
                print "Zero energy at n element =  %d" %n
                print "m = %d"  %mlist[m]
                print "Eigenvalue:", Emat[n,m]
                rhoo = cdtens[m,:,n]
                figure()
                plot(r, rhoo, label='charge density energy=%f' %Emat[n,m])
                title('Momentum Channel %d' %mlist[m])
                figtext(0.2, 0.85, 'Energy state %d' %n)
                legend()
                show()
            for i in range (0,N): 
                dostens[m,n,i]+=A/(gam**2+(E[i]-Ematp[n,m])**2)**2
                if info[1] != 0:
                    dostens[m,n,i]+=A/(gam**2+(E[i]-Ematn[n,m])**2)**2
 
#    for m in range (0,c):
#        for n in range (N0min, N0max):
#            dos[:] += dostens[m,n,:]
    if info[1] == 0.0:
        dostens = 2.0 * dostens
    np.save("dostens-cn-U=%g-B=%g-ms=%d-N=%d-dm6" 
            %(info[0], info[1], c, len(r)), dostens)
    np.save("globdos(pos)", dos)
    timeend2 = time.time()
    print "Total time:", timeend2 - timestart2 
    return E, dostens
def DOS2(Ematp, Ematn, mlist ,r):
    info = np.load("EMinfo.npy")
    N0 = len(Ematp)
    c = len(mlist)
    N = 4 * N0
    E = np.zeros((N))
    N0min = 0
    N0max = N0
    dos = np.zeros((N))
    dostens = np.zeros((c,N0,N))
    doschan = np.zeros((N))
    rmax = r[N0/2 -1.0]
    gam = np.pi * 0.8 / rmax
    Emax = 24.0
    Emin = -Emax
    wf = np.load("cdtens-cn-U=%g-B=%g-ms=%d-N=%d-dm6.npy" %(info[0],B0,len(mlist),len(r)))
    timestart2 = time.time()
    A = 2.0/np.pi*gam**3 
    for i in range (0,N):
        E[i] = Emin + float(i) * (Emax - Emin)/N
    print 'Filling Global Density of States'
    for i in range(0, N):
        if i % (N/100) == 0:
            print (100*(i+1))/N, "% complete"
        X = A / (gam**2 + (E[i] - Ematp)**2)**2
        if (info[1] != 0):
            X += A / (gam**2 + (E[i] - Ematn)**2)**2        
        dostens[:, :, i] += np.transpose(X)

 
#    for m in range (0,c):
#        for n in range (N0min, N0max):
#            dos[:] += dostens[m,n,:]
    if info[1] == 0.0:
        dostens = 2.0 * dostens
    np.save("dostens-cn-U=%g-B=%g-ms=%d-N=%d-dm6" 
            %(info[0], info[1], c, len(r)), dostens)
    np.save("globdos(pos)", dos)
    timeend2 = time.time()
    print "Total time:", timeend2 - timestart2 
    return E, dostens

if __name__ == '__main__':
   ### dm terminates at element 6
   print "Running Coulomb potential."
   N = 800
   rmin = 0.01
   rmax = 25.0
   alpha = 0.2
   B0 = 0.0 #2.370
   r = zeros((N))
   pot = zeros((N))
   a = 15
   Ustr = 0.0
   info = np.zeros((2))
   info[0] = Ustr
   info[1] = B0
   r0 =1.0
   np.save("EMinfo", info)
#   mlist = zeros((2*a + 1))
   mlist = np.array(range(0,a))
#   mlist[0] = 1
   for i in range (0,N):
       r[i] = rmin +  i*(rmax-rmin) / N
   pot = -Ustr /np.sqrt(r**2 + r0**2) #Coulomb
   #pot = -Ustr * np.exp(-alpha * r) #Exponential
   print "Momentum Channels:",  mlist
   Ematp, Ematn, cdtens = diracham(r, pot, mlist,B0)
   E, dostens = DOS2(Ematp, Ematn, mlist ,r)
   np.save("rvec",r)
   np.save("mlist",mlist)
   np.save("Evec",E)
   np.save("potvec-cn-U=%g-grid=%g-r0=%g" %(Ustr, N, r0),pot)
