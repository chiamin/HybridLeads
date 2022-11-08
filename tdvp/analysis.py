import pylab as pl
from collections import OrderedDict
import sys
from math import pi, acos
import plotsetting as ps
import numpy as np
import fitfun as ff
import matplotlib.colors as colors
import cmasher as cmr

def exactG (V, tp):
    if V != 0:
        return 0.5 * pi / acos(-0.5*V)
    else:
        return 4*tp**2/(1+tp**2)**2

def get_para (fname, key, typ):
    with open(fname) as f:
        for line in f:
            if key in line:
                return typ(line.split()[-1])

def get_data (fname):
    dt = get_para (fname, 'dt', float)
    t_contact = get_para (fname, 't_contact', float)
    muL = get_para (fname, 'mu_leadL', float)
    muR = get_para (fname, 'mu_leadR', float)
    Vg = muR - muL

    I_tx = OrderedDict()
    dens = OrderedDict()
    iL_dev = OrderedDict()
    iR_dev = OrderedDict()
    with open(fname) as f:
        for line in f:
            if 'step =' in line:
                step = int(line.split()[-1])
                I_tx[step] = OrderedDict()
                I_t = I_tx[step]
                dens[step] = OrderedDict()
                ns = dens[step]
            elif 'device site =' in line:
                idevL = int(line.split()[-2])
                idevR = int(line.split()[-1])
                iL_dev[step] = idevL
                iR_dev[step] = idevR
            elif 'current ' in line:
                tmp = line.split()
                ilink = int(tmp[1])

                I = float(tmp[-1].strip('()').split(',')[0])
                cc = 2*pi/Vg
                if ilink == idevL-1 or ilink == idevR:
                    cc *= t_contact
                I_t[ilink] = I * cc
            elif '\tn ' in line:
                tmp = line.split()
                i = int(tmp[1])
                n = float(tmp[-1])
                ns[i] = n

    xss, tss, Iss, nss = [],[],[],[]
    I_L, I_R, I_mean, ts = [],[],[],[]
    w = 1
    for step in I_tx:
        t = step*dt
        Is = list(I_tx[step].values())
        xs = list(I_tx[step].keys())
        central_site = int((iL_dev[step] + iR_dev[step])/2)
        xss += [i - central_site for i in xs]
        tss += [t for i in range(len(Is))]
        Iss += Is
        nss += dens[step].values()

        iL, iR = iL_dev[step]-1, iR_dev[step]
        if iL < len(Is) and iR < len(Is):
            I_L.append (Is[iL])
            I_R.append (Is[iR])
            whf = int(w/2)
            I_mean.append (np.mean(Is))#[imp_site-whf:imp_site+whf+1]))
            ts.append (t)

    return xss, tss, Iss, nss, ts, I_L, I_R, I_mean

def get_average_current (ts, js, ax=None, tbeg=0., tend=float('inf'), **args):
    ibeg, iend = 0, len(js)
    for i in range(len(js)):
        if ts[i] > tbeg:
            ibeg = i
            break
    for i in range(len(js)):
        if ts[i] > tend:
            iend = i
            break
    G = np.mean (js[ibeg:iend])
    err = np.std (js[ibeg:iend])

    if ax != None:
        ax.plot (ts, js, **args)
        ax.axhline (G, ls='--', c='k')
        ax.axvline (tbeg, ls=':', c='gray')
    return G, err

def to_imag_data (xss, tss, Iss):
    tmap = dict()
    ti = 0
    tpre = np.nan
    for t in tss:
        if t != tpre:
            tmap[t] = ti
            tpre = t
            ti += 1
    Nt = len(tmap)
    tmin = min(tss)
    tmax = max(tss)

    xmin = min(xss)
    xmax = max(xss)
    Nx = xmax-xmin+1

    Idata = np.empty ((Nt,Nx))
    Idata.fill (np.nan)
    for x,t,I in zip(xss,tss,Iss):
        ti = tmap[t]
        xi = x-xmin
        Idata[ti,xi] = I
    return Idata, tmin, tmax, xmin, xmax

def plot_data (fname, ax2, tbeg):
    with open(fname) as f:
        for line in f:
            if 'Largest link dim' in line:
                m = int(line.split()[-1])

    xss, tss, Iss, nss, ts, I_Lt, I_Rt, I_mean = get_data (fname)

    print (fname, Iss[-1])

    f1,ax1 = pl.subplots()
    Idata, tmin, tmax, xmin, xmax = to_imag_data (xss, tss, Iss)
    Ilim = max([abs(i) for i in Iss])
    midnorm = ps.MidpointNormalize (vmin=-Ilim, vmax=Ilim, vzero=0.)
    sc = ax1.imshow (Idata, origin='lower', extent=[xmin, xmax, tmin, tmax], aspect='auto', cmap='cmr.pride_r',norm=midnorm)
    #sc = ax1.scatter (xss, tss, c=Iss)
    cb = pl.colorbar (sc)
    ax1.set_title ('$m='+str(m)+'$')
    ax1.set_xlabel ('site')
    ax1.set_ylabel ('time')
    cb.ax.set_ylabel ('current')

    f,ax = pl.subplots()
    Idata, tmin, tmax, xmin, xmax = to_imag_data (xss, tss, nss)
    sc = ax.imshow (Idata, origin='lower', extent=[xmin, xmax, tmin, tmax], aspect='auto')
    #sc = ax.scatter (xss, tss, c=nss)
    cb = pl.colorbar (sc)
    ax.set_title ('$m='+str(m)+'$')
    ax.set_xlabel ('site')
    ax.set_ylabel ('time')
    cb.ax.set_ylabel ('density')

    ps.set([ax, ax1])

    IL, errL = get_average_current (ts, I_Lt, ax=ax2, tbeg=tbeg, marker='.')
    IR, errR = get_average_current (ts, I_Rt, ax=ax2, tbeg=tbeg, marker='.')

    muL = get_para (fname, 'mu_leadL', float)
    muR = get_para (fname, 'mu_leadR', float)
    Vg = muR - muL
    GL = IL / Vg
    GR = IR / Vg
    GerrL = errL / Vg
    GerrR = errR / Vg

if __name__ == '__main__':
    fc,axc = pl.subplots()
    tbeg = 20

    files = [i for i in sys.argv[1:] if i[0] != '-']
    for fname in files:
        print (fname)

        plot_data (fname, axc, tbeg)

        if '-pdf' in sys.argv:
            f1.savefig (fname+'_I.pdf')
    #axc.set_xlim (0, 30)
    axc.set_xlabel ('time')
    axc.set_ylabel ('current')
    axc.legend()
    ps.set(axc)
    if '-pdf' in sys.argv:
        fc.savefig ('G.pdf')
    pl.show()
