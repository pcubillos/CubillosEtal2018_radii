#! /usr/bin/env python
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.constants   as sc
plt.ioff()

sys.path.append("./code")
import colormaps as cm

sys.path.append("./pyratbay")
import pyratbay.constants as pc

mh = sc.m_p + sc.m_e


# Load grid data:
grid = np.load("run02_grid/MRT_solar.npz")
Rp, Mp, Teq = grid["Rp"], grid["Mp"], grid["Teq"]
ptransit = np.load("run02_grid/rtransit_solar.npz")["pt"]/pc.bar

cont = np.copy(ptransit[:,:,0,0])
cont[cont>0] = 1.0

dr = Rp[1] - Rp[0]
Rplot = np.concatenate((Rp, [Rp[-1]+dr])) - dr/2.0
logM = np.log(Mp)
dlm = logM[1] - logM[0]
Mplot = np.exp(np.concatenate((logM[0:1]-dlm, logM)) + dlm/2)

palette = cm.inferno_r
palette.set_under(color=(1,1,1,0))
palette.set_bad(color='w')

lw = 1.25
fs = 12
matplotlib.rcParams.update({'ytick.labelsize':10})
matplotlib.rcParams.update({'xtick.labelsize':10})

itemp = [0, 4, 8, 12, 16, 20, 24, 27]

plt.figure(15, (8.5, 4))
plt.clf()
plt.subplots_adjust(0.06, 0.13, 0.865, 0.97, hspace=0.04, wspace=0.04)
for i in np.arange(len(itemp)):
  it = itemp[i]
  Z = np.copy(ptransit[:,:,it,0])
  Zmin = 0.003
  Zmax = 0.3
  Z[(Z>0) & (Z<Zmin)] = Zmin
  ax = plt.subplot(2,4,i+1)
  plt.contour(Mp, Rp, cont, 1, colors="0.7", linewidths=0.75, zorder=0)
  plt.pcolor(Mplot, Rplot, np.log10(Z), edgecolors="face",
             vmin=np.log10(Zmin), vmax=np.log10(Zmax), cmap=palette)
  rlambda10 = 0.1*sc.G*(Mp*pc.mearth)*mh/(10.0*sc.k*Teq[it])/pc.rearth
  rlambda20 = 0.1*sc.G*(Mp*pc.mearth)*mh/(20.0*sc.k*Teq[it])/pc.rearth
  plt.plot(Mp, rlambda10, "green",     lw=lw)
  plt.plot(Mp, rlambda20, "limegreen", lw=lw)
  #Rh = Mp**(1/3.) * (Teq[it]/(1.75*3000))**-2
  #plt.plot(Mp, Rh, "-", color="orangered", lw=lw)
  ax.set_xscale('log')
  if i%4 == 0:
    plt.ylabel(r"${\rm Radius}\ (R_{\oplus})$", fontsize=fs)
  else:
    ax.set_yticklabels([])
  if i>=4:
    plt.xlabel(r"${\rm Mass}\ (M_{\oplus})$",   fontsize=fs)
  else:
    ax.set_xticklabels([])
  plt.xlim(Mplot[0], Mplot[-1])
  plt.ylim(Rplot[0], Rplot[-1])
  plt.text(1.1, 36, r"$T=\ {:4.0f}\ {{\rm K}}$".format(Teq[it]),
           fontsize=fs)
cax = plt.axes([0.87, 0.13, 0.016,  0.84])
cax = plt.colorbar(cax=cax)
cax.set_label(r"${\rm Transit\ pressure}\ \ ({\rm bar})$", fontsize=fs)
cbticks = [3e-1, 1e-1, 3e-2, 1e-2, 3e-3]
cax.set_ticks(np.log10(cbticks))
cax.set_ticklabels([r"$3\times10^{-1}$", r"$1\times10^{-1}$",
          r"$3\times10^{-2}$", r"$1\times10^{-2}$", r"$<3\times10^{-3}$"])
cax.ax.invert_yaxis()
cax.solids.set_edgecolor("face")
plt.savefig("./figs/transit_pressure_lambda.ps")
