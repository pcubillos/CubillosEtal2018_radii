#! /usr/bin/env python
import sys
import numpy as np
import scipy.constants as sc
import matplotlib as mpl
import matplotlib.pyplot as plt
plt.ioff()

sys.path.append("./code")
import colormaps as cm

sys.path.append("./pyratbay")
import pyratbay.constants as pc


# Load grid:
metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
grid = np.load("run02_grid/MRT_{:s}.npz".format(metal[0]))
Rp, Mp, Teq = grid["Rp"], grid["Mp"], grid["Teq"]

mh = sc.m_p + sc.m_e
R = Rp [:,np.newaxis,np.newaxis]
M = Mp [np.newaxis,:,np.newaxis]
T = Teq[np.newaxis,np.newaxis,:]
Lambda = 0.1 * M*pc.mearth * sc.G*mh/(sc.k*T*R*pc.rearth) < 10

# Xaxis
xlabel = [r"${\rm Radius}\ (R_{\oplus})$",
          r"${\rm Mass}\ (M_{\oplus})$",
          r"${\rm Temperature\ (K)}$"]
X  = [Rp, Mp, Teq]
Xl = [[0,41], [0.9, 650], [250, 3050]]

# Colorbar
irad  = [0, 2, 6, 18, 38, 78]
imass = [0,  8, 18, 24, 31, 39, 44]
itemp = [0,  5, 11, 16, 21, 27]
clabel = [r"${\rm Mass}\ (M_{\oplus})$",
          r"${\rm Radius}\ (R_{\oplus})$",
          r"${\rm Mass}\ (M_{\oplus})$"]
palette = [cm.viridis_r, cm.magma_r, cm.viridis_r]
index = [imass, irad, imass]
W     = [Mp, Rp, Mp, Teq]
Wlen  = np.array([len(Mp), len(Rp), len(Mp)])
c0 = np.array([5.0, 30.0, 5.0])
c1 = (255-c0)/(Wlen-1.0)
def col(i,t):
    return int(c0[t] + c1[t]*i)

# Slice
llabel = [r"$T = {}\ {{\rm K}}$",
          r"$T = {}\ {{\rm K}}$",
          r"$R_{{\rm p}}={}\ R_{{\oplus}}$"]
R0, T0 = 8, 7
Z  = [Teq[T0], Teq[T0], Rp[R0]]

# Plot setup:
fs = 15
lw = 1.5
fname = ["photo", "transit"]
ylabel = [r"$\rm Photospheric\ pressure\ \ (bar)$",
          r"$\rm Transit\ pressure\ \ (bar)$"]
yran = [10, 1e-3]


for z in np.arange(len(metal)):
  grid1 = np.load("run02_grid/MRT_{:s}.npz".format(metal[z]))
  grid2 = np.load("run02_grid/rtransit_{:s}.npz".format(metal[z]))
  p0 = [grid1["p0"]/pc.bar, grid2["pt"]/pc.bar]

  npanels = len(X)
  for r in [0,1]:
    # Slice the grid pressures:
    p1 = p0[r][:,:,:,0]
    Y = [p1[:,:,T0].swapaxes(0,1), p1[:,:,T0], p1[R0,:,:]]
    M = Lambda[:,:,T0].swapaxes(0,1), Lambda[:,:,T0], Lambda[R0,:,:]
    # The plot:
    fig = plt.figure(-6+r, (7,11))
    plt.clf()
    plt.subplots_adjust(0.12, 0.08, 0.84, 0.97, hspace=0.27)
    for t in np.arange(npanels):
      ax = plt.subplot(npanels, 1, 1+t)
      for k in np.arange(0, Wlen[t], 3):
        plt.plot(X[t], Y[t][k], lw=lw, color=palette[t](col(k,t)))
      # Boil-off models:
      #for k in np.arange(0, Wlen[t], 3):
      #  m = M[t][k]
      #  if np.sum(m) > 0:
      #    plt.plot(X[t][m], Y[t][k][m], "x", mew=1.25, color="0.6", zorder=0)
      plt.xlim(Xl[t])
      plt.ylim(yran)
      if t == 1:
        plt.ylabel(ylabel[r], fontsize=fs+2)
      if len(X[t])==45:
        ax.set_xscale("log")
      ax.set_yscale("log")
      plt.xlabel(xlabel[t], fontsize=fs)
      # Legend
      plt.plot([],[], color="k", lw=lw, label=llabel[t].format(Z[t]))
      plt.legend(loc="upper right", fontsize=fs-2)
      # Colorbar
      bb = ax.get_position()
      ax2 = fig.add_axes([0.86, bb.y0, 0.025, bb.height])
      norm = mpl.colors.Normalize(vmin=0, vmax=255)
      c = c0[t] + (255-c0[t])*np.arange(len(index[t]))/(len(index[t])-1)
      ticks = np.array(c, np.double)
      tl = ticks - (ticks[2]-ticks[1])/2
      tl[0] = ticks[0]/2
      cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=palette[t],
            boundaries=[0]+list(ticks),
            norm=norm, orientation='vertical')
      cb1.set_label(clabel[t], fontsize=fs)
      cb1.set_ticks(tl, True)
      if len(W[t]) == 45:
        cb1.ax.invert_yaxis()
      tlabs = np.array(np.round(W[t][index[t]]), int)
      cb1.set_ticklabels(tlabs, True)
    plt.savefig("figs/pressure_{:s}_{:s}.ps".format(fname[r], metal[z]))
