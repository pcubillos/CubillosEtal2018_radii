#! /usr/bin/env python
import sys
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
import scipy.optimize  as so
plt.ioff()

sys.path.append("./code")
import colormaps as cm

sys.path.append("./pyratbay")
import pyratbay.constants as pc


# Create the grid:
nt = 28 # Number of temperature values
nr = 79 # Number of planetary radius values
nm = 45 # Number of planetary mass values
Teq = np.linspace(300, 3000, nt)  # K
Rp  = np.linspace(1.0, 40.0, nr)  # Rearth
Mp  = np.logspace(np.log10(1.0), np.log10(630), 45)  # Mearth
Mp[-1]  = 600.0

# Density:
dens = Mp*pc.mearth / (4./3*np.pi*np.expand_dims(Rp*pc.rearth,1)**3)
ldens = np.log10(dens)

# Density plot:
palette = cm.magma_r
palette.set_under(color='w')
palette.set_over(color='w')

# Plot edges:
dr = Rp[1] - Rp[0]
Rplot = np.concatenate((Rp, [Rp[-1]+dr])) - dr/2.0
logM = np.log(Mp)
dlm = logM[1] - logM[0]
Mplot = np.exp(np.concatenate((logM[0:1]-dlm, logM)) + dlm/2)

# The plot:
plt.figure(0, (7, 4.5))
plt.clf()
plt.subplots_adjust(0.1, 0.15, 0.95, 0.95)
ax = plt.subplot(111)
a=plt.pcolor(Mplot, Rplot, ldens, cmap=palette, vmin=np.log10(0.03),
             vmax=np.log10(30.0), edgecolors="face")
ax.set_xscale('log')
plt.xlabel(r"${\rm Mass}\ (M_{\oplus})$", fontsize=14)
plt.ylabel(r"${\rm Radius}\ (R_{\oplus})$", fontsize=14)
# Solar system:
plt.plot([1],     [1],     "o", ms=8, mew=1.25, mec="w", mfc="none")
plt.plot([17.15], [3.883], "o", ms=8, mew=1.25, mec="w", mfc="none")
plt.plot([317],   [11.2],  "o", ms=8, mew=1.25, mec="w", mfc="none")
plt.xlim(Mplot[0], Mplot[-1])
plt.ylim(Rplot[0], Rplot[-1])
# Color bar:
cax = plt.colorbar()
cax.set_label(r"${\rm Bulk\ density}\ \ ({\rm g\ cm}^{-3})$", fontsize=14)
cbticks = [0.03, 0.1, 0.3, 1.0, 3.0, 10.0, 30.0]
cax.set_ticks(np.log10(cbticks))
cax.set_ticklabels(cbticks)
cax.solids.set_edgecolor("face")
cax.ax.invert_yaxis()
plt.savefig("./figs/density.ps")
