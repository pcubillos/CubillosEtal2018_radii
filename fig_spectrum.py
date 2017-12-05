#! /usr/bin/env python
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.interpolate as si
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.tools     as pt
import pyratbay.wine      as pw


# Emission contribution function:
pyrat = pb.pbay.run("hot_neptune.cfg")
cf  = pt.cf(pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
bcf = pt.bandcf(cf, pyrat.obs.bandtrans, pyrat.spec.wn, pyrat.obs.bandidx)

# Pressure-radius interpolation:
rint = si.interp1d(pyrat.atm.press/pc.bar, pyrat.atm.radius/pc.rearth)

# Transmission spectrum and transmittance:
pyrat.od.path = "transit"
pyrat = pb.pyrat.run(pyrat)
bandflux = np.sqrt(pw.bandintegrate(pyrat=pyrat))*pyrat.phy.rstar/pc.rearth
rtrans = gaussf(np.sqrt(pyrat.spec.spectrum)*pyrat.phy.rstar/pc.rearth, 1.0)

transmittance = pt.transmittance(pyrat.od.depth, pyrat.od.ideep)
bcft = pt.bandcf(transmittance, pyrat.obs.bandtrans, pyrat.spec.wn,
                 pyrat.obs.bandidx)

# Plot setup:
matplotlib.rcParams.update({'ytick.labelsize':10})
matplotlib.rcParams.update({'xtick.labelsize':10})
yran = 3.75, 4.83
fs = 12
lw = 1.5
fcol = ["mediumblue", "limegreen", "darkorange", "r"]
h = [np.amax(pyrat.obs.bandtrans[0]),
     np.amax(pyrat.obs.bandtrans[1]),
     np.amax(pyrat.obs.bandtrans[2]),
     np.amax(pyrat.obs.bandtrans[3])]

# The plot:
plt.figure(-10, (8.5, 4))
plt.clf()
# Transmission
ax=plt.axes([0.075, 0.25, 0.47, 0.65])
plt.plot(1e4/pyrat.spec.wn, rtrans, lw=1.0, color="orange")
for j in np.arange(pyrat.obs.nfilters):
  plt.plot(1e4/pyrat.spec.wn[pyrat.obs.bandidx[j]],
           0.16*pyrat.obs.bandtrans[j]/h[j] + yran[0], lw=lw, color=fcol[j])
  plt.plot(1e4/pyrat.obs.bandwn[j], bandflux[j], "o", ms=6,
           color=fcol[j], mec="k", mew=1.0)
plt.ylabel(r"${\rm Transmission\ radius}\ (R_{\oplus})$", fontsize=fs)
plt.xlabel(r"$\rm Wavelength\ (um)$",            fontsize=fs)
plt.ylim(yran)
plt.xlim(0.3, 1.2)
# Transmittance
ax=plt.axes([0.61, 0.25, 0.12, 0.65])
plt.plot(bcft[0], pyrat.atm.radius/pc.rearth, lw=lw, color=fcol[0])
ax.set_xticks([0, 0.5, 1.0])
plt.ylim(yran)
plt.xlim(-0.02, 1.02)
plt.ylabel(r"${\rm Impact\ parameter}\ (R_{\oplus})$", fontsize=fs)
plt.xlabel(r"${\rm Transmittance}$", fontsize=fs)
# Contribution function:
ax=plt.axes([0.81, 0.25, 0.12, 0.65])
ax2 = plt.twinx()
ax2.plot(bcf[0]/np.amax(bcf[0]), pyrat.atm.radius/pc.rearth,
         lw=lw, color=fcol[0], ls="-")
ax2.set_ylim(yran)
ax2.set_xlim(-0.02, 1.12)
ax2.set_xticks([0, 0.5, 1.0])
ax2.set_ylabel(r"${\rm Radius}\ (R_{\oplus})$", fontsize=fs)
# Pressure scale:
ax.set_ylim(yran)
yticks = np.logspace(1, -6, 8)
syticks = [r"$10^{{{:.0f}}}$".format(np.log10(x)) for x in yticks]
ax.set_yticks(rint(yticks))
ax.set_yticklabels(syticks)
ax.set_ylabel(r"$\rm Pressure\ (bar)$", fontsize=fs)
ax.set_xlabel(r"$\rm Contrib.\ function$", fontsize=fs)
plt.savefig("../figs/neptune_spectrum.ps")
