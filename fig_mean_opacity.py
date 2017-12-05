#! /usr/bin/env python

import sys
import os
import numpy as np
import scipy.interpolate as si
import matplotlib as mpl
import matplotlib.lines as mlines
import scipy.constants as sc
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

sys.path.append("../code")
import colormaps as cm
import sma       as sma

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.atmosphere as pa
import pyratbay.constants  as pc
import pyratbay.tools      as pt
import pyratbay.wine       as pw


def Bwn(wn, temp):
  """
  Planck function as function of wavenumber in mks units.
 
  Parameters:
  -----------
  nu: [Type] 1D ndarray, [Units] m-1.
        Frequencies to sample the Planck function.
  temp: [Type] Scalar, [Units] Kelvin degrees.
        Temperature of the blackbody.
        
  Returns:
  --------
  Bnu: [Type] 1D ndarray, [Units] W m^-2 sr^-2 m
       The planck function for temperature temp evaluated at wavenumber wn.
  """
  return 2 * sc.h * sc.c**2 * wn**3 / (np.exp(sc.h*sc.c*wn/(sc.k*temp)) - 1)


def Uwn(wn, temp):
  """
  Temperature derivative of the Plank function.

  Parameters
  ----------
  nu: 1D float ndarray
     Frequencies to sample the Planck function (m-1).
  temp: Float
     Temperature of the blackbody (K).

  Returns
  -------
  dB/dT: 1D float ndarray
     Temperature derivative of the Planck function (W sr-1 m-2 m K-1).
  """
  return 2*sc.h**2 * sc.c**3 * wn**4 / (sc.k * temp**2) * \
     np.exp(sc.h*sc.c*wn/(sc.k*temp)) / (np.exp(sc.h*sc.c*wn/(sc.k*temp))-1)**2



# Load model with LBL opacities from all molecules:
pyrat = pb.pyrat.init("./hot_neptune.cfg")
press   = pyrat.atm.press
wn      = pyrat.spec.wn
wl      = 1e4/pyrat.spec.wn
nlayers = pyrat.atm.nlayers

metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Load grid:
z = 1
grid1 = np.load("../run02_grid/MRT_{:s}.npz".format(metal[z]))
Rp, Mp, Teq = grid1["Rp"], grid1["Mp"], grid1["Teq"]
p0 = grid1["p0"]/pc.bar

grid2 = np.load("../run02_grid/rtransit_{:s}.npz".format(metal[z]))
pt = grid2["pt"]/pc.bar

nr, nm, nt, nf = np.shape(p0)

smaxis = sma.sma(Teq, pyrat.phy.mstar/pc.msun) * pc.au
pyrat.od.maxdepth = np.inf  # Make sure to go all the way down
idx  = pyrat.obs.bandidx[0]
band = pyrat.obs.bandtrans[0]

# Set mass and radius (though these are irrelevant):
imass = 39
irad  = 20

pyrat.phy.mplanet = Mp[imass] * pc.mearth
pyrat.phy.rplanet = Rp[irad]  * pc.rearth


kp1 = np.zeros((nt, nlayers))
kp2 = np.zeros((nt, nlayers))
kr1 = np.zeros((nt, nlayers))
kr2 = np.zeros((nt, nlayers))

for i in np.arange(nt):
  spec, p, temp, q = pa.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                    format(metal[z], Teq[i]))
  spec = list(spec)
  pyrat.phy.smaxis = smaxis[i]
  pyrat.phy.rhill  = pyrat.phy.smaxis * (pyrat.phy.mplanet /
                                         (3*pyrat.phy.mstar))**(1.0/3.0)
  pyrat.refpressure = p0[irad,imass,i,0] * pc.bar
  pyrat = pb.pyrat.run(pyrat, [temp, q])
  # Mass-density profile (g cm-3):
  rho = pyrat.atm.mm*pc.amu * press/(pc.k*Teq[i])
  # B_wn(T):
  B = Bwn(wn*pc.m, Teq[i]) * 1e5
  # dB/dT:
  u = Uwn(wn*pc.m, Teq[i]) * 1e5
  for j in np.arange(nlayers):
    kappa = pyrat.od.ec[j] / rho[j]
    # Integrate:
    Bint1 = np.trapz(B, wn)
    Bint2 = np.trapz(B[idx]*band, wn[idx])
    uint1 = np.trapz(u, wn)
    uint2 = np.trapz(u[idx]*band, wn[idx])
    # Full range:
    kp1[i,j] =    np.trapz(B*kappa, wn) / Bint1
    kr1[i,j] = 1/(np.trapz(u/kappa, wn) / uint1)
    # Kepler band:
    kp2[i,j] =    np.trapz(B[idx]*band*kappa[idx], wn[idx]) / Bint2
    kr2[i,j] = 1/(np.trapz(u[idx]*band/kappa[idx], wn[idx]) / uint2)


palette = cm.viridis_r
lw = 1.5
fs = 16
fig = plt.figure(4)
plt.clf()
plt.subplots_adjust(0.12, 0.15, 0.82, 0.95)
ax = plt.subplot(111)
n = [0, 12, 24, 36, 48, 60, 72, 84, 96, 108, 119]
for j in n:
  plt.semilogy(Teq, kp2[:,j], ".:", color=palette(int(255*j/nlayers)), lw=1.5)
  plt.semilogy(Teq, kr2[:,j], ".-", color=palette(int(255*j/nlayers)), lw=1.5)
ax.set_yticks([1e-10, 1e-8, 1e-6, 1e-4, 1e-2,1e0, 1e2, 1e4, 1e6, 1e8])
plt.ylim(1e-11, 1e2)
plt.ylim(1e-5, 1e3)
plt.xlim(250, 3050)
plt.xlabel(r"$\rm Temperature\ \ (K)$",     fontsize=fs)
plt.ylabel(r"${\rm Mean\ opacity\ \ (cm^{2}g^{-1})}$", fontsize=fs)
# Legend:
plt.plot([],[], ".:", color="k", lw=2, label=r"$\kappa_{\rm Planck}$")
plt.plot([],[], ".-", color="k", lw=2, label=r"$\kappa_{\rm Rosseland}$")
plt.legend(loc="upper left", fontsize=fs+2)
# Pressure color bar:
tlabs = [r"$10^{-8}$", r"$10^{-7}$", r"$10^{-6}$",
         r"$10^{-5}$", r"$10^{-4}$",  r"$10^{-3}$",
         r"$10^{-2}$", r"$10^{-1}$",  r"$10^{0}$",
         r"$10^{1}$", r"$10^{2}$"]
ax2 = fig.add_axes([0.85, 0.15, 0.025, 0.8])
norm = mpl.colors.Normalize(vmin=0, vmax=255)
ticks = np.array(n)*255.0/nlayers
tl = ticks - ticks[1]/2
ticks[0] = 1e-4
tl[0] = 5e-5
cb1 = mpl.colorbar.ColorbarBase(ax2, cmap=palette,
      boundaries=[0]+list(ticks), norm=norm, orientation='vertical')
cb1.ax.tick_params(labelsize=16)
cb1.set_label(r"$\rm Pressure\ \ (bar)$", fontsize=fs)
cb1.set_ticks(tl, True)
cb1.set_ticklabels(tlabs, True)
cb1.ax.invert_yaxis()
plt.savefig("../figs/mean_opacity.ps")
