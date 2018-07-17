#! /usr/bin/env python
import sys
import matplotlib
import numpy as np
import scipy.interpolate as si
import scipy.constants   as sc
import matplotlib.pyplot as plt
plt.ioff()

sys.path.append("./code")
import colormaps as cm
import sma       as sma

sys.path.append("./pyratbay")
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa


# The grid:
grid  = np.load("run02_grid/MRT_solar.npz")
Rp, Mp, Teq = grid["Rp"], grid["Mp"], grid["Teq"]
nr, nm, nt  = len(Rp), len(Mp), len(Teq)
smaxis = sma.sma(Teq,1.3) * pc.au

metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Lambda:
mh = sc.m_p + sc.m_e
R = Rp [:,np.newaxis,np.newaxis]
M = Mp [np.newaxis,:,np.newaxis]
T = Teq[np.newaxis,np.newaxis,:]
Lambda = 0.1 * M*pc.mearth * sc.G*mh/(sc.k*T*R*pc.rearth)

# Plot setup:
dens = Mp*pc.mearth / (4./3*np.pi*np.expand_dims(Rp*pc.rearth,1)**3)
mask = (dens>0.03) & (dens <=30.0)

matplotlib.rcParams.update({'ytick.labelsize':14})
matplotlib.rcParams.update({'xtick.labelsize':14})

lab = [r"$0.1\times\ {\rm solar}$",
       r"$1.0\times\ {\rm solar}$",
       r"$10.0\times\ {\rm solar}$",
       r"$100.0\times\ {\rm solar}$"]

col = ["navy", "red", "limegreen", "orange"]

lw = 1.5
fs = 18

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Collect (Rt-Rp)/H:
fmean   = np.zeros((nz,nt))
fmedian = np.zeros((nz,nt))
fstd    = np.zeros((nz,nt))

for z in np.arange(nz):
  grid  = np.load("run02_grid/MRT_{:s}.npz".format(metal[z]))
  tgrid = np.load("run02_grid/rtransit_{:s}.npz".format(metal[z]))
  p0       = grid ["p0"]/pc.bar
  ptransit = tgrid["pt"]/pc.bar
  # Delta-R factor:  D(radius) = H * fdr:
  fdr = np.log(p0/ptransit)[:,:,:,0]

  Z = np.copy(fdr)
  Z[~np.isfinite(Z)] = np.nan
  Z[Lambda<10] = np.nan
  vmin, vmax = np.nanmin(Z), np.nanmax(Z)
  print(vmin, vmax)

  for t in np.arange(nt):
    Z = np.copy(fdr[:,:,t])
    Z[~np.isfinite(Z)] = np.nan
    Z[Lambda[:,:,t]<10] = np.nan
    fmean[z,t]   = np.nanmean(Z)
    fmedian[z,t] = np.nanmedian(Z)
    fstd[z,t]    = np.nanstd(Z)
    print("{:.2f} {:.2f} {:4.0f}".format(np.nanmean(Z), np.nanmedian(Z),
                                         Teq[t]))

fig=plt.figure(-4, (7,4.5))
plt.clf()
plt.subplots_adjust(0.12, 0.15, 0.95, 0.95)
ax = plt.subplot(111)
for z in np.arange(nz):
  plt.plot(Teq, fmean[z], ".-", color=col[z], lw=lw, label=lab[z], ms=6)
  #plt.plot(Teq, fmedian[z], "--", color=col[z], lw=lw)
plt.xlabel(r"$\rm Temperature\ (K)$",    fontsize=fs)
plt.ylabel(r"$(R_{\rm T}-R_{\rm p})/H$", fontsize=fs)
plt.xlim(200, 3100)
l = plt.legend(loc="best")

renderer = fig.canvas.get_renderer()
shift = max([t.get_window_extent(renderer).width for t in l.get_texts()])
shift = 80.0
for t in l.get_texts():
  t.set_ha('right')
  t.set_position((shift,0))
plt.savefig("figs/scaleheight.ps")

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
thin = 10

fig=plt.figure(0, (7, 4.5))
plt.clf()
plt.subplots_adjust(0.12, 0.15, 0.95, 0.95)
ax = plt.subplot(111)
for z in np.arange(nz):
  grid  = np.load("run02_grid/MRT_{:s}.npz".format(metal[z]))
  tgrid = np.load("run02_grid/rtransit_{:s}.npz".format(metal[z]))
  p0 = grid ["p0"][:,:,:,0]
  Rt = tgrid["Rt"][:,:,:,0]

  Rcloud = np.zeros((nr,nm,nt)) # Radius at 1e-5 bar
  Rhill  = np.zeros((nr,nm,nt))

  # Get the altitude at 10**-5 bars (else, Hill radius):
  for itemp in np.arange(nt):
    spec, press, temp, q = pa.readatm("run01_atm/{:s}_{:04.0f}K.atm".
                                      format(metal[z], Teq[itemp]))
    press *= pc.bar                  # in CGS units
    mm     = pa.meanweight(q, spec)  # mean molecular mass
    for imass in np.arange(nm):
      for irad in np.arange(nr):
        if Rt[irad,imass,itemp] == 0.0:
          continue
        radius = pa.hydro_m(press, temp, mm, Mp[imass]*pc.mearth,
                     p0[irad,imass,itemp], Rp[irad]*pc.rearth)
        rhill  = smaxis[itemp] * (Mp[imass]*pc.mearth /
                                  (3*1.3*pc.msun))**(1.0/3.0)
        Rhill[irad,imass,itemp] = rhill
        radius[radius<0]     = np.nan
        radius[radius>rhill] = np.nan
        good = radius > 0
        # Radius at 1e-5 bar:
        if press[good][0] < 1e-5 * pc.bar:
          rinterp = si.interp1d(press[::-1], radius[::-1])
          Rcloud[irad,imass,itemp] = rinterp(1e-5 * pc.bar)
        # Else, Hill radius:
        elif radius[good][0] > Rt[irad,imass,itemp]:
          Rcloud[irad,imass,itemp] = rhill

  Z = 100.0 * (Rcloud-Rt) / Rt
  m = np.expand_dims(mask,2) & (np.isfinite(Z)) & (Lambda>=10) & (Rcloud!=Rhill)
  plt.loglog(Lambda[m][::thin], Z[m][::thin], ".", ms=6,
             color=col[z], label=lab[z])

plt.xlabel(r"$\Lambda$", fontsize=fs)
plt.ylabel(r"$\Delta R_{\rm c}/R\ \ (\%)$", fontsize=fs)
plt.xlim(9, 4000)
plt.ylim(0.05, 300)
plt.savefig("figs/Lambda_rcloud_metal.ps")
