#! /usr/bin/env python
import sys
import numpy as np
import scipy.interpolate as si
import scipy.constants   as sc
import matplotlib
import matplotlib.pyplot as plt
plt.ioff()

sys.path.append("./code")
import colormaps as cm
import sma       as sma

sys.path.append("./pyratbay")
import pyratbay.constants  as pc
import pyratbay.atmosphere as pa


# The grid:
grid  = np.load("run02_grid/MRT_solarb.npz")
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
dr = Rp[1] - Rp[0]
Rplot = np.concatenate((Rp, [Rp[-1]+dr])) - dr/2.0
logM = np.log(Mp)
dlm = logM[1] - logM[0]
Mplot = np.exp(np.concatenate((logM[0:1]-dlm, logM)) + dlm/2)

dens = Mp*pc.mearth / (4./3*np.pi*np.expand_dims(Rp*pc.rearth,1)**3)
mask = (dens>0.03) & (dens <=30.0)

palette = cm.inferno_r
palette.set_under(color=(1,1,1,0))
palette.set_bad(color='w')

lw = 1.25
fs = 12


for z in np.arange(nz):
  grid  = np.load("run02_grid/MRT_{:s}.npz".format(metal[z]))
  tgrid = np.load("run02_grid/rtransit_{:s}.npz".format(metal[z]))
  p0 = grid ["p0"][:,:,:,0]
  Rt = tgrid["Rt"][:,:,:,0]

  Rcloud = np.zeros((nr,nm,nt)) # Radius at 1e-5 bar
  Rhill  = np.zeros((nr,nm,nt)) # Radius at 1e-5 bar

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

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  Z = 100*(Rcloud-Rt) / Rt
  Z[~np.isfinite(Z)] = np.nan
  Z[p0==0]     = np.nan
  Z[Lambda<10] = np.nan
  print(np.nanmin(Z), np.nanmax(Z))

  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  # The plot:
  matplotlib.rcParams.update({'ytick.labelsize':10})
  matplotlib.rcParams.update({'xtick.labelsize':10})
  plt.figure(11, (8.5, 4))
  plt.clf()
  plt.subplots_adjust(0.07, 0.13, 0.88, 0.97, hspace=0.05, wspace=0.05)
  j=0
  for i in [0, 4, 8, 12, 16, 20, 24, 27]:
    Z = 100.0 * (Rcloud[:,:,i]-Rt[:,:,i]) / Rt[:,:,i]
    Zmin = 0.0
    Zmax = 200.0
    ax = plt.subplot(2, 4, j+1)
    plt.contour(Mp, Rp, mask, 1, colors="0.7", linewidths=0.75, zorder=0)
    plt.pcolor(Mplot, Rplot, Z, edgecolors="face",
               vmin=Zmin, vmax=Zmax, cmap=palette)
    rlambda10 = 0.1*sc.G*(Mp*pc.mearth)*mh/(10.0*sc.k*Teq[i])/pc.rearth
    plt.plot(Mp, rlambda10, "limegreen", lw=lw)
    ax.set_xscale('log')
    if j%4 == 0:
      plt.ylabel(r"${\rm Radius}\ (R_{\oplus})$", fontsize=fs)
    else:
      ax.set_yticklabels([])
    if j >= 4:
      plt.xlabel(r"${\rm Mass}\ (M_{\oplus})$",   fontsize=fs)
    else:
      ax.set_xticklabels([])
    plt.xlim(Mplot[0], Mplot[-1])
    plt.ylim(Rplot[0], Rplot[-1])
    plt.text(1.1, 36, r"$T=\ {:4.0f}\ {{\rm K}}$".format(Teq[i]),
             fontsize=fs)
    j += 1
  cax = plt.axes([0.89, 0.13, 0.016, 0.84])
  cax = plt.colorbar(cax=cax)
  cax.set_label(r"$\Delta R_{\rm c}/R\ \ (\%)$", fontsize=fs)

  cbticks = [0, 50, 100, 150, 200]
  cax.set_ticks(cbticks)
  cax.set_ticklabels([r"$0$", r"$50$", r"$100$", r"$150$", r"$>200$"])
  cax.solids.set_edgecolor("face")
  plt.savefig("./figs/cloud_transit_{:s}.ps".format(metal[z]))


  # ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
  thin = 4
  fs2 = 18
  matplotlib.rcParams.update({'ytick.labelsize':14})
  matplotlib.rcParams.update({'xtick.labelsize':14})

  fig=plt.figure(0, (7, 4.5))
  plt.clf()
  plt.subplots_adjust(0.12, 0.15, 0.84, 0.95)
  ax = plt.subplot(111)
  for i in np.arange(nt):
    rc  = Rcloud[:,:,i]
    rt  = Rt    [:,:,i]
    rh  = Rhill [:,:,i]
    lam = Lambda[:,:,i]
    Z = 100.0 * (rc-rt) / rt
    m1 = mask & (np.isfinite(Z)) & (rc!=rh)
    m2 = mask & (np.isfinite(Z)) & (lam >= 10) & (rc!=rh)

    l1 = lam[m1]
    l2 = lam[m2]
    dR1 = Z[m1]
    dR2 = Z[m2]
    plt.plot(l2[::thin], dR2[::thin], ".", ms=4,
             color=cm.plasma(int(254.0*i/(nt-1))))
  plt.xscale("log")
  plt.yscale("log")
  plt.xlabel(r"$\Lambda$", fontsize=fs2)
  plt.ylabel(r"$\Delta R_{\rm c}/R\ \ (\%)$", fontsize=fs2)
  plt.xlim(9, 4000)
  plt.ylim(0.1, 300)

  if z == 1:
    a = np.logspace(0.1, 4.0, 100)
    c= -1
    b = 400*a**c
    plt.plot(a, b, color="limegreen", lw=2, linestyle="--")
  # Colorbar
  bb = ax.get_position()
  ax2 = fig.add_axes([0.86, bb.y0, 0.025, bb.height])
  norm = matplotlib.colors.Normalize(vmin=0, vmax=255)
  cb1 = matplotlib.colorbar.ColorbarBase(ax2, cmap=cm.plasma, norm=norm,
                                  orientation='vertical')
  cb1.set_label(r"${\rm Temperature\ (K)}$", fontsize=fs2)
  tl = np.linspace(300, 3000, 10, dtype=int)
  ticks = (tl-300) * 255.0 / 2700.0
  cb1.set_ticks(ticks, True)
  cb1.set_ticklabels(tl, True)
  cb1.solids.set_edgecolor("face")
  plt.savefig("figs/Lambda_rcloud_{:s}.ps".format(metal[z]))
