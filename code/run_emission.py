#! /usr/bin/env python
import sys, os
import numpy as np

sys.path.append("../code")
import sma as sma

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants as pc
import pyratbay.tools     as pt


# Load Pyrat object:
pyrat   = pb.pyrat.init("MRTz_grid.cfg")
press   = pyrat.atm.press    # Atmospheric pressure array (barye)
nlayers = pyrat.atm.nlayers
pyrat.verb = 0  # Mute it

# Create new grid:
nf = pyrat.obs.nfilters  # Number of observing filter
nt = 28 # Number of temperature values
nr = 79 # Number of planetary radius values
nm = 45 # Number of planetary mass values
# The arrays:
Teq = np.linspace(300, 3000, nt)  # K
Rp  = np.linspace(1.0, 40.0, nr)  # Rearth
Mp  = np.logspace(np.log10(1.0), np.log10(630), 45)  # Mearth
Mp[-1]  = 600.0

# Reference pressure at the photosphere:
p0    = np.zeros((nr, nm, nt, nf), np.double)
# Has-this-run-been-computed flag:
flag  = np.zeros((nr, nm, nt), bool)
# Number of iterations:
niter = np.zeros((nr, nm, nt), int)
# Height of Hill radius over planetary radius in scale heights:
Hhill = np.zeros((nr, nm, nt), np.double)

# Max number of iterations per model:
nmax = 7
ph = np.zeros(nmax,   np.double) # Hill radius pressure
rp = np.zeros(nmax+1, np.double) # Reference pressure at each iteration

# Density:
dens = Mp*pc.mearth/(4./3*np.pi*np.expand_dims(Rp*pc.rearth,1)**3)
# Density flag:
dflag = np.ones((nr, nm), bool)
dflag[np.where(dens <  0.03)] = 0
dflag[np.where(dens > 30.00)] = 0

# Compute orbital semi-major axis from Teq and stellar mass:
smaxis = sma.sma(Teq, pyrat.phy.mstar/pc.msun) * pc.au

# Metallicity:
metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
nz = len(metal)

# Find p0, the pressure corresponding to Rp for each model in the grid:
for z in np.arange(nz):
  # Reset arrays for each metallicity
  p0   [:] = 0.0
  flag [:] = False
  niter[:] = 0
  Hhill[:] = 0.0
  for itemp in np.flipud(np.arange(nt)):
    spec, p, temp, q = pb.atmosphere.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                           format(metal[z], Teq[itemp]))
    # Set orbital semi-major axis:
    pyrat.phy.smaxis = smaxis[itemp]
    # Compute mean molecular mass:
    mm = np.sum(q*pyrat.mol.mass, axis=1)
    for imass in np.flipud(np.arange(nm)):
      pyrat.phy.mplanet = Mp[imass] * pc.mearth
      pyrat.phy.rhill   = pyrat.phy.smaxis * (pyrat.phy.mplanet /
                                             (3*pyrat.phy.mstar))**(1.0/3.0)
      for irad in np.arange(nr):
        # Invalid density:
        if not dflag[irad, imass]:
          continue
        print("\n::::::::::::::::::::::::::::::::::::\n"
                "r={:02d}/{:02d}  m={:03d}/{:03d}  t={:02d}/{:02d}  {:s}\n"
                "::::::::::::::::::::::::::::::::::::".
                format(irad, nr-1, imass, nm-1, itemp, nt-1, metal[z]))
        pyrat.phy.rplanet = Rp[irad] * pc.rearth
        if pyrat.phy.rplanet > pyrat.phy.rhill:
          print("Planet radius is greater than Hill radius.")
          break
        # Initial guess for the pressure reference level:
        rp[:] = 0.1*pc.bar  # Default
        if   irad != 0     and flag[irad-1,imass,itemp]:  # From prev Rp
          rp[:] = p0[irad-1,imass,itemp,0]
        elif imass != nm-1 and flag[irad,imass+1,itemp]:  # From prev Mp
          rp[:] = p0[irad,imass+1,itemp,0]
        # Iterate (nmax times) untill ref. pressure converges:
        for i in np.arange(nmax):
          # Compute spectra:
          pyrat.refpressure = rp[i]
          pyrat = pb.pyrat.run(pyrat, [temp, q])
          cf  = pt.cf(pyrat.od.depth, pyrat.atm.press, pyrat.od.B)
          bcf = pt.bandcf(cf, pyrat.obs.bandtrans, pyrat.spec.wn,
                          pyrat.obs.bandidx)
          rp[i+1] = 10**np.average(np.log10(press), weights=bcf[0])
          ph[i]   = press[pyrat.atm.rtop]
          print("Photosphere: {:.4f} bar  --  dH={:6.3f}  --  update={:6.2f}%".
                format(rp[i+1]/pc.bar, np.log(rp[i+1]/ph[i]),
                       100.0*np.abs(rp[i+1]-rp[i])/rp[i+1]))
          # Ref. pressure converged (improvement smaller than 2% in pressure):
          if i>=2 and np.abs(rp[i+1]-rp[i])/rp[i+1] < 0.02:
            break
        # Save results:
        for j in np.arange(pyrat.obs.nfilters):
          p0 [irad,imass,itemp,j] = 10**np.average(np.log10(press),
                                                   weights=bcf[j])
        Hhill[irad,imass,itemp] = np.log(rp[i+1]/ph[i])
        flag [irad,imass,itemp] = True
        niter[irad,imass,itemp] = i
        if np.log(rp[i+1]/ph[i]) < 2.0: # Break radii loop if < 2 scale heights
          break
      np.savez("MRT_{:s}.npz".format(metal[z]), Rp=Rp, Mp=Mp, Teq=Teq,
               flag=flag, p0=p0, niter=niter, Hhill=Hhill)
