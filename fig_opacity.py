#! /usr/bin/env python
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter1d as gaussf
plt.ioff()

sys.path.append("../pyratbay")
import pyratbay as pb
import pyratbay.constants  as pc
import pyratbay.tools      as pt
import pyratbay.atmosphere as pa


ntemp = 28
metal = ["0.1xsolar", "solar", "10xsolar", "100xsolar"]
Teq = np.linspace(300, 3000, ntemp)

sigma = 2.0
lw    = 1.5
yticks = np.logspace(-21, 9, 11)


# Initialize pyrat object:
pyrat = pb.pyrat.init("spectrum_lbl.cfg")
# Wavelength:
wl = 1.0/(pyrat.spec.wn*pc.um)
# Atmospheric pressure profile (bar):
press = pyrat.atm.press/pc.bar
#  60  0.001 bar
#  72  0.01  bar
#  84  0.1   bar
#  95  1.0   bar
# 107 10.0   bar
sp = ["0.001", "0.1", "10"]
ipress = [60, 84, 107]
nheights = len(ipress)

# The samples:
t = 0, 9, 27
z = 1, 1, 1
l = 60, 107
sp = ["0.001", "10"]

exc, excl = [], []
for i in np.arange(3):
  spec, p, temp, q = pb.atmosphere.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                       format(metal[z[i]], Teq[t[i]]))
  pyrat.atm.temp = temp
  pyrat.atm.q    = q
  for j in np.arange(pyrat.mol.nmol):
    pyrat.atm.d[:,j] = pa.IGLdensity(pyrat.atm.q[:,j], pyrat.mol.mass[j],
                                   pyrat.atm.press,  pyrat.atm.temp)
  exc.append([])
  excl.append([])
  for j in np.arange(2):
    ext = pyrat.get_ec(l[j])
    exc[i].append(ext[0])
    excl[i].append(ext[1])

ecl = excl[0][0]
imol = np.where(np.in1d(ecl, ['H2O', 'CH4', 'CO', 'CO2', 'HCN', 'NH3']))[0]


dash = [(),          # Solid
        (5,2),       # Dash
        (7,2,2,2), ] # Dash-dot
spec, press, T, q  = pa.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                        format(metal[0], Teq[0]))
# Species to show:
show = ["H2O", "CO", "CO2", "HCN", "NH3", "CH4", "H2", "He", "H", "Na", "K"]
ishow = np.where(np.in1d(spec, show))[0]
# Species with solid dash:
mol = ["H2O", "CO", "CO2", "HCN", "NH3", "CH4"]
ls = [dash[1],]*len(spec)
for i in np.where(np.in1d(spec, mol))[0]:
 ls[i] = dash[0]
for i in np.where(np.in1d(spec, ['Na','K']))[0]:
 ls[i] = dash[2]

# Latex species names:
lmol = [r'$\rm H_2$', r'$\rm He$', r'$\rm Na$', r'$\rm K$', r'$\rm H$',
        r'$\rm H_2O$', r'$\rm CH_4$', r'$\rm CO$', r'$\rm CO_2$',
        r'$\rm HCN$', r'$\rm NH_3$',
        r'$\rm N_2$', r'$\rm C_2H_2$', r'$\rm C_2H_4$']

# Colors:
ac = ['orange', 'k', 'limegreen', 'g', 'teal', 'b', 'peru', 'r',
      'gold', 'm', 'dodgerblue', '', '', '', '', '']
oc = ['blue', 'peru', 'r', 'gold', 'm', 'dodgerblue']

ap = 1/(np.log10(np.amin(press))-np.log10(np.amax(press)))
bp = -ap*np.log10(np.amax(press))
yy = [[np.log10(1.5 *press[l[0]])*ap +bp,
       np.log10(0.75*press[l[0]])*ap +bp],
      [np.log10(1.5 *press[l[1]])*ap +bp,
       np.log10(0.75*press[l[1]])*ap +bp]]

shade = "0.85"
fs   = 12
thin =  6
matplotlib.rcParams.update({'xtick.labelsize':fs-1})
matplotlib.rcParams.update({'ytick.labelsize':fs-1})
al, aw = 0.085, 0.23
ah = 0.284
bot = 0.06 + np.linspace(2,0,3)*(ah+0.035)
ol, ow = 0.41, 0.562
oh = ah/2.2
dh = ah - 2*oh

plt.figure(0, (8.5, 10.5))
plt.clf()
for i in np.arange(3):
  for j in np.arange(2):
    ax0 = plt.axes([al, bot[i], ow+ol-al, ah])
    ax0.set_xlim(al, ow+ol)
    ax0.set_ylim(0,1)
    ax0.set_axis_off()
    verts = [(al+aw, yy[j][0])] + [(al+aw, yy[j][1])] + \
            [(ol, (dh+(2-j)*oh)/(dh+2*oh))] + [(ol, (dh+(1-j)*oh)/(dh+2*oh))]
    poly = matplotlib.patches.Polygon(verts, facecolor=shade, edgecolor='none')
    ax0.add_patch(poly)
  # Abundance panels:
  lw = 1.5
  zz, tt = z[i], t[i]
  ax = plt.axes([al, bot[i], aw, ah])
  spec, p, temp, q = pa.readatm("../run01_atm/{:s}_{:04.0f}K.atm".
                                           format(metal[zz], Teq[tt]))
  plt.text(0.05, 0.93, r"$T\, =\, {:4.0f}\ \rm K$".format(temp[tt]), ha="left",
             transform=ax.transAxes, fontsize=fs)
  #plt.text(0.05, 0.84, "[M/H]={:d}".format(zz-1), ha="left",
  #           transform=ax.transAxes, fontsize=fs)
  plt.axhspan(0.75*press[l[0]], 1.5*press[l[0]], color="0.8", zorder=0)
  plt.axhspan(0.75*press[l[1]], 1.5*press[l[1]], color="0.8", zorder=0)
  for j in ishow:
    plt.loglog(q[:,j], press, label=lmol[j], lw=lw, dashes=ls[j], color=ac[j])
    ax.set_yticks([1e2, 1, 1e-2, 1e-4, 1e-6, 1e-8])
    ax.set_xticks([1, 1e-4, 1e-8, 1e-12, 1e-16])
  if i==0:
    plt.legend(loc='lower left', fontsize=fs-2, labelspacing=0.2)
  if i==2:
    plt.xlabel(r"$\rm Mole\ fraction$", fontsize=fs)
  plt.xlim(1e-15, 3.0)
  plt.ylim(np.amax(press), np.amin(press))
  plt.ylabel(r"$\rm Pressure\ (bar)$", fontsize=fs)
  # Opacity panels:
  lw = 1.25
  for j in np.arange(2):
    ec, ecl = exc[i][j], excl[i][j]
    alkali = np.sum(ec[np.in1d(ecl, ['Na', 'K'])],        axis=0)
    cia    = np.sum(ec[np.in1d(ecl, ['H2-H2', 'H2-He'])], axis=0)
    ray = np.sum(ec[np.in1d(ecl,['lecavelier','dalgarno_He','dalgarno_H'])],0)
    ax = plt.axes([ol, bot[i]+oh*(1-j)+dh, ow, oh])
    for k in imol:
      plt.plot(wl[::thin], gaussf(ec[k], sigma)[::thin], color=oc[k], lw=lw)
    plt.plot(wl[::thin],   gaussf(alkali,sigma)[::thin], color="limegreen",
             lw=lw, label=r'$\rm Alkali$')
    plt.plot(wl[::thin],   gaussf(ray,   sigma)[::thin], color="0.5", lw=lw,
             label=r'$\rm Rayleigh$', dashes=(8,2))
    plt.plot(wl[::thin],   gaussf(cia,   sigma)[::thin], color="k",   lw=lw,
             label=r'$\rm CIA$', dashes=(8,2))
    plt.yscale("log")
    ax.set_xticks([0.4, 0.6, 0.8, 1.0, 1.2])
    plt.xlim(np.amin(wl), np.amax(wl))
    ax.set_yticks(np.logspace(-30, 9, 14))
    ymax = 2*np.amax(ec)
    plt.ylim(ymax*1e-18, ymax)
    plt.text(0.97, 0.87, r"${:s}\ \rm bar$".format(sp[j]), ha="right",
             transform=ax.transAxes, fontsize=fs)
    if i == 0 and j==0:
      plt.legend(loc="upper left", fontsize=fs-2,  labelspacing=0.15)
    if j == 1:
      plt.xlabel(r"$\rm Wavelength\ (um)$", fontsize=fs)
    else:
      plt.ylabel(r"${\rm Extinction\ coefficient\ (cm}^{-1})$" + 25*" ",
                 fontsize=fs)
      ax.set_xticklabels([])
plt.savefig("../figs/opacity.ps")
