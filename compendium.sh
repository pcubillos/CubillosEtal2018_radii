# Clone code
# From the directory where this file is located, execute:
topdir=`pwd`
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 80fd85d
make


# Make TEA atmospheric models:
cd $topdir/run01_atm
python $topdir/code/run_atm.py

# Download H2O HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp-H2O_0.3-1.2um.txt
unzip '*.zip'
rm -f *.zip

# Compile HITEMP H2O TLI file:
cd $topdir/run01_atm
python $topdir/pyratbay/pbay.py -c tli_H2O.cfg
# Compile opacity grid:
python $topdir/pyratbay/pbay.py -c opacity_H2O.cfg


# Compute optical photospheric pressure over the MRTz grid:
cd $topdir/run02_grid
python $topdir/code/run_emission.py

# Compute transmission radii for grid:
cd $topdir/run02_grid
python $topdir/code/run_transmission.py


# Figure 1:
cd $topdir
python $topdir/fig_density.py

# Figure 2:
cd $topdir/run01_atm
python $topdir/fig_spectrum.py

# Figure 3:
cd $topdir
python $topdir/fig_slices.py

# Figure 4:
cd $topdir
python $topdir/fig_Lambda.py

# Figures 5 and 6:
cd $topdir
python $topdir/fig_clouds.py

# Figure 7:
cd $topdir
python $topdir/fig_metal.py
