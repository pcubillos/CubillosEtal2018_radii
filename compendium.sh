# Clone code
# From the directory where this file is located, execute:
topdir=`pwd`
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 0ddd082
make

cd $topdir
git clone https://github.com/pcubillos/repack
cd $topdir/repack
git checkout 906f885
make

# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Make TEA atmospheric models:
cd $topdir/run01_atm
python $topdir/code/run_atm.py


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Download H2O HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp-H2O_0.3-1.2um.txt
wget --user=HITRAN --password=getdata -N -i wget_hitemp_CO2_0.3-1.2um.txt
unzip '*.zip'
rm -f *.zip
# Download Exomol data:
wget -i wget_exomol_HCN.txt
bzip2 -d 1H-*.bz2  # Unzip only the HCN data
wget -i wget_exomol_NH3_0.3-1.2um.txt
wget -i wget_exomol_CH4_0.3-1.2um.txt
# Download Li2015 CO data:
wget -i wget_Li2015_CO.txt
tar -xvzf apjs504015_data.tar.gz
rm -f apjs504015_data.tar.gz ReadMe Table_S1.txt Table_S2.txt \
      Table_S3.txt Table_S6.par

# Generate partition-function files for NH3, HCN, and CH4:
cd $topdir/run01_atm
python $topdir/code/pf_tips_NH3.py
python $topdir/pyratbay/scripts/PFformat_Exomol.py  \
       $topdir/inputs/opacity/1H-12C-14N__Harris.pf \
       $topdir/inputs/opacity/1H-13C-14N__Larner.pf
python $topdir/pyratbay/scripts/PFformat_Exomol.py \
       $topdir/inputs/opacity/12C-1H4__YT10to10.pf

# Repack large line-transition datasets:
cd $topdir/run01_atm
$topdir/repack/repack.py repack_NH3.cfg
$topdir/repack/repack.py repack_CH4.cfg

# Compile TLI files:
cd $topdir/run01_atm
python $topdir/pyratbay/pbay.py -c tli_H2O.cfg
python $topdir/pyratbay/pbay.py -c tli_CO2.cfg
python $topdir/pyratbay/pbay.py -c tli_CO.cfg
python $topdir/pyratbay/pbay.py -c tli_CH4.cfg
python $topdir/pyratbay/pbay.py -c tli_HCN.cfg
python $topdir/pyratbay/pbay.py -c tli_NH3.cfg


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compile opacity grid:
cd $topdir/run01_atm
python $topdir/pyratbay/pbay.py -c opacity_H2O-CH4-NH3.cfg


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Compute optical photospheric pressure over the MRTz grid:
cd $topdir/run02_grid
python $topdir/code/run_emission.py

# Compute transmission radii for grid:
cd $topdir/run02_grid
python $topdir/code/run_transmission.py


# ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
# Figure 1:
cd $topdir
python $topdir/fig_density.py

# Figure 2:
cd $topdir/run01_atm
python $topdir/fig_opacity.py

# Figure 3:
cd $topdir/run01_atm
python $topdir/fig_spectrum.py

# Figure 4:
cd $topdir
python $topdir/fig_slices.py

# Figure 5:
cd $topdir
python $topdir/fig_Lambda.py

# Figure 6:
cd $topdir/run01_atm
python $topdir/fig_mean_opacity.py

# Figures 7 and 8:
cd $topdir
python $topdir/fig_clouds.py

# Figure 9:
cd $topdir
python $topdir/fig_metal.py
