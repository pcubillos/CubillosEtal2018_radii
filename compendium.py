# Clone code
# From the directory where this file is located, execute:
topdir=`pwd`
git clone --recursive https://github.com/pcubillos/pyratbay
cd $topdir/pyratbay
git checkout 8920b91
make

# TEA patch: with a text editor open:
$topdir/pyratbay/modules/TEA/tea/balance.py
# And change line 147 to:
            free_id.append(n + m)


# Make TEA atmospheric models:
cd $topdir/run01_atm
python $topdir/code/runatm.py

# Download H2O HITEMP data:
cd $topdir/inputs/opacity
wget --user=HITRAN --password=getdata -N -i wget_hitemp-H2O_0.3-1.2um.txt
unzip '*.zip'
rm -f *.zip

# Compile HITEMP H2O TLI file:
cd $topdir/run01_atm
$topdir/pyratbay/pbay.py -c tli_H2O.cfg
# Compile opacity grid:
$topdir/pyratbay/pbay.py -c opacity_H2O.cfg

# :::  OK!  :::


