#!/bin/bash

# Scripts to run NICER code tests

base="/Users/sharon/Dropbox/code"
exe_dir="$base/nicer"
#pwd
make spot

out_dir="$exe_dir/LINE"

normalizeFlux=false  # false = no, true = yes; still need to add or remove -N flag in command line (-N means yes)

gray=0    # 0 = isotropic emission, 1 = graybody

# integers
numtheta=1 # number of theta bins; for a small spot, only need one
NS_model=3 # 1 (oblate) or 3 (spherical); don't use 2
numbins=128

# doubles -- using "e" notation
spin=1 # in Hz
mass=1400000e-6 # in Msun
radius=12000000e-6 # in km
inclination=90 # in degrees
emission=90 # in degrees
rho=1e-2  # in radians
phaseshift=0.0 # this is used when comparing with data
temp=35e-2 # in keV (Surface temperature)
#distance=3.08567758e20 # 10 kpc in meters
distance=6.1713552e18 # 200 pc in meters


## MAKING THE DATA FILE
if test ! -d "$out_dir"
   	then mkdir -p "$out_dir"
fi

# TEST 1 = 1 Hz, tiny spot; narrow line
numtheta=1 # number of theta bins; for a small spot, only need one
out_file="$out_dir/line1.txt"
spin=1
# NICER line emission
spectraltype=1
numbands=200
elow=0.99998
ehi=1.00002
e0=0.809
deltae=1e-5
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$gray" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -x "$e0" -X "$deltae"

# TEST 2 = 1 Hz, big spot; narrow line
numtheta=16
numphi=16
rho=1
out_file="$out_dir/line2.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$gray" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -x "$e0" -X "$deltae"

times
# TEST 3 = 400 Hz, tiny spot; wide line
numtheta=1 # number of theta bins; for a small spot, only need one
spin=400
rho=1e-2
# NICER line emission
spectraltype=1
numbands=300
elow=0.995
ehi=1.005
e0=0.7
deltae=1e-3
out_file="$out_dir/line3.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$gray" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -x "$e0" -X "$deltae"
times

# TEST 4 = 400 Hz, big spot; wide line
numtheta=16 # number of theta bins; for a small spot, only need one
numphi=16
spin=400
rho=1
out_file="$out_dir/line4.txt"
## RUNNING THE CODE
./spot -m "$mass" -r "$radius" -f "$spin" -i "$inclination" -e "$emission" -l "$phaseshift" -n "$numbins" -q "$NS_model" -o "$out_file" -p "$rho" -T "$temp" -D "$distance" -t "$numtheta" -g "$gray" -s "$spectraltype" -S "$numbands" -u "$elow" -U "$ehi" -x "$e0" -X "$deltae"

