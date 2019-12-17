#set -e
#module purge 

#module load python-dev 
#module load py-bglibpy
#module load brion
#module load neurodamus-neocortex

. ~/imped/bin/activate

export HOC_LIBRARY_PATH=/home/arnaudon/codes/bbp/lib/hoclib

export OMP_NUM_THREADS=1

python3 compute_impedance.py
