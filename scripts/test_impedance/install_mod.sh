rm -rf neocortex
rm -rf x86_64

git clone --recursive ssh://bbpcode.epfl.ch/sim/models/neocortex

cp neocortex/mod/v5/Ca_HVA.mod neocortex/mod/v6

nrnivmodl neocortex/mod/v6
