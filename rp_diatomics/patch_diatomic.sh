#!/bin/bash

cat diatomic_atomtype.atp >> atomtypes.atp
cat diatomic_ffnonbonded.itp >> ffnonbonded.itp
rm ffbonded.itp
cat ffbonded_top.itp diatomic_ffbonded.itp ffbonded_bottom.itp >> ffbonded.itp
