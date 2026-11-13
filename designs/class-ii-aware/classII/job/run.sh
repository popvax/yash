#!/bin/sh
# Class II binding for Yash's typed alleles. Needs netMHCIIpan on PATH
# (academic licence, obtained by the user, not bundled here).
set -eu
ALLELES=$(grep -v '^#' alleles.txt | paste -sd, -)

# Per-region: does each 25-aa antigen block contain a CD4 epitope at all?
netMHCIIpan -f regions.fasta -inptype 0 -a "$ALLELES" -length 15 \
    -BA -xls -xlsfile regions.classII.xls > regions.classII.txt

# Whole cassette: also covers epitopes formed across antigen junctions,
# which per-region scoring cannot see.
netMHCIIpan -f cassette.fasta -inptype 0 -a "$ALLELES" -length 15 \
    -BA -xls -xlsfile cassette.classII.xls > cassette.classII.txt
