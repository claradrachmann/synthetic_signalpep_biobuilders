#!/usr/bin/env python3

import sys
import os
import subprocess
from pathlib import Path
import numpy as np

n = len(sys.argv)
if n < 4:
    sys.exit("ERROR: must supply, in that order, template genome file, desired locus ORF id, number of noisy samples from each noise level, and number of bases upstream to consider")

template_genome = sys.argv[1]
locus = sys.argv[2]
n_samples = sys.argv[3]
bases_upstream = sys.argv[4]

template_genome = Path(template_genome).resolve()

if not template_genome.exists():
    sys.exit("Template genome not found")

template_genome = template_genome.as_posix()

try:
    Path(locus).mkdir()
except FileExistsError:
    sys.exit("ERROR: Output directory already exists. Not overwriting.")
os.chdir(locus)

print("Picking out sequences")
subprocess.run("python ../complete.py "+template_genome+" "+locus+" "+bases_upstream+" 0",shell=True)

print("Generating homology model")
subprocess.run("python ../em5.py "+locus,shell=True)

with open(locus+"-hom-upseqs.final.hmm") as f:
    np.savetxt(
        locus+"-hom-upseqs.final.csv",
        np.array(list(map(lambda l: list(map(float,l.split()[1:5])),list(f)[19:-1:3]))),
        fmt="%.5f",
        delimiter=",",
        header="A,C,G,T",
        comments="")

print("Generating synthetic promoter samples")
subprocess.run("julia ../hmm-match-state-emission.jl "+locus+" "+n_samples,shell=True)

# Genomes/Aspni_DSM_1/Aspni_DSM_1_GeneCatalog_proteins_20130526.aa.fasta
