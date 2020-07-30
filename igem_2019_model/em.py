#!/usr/bin/env python3
import sys
import glob
import subprocess
import os
import shutil
import itertools
from Bio import SeqIO

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

locus = sys.argv[1]

os.mkdir("em")
subprocess.run("cp "+locus+"-hom-upseqs.clu em/"+locus+"-hom-upseqs.0.clu",shell=True)

with cd("em"):
    initialize_command = "hmmbuild --enone --fragthresh 0 "+locus+"-hom-upseqs.0.hmm "+locus+"-hom-upseqs.0.clu"
    subprocess.run(initialize_command,shell=True)
    subprocess.run("hmmemit -o "+locus+"-hom-upseqs.consensus.0.nt.fasta -C "+locus+"-hom-upseqs.0.hmm",shell=True)
    i = 0
    while True:
        i += 1
        command_line = "hmmalign --dna --trim -o "+locus+"-hom-upseqs."+str(i)+".sto "+locus+"-hom-upseqs."+str(i-1)+".hmm "+locus+"-hom-upseqs.0.clu && hmmbuild --enone "+locus+"-hom-upseqs."+str(i)+".hmm "+locus+"-hom-upseqs."+str(i)+".sto"
        subprocess.run(command_line,shell=True)
        subprocess.run("hmmemit -o "+locus+"-hom-upseqs.consensus."+str(i)+".nt.fasta -C "+locus+"-hom-upseqs."+str(i)+".hmm",shell=True)
        s_before = next(SeqIO.parse(locus+"-hom-upseqs.consensus."+str(i-1)+".nt.fasta","fasta")).seq
        s_now = next(SeqIO.parse(locus+"-hom-upseqs.consensus."+str(i)+".nt.fasta","fasta")).seq
        if s_before == s_now:
            break
    shutil.copy(locus+"-hom-upseqs."+str(i)+".hmm","../"+locus+"-hom-upseqs.final.hmm")
