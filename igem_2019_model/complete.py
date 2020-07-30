#!/usr/bin/env python3
import sys
import glob
import subprocess
import pathlib

from Bio import SearchIO
from Bio import SeqIO

template_genome = sys.argv[1]
locus = sys.argv[2]
bases_upstream = int(sys.argv[3])
bases_downstream = int(sys.argv[4])

genomes_dir = pathlib.Path(template_genome).parent.parent.as_posix()

proteins = SeqIO.index(template_genome,"fasta",key_function = lambda s: s.split("|")[3][0:-4])
protein = proteins[locus]
print(str(len(protein)),"AA long protein")
SeqIO.write(protein, locus+".aa.fasta", "fasta")

# TODO Look only for matches with a single domain
command_line = "cat "+genomes_dir+"/*/*.aa.fasta | phmmer --noali -o "+locus+"-matches.txt "+locus+".aa.fasta -"
subprocess.run(command_line,shell=True)

res = SearchIO.read(locus+"-matches.txt","hmmer3-text")

match_list = list(map(lambda h: h.id.split("|"), filter(lambda h: h.is_included, res.hits)))


matches = {g: [] for g in set(map(lambda m: m[1], match_list))}
for match in match_list:
    matches[match[1]].append(match[3])

print(str(sum(map(len,matches.values()))),"putative homologs found")

for genome, ms in matches.items():
  if len(ms) > 0:
    matches[genome] = [ms[0]]

print(str(sum(map(len,matches.values()))),"putative orthologs found")

hits = {g: [] for g in matches.keys()}
for genome_name, genes in matches.items():
  # There should only be 1 such file
  files = glob.glob(genomes_dir+"/"+genome_name+"/"+genome_name+"_GeneCatalog_genes_*.gff")
  with open(files[0],"r") as f:
    for line in f:
      if "start_codon" in line:
        for gene in genes:
          if "name \""+gene+"\"" in line:
            hit_line = line.split("\t")
            hits[genome_name].append((gene,hit_line[0],hit_line[6],int(hit_line[3]),int(hit_line[4])))

print(str(sum(map(len,hits.values()))),"hits recovered from",len(hits),"genomes")

def find_upseq(loc, genome, n_up=1000, n_down=200):
    contig = genome[loc[1]]
    print(loc)
    if loc[2] == "+":
        upseq = contig.seq[loc[3]-1-n_up:loc[4]+n_down]
    else:
        upseq = contig.seq[loc[3]-1-n_down:loc[4]+n_up].reverse_complement()
    return SeqIO.SeqRecord(upseq.upper(), loc[0], description=str(n_up)+"_bases_upstream")

with open(locus+"-hom-upseqs.nt.fasta","w") as f:
    for genome_name, locs in hits.items():
        print(genome_name)
        genome = SeqIO.index(genomes_dir+"/"+genome_name+"/"+genome_name+"_AssemblyScaffolds_Repeatmasked.fasta","fasta")
        for loc in locs:
            sr = find_upseq(loc, genome, bases_upstream, bases_downstream)
            sr.id = genome_name+"|"+sr.id
            if len(sr.seq) == bases_upstream + 3 + bases_downstream and 'N' not in sr.seq:
                SeqIO.write(sr, f, "fasta")

hom_upseqs = SeqIO.index(locus+"-hom-upseqs.nt.fasta","fasta")
print(str(len(hom_upseqs)),"proper upstream regions found")

command_line_align = "muscle -clwstrict -in "+locus+"-hom-upseqs.nt.fasta -out "+locus+"-hom-upseqs.clu"
subprocess.run(command_line_align,shell=True)
