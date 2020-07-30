# Synthetic promoter generation from conserved regions upstream of homologous genes

## Dependencies
* Linux or Mac OS X
* [MUSCLE](http://www.drive5.com/muscle/)
* [HMMER](http://hmmer.org/)
* Python 3 with the following packages:
  * Biopython
  * numpy
* Julia 1 with the following packages:
  * JuMP
  * GLPK
  * CSV
  * StatsBase

## Invocation
`python prohmmoter.py <TemplateGenome> <Locus> <NNoisy> <NUpstream>`

Where:

* `TemplateGenome` is the `GeneCatalog_proteins` file from the genome
you want to use as the reference
when selecting the genes whose promoters you want.
* `Locus` is the identifier of the gene whose promoter you want
* `NNoisy` is the number of noisy samples you want from each defined noise level
(1, 1.5, and 2 times as high entropy as picking the promoter from one of the genomes)
* `NUpstream` is the number of bases upstream of the start codon you wish to use
as the basis for what you consider the promoter

The script will then create a folder named `<Locus>`
in the current working directory,
where the scripts used for each sub-task will output its files.

The relevant final output files are:

* `<Locus>-promoter-samples.nt.fasta` containing the generated synthetic promoters
* `<Locus>-hom-upseqs.final.hmm` containing the HMM of the specific promoter
* `<Locus>-hom-upseqs.nt.fasta` containing the sequences upstream of each identified ortholog

## Files required
For each genome in a folder `<GenomeName>`:

* `<GenomeName>_AssemblyScaffolds_Repeatmasked.fasta`
* `<GenomeName>_GeneCatalog_genes_<optional? date>.gff`

For the template genome, additionally:

* `<GenomeName>_GeneCatalog_proteins_<optional? date>.aa.fasta`
