#!/usr/bin/env python3
# coding=utf-8
import sys
from datetime import datetime
import argparse
import pdb
import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import subprocess
from io import StringIO
from random import choices
import seaborn as sns
sns.set(style="darkgrid")


def get_parser():
    # Build commandline parser
    parser = argparse.ArgumentParser(
        description="Take in modified HMM-output from signalP-3.0, extract N, H and C regions, define amino acid distributions for"
                    "each region and make a fasta file with proposed synthetic signal peptides")
    # Arguments
    parser.add_argument("-in", "--infiles", type=str, nargs='+', dest="infiles", help="Input files", required=True)
    parser.add_argument("-out", "--outfile", type=str, dest="outfile",
                        help="Name of output file. Default: {infile}_synthetic_signalpeps.fa")
    parser.add_argument("-n_pep", "--n_peptides", type=int, dest="n_pep", default=10000,
                        help="Number of synthetic signal peptides to compute. Default: 1000")

    return parser


def get_args():
    parser = get_parser()
    args = parser.parse_args()
    if not args.outfile:
        args.outfile = os.path.splitext(args.infiles[0])[0] + '_synthetic_signalpeps.fa'

    return args


def process_signalp3(filename, aa_count_dict, region_lengths, critical_creg_positions):
    """
    :param filename: name of a single signalp3 modified output file
    :param aa_count_dict: dict with key:region, val:(dict with key:amino_acid, val:count)
    :param region_lengths: dict with key: region, val:list of observed length
    :param critical_creg_positions: dict of observed amino acids on positions -3, -2 and -1
    :return: dict of updated region counts, updated lists of observed region lengths and critical_creg_positions
    """

    try:
        infile = open(filename, 'r')
    except IOError as err:
        print("Cannot open file:", str(err));
        sys.exit(1)

    all_lines = infile.readlines()

    # we have cut off the data so we know that each sequence takes up 42 lines
    for i in range(0,len(all_lines),42):
        header = all_lines[i]
        seq_data = pd.read_csv(StringIO(''.join(all_lines[i+1:i+42])), sep='\t')
        try:
            # cut off data at most probable cleavage site - we do not want to include any positions after this
            signalpep_data = seq_data.iloc[:seq_data['C'].idxmax()].drop(columns=['C','S','pos'])
        except:
            pdb.set_trace()
        # extract region of highest likelihood
        if not signalpep_data.empty:
            # get amino acids on critical C region positions (cleavage site)
            critical_creg_aa = signalpep_data.iloc[-3:, :]['aa'].to_list()
            for pos in critical_creg_positions:
                critical_creg_positions[pos].append(critical_creg_aa[pos])

            signalpep_data['region'] = signalpep_data[['n-reg','h-reg','c-reg']].idxmax(axis=1).str.split('-').str[0]
            # get lengths of regions 
            reg_count = signalpep_data['region'].value_counts()
            for reg in reg_count.keys():
                region_lengths[reg].append(reg_count[reg])
            # add number of each amino acid to count dicts for each region
            region_dfs = {k: v for k, v in signalpep_data.groupby(by='region')}
            for reg, df in region_dfs.items():
                val_counts = df['aa'].value_counts()
                try:
                    for aa in val_counts.keys():
                        aa_count_dict[reg][aa.upper()] += val_counts[aa]
                except KeyError as err:
                    print(f'Obs! Non-recognized amino acid: {aa}')

    return aa_count_dict, region_lengths, critical_creg_positions


def synthetic_signalpep(region_aa_count_dict, region_lengths, optimal_creg_cleavage, n_peptides=1000):
    """
    Create n_peptides synthetic signal peptides by sampling the optimal length
    from the region_aa_count_dict distribution for each region
    :param region_aa_count_dict: key:region, value:(dict with key:amino_acid, val:count)
    :param region_lengths: dict: key: region: value: optimal length
    :param optimal_creg_cleavage: str: amino acids which should be the end of C region
    :param n_peptides: number of synthetic signal peptides to compute
    :return: signal_peptides: list of computed synthetic signal peptides
    """
    signal_peptides = list()
    for i in range(n_peptides):
        signalpep = ''
        for reg in ('n','h','c'):
            aa_population = [k for k, v in sorted(region_aa_count_dict[reg].items())]
            aa_weights = [v for k, v in sorted(region_aa_count_dict[reg].items())]
            if reg != 'c':
                # sample from the amino acid distribution of the given region and add to the signal peptide sequence
                signalpep += ''.join(choices(aa_population,aa_weights,k=region_lengths[reg]))
            else:
                # if C region, sequence should end on optimal_creg_cleavage
                signalpep += ''.join(choices(aa_population, aa_weights, k=region_lengths[reg]-len(optimal_creg_cleavage)))
                signalpep += optimal_creg_cleavage
        signal_peptides.append(signalpep)

    return signal_peptides



def main(args):
    # initialize
    amino_acids = 'ARNDCQEGHILKMFPSTWYV'
    region_aa_count_dict = {'n':{aa:0 for aa in amino_acids}, 'h':{aa:0 for aa in amino_acids}, 'c':{aa:0 for aa in amino_acids}}
    region_lengths = {'n':[],'h':[],'c':[]}
    critical_creg_positions = {-3:[],-2:[],-1:[]}

    # create one dict of collected SignalP 3.0 data for all sequences
    for file in args.infiles:
        region_aa_count_dict, region_lengths, critical_creg_positions = process_signalp3(file, region_aa_count_dict, region_lengths, critical_creg_positions)
    pdb.set_trace()
    df = pd.DataFrame(region_lengths)


    # plotting
    aa_alphabet = 'ARNDCQEGHILKMFPSTWYV'
    cmap = plt.get_cmap('tab20')
    aa_colors = {aa_alphabet[idx]:cmap(col_i) for idx, col_i in enumerate(np.linspace(0, 1, len(aa_alphabet)))}

    for reg in ('n', 'h', 'c'):
        aa_labels = sorted(list(region_aa_count_dict[reg].keys()))
        vals = [region_aa_count_dict[reg][x]/sum(region_aa_count_dict[reg].values()) for x in aa_labels]
        cols = [aa_colors[aa] for aa in aa_labels]

        cmap = plt.get_cmap('tab20')
        fig, ax = plt.subplots(figsize=(10, 6))
        positions = np.arange(len(region_aa_count_dict[reg]))
        ax.bar(positions, vals, color=cols)
        ax.set_xticks(positions)
        ax.set_xticklabels(aa_labels)
        ax.set_title(f"Amino acid frequencies in {reg.upper()} regions\nof native A. niger signal peptides", fontsize=18)
        #ax.set_ylabel(f"Frequency")
        plt.savefig(f"aa_freq_{reg}_reg.svg")


    sys.exit(1)
    # find most frequent length for each region
    optimal_legths = {reg:max(set(region_lengths[reg]), key=region_lengths[reg].count) for reg in ('n','h','c')}
    
    # find most frequent amino acids on critical C region positions
    optimal_creg_cleavage = ''.join([max(set(critical_creg_positions[pos]), key=critical_creg_positions[pos].count) for pos in sorted(critical_creg_positions)])

    good_sigpeps = list()
    n_rounds = 0
    glycoamylase = 'MSFRSLLALSGLVCTGLANVISKRATLDSWLSNEATVARTAILNNIGADGAWVSGADSGIVVASPSTDNPDYFYTWTRD' \
                   'SGLVLKTLVDLFRNGDTSLLSTIENYISAQAIVQGISNPSGDLSSGAGLGEPKFNVDETAYTGSWGRPQRDGPALRATA' \
                   'MIGFGQWLLDNGYTSTATDIVWPLVRNDLSYVAQYWNQTGYDLWEEVNGSSFFTIAVQHRALVEGSAFATAVGSSCSWC' \
                   'DSQAPEILCYLQSFWTGSFILANFDSSRSGKDANTLLGSIHTFDPEAACDDSTFQPCSPRALANHKEVVDSFRSIYTLN' \
                   'DGLSDSEAVAVGRYPEDTYYNGNPWFLCTLAAAEQLYDALYQWDKQGSLEVTDVSLDFFKALYSDAATGTYSSSSSTYS' \
                   'SIVDAVKTFADGFVSIVETHAASNGSMSEQYDKSDGEQLSARDLTWSYAALLTANNRRNSVVPASWGETSASSVPGTCA' \
                   'ATSAIGTYSSVTVTSWPSIVATGGTTTTATPTGSGSVTSTSKTTATASKTSTSTSSTSCTTPTAVAVTFDLTATTTYGE' \
                   'NIYLVGSISQLGDWETSDGIALSADKYTSSDPLWYVTVTLPAGESFEYKFIRIESDDSVEWESDPNREYTVPQACGTST' \
                   'ATVTDTWR'

    synthetic_sigpeps = synthetic_signalpep(region_aa_count_dict, optimal_legths, optimal_creg_cleavage, args.n_pep)
    sigpeps_gla = [pep + glycoamylase for pep in synthetic_sigpeps]

    # write signal peptides to fasta files of 10000 entries per file (for signalP 5.0)
    for chunk_no in range(0, args.n_pep, 10000):
        fasta = open(args.outfile+f'_{int(chunk_no/10000)}', 'w')
        fasta.write('\n'.join([f'>seq_{i}_\n' + sigpeps_gla[i] for i in range(chunk_no,chunk_no+10000)]))
        fasta.close()

    """
    # compute possible signal peptides as long as there is less than 10 good ones
    while len(good_sigpeps) < 10 and n_rounds < args.max_rounds:
        synthetic_sigpeps = synthetic_signalpep(region_aa_count_dict, optimal_legths, optimal_creg_cleavage, args.n_pep)
        sigpeps_gla = [pep+glycoamylase for pep in synthetic_sigpeps]

        # write 4 temporary fastas with 250 synthetic signal peptides + GlaA
        # run signalP 5.0 on these and extract signal peptides of >0.9 certainty
        #if not os.path.exists('tmp/'): os.mkdir('tmp')
        fasta = open(outfile, 'w')
        fasta.write('\n'.join([f'>seq_{i}\n'+sigpeps_gla[i] for i in range(len(sigpeps_gla))]))
        fasta.close()

        # run signalP 5.0 to extract synthetic signal peptides with >0.9 accuracy
        #process = subprocess.run(['signalp-5.0', '-t', 'euk', '-format', 'short', 'tmp/synthetic_sigpeps.fa'], capture_output=True)
        #pdb.set_trace()
        #signalp5_output = StringIO(input.communicate()[0].decode('utf-8'))
        #signalp5_output = process.stdout.decode('utf-8')


        #os.remove('tmp/synthetic_sigpeps.fa')

        n_rounds += 1
    """

    return



if __name__ == '__main__':
    print("# Running model for predicting synthetic signal peptides...")
    start_time = datetime.now()
    args = get_args()
    print("# args:", args)
    main(args)
    print(f"Done! Wrote fasta: {args.outfile}")
    end_time = datetime.now()
    print('# Duration: {}'.format(end_time - start_time))

