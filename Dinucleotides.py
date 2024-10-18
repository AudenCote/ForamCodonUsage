import os, sys
from Bio import SeqIO
import dinuq
import numpy as np
import pandas as pd
import random
from tqdm import tqdm
from copy import deepcopy
from statistics import mean
import math
import argparse


universal_no6fold = {
	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
	'GCG': ['A', 'four', 0], 'CGT': ['R', 'four', 0], 'CGC': ['R', 'four', 0],
	'CGG': ['R', 'four', 0], 'CGA': ['R', 'four', 0], 'AGA': ['R', 'two', 0],
	'AGG': ['R', 'two', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
	'TTA': ['L', 'two', 0], 'TTG': ['L', 'two', 0], 'CTT': ['L', 'four', 0],
	'CTC': ['L', 'four', 0], 'CTA': ['L', 'four', 0], 'CTG': ['L', 'four', 0],
	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'four', 0],
	'TCC': ['S', 'four', 0], 'TCA': ['S', 'four', 0], 'TCG': ['S', 'four', 0],
	'AGT': ['S', 'two', 0], 'AGC': ['S', 'two', 0], 'ACT': ['T', 'four', 0],
	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}

aas = list(dict.fromkeys([universal_no6fold[codon][0] for codon in universal_no6fold]))
codons_per_aa = { aa : [codon for codon in universal_no6fold if universal_no6fold[codon][0] == aa] for aa in aas }

codon_folds = { codon : { 0 : len(list(dict.fromkeys([c[0] for c in codons_per_aa[aa] if c[1:3] == codon[1:3]]))), 1 : len(list(dict.fromkeys([c[1] for c in codons_per_aa[aa] if c[0] == codon[0] and c[2] == codon[2]]))), 2 : len(list(dict.fromkeys([c[2] for c in codons_per_aa[aa] if c[0:2] == codon[0:2]])))} for aa in codons_per_aa for codon in codons_per_aa[aa] }
codon_folds = { codon: { pos : codon_folds[codon][pos] if codon_folds[codon][pos] != 2 else ('2TC' if codon[pos] in 'TC' else '2AG') for pos in codon_folds[codon] } for codon in codon_folds }

def percentile(li, percentile):
    n = len(li)
    idx = n * percentile / 100
    return sorted(li)[math.floor(idx)]

def persite_sims(fasta_file):

	alldns = ['AA', 'AT', 'AC', 'AG', 'TT', 'TA', 'TC', 'TG', 'CC', 'CG', 'CA', 'CT', 'GG', 'GC', 'GT', 'GA']

	dt = { d : 0 for d in alldns }
	nt = { 1 : { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 }, '2TC' : { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 }, '2AG' : { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 }, 3 : { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 }, 4 : { 'A' : 0, 'T' : 0, 'G' : 0, 'C' : 0 } }

	with open('Databases/SSDU.csv', 'w') as o:
		o.write('Seq,Taxon,TPM,Pos,Dn,ObsFreq,ExpFreq\n')#Mean,CI95Low,CI95High\n')
		for rec in tqdm([r for r in SeqIO.parse(fasta_file, 'fasta')]):
			tax = rec.id[:4] + rec.id[5:10]
			tpm = rec.id.split('TPM')[-1].split('_')[0]

			pos1_sn = deepcopy(nt); pos2_sn = deepcopy(nt); pos3_sn = deepcopy(nt)
			pos1_dn = deepcopy(dt); pos2_dn = deepcopy(dt); pos3_dn = deepcopy(dt)

			aa_seq = ''

			for i, nuc in enumerate(rec.seq):
				if i % 3 == 0:

					codon = str(rec.seq[i:i+3])
					if 'N' not in codon:
						aa_seq = aa_seq + universal_no6fold[codon][0]

						pos1_sn[codon_folds[codon][0]][str(rec.seq[i])] += 1
						pos2_sn[codon_folds[codon][1]][str(rec.seq[i + 1])] += 1
						pos3_sn[codon_folds[codon][2]][str(rec.seq[i + 2])] += 1

						pos1_dn[str(rec.seq[i:i+2])] += 1
						pos2_dn[str(rec.seq[i + 1:i+3])] += 1

						try:
							pos3_dn[str(rec.seq[i + 2:i+4])] += 1
						except:
							continue

			pos1_sn = { s : { k : pos1_sn[s][k]/sum(pos1_sn[s].values()) for k in pos1_sn[s] } if sum(pos1_sn[s].values()) != 0 else { k : 0 for k in pos1_sn[s] } for s in pos1_sn  }
			pos2_sn = { s : { k : pos2_sn[s][k]/sum(pos2_sn[s].values()) for k in pos2_sn[s] } if sum(pos2_sn[s].values()) != 0 else { k : 0 for k in pos2_sn[s] } for s in pos2_sn  }
			pos3_sn = { s : { k : pos3_sn[s][k]/sum(pos3_sn[s].values()) for k in pos3_sn[s] } if sum(pos3_sn[s].values()) != 0 else { k : 0 for k in pos3_sn[s] } for s in pos3_sn  }
			pos1_dn = { k : pos1_dn[k]/sum([pos1_dn[k2] for k2 in pos1_dn]) for k in pos1_dn }
			pos2_dn = { k : pos2_dn[k]/sum([pos2_dn[k2] for k2 in pos2_dn]) for k in pos2_dn }
			pos3_dn = { k : pos3_dn[k]/sum([pos3_dn[k2] for k2 in pos3_dn]) for k in pos3_dn }

			sims = { 'pos1' : [], 'pos2' : [], 'pos3' : [] }
			for i in range(100):
				sim_seq = ''

				for aa in aa_seq:

					codon_options = codons_per_aa[aa]
					codon_freqs = { codon : pos1_sn[codon_folds[codon][0]][codon[0]] * pos2_sn[codon_folds[codon][1]][codon[1]] * pos3_sn[codon_folds[codon][2]][codon[2]] for codon in codon_options}

					sim_seq = sim_seq + random.choices(list(codon_freqs.keys()), weights=codon_freqs.values())[0]

				sim_pos1_dn = deepcopy(dt); sim_pos2_dn = deepcopy(dt); sim_pos3_dn = deepcopy(dt)

				for i, nuc in enumerate(sim_seq):
					if i % 3 == 0:

						sim_pos1_dn[str(sim_seq[i:i+2])] += 1
						sim_pos2_dn[str(sim_seq[i + 1:i+3])] += 1

						try:
							sim_pos3_dn[str(sim_seq[i + 2:i+4])] += 1
						except:
							continue

				sim_pos1_dn = { k : sim_pos1_dn[k]/sum([sim_pos1_dn[k2] for k2 in sim_pos1_dn]) for k in sim_pos1_dn }
				sim_pos2_dn = { k : sim_pos2_dn[k]/sum([sim_pos2_dn[k2] for k2 in sim_pos2_dn]) for k in sim_pos2_dn }
				sim_pos3_dn = { k : sim_pos3_dn[k]/sum([sim_pos3_dn[k2] for k2 in sim_pos3_dn]) for k in sim_pos3_dn }

				sims['pos1'].append(sim_pos1_dn)
				sims['pos2'].append(sim_pos2_dn)
				sims['pos3'].append(sim_pos3_dn)

			sim_pos1_dn = { dn : [subdict[dn] for subdict in sims['pos1']] for dn in alldns }
			sim_pos2_dn = { dn : [subdict[dn] for subdict in sims['pos2']] for dn in alldns }
			sim_pos3_dn = { dn : [subdict[dn] for subdict in sims['pos3']] for dn in alldns }

			

			sim_pos1_dn = { dn : [mean(sim_pos1_dn[dn]), percentile(sim_pos1_dn[dn], 2.5), percentile(sim_pos1_dn[dn], 97.5)] for dn in sim_pos1_dn }
			sim_pos2_dn = { dn : [mean(sim_pos2_dn[dn]), percentile(sim_pos2_dn[dn], 2.5), percentile(sim_pos2_dn[dn], 97.5)] for dn in sim_pos1_dn }
			sim_pos3_dn = { dn : [mean(sim_pos3_dn[dn]), percentile(sim_pos3_dn[dn], 2.5), percentile(sim_pos3_dn[dn], 97.5)] for dn in sim_pos1_dn }

			for dn in alldns:
				o.write(rec.id + ',' + tax + ',' + tpm + ',Pos1,' + dn + ',' + str(pos1_dn[dn]) + ',' + str(sim_pos1_dn[dn][0]) + ',' + str(sim_pos1_dn[dn][1]) + ',' + str(sim_pos1_dn[dn][2]) + '\n')
				o.write(rec.id + ',' + tax + ',' + tpm + ',Pos2,' + dn + ',' + str(pos2_dn[dn]) + ',' + str(sim_pos2_dn[dn][0]) + ',' + str(sim_pos2_dn[dn][1]) + ',' + str(sim_pos2_dn[dn][2]) + '\n')
				o.write(rec.id + ',' + tax + ',' + tpm + ',Pos3,' + dn + ',' + str(pos3_dn[dn]) + ',' + str(sim_pos3_dn[dn][0]) + ',' + str(sim_pos3_dn[dn][1]) + ',' + str(sim_pos3_dn[dn][2]) + '\n')


if __name__ == '__main__':
	persite_sims('Databases/AllCuratedNTDSeqs_Seq70_ATGTrimmed_TPMV2Raw.fasta')



