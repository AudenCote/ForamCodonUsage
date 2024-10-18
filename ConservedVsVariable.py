import os
import sys
from Bio import SeqIO
from statistics import mean
from tqdm import tqdm


genetic_code = {
        	'GCT': ['A', 'four', 0], 'GCC': ['A', 'four', 0], 'GCA': ['A', 'four', 0],
        	'GCG': ['A', 'four', 0], 'CGT': ['R', 'four', 0], 'CGC': ['R', 'four', 0],
        	'CGG': ['R', 'four', 0], 'CGA': ['R', 'four', 0], 'AGA': ['R_', 'two', 0],
        	'AGG': ['R_', 'two', 0], 'AAT': ['N', 'two', 0], 'AAC': ['N', 'two', 0],
        	'GAT': ['D', 'two', 0], 'GAC': ['D', 'two', 0], 'TGT': ['C', 'two', 0],
        	'TGC': ['C', 'two', 0], 'CAA': ['Q', 'two', 0], 'CAG': ['Q', 'two', 0],
        	'GAA': ['E', 'two', 0], 'GAG': ['E', 'two', 0], 'GGT': ['G', 'four', 0],
        	'GGC': ['G', 'four', 0], 'GGA': ['G', 'four', 0], 'GGG': ['G', 'four', 0],
        	'CAT': ['H', 'two', 0], 'CAC': ['H', 'two', 0], 'ATT': ['I', 'three', 0],
        	'ATC': ['I', 'three', 0], 'ATA': ['I', 'three', 0], 'ATG': ['M', 'one', 0],
        	'TTA': ['L_', 'two', 0], 'TTG': ['L_', 'two', 0], 'CTT': ['L', 'four', 0],
        	'CTC': ['L', 'four', 0], 'CTA': ['L', 'four', 0], 'CTG': ['L', 'four', 0],
        	'AAA': ['K', 'two', 0], 'AAG': ['K', 'two', 0], 'TTT': ['F', 'two', 0],
        	'TTC': ['F', 'two', 0], 'CCT': ['P', 'four', 0], 'CCC': ['P', 'four', 0],
        	'CCA': ['P', 'four', 0], 'CCG': ['P', 'four', 0], 'TCT': ['S', 'four', 0],
        	'TCC': ['S', 'four', 0], 'TCA': ['S', 'four', 0], 'TCG': ['S', 'four', 0],
        	'AGT': ['S_', 'two', 0], 'AGC': ['S_', 'two', 0], 'ACT': ['T', 'four', 0],
        	'ACC': ['T', 'four', 0], 'ACA': ['T', 'four', 0], 'ACG': ['T', 'four', 0],
        	'TGG': ['W', 'one', 0], 'TAT': ['Y', 'two', 0], 'TAC': ['Y', 'two', 0],
        	'GTT': ['V', 'four', 0], 'GTC': ['V', 'four', 0], 'GTA': ['V', 'four', 0],
        	'GTG': ['V', 'four', 0], 'TAA': ['*', 'none', 0], 'TGA': ['*', 'none', 0],
        	'TAG': ['*', 'none', 0], 'XXX': ['_missing', 'none', 0]}


nucls = { rec.id : str(rec.seq) for rec in SeqIO.parse('Databases/AllCuratedNTDSeqs_Seq70_ATGTrimmed_TPMV2Raw.fasta', 'fasta') }
rec_corr = { line.split(',')[1].strip() : line.split(',')[0].strip() for line in open('UTR_RecID_Corr.csv') }

def group_sites():

	with open('Databases/Forams_Akashi.csv', 'w') as o:
		o.write('Seq,Conserved.GC,Conserved.Total,Variable.GC,Variable.Total\n')
		for file in tqdm(os.listdir('PostCurationAlignments')):
			if('.fas' in file):
				sites = { }
				for rec in SeqIO.parse('PostCurationAlignments/' + file, 'fasta'):
					for c, char in enumerate(rec.seq):
						if(c not in sites):
							sites.update({ c : { } })

						if(char not in sites[c]):
							sites[c].update({ char : 0 })

						sites[c][char] += 1

				conserved = []; variable = []
				for site in sites:
					if(max([sites[site][char]/sum(sites[site].values()) for char in sites[site] if char != '-']) >= .8):
						conserved.append(site)
					elif(max([sites[site][char]/sum(sites[site].values()) for char in sites[site] if char != '-']) <= .2):
						variable.append(site)


				if(len(conserved) > 1 and len(variable) > 1):

					for rec in SeqIO.parse('PostCurationAlignments/' + file, 'fasta'):
						try:
							rec.id = rec_corr[rec.id]
							if rec.id in nucls:
								gc_cons = 0; total_cons = 0; gc_var = 0; total_var = 0
								c = 0
								for i, char in enumerate(rec.seq):
									if i in conserved:
										while nucls[rec.id][c:c+3] in genetic_code and char != genetic_code[nucls[rec.id][c:c+3]][0]:
											c += 3
										if nucls[rec.id][c:c+3] in genetic_code and genetic_code[nucls[rec.id][c:c+3]][1] == 'four':
											total_cons += 1
											if nucls[rec.id][c:c+3][2] in 'GCgc':
												gc_cons += 1
									elif(i in variable):
										while nucls[rec.id][c:c+3] in genetic_code and char != genetic_code[nucls[rec.id][c:c+3]][0]:
											c += 3
										if nucls[rec.id][c:c+3] in genetic_code and genetic_code[nucls[rec.id][c:c+3]][1] == 'four':
											total_var += 1
											if nucls[rec.id][c:c+3][2] in 'GCgc':
												gc_var += 1
									
									if(char != '-'):
										c += 3


							o.write(rec.id + ',' + str(gc_cons) + ',' + str(total_cons) + ',' + str(gc_var) + ',' + str(total_var) + '\n')
						except Exception as e:
							print(e) 

group_sites()















