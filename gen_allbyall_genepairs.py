'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
'''

import pandas as pd
import csv, argparse
from itertools import combinations


def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--gene-list', type=str, help='.txt file with list of genes to make all by all combos')
	return parser


def read_args():
	args = get_parser().parse_args()
	gene_df = pd.read_csv(args.gene_list, sep='\t')
	return gene_df


def gen_gene_pairs(df):
	genes = list(df['Gene Symbol'])
	gene_pairs = list(combinations(genes, 2))
	return gene_pairs


def write_output(outputfile, gene_pairs, gene_df):
	with open(outputfile, 'w') as o:
		w = csv.writer(o, delimiter='\t')
		w.writerow(['Gene Symbol', 'Annotation'])
		for i, r in enumerate(gene_pairs):
			w.writerow([r[0] + ':' + r[1], 'GP'])
		for i, r in gene_df.iterrows():  ###Making same gene pairs
			w.writerow([r['Gene Symbol'] + ':' + r['Gene Symbol'], 'GP'])
	return


if __name__ == '__main__':
	gene_df = read_args()
	outputfile = 'gene_pairs.txt'
	write_output(outputfile, gen_gene_pairs(gene_df), gene_df)
