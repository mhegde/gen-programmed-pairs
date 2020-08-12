'''
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
'''

import pandas as pd
import csv, argparse
from itertools import product
from matplotlib import rc

rc("pdf", fonttype=42)


def get_parser():
	parser = argparse.ArgumentParser()
	parser.add_argument('--gene-pairs', type=str, help='File with required gene pairs')
	parser.add_argument('--sg-picks', type=str, help='File with guide picks for all genes')
	parser.add_argument('--no-site', type=str, help='File with list of no-site controls')
	parser.add_argument('--intergenic', type=str, help='File with list of intergenic controls')
	parser.add_argument('--num-ctrls', type=int, default=5, help='Number of no-site and intergenic controls to pair each guide with')
	parser.add_argument('--ctrl-ctrl', type=int, default=100, help='Number of control-control combinations')
	parser.add_argument('--outputfile', type=str, help='Output file name')
	return parser


def read_args():
	args = get_parser().parse_args()
	gene_pairs = pd.read_csv(args.gene_pairs, sep='\t')
	sg_picks = pd.read_csv(args.sg_picks)
	no_site = pd.read_csv(args.no_site, sep='\t')
	intergenic = pd.read_csv(args.intergenic, sep='\t')
	num_ctrls = args.num_ctrls
	ctrl_ctrl = args.ctrl_ctrl
	outputfile = args. outputfile
	return gene_pairs, sg_picks, no_site, intergenic, num_ctrls, ctrl_ctrl, outputfile


def get_list(df):
	g_list = []
	for i,r in df.iterrows():
		g_list.append(r['sgRNA Sequence']+':'+r['Target Gene Symbol'])
	return g_list


def single_combos(df1, df2, dr, w):
	glist1 = get_list(df1)
	glist2 = get_list(df2)
	req_comb = list(product(glist1, glist2))
	write_output(req_comb, dr, w)
	return


def sg_ctrls_combos(df, ctrls1, ctrls2, dr, w):
	glist = get_list(df)
	ctrls1 = get_list(ctrls1)
	ctrls2 = get_list(ctrls2)
	for i, g in enumerate(glist):
		req_comb = list(product([g],ctrls1))
		write_output(req_comb, dr, w)
		req_comb = list(product([g],ctrls2))
		write_output(req_comb, dr, w)
	return


def ctrls_combos(ctrls1, ctrls2, dr, w):
	ctrls1 = get_list(ctrls1)
	ctrls2 = get_list(ctrls2)
	for (c1,c2) in zip(ctrls1, ctrls2):
		sg1, g1 = c1.split(':')
		sg2, g2 = c2.split(':')
		w.writerow([sg1+':'+dr+':'+sg2, g1+':'+g2])
	return


def write_output(comb_list, dr, w):
	for c in comb_list:
		sg1, g1 = c[0].split(':')
		sg2, g2 = c[1].split(':')
		w.writerow([sg1+':'+dr+':'+sg2, g1+':'+g2])
	return

if __name__ == '__main__':
	gene_pairs, sg_picks, no_site, intergenic, num_ctrls, ctrl_ctrl, outputfile = read_args()
	genes = {}
	with open(outputfile, 'w') as o:
		w = csv.writer(o, delimiter='\t')
		w.writerow(['sgRNA combinations', 'Gene combinations'])
		for i,r in gene_pairs.iterrows():
			if r['Annotation'] == 'GP':
				g1, g2 = r['Gene Symbol'].split(':')
				g1_df = sg_picks[sg_picks['Target Gene Symbol'] == g1]
				g2_df = sg_picks[sg_picks['Target Gene Symbol'] == g2]
				single_combos(g1_df, g2_df, 'TAATTTCTACTATCGTAGAT', w)
				if g1 not in genes.keys():
					sg_ctrls_combos(g1_df, no_site.sample(num_ctrls), intergenic.sample(num_ctrls), 'TAATTTCTACTATCGTAGAT', w)
					genes[g1] = 1
				if g2 not in genes.keys():
					sg_ctrls_combos(g2_df, no_site.sample(num_ctrls), intergenic.sample(num_ctrls), 'TAATTTCTACTATCGTAGAT', w)
					genes[g2] = 1
			elif r['Annotation'] == 'Essential':
				g1, g2 = r['Gene Symbol'].split(':')
				g1_df = sg_picks[sg_picks['Target Gene Symbol'] == g1]
				sg_ctrls_combos(g1_df, no_site.sample(num_ctrls), intergenic.sample(num_ctrls), 'TAATTTCTACTATCGTAGAT', w)
		nt_ctrls = no_site.sample(ctrl_ctrl)
		os_ctrls = intergenic.sample(ctrl_ctrl)
		ctrls_combos(nt_ctrls[:int(ctrl_ctrl/2)], nt_ctrls[int(ctrl_ctrl/2):ctrl_ctrl], 'TAATTTCTACTATCGTAGAT', w)
		ctrls_combos(os_ctrls[:int(ctrl_ctrl/2)], os_ctrls[int(ctrl_ctrl/2):ctrl_ctrl], 'TAATTTCTACTATCGTAGAT', w)

