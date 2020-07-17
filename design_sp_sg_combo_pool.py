'''
This script helps design combination pools with specific guide pairs. It also checks for sequence similarity between
guides
Author: Mudra Hegde
Email: mhegde@broadinstitute.org
'''
import pandas as pd
import csv, argparse, random
from difflib import SequenceMatcher
from itertools import product


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sp-input-file',
        type=str,
        help='File with S.pyogenes guide designs')
    parser.add_argument('--sa-input-file',
        type=str,
        help='File with S.aureus guide designs')
    parser.add_argument('--sp-nosite-file',
        type=str,
        help='File with S.pyogenes NO_SITE control guides')
    parser.add_argument('--sp-onesite-file',
        type=str,
        help='File with S.pyogenes ONE_NON-GENE_SITE control guides')
    parser.add_argument('--sa-nosite-file',
        type=str,
        help='File with S.aureus NO_SITE control guides')
    parser.add_argument('--sa-onesite-file',
        type=str,
        help='File with S.aureus ONE_NON-GENE_SITE control guides')
    parser.add_argument('--gene-pairs',
        type=str,
        help='File with required gene pairs')
    parser.add_argument('--overlap-check',
        type=int,
        default=12, 
        help='Length of overlap to be avoided')
    parser.add_argument('--num-ctrls',
        type=int,
        default=5,
        help='Number of control guides to be paired with targeting guides')
    parser.add_argument('--ctrl-ctrl',
        type=int,
        default=50,
        help='Number of control by control guide combinations'
    parser.add_argument('--outputfile',
        type=str,
        help='Outputfile')
    return parser


def revcom(s):
    basecomp = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    letters = list(s[::-1])
    letters = [basecomp[base] for base in letters]
    return ''.join(letters)


def get_comb(df1, df2):
    list1, list2 = [], []
    for i, r in df1.iterrows():
        list1.append(r['sgRNA Sequence'] + ':' + r['Target Gene Symbol'])
    for i, r in df2.iterrows():
        list2.append(r['sgRNA Sequence'] + ':' + r['Target Gene Symbol'])
    req_comb = list(product(list1, list2))
    return req_comb


def write_output(comb_list, o_check, w):
    for c in comb_list:
        sg1, g1 = c[0].split(':')
        sg2, g2 = c[1].split(':')
        s = SequenceMatcher(None, sg1, sg2)  # Sequence overlap check
        s_rev = SequenceMatcher(None, revcom(sg1), sg2)  # Sequence overlap check of reverse complement
        if (s.find_longest_match(0, 20, 0, 21).size < o_check) and (
                s_rev.find_longest_match(0, 20, 0, 21).size < o_check):
            w.writerow([sg1+':'+sg2, g1+':'+g2])
        else:
            print(sg1 + ':' + sg2)
    return w


def get_sg_pairs(sp, sa, sp_nosite, sp_onesite, sa_nosite, sa_onesite, gp, o_check, num_ctrls, ctrl_ctrl, outputfile):
    pair_check = {}
    with open(outputfile,'w') as o:
        w = csv.writer(o,delimiter='\t')
        w.writerow(['sgRNA Combination','Gene pair'])
        for i,g in gp.iterrows():
            if g['Annotation'] != 'Essential':
                print(g['Gene Symbol'])
                g1,g2 = g['Gene Symbol'].split(':')
                g1_sp = sp[sp['Target Gene Symbol']==g1]
                g2_sa = sa[sa['Target Gene Symbol']==g2]
                ori1_combo = get_comb(g1_sp, g2_sa)
                write_output(ori1_combo, o_check, w)
                if 'control:'+g2 not in pair_check.keys():
                    g2_sa_nosite_combo = get_comb(sp_nosite.sample(num_ctrls), g2_sa)
                    write_output(g2_sa_nosite_combo, o_check, w)
                    g2_sa_onesite_combo = get_comb(sp_onesite.sample(num_ctrls), g2_sa)
                    write_output(g2_sa_onesite_combo, o_check, w)
                    pair_check['control:'+g2] = 1
                if g1 + ':control' not in pair_check.keys():
                    g1_sp_nosite_combo = get_comb(g1_sp, sa_nosite.sample(num_ctrls))
                    write_output(g1_sp_nosite_combo, o_check, w)
                    g1_sp_onesite_combo = get_comb(g1_sp, sa_onesite.sample(num_ctrls))
                    write_output(g1_sp_onesite_combo, o_check, w)
                    pair_check[g1 + ':control'] = 1
                g2_sp = sp[sp['Target Gene Symbol']==g2]
                g1_sa = sa[sa['Target Gene Symbol']==g1]
                ori2_combo = get_comb(g2_sp, g1_sa)
                write_output(ori2_combo, o_check, w)
                if 'control:'+g1 not in pair_check.keys():
                    g1_sa_nosite_combo = get_comb(sp_nosite.sample(num_ctrls), g1_sa)
                    write_output(g1_sa_nosite_combo, o_check, w)
                    g1_sa_onesite_combo = get_comb(sp_onesite.sample(num_ctrls), g1_sa)
                    write_output(g1_sa_onesite_combo, o_check, w)
                    pair_check['control:'+g1] = 1
                if g2+':control' not in pair_check.keys():
                    g2_sp_nosite_combo = get_comb(g2_sp, sa_nosite.sample(num_ctrls))
                    write_output(g2_sp_nosite_combo, o_check, w)
                    g2_sp_onesite_combo = get_comb(g2_sp, sa_onesite.sample(num_ctrls))
                    write_output(g2_sp_onesite_combo, o_check, w)
                    pair_check[g2 + ':control'] = 1
            else:
                g1, g2 = g['Gene Symbol'].split(':')
                g1_sp = sp[sp['Target Gene Symbol'] == g1]
                if g1+':control' not in pair_check.keys():
                    g1_sp_nosite_combo = get_comb(g1_sp, sa_nosite.sample(num_ctrls))
                    write_output(g1_sp_nosite_combo, o_check, w)
                    g1_sp_onesite_combo = get_comb(g1_sp, sa_onesite.sample(num_ctrls))
                    write_output(g1_sp_onesite_combo, o_check, w)
                    pair_check[g1 + ':control'] = 1
                g1_sa = sa[sa['Target Gene Symbol'] == g1]
                if 'control:'+g1 not in pair_check.keys():
                    g1_sa_nosite_combo = get_comb(sp_nosite.sample(num_ctrls), g1_sa)
                    write_output(g1_sa_nosite_combo, o_check, w)
                    g1_sa_onesite_combo = get_comb(sp_onesite.sample(num_ctrls), g1_sa)
                    write_output(g1_sa_onesite_combo, o_check, w)
                    pair_check['control:'+g1] = 1
        #get_control_by_control(sp_nosite.sample(50), sa_nosite.sample(50), o_check, w)
        ctrl_comb = get_comb(sp_nosite.sample(ctrl_ctrl), sa_nosite.sample(ctrl_ctrl))
        write_output(ctrl_comb, o_check, w)
        ctrl_comb = get_comb(sp_onesite.sample(ctrl_ctrl), sa_onesite.sample(ctrl_ctrl))
        write_output(ctrl_comb, o_check, w)
        #get_control_by_control(sp_onesite.sample(50), sa_onesite.sample(50), o_check, w)
    return w


def get_control_by_control(ctrls1, ctrls2, o_check, w):
    list1, list2 = [], []
    for i, r in ctrls1.iterrows():
        list1.append(r['sgRNA Sequence'] + ':' + r['Target Gene Symbol'])
    for i, r in ctrls2.iterrows():
        list2.append(r['sgRNA Sequence'] + ':' + r['Target Gene Symbol'])
    for (c1,c2) in zip(list1, list2):
        sg1, g1 = c1.split(':')
        sg2, g2 = c2.split(':')
        s = SequenceMatcher(None, sg1, sg2)  # Sequence overlap check
        s_rev = SequenceMatcher(None, revcom(sg1), sg2)  # Sequence overlap check of reverse complement
        if (s.find_longest_match(0, 20, 0, 21).size < o_check) and (
                s_rev.find_longest_match(0, 20, 0, 21).size < o_check):
            w.writerow([sg1+':'+sg2, g1+':'+g2])
    return


def read_args(args):
    sp_input = pd.read_csv(args.sp_input_file)
    sa_input = pd.read_csv(args.sa_input_file)
    sp_nosite = pd.read_csv(args.sp_nosite_file, sep='\t')
    sp_onesite = pd.read_csv(args.sp_onesite_file, sep='\t')
    sa_nosite = pd.read_csv(args.sa_nosite_file, sep='\t')
    sa_onesite = pd.read_csv(args.sa_onesite_file, sep='\t')
    gene_pairs = pd.read_table(args.gene_pairs)
    overlap_check = args.overlap_check
    num_ctrls = args.num_ctrls
    ctrl_ctrl = args.ctrl_ctrl
    outputfile = args.outputfile
    return sp_input, sa_input, sp_nosite, sp_onesite, sa_nosite, sa_onesite, gene_pairs, overlap_check, num_ctrls, ctrl_ctrl, outputfile


if __name__ == '__main__':
    args = get_parser().parse_args()
    sp_input, sa_input, sp_nosite, sp_onesite, sa_nosite, sa_onesite, gene_pairs, overlap_check, num_ctrls, ctrl_ctrl, outputfile = read_args(args)
    val = get_sg_pairs(sp_input, sa_input, sp_nosite, sp_onesite, sa_nosite, sa_onesite, gene_pairs, overlap_check, num_ctrls, ctrl_ctrl, outputfile)



