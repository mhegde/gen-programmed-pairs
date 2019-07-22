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
    parser.add_argument('--sp-ctrls-file',
        type=str,
        help='File with S.pyogenes control guides')
    parser.add_argument('--sa-ctrls-file',
        type=str,
        help='File with S.aureus control guides')
    parser.add_argument('--gene-pairs',
        type=str,
        help='File with required gene pairs')
    parser.add_argument('--overlap-check',
        type=int,
        default=12, 
        help='Length of overlap to be avoided')
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


def get_sg_pairs(sp, sa, sp_ctrl, sa_ctrl, gp, o_check, outputfile):
    pair_check = {}
    sp_orig = sp_ctrl
    sa_orig = sa_ctrl
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
                if 'non-targeting:'+g2 not in pair_check.keys():
                    if len(sp_ctrl) <= len(g2_sa):
                        sp_ctrl = sp_orig
                    w, sp_ctrl = get_control_combos(g2_sa, sp_ctrl, 'sa', o_check, w)
                    pair_check['non-targeting:'+g2] = 1
                if g1 + ':non-targeting' not in pair_check.keys():
                    if len(sa_ctrl) == 0:
                        sa_ctrl = sa_orig
                    w, sa_ctrl = get_control_combos(g1_sp, sa_ctrl, 'sp', o_check, w)
                    pair_check[g1 + ':non-targeting'] = 1
                g2_sp = sp[sp['Target Gene Symbol']==g2]
                g1_sa = sa[sa['Target Gene Symbol']==g1]
                ori2_combo = get_comb(g2_sp, g1_sa)
                write_output(ori2_combo, o_check, w)
                if 'non-targeting:'+g1 not in pair_check.keys():
                    if len(sp_ctrl) <= len(g1_sa):
                        sp_ctrl = sp_orig
                    w, sp_ctrl = get_control_combos(g1_sa, sp_ctrl, 'sa', o_check, w)
                    pair_check['non-targeting:'+g1] = 1
                if g2+':non-targeting' not in pair_check.keys():
                    if len(sa_ctrl) <= len(g2_sp):
                        sa_ctrl = sa_orig
                    w, sa_ctrl = get_control_combos(g2_sp, sa_ctrl, 'sp', o_check, w)
                    pair_check[g2 + ':non-targeting'] = 1
            else:
                g1, g2 = g['Gene Symbol'].split(':')
                g1_sp = sp[sp['Target Gene Symbol'] == g1]
                if g1+':non-targeting' not in pair_check.keys():
                    if len(sa_ctrl) <= len(g1_sp):
                        sa_ctrl = sa_orig
                    w, sa_ctrl = get_control_combos(g1_sp, sa_ctrl, 'sp', o_check, w)
                    pair_check[g1 + ':non-targeting'] = 1
                g1_sa = sa[sa['Target Gene Symbol'] == g1]
                if 'non-targeting:'+g1 not in pair_check.keys():
                    if len(sp_ctrl) <= len(g1_sa):
                        sa_ctrl = sp_ctrl
                    w, sp_ctrl = get_control_combos(g1_sa, sp_ctrl, 'sa', o_check, w)
                    pair_check['non-targeting:'+g1] = 1
        if len(sp_ctrl) <= 100:
            sp_ctrl = sp_orig
        if len(sa_ctrl) <= 100:
            sa_ctrl = sa_orig
        get_control_by_control(sp_ctrl, sa_ctrl, o_check, w)
    return w


def get_control_combos(df1, ctrls, sg_type, o_check, w):
    df1.index = list(range(0,len(df1)))
    ran_ctrls = random.sample(ctrls, len(df1))
    new_ctrls = [x for x in ctrls if x not in ran_ctrls]
    req_comb = []
    if sg_type == 'sp':
        for i,r in df1.iterrows():
            req_comb.append((r['sgRNA Sequence']+':'+r['Target Gene Symbol'], ran_ctrls[i]+':non-targeting'))
    else:
        for i, r in df1.iterrows():
            req_comb.append((ran_ctrls[i] + ':non-targeting', r['sgRNA Sequence'] + ':' + r['Target Gene Symbol']))
    w = write_output(req_comb, o_check, w)
    return w, new_ctrls


def get_control_by_control(ctrls1, ctrls2, o_check, w):
    ran_ctrls_1 = random.sample(ctrls1, 100)
    ran_ctrls_2 = random.sample(ctrls2, 100)
    for i, r in enumerate(ran_ctrls_1):
        control_2 = ran_ctrls_2[i]
        s = SequenceMatcher(None, r, control_2)  # Sequence overlap check
        s_rev = SequenceMatcher(None, revcom(r), control_2)  # Sequence overlap check of reverse complement
        if (s.find_longest_match(0, 20, 0, 21).size < o_check) and (
                s_rev.find_longest_match(0, 20, 0, 21).size < o_check):
            w.writerow([r+':'+control_2, 'non-targeting:non-targeting'])
    return


def read_args(args):
    sp_input = pd.read_csv(args.sp_input_file)
    sa_input = pd.read_csv(args.sa_input_file)
    sp_ctrls = pd.read_csv(args.sp_ctrls_file, sep='\t')
    sp_ctrls = list(sp_ctrls['sgRNA Seq'])
    sa_ctrls = pd.read_csv(args.sa_ctrls_file, sep='\t',header=1)
    sa_ctrls = list(sa_ctrls['sgRNA Seq'])
    gene_pairs = pd.read_table(args.gene_pairs)
    overlap_check = args.overlap_check
    outputfile = args.outputfile
    return sp_input, sa_input, sp_ctrls, sa_ctrls, gene_pairs, overlap_check, outputfile


if __name__ == '__main__':
    args = get_parser().parse_args()
    sp_input, sa_input, sp_ctrls, sa_ctrls, gene_pairs, overlap_check, outputfile = read_args(args)
    val = get_sg_pairs(sp_input, sa_input, sp_ctrls, sa_ctrls, gene_pairs, overlap_check, outputfile)


