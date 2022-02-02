#!/usr/bin/python3
# _*_ coding: utf-8 _*_

import os
import sys
import time
import pandas as pd
from argparse import ArgumentParser

from mirsnp.mkaltfa import (fetch_transcript_fasta,
                            run_sub_variant_bases,
                            _fa_groups, _direct_groups)
from mirsnp.formatpt import merge_refalt
from mirsnp.utils import check_outputf, run_subshell
from mirsnp.trans_input import trans_input
from mirsnp.parse_gtf import parse_gtf
from mirsnp.config import DEFAULT_FLANK, DEFAULT_ENERGY
from mirsnp._version import __version__


def mir_snp(fpath, mirna_fa, out_dir, transcript_fa, gtf,
            flank=DEFAULT_FLANK, run_rnahybrid=True):
    """
    fpath: variants file
    mirna_fa: miRNA sequences, FASTA
    out_dir: main output directory
    transcript_fa: human genome exon region fasta，created by gffread
    gtf: gtf
    flank: flank sequence length from vairant position
    """
    fpath = os.path.abspath(fpath)
    out_dir = os.path.abspath(out_dir)

    seq_dir = os.path.join(out_dir, 'seq')
    if not os.path.exists(seq_dir):
        os.makedirs(seq_dir)

    print('##start {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))

    print('='*40)
    print('##input file: {}'.format(fpath))
    print('##ref fasta: {}'.format(transcript_fa))
    print('##gtf: {}'.format(gtf))
    print('##miRNA fasta: {}'.format(mirna_fa))
    print('##flank length: {}'.format(flank))
    print('##output directory: {}'.format(out_dir))
    print('='*40)

    print('transform variants file: ')
    variants_file = os.path.join(out_dir, 'variants.txt')
    if not os.path.exists(variants_file):
        trans_input(fpath, variants_file, gtf)

    print('flank length: {}'.format(flank))

    print('extract transcript sequence - ref: ')
    ref_fa = os.path.join(seq_dir, "ref.fa")
    fetch_transcript_fasta(variants_file, transcript_fa, ref_fa)
    check_outputf(ref_fa)

    print('extract transcript sequence - alt: ')
    alt_log = os.path.join(seq_dir, "alt.log")
    if not os.path.exists(alt_log):
        run_sub_variant_bases(variants_file, ref_fa, seq_dir)
        check_outputf(alt_log)

    # run RNAhybrid
    if run_rnahybrid:
        if not os.path.exists(mirna_fa):
            sys.exit('error: query miRNA fasta not found: {}'.format(mirna_fa))
        for gp in _fa_groups:
            for tp in _direct_groups:
                target_fa = os.path.join(seq_dir, f"{gp}-flank{flank}bp-{tp}.fa")
                output_pt = os.path.join(seq_dir, f"{gp}-flank{flank}bp-{tp}-pattern.txt")
                if not os.path.exists(output_pt):
                    print('run RNAhybrid {} - {}'.format(gp, tp))
                    cmd = 'RNAhybrid -f 2,7 -e {} -s 3utr_human -t "{}" -q {} > "{}"'.format(DEFAULT_ENERGY, target_fa, mirna_fa, output_pt)
                    run_subshell(cmd)

        trans_info = parse_gtf(gtf)
        stats_dir = os.path.join(out_dir, 'stats')
        if not os.path.exists(stats_dir):
            os.mkdir(stats_dir)

        df_lst = []
        for tp, tpnm in _direct_groups.items():
            ref_pt = os.path.join(seq_dir, f"ref-flank{flank}bp-{tp}-pattern.txt")
            alt_pt = os.path.join(seq_dir, f"alt-flank{flank}bp-{tp}-pattern.txt")
            print(f"parse RNAhybrid pattern, compare ref and alt - {tp}")
            dfm = merge_refalt(ref_pt, alt_pt, stats_dir, trans_info, name=tp, td=tpnm)
            df_lst.append(dfm)

        result_file = os.path.join(out_dir, 'result.txt')
        df_cat = pd.concat(df_lst, sort=False)
        df_f1 = df_cat.index.to_frame()['flag'].str.split(':', expand=True)
        df_f1.columns = ['target_name', 'miRNA', 'gene', 'strand', 'target_position']
        df_f2 = df_f1['target_name'].str.split('_', expand=True)
        df_f2.columns = ['transcript', 'rsid', 'chr', 'start_position', 'ref', 'alt', 'variant_position']
        df_res = df_f2.join(df_f1).join(df_cat)

        df_res.to_csv(result_file, sep='\t', index=False)
        check_outputf(result_file)

    print('##end {}'.format(time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())))


def abspath(p):
    return os.path.abspath(p)


def main(args):
    parser = ArgumentParser(description=f'Detecting the effect of gene mutation (SNP/InDel) on microRNA binding target genes')
    parser.add_argument('-v', '--version', dest='version', action='version', version=__version__)
    parser.add_argument('-i', dest='fpath', metavar='FILE', help='variants file')
    parser.add_argument('-f', dest='mirna', metavar='FILE', help='miRNA fasta')
    parser.add_argument('-r', dest='reffa', metavar='FILE', help='ref fasta(exon regions extracted by gffread)')
    parser.add_argument('-g', dest='gtf', metavar='FILE', help='a standard Ensembl Genomes gtf file')
    parser.add_argument('-o', dest='outdir', metavar='PATH', help='output directory')

    _args = parser.parse_args(args)

    mir_snp(fpath=abspath(_args.fpath),
            mirna_fa=abspath(_args.mirna),
            gtf=abspath(_args.gtf),
            out_dir=abspath(_args.outdir),
            transcript_fa=abspath(_args.reffa),
            run_rnahybrid=True,
            flank=DEFAULT_FLANK)
