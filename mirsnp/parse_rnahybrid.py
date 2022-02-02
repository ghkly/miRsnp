#!/usr/bin/python3
# _*_ coding: utf-8 _*_


import os
import re
import sys
import gzip
import pandas as pd

from mirsnp.utils import check_outputf
from mirsnp.parse_gtf import parse_gtf


PATTERN_BASES = "AGCTUN"


def get_trans_name(al):
    trans_name = re.search(r"target: (.*?)\n", al).group(1)
    return trans_name


def get_mirna_name(al):
    mirna_name = re.search(r"miRNA : (.*?)\n", al).group(1)
    return mirna_name


def get_position(al):
    position = re.search(r"position  (\d+)\n", al).group(1)
    return position


def get_pvalue(al):
    pvalue = re.search(r"p-value: (.*?)\n", al).group(1)
    return pvalue


def get_mfe(al):
    mfe = re.search(r"mfe: (.*?) kcal/mol", al).group(1)
    return mfe


def get_length(al):
    """stat target seq legth and miRNA seq length"""
    lens = re.findall(r"length: (\d+)\n", al)
    target_len = int(lens[0])
    mirna_len = int(lens[1])
    return target_len, mirna_len


def str2list(seq):
    return list(seq)


def get_start_pos(seq, start=0):
    """position of the first base of target seq
    seq: target 5'   G  UGC  GG ... 3'
         miRNA  3' ACA  U    UU ... 5'
    """
    pos = 10
    for n, v in enumerate(seq):
        if v in PATTERN_BASES and n >= start:
            pos = n
            break
    return pos


def repair_end(seq1, seq2):
    """
    seq1:    target 5'   G  UGC  GG     CCAUCCAAGA       A 3'
    seq2:                 AG   UG  GGAGG          UGGCAUG
    seq2:                 UC   AC  CUUCC          AUUGUAC
    seq1:    miRNA  3' ACA  U    UU                        5'
    """
    seq1 = str2list(seq1)
    seq2 = str2list(seq2)
    vseq1 = list(reversed(seq1))
    vseq2 = list(reversed(seq2))

    p1 = get_start_pos(vseq1)
    p2 = get_start_pos(vseq2)

    if p1 > 3 and p2 > 3:
        for i in range(3, min([p1, p2])):
            vseq1[i] = '-'
    return ''.join(reversed(vseq1))


def repair_start(seq1, seq2):
    """
    seq1:    target 5'   G  UGC  GG     CCAUCCAAGA       A 3'
    seq2:                 AG   UG  GGAGG          UGGCAUG
    seq2:                 UC   AC  CUUCC          AUUGUAC
    seq1:    miRNA  3' ACA  U    UU                        5'
    """
    seq1 = str2list(seq1)
    seq2 = str2list(seq2)
    p1 = get_start_pos(seq1, start=8)
    p2 = get_start_pos(seq2, start=8)
    if p1 > 10 and p2 > 10:
        for i in range(10, min([p1, p2])):
            seq1[i] = '-'
    return ''.join(seq1)


def get_rna_seq(seq):
    """
    target 5' --G  UGC  GG     CCAUCCAAGA       A 3'
    ->
    --G  UGC  GG     CCAUCCAAGA       A
    """
    seq = str2list(seq)[10:-3]
    return ''.join(seq)


def merge_rna_seq(seq1, seq2):
    seq1 = str2list(seq1)
    seq2 = str2list(seq2)
    zp = list(zip(seq1, seq2))
    seq = ''
    for z in zp:
        zs = set(list(filter(lambda x: x != " ", z)))
        if len(zs) == 0:
            seq = seq + '-'
        elif len(zs) == 1:
            seq = seq + list(zs)[0]
        else:
            print("提示，序列重叠!")
            sys.exit()
    return seq


def get_merged_seq(al):
    """
    al:
    target 5'   G  UGC  GG     CCAUCCAAGA       A 3'
             AG   UG  GGAGG          UGGCAUG
             UC   AC  CUUCC          AUUGUAC
    miRNA  3' ACA  U    UU                        5'
    """
    para = re.search("target 5' (.*?) 5'\n", al, re.S).group(0)
    arrs = para.strip('\n').split('\n')
    arrs = list(filter(lambda x: x != '', arrs))
    arr1, arr2, arr3, arr4 = arrs

    mod_seq1 = repair_start(arr1, arr2)
    repair_seq1 = repair_end(mod_seq1, arr2)
    mod_seq2 = repair_start(arr4, arr3)
    repair_seq2 = repair_end(mod_seq2, arr3)

    s1 = get_rna_seq(repair_seq1)
    s2 = get_rna_seq(arr2)
    s3 = get_rna_seq(arr3)
    s4 = get_rna_seq(repair_seq2)

    target_seq = merge_rna_seq(s1, s2)
    mirna_seq = merge_rna_seq(s4, s3)
    return target_seq, mirna_seq


def get_seq_pos(start, target_length, seq_length):
    end_pos = start + seq_length - 1
    tail_dist = target_length - end_pos

    rev_start = tail_dist + 1
    rev_end = rev_start + seq_length - 1

    return [rev_start, rev_end]


def get_segs(exons):
    dist = [0]
    for ex in exons:
        s, e = ex.split(':')
        exlen = int(e) - int(s) + 1
        dist.append(exlen)

    segs = []
    for i in range(len(dist) - 1):
        segs.append([sum(dist[: i + 1]) + 1, sum(dist[:i + 2])])

    segs = [f"{i[0]}-{i[1]}" for i in segs]

    return segs


def get_chr_regions(seq_pos, segs, exons):
    regions = []
    for n, seg in enumerate(segs):
        exon = exons[n]
        seg_start, seg_end = seg.split('-')
        exon_start, exon_end = exon.split(':')

        seq_pos_range = range(seq_pos[0], seq_pos[1]+1)
        seg_pos_range = range(int(seg_start), int(seg_end)+1)
        exon_pos_range = range(int(exon_start), int(exon_end)+1)

        inters = (set(seg_pos_range) & set(seq_pos_range))
        map_pos = dict(zip(seg_pos_range, exon_pos_range))
        if inters != set():
            pos_s = min(inters)
            pos_e = max(inters)
            range_s = map_pos.get(pos_s)
            range_e = map_pos.get(pos_e)
            regions.append(f"{pos_s}:{pos_e}|{range_s}:{range_e}|{exon}")
    return regions


def parse_rnahybrid(fname, trans_info, outdir=None,
                    continuous=False, u2t=False):
    """
    fname: RNAhybrid_pattern.txt
    """
    print(f"process RNAhybrid result: {fname}")
    if outdir is None:
        outdir = os.path.split(fname)[0]

    fn = os.path.split(fname)[1]
    stem, suffix = os.path.splitext(fn)
    if fname.endswith('.gz'):
        stem, suffix = os.path.splitext(stem)
    output = os.path.join(outdir, f"{stem}-parse{suffix}")

    if not os.path.exists(output):
        if fname.endswith('.gz'):
            fo = gzip.open(fname, 'rb')
            cont = fo.read().decode("gbk")
        else:
            fo = open(fname, 'r')
            cont = fo.read()

        aligns = re.findall(f"target: .*?5'\n\n\n", cont, re.S)

        n = 0
        rows = []
        for al in aligns:
            position = int(get_position(al))
            target_length, mirna_seq_length = get_length(al)
            trans = get_trans_name(al)
            mirna = get_mirna_name(al)
            pvalue = get_pvalue(al)
            energy = get_mfe(al)
            target_seq, mirna_seq = get_merged_seq(al)

            ttag = trans
            if '_' in ttag:
                ttag = trans.split('_')[0]

            dft = [""]*4
            chrom = trans_info.get(ttag, dft)[1]
            strand = trans_info.get(ttag, dft)[2]
            exons = trans_info.get(ttag, dft)[3].split('|')
            gene = trans_info.get(ttag, dft)[0]

            if continuous:
                target_seq = target_seq.replace('-', '')
                mirna_seq = mirna_seq.replace('-', '')
            if u2t:
                target_seq = target_seq.replace('U', 'T')
                mirna_seq = mirna_seq.replace('U', 'T')

            target_seq_length = len(target_seq.replace('-', ''))
            target_start = position
            target_end = position + target_seq_length - 1

            seq_pos = [target_start, target_end]
            if strand == '-':
                exons = list(reversed(exons))
                seq_pos = get_seq_pos(target_start, int(target_length), target_seq_length)
            segs = get_segs(exons)
            regions = get_chr_regions(seq_pos, segs, exons)

            region = 1
            for rg in regions:
                p_range, r_range, exon = rg.split('|')
                chr_start, chr_end = r_range.split(':')
                row = [chrom, chr_start, chr_end, gene, strand,
                       mirna, trans, position,
                       target_length, pvalue, energy,
                       ''.join(reversed(mirna_seq)), target_seq,
                       mirna_seq_length, target_seq_length,
                       region,
                       ]
                region += 1
                rows.append(row)
            n += 1

        fo.close()

        columns = ["chr", "start", "end", "gene", "strand",
                   "miRNA", "transcript", "position",
                   "target_length", "pvalue", "energy",
                   "miRNA_seq(5'-3')", "target_seq(5'-3')",
                   "miRNA_seq_length", "target_seq_length",
                   "region"]
        df = pd.DataFrame(rows, columns=columns)
        df.to_csv(output, sep="\t", index=False)
        check_outputf(output)
    return output


def parse_pattern(pattern_file, gtf, out_dir, continuous=True,
                  u2t=True):
    """
    pattern_file: RNAhybrid_hsa_Pattern.txt
    gtf: Homo_sapiens.GRCh38.101.gtf
    out_dir：output directory
    continuous：remove `-` in seq
    u2t: U -> T
    out_bed: output a bed format file
    """

    fn = os.path.split(pattern_file)[1]
    stem, suffix = os.path.splitext(fn)
    if pattern_file.endswith('.gz'):
        stem, suffix = os.path.splitext(stem)
    output = os.path.join(out_dir, f"{stem}-parse{suffix}")

    if not os.path.exists(output):
        trans_info = parse_gtf(gtf)

        parse_rnahybrid(pattern_file, out_dir, trans_info,
                        continuous=continuous, u2t=u2t)


if __name__ == "__main__":
    pattern_file_ = os.path.abspath(sys.argv[1])
    gtf_ = os.path.abspath(sys.argv[2])
    out_dir_ = os.path.abspath(sys.argv[3])
    parse_pattern(pattern_file=pattern_file_,
                  gtf=gtf_,
                  out_dir=out_dir_,
                  )
