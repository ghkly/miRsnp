#!/usr/bin/python3
# _*_ coding: utf-8 _*_

import os
import pandas as pd

from mirsnp.utils import check_outputf
from mirsnp.parse_gtf import parse_gtf


def get_trans_info(gtf):
    info = parse_gtf(gtf)
    res = {}
    for trans, data in info.items():
        gene = data[0]
        strand = data[2]
        exons = data[3]
        if strand == '-':
            exons = '|'.join(list(reversed(exons.split('|'))))
        if gene not in res:
            res[gene] = {}
        tag = '{}:{}'.format(trans, strand)
        res[gene][tag] = exons
    return res


def in_seg(ps, seg):
    st, ed = seg.split(':')
    st = int(st)
    ed = int(ed)
    if st > ed:
        st, ed = ed, st
    if int(st) < int(ps) < int(ed):
        return True
    return False


def in_exons(ps, exons):
    segs = exons.split('|')
    sts = [in_seg(ps, seg) for seg in segs]
    return any(sts)


def get_trans(gene, pos, trans_info):
    res = []
    mp = trans_info.get(gene, {})
    for tag, exons in mp.items():
        trans, strand = tag.split(':')
        if in_exons(pos, exons):
            res.append((trans, strand))
    return res


def trans_input(fp, output, gtf):
    suffix = os.path.splitext(fp)[1].lower()
    if suffix == '.xlsx':
        df = pd.read_excel(fp, dtype=str)
    elif suffix == '.csv':
        df = pd.read_csv(fp, dtype=str)
    else:
        df = pd.read_csv(fp, dtype=str, sep='\t')
    df = df.dropna(subset=['chr', 'pos'], axis=0)

    trans_info = get_trans_info(gtf)
    fw = open(output, 'w', encoding='utf-8')
    header = ['rs_snp', 'rstag', 'gene', 'strand', 'miRNA', 'transcript']
    fw.write('\t'.join(header) + '\n')

    for idx, row in df.to_dict(orient='index').items():
        rsid = row['rsid']
        ref = row['ref']
        alt = row['alt']
        chrom = row['chr']
        pos = row['pos']
        gene = row.get('gene', '').strip()
        rs_snp = '{}|{}:{}:{}:{}'.format(rsid, chrom, pos, ref, alt)
        rstag = '{}_{},{}'.format(rsid, ref, alt)
        mirna = 'none'

        tsinfo = get_trans(gene, pos, trans_info)
        if tsinfo:
            for transcript, strand in tsinfo:
                line = [rs_snp, rstag, gene, strand, mirna, transcript]
                fw.write('\t'.join(line) + '\n')
    fw.close()
    check_outputf(output)
    return output
