#!/usr/bin/python3
# _*_ coding: utf-8 _*_
# @Time     : 2022/1/31 19:21

import re
import os
import pandas as pd

from mirsnp.parse_rnahybrid import parse_rnahybrid, get_position, get_mirna_name, get_trans_name
from mirsnp.utils import check_outputf, DEFAULT_ENERGY, INTERVAL


def get_pt_map(pt):
    pt_map = {}
    with open(pt) as fr:
        aligns = re.findall(f"target: .*?5'\n\n\n", fr.read(), re.S)
        for al in aligns:
            position = get_position(al)
            trans = get_trans_name(al)
            mirna = get_mirna_name(al)
            flag = f"{trans}:{mirna}:{position}"
            pt_map[flag] = al.replace('\n', '#')
    return pt_map


def get_target_pt(flag, pt_map):
    fs = flag.split(':')
    k = f"{fs[0]}:{fs[1]}:{fs[-1]}"
    return pt_map.get(k, "")


def parse_common_variants(variants_file):
    df = pd.read_csv(variants_file, sep='\t')
    gene_map = dict(zip(df['transcript'], df['gene']))
    strand_map = dict(zip(df['transcript'], df['strand']))
    return gene_map, strand_map


def isin_interval(tpos, seq_pos):
    tpos = int(tpos)
    start, end = seq_pos.split(' ')[1].split('-')
    if int(start) <= tpos <= int(end):
        return 'Yes'
    return 'No'


def get_target_pos(mirna53, target53, flag, interval):
    pos = int(flag.split(':')[-1])
    tpos = int(flag.split(':')[0].split('_')[6])
    target35 = ''.join(reversed(target53))

    pos_range = []
    m = 0
    for i in target53:
        if i != '-':
            pos_range.append(pos+m)
            m += 1
        else:
            pos_range.append(-1)
    rev_range = list(reversed(pos_range))

    spos = 1
    mirna_seqs = []
    target_seqs = []
    target_poss = []
    for n, b in enumerate(mirna53):
        tb = target35[n]
        if b != '-':
            if interval[0] <= spos <= interval[1]:
                mirna_seqs.append(b)
                target_seqs.append(tb)
                target_poss.append(rev_range[n])
            spos += 1
        else:
            if interval[0] <= spos <= interval[1]:
                target_seqs.append(tb)
                target_poss.append(rev_range[n])

    mirna_seq = ''.join(mirna_seqs)
    target_seq = ''.join(reversed(target_seqs)).replace('-', '')

    fpos = list(filter(lambda x: x != -1, target_poss))
    target_poss = list(reversed(fpos))

    target_seq_pos = f"{target_seq} {target_poss[0]}-{target_poss[-1]}"
    intgt = isin_interval(tpos, target_seq_pos)
    # "{}|{}|{}".format(mirna_seq, target_seq_pos, intgt)
    return intgt


def fmt_pt(pt, tag, trans_info, interval=INTERVAL):
    pp = parse_rnahybrid(pt, trans_info)
    pt_map = get_pt_map(pt)
    df = pd.read_csv(pp, sep='\t', low_memory=False)
    if not df.empty:
        if df.shape[0] >= 1 and df.shape[1] > 3:
            df = df[df['region'] == 1].drop_duplicates()
            df['flag'] = df["transcript"] + ":" + df["miRNA"].astype(str) + \
                          ":" + df['gene'].astype(str) + \
                          ":" + df['strand'].astype(str) + \
                          ":" + df['position'].astype(str)
            df = df[["flag", "pvalue", "energy", "miRNA_seq(5'-3')", "target_seq(5'-3')"]]

            # itv = f"{interval[0]}-{interval[1]}"
            df['variant_in_2-8'] = df.apply(lambda row: get_target_pos(
                row["miRNA_seq(5'-3')"], row["target_seq(5'-3')"], row["flag"], interval), axis=1)
            df['pattern'] = df['flag'].apply(lambda x: get_target_pt(x, pt_map))
            df = df.set_index('flag', drop=True)
            df = df.add_prefix('{}_'.format(tag))
            return df


def cal_dist(flags):
    # ENST00000551241_rs1859333_12_112933161_T_C_36:1364:OAS1:+:29
    res = []
    for flag in flags:
        ar = flag.split(':')
        variant_pos = ar[0].split('_')[-1]
        tpos = ar[-1]
        if variant_pos.isdigit() and tpos.isdigit():
            res.append(str(int(variant_pos) - int(tpos)))
        else:
            res.append('NA')
    return res


def get_diff(x, y):
    if "NA" not in [x, y]:
        diff = round(float(y) - float(x), 2)
        if diff < 0:
            return [diff, "More stable"]
        elif diff > 0:
            return [diff, "Less stable"]
        else:
            return [diff, "No change"]
    return ['NA', 'NA']


def get_energy_classify(x, y):
    if "NA" in [x, y]:
        if x == "NA":
            return "new"
        elif y == "NA":
            return "off-target"
    return 'consistent'


def classify_and_sort(df, default_energy=DEFAULT_ENERGY):
    """根据能量值差异分类"""
    if not df.empty:
        df = df.copy()
        df['classify'] = df.apply(lambda row: get_energy_classify(row["ref_target_seq(5'-3')"], row["alt_target_seq(5'-3')"]), axis=1)

        df.loc[df.classify == 'new', 'ref_energy'] = default_energy
        df.loc[df.classify == 'off-target', 'alt_energy'] = default_energy

        df['diff_energy'] = df.apply(lambda row: get_diff(row['ref_energy'], row['alt_energy'])[0], axis=1)
        df['state'] = df.apply(lambda row: get_diff(row['ref_energy'], row['alt_energy'])[1], axis=1)
        return df


def merge_refalt(refpt, altpt, outdir, trans_info,
                 name='f', td='raw', interval=INTERVAL):
    df_ref = fmt_pt(refpt, 'ref', trans_info, interval=interval)
    df_alt = fmt_pt(altpt, 'alt', trans_info, interval=interval)

    output_all = os.path.join(outdir, f"energy_stats-{name}-result.csv")
    dfm = None
    if df_ref is not None and df_alt is not None:
        dfm = pd.merge(df_ref, df_alt, left_index=True, right_index=True, how='outer').fillna('NA')
        dfm = classify_and_sort(dfm)
        dfm['transcript_direction'] = td

        dfm = dfm[dfm['state'] != 'No change']
        dfm = dfm.drop_duplicates()
        dfm['target_position_to_variant_position'] = cal_dist(dfm.index)
        dfm.to_csv(output_all)
        check_outputf(output_all)
    else:
        print('either {} or {} is null, no output!'.format(refpt, altpt))
    return dfm

