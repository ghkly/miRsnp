#!/usr/bin/python3
# _*_ coding: utf-8 _*_

# compare energy between ref and alt
import os
import re
import pandas as pd

from mirsnp.utils import (base_match, rev_match, group_fasta,
                          format_seq, check_outputf)
from mirsnp.config import _fa_groups, _direct_groups, BASES, DEFAULT_FLANK


def fetch_transcript_fasta(variants_file, total_fa, out_fa):
    df = pd.read_csv(variants_file, sep='\t')
    target_trans = df['transcript'].unique().tolist()
    print('total transcripts: {}'.format(len(target_trans)))
    if not os.path.exists(out_fa):
        group_fa = group_fasta(total_fa)
        with open(out_fa, "w") as rfw:
            n = 1
            for k, v in group_fa.items():
                ktrans = re.search(r">(.*?)[ |\t]", k, re.S).group(1)
                if ktrans in target_trans:
                    # print('fetch ref transcript sequence {}: {}'.format(n, ktrans))
                    rfw.write("{}\n".format(k))
                    rfw.write("{}\n".format(format_seq(v, 1, 70)))
                    n += 1
    return out_fa


def flatten_list(arr):
    res = []
    for i in arr:
        if len(i) >= 1:
            res.extend([[j] for j in i])
        else:
            res.append([i])
    return res


def parse_variants(fpath):
    print('parse variants: ')
    res = {}
    with open(fpath) as fr:
        fr.readline()
        for line in fr:
            lst = line.strip().split('\t')
            transcript = lst[-1]
            if transcript not in res:
                res[transcript] = []
            res[transcript].append((lst[0], ))
    return res


def mod_segs(exons):
    """生成负链segs"""
    dist = [0]
    for ex in exons:
        s, e = ex.split('-')
        exlen = int(e) - int(s) + 1
        dist.append(exlen)

    segs = []
    for i in range(len(dist) - 1):
        segs.append([sum(dist[: i + 1]) + 1, sum(dist[:i + 2])])

    segs = [f"{i[0]}-{i[1]}" for i in segs]

    return segs


def split_transcript_segs(fa):
    """
    fa: fasta

    return:
        {"1": {transcript1: {"1-10": ATCG...}, transcript2: {}...}}
    """
    grouped_fa = group_fasta(fa)
    transcript_segs = {}
    tpos_map = {}
    exons_map = {}
    print('split transcript segs: ')
    # 转录本计数位置
    for idx, seq in grouped_fa.items():
        tpos_start = 1
        transcript = re.search(r">(.*?) ", idx, re.S).group(1)
        strand = re.search(r"-\d+\|(.*?) exons", idx, re.S).group(1)
        chrom = re.search(r"loc:(\w+)\|", idx, re.S).group(1)
        exons = re.search(r"exons:(.*?) ", idx, re.S).group(1).split(',')
        segs = re.search(" segs:(.*?)$", idx, re.S).group(1).split(',')

        transcript = f"{transcript}:{strand}"

        if chrom not in transcript_segs:
            transcript_segs[chrom] = {}

        if strand == '+':
            pos_map = dict(zip(segs, exons))
        else:
            # 负链，首先序列反向互补，然后根据染色体位置生成新的segs
            seq = base_match(seq, reverse=True, compl=True)
            pos_map = dict(zip(mod_segs(exons), exons))

        rg_map = {}
        spos = {}
        for pos, region in pos_map.items():
            pos_start, pos_end = pos.split('-')
            pos_start, pos_end = int(pos_start), int(pos_end)
            seg_seq = seq[pos_start-1: pos_end]
            rg_map[region] = seg_seq

            tpos_end = tpos_start + len(seg_seq) - 1
            spos[region] = f"{tpos_start}-{tpos_end}"
            tpos_start = tpos_end + 1

        transcript_segs[chrom][transcript] = rg_map
        tpos_map[transcript] = spos

        exons2 = re.search("exons:(.*?) segs:", idx).group(1)
        exons_map[transcript] = exons2

    return exons_map, transcript_segs, tpos_map


def fmt_variant_name(v):
    """
    v:
        ('rs1045288|11:237087:A:G',
        'rs1128320|11:244167:C:T',)
    """
    res = ''
    for i in v:
        rs = i.split('|')[0]
        ref = i.split(':')[-2]
        alt = i.split(':')[-1]
        res = res + f"_{rs}({ref},{alt})"
    return res


def fmt_identify_name(v):
    """
    v:
        ('rs1045288|11:237087:A:G',
        'rs1128320|11:244167:C:T',)
    """
    return ','.join(v)


def fmt_identify_name1(v):
    """
    v:
        ('rs1045288|11:237087:A:G',
        'rs1128320|11:244167:C:T',)
    """
    v = [i.replace('|', '_').replace(':', '_') for i in v]
    return ','.join(v)


def get_abs_tpos(chrom_pos, t_pos_map):
    """
    chrom_pos: 43225635
    t_pos_map： {'43225635-43226143': '1-509',
    '43226327-43226395': '510-578',
    '43229261-43229478': '579-796'}
    """
    for k, v in t_pos_map.items():
        c = int(chrom_pos)
        s, e = k.split('-')
        ps, pe = v.split('-')
        if int(s) <= c <= int(e):
            dist = c - int(s)
            return int(ps) + dist
        else:
            continue


def get_overlap(region, variants_pos, chrom_variants, t_pos_map):
    start, end = [int(i) for i in region.split('-')]
    pos = variants_pos & set(range(start, end+1))
    res = []
    for i in reversed(sorted(list(pos))):
        tpos = get_abs_tpos(i, t_pos_map)
        res.append([i, chrom_variants.get(i), tpos])
    return res


def get_seg_map(seg_range, seg_seq):
    """
    return:
        染色体位置对应序列计数位置
        {123456: 1, 123457: 2, ...}
    """
    start, end = [int(i) for i in seg_range.split('-')]
    seg_map = dict(zip(range(start, end+1), range(0, len(seg_seq))))
    return seg_map


def map_trans(df):
    df['transcript'] = df['flag'].apply(lambda x: x.split(':')[0].split('_')[0])
    return df


def map_mirna_name(df):
    df['miRNA'] = df['flag'].apply(lambda x: x.split(':')[1])
    return df


def map_gene(df, gene_map):
    df['gene'] = df['transcript'].apply(lambda x: gene_map.get(x.split(':')[0], 'NA'))
    return df


def map_strand(df, strand_map):
    df['strand'] = df['transcript'].apply(lambda x: strand_map.get(x.split(':')[0], 'NA'))
    return df


def get_log_info(fname):
    res = {}
    with open(fname, 'r', encoding='utf-8') as fo:
        for line in fo:
            if 'not_replaced' not in line:
                arr = line.strip().split(' ')
                trans = arr[1].split(':')[2]
                strand = arr[1].split(':')[-1]
                rs_pos = arr[3].split(':')[-1]
                rs_tag = arr[2].replace('|', '_').replace(':', '_')
                variant_tag = f"{trans}_{rs_tag}"
                res[variant_tag] = [rs_pos, strand]
    return res


def get_exons(fa, t):
    with open(fa, 'r') as fo:
        cont = fo.read()
    t1 = re.search(f"{t}(.*?)\n", cont, re.S).group(1)
    exons = re.search("exons:(.*?) segs:", t1).group(1)
    return exons


def sum_segs(exons):
    res = []
    for i in exons.split(','):
        start, end = i.split('-')
        res.append(int(end) - int(start) + 1)
    return sum(res)


def get_seq_pos(start, target_length, seq_length):
    end_pos = start + seq_length - 1
    tail_dist = target_length - end_pos

    rev_start = tail_dist + 1
    rev_end = rev_start + seq_length - 1

    return [rev_start, rev_end]


def sub_dna_seq(chrom, transcript, seg_seq, overlap_variants, seg_map, flog):
    seq = list(seg_seq)
    pos_list = []
    for i in overlap_variants:
        chrom_pos, ref_alt_list, tpos = i
        ref = ref_alt_list[0]
        alt = ref_alt_list[1]
        rs = ref_alt_list[2]
        chrom_pos = int(chrom_pos)
        ref_len = len(ref)
        alt_len = len(alt)

        seg_start = seg_map.get(chrom_pos)
        last_chrom_pos = list(seg_map)[-1]
        last_seg_pos = seg_map.get(last_chrom_pos)

        if chrom_pos == last_chrom_pos:
            if ref_len == 1:
                seg_end = last_seg_pos + 1
            else:
                seg_end = None
        else:
            seg_end = seg_map.get(chrom_pos + len(ref))

        fasta_ref = ''.join(seg_seq[seg_start: seg_end])
        info = f"transcript:{transcript} {rs}|{chrom}:{chrom_pos}:{ref}:{alt} pos+:{tpos} ref_fasta_allele:{fasta_ref} allele:{ref} -> {alt}"
        pos_list.append(tpos)

        if seg_end is not None:
            if ref_len == alt_len:
                if ref_len == 1:
                    if ref == fasta_ref:
                        flog.write(f"replace snp:{info}\n")
                        seq[seg_start] = alt
                    elif ref == BASES.get(fasta_ref):
                        flog.write(f"replace snp(reverse complementary):{info}\n")
                        seq[seg_start] = BASES.get(alt)
                    else:
                        flog.write(f"replace snp error:{info}, not_replaced\n")
                else:
                    if ref == fasta_ref:
                        flog.write(f"replace complex:{info}\n")
                        seq[seg_start: seg_end] = alt
                    elif ref == rev_match(fasta_ref):
                        flog.write(f"replace complex(reverse complementary):{info}\n")
                        seq[seg_start: seg_end] = rev_match(alt)
                    else:
                        flog.write(f"replace complex error:{info}, not_replaced\n")

            elif ref_len > alt_len:
                diff = ref_len - alt_len
                cmp = "-"*diff
                if ref == fasta_ref:
                    flog.write(f"replace del:{info}\n")
                    seq[seg_start: seg_end] = list(alt + cmp)

                elif ref == rev_match(fasta_ref):
                    flog.write(f"replace del(reverse complementary):{info}\n")
                    seq[seg_start: seg_end] = list(rev_match(alt)+cmp)
                else:
                    flog.write(f"replace del error:{info}, not_replaced\n")

            elif ref_len < alt_len:
                if ref_len == 1:
                    if ref == fasta_ref:
                        flog.write(f"replace ins:{info}\n")
                        seq[seg_start] = alt
                    elif ref == BASES.get(fasta_ref):
                        flog.write(f"replace ins(reverse complementary):{info}\n")
                        seq[seg_start] = rev_match(alt)
                    else:
                        flog.write(f"replace ins error:{info}, not_replaced\n")
                else:
                    seq[seg_start: seg_end] = ['-']*ref_len
                    if ref == fasta_ref:
                        flog.write(f"replace ins:{info}\n")
                        seq[seg_start] = alt
                    elif ref == rev_match(fasta_ref):
                        flog.write(f"replace ins(reverse complementary):{info}\n")
                        seq[seg_start] = rev_match(alt)
                    else:
                        flog.write(f"replace ins error:{info}, not_replaced\n")
            else:
                pass
        else:
            flog.write(f"warning, {info} out of range of reference sequence, end location of transcript:{last_chrom_pos} not_replaced！\n")

    return ''.join(seq).replace('-', ''), pos_list


def sub_variant_bases(transcript_segs, tpos_map, variants, outdir, exons_map, flank=DEFAULT_FLANK):
    alt_log = os.path.join(outdir, "alt.log")
    flog = open(alt_log, 'w', encoding='utf-8')

    fw_mp = {}
    for gp in _fa_groups:
        if gp not in fw_mp:
            fw_mp[gp] = {}
        for tp in _direct_groups:
            out_fa = os.path.join(outdir, f"{gp}-flank{flank}bp-{tp}.fa")
            fw_mp[gp][tp] = open(out_fa, 'w', encoding='utf-8')

    print('sub variant bases: ')
    existed_trans = []

    for chrom, chrom_transcript_segs in transcript_segs.items():
        for transcript, segs in chrom_transcript_segs.items():
            # transcript: ENST00000496421:+
            transname, strand = transcript.split(':')
            t_variants = variants.get(transname, [])
            # t_variants: [('rs773840258|1:113900018:CCACTTTCT:C',), ('rs781125916|1:113900022:T:C',), ...]
            t_pos_map = tpos_map.get(transcript)

            if t_variants:
                for variant in t_variants:
                    # print(f"process: chr:{chrom} {transcript} {variant}")
                    chrom_positions = set([int(i.split(':')[1]) for i in variant])
                    chrom_variants = {int(i.split(':')[1]): [i.split(':')[2],
                                                             i.split(':')[3],
                                                             i.split('|')[0]] for i in variant}
                    alt_seq_list = []
                    ref_seq_list = []
                    pos_list_all = []
                    for seg_range, seg_seq in segs.items():
                        seg_map = get_seg_map(seg_range, seg_seq)
                        overlap_variants = get_overlap(seg_range, chrom_positions, chrom_variants, t_pos_map)
                        if overlap_variants:
                            new_seg_seq, pos_list = sub_dna_seq(chrom, transcript, seg_seq, overlap_variants, seg_map, flog)
                            pos_list_all.extend(pos_list)
                        else:
                            new_seg_seq = seg_seq
                        alt_seq_list.append(new_seg_seq)
                        ref_seq_list.append(seg_seq)

                    alt_seq = ''.join(alt_seq_list)
                    ref_seq = ''.join(ref_seq_list)
                    if strand == '-':
                        alt_seq = rev_match(alt_seq)
                        ref_seq = rev_match(ref_seq)

                    # get flank seq
                    if len(pos_list_all) == 1:
                        tpos = int(pos_list_all[0])
                        exons = exons_map.get(transcript)
                        if exons is not None:
                            target_len = sum_segs(exons)
                            if strand != '+':
                                tpos = get_seq_pos(tpos, target_len, 1)[0]

                            # pos_range = [1, tpos + flank]
                            # if tpos > flank:
                            #     pos_range = [tpos - flank, tpos + flank]

                            if tpos <= flank:
                                pos_range = [1, tpos + flank + (flank - tpos + 1)]
                            else:
                                tail_dist = target_len - tpos
                                if tail_dist >= flank:
                                    pos_range = [tpos - flank, tpos + flank]
                                else:
                                    pos_range = [tpos - flank - (flank - tail_dist), tpos + flank]
                            variant_pos = tpos - pos_range[0] + 1
                            ref_raw = ref_seq[pos_range[0] - 1: pos_range[1]]
                            alt_raw = alt_seq[pos_range[0] - 1: pos_range[1]]

                            if ref_raw != '' and alt_raw != '':
                                ref_rev = base_match(ref_raw, reverse=True, compl=False)
                                ref_cmp = base_match(ref_raw, reverse=False, compl=True)
                                ref_rmp = base_match(ref_raw, reverse=True, compl=True)

                                alt_rev = base_match(alt_raw, reverse=True, compl=False)
                                alt_cmp = base_match(alt_raw, reverse=False, compl=True)
                                alt_rmp = base_match(alt_raw, reverse=True, compl=True)

                                for gp in _fa_groups:
                                    fa_idn = f'>{transname}_{variant[0].replace("|", "_").replace(":", "_")}_{variant_pos}'
                                    # >ENST00000445409_rs2660_12_112919637_G_A_36
                                    for tp in _direct_groups:
                                        fw_mp[gp][tp].write(f'{fa_idn}\n' + eval(f'{gp}_{tp}') + '\n')

                                existed_trans.append(transname)
    flog.close()
    for gp in _fa_groups:
        for tp in _direct_groups:
            fw_mp[gp][tp].close()


def run_sub_variant_bases(variant_file, fasta, outdir, flank=DEFAULT_FLANK):
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    variants = parse_variants(variant_file)
    exons_map, transcript_segs, tpos_map = split_transcript_segs(fasta)
    sub_variant_bases(transcript_segs, tpos_map, variants, outdir, exons_map, flank=flank)

