#!/usr/bin/python
# _*_ coding: utf-8 _*_
# @Time     : 2021/12/3 0003 11:15
# @Author   : yangliuke

import os
import sys
import pandas as pd


def fmt_anno_tbl(fpath):
    res = []
    with open(fpath, 'r') as fo:
        for line in fo:
            lst = line.strip('\n').split('\t')
            c, s, e, ref, alts, *o_ = lst
            for alt in alts.strip('"').split(","):
                rl, al = len(ref), len(alt)
                if rl != al:
                    if rl > al:
                        ref = ref[1:]
                        alt = '-'
                        s = int(s) + 1
                        e = s + len(ref) - 1
                    else:
                        ref = '-'
                        alt = alt[1:]

                new_lst = [c, int(s), int(e), ref, alt, *o_]
                res.append(new_lst)

    with open(fpath, 'w') as fw:
        for lne in res:
            line = [str(i) for i in lne]
            fw.write('\t'.join(line) + '\n')


def trans_fmt(fp):
    prefix, suffix = os.path.splitext(fp)
    if suffix.lower() == '.xlsx':
        df = pd.read_excel(fp, dtype=str)
    elif suffix.lower() == '.csv':
        df = pd.read_csv(fp, dtype=str)
    else:
        df = pd.read_csv(fp, dtype=str, sep='\t')
    df = df.dropna(subset=['chr', 'pos'])
    return df


def run_annovar(fpath, humandb):
    df = trans_fmt(fpath)

    prefix, _ = os.path.splitext(fpath)
    output = prefix + '-annovar' + '.tsv'

    raw_cols = df.columns.tolist()

    vdf = pd.DataFrame()
    vdf['#CHROM'] = 'chr' + df['chr']
    vdf['START'] = df['pos'].astype(int)
    vdf['END'] = df['pos'].astype(int)
    vdf['REF'] = df['ref']
    vdf['ALT'] = df['alt']

    raw_cols.remove('chr')
    raw_cols.remove('pos')
    raw_cols.remove('ref')
    raw_cols.remove('alt')

    for i in raw_cols:
        vdf[i] = df[i]

    vdf['ref_len'] = vdf['REF'].apply(lambda x: len(x)-1)
    vdf['END'] = vdf['START'] + vdf['ref_len']
    vdf = vdf.drop(['ref_len'], axis=1)

    vdf = vdf.drop_duplicates().sort_values(['#CHROM', 'START'])
    vdf.to_csv(output, sep='\t', index=False, header=0)
    print(f'output：{output}')
    fmt_anno_tbl(output)

    fn = os.path.basename(output)
    fb = os.path.splitext(fn)[0]
    annovardir = './annotate'
    if not os.path.exists(annovardir):
        os.mkdir(annovardir)
    anno_sh = prefix + '-annovar' + '.sh'
    with open(anno_sh, 'w') as fw:
        fw.write(f"""#!/usr/bin/bash

table_annovar.pl \\
    ./{fn} \\
    {humandb} \\
    --buildver hg38 \\
    --protocol refGene,avsnp150,1000g2015aug_eur,1000g2015aug_afr,1000g2015aug_amr,1000g2015aug_eas,1000g2015aug_sas \\
    --operation g,f,f,f,f,f,f \\
    --nastring . \\
    --thread 2 \\
    --remove \\
    --csvout \\
    --polish \\
    --otherinfo \\
    --outfile {annovardir}/{fb}
            """)
    print(f'输出文件：{anno_sh}\n请执行 bash {anno_sh}')
    return anno_sh


if __name__ == "__main__":
    fpath_ = sys.argv[1]
    humandb_ = sys.argv[2]

    sh = run_annovar(fpath_, humandb_)
    os.system('bash {}'.format(sh))

    # fin = '../example/variants_input.xlsx'
    # run_annovar(fin, humandb_)

