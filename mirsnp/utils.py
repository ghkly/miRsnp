#!/usr/bin/python3
# _*_ coding: utf-8 _*_

import os
import re
import time
import subprocess
from numpy import reshape

from mirsnp.config import *


def group_fasta(fa):
    """
    return:
        {>id: ATCG..., }
    """
    ids = []
    seqs = []
    seq = ''
    with open(fa, 'r') as fo:
        n = 0
        while True:
            line = fo.readline().strip('\n')

            if line.startswith('>'):
                ids.append(line)
                if seq:
                    seqs.append(seq)
                seq = ''
            else:
                if line:
                    seq += line

            if line == '':
                seqs.append(seq)
                break

            n += 1
    seqmap = dict(zip(ids, seqs))
    return seqmap


def format_seq(seqs, cols=1, nums=10, sep='\n'):
    """
    按照指定间隔长度输出字符串(连续的单个字符串):
    hfjksdffhieowfkdsjfs -> hfjksdff hieowfkd sjfs
    cols：每行列数
    nums：每列连续字符串长度
    """
    if cols <= 0 or nums <= 0:
        print('cols and nums must be integer(>0)')
        return seqs

    length = len(seqs)
    heads = seqs[:(length - length % nums)]
    tails = seqs[(length - length % nums):]
    seqs_list = re.findall(r'.{%s}' % nums, heads)
    if tails:
        seqs_list.append(tails)

    if cols == 1:
        formatseq = seqs_list
    elif cols > 1:
        list_length = len(seqs_list)
        tail_length = list_length % cols
        if tail_length == 0:
            seqmatrix = reshape(seqs_list, [list_length // cols, cols])
            formatseq = [' '.join(lst) for lst in seqmatrix]
        else:
            head_list = seqs_list[:(list_length - tail_length)]
            tail_list = seqs_list[(list_length - tail_length):]
            seqmatrix = reshape(head_list, [len(head_list) // cols, cols])
            formatseq = [' '.join(lst) for lst in seqmatrix]
            formatseq.append(' '.join(tail_list))
    else:
        formatseq = [seqs]
    return sep.join(formatseq)


def base_match(seq, reverse=False, compl=False):
    """获取DNA序列的反向/互补序列"""
    bases = list(seq)
    if reverse:
        if compl:
            reverse_seq = ''.join(reversed(bases))
            return ''.join([BASES.get(i, i) for i in list(reverse_seq)])
        else:
            return ''.join(reversed(bases))
    else:
        if compl:
            return ''.join([BASES.get(i, i) for i in bases])
        else:
            return seq


def rev_match(seq):
    return base_match(seq, reverse=True, compl=True)


def split_row(df, col_name, sep=','):
    tmp_col = f"{col_name}_tmp"
    df[tmp_col] = df[col_name]
    df = df.drop(tmp_col, axis=1).join(df[tmp_col].str.split(
        sep, expand=True).stack().reset_index(
        level=1, drop=True).rename(tmp_col))
    df[col_name] = df[tmp_col]
    df = df.drop([tmp_col], axis=1)
    return df


def strip_keys(dt):
    return {k.split(' ')[0]: v for k, v in dt.items()}


def check_outputf(outputf):
    if os.path.exists(outputf):
        print("output：{}".format(outputf))
    else:
        print("output failed：{}".format(outputf))


def run_subshell(cmd):
    try:
        s = subprocess.Popen(cmd, stderr=subprocess.PIPE,
                             stdout=subprocess.PIPE,
                             shell=True, preexec_fn=os.setpgrp)
        s.wait()
        time.sleep(1)
    except Exception as e:
        print(str(e))
