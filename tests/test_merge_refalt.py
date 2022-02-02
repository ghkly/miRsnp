#!/usr/bin/python3
# _*_ coding: utf-8 _*_
# @Time     : 2022/1/31 11:11

from mirsnp.formatpt import merge_refalt, fmt_pt
from mirsnp.parse_gtf import parse_gtf

ref_pt = '../data/ref-flank35bp-raw-pattern.txt'
alt_pt = '../data/alt-flank35bp-raw-pattern.txt'
gtf = '../data/Homo_sapiens.GRCh38.101-example.gtf'


def test():
    trans_info = parse_gtf(gtf)
    merge_refalt(ref_pt, alt_pt, './data', trans_info)

