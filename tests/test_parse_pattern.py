#!/usr/bin/python3
# _*_ coding: utf-8 _*_


from mirsnp.parse_rnahybrid import parse_pattern

ref_pt = '../data/ref-flank35bp-raw-pattern.txt'
alt_pt = '../data/alt-flank35bp-raw-pattern.txt'
gtf = '../data/Homo_sapiens.GRCh38.101-example.gtf'


def test():
    parse_pattern(ref_pt, gtf, '../data', continuous=False, u2t=False)
    parse_pattern(alt_pt, gtf, '../data', continuous=False, u2t=False)




