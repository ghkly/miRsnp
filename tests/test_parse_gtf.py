#!/usr/bin/python3
# _*_ coding: utf-8 _*_
# @Time     : 2022/1/30 18:44


from mirsnp.parse_gtf import parse_gtf
gtf = '../data/Homo_sapiens.GRCh38.101-example.gtf'


def test():
    print(parse_gtf(gtf))

