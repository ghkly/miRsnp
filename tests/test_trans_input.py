#!/usr/bin/python
# _*_ coding: utf-8 _*_

from mirsnp.trans_input import trans_input


def test():
    fin = '../data/variants_input.xlsx'
    gtf = r'../data/Homo_sapiens.GRCh38.101-example.gtf'
    trans_input(fin, '../data/variants.txt', gtf)
