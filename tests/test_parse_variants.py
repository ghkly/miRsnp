#!/usr/bin/python
# _*_ coding: utf-8 _*_
# @Time     : 2021/12/20 0020 11:36
# @Author   : yangliuke

from mirsnp.mkaltfa import parse_variants

# {'ENST00000171111': [('rs45524632|19:10486312:C:A',), ('rs45524632|19:10486312:C:G',)],
# 'ENST00000202917': [('rs1051042|12:112919432:G:C',), ('rs1131476|12:112919404:G:A',),
# ('rs1131476|12:112919404:G:C',), ('rs1131476|12:112919404:G:T',), ('rs2660|12:112919637:G:A',)], ...}

variants_file = '../data/variants.txt'


def test():
    print(parse_variants(variants_file))

