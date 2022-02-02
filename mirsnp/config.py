#!/usr/bin/python3
# _*_ coding: utf-8 _*_
# @Time     : 2022/2/2 17:17


BASES = {'A': 'T', 'C': 'G', 'T': 'A', 'G': 'C'}
DEFAULT_FLANK = 35
DEFAULT_ENERGY = -10
INTERVAL = (2, 8)
OFFSET = 8

_fa_groups = ['ref', 'alt']
_direct_groups = {
    'raw': 'raw',
    # 'rev': 'reverse',
    # 'cmp': 'complementary',
    'rmp': 'reverse complementary'
}
# raw: raw direction of transcript (from left to right)
# rev: reverse of raw
# cmp: complementary of raw
# rmp: reverse complementary of raw