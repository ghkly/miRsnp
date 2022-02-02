#!/usr/bin/python3
# _*_ coding: utf-8 _*_
# @Time     : 2022/1/30 18:44


import re
import gzip


def fmt_gtf(fp, trans_type):
    if fp.endswith('.gz'):
        fr = gzip.open(fp, 'rb')
    else:
        fr = open(fp, 'rb')
    for line in fr:
        line = line.decode()
        if not line.startswith('#'):
            arr = line.split('\t')
            desp = arr[8]
            if arr[2] == trans_type:
                trans = re.search("transcript_id (.*?);", desp, re.S).group(1).strip('"')
                try:
                    gene = re.search("gene_name (.*?);", desp, re.S).group(1).strip('"')
                except AttributeError:
                    gene = re.search("gene_id (.*?);", desp, re.S).group(1).strip('"')
                exon_number = re.search("exon_number (.*?);", desp, re.S).group(1).strip('"')
                row = [
                    trans,
                    gene,
                    arr[0],
                    arr[6],
                    '{}:{}'.format(arr[3], arr[4]),
                    str(int(arr[4]) - int(arr[3]) + 1),
                    exon_number
                ]
                yield row
    fr.close()


def parse_gtf(fp, trans_type='exon'):
    """
    return: {'ENST00000497304': ['INTS11', '1', '-', '1312226:1312700|1312018:1312147|1311598:1311924'], ...}
    """
    data = fmt_gtf(fp, trans_type)
    res = {}
    for row in data:
        transcript, interval, exon_number = row[0], row[4], row[-1]
        if transcript not in res:
            res[transcript] = {}
        if exon_number == '1':
            res[transcript] = row[1:5]
        else:
            res[transcript][3] += "|{}".format(interval)
    return res


if __name__ == "__main__":
    gtf = r'H:\resource\test.gtf.gz'
    parse_gtf(gtf)
