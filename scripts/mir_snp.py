#!/usr/bin/python3

import os
import sys

cdir = os.path.abspath(os.path.dirname(__file__))
topdir = os.path.dirname(cdir)
sys.path.append(topdir)
from mirsnp.mrsp import main

if __name__ == "__main__":
    args = None
    if len(sys.argv) == 1:
        args = ["--help"]
    main(sys.argv[1:])
