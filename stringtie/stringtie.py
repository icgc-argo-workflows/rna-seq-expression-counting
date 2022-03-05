#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (c) 2022, icgc-argo-workflows

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  Authors:
    DabinJeong
"""

import os
import sys
import argparse
import pandas as pd
import subprocess

def stringtie(bam, annotation, output_pattern):
    cmd = ['stringtie',bam,'-e','-G',annotation,'-o',"{}.stringtie.gtf".format(output_pattern),'-A',"{}.stringtie".format(output_pattern)]
    try:
        subprocess.run(' '.join(cmd), check=True, shell=True)
    except Exception as e:
        sys.exit("Error: %s. stringtie failed: \n" %(e))

def main():
    """
    Python implementation of tool: stringtie

    This is auto-generated Python code, please update as needed!
    """

    parser = argparse.ArgumentParser(description='Tool: stringtie')
    parser.add_argument('-a','--annotation', dest='annotation', type=str, help='annotation file in gtf format')
    parser.add_argument('-bam', type=str, help='input bam file')
    parser.add_argument('-o','--output_pattern', dest='output_pattern', type=str)
    #parser.add_argument('-i', '--input-file', dest='input_file', type=str, help='Input file', required=True) 
    #parser.add_argument('-o', '--output-dir', dest='output_dir', type=str, help='Output directory', required=True)
    args = parser.parse_args()

    if not os.path.isfile(args.bam):
        sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.bam)

    stringtie(args.bam, args.annotation, args.output_pattern)

if __name__ == "__main__":
    main()
