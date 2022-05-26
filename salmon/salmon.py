#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
  Copyright (c) 2021, icgc-argo-rna-wg

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
import subprocess

def salmon_index(refseq, index):
    cmd = ['salmon','index','-t',refseq,'-i',index]
    try:
        subprocess.run(' '.join(cmd), check=True, shell=True)
    except Exception as e:
        sys.exit("Error: %s. salmon index failed: %s\n" %(e,refseq))

def salmon_quant(index, libType, read1, read2, outdir):
    cmd = ['salmon', 'quant', '-i', index, '--libType', libType, '-1', read1, '-2', read2, '-o', outdir]
    try:
        subprocess.run(' '.join(cmd), check=True, shell=True)
    except Exception as e:
        sys.exit("Error: %s. salmon index failed: %s, %s\n" %(e,read1,read2))

def main():
    """
    Python implementation of tool: salmon

    This is auto-generated Python code, please update as needed!
    """

    parser = argparse.ArgumentParser(description='Tool: salmon')
    parser.add_argument('--referenceSeq',required=True)
    parser.add_argument('--index',required=True)
    parser.add_argument('-r1','--read1', dest='read1',type=str, help='Input file: read1', required=True)
    parser.add_argument('-r2','--read2', dest='read2',type=str, help='Input file: read2', required=True)
    parser.add_argument('-o', '--output-dir', dest='output_dir', type=str,
                        help='Output directory', required=True)
    parser.add_argument('--threads')
    parser.add_argument('--mem')
    args = parser.parse_args()

    if not os.path.isfile(args.read1):
        sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.input_file)
    if not os.path.isfile(args.read2):
        sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.input_file)

    #if not os.path.isdir(args.index):
    #    sys.exit('Error: index directory of salmon %s does not exist or is not accessible!' % args.index)
    #if not os.path.isdir(args.output_dir):
    #    sys.exit('Error: specified output dir %s does not exist or is not accessible!' % args.output_dir)

    subprocess.run("mkdir -p {}".format(args.output_dir), check=True, shell=True)
    
    salmon_index(args.referenceSeq, args.index)

    salmon_quant(args.index, "A", args.read1, args.read2, args.output_dir)

if __name__ == "__main__":
    main()

