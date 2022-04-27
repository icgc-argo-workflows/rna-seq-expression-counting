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
import numpy as np 
import subprocess

def quantile_normalize(df):
    """
    input: dataframe with numerical columns (sample X feature)
    output: dataframe with quantile normalized values
    """
    df_sorted = pd.DataFrame(np.sort(df.values,
                                     axis=1),
                             index=df.index,
                             columns=df.columns)
    df_mean = df_sorted.mean(axis=0)
    df_mean.index = np.arange(1, len(df_mean) + 1)
    df_qn =df.T.rank(method="min").stack().astype(int).map(df_mean).unstack()
    return(df_qn)

def upper_quartile_normalize(df):
    df_dropZeros = df.loc[(df!=0).any(axis=1)]
    df_uq = df_dropZeros/np.percentile(df,75)
    return df_uq

def main():
    """
    Python implementation of tool: normalization

    This is auto-generated Python code, please update as needed!
    """

    parser = argparse.ArgumentParser(description='Tool: normalization')
    parser.add_argument('-i', '--input-file', dest='input', type=str,
                        help='Input file', required=True)
    parser.add_argument('-geneLength',required=True,help='gene length file')
    parser.add_argument('-o', '--output-file', dest='output', type=str,
                        help='Output file', required=True)
    args = parser.parse_args()

    if not os.path.isfile(args.input):
        sys.exit('Error: specified input file %s does not exist or is not accessible!' % args.input)
    if not os.path.isfile(args.geneLength):
        sys.exit('Error: specified geneLength file %s does not exist or is not accessible!' % args.geneLength)

    df_inp = pd.read_csv(args.input,sep='\t').dropna()
    
    df_length = pd.read_csv(args.geneLength,sep='\t')
    dict_length = dict(zip(df_length['gene'],df_length['merged']))

    df_readCounts = df_inp.loc[:,['gene','readCounts']].set_index('gene')

    # FPKM
    df_FPKM = df_readCounts.apply(lambda x:x/sum(x)).apply(lambda x:x/df_readCounts.index.map(dict_length)).applymap(lambda x:x*10**9).dropna().rename(columns = {"readCounts":'FPKM'})
    df_FPKM_UQ = upper_quartile_normalize(df_FPKM)
    
   # TPM
    df_TPM = df_readCounts.apply(lambda x:x/df_readCounts.index.map(dict_length)).dropna().apply(lambda x:x/sum(x)).applymap(lambda x:x*10**6).rename(columns = {"readCounts":'TPM'})
    df_TPM_UQ = upper_quartile_normalize(df_TPM)

    df_readCounts.join(df_FPKM).join(df_FPKM_UQ,rsuffix="_UQ").join(df_TPM).join(df_TPM_UQ,rsuffix="_UQ").reset_index().to_csv(args.output,sep='\t',index=False,header=True) 
    

if __name__ == "__main__":
    main()
