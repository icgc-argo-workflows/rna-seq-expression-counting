#!/usr/bin/env python
# coding: utf-8
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='python normalize.py input output')
parser.add_argument('input',help="input tsv with gene, gene length, and read counts. Column names should be 'gene','length', and 'readCounts'.")
parser.add_argument('output',help='output tsv')

args = parser.parse_args()

df_inp = pd.read_csv(args.input,sep='\t').dropna()
dict_length = dict(zip(df_inp['gene'],df_inp['length']))

df_readCounts = df_inp.loc[:,['gene','readCounts']].set_index('gene') 

# FPKM
df_FPKM = df_readCounts.apply(lambda x:x/sum(x)).apply(lambda x:x/df_readCounts.index.map(dict_length)).applymap(lambda x:x*10**9).dropna().rename(columns = {"readCounts":'FPKM'})


# TPM
df_TPM = df_readCounts.apply(lambda x:x/df_readCounts.index.map(dict_length)).dropna().apply(lambda x:x/sum(x)).applymap(lambda x:x*10**6).rename(columns = {"readCounts":'TPM'})


df_readCounts.join(df_FPKM).join(df_TPM).reset_index().to_csv(args.output,sep='\t',index=False,header=True)


