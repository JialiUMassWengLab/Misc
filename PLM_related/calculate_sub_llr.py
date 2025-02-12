#! /usr/bin/python

import re
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Data.IUPACData import protein_letters_3to1
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

def processBatch(df,input_seqs,input_muts,tokenizer,model):
    inputs = tokenizer([str(s) for s in input_seqs], return_tensors="pt",padding=True)
    with torch.no_grad():
        logits = model(**inputs).logits

    data = {'ref':{},'alt':{},'ref_logit':{},'alt_logit':{},'max_logit_aa':{}}
    for i,index in enumerate(input_muts):
        S = df.loc[index]
        if not S['sub'][:3] in protein_letters_3to1 or not S['sub'][-3:] in protein_letters_3to1:
            #print(S)
            continue
        ref_id = tokenizer.convert_tokens_to_ids(protein_letters_3to1[S['sub'][:3]])
        alt_id = tokenizer.convert_tokens_to_ids(protein_letters_3to1[S['sub'][-3:]])
        pos = (inputs.input_ids == tokenizer.mask_token_id)[i].nonzero(as_tuple=True)[0]
        data['ref'].update({index: protein_letters_3to1[S['sub'][:3]]})
        data['alt'].update({index: protein_letters_3to1[S['sub'][-3:]]})
        data['ref_logit'].update({index: float(logits[i,int(pos),ref_id])})
        data['alt_logit'].update({index: float(logits[i,int(pos),alt_id])})
        data['max_logit_aa'].update({index: tokenizer.decode(logits[i,int(pos)].argmax(axis=-1))})
        
    res = pd.DataFrame.from_dict(data)
    res['llr'] = res['alt_logit'] - res['ref_logit']
    return res


def main():
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    model = AutoModelForMaskedLM.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    record_dict = SeqIO.index("mRNA_translated_seqs.fa", "fasta")
    
    df1 = pd.read_csv('Orphanet_LoF_GoF_clinvarVariants.csv',index_col=0)
    df2 = pd.read_csv('Benign_clinvarVariants.csv',index_col=0)
    df = pd.concat([df1,df2],join='inner',axis=0)
    df['sub'] = df['Name'].map(lambda x: re.search(r'\(p\.\w+\)$',x).group()[3:-1] if not re.search(r'\(p\.\w+\)$',x)==None else '')
    df = df.query('sub!=""')

    input_seqs = []
    input_muts = []
    for index,row in df.iterrows():
        match = re.search(r'\d+',row['sub'])
        if match != None and re.search(r'Ter$',row['sub'])==None and re.search(r'fs$',row['sub'])==None and re.search(r'_',row['sub'])==None:
            txID = re.search(r'^N\w+\.\d+\(',row['Name'])
            if not txID == None:
                txID = txID.group()[:-1]
                aa_seq = record_dict[txID].seq
                if len(aa_seq) < 1023:
                    pos = int(match.group())-1
                    input_muts.append(index)
                    input_seqs.append(aa_seq[:pos]+'<mask>'+aa_seq[pos+1:])
    print(len(input_seqs))
    
    array = []
    batchSize = 1000
    for i in range(0, len(input_seqs), batchSize):
        res = processBatch(df,input_seqs[i:i+batchSize],input_muts[i:i+batchSize],tokenizer,model)
        array.append(res)

    out = pd.concat(array,join='inner',axis=0)
    out = pd.concat([out,df['ClinicalSignificance']],join='inner',axis=1)
    out.to_csv('Orphanet_llr.csv')


if __name__=='__main__':
    main()
