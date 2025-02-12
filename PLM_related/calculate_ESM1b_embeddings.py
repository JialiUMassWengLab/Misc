#! /usr/bin/python

import re
import os
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Data.IUPACData import protein_letters_3to1
import torch
from transformers import AutoTokenizer, AutoModelForMaskedLM

def transformSeq(seq,sub):
    new_seq = None
    if not re.search(r'^Ter',sub)==None: #skip extension cases for now since number is small; take care of it later
        return new_seq
    if re.search(r'_',sub) == None: #Single aa change
        match = re.search(r'\d+',sub)
        pos = int(match.group())-1
        new = re.sub(r'^[A-Z][a-z][a-z]\d+','',sub)
        if new == 'del':
            new_seq = seq[:pos] + seq[pos+1:]
        elif new == 'dup':
            new_seq = seq[:pos] + seq[pos] + seq[pos:]
        else: #delins and single aa substitution
            new_aa = ''
            if not re.search(r'^delins',new)==None:
                new = new[6:]
            for i in range(0,len(new),3):
                if new[i:i+3] in protein_letters_3to1:
                    new_aa += protein_letters_3to1[new[i:i+3]]
            new_seq = seq[:pos] + new_aa
            if new[-3:] != 'Ter':
                new_seq += seq[pos+1:]
    else:
        match = re.search(r'^[A-Z][a-z][a-z]\d+_[A-Z][a-z][a-z]\d+',sub)
        pos1 = int(re.search(r'\d+',match.group().split('_')[0]).group())-1
        pos2 = int(re.search(r'\d+',match.group().split('_')[1]).group())-1
        new = re.sub(r'^[A-Z][a-z][a-z]\d+_[A-Z][a-z][a-z]\d+','',sub)
        if new == 'del':
            new_seq = seq[:pos1] + seq[pos2+1:]
        elif new == 'dup':
            new_seq = seq[:pos1] + seq[pos1:pos2+1] + seq[pos1:]
        elif new[:3] == 'ins':
            new = new[3:]
            new_aa = ''
            for i in range(0,len(new),3):
                if new[i:i+3] in protein_letters_3to1:
                    new_aa += protein_letters_3to1[new[i:i+3]]
            new_seq = seq[:pos1+1] + new_aa
            if new[-3:] != 'Ter':
                new_seq += seq[pos2:]
        elif not re.search(r'delins',new)==None:
            new = new[6:]
            new_aa = ''
            for i in range(0,len(new),3):
                if new[i:i+3] in protein_letters_3to1:
                    new_aa += protein_letters_3to1[new[i:i+3]]
            new_seq = seq[:pos1] + new_aa
            if new[-3:] != 'Ter':
                new_seq += seq[pos2+1:]
                
    return new_seq
                        

def WT_embeddings():
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    model = AutoModelForMaskedLM.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    record_dict = SeqIO.index("mRNA_translated_seqs.fa", "fasta")
    
    input_seqs = []
    input_ids = []
    for txID,record in record_dict.items():
        if len(record.seq) < 1023:
            input_seqs.append(record.seq)
            input_ids.append(record.id)
            
    embedding_tensors = {}
    batchSize = 1000
    for i in range(0, len(input_seqs), batchSize):
        inputs = tokenizer([str(s) for s in input_seqs[i:i+batchSize]], return_tensors="pt", padding=True)
        with torch.no_grad():
            embeddings = model(**inputs, output_hidden_states=True).hidden_states
            for j in range(len(inputs['attention_mask'])):
                embedding_tensors.update({input_ids[j+i]+'_WT':torch.mean(embeddings[-1][j][torch.nonzero(inputs['attention_mask'][j],as_tuple=True)],0)})

    torch.save(embedding_tensors,'embedding_tensors.wildtype.pt')

def MT_embeddings():
    tokenizer = AutoTokenizer.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    model = AutoModelForMaskedLM.from_pretrained("facebook/esm1b_t33_650M_UR50S")
    record_dict = SeqIO.index("mRNA_translated_seqs.fa", "fasta")

    #df1 = pd.read_csv('Orphanet_LoF_GoF_clinvarVariants.csv',index_col=0)
    #df2 = pd.read_csv('Benign_clinvarVariants.csv',index_col=0)
    df1 = pd.read_csv('OMIM_LoF_GoF_clinvarVariants.csv',index_col=0)
    df2 = pd.read_csv('Benign_clinvarVariants_OMIM.csv',index_col=0)
    df2['direction'] = ['Benign'] * df2.shape[0]
    df = pd.concat([df1,df2],join='inner',axis=0)
    df['sub'] = df['Name'].map(lambda x: re.search(r'\(p\.\w+\)$',x).group()[3:-1] if not re.search(r'\(p\.\w+\)$',x)==None else '')
    df = df.query('sub!=""')
    
    input_seqs = []
    input_ids = []
    for index,row in df.iterrows():
        if re.search(r'fs$',row['sub'])==None:
            txID = re.search(r'^N\w+\.\d+\(',row['Name'])
            if not txID == None:
                txID = txID.group()[:-1]
                record = record_dict[txID]
                if len(record.seq) < 1023:
                    mu_seq = transformSeq(record.seq,row['sub'])
                    if mu_seq == None:
                        #print(str(index),txID,row['sub'])
                        continue
                    elif len(mu_seq) < 1023:
                        input_seqs.append(mu_seq)
                        input_ids.append(':'.join([str(index),txID,row['sub']]))
                        
    embedding_tensors = {}
    batchSize = 1000
    for i in range(0, len(input_seqs), batchSize):
        print('Processing sequences %d to %d...' % (i,i+batchSize))
        inputs = tokenizer([str(s) for s in input_seqs[i:i+batchSize]], return_tensors="pt", padding=True)
        with torch.no_grad():
            embeddings = model(**inputs, output_hidden_states=True).hidden_states
            for j in range(len(inputs['attention_mask'])):
                embedding_tensors.update({input_ids[j+i]:torch.mean(embeddings[-1][j][torch.nonzero(inputs['attention_mask'][j],as_tuple=True)],0)})

    #torch.save(embedding_tensors,'embedding_tensors.training.pt')
    torch.save(embedding_tensors,'embedding_tensors.testing.pt')


def main():
    #WT_embeddings()
    MT_embeddings()


if __name__=='__main__':
    main()
