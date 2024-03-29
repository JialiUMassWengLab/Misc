import pandas as pd
import requests as rs
import re
import os
import boto3
from botocore.client import Config
from transformers import AutoTokenizer, AutoModelForMaskedLM
from torch.nn.functional import cosine_similarity
import torch
import pronto
#from disqover import *

def getEFOterms(efo_term):
    response = rs.get('https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms?id=%s' % efo_term)
    if response.status_code == 200:
        results = response.json()['_embedded']['terms'][0]
        return results['label'],results['description']
    else:
        print('Error when querying EFO term %s; code %d' % (efo_term,response.status_code))
        
def construct_efo_text(efo_term):
    response = rs.get('https://www.ebi.ac.uk/ols4/api/ontologies/efo/terms?id=%s' % efo_term)
    if response.status_code == 200:
        results = response.json()['_embedded']['terms'][0]
        #print(efo_term,results['description'])
        s = "[CLS] " + results['label']
        if len(results['description'])>0:
            s = s + " is " + results['description'][0].split('.')[0]
        return s
    else:
        print('Error when querying EFO term %s; code %d' % (efo_term,response.status_code))

        
def get_efo_list():
    # load AWS S3 credentials
    session = boto3.Session(profile_name='insight')
    bucket_name = 'tak-insight-priv-compbio-silv'
    path_prefix = 'home/chris_deboever/gsp_ie/20231203/'
    # connect to S3 using boto
    client = session.client("s3")

    # get a list of all genes
    paginator = client.get_paginator('list_objects_v2')
    pages = paginator.paginate(Bucket=bucket_name, Prefix=path_prefix)
    geneList = []
    for page in pages:
        geneList += [os.path.basename(x['Key']).replace(' GSP.xlsx','') for x in page['Contents'] if '.xlsx' in x['Key']]
    #print(geneList)

    #read GWAS summary dataframe
    array = []
    for gene in geneList:
        gwas_path = os.path.join(path_prefix,gene+' GSP','GWAS Summary.tsv.gz')
        gwas_df = pd.read_csv(client.get_object(Bucket=bucket_name, Key=gwas_path)['Body'],sep='\t',compression='gzip')
        array.append(gwas_df)
        
    summary_df = pd.concat(array,join='inner',axis=0)
    summary_df['study_traitEfos'] = summary_df['study_traitEfos'].fillna('')
    efoList = set(summary_df['study_traitEfos'])
    efoList = [i for i in efoList if not re.search(r'^EFO',i)==None and re.search(r',',i)==None]
    return efoList


def get_embeddings(text_list,tokenizer,model):
    mod_text_list = [x for x in text_list]
    inputs = tokenizer(mod_text_list, return_tensors="pt", truncation=True, max_length=50, padding=True)
    print(inputs['input_ids'].shape)
    with torch.no_grad():
        embeddings = model(**inputs, output_hidden_states=True).hidden_states

    return embeddings[-1][:,0,:]
    #return torch.mean(embeddings[-1],1)

def get_HPO_embeddings(all_hpo_terms,tokenizer,model):   
    #all_hpo_names = [x for x in all_hpo_terms.keys()]
    all_hpo_names = [x for x in all_hpo_terms]
    inputs = tokenizer(all_hpo_names, return_tensors="pt", truncation=True, max_length=50, padding=True)
    print(inputs['input_ids'].shape)
    with torch.no_grad():
        embeddings = model(**inputs, output_hidden_states=True).hidden_states

    return embeddings[-1][:,0,:]

def get_EFO_embeddings(efo_list,tokenizer,model):
    #all_efo_names = [getEFOterms(i) for i in efo_list]
    #all_efo_names = [construct_efo_text(i) for i in efo_list]
    all_efo_names = [i for i in efo_list]
    inputs = tokenizer(all_efo_names, return_tensors="pt", truncation=True, max_length=50, padding=True)
    with torch.no_grad():
        embeddings = model(**inputs, output_hidden_states=True).hidden_states

    return embeddings[-1][:,0,:]
            

def checkOneTerm(efo_index,hpo_def_embed,hpo_name_embed,efo_def_embed,efo_name_embed):
    s_name = cosine_similarity(efo_name_embed[efo_index],hpo_name_embed)
    s_def = cosine_similarity(efo_def_embed[efo_index],hpo_def_embed)
    s_both = torch.mul(s_name,s_def)
    return sa,sb,sab


def main():
    efo_list = get_efo_list()
    print(len(efo_list))
    efo_name_list = []
    efo_def_list = []
    efo_indices = {}
    i = 0
    for term in efo_list:
        name,desc = getEFOterms(term)
        name = name.replace('obsolete_','')
        if len(desc) > 0:
            desc = re.sub(r'\s\([\w\.\s,;/-]+\)','',desc[0]).split('.')[0]
        else:
            desc = ""
        efo_name_list.append(name)
        efo_def_list.append(desc)
        efo_indices.update({i: term})
        i += 1
    
    hpo = pronto.Ontology.from_obo_library('hp.obo')
    pheno_abnormality_terms = set(hpo['HP:0000118'].subclasses())
    hpo_indices = {}
    hpo_def_list = []
    hpo_name_list = []
    i = 0
    for entity in hpo.terms():
        if entity in pheno_abnormality_terms:
            name = entity.name
            definition = ""
            if not entity.definition == None:
                definition = re.sub(r'\s\([\w\.\s,;/-]+\)','',entity.definition).split('.')[0]
            hpo_name_list.append(name)
            hpo_def_list.append(definition)
            hpo_indices.update({i: entity.id})
            i += 1
        
    tokenizer = AutoTokenizer.from_pretrained("microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
    model = AutoModelForMaskedLM.from_pretrained("microsoft/BiomedNLP-BiomedBERT-base-uncased-abstract-fulltext")
    A = get_embeddings(hpo_def_list,tokenizer,model)
    B = get_embeddings(hpo_name_list,tokenizer,model)
    C = get_embeddings(efo_def_list,tokenizer,model)
    D = get_embeddings(efo_name_list,tokenizer,model)

    pheno_matches = {'efo_name':{},'hpo_hit_id':{},'hpo_hit_name':{},'score':{}}
    for i,efo in enumerate(efo_name_list):
        efo_id = efo_indices[i]
        s_name = cosine_similarity(D[i],B)
        s = s_name
        if not efo_def_list[i] == "" and torch.max(s_name) < .999:
            s_def = cosine_similarity(C[i],A)
            s_both = torch.mul(s_name,s_def)
            s = s_both

        best_index = int(s.argmax())
        pheno_matches['efo_name'].update({efo_id: efo})
        pheno_matches['hpo_hit_id'].update({efo_id: hpo_indices[best_index]})
        pheno_matches['hpo_hit_name'].update({efo_id: hpo_name_list[best_index]})
        pheno_matches['score'].update({efo_id: s[best_index]})
        print('%s\t%s' % (efo,hpo_name_list[best_index]))
        
    df = pd.DataFrame.from_dict(pheno_matches)
    df.to_csv('LLM_both_matches.csv')
    print(df)

    
if __name__ == '__main__':
    main()
