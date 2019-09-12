#! /usr/bin/python

import re
import os
import sys
import json
import requests
import numpy as np
import pandas as pd
import xml.etree.ElementTree as ET

def APIget(gene):
    #r = requests.get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=10&sort=relevance&term='+gene+'%20AND%20senescence')
    r = requests.get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&retmode=json&retmax=20&sort=relevance&term=(lipopolysaccharide%5BTitle%2FAbstract%5D%20OR%20LPS%5BTitle%2FAbstract%5D)%20AND%20'+gene+'%5BTitle%2FAbstract%5D')
    idList = json.loads(r.text)["esearchresult"]["idlist"]
    s = ','.join(idList)
    r = requests.get('http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&rettype=abstract&id='+s)
    root = ET.fromstring(r.text.encode('utf-8'))
    #for pubmedID in root.iter('PMID'):
    #    print pubmedID.text
    data = {'title':{}, 'abstract':{}, 'journal':{}}
    for article in root.findall('PubmedArticle'):
        info = article.find('MedlineCitation')
        pubmedid = info.find('PMID').text
        journal = info.find('Article/Journal/Title').text
        abstract = info.find('Article/Abstract')
        title = info.find('Article/ArticleTitle').text
        if pubmedid is not None and abstract is not None and journal is not None and title is not None:
            text = ''
            for seg in abstract.iter('AbstractText'):
                label = seg.get('Label')+': ' if seg.get('Label') else ''
                if seg.text is not None:
                    text += label+' '+seg.text
            data['title'].update({pubmedid: title})
            data['abstract'].update({pubmedid: text})
            data['journal'].update({pubmedid: journal})

    df = pd.DataFrame.from_dict(data)
    #print df
    return df


def main():
    genes = pd.read_csv('temp_pattern_nonBlood/nonBlood.sig.csv',header=None,index_col=1)
    for gene in genes.index:
        print gene
        Df = APIget(gene)
        if Df.shape[0]==0:
            continue
        with open('papers.%s.txt' % gene,'w') as out:
            for pmid,row in Df.iterrows():
                out.write('Pubmed ID:\n'+pmid+'\n')
                for index,text in row.iteritems():
                    out.write(index+':\n'+text.encode('utf-8')+'\n')
                out.write('\n\n')
            


if __name__=='__main__':
    main()
