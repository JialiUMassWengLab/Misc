#! /usr/bin/python

import re
import os
import sys
import glob
import pickle
import numpy as np
import pandas as pd
import networkx as nx

def const_network(prefix):
    df1 = pd.read_csv(prefix+'.partners.csv',index_col=0)
    df2 = pd.read_csv(prefix+'.partners.degrees.csv',index_col=0)
    G = nx.Graph()
    for index,row in df1.iterrows():
        partners = row['partner'].split('|')
        db = row['database']
        goi = row['GOI']
        if len(partners) > 50:
            continue
        if not goi in G.nodes():
            G.add_node(goi)
        for partner in partners:
            if not re.search(r'^[A-Z][A-Z0-9]+',partner):
                continue
            partner_info = df2.query('partner==@partner & database==@db')
            if partner_info.shape[0] == 0:
                continue
            degree = partner_info.iloc[0,2]
            if degree > 50:
                continue
            wt = 1./float(len(partners) * degree)
            if not (goi,partner) in G.edges():
                G.add_edge(goi,partner,weight=wt)
            else:
                G[goi][partner]['weight'] += wt
    return G
                
def run_pagerank(prefix,G):
    with open(prefix+'_list.txt','r') as f:
        GOI_list = [i[:-1] for i in f.readlines()]

    pr = pd.Series(nx.pagerank(G,weight='weight',max_iter=1000))
    res = pd.DataFrame(pr,columns=['pageRank'])
    res['is_goi'] = res.index.map(lambda x: x in GOI_list)
    nx.set_node_attributes(G,res[['is_goi']].transpose().to_dict())
    nx.set_node_attributes(G,res[['pageRank']].transpose().to_dict())
    return G,res

def annotate_partners(prefix):
    res = pd.read_csv(prefix+'.graph.pageRank.csv',index_col=0).iloc[:,:3]
    not_supported_list = res[res['is_goi']==False].index.to_list()
    df = pd.read_csv(prefix+'.partners.csv',index_col=0)
    data = {}
    for db in list(df['database'].unique()):
        data.update({'genetically_supported_partners_'+db:{}})
    for pair,row in df.iterrows():
        name = 'genetically_supported_partners_'+row['database']
        #for partner in [i for i in row['partner'].split('|') if i in not_supported_list]:
        for partner in [i for i in row['partner'].split('|') if i in res.index.to_list()]:
            if partner == row['GOI']:
                continue
            if partner in data[name]:
                data[name][partner].append(row['GOI'])
            else:
                data[name].update({partner:[row['GOI']]})

    for db,series in data.items():
        for partner in series.keys():
            data[db][partner] = '|'.join(data[db][partner])
            
    out = pd.DataFrame.from_dict(data).fillna('')
    res = res.merge(out,left_index=True,right_index=True,how='left').fillna('')
    return res
    

def main():
    prefix = sys.argv[1]
    #G = const_network(prefix)
    #G,res = run_pagerank(prefix,G)
    #pickle.dump(G, open('%s.graph.pickle' % prefix, 'wb'))
    
    #Bcell = pd.read_csv('../B_cellmarker_1060.txt',sep='\t',index_col=0)
    #res = res.merge(Bcell,left_index=True,right_index=True,how='left').fillna(0)
    res = annotate_partners(prefix)
    res = res.sort_values('pageRank',ascending=False)
    res.to_csv(prefix+'.graph.pageRank.csv')


if __name__=='__main__':
    main()
