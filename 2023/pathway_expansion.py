#! /usr/bin/python

import re
import os
import sys
import time
import numpy as np
import pandas as pd
import requests as rs

def get_Reactome_components(stId):
    url = 'https://reactome.org/ContentService/data/complex/%s/subunits?excludeStructures=false' % stId
    r = rs.get(url)
    components = []
    if r.status_code == 200:
        for comp in r.json():
            if comp['className']=='Protein':
                #components.append(comp['name'][0])
                components.append(comp['stId'])
        return components

def expand_complex_into_pairwise_interactions(df,comp_col,geneList):
    array = []
    for index,row in df.iterrows():
        components = list(row[comp_col])
        for i in [i for i in components if i in geneList]:
            for j in [j for j in components if not j==i]:
                new_row = row.copy()
                new_row['GOI'] = i
                new_row['partner'] = j
                array.append(new_row)
    new_df = pd.concat(array,join='inner',axis=1).transpose()
    return new_df

def convert_query_output_to_df_IntAct(results,fields2keep):
    array = []
    for item in results:
        data = {}
        for field in fields2keep:
            data.update({field:item[field]})
        array.append(pd.Series(data,name=item['ac']))
    if len(array)> 0:
        return pd.concat(array,axis=1).transpose()
    else:
        return pd.DataFrame(columns=fields2keep)


def query_BioGRID(geneList):
    url = 'https://webservice.thebiogrid.org/interactions/'
    params = {
        'accessKey': '476c99426ff0fb4d17a7fdc52d702afa', 'interSpeciesExcluded': 'true', \
        'searchNames': 'true', 'includeInteractors': 'true', 'taxId': '9606', 'format': 'json', \
        'geneList': '|'.join(geneList),
    }
    
    r = rs.get(url,params=params)
    if r.status_code==200:
        df = pd.read_json(r.content).transpose()
        #filtering criteria
        df = df.query('EXPERIMENTAL_SYSTEM_TYPE == "physical"')
        df.loc[:,'GOI'] = df.apply(lambda x: x['OFFICIAL_SYMBOL_A'] if x['OFFICIAL_SYMBOL_A'] in geneList else x['OFFICIAL_SYMBOL_B'],axis=1)
        df.loc[:,'partner'] = df.apply(lambda x: x['OFFICIAL_SYMBOL_B'] if x['OFFICIAL_SYMBOL_A'] in geneList else x['OFFICIAL_SYMBOL_A'],axis=1)
        df.loc[:,'database'] = ['BioGRID'] * df.shape[0]
        return df
    else:
        print('Query failed with status code: %d' % r.status_code)


def query_IntAct(geneList):
    #first get interactorACs
    base_url = 'https://www.ebi.ac.uk/intact/ws/interactor/findInteractor/'
    params = {'pageSize': 50}
    interactor_ids = {}    
    for goi in geneList:
        url = base_url + goi
        r = rs.get(url,params=params)
        if r.status_code==200:
            results = r.json()['content']
            for item in results:
                #keep only human proteins and exact match to GOI
                if item['interactorTaxId']==9606 and item['interactorName']==goi:
                    if goi in interactor_ids:
                        interactor_ids[goi].append(item['interactorAc'])
                    else:
                        interactor_ids.update({goi: [item['interactorAc']]})
        else:
            print('Query failed with status code: %d' % r.status_code)
    
    #query interactions
    fields2keep = ['uniqueIdA','uniqueIdB','moleculeA','moleculeB','detectionMethod','type','hostOrganism','intactMiscore','publicationPubmedIdentifier']
    url = 'https://www.ebi.ac.uk/intact/ws/interaction/findInteractionWithFacet'
    params = {
        'batchSearch': 'true', 'intraSpeciesFilter': 'true', 'maxMIScore': 1,'minMIScore': 0, \
        'interactionTypesFilter': ['association', 'physical association', 'direct interaction'], \
        'interactorSpeciesFilter': 'Homo sapiens', 'interactorTypesFilter': 'protein', \
        'pageSize': 200, 'page': 0, 'negativeFilter': 'POSITIVE_ONLY',
    }
    df_array = []
    '''
    #gene by gene query
    for goi in geneList:
        #skip if no IDs found for GOI
        if not goi in interactor_ids:
            continue
        
        params.update({'query': ','.join(interactor_ids[goi])})
        r = rs.post(url,params=params)
        if r.status_code==200:
            results = r.json()['data']['content']
            df_array.append(convert_query_output_to_df_IntAct(results,fields2keep))
        else:
            print('Query failed with status code: %d' % r.status_code)
    '''
    #batch query
    import itertools
    idList = list(itertools.chain.from_iterable(interactor_ids.values()))
    params.update({'query': ','.join(idList)})
    r = rs.post(url,params=params)
    if r.status_code==200:
        output = r.json()
        page_num = output['data']['totalPages']
        total_num = output['data']['totalElements']
        print('Query IntAct with %d hits in %d pages' % (total_num,page_num))
        results = output['data']['content']
        df_array.append(convert_query_output_to_df_IntAct(results,fields2keep))
        for p in range(1,(page_num+1)):
            params.update({'page':p})
            r = rs.post(url,params=params)
            if r.status_code==200:
                results = r.json()['data']['content']
                df_array.append(convert_query_output_to_df_IntAct(results,fields2keep))
            else:
                print('Query failed with status code: %d' % r.status_code)
                break
    else:
        print('Query failed with status code: %d' % r.status_code)
    
    df = pd.concat(df_array,join='inner',axis=0)
    #filtering criteria
    df = df.query('type == "physical association" | type == "association"')
    df.loc[:,'GOI'] = df.apply(lambda x: x['moleculeA'] if x['moleculeA'] in geneList else x['moleculeB'],axis=1)
    df.loc[:,'partner'] = df.apply(lambda x: x['moleculeB'] if x['moleculeA'] in geneList else x['moleculeA'],axis=1)
    df.loc[:,'database'] = ['IntAct'] * df.shape[0]
    return df


def query_STRING(geneList):
    url = 'https://string-db.org/api/json/interaction_partners'
    params = {'identifiers': '%0d'.join(geneList), 'species': 9606, \
              'network_type': 'physical', 'limit': 500}
    r = rs.get(url,params=params)
    if r.status_code==200:
        df = pd.read_json(r.content)
        #filtering criteria
        df = df.query('escore > 0 | dscore > 0')
        df.loc[:,'GOI'] = df.apply(lambda x: x['preferredName_A'] if x['preferredName_A'] in geneList else x['preferredName_B'],axis=1)
        df.loc[:,'partner'] = df.apply(lambda x: x['preferredName_B'] if x['preferredName_A'] in geneList else x['preferredName_A'],axis=1)
        df.loc[:,'database'] = ['STRING'] * df.shape[0]
        return df
    else:
        print('Query failed with status code: %d' % r.status_code)


def query_ReactomeComplex(geneList):
    url = 'https://reactome.org/ContentService/search/query'
    params = {
        'query': ','.join(geneList), 'species': 'Homo sapiens', 'types': 'Complex', \
        'cluster': 'true', 'parserType': 'STD', 'rows': 500,
    }
    r = rs.get(url,params=params)
    if r.status_code==200:
        results = r.json()["results"][0]['entries']
        array = []
        for item in results:
            array.append(pd.Series(item))
        df = pd.concat(array,axis=1).transpose()
        #filtered criteria
        df = df.loc[df.apply(lambda x: x['species'][0]=='Homo sapiens',axis=1),:]
        df.loc[:,'Components'] = df.apply(lambda x: get_Reactome_components(x['stId']),axis=1)
        
        import itertools
        idList = list(itertools.chain.from_iterable(df['Components']))
        geneNames = {}
        for ID in idList:
            r = rs.get('https://reactome.org/ContentService/data/query/%s/referenceEntity' % ID)
            if r.status_code == 200:
                #print(r.text.split('\t')[1].split(' ')[1])
                geneNames.update({ID: r.text.split('\t')[1].split(' ')[1]})
                
        df.loc[:,'Components2'] = df.apply(lambda x: list(map(lambda y: geneNames[y], x['Components'])),axis=1)
        new_df = expand_complex_into_pairwise_interactions(df,'Components2',geneList)
        new_df.loc[:,'database'] = ['Reactome'] * new_df.shape[0]
        return new_df
    
    else:
        print('Query failed with status code: %d' % r.status_code)


def query_LIANA(geneList):
    df = pd.read_csv('/mnt/efs/home/iet5740/Projects/GSP/pathway_expansion/LIANA.ligand_receptor.csv')
    #filtering for interactions containing GOIs
    df = df.query('source_genesymbol == @geneList | target_genesymbol == @geneList')
    df.loc[:,'GOI'] = df.apply(lambda x: x['source_genesymbol'] if x['source_genesymbol'] in geneList else x['target_genesymbol'],axis=1)
    df.loc[:,'partner'] = df.apply(lambda x: x['target_genesymbol'] if x['source_genesymbol'] in geneList else x['source_genesymbol'],axis=1)
    df.loc[:,'database'] = ['LIANA'] * df.shape[0]
    return df


def query_complexPortal(geneList):
    df = pd.read_csv('/mnt/efs/home/iet5740/Projects/GSP/pathway_expansion/complex-portal-9606.map2geneNames.tsv',sep='\t')
    df.loc[:,'Components'] = df.apply(lambda x: str(x['Components']).split('|'),axis=1)
    #filter for complexes containing GOIs
    df = df.loc[df.apply(lambda x: len(x['Components'])>1 and len([i for i in x['Components'] if i in geneList])>0,axis=1),:]    
    new_df = expand_complex_into_pairwise_interactions(df,'Components',geneList)
    new_df.loc[:,'database'] = ['complexPortal'] * new_df.shape[0]
    return new_df
    

def compile_partners(geneList):
    biogrid_df = query_BioGRID(geneList)
    intact_df = query_IntAct(geneList)
    string_df = query_STRING(geneList)
    reactome_df = query_ReactomeComplex(geneList)
    liana_df = query_LIANA(geneList)
    complexportal_df = query_complexPortal(geneList)
    array = []
    for df in [biogrid_df, intact_df, string_df, reactome_df, liana_df, complexportal_df]:
        df2 = df[['GOI','partner','database']].groupby('GOI').agg(lambda x: set(x)).reset_index()
        df2.loc[:,'database'] = df2.apply(lambda x: list(x['database'])[0],axis=1)
        df2.index = df2.apply(lambda x: ' - '.join([x['database'],x['GOI']]),axis=1)
        array.append(df2)
    df = pd.concat(array,join='inner',axis=0)
    return df

    
def main():
    geneList = []
    with open('../indication_expansion_target_list.txt','r') as f:
        geneList = list(map(lambda x: x[:-1],f.readlines()))
    print(geneList)

    #df = query_complexPortal(geneList)
    #df = query_ReactomeComplex(geneList)
    #df = query_IntAct(geneList)    
    df = compile_partners(geneList)
    df.loc[:,'partner'] = df.apply(lambda x: '|'.join(list(x['partner'])),axis=1)
    df.to_csv('indication_expansion_partners.csv')
    print(df)
    

if __name__=='__main__':
    main()
