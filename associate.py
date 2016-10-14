import pandas as pd
import csv
import re
import sys
import numpy
import functools
import scipy.stats as st

def cleanFloat(num):
    if isinstance(num,str):
        num = num.strip()
        num = re.sub(r'\`','',num)
        num = re.sub(r'\.\.','.',num)
        if not re.search(r'\d',num):
            num = None
    return num

def GT2num(genotype,minor_allele):
    if genotype==None:
        return numpy.nan
    else:
        match = re.search(r'([\w-]+)/([\w-]+)',genotype)
        minor_num = 0
        if match:
            if match.group(1) == minor_allele:
                minor_num += 1
            if match.group(2) == minor_allele:
                minor_num += 1

        return minor_num

def getSNPs():
    data = pd.read_csv("../all_SNP_calls.csv",header=None,names=['rs','gt','sid','pid'])
    data = data[map(lambda x: not re.search(r'[\w-]+/[\w-]+',x)==None, data['gt'])][['rs','gt','pid']]
    data = data.drop_duplicates()
    data = data[~data[['pid','rs']].duplicated(keep=False)]
    data = data.pivot(index='pid',columns='rs')
    data.columns = map(lambda x: x[1], data.columns)
    return data

def getMetabolites():
    data = pd.read_csv("EPIDEMIC.concentration.csv",header=None,names=['pid','age','height','weight','metabolites','concentration'])
    data = data[~(data['pid']=='SUBJECT ID')]
    data = data.drop_duplicates()
    data = data[~data[['pid','metabolites']].duplicated()]
    data['height'][data['height']=='N/a'] = numpy.nan
    data['height'][data['height']=='Not Listed'] = numpy.nan
    data['weight'][data['weight']=='N/a'] = numpy.nan
    data['weight'][data['weight']=='Not Listed'] = numpy.nan
    data['age'] = data['age'].astype('float')
    data['height'] = data['height'].astype('float')
    data['weight'] = data['weight'].astype('float')

    metas_to_filter = []
    for key,value in data['metabolites'].value_counts().iteritems():
        if value < 100:
            metas_to_filter.append(key)

    wide = data[['pid','metabolites','concentration']].pivot(index='pid',columns='metabolites')
    wide.columns = map(lambda x: x[1], wide.columns)
    for col in metas_to_filter:
        wide = wide.drop(col,1)

    wide = wide.iloc[1:-1,1:]
    merged = pd.merge(data[['pid','age','height','weight']],wide,right_index=True,left_on='pid').drop_duplicates()
    merged.index = merged['pid']
    merged = merged.drop('pid',1)
    merged = merged[~merged.index.duplicated(keep=False) | ~numpy.isnan(merged['age'])]
    for col in merged.columns[4:]:
        merged[col] = map(cleanFloat,merged[col])
        merged[col] = merged[col].astype('float')

    return merged
    #print data.shape
    #print data.columns

def getMergedTable():
    met = getMetabolites()
    print met.shape

    snps = getSNPs()
    for column in snps.columns:
        if (numpy.sum(map(lambda x: x==None, snps[column]))/float(snps.shape[0]) > 0.4):
            snps = snps.drop(column,1)
    print snps.shape
    table = pd.merge(met,snps,left_index=True,right_index=True)
    print table.shape

    minor = {}
    with open('../reference_files/Allele_class.csv','rU') as csvfile:
        reader = csv.DictReader(csvfile)
        row = next(reader)
        row = next(reader)
        minor = row
    for column in table.columns[met.shape[1]:]:
        table[column] = map(functools.partial(GT2num,minor_allele=minor[column]), table[column])

    table.to_csv('Metabolites_SNP_joinedTable.csv',header=True,index=True)
    return table


def main():
    #data = getMergedTable()
    data = pd.read_csv('Metabolites_SNP_joinedTable.csv',index_col=0)

    snp_cols = []
    for col in data.columns:
        if re.search(r'rs\d+',col) or re.search(r'Hcv\d+',col):
            snp_cols.append(col)

    pv_dict = {}
    dict_rows = []
    for col in data.columns:
        if not col in snp_cols+['age','height','weight']:
            #print col
            dict_rows.append(col)
            for snp in snp_cols:
                #print snp
                if not snp in pv_dict:
                    pv_dict[snp] = []

                curr = data[[col,snp]]
                curr[col].astype('float',error='coerce',copy=False)
                curr = curr[~numpy.isnan(curr[col]) & ~numpy.isnan(curr[snp])]
                #print curr
                curr['type'] = map(lambda x: x>0,curr[snp])
                a = curr[curr['type']][col]
                b = curr[~curr['type']][col]
                if len(a)>10 and len(b)>10:
                    print st.ranksums(a,b)[1]
                    pv_dict[snp].append(st.ranksums(a,b)[1])
                else:
                    pv_dict[snp].append(numpy.nan)
            
    pvalues = pd.DataFrame(pv_dict, index=dict_rows)
    #return pvalues
    #print pvalues

    from statsmodels.sandbox.stats.multicomp import multipletests
    pv = pvalues.values.reshape(pvalues.shape[0]*pvalues.shape[1])
    mask = numpy.isfinite(pv)
    pc = numpy.empty(pv.shape)
    pc.fill(numpy.nan)
    pc[mask] = multipletests(pv[mask],method='fdr_bh')[1]
    print pc[pc < 0.1]
    corrected_pvalues = pd.DataFrame(pc.reshape(pvalues.shape[0],pvalues.shape[1]),index=pvalues.index,columns=pvalues.columns)
    #print corrected_pvalues
    #print corrected_pvalues[corrected_pvalues < 0.1]
    corrected_pvalues.to_csv('Metabolites_SNP_pvalues.csv',index=True,header=True)


if __name__=='__main__':
    main()
