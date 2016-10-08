#! /usr/bin/python

import pandas as pd
import numpy
import scipy
import sys
import csv
import re

drug_names = {
    'gabapentin': ['GABAPENTIN','FANATREX','GRAILSE','NEUROTIN','NUPENTIN','NEOGAB'],
    'ibuprofen': ['IBUPROFEN','ADVIL','BRUFEN','MORTIN','NUROFEN'],
    'duloxetine': ['DELOXETINE','CYMBALTA'],
    'acetaminophen': ['ACETAMINOPHEN','TYLENOL','PANADOL'],
    'alprazolam': ['ALPRAZOLAM','XANAX']
}

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

def resultCategory(a,b):
    if not len(a) == len(b):
        print "Length not equal!"
        return None

    goodResponder = {'Good':0, 'Average':0, 'Poor':0}
    badResponder = {'Good':0, 'Average':0, 'Poor':0}
    for i in range(len(a)):
        if a[i]:
            if b[i] < 2: 
                goodResponder['Poor'] += 1
            elif b[i] < 4:
                goodResponder['Average'] += 1
            else:
                goodResponder['Good'] += 1
        else:
            if b[i] < 2: 
                badResponder['Poor'] += 1
            elif b[i] < 4:
                badResponder['Average'] += 1
            else:
                badResponder['Good'] += 1
    return (goodResponder, badResponder)            

def resultCategory1(a,b,le):
    a = le.inverse_transform(a)
    if not len(a) == len(b):
        print "Length not equal!"
        return None

    goodResponder = {'Good':0, 'Average':0, 'Poor':0}
    badResponder = {'Good':0, 'Average':0, 'Poor':0}
    for i in range(len(a)):
        if a[i] == 'Good':
            goodResponder[b[i]] += 1
        else:
            badResponder[b[i]] += 1

    return (goodResponder, badResponder)

def leafCategories(dt):
    cats = []
    for i in range(len(list(dt.tree_.value))):
        if dt.tree_.feature[i]==-2:
            cats.append(numpy.argmax(dt.tree_.value[i][0]))
    return cats

def scorer(estimator, X, y):
    good, bad = resultCategory(estimator.predict(X), y)

    if sum(good.values()) > 0 and sum(bad.values()) > 0:
        obs = numpy.array([[good['Poor'], good['Average'], good['Good']],
                           [bad['Poor'], bad['Average'], bad['Good']]])
        return scipy.stats.chi2_contingency(obs, correction=False)[1]
    else:
        return 1

    #return sum(bad.values())

def getMED(drug_list):
    MED = pd.read_csv('all_r_m_aug1.MED.csv',header=None,names=['pid','visit','drug','p.rate'])
    MED['p.rate'] = pd.to_numeric(MED['p.rate'], errors='coerce')
    MED = MED[~numpy.isnan(MED['p.rate'])]
    MED_filter = MED[map(lambda x: not re.search(drug_list[0],x)==None, MED['drug'])][['pid','visit','p.rate']]
    for i in range(1,len(drug_list)):
        MED_filter = MED_filter.append(MED[map(lambda x: not re.search(drug_list[i],x)==None, MED['drug'])][['pid','visit','p.rate']])
    MED_filter = MED_filter.sort_values(by=['pid','p.rate'],ascending=[True,False],axis=0)
    MED_filter = MED_filter.drop_duplicates()
    MED_filter = MED_filter[(MED_filter['p.rate'] <= 5) & (MED_filter['p.rate'] >= 0)]
    #print MED_filter

    resp2 = MED_filter.loc[MED_filter.visit=='2',['pid','p.rate']]
    resp2 = resp2[~resp2.pid.duplicated()]
    resp1 = MED_filter.loc[(MED_filter.visit=='1') & (~MED_filter.pid.isin(resp2.pid)), ['pid','p.rate']]
    resp1 = resp1[~resp1.pid.duplicated()]
    resp3 = MED_filter.loc[(MED_filter.visit=='3') & (~MED_filter.pid.isin(resp2.pid))
                           & (~MED_filter.pid.isin(resp1.pid)), ['pid','p.rate']]
    resp3 = resp3[~resp3.pid.duplicated()]
    resp4 = MED_filter.loc[(MED_filter.visit=='4') & (~MED_filter.pid.isin(resp2.pid))
                           & (~MED_filter.pid.isin(resp1.pid)) & (~MED_filter.pid.isin(resp3.pid)), ['pid','p.rate']]
    resp4 = resp4[~resp4.pid.duplicated()]

    #print resp2.shape,resp1.shape,resp3.shape,resp4.shape
    merged = pd.concat([resp2,resp1,resp3,resp4])
    merged['grade'] = pd.cut(merged['p.rate'], 3, labels=['Poor','Average','Good']).astype('category')
    return merged

def getReport(drug):
    data = pd.read_csv(drug+"Response.csv",header=None,names=['sid','pid','Rank','maxRank','isGoodResponder','metabolism'])
    data.metabolism = pd.to_numeric(data.metabolism, errors="coerce")
    data = data[~numpy.isnan(data.metabolism)]
    data = data[data['metabolism'] > 0]
    data['metabolism'] = data['metabolism'].astype('category', categories=[1,2,3])
    #print data.metabolism.unique()
    data['metabolism'].cat.categories = ['High','Normal','Low']
    return data
    
def getDemo():
    demo = pd.read_csv("reference_files/patientInfo.csv", header=None, names=['pid','gender','race'])
    demo.gender = demo.gender.astype('category')
    demo.gender.cat.categories = [0,1]
    demo.race = demo.race.astype('category')
    return demo

def getSNPs():
    data = pd.read_csv("all_SNP_calls.csv",header=None,names=['rs','gt','sid','pid'])
    data = data[map(lambda x: not re.search(r'[\w-]+/[\w-]+',x)==None, data['gt'])][['rs','gt','pid']]
    data = data[~data[['pid','rs']].duplicated(keep=False)]
    data = data.pivot(index='pid',columns='rs')
    data.columns = map(lambda x: x[1], data.columns)
    return data

def logisticMED_vs_report(drug):
    drug_list = drug_names[drug]

    resp = getReport(drug)
    demo = getDemo()
    MED = getMED(drug_list)
    #print MED

    merged = pd.merge(demo,MED,on='pid')
    merged = pd.merge(merged,resp[['pid','isGoodResponder','metabolism']].drop_duplicates(),on='pid')
    merged = merged.sort_values(by='pid',ascending=True,axis=0)
    merged.index = merged['pid']
    merged = merged.drop('pid',1)
    #print merged

    cats = pd.get_dummies(merged)
    #from sklearn.feature_selection import chi2
    #scores, pvalues = chi2(processed.filter(regex=('race_*|gender|isGoodResponder|metabolism')).values, processed.filter(regex=('grade_*')).values)
    #print pvalues
    import statsmodels.api as sm
    x = cats.filter(items=['gender','isGoodResponder','metabolism_Low','metabolism_High','race_W','race_B/AA','race_H/L']).astype('float')
    x = sm.add_constant(x, prepend = False)
    logit = sm.MNLogit(merged['grade'], x)
    print logit.fit(maxiter=1000,method='powell').summary()

def getJoinedGenotypeTable(drug):
    import functools

    MED = getMED(drug_names[drug])
    snps = getSNPs()
    snps = pd.merge(MED,snps,left_on='pid',right_index=True)
    snps.index = snps['pid']
    snps = snps.drop('pid',1)
    print(snps.shape)
    for column in snps.columns:
        if (numpy.sum(map(lambda x: x==None, snps[column]))/float(snps.shape[0]) > 0.4):
            snps = snps.drop(column,1)
    for row in snps.index:
        if (numpy.sum(map(lambda x: x==None, snps.loc[row]))/float(snps.shape[1]-2) > 0.5):
            snps = snps.drop(row,0)

    minor = {}
    with open('reference_files/Allele_class.csv','rU') as csvfile:
        reader = csv.DictReader(csvfile)
        row = next(reader)
        row = next(reader)
        minor = row

    for column in snps.columns[2:]:
        snps[column] = map(functools.partial(GT2num,minor_allele=minor[column]), snps[column])

    demo = getDemo()
    snps = pd.merge(snps,demo,left_index=True,right_on='pid')
    snps.index = snps['pid']
    snps = snps.drop('pid',1)
    snps = pd.concat([snps.iloc[:,:-1], pd.get_dummies(snps.iloc[:,-1:]).filter(items=['race_W','race_B/AA','race_H/L','race_A'])], axis=1)
    snps['grade1'] = map(lambda x: x>=3.0, snps['p.rate'])
    print(snps.shape)
    return snps
    '''
    snps = pd.concat([snps.iloc[:,0:2],pd.get_dummies(snps.iloc[:,2:len(snps.columns)])],axis=1)
    
    import sklearn.tree as sktr
    dt = sktr.DecisionTreeClassifier(min_samples_split=20, random_state=99)
    dt.fit(snps.iloc[:,2:len(snps.columns)], snps['p.rate'])
    print dt.predict(snps.iloc[30:80,2:len(snps.columns)])
    '''

def learningMED(drug):
    import sklearn.preprocessing as skpre
    import sklearn.ensemble as sken
    from sklearn import cross_validation
    from sklearn import utils

    snps = getJoinedGenotypeTable(drug)
    #X_train,X_test,Y_train,Y_test = cross_validation.train_test_split(skpre.Imputer(strategy='most_frequent',axis=1).fit_transform(snps.iloc[:,2:]),
    #                                                                  snps['p.rate'],test_size=0.4,random_state=199)
    X_train,X_test,Y_train,Y_test = cross_validation.train_test_split(skpre.Imputer(strategy='most_frequent',axis=1).fit_transform(snps.iloc[:,2:]),
                                                                      snps['grade1'],test_size=0.4,random_state=199)
    X_train1,X_test1,Y_train1,Y_test1 = cross_validation.train_test_split(skpre.Imputer(strategy='most_frequent',axis=1).fit_transform(snps.iloc[:,2:]),
                                                                          snps['p.rate'],test_size=0.4,random_state=199)
    '''
    train_df = pd.concat([pd.DataFrame(X_train,columns=snps.columns[2:]), pd.Series(Y_train,name='p.rate')], axis=1)
    target = min(train_df['p.rate'].value_counts())
    train_resamp=[]
    for i in range(6):
        if len(train_df[train_df['p.rate']==i]) > target:
            train_resamp.append(utils.resample(train_df[train_df['p.rate'] == i], replace=True, n_samples=target, random_state=73))
        else:
            train_resamp.append(numpy.array(train_df[train_df['p.rate'] == i]))
    train_resamp = numpy.array(train_resamp)
    print train_resamp.shape
    train_resamp = numpy.reshape(train_resamp,(train_resamp.shape[0]*train_resamp.shape[1],train_resamp.shape[2]))
    X_train = train_resamp[:,:-1]
    Y_train = train_resamp[:,-1]
    '''
    ntrees = [10,15,20,25,30,40]
    leafSampleSizes = [2,5,10,15,20,25,30,40]
    depths = [8,10,12,15,20]

    for n in ntrees:
        for m in leafSampleSizes:
        #for d in depths:
            rf = sken.RandomForestClassifier(min_samples_split=m,n_estimators=n,random_state=99)
            #rf = sken.RandomForestClassifier(max_depth=d,n_estimators=n,random_state=99)
            #rf = sken.RandomForestRegressor(min_samples_split=m,n_estimators=n,random_state=99)
            scores = cross_validation.cross_val_score(rf,X_train,Y_train,cv=5)
            print n,m,scores

    rf = sken.RandomForestClassifier(min_samples_split=15,n_estimators=25,random_state=119)
    #rf = sken.RandomForestClassifier(max_depth=12,n_estimators=40,random_state=119)
    rf.fit(X_train,Y_train)
    cats = resultCategory(rf.predict(X_test),Y_test1)
    #print list(rf.predict(X_test))
    #print list(Y_test)
    print cats
    print scorer(rf,X_test,Y_test1)
    '''
    le = skpre.LabelEncoder()
    le.fit(Y_train)
    for n in ntrees:
        for m in leafSampleSizes:            
            rf = sken.RandomForestClassifier(min_samples_split=m,n_estimators=n,random_state=99)
            scores = cross_validation.cross_val_score(rf,X_train,le.transform(Y_train),cv=5)
            print n,m,scores

    rf = sken.RandomForestClassifier(min_samples_split=25,n_estimators=20,random_state=99)
    rf.fit(X_train,le.transform(Y_train))
    #print list(rf.predict(X_test))
    #print list(Y_test)
    cats = resultCategory1(rf.predict(X_test),Y_test,le)
    print cats

    nest = [100,150,200,250,300,350,400]
    maxdepth = [3,4,5,6]
    for n in nest:
        for d in maxdepth:
            gb = sken.GradientBoostingClassifier(random_state=91,max_depth=d,n_estimators=n)
            scores = cross_validation.cross_val_score(gb,X_train,Y_train,cv=5,scoring='mean_squared_error')
            print n,d,scores

    gb = sken.GradientBoostingClassifier(random_state=91,max_depth=4,n_estimators=200)
    gb.fit(X_train,Y_train)
    print list(gb.predict(X_test))
    cats = resultCategory(gb.predict(X_test),Y_test)
    print cats
    '''

def main():
    #logisticMED_vs_report(sys.argv[1])
    learningMED(sys.argv[1])

if __name__=="__main__":
    main()
