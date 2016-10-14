#! /usr/bin/python 

import sys
import csv
import re
import pandas as pd
import datetime as dt

def getClinics():
    outfile = open('ClinicsNonOpioid.txt','w')
    with open('reference_files/Non-Opioid Response.csv','rU') as infile:
        rows = csv.reader(infile,delimiter=',')
        for row in rows:
            visit = 0
            if row[13] == '2': 
                visit = 2
            elif row[14] == '3':
                visit = 3

            if visit == 0 or row[9]=='pid': 
                continue

            newrow = [row[9].upper().strip(),visit,row[10].upper().strip(),row[11]]
            writer = csv.writer(outfile,delimiter=',')
            writer.writerow(newrow)

    outfile.close()

def getDemo():
    demo = pd.read_csv("patientInfo.csv",header=None,names=['pid','gender','race','dob'])
    demo['dob'] = pd.to_datetime(demo['dob'])
    demo['age'] = (pd.datetime.now().date()-demo['dob']).astype('timedelta64[Y]')
    demo = demo.drop_duplicates()
    demo = demo[~demo['pid'].duplicated()]
    demo.index = demo['pid']
    demo.drop('pid',1,inplace=True)
    return demo

def getJoinedDemo():
    demo = getDemo()
    clinics = pd.read_csv('ClinicsNonOpioid.txt',header=None,names=['pid','visit','clinic','date'])
    clinics['date'] = pd.to_datetime(clinics['date'])
    clinics = clinics[~clinics[['pid','visit']].duplicated()]
    clinics = clinics.pivot(index='pid',columns='visit')
    clinics.columns = map(lambda x: x[0]+str(x[1]), clinics.columns)

    NOpatients = []
    with open('Non-Opioid-physician-rating.csv','rU') as infile:
        rows = csv.reader(infile,delimiter=',')
        for row in rows:
            row[0] = row[0].strip()
            if not row[0].upper() in NOpatients:
                NOpatients.append(row[0].upper())

    merged = pd.merge(demo,clinics,left_index=True,right_index=True)
    merged = merged[merged.index.isin(NOpatients)]

    labels = ['{0} - {1}'.format(i,i+4) for i in range(0,125,5)]
    merged['ageGroup'] = pd.cut(merged['age'], range(0,126,5), right=False, labels=labels)
    #print merged
    print merged.shape
    return merged

def main():
    getClinics()
    merged = getJoinedDemo()
    print merged['gender'].value_counts()
    print merged['race'].value_counts()
    print merged['ageGroup'].value_counts()
    print merged['age'].min()
    print merged['age'].max()
    print merged[(merged['date2']>'2013-01-01') & (merged['date2']<'2016-10-01')]['date2'].min()
    print merged[(merged['date2']>'2013-01-01') & (merged['date2']<'2016-10-01')]['date2'].max()
    print merged[(merged['date3']>'2013-01-01') & (merged['date3']<'2016-10-01')]['date3'].min()
    print merged[(merged['date3']>'2013-01-01') & (merged['date3']<'2016-10-01')]['date3'].max()

    clinics = []
    for key,value in merged['clinic2'].value_counts().iteritems():
        if value > 5 and not key in clinics:
            clinics.append(key)
    for key,value in merged['clinic3'].value_counts().iteritems():
        if value > 5 and not key in clinics:
            clinics.append(key)
    print len(clinics)


if __name__ == "__main__":
    main()
