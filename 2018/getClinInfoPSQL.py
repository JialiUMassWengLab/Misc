#! /usr/bin/python

import re
import os
import sys
import psycopg2
import numpy as np
import pandas as pd

def connect():
    #su postgres
    #password: Welcome@ms
    #cd /usr/local/pgsql/data
    #psql
    conn = psycopg2.connect(dbname="test_db",user="postgres")
    return conn

def getClinicInfo(tbl,conn):
    df = pd.read_csv(tbl,sep='\t',index_col=0,low_memory=False)
    cols = set(map(lambda x: x.split('-')[0],df.columns[:-2]))
    alqIDs = ','.join([i for i in cols if not re.search(r'^\d',i)==None])
    clinInfo = {'PatientID':{},'Disease':{},'PathFibrosis':{},'APOE':{},'MMSE':{},'CDR':{},'Center':{}}
    cur = conn.cursor()
    cur.execute("SELECT tbaliquot.\"AliquotID\",tbclinicalhx.\"PatientID\",tbclinicalhx.\"Disease\",tbclinicalhx.\"PathFibrosis\",tbclinicalhx.\"TestAPOE\",tbclinicalhx.\"TestMMSE\",tbclinicalhx.\"TestCDR\",tbpracticeinfo.\"PracticeName\" FROM tbspecimeninfo \
    INNER JOIN tbaliquot ON (tbspecimeninfo.specimenindexid = tbaliquot.\"Parent\") \
    INNER JOIN tbclinicalhx ON (tbspecimeninfo.caseid = tbclinicalhx.\"CaseID\") \
    INNER JOIN tbcases ON (tbspecimeninfo.caseid = tbcases.\"ID\") \
    INNER JOIN tbpracticeinfo ON (tbcases.\"Practice\" = tbpracticeinfo.\"PracticeIndexID\") \
    WHERE tbaliquot.\"AliquotID\" IN (%s)" % (alqIDs))
    rows = cur.fetchall()
    for row in rows:
        aliquot = str(row[0]).replace('.0','')
        clinInfo['PatientID'].update({aliquot: str(row[1]).replace('.0','')})
        clinInfo['Disease'].update({aliquot: str(row[2])})
        clinInfo['PathFibrosis'].update({aliquot: str(row[3])})
        clinInfo['APOE'].update({aliquot: str(row[4])})
        clinInfo['MMSE'].update({aliquot: str(row[5])})
        clinInfo['CDR'].update({aliquot: str(row[6])})
        clinInfo['Center'].update({aliquot: str(row[7])})
    cur.close()
    del cur
    return pd.DataFrame.from_dict(clinInfo)


def getDemographicInfo(tbl,conn):
    df = pd.read_csv(tbl,sep='\t',index_col=0,low_memory=False)
    cols = set(map(lambda x: x.split('-')[0],df.columns[:-2]))
    alqIDs = ','.join([i for i in cols if not re.search(r'^\d',i)==None])
    clinInfo = {'BMI':{},'Ethnicity':{},'Gender':{}}
    cur = conn.cursor()
    cur.execute("SELECT tbaliquot.\"AliquotID\",tbclinicalhx.\"BMI\",tbclinicalhx.\"Ethnicity\",tbpatientinfo.\"PatientGender\" FROM tbspecimeninfo \
    INNER JOIN tbaliquot ON (tbspecimeninfo.specimenindexid = tbaliquot.\"Parent\") \
    INNER JOIN tbclinicalhx ON (tbspecimeninfo.caseid = tbclinicalhx.\"CaseID\") \
    INNER JOIN tbcases ON (tbspecimeninfo.caseid = tbcases.\"ID\") \
    INNER JOIN tbpracticeinfo ON (tbcases.\"Practice\" = tbpracticeinfo.\"PracticeIndexID\") \
    INNER JOIN tbpatientinfo ON (tbpatientinfo.\"PatientIndexID\" = tbclinicalhx.\"PatientID\") \
    WHERE tbaliquot.\"AliquotID\" IN (%s)" % (alqIDs))
    rows = cur.fetchall()
    for row in rows:
        aliquot = str(row[0]).replace('.0','')
        clinInfo['BMI'].update({aliquot: row[1]})
        clinInfo['Ethnicity'].update({aliquot: str(row[2])})
        clinInfo['Gender'].update({aliquot: str(row[3])})
    cur.close()
    del cur
    return pd.DataFrame.from_dict(clinInfo)


def main():
    conn = connect()
    demo = getDemographicInfo(sys.argv[1],conn)
    info = getClinicInfo(sys.argv[1],conn)
    info['Disease'] = map(lambda x: 'NCI' if x=='Normal Control' else x,info['Disease'])
    info['Disease'] = map(lambda x: 'AD' if x=='Alzheimer\'s' else x,info['Disease'])
    #info[['PatientID','Center','Disease','PathFibrosis']].to_csv('sampleInfo.csv')
    info = pd.concat([info,demo],join='inner',axis=1)
    info[['PatientID','BMI','Ethnicity','Gender','Center','Disease','APOE','MMSE','CDR']].to_csv('sampleInfo2.csv')
    

if __name__=='__main__':
    main()
