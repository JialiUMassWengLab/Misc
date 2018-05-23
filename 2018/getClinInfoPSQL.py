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
    clinInfo = {'PatientID':{},'Disease':{},'PathFibrosis':{},'PathSteatosis':{},'PathInflammation':{},'PathBallooning':{},'Center':{}}
    cur = conn.cursor()
    cur.execute("SELECT tbaliquot.\"AliquotID\",tbclinicalhx.\"PatientID\",tbclinicalhx.\"Disease\",tbclinicalhx.\"PathFibrosis\",tbclinicalhx.\"PathSteatosis\",tbclinicalhx.\"PathInflammation\",tbclinicalhx.\"PathBallooning\",tbpracticeinfo.\"PracticeName\" FROM tbspecimeninfo \
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
        clinInfo['PathSteatosis'].update({aliquot: str(row[4])})
        clinInfo['PathInflammation'].update({aliquot: str(row[5])})
        clinInfo['PathBallooning'].update({aliquot: str(row[6])})
        clinInfo['Center'].update({aliquot: str(row[7])})
    cur.close()
    del cur
    return pd.DataFrame.from_dict(clinInfo)


def main():
    conn = connect()
    info = getClinicInfo(sys.argv[1],conn)
    info[['PatientID','Center','Disease','PathFibrosis']].to_csv('sampleInfo.csv')
    

if __name__=='__main__':
    main()
