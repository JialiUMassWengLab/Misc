#! /usr/bin/python

import psycopg2
import os
import re
import sys
import glob

def connect():
    conn = psycopg2.connect(dbname="test_db",user="postgres")
    return conn

def getDict(fn):
    data = {}
    with open(fn,'rU') as infile:
        for line in infile:
            a = line.split()
            if re.search(r'\.',a[2]):
                data[a[0]] = float(a[2])
            else:
                data[a[0]] = int(a[2])

    return data
    
def getIndex(id_field,db,conn):
    cur = conn.cursor()
    cur.execute("SELECT DISTINCT %s FROM %s" % (id_field,db))
    rows = cur.fetchall()
    nums = [0]
    for row in rows:
        nums.append(re.sub(r'^[A-Za-z]+0+','',row[0]))

    cur.close()
    del cur
    return max(map(int,nums))


def populateGeneTbl(conn):
    names = {}
    annos = {}
    with open('/mnt/shares2/annotations/hg38/gencode_basic_gene_names.tsv','rU') as infile:
        cur = conn.cursor()
        added = {}
        id_repeat = set([])
        for line in infile:
            a = line.split()
            if not a[0] in id_repeat:
                if a[1] in added:
                    added[a[1]] += 1
                    a[1] += '_'+str(added[a[1]])
                else:
                    added[a[1]] = 1

                #a[1] = a[1].replace('.','_')
                #a[1] = a[1].replace('-','_')
                names[a[0]] = a[1]
                annos[a[0]] = a[2]
                print a[0],a[1],a[2]
                cur.execute("""INSERT INTO gene_info(ensembl_id, common_name, gencode_anno) VALUES ('%s', '%s', '%s')""" % (a[0],a[1],a[2]))
                id_repeat.add(a[0])

        cur.close()
        conn.commit()


def populateExprTbl(conn,sampleID,runID,prefix):
    rawCounts = dedupCounts = {}
    PCrawTPM = PCdedupTPM = {}
    rawCounts = getDict(prefix+'.gtf24.htseq-count.reads.tab')
    dedupCounts = getDict(prefix+'.dedup.gtf24.htseq-count.reads.tab')
    PCrawTPM = getDict(prefix+'.gtf24.htseq-count.proteinCoding.TPM.tab')
    PCdedupTPM = getDict(prefix+'.dedup.gtf24.htseq-count.proteinCoding.TPM.tab')
    if rawCounts and dedupCounts and PCrawTPM and PCdedupTPM:
        cur = conn.cursor()
        cur.execute("""INSERT INTO sample_data(sample_id,sample_name,run_id) VALUES ('%s','%s','%s')""" \
                    % (sampleID,os.path.basename(prefix),runID))

        for gene in rawCounts.keys():
            if gene in PCrawTPM and gene in PCdedupTPM:
                cur.execute("""INSERT INTO expression(sample_id,ensembl_id,counts_raw,counts_dedup,pc_tpm_raw,pc_tpm_dedup) VALUES ('%s','%s','%d','%d','%f','%f')""" % (sampleID,gene,rawCounts[gene],dedupCounts[gene],PCrawTPM[gene],PCdedupTPM[gene]))
            else:
                cur.execute("""INSERT INTO expression(sample_id,ensembl_id,counts_raw,counts_dedup) VALUES ('%s','%s','%d','%d')""" % (sampleID,gene,rawCounts[gene],dedupCounts[gene]))  
        
        cur.close()
        conn.commit()
                

def populateRunTbl(conn,runID,runName):
    cur = conn.cursor()
    cur.execute("""INSERT INTO run_data(run_id,run_name) VALUES ('%s','%s')""" % (runID,runName))
    cur.close()
    conn.commit()


def main():
    conn = connect()
    #populateGeneTbl(conn)

    sample_index,run_index,experiment_index = getIndex('sample_id','sample_data',conn),getIndex('run_id','run_data',conn),getIndex('exp_id','experiment_data',conn)
    runName = os.path.basename(sys.argv[1])
    if runName == '':
        runName = os.path.basename(os.path.dirname(sys.argv[1]))
    run_index += 1
    runID = ''.join(['Run']+['0']*(6-len(str(run_index)))+[str(run_index)])
    populateRunTbl(conn,runID,runName)

    files = glob.glob(os.path.join(sys.argv[1],'*.dedup.htseq.out'))
    for fn in files:
        fn = fn.replace('.dedup.htseq.out','')
        print fn
        sampleName = os.path.basename(fn)
        sample_index += 1
        sampleID = ''.join(['RD']+['0']*(7-len(str(sample_index)))+[str(sample_index)])
        print sampleID
        populateExprTbl(conn,sampleID,runID,fn)
    conn.close()


if __name__=='__main__':
    main()

