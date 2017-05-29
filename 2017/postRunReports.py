#! /usr/bin/python

import sys
import os
import re
import json
import glob
import requests
import subprocess
import renderPage
from lxml import etree, html

def pull(outFileName):
    r = requests.get('http://10.0.1.9:8090/rest/api/content/3211305?expand=body.storage',auth=('jzhuang','Xiamen533211@@'))
    j = json.loads(r.text)
    document_root = html.fromstring(j['body']['storage']['value'])
    #print etree.tostring(document_root, encoding='unicode', pretty_print=True).encode("UTF-8")
    f = open(outFileName,'wb')
    f.write(etree.tostring(document_root, encoding='unicode', pretty_print=True).encode("UTF-8"))
    #f.write(etree.tostring(document_root, pretty_print=True))
    f.close()

def postSampleReport(sourceDir,sampleName,parentID,auth):
    s = renderPage.renderSamplePage(sourceDir,sampleName)
    #b = os.path.basename(sourceDir).split('_')
    #title = '_'.join(b[1:]+[b[0]])
    title = os.path.basename(sourceDir)
    if title == '':
        title = os.path.basename(os.path.dirname(sourceDir))

    pageData = {"type":"page",
                "title":title+'__'+sampleName,
                "space":{"key":"RR"},
                "ancestors":[{"id":parentID}],
                "body":{"storage":{"value":s,"representation":"storage"}},
    }
    #print pageData
    r = requests.post('http://10.0.1.9:8090/rest/api/content/',data=json.dumps(pageData),auth=auth,headers=({'Content-Type':'application/json'}))
    if 'message' in r.json():
        print r.json()['message']
        if re.search(r'A page with this title already exists',r.json()['message']):
            pageTitle = re.sub(r' ','%20',pageData['title'])
            r1 = requests.get('http://10.0.1.9:8090/rest/api/content?title=%s' % pageTitle,auth=auth)
            pageID = str(r1.json()['results'][0]['id'])
            version = r1.json()['results'][0]['_expandable']['version']
            if version:
                version = int(version)+1
            else:
                version = 2
            pageData.update({"version":{"number": version}})
            r = requests.put('http://10.0.1.9:8090/rest/api/content/'+str(r1.json()['results'][0]['id']),data=json.dumps(pageData),auth=auth,headers=({'Content-Type':'application/json'}))
    return r.json()['id']

def postRunReport(sourceDir,auth):
    s = renderPage.renderRunPage(sourceDir)
    #b = os.path.basename(sourceDir).split('_')
    #title = '_'.join(b[1:]+[b[0]])
    title = os.path.basename(sourceDir)
    if title == '':
        title = os.path.basename(os.path.dirname(sourceDir))

    pageData = {"type":"page",
                "title":title+'__Run Summary',
                "space":{"key":"RR"},
                "ancestors":[{"id":"3211370"}],
                "body":{"storage":{"value":s,"representation":"storage"}},
    }
    r = requests.post('http://10.0.1.9:8090/rest/api/content/',data=json.dumps(pageData),auth=auth,headers=({'Content-Type':'application/json'}))
    if 'message' in r.json():
        if re.search(r'A page with this title already exists',r.json()['message']):
            pageTitle = re.sub(r' ','%20',pageData['title'])
            r1 = requests.get('http://10.0.1.9:8090/rest/api/content?title=%s' % pageTitle,auth=auth)
            pageID = str(r1.json()['results'][0]['id'])
            version = r1.json()['results'][0]['_expandable']['version']
            if version:
                version = int(version)+1
            else:
                version = 2
            pageData.update({"version":{"number": version}})
            r = requests.put('http://10.0.1.9:8090/rest/api/content/'+pageID,data=json.dumps(pageData),auth=auth,headers=({'Content-Type':'application/json'}))
    return r.json()['id']

def postChildPage(html,parentID,auth):
    with open(html,'rb') as htmlFile:
        s = htmlFile.read()
        print s
        pageData = {"type":"page",
                    "title":html,
                    "space":{"key":"BR"},
                    "ancestors":[{"id":parentID}],
                    "body":{"storage":{"value":s,"representation":"storage"}},
                }
        print pageData
        r = requests.post('http://10.0.1.9:8090/rest/api/content/',data=json.dumps(pageData),auth=auth,headers=({'Content-Type':'application/json'}))
        #print r.json()
        print r.json()['id']

def attach(prefix,pageID,auth):
    url = 'http://10.0.1.9:8090/rest/api/content/%s/child/attachment' % pageID
    pdf_files = glob.glob(prefix+'[\._]*.pdf')
    png_files = glob.glob(prefix+'[\._]*.png')
    tsv_files = glob.glob(prefix+'[\._]*.tsv')
    xlsx_files = glob.glob(prefix+'[\._]*.xlsx')
    rsem_files = glob.glob(prefix+'[\._]*.results')
    html_files = glob.glob(os.path.join(prefix+'_fastqc','*.html'))
    #fpkm_files = glob.glob(os.path.join(prefix+'_cufflinks_out','Protein_coding*_tracking'))

    fileDict = {'pdf': pdf_files,
                'png': png_files,
                'tsv': tsv_files,
                'xlsx': xlsx_files,
                'html': html_files,
                'rsem': rsem_files,
                #'fkpm': fpkm_files,
            }

    for fileType in fileDict.values():
        for File in fileType:
            files = {'file': open(File,'rb'),
                     'comment': 'this is %s' % os.path.basename(File),
                 }
            r = requests.post(url,files=files,headers=({'X-Atlassian-Token': 'no-check'}),auth=auth)
            #print r.json()

def put():
    url = 'http://10.0.1.9:8090/rest/api/content/3211305'
    with open('samplePage.html','rU') as template:
        s = template.read()
        print s

    pageData = {"id": "3211305",
                "type":"page",
                "title":"Sample page mirror",
                "space":{"key":"BR"},
                "body":{"storage":{"value":s,"representation":"storage"}},
                "version":{"number":6},
            }
    r = requests.put(url,data=json.dumps(pageData),auth=('jzhuang','Xiamen533211!!'),headers=({'Content-Type':'application/json'}))
    print r.text

def deleteRun(runName,auth):
    url = 'http://10.0.1.9:8090/rest/api/content/search?cql=space=RR+and+type=page+and+title~"%s"' % runName
    r = requests.get(url,auth=auth)
    for result in r.json()['results']:
        pageID = result['id']
        print pageID
        url1 = 'http://10.0.1.9:8090/rest/api/content/%s' % pageID
        r1 = requests.delete(url1,auth=auth)

def main():
    if sys.argv[1]=='del':
        auth = (sys.argv[3],sys.argv[4])
        deleteRun(sys.argv[2],auth)
    else:
        runName = os.path.basename(sys.argv[1])
        if runName == '':
            runName = os.path.basename(os.path.dirname(sys.argv[1]))

        print runName
        auth = (sys.argv[2],sys.argv[3])
        subprocess.check_output('bash run_samples_summary_stats.bash %s' % sys.argv[1],shell=True,stderr=subprocess.STDOUT)
        #subprocess.check_output('python correlation.py %s %s' % (sys.argv[1],runName+'.exprPCC'),shell=True,stderr=subprocess.STDOUT)
        subprocess.check_output('python countGenesDetected.py %s' % sys.argv[1],shell=True,stderr=subprocess.STDOUT)
        subprocess.check_output('Rscript base_mapping_distribution.R %s' % sys.argv[1],shell=True,stderr=subprocess.STDOUT)
        runPageID = postRunReport(sys.argv[1],auth)
        print runPageID
        attach(os.path.join(sys.argv[1],runName),runPageID,auth)

        samples = glob.glob(os.path.join(sys.argv[1],'*Aligned.sortedByCoord.out.bam.stats.out'))
        for sample in samples:
            name = re.sub(r'Aligned.sortedByCoord.out.bam.stats.out','',os.path.basename(sample))
            pageID = postSampleReport(sys.argv[1],name,runPageID,auth)
            print pageID
            attach(os.path.join(sys.argv[1],name),pageID,auth)
            #postChildPage('/mnt/shares2/analysis/161110_TH024_JRP050/Shalfr2_fastqc/Shalfr2_R1_val_1_fastqc.html','3211336')


if __name__=='__main__':
    main()
