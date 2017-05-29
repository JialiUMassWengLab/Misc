#! /usr/bin/python

import os
import sys
import re
import json
import glob
import numpy
import jinja2
#from lxml import etree, html
#from django.template import Context, Template, loader

def loadStatsData(sourceDir,sampleName):
    picardStats = []
    with open(os.path.join(sourceDir,sampleName+'Aligned.sortedByCoord.out.bam.picard.stats.out'),'rU') as infile:
        for line in infile:
            a = line.split()
            if len(a)>1:
                '''
                if not re.search(r'\.',a[1]):
                    a[1] = int(a[1])
                else:
                    a[1] = float(a[1])
                '''
                picardStats.append({'measure':a[0],
                                    'value':a[1]
                                })

    bamStats = []
    with open(os.path.join(sourceDir,sampleName+'Aligned.sortedByCoord.out.bam.stats.out'),'rU') as infile:
        for line in infile:
            if line=='\n' or re.search(r'^#',line):
                continue
            a = line.split(':')
            if len(a)>1:
                bamStats.append({'measure':a[0],
                                 'value':int(a[1].strip())
                             })
            elif re.search(r'Non primary hits',line):
                bamStats.append({'measure':'Non primary hits',
                                 'value':int(re.sub(r'Non primary hits','',line).strip())
                             })

    readDist = {}
    with open(os.path.join(sourceDir,sampleName+'Aligned.sortedByCoord.out.bam.readDist.out'),'rU') as infile:
        lines = infile.readlines()
        readDist['total_reads'] = int(re.sub(r'Total Reads','',lines[0]).strip())
        readDist['total_tags'] = int(re.sub(r'Total Tags','',lines[1]).strip())
        readDist['total_assigned_tags'] = int(re.sub(r'Total Assigned Tags','',lines[2]).strip())
        
        readDist['stats'] = []
        for line in lines[5:15]:
            a = line.split()
            readDist['stats'].append({'group':a[0],
                                      'base':int(a[1]),
                                      'tag_count':int(a[2]),
                                      'tag_density':float(a[3]),
                                  })
            
    starStats = []
    with open(os.path.join(sourceDir,sampleName+'Log.final.out'),'rU') as infile:
        for i in range(5):
            next(infile)
        for line in infile:
            if re.search(r'\|',line):
                b = line.split('|')
                measure = b[0].strip()
                value = b[1].strip()
                value = value.strip('%')
                starStats.append({'measure': measure,
                                  'value': value,
                              })
        
    rsemStats = []
    with open(os.path.join(sourceDir,sampleName+'_geneTypes.tsv'),'rU') as infile:
        next(infile)
        for line in infile:
            rsemStats.append(line[:-1].split('\t'))

    geneCounts = []
    with open(os.path.join(sourceDir,sampleName+'.genes.detected.tsv'),'rU') as infile:
        for line in infile:
            a = line.split()
            geneCounts.append({'TPM':a[0], 'counts':a[1]})

    
    trace = {}
    with open(os.path.join(sourceDir,sampleName+'.ERCC.comp.txt'),'rU') as infile:
        x = []
        y = []
        text = []
        for line in infile:
            a = line.split(',')
            text.append(a[0])
            y.append(numpy.log10(float(a[1])+0.1))
            x.append(numpy.log10(float(a[2])))
            
        trace.update({'x': x, 'y': y, 'name': 'ERCC', 'type': 'scatter',
                      'text': text, 'mode': 'markers', 'marker': {'size': 10},
                      'hoverinfo': 'x+y+text',
                  })
    ERCCjson = json.dumps([trace])


    graphs = []
    #graphs = {'pdf':[], 'png':[]}
    for graph in glob.glob(os.path.join(sourceDir,sampleName+'[\._]*.pdf')):
        if not re.search(r'_RPKM_saturation',graph):
            #graphs['pdf'].append(os.path.basename(graph))
            graphs.append(os.path.basename(graph))            
    for graph in glob.glob(os.path.join(sourceDir,sampleName+'[\._]*.png')):
        #graphs['png'].append(os.path.basename(graph))
        graphs.append(os.path.basename(graph))
    tables = []
    for table in glob.glob(os.path.join(sourceDir,sampleName+'[\._]*.tsv')):
        tables.append(os.path.basename(table))
    for table in glob.glob(os.path.join(sourceDir,sampleName+'[\._]*.results')):
        tables.append(os.path.basename(table))

    contextDict = {'picardStats': picardStats,
                   'bamStats': bamStats,
                   'readDist': readDist,
                   'starStats': starStats,
                   'rsemStats': rsemStats,
                   'geneCounts': geneCounts,
                   'ERCCjson': ERCCjson,
                   'graphs': graphs,
                   'tables': tables,
               }
    return contextDict


def loadStarSummary(sourceDir,sampleName):
    starSummary = {'Total': 0, 'Unmapped': 0.0, 'Chimeric': 0.0, 'Unannotated': 0.0}
    with open(os.path.join(sourceDir,sampleName+'Log.final.out'),'rU') as infile:
        for i in range(5):
            next(infile)
        for line in infile:
            if re.search(r'\|',line):
                b = line.split('|')
                measure = b[0].strip()
                value = b[1].strip()
                value = value.strip('%')
                if measure=='Number of input reads':
                    starSummary['Total'] = int(value)
                elif re.search(r'% of reads unmapped',measure):
                    starSummary['Unmapped'] += float(value)
                elif measure=='% of chimeric reads':
                    starSummary['Chimeric'] = float(value)

    return starSummary    


def loadCompStats(sourceDir):
    files = glob.glob(os.path.join(sourceDir,'*_geneTypes.tsv'))
    categories = ['protein_coding','rRNA','pseudogene','misc_RNA','Mt_protein_coding','Mt_rRNA','Mt_tRNA','other']
    readCount = {}
    readPct = {}
    tpmPct = {}
    samples = []
    for cat in categories:
        readCount[cat] = []
        readPct[cat] = []
        tpmPct[cat] = []
    for cat in ['Chimeric','Unannotated','Unmapped']:
        readCount[cat] = []

    for fn in files:
        sample = re.sub(r'_geneTypes.tsv','',os.path.basename(fn))
        samples.append(sample)
        starSummary = loadStarSummary(sourceDir,sample)
        annotated = 0.0
        with open(fn,'rU') as infile:
            next(infile)
            readCountOther = 0.0
            readOther = 0.0
            tpmOther = 0.0
            for line in infile:
                a = line.split('\t')
                count = float(a[1].split('(')[0].strip())
                match1 = re.search(r'\(([\d\.]+)%\)',a[1])
                match2 = re.search(r'\(([\d\.]+)%\)',a[2])
                if match1 and match2:
                    value = count/float(starSummary['Total'])*100
                    if a[0] in categories:
                        readCount[a[0]].append( '{0:.2f}'.format(value) )
                        readPct[a[0]].append(match1.group(1))
                        tpmPct[a[0]].append(match2.group(1))
                    else:
                        readCountOther += value
                        readOther += float(match1.group(1))
                        tpmOther += float(match2.group(1))

                    annotated += value

            readCount['other'].append( '{0:.2f}'.format(readCountOther) )
            readCount['Unmapped'].append(starSummary['Unmapped'])
            readCount['Chimeric'].append(starSummary['Chimeric'])
            readCount['Unannotated'].append( '{0:.2f}'.format(100-annotated-starSummary['Unmapped']) )
            readPct['other'].append(str(readOther))
            tpmPct['other'].append(str(tpmOther))
            
    readData = []
    tpmData = []
    for cat in categories:
        trace1 = {'x': readPct[cat],
                  'y': samples,
                  'name': cat,
                  'orientation': 'h',
                  'type': 'bar',
                  'marker': {'width': 1}
              }
        trace2 = {'x': tpmPct[cat],
                  'y': samples,
                  'name': cat,
                  'orientation': 'h',
                  'type': 'bar',
                  'marker': {'width': 1}
              }
        trace3 = {'x': readCount[cat],
                  'y': samples,
                  'name': cat,
                  'orientation': 'h',
                  'type': 'bar',
                  'marker': {'width': 1}
              }
        readData.append(trace3)
        tpmData.append(trace2)

    for cat in ['Unannotated','Unmapped']:
        readData.append({'x': readCount[cat],
                         'y': samples,
                         'name': cat,
                         'orientation': 'h',
                         'type': 'bar',
                         'marker': {'width': 1}
                     })

    readCompJSON = json.dumps(readData)
    tpmCompJSON = json.dumps(tpmData)
    return readCompJSON,tpmCompJSON

def loadRunStats(sourceDir):
    runName = os.path.basename(sourceDir)
    if runName == '':
        runName = os.path.basename(os.path.dirname(sourceDir))

    statsHeader = []
    statsSamples = []
    with open(os.path.join(sourceDir,runName+'_summary.tsv'),'rU') as infile:
        header = next(infile)
        statsHeader = header.split()
        for line in infile:
            sample = line.split()
            statsSamples.append(sample)

    erccHeader = []
    erccSamples = []
    with open(os.path.join(sourceDir,runName+'_ercc.summary.tsv'),'rU') as infile:
        header = next(infile)
        erccHeader = header.split('\t')
        for line in infile:
            sample = line.split('\t')
            erccSamples.append(sample)

    genesHeader = ['Sample','Liver(>=1)','Liver(>=5)','Liver(>=10)','ProteinCoding(>=1)','ProteinCoding(>=5)','ProteinCoding(>=10)','Tissue(>=1)','Tissue(>=5)','Tissue(>=10)','Total']
    genesSamples = []
    with open(os.path.join(sourceDir,runName+'_genesDetected.summary.tsv'),'rU') as infile:
        next(infile)
        for line in infile:
            sample = re.sub(r'"','',line).split(',')
            genesSamples.append(sample)
    

    graphs = {'pdf':[], 'png':[]}
    for graph in glob.glob(os.path.join(sourceDir,runName+'*.pdf')):
        graphs['pdf'].append(os.path.basename(graph))
    for graph in glob.glob(os.path.join(sourceDir,runName+'*.png')):
        graphs['png'].append(os.path.basename(graph))
    tables = []
    for table in glob.glob(os.path.join(sourceDir,runName+'*.tsv')):
        tables.append(os.path.basename(table))
    for table in glob.glob(os.path.join(sourceDir,runName+'*.xlsx')):
        tables.append(os.path.basename(table))

    readCompJSON, tpmCompJSON = loadCompStats(sourceDir)

    contextDict = {'statsHeader': statsHeader,
                   'statsSamples': statsSamples,
                   'genesHeader': genesHeader,
                   'genesSamples': genesSamples,
                   'erccHeader': erccHeader,
                   'erccSamples': erccSamples,
                   'readCompJSON': readCompJSON,
                   'tpmCompJSON': tpmCompJSON,
                   'graphs': graphs,
                   'tables': tables,
               }

    return contextDict

def renderSamplePage(sourceDir,sampleName):
    #from django.conf import settings
    #settings.configure()
    #sys.path.append('/home/jzhuang@ms.local/run_reports/mysite')
    #os.environ['DJANGO_SETTINGS_MODULE'] = 'mysite.settings'
    #contextDict = loadStatsData('/mnt/shares2/analysis/161110_TH024_JRP050/','Shalfr2')
    #template = loader.get_template('templates/rnaseq_pipeline_report_template.html')
    #with open('rnaseq_pipeline_report_template.html','rU') as html_template:
        #template = Template(html_template.read())
    #print template.render(Context(contextDict))

    contextDict = loadStatsData(sourceDir,sampleName)
    templateLoader = jinja2.FileSystemLoader( searchpath="/" )
    templateEnv = jinja2.Environment( loader=templateLoader, autoescape=True )
    TEMPLATE_FILE = "/home/jzhuang@ms.local/run_reports/templates/rnaseq_pipeline_sample_report_template.html"
    template = templateEnv.get_template( TEMPLATE_FILE )
    return template.render( contextDict )

def renderRunPage(sourceDir):
    contextDict = loadRunStats(sourceDir)
    templateLoader = jinja2.FileSystemLoader( searchpath="/" )
    templateEnv = jinja2.Environment( loader=templateLoader, autoescape=True )
    TEMPLATE_FILE = "/home/jzhuang@ms.local/run_reports/templates/rnaseq_pipeline_run_report_template.html"
    template = templateEnv.get_template( TEMPLATE_FILE )
    return template.render( contextDict )
    

def main():
    #print renderRunPage('/home/jzhuang@ms.local/biodatashare/analysis/161209_LJI')
    print renderSamplePage('/home/jzhuang@ms.local/biodatashare/analysis/161129_TH021_VH030_new/','VH02')

if __name__=='__main__':
    main()
