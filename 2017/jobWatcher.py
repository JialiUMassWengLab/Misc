#! /usr/bin/python

import re
import os
import csv
import sys
import time
import json
import getpass
import argparse
import subprocess
import jobController

def main():
    parse = argparse.ArgumentParser(description='Monitoring jobs and automatically terminate instances where the job have finished',formatter_class=argparse.RawTextHelpFormatter)

    parse.add_argument('InfoFile',nargs=1,help='Info file name')
    args = parse.parse_args()
    print args

    instance_status = jobController.checkInstances()
    sampleList = []
    with open(args.InfoFile[0],'rU') as infile:
        reader = csv.DictReader(infile,delimiter='\t')
        for row in reader:
            if not row['instance_id'] in instance_status or not instance_status[row['instance_id']]=='running':
                print 'Instance %s for sample %s is not running!' % (row['instance_id'],row['sample_name'])
                continue
            sampleList.append(row)

    n_samples = len(sampleList)
    terminated = set([])
    while len(terminated) < n_samples:
        print sampleList
        terminateList = []
        keepList = []
        jobController.instanceCommand(sampleList,'docker ps -a')
        for sample in sampleList:
            for line in sample['command_output'].split('\n'):
                if not re.search(r'molestethdockerid',line):
                    continue
                a = line.split()
                if a[0]==sample['docker_id'] and re.search(r'Exited',line):
                    match = re.search(r'Exited \((\d+)\)',line)
                    if match.group(1)=='0':
                        terminateList.append(sample)
                        terminated.add(sample['sample_name'])
                else:
                    keepList.append(sample)

            sample.pop('command_output',None)
                                            
        print terminated
        print terminateList
        print keepList
        if len(terminateList) > 0:
            jobController.instanceCommand(terminateList,'docker logs')
            for sample in terminateList:
                sample.pop('command_output',None)
            jobController.copyOutputFromS3(terminateList)
            jobController.terminate(terminateList)

        jobController.instanceCommand(keepList,'docker logs')
        for sample in keepList:
            sample.pop('command_output',None)
        jobController.copyOutputFromS3(keepList)
        sampleList = keepList
        time.sleep(120)


if __name__=='__main__':
    main()

