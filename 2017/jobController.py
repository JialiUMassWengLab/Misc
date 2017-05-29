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

def checkInstances():
    output = subprocess.check_output( 'aws ec2 describe-instances --output text --query "Reservations[*].Instances[*].[InstanceId,State.Name,Tags[0].Value]"', shell=True, stderr=subprocess.STDOUT )
    #print output

    instances_status = {}
    instances = output.split('\n')[:-1]
    for instance in instances:
        a = instance.split()
        instances_status.update({a[0]: a[1]})
        if len(a) > 2 and not re.search(r'\'s',a[2]):
            print '\t'.join(a)

    return instances_status


def printCommandOutput(sampleList):
    for sample in sampleList:        
        print '\t'.join(['Sample','Instance ID','Docker container ID'])
        print '\t'.join([sample['sample_name'],sample['instance_id'],sample['docker_id']])
        print '----------------------------------------------'
        print sample['command_output']
        print '=============================================='
        print '==============================================\n'


def copyOutputFromS3(sampleList):
    for sample in sampleList:
        try:
            subprocess.check_output( 'aws s3 cp s3://ms-sequencing-data/AWS\ Logs/%s/%s/awsrunShellScript/0.awsrunShellScript/stdout ./%s' % (sample['command_id'],sample['instance_id'],'_'.join([sample['sample_name'],sample['instance_id']])+'.out'),  shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print e.output

        try:
            subprocess.check_output( 'aws s3 cp s3://ms-sequencing-data/AWS\ Logs/%s/%s/awsrunShellScript/0.awsrunShellScript/stderr ./%s' % (sample['command_id'],sample['instance_id'],'_'.join([sample['sample_name'],sample['instance_id']])+'.err'),  shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print e.output

        try:
            subprocess.check_output( 'aws s3 rm --recursive s3://ms-sequencing-data/AWS\ Logs/%s/' % sample['command_id'],  shell=True, stderr=subprocess.STDOUT)
        except subprocess.CalledProcessError as e:
            print e.output


def submit(samples,InfoFile,configFile):
    if len(samples) == 0:
        print 'Sample list cannot be empty!'
        return

    config = {
        'imageID': 'ami-5ec1673e',
        'instanceType': 'c3.8xlarge',
        'dockerImage': 'molestethdockerid/rnaseq-pipeline-image1',
        'dockerCommand': 'docker run -d -v /mnt/efs:/mnt/efs -v /home/ec2-user:/host/home',
        'commandInContainer': 'bash /mnt/efs/scripts/RNASeq_pipeline_AWS_st.bash /mnt/efs/161110_TH024_JRP050 /mnt/efs/test/ ',
    }
    if not configFile=='' and os.path.isfile(configFile):
        with open(configFile,'rU') as infile:
            reader = csv.reader(infile,delimiter=':')
            for row in reader:
                config.update({row[0]:row[1]})
        print config
                
    jobsInfo = {}
    #launching instances and keep track of instance IDs
    for sample in samples:
        jobsInfo[sample] = {}
        instance_id = subprocess.check_output('aws ec2 run-instances --iam-instance-profile Name=SSMaccess-profile --image-id %s --count 1 --instance-type %s --key-name us_west_oregon --security-groups launch-wizard-3 --user-data file://setup_aws_instance.sh --block-device-mappings file://mapping.json --output text --query "Instances[0].InstanceId"' % (config['imageID'],config['instanceType']), shell=True, stderr=subprocess.STDOUT)
        instance_id = instance_id[:-1]
        name_tag = '_'.join([getpass.getuser(),InfoFile,sample])
        subprocess.check_output( 'aws ec2 create-tags --resources %s --tags Key=Name,Value=%s' % (instance_id,name_tag), shell=True, stderr=subprocess.STDOUT )
        if re.search(r'^i-\w+',instance_id):
            jobsInfo[sample].update({'instance_id':instance_id})
        else:
            print 'Launching failed for sample %s:\n%s' % (sample,instance_id)

    time.sleep(120)

    print jobsInfo
    #run docker
    for sample in samples:
        command_id = subprocess.check_output( 'aws ssm send-command --instance-ids %s --document-name "AWS-RunShellScript" --comment "Run docker command" --parameters commands=\'%s %s %s\' --output-s3-bucket-name ms-sequencing-data --output-s3-key-prefix "AWS Logs" --output text --query "Command.CommandId"' % (jobsInfo[sample]['instance_id'],config['dockerCommand'],config['dockerImage'],config['commandInContainer']+sample), shell=True, stderr=subprocess.STDOUT )[:-1]
        jobsInfo[sample].update({'command_id':command_id})

    print jobsInfo
    #get docker container IDs
    finishedCommands = 0
    timepassed = 0
    while finishedCommands < len(samples) and timepassed < 600:
        for sample in samples:
            if not 'docker_id' in jobsInfo[sample]:
                command_output = json.loads(subprocess.check_output( 'aws ssm list-command-invocations --instance-id %s --command-id %s --details' % (jobsInfo[sample]['instance_id'],jobsInfo[sample]['command_id']), shell=True, stderr=subprocess.STDOUT ))
                if len(command_output['CommandInvocations'])>0:
                    if len(command_output['CommandInvocations'][0]['CommandPlugins'])>0 and not command_output['CommandInvocations'][0]['CommandPlugins'][0]['Status']=='InProgress':
                        a = command_output['CommandInvocations'][0]['CommandPlugins'][0]['Output'].split('\n')
                        if re.search(r'[0-9a-f]+',a[0]):
                            jobsInfo[sample].update({'docker_id':a[0][:12]})
                            finishedCommands += 1
                        elif re.search(r'failed to run commands',a[2]):
                            command_id = subprocess.check_output( 'aws ssm send-command --instance-ids %s --document-name "AWS-RunShellScript" --comment "Run docker command" --parameters commands=\'%s %s %s\' --output text --query "Command.CommandId"' % (jobsInfo[sample]['instance_id'],config['dockerCommand'],config['dockerImage'],config['commandInContainer']+sample), shell=True, stderr=subprocess.STDOUT )[:-1]
                            jobsInfo[sample].update({'command_id':command_id})

        timepassed += 1
        time.sleep(1)

    #print jobs information to Info file
    print jobsInfo
    with open(InfoFile,'w') as outfile:
        outfile.write( '\t'.join(['sample_name','instance_id','command_id','docker_id'])+'\n' )
        for key,value in jobsInfo.iteritems():
            if 'docker_id' in value:
                outfile.write( '\t'.join([key,value['instance_id'],value['command_id'],value['docker_id']])+'\n' )
            else:
                print 'Submission of job %s (%s) failed!' % (key, value['instance_id'])


def instanceCommand(sampleList, cmd):
    for sample in sampleList:
        command = cmd
        if command == 'docker top' or command == 'docker logs' or command == 'docker stop' or command == 'docker start':
            command += ' %s' % sample['docker_id']

        if cmd == 'docker logs':
            command_id = subprocess.check_output( 'aws ssm send-command --instance-ids %s --document-name "AWS-RunShellScript" --comment "Run docker command" --parameters commands=\'%s\' --output-s3-bucket-name ms-sequencing-data --output-s3-key-prefix "AWS Logs" --output text --query "Command.CommandId"' % (sample['instance_id'],command), shell=True, stderr=subprocess.STDOUT )[:-1]
        else:
            command_id = subprocess.check_output( 'aws ssm send-command --instance-ids %s --document-name "AWS-RunShellScript" --comment "Run docker command" --parameters commands=\'%s\' --output text --query "Command.CommandId"' % (sample['instance_id'],command), shell=True, stderr=subprocess.STDOUT )[:-1]            
        sample.update({'command_id':command_id})

    finishedCommands = 0
    while finishedCommands < len(sampleList):
        for sample in sampleList:
            if not 'command_output' in sample:
                command_output = json.loads(subprocess.check_output( 'aws ssm list-command-invocations --instance-id %s --command-id %s --details' % (sample['instance_id'],sample['command_id']), shell=True, stderr=subprocess.STDOUT ))
                if len(command_output['CommandInvocations'])>0:
                    if len(command_output['CommandInvocations'][0]['CommandPlugins'])>0 and not command_output['CommandInvocations'][0]['CommandPlugins'][0]['Status']=='InProgress':
                        sample.update({'command_output': command_output['CommandInvocations'][0]['CommandPlugins'][0]['Output']})
                        sample.update({'status': command_output['CommandInvocations'][0]['CommandPlugins'][0]['Status']})
                        finishedCommands += 1    
        time.sleep(1)

    return sampleList


def terminate(sampleList):
    samples = instanceCommand(sampleList, 'docker stop')
    for sample in samples:
        if sample['status'] != 'Success':
            print 'Stop docker container failed for sample %s on instance %s!' % (sample['sample_name'],sample['instance_id'])
            del samples[samples.index(sample)]
    '''
    samples = instanceCommand(samples, 'sudo umount /mnt/efs')
    for sample in samples:
        if sample['status'] != 'Success':
            print 'Unmount EFS disk failed for sample %s!' % sample['sample_name']
            del samples[samples.index(sample)]
    '''
    for sample in samples:
        print subprocess.check_output( 'aws ec2 terminate-instances --instance-ids %s' % sample['instance_id'], shell=True, stderr=subprocess.STDOUT )


def watcher(InfoFile):
    instance_status = checkInstances()
    sampleList = []
    with open(InfoFile,'rU') as infile:
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
        instanceCommand(sampleList,'docker ps -a')
        for sample in sampleList:
            for line in sample['command_output'].split('\n'):
                if not re.search(r'molestethdockerid',line):
                    continue
                a = line.split()
                if a[0]==sample['docker_id'] and re.search(r'Exited',line):
                    match = re.search(r'Exited \((\d+)\)',line)
                    #if match.group(1)=='0':
                    if match:
                        terminateList.append(sample)
                        terminated.add(sample['sample_name'])
                else:
                    keepList.append(sample)

            sample.pop('command_output',None)

        print terminated
        print terminateList
        print keepList
        if len(terminateList) > 0:
            instanceCommand(terminateList,'docker logs')
            for sample in terminateList:
                sample.pop('command_output',None)
            copyOutputFromS3(terminateList)
            terminate(terminateList)

        instanceCommand(keepList,'docker logs')
        for sample in keepList:
            sample.pop('command_output',None)
        copyOutputFromS3(keepList)
        sampleList = keepList
        time.sleep(120)

    
def main():
    parse = argparse.ArgumentParser(description='''Predefined commands are:
    submit: launching instances and running docker images there
    top: running "docker top" to check the specified containers
    ps: listing all the running docker containers
    ls: listing all the files in the /home/ec2-user directory
    log: chekcing outputs/errors on the containers
    terminate: terminating specified instances''',formatter_class=argparse.RawTextHelpFormatter)

    parse.add_argument('Command',nargs=1,help='Commands to be executed on the instance(s)')
    parse.add_argument('InfoFile',nargs=1,help='Info file name')
    #parse.add_argument('--command-to-run','-c',nargs=1,dest='command',help='Commands to be executed on the instance(s)',default='top')
    parse.add_argument('--sample-list','-s',nargs='?',dest='sample_list',help='The list of samples (separated by ",") to be checked. All samples are checked if omitted')
    parse.add_argument('--instance-list','-i',nargs='?',dest='instance_list',help='The list of instances (separated by ",") to be checked. All instances are checked if omitted')
    parse.add_argument('--configure-file','-c',nargs='?',default='',dest='configFile',help='Configure file for submitting jobs')
    args = parse.parse_args()
    print args

    if args.Command[0] == 'submit':
        submit(args.sample_list.split(','),args.InfoFile[0],args.configFile)
        print '\n-----------Submissions completed------------\n'
        watcher(args.InfoFile[0])
    else:
        cmdDict = {
            'top': 'docker top',
            'ps':  'docker ps -a',
            'ls':  'ls -l /home/ec2-user',
            'log': 'docker logs',
        }
        instance_status = checkInstances()

        #l = ['VH01','VH02']
        sampleList = []
        with open(args.InfoFile[0],'rU') as infile:
            reader = csv.DictReader(infile,delimiter='\t')
            for row in reader:
                if not row['instance_id'] in instance_status or not instance_status[row['instance_id']]=='running':
                    print 'Instance %s for sample %s is not running!' % (row['instance_id'],row['sample_name'])
                    continue

                if args.sample_list == None and args.instance_list == None:
                    sampleList.append(row)
                elif not args.sample_list == None and row['sample_name'] in args.sample_list.split(','):
                    sampleList.append(row)
                elif not args.instance_list == None and row['instance_id'] in args.instance_list.split(','):
                    sampleList.append(row)

            if args.Command[0] == 'terminate':
                terminate(sampleList)
            else:
                command = args.Command[0]
                if command in cmdDict:
                    command = cmdDict[command]

                samples = instanceCommand(sampleList,command)
                printCommandOutput(samples)
                if args.Command[0] == 'log':
                    copyOutputFromS3(samples)



if __name__=='__main__':
    main()

