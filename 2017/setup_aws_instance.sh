#!/bin/bash
cd /tmp
curl https://amazon-ssm-us-west-2.s3.amazonaws.com/latest/linux_amd64/amazon-ssm-agent.rpm -o amazon-ssm-agent.rpm
yum install -y amazon-ssm-agent.rpm

mkdir /mnt/efs
mount -t nfs4 -o nfsvers=4.1,rsize=1048576,wsize=1048576,hard,timeo=600,retrans=2 172.31.16.174:/ /mnt/efs
yum update -y
yum install -y docker
service docker start
usermod -a -G docker ec2-user
