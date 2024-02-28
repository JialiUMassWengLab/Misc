import os
import re
import sys
import glob
import tiledb
import tiledb.cloud
import tiledb.cloud.utilities
import tiledbvcf
import numpy as np

def get_credentials(profile, credentials_fn='~/.aws/credentials'):
    """
    Copied from Chris' code: https://github.com/oneTakeda/Computational-Biology-cdtools.git
    Computational-Biology-cdtools/cdtools/aws.py

    Get AWS access key ID and secret access key for a specified profile.
    """
    from os.path import expanduser
    home = expanduser("~")
    if credentials_fn[0] == '~':
        credentials_fn = expanduser(credentials_fn)
    with open(credentials_fn) as f:
        lines = f.read()
        profiles = [x.split('\n')[0][1:-1] for x in lines.strip().split('\n\n')]
        ids = dict(zip(profiles, [x.split('\n')[1].split(' = ')[1] for x in lines.strip().split('\n\n')]))
        keys = dict(zip(profiles, [x.split('\n')[2].split(' = ')[1] for x in lines.strip().split('\n\n')]))
    assert profile in profiles, 'Profile not found\n'
    access_key_id = ids[profile]
    secret_access_key = keys[profile]
    return(access_key_id, secret_access_key)


def ingestion(dataset,batchSize):
    key_id,secret_key = get_credentials('insight')
    config = tiledb.Config()
    config["vfs.s3.aws_access_key_id"] = key_id
    config["vfs.s3.aws_secret_access_key"] = secret_key
    config["vfs.s3.region"] = "us-east-1"
    s3_cfg = tiledbvcf.ReadConfig(tiledb_config=config)
    vfs = tiledb.VFS(config=config)

    array_uri = "s3://tak-insight-priv-tiledb-plat/groups/%s" % dataset
    if (vfs.is_dir(array_uri)):
        print(f"Deleting existing array '{array_uri}'")
        vfs.remove_dir(array_uri)
        print("Done.")
        
    all_samples = glob.glob(os.path.join(dataset,"*.vcf.gz"))
    print("%d samples in total to be ingested into dataset %s" % (len(all_samples),dataset))
    batches = []
    for i in range(0, len(all_samples), batchSize):
        batches.append(all_samples[i:i+batchSize])
    print("Divide into %d batches each of %d samples" % (len(batches),batchSize))

    ds = tiledbvcf.Dataset(uri=array_uri, mode="w", cfg=s3_cfg)
    ds.create_dataset(vcf_attrs=batches[0][0], enable_allele_count=False, enable_variant_stats=False)
    for i in range(0,len(batches)):
        print("Ingesting batch %d..." % (i+1))
        ds.ingest_samples(sample_uris = batches[i])
        
    
def main():
    print(
        f"tiledb v{tiledb.version.version}\n"
        f"tiledb-vcf v{tiledbvcf.version}\n"
        f"tiledb-cloud v{tiledb.cloud.version.version}\n"
    )
    dataset = sys.argv[1]
    batchSize = int(sys.argv[2])
    ingestion(dataset,batchSize)


if __name__=='__main__':
    main()
