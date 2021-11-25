#! /usr/bin/python

import os
import sys
import glob
import pandas as pd
#import pyreadstat
from sas7bdat import SAS7BDAT

def main():
    path = sys.argv[1]
    files = glob.glob(os.path.join(path,'*.sas7bdat'))
    for fn in files:
        if not os.path.exists(fn.replace('.sas7bdat','.csv')):
            print(fn)
            with SAS7BDAT(fn,encoding="GBK") as f:
                header = f.columns
                meta = {}
                for column in header:
                    #meta.update({column.name.decode("utf-8"): column.label.decode("utf-8")})
                    meta.update({column.name: column.label})
                pd.Series(meta).to_csv(fn.replace('.sas7bdat','.metadata'),header=None)
                
                df = f.to_data_frame()
                df.to_csv(fn.replace('.sas7bdat','.csv'),index=False)


if __name__=='__main__':
    main()
