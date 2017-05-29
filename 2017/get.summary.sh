#! /bin/bash

for i in /mnt/shares/Users/jzhuang/*Chimeric.out.junction
do
    j=`basename $i`
    j=${j/junction/summary}
    cat $i | awk '$1!="chrM" && $4!="chrM" && $7>0 && $8+$9<=5 {print $1,$2,$3,$4,$5,$6,$7,$8,$9}' | sort | uniq -c | sort -k1,1rn > $j
done

