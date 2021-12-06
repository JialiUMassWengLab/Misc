#! /usr/bin/bash

while read -r line;
do
    echo $line;
    synapse -u jzhuang_denovo -p Welc0me@denovo get $line --downloadLocation Source;
done < Source/wgs_list.txt
