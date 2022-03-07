#!/bin/bash
# FOR PLANTS
# Author : qinz001@qq.com 

# Dependence : mmseq  signalp hmmsearch mcl mcxload and  Addann.py that like a join function

# Input Genome pep_less200.fa  and a hmm file download from MtSSPdb   ./SSPplant_HMM_v5Byqz.hmm

Infa=$1
hmmfile=$2
name=`basename ${Infa}`
# Step 0 

cat ${i} | paste - - | awk '{print ">"$1"\t"$2}' > ${name}.qzfa

./filterpep.py  ${name}.qzfa | sponge ${name}.qzfa

# Step 1 SignalP 

signalp -batch 1000000 -stdout -fasta ${Infa} > ./${name}.signalout

awk -F"\t" '{if($2!=OTHER && $3>=0.85 )print $1"\t"$NF ;else print $1"\t""Non-secretory protein"}' ./${name}.signalout \
    | grep -v "Non-secretory protein" > ./${name}.signalout.filter

Addann.py -i ./${name}.signalout.filter -c 1 -a ./${name}.qzfa |  grep -v "#" | sponge ./${name}.signalout.filter 


rm ./${name}.qzfa 

awk -F"\t" '{print ">"$1"\n"$3 }'  ./${name}.signalout.filter  > ./${name}.signalout.filter.fa

# Step 2 Hmmsearch 

hmmsearch --tblout  ./${name}.hmmout -o ./${name}.hmmoutall  ./${2} 

grep  -v "#" ./${name}.hmmout | awk '{print $1"\t"$3}' \
    | sort -k1,1 -u  > ${name}.family

# Step 3 MergeResult ${name}.family ${name}.signalout.filter 

Addann.py -i ./${name}.signalout.filter -c 1 -a ./${name}.family \
    |  grep -v "#" | sponge ./${name}.signalout.filter

cp ./${name}.signalout.filter  ./${name}.finalout

# Step 4 mmseq 

grep -w "-" ./${name}.finalout | awk -F"\t"  '{print ">"$1"\n"$3 }' > unknown.fa

mmseqs easy-cluster unknown.fa   clusterRes ./  --min-seq-id 0.5 -c 0.65 --cov-mode 1 

mcxload -abc ./clusterRes_cluster.tsv --stream-mirror -write-tab data.tab -o data.mci

mcl ./data.mci  -I 1.4  -use-tab data.tab

awk '{print "cluster"NR"\t"$0}' out.data.mci.I14 > ${name}.cluster

while read line ;do                                                                                                                  11:36PM
    gene=`echo ${line} | awk '{print $1}'  ` ;
    res=`grep -w ${gene} tmp1 | awk '{print $1}' ` 
    echo  $res $line
done <  unknown.fa >  unknown.cluster.fa


