#!/bin/bash
input=$1
output1=$2
output2=$3
output3=$4
output4=$5

touch $output1 $output2 $output3 $output4

while read -r line; do
set $line
x='>'$1'\t'$2
if  [ $5 = 'AMP' ]; then 
printf '%b\n' "$x" >> $output1
fi
if  [ $7 = 'ABP' ]; then 
printf '%b\n' "$x" >> $output2
fi
if  [ $9 = 'AFP' ]; then 
printf '%b\n' "$x" >> $output3
fi
if  [ $11 = 'AVP' ]; then 
printf '%b\n' "$x" >> $output4
fi

done <$input

sed -i 's/\t/\n/g' $output1 $output2 $output3 $output4

##run blastp
blastp -query amp.fa -db ../seq_alignment/antimicrobial -evalue 0.00001 -out amp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70