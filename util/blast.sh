#!/bin/bash
input=$1
output1=$2
output2=$3
output3=$4
output4=$5
path=$6

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
blastp -query $output1 -db $path/blastp_alignment/amp/antimicrobial -evalue 0.00001 -out $path/data/results/amp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query $output2 -db $path/blastp_alignment/abp/antibacterial -evalue 0.00001 -out $path/data/results/abp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query $output3 -db $path/blastp_alignment/afp/antifungal -evalue 0.00001 -out $path/data/results/afp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query $output4 -db $path/blastp_alignment/avp/antiviral -evalue 0.00001 -out $path/data/results/avp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70

##convert non_aligned fasta to tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output1 > $path/data/results/amp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output2 > $path/data/results/abp.tab 
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output3 > $path/data/results/afp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output4 > $path/data/results/avp.tab

##extract peptide sequence having BLAST positive result in fasta format
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $path/data/results/amp_align $path/data/results/amp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > $path/data/results/amp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $path/data/results/abp_align $path/data/results/abp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > $path/data/results/abp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $path/data/results/afp_align $path/data/results/afp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > $path/data/results/afp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' $path/data/results/avp_align $path/data/results/avp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > $path/data/results/avp_align.fa

rm $path/data/results/amp.fa $path/data/results/abp.fa $path/data/results/afp.fa $path/data/results/avp.fa
rm $path/data/results/amp_align $path/data/results/abp_align $path/data/results/afp_align $path/data/results/avp_align
rm $path/data/results/amp.tab $path/data/results/abp.tab $path/data/results/afp.tab $path/data/results/avp.tab