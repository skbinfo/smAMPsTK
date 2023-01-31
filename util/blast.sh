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
blastp -query ../data/results/amp.fa -db ../blastp_alignment/amp/antimicrobial -evalue 0.00001 -out ../data/results/amp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query ../data/results/abp.fa -db ../blastp_alignment/abp/antimicrobial -evalue 0.00001 -out ../data/results/abp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query ../data/results/afp.fa -db ../blastp_alignment/afp/antimicrobial -evalue 0.00001 -out ../data/results/afp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query ../data/results/avp.fa -db ../blastp_alignment/avp/antimicrobial -evalue 0.00001 -out ../data/results/avp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70

##convert non_aligned fasta to tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output1 > ../data/results/amp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output2 > ../data/results/abp.tab 
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output3 > ../data/results/afp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output4 > ../data/results/avp.tab

##extract peptide sequence having BLAST positive result in fasta format
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ../data/results/amp_align ../data/results/amp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > ../data/results/amp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ../data/results/abp_align ../data/results/abp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > ../data/results/abp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ../data/results/afp_align ../data/results/afp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > ../data/results/afp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' ../data/results/avp_align ../data/results/avp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > ../data/results/avp_align.fa

rm ../data/results/amp.fa ../data/results/abp.fa ../data/results/afp.fa ../data/results/avp.fa
rm ../data/results/amp_align ../data/results/abp_align ../data/results/afp_align ../data/results/avp_align
rm ../data/results/amp.tab ../data/results/abp.tab ../data/results/afp.tab ../data/results/avp.tab