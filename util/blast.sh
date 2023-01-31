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
blastp -query amp.fa -db ../seq_alignment/amp/antimicrobial -evalue 0.00001 -out amp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query abp.fa -db ../seq_alignment/abp/antimicrobial -evalue 0.00001 -out abp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query afp.fa -db ../seq_alignment/afp/antimicrobial -evalue 0.00001 -out afp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70
blastp -query avp.fa -db ../seq_alignment/avp/antimicrobial -evalue 0.00001 -out avp_align -outfmt "6 qaccver qlen saccver pident length mismatch gapopen gaps qstart qend sstart send evalue bitscore qcovs" -num_threads 70

##convert non_aligned fasta to tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output1 > amp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output2 > abp.tab 
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output3 > afp.tab
awk 'BEGIN{RS=">";OFS="\t"}NR>1{print $1,$2}' $output4 > avp.tab

##extract peptide sequence having BLAST positive result in fasta format
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' amp_align amp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > amp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' abp_align abp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > abp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' afp_align afp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > afp_align.fa
awk 'NR==FNR{c[$1]++;next};c[$1] > 0' avp_align avp.tab|awk 'BEGIN{OFS='\t'}{print ">"$1"\n"$2; col1+=1}' > avp_align.fa