#! /usr/bin/env python3
'''
This program is intended to extract a number of antimicrobial peptides from the given 
set of transcriptomic data (in FASTA format). 

smAMPs.py: A pipeline for identification of smORF encoding Antimicrobal peptides
v1.0 19/10/2022
Author: Mohini Jaiswal
email: mohini12jaiswal@nipgr.ac.in
'''

import orfipy_core
import sys, getopt, os, subprocess
import argparse
import pandas as pd
import re
import stat
import configparser
from distutils import spawn
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import SimpleFastaParser
from orffinder import orffinder

def makeArgParser():
    parser = argparse.ArgumentParser(description='smAMPs.py: extract antimicrobial peptides from plant transcriptome data')
    parser.add_argument('-domain', help="Include domain search for the translated peptides, if found then represented as a column in the result output table", dest='dmsearch', action='store_true')
    parser.add_argument('configfile', help="having path of input file i.e. genomic file (in FASTA format)", action='store')
    parser.add_argument('-th', help="provide a threshold value (range lies between 0 to 1 OR 0 to -1) to classify result of SVM models.\nDefault: 0.0", dest= 'threshold', default='0.0', type=float, required=True)
    args = parser.parse_args()
    return args

def parse_fasta (lines):
    descs = []
    seqs = []
    data = ''
    for line in lines:
        if line.startswith('>'):
            if data:   # have collected a sequence, push to seqs
                seqs.append(data)
                data = ''
            descs.append(line[1:])  # Trim '>' from beginning
        else:
            data += line.rstrip('\r\n')
    # there will be yet one more to push when we run out
    seqs.append(data)
    return descs, seqs

def thres_result (filee, pos, neg):
    thres_val = []
    with open(filee, 'rb') as o:
        for line in o:
            if (float(line) >= float(args.threshold)):
                thres_val.append(pos)
            else:
                thres_val.append(neg)
    return thres_val

def parse_file (filee):
    svm_sc = []
    with open(filee, 'r') as o:
        for line in o:
            svm_sc.append(line.replace('\n', ''))
    return svm_sc

def orf_len_filt(infile, outfile):
    with open(outfile, 'w') as ouf:
        with open(infile) as handle:
            for record in SimpleFastaParser(handle):
                if ((len(str(record[1])) <=303) and (len(str(record[1])) % 3 == 0)):
                    ouf.write('%s\n' %('>' + record[0] + '\t' + record[1])) 
             
def pep_len_filt(infile, outfile):    
    with open(outfile, 'w') as ouf:
        with open(infile) as handle:
            for record in SimpleFastaParser(handle):
                if record[1].count('*') == 0:
                    if ((len(record[1]) >= 5) and (len(record[1]) <=101)):
                        ouf.write('%s\n' %('>' + record[0] + '\t' + record[1]))

def remove_d(infile, outfile):
    df = pd.read_csv(infile, header=None, sep='\t')
    df.columns=['header', 'Sequence']
    df.drop_duplicates(subset='Sequence', keep="first", inplace=True)
    with open(outfile, 'w') as ouf:
        for i, row in df.iterrows():
            ouf.write('%s\n' %(row['header'] + '\n' + row['Sequence']))

def dom_search (file1, file2, var1, var2):
    domid = []
    domdes = []
    mj1 = []
    did_n = []
    des_n = []
    sid = []
    silen = []
    with open(file1, 'r') as dl:
        for line in dl:
            domid.append(line.split('\t')[0])
            domdes.append(line.split('\t')[1].replace('\n', ''))
    with open(file2, 'r') as ipf:
        for line in ipf:
            for x, d in zip(domid, domdes):
                if line.split('\t')[11] == x:
                    mj1.append((line.split('\t')))
                    mj1.extend([d])
    for x, y in zip(var1, var2):
        for i in range(len(mj1)):
            if mj1[i][0] == x and int(mj1[i][2]) == int(y):
                did_n.append(mj1[i][11])
                des_n.append(mj1[i+1])
                sid.append(x)
                silen.append(y)
            else:
                pass

    return sid, silen, did_n, des_n

if __name__ == '__main__':
    args = makeArgParser()

    confile = args.configfile
    if not confile:
        parser.print_help()
        sys.exit(1)

    #to read config file
    config = configparser.ConfigParser()
    config.read(confile)
    in_file = config.get("input_req", "path1")
    mipe_tool = config.get("input_req", "Mipep_path2")
    interpro = config.get("input_req", "Interpro_path")
    out_path = config.get("output_folder", "final_result")

    ###
    MMPATH = sys.argv[0].replace('smAMPs.py','')
    if MMPATH == "":
        MMPATH = os.getcwd()
    else:
        pass
    MPATH=MMPATH.replace("//","/")

    os.chmod(MPATH+'/util/svm_classify', stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH|stat.S_IXOTH)
    os.chmod(MPATH+'/util/blast.sh', stat.S_IRWXU|stat.S_IRWXG|stat.S_IROTH|stat.S_IXOTH)

    ##to remove 'N' from input file, if any and rewrite it
    with open(out_path+'/input.fa', 'w') as ouf:
        for record in SeqIO.parse(in_file, 'fasta'):
            ouf.write('%s\n' % ('>' + record.description + '\n' + record.seq.replace('B', '').replace('D', '').replace('E', '').replace('F', '').replace('H', '').replace('I', '').replace('J', '').replace('K', '').replace('L', '').replace('M', '').replace('N', '').replace('O', '').replace('P', '').replace('Q', '').replace('R', '').replace('S', '').replace('U', 'T').replace('V', '').replace('W', '').replace('X', '').replace('Y', '').replace('Z', '')))

    ##code for MiPepid
    os.chdir(os.path.abspath(mipe_tool))
    args_three = ['python3', './src/mipepid.py', out_path+'/input.fa', out_path+'/out.txt']
    try:
        subprocess.run(args_three, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)

    ##code for orfipy
    descriptions, sequences = parse_fasta(open(out_path+'/input.fa', 'r').read().split('\n'))
    out_file = open(out_path+'/out.txt', 'a')
    out_file.write("Result of Orfipy\n")
    for faseq, head in zip(sequences, descriptions):
        for start,stop,strand,description in orfipy_core.orfs(faseq,minlen=15,maxlen=300):
            #print(head,start,stop,strand,description,faseq[start:stop])
            out_file.write("%s\n" % (head + ' ' + str(start) + ' ' + str(stop) + ' ' + strand + ' ' + description + ' ' + faseq[start:stop+3]))

    out_file.close()

    ##code for orffinder
    with open(out_path+'/out.txt', 'a') as o:
        o.write("Result of ORffinder\n")
        for rec in SeqIO.parse(out_path+'/input.fa', "fasta"):
            a= orffinder.getORFNucleotides(rec, return_loci=True, minimum_length=15, remove_nested=True)
            for orf in a:
                o.write("%s" %(rec.id + '\t' + 'the orfs is: ' + str(orf) + '\t' + orf["nucleotide"] + '\n'))
        o.write("%s" %('End of ORFfinder result'))   


    ###to extract seqID and smORFs from MiPepid_result 
    Id = []
    sequen = []
    with open(out_path+'/out.txt', 'r') as fp:
        with open (out_path+'/all_orf.txt', 'w') as jo:
            for line in iter(fp.readline, 'sORF_ID,sORF_seq,transcript_DNA_sequence_ID,start_at,end_at,classification,probability\n'):
                pass
            for line in iter(fp.readline, 'Result of Orfipy\n'):
                b = line.split(',')
                if (b[5] == 'coding'):
                    Id.append(b[2])
                    sequen.append(b[1])
                    jo.write("%s\n" % ('>' + Id[-1] + '\n' + sequen[-1]))
                
    ###to extract seqID and smORFs from orfipy_result 
    Id2 = []
    sequen2 = []
    with open(out_path+'/out.txt', 'r') as fp:
        with open (out_path+'/all_orf.txt', 'a') as joo:
            for line2 in iter(fp.readline, 'Result of Orfipy\n'):
                pass
            for line2 in iter(fp.readline, 'Result of ORffinder\n'):
                c = line2.split()
                Id2.append(c[0])
                sequen2.append(c[-1])
                joo.write("%s\n" % ('>' + Id2[-1] + '\n' + sequen2[-1]))

    ###to extract seqID and smORFs from ORFfinder_result
    Id3 = []
    sequen3 = []
    with open(out_path+'/out.txt', 'r') as fp:
        with open (out_path+'/all_orf.txt', 'a') as jooo:
            for line3 in iter(fp.readline, 'Result of ORffinder\n'):
                pass
            for line3 in iter(fp.readline, 'End of ORFfinder result'):
                d = line3.split('\t')
                Id3.append(d[0])
                sequen3.append(d[-1])
                jooo.write("%s" % ('>' + Id3[-1] + '\n' + sequen3[-1]))

    ####to remove redundant ORFs (with condition of length (303bp) and % 3 ==0)
    orf_len_filt(out_path+'/all_orf.txt', out_path+'/all_orf.tab')
    ####remove duplicates of tab separated file
    remove_d(out_path+'/all_orf.tab', out_path+'/all_forf.txt')
       
    ##translating ORFs by using EMBOSS transeq
    args_a = ['transeq', '-sequence',out_path+'/all_forf.txt', '-outseq',out_path+'/pep_seq.txt', '-frame','F', '-trim','Y']
    try:
        subprocess.run(args_a, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)

    ##prepare a peptide fasta file (no redundant sequence) of length range between 5-100aa
    pep_len_filt(out_path+'/pep_seq.txt', out_path+'/pep_seq.tab')
    remove_d(out_path+'/pep_seq.tab', out_path+'/pep_seq_f.txt') 
      
    os.chdir(os.path.abspath(MPATH))

    ##to prepare files for model
    cmd1 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/fasta2sfasta', '-i',out_path+'/pep_seq_f.txt', '-o',out_path+'/pep_f.sfa']
    cmd2 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/pro2aac', '-i',out_path+'/pep_f.sfa', '-o',out_path+'/pep_f.mono']
    cmd3 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/pro2dpc', '-i',out_path+'/pep_f.sfa', '-o',out_path+'/pep_f.dpc']
    cmd4 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/col2svm', '-i',out_path+'/pep_f.mono', '-o',out_path+'/pep_f_col.aac', '-s',str(0)]
    cmd5 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/col2svm', '-i',out_path+'/pep_f.dpc', '-o',out_path+'/pep_f_col.dpc', '-s',str(0)]

    commands = [cmd1, cmd2, cmd3, cmd4, cmd5]
    for i in commands:
        subprocess.check_call(i, stdin=None, stdout=None, stderr=None, shell=False)

    ##run svm_classify to get the results
    cmd_sv1 = [MPATH+'/util/svm_classify', out_path+'/pep_f_col.dpc', MPATH+'/models/model_amp', out_path+'/amp_res']
    cmd_sv2 = [MPATH+'/util/svm_classify', out_path+'/pep_f_col.aac', MPATH+'/models/model_abp', out_path+'/abp_res']
    cmd_sv3 = [MPATH+'/util/svm_classify', out_path+'/pep_f_col.aac', MPATH+'/models/model_afp', out_path+'/afp_res']
    cmd_sv4 = [MPATH+'/util/svm_classify', out_path+'/pep_f_col.aac', MPATH+'/models/model_avp', out_path+'/avp_res']
    command_sv = [cmd_sv1, cmd_sv2, cmd_sv3, cmd_sv4]
    for i in command_sv:
        subprocess.check_call(i, stdin=None, stdout=None, stderr=None, shell=False)

    amp_val = thres_result(out_path+'/amp_res', 'AMP', 'Non-AMP')
    abp_val = thres_result(out_path+'/abp_res', 'ABP', 'Non-ABP')
    afp_val = thres_result(out_path+'/afp_res', 'AFP', 'Non-AFP')
    avp_val = thres_result(out_path+'/avp_res', 'AVP', 'Non-AVP')

    ######
    seq_id, pept = parse_fasta (open(out_path+'/pep_seq_f.txt', 'r').read().split('\n'))
    length_sq = []
    for sq in pept:
        length_sq.append(len(sq))

    amp_sc = parse_file (out_path+'/amp_res')
    abp_sc = parse_file (out_path+'/abp_res')
    afp_sc = parse_file (out_path+'/afp_res')
    avp_sc = parse_file (out_path+'/avp_res')

    ##to write final output
    with open(out_path+'/output.tsv', 'w') as ouw:
        for idseq, sequ, length, am_score, am_res, ab_score, ab_res, af_score, af_res, av_score, av_res in zip(seq_id, pept, length_sq, amp_sc, amp_val, abp_sc, abp_val, afp_sc, afp_val, avp_sc, avp_val):
            ouw.write('%s\n' % (idseq + '\t' + sequ + '\t' + str(length) + '\t' + am_score + '\t' + am_res + '\t' + ab_score + '\t' + ab_res + '\t' + af_score + '\t' + af_res + '\t' + av_score + '\t' + av_res))

    ##for BLASTP alignment
    align_arg = ['sh', MPATH+'/util/blast.sh', out_path+'/output.tsv', out_path+'/amp.fa', out_path+'/abp.fa', out_path+'/afp.fa', out_path+'/avp.fa', MPATH]
    try:
        subprocess.run(align_arg, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)

    if args.dmsearch:
        os.chdir(os.path.abspath(interpro))
        ##for amp
        args_ipro = ['./interproscan.sh', '-i', out_path+'/amp_align.fa', '-f', 'TSV', '-d', out_path]
        try:
            subprocess.run(args_ipro, check=True)
            sys.exit
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)
        sid1, silen1, did_n1, des_n1 = dom_search (MPATH+'/util/domain_list', out_path+'/amp_align.fa.tsv', seq_id, length_sq)
        with open(out_path+'/domain_result_amp', 'w') as ow:
            for w, x, y, z in zip(sid1, silen1, did_n1, des_n1):
                ow.write(w + '\t' + str(x) + '\t' + y + '\t' + z + '\n')

        ##for abp
        args_ipro = ['./interproscan.sh', '-i', out_path+'/abp_align.fa', '-f', 'TSV', '-d', out_path]
        try:
            subprocess.run(args_ipro, check=True)
            sys.exit
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)
        sid1, silen1, did_n1, des_n1 = dom_search (MPATH+'/util/domain_list', out_path+'/abp_align.fa.tsv', seq_id, length_sq)
        with open(out_path+'/domain_result_abp', 'w') as ow:
            for w, x, y, z in zip(sid1, silen1, did_n1, des_n1):
                ow.write(w + '\t' + str(x) + '\t' + y + '\t' + z + '\n')

        ##for afp
        args_ipro = ['./interproscan.sh', '-i', out_path+'/afp_align.fa', '-f', 'TSV', '-d', out_path]
        try:
            subprocess.run(args_ipro, check=True)
            sys.exit
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)
        sid1, silen1, did_n1, des_n1 = dom_search (MPATH+'/util/domain_list', out_path+'/afp_align.fa.tsv', seq_id, length_sq)
        with open(out_path+'/domain_result_afp', 'w') as ow:
            for w, x, y, z in zip(sid1, silen1, did_n1, des_n1):
                ow.write(w + '\t' + str(x) + '\t' + y + '\t' + z + '\n')

        ##for avp
        args_ipro = ['./interproscan.sh', '-i', out_path+'/avp_align.fa', '-f', 'TSV', '-d', out_path]
        try:
            subprocess.run(args_ipro, check=True)
            sys.exit
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)
        sid1, silen1, did_n1, des_n1 = dom_search (MPATH+'/util/domain_list', out_path+'/avp_align.fa.tsv', seq_id, length_sq)
        with open(out_path+'/domain_result_avp', 'w') as ow:
            for w, x, y, z in zip(sid1, silen1, did_n1, des_n1):
                ow.write(w + '\t' + str(x) + '\t' + y + '\t' + z + '\n')
