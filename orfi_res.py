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
import re
import configparser
from distutils import spawn
from Bio import SeqIO
from Bio.Seq import Seq
from orffinder import orffinder
from subprocess import Popen, PIPE

def makeArgParser():
    parser = argparse.ArgumentParser(description='smAMP.py: extract antimicrobial peptides from plant genomic data')
    parser.add_argument('-out_dir', help="Output Files Directory; Default: current diretory")
    parser.add_argument('-domain', help="Include domain search for the translated peptides, if found then represented as a column in the result output table", dest='dmsearch', action='store_true')
    parser.add_argument('configfile', help="having path of input file i.e. genomic file (in FASTA format)", action='store')
    parser.add_argument('-th', help="provide a threshold value (range lies between 0 to 1 OR 0 to -1) to classify result of SVM models.\nDefault: 0.0", dest= 'threshold', default='0.0', type=float, required=True)
    parser.add_argument('-t', help="number of threads [default: 1]", dest='threads', default='1', type=int)
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

    #create directory for storing SignalP result
    if not os.path.exists(out_path + '/signal_res'):
        os.makedirs(out_path + '/signal_res')

    #to create output directory if not provided
    outdir = args.out_dir
    if not outdir:
        outdir="smAMP_"+os.path.basename(in_file)+'_out'

    ##to remove 'N' from input file, if any and rewrite it
    with open(out_path+'/input.fa', 'w') as ouf:
        for record in SeqIO.parse(in_file, 'fasta'):
            ouf.write('%s\n' % ('>' + record.description + '\n' + record.seq.replace('B', '').replace('D', '').replace('E', '').replace('F', '').replace('H', '').replace('I', '').replace('J', '').replace('K', '').replace('L', '').replace('M', '').replace('N', '').replace('O', '').replace('P', '').replace('Q', '').replace('R', '').replace('S', '').replace('U', '').replace('V', '').replace('W', '').replace('X', '').replace('Y', '').replace('Z', '')))

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
                    jo.write("%s\n" % (Id[-1] + '\t' + sequen[-1]))
                
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
                joo.write("%s\n" % (Id2[-1] + '\t' + sequen2[-1]))

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
                jooo.write("%s" % (Id3[-1] + '\t' + sequen3[-1]))

    ##remove duplicate smORF from all ORFs file
    with open(out_path+'/all_orf.txt', "r") as txt_file:
        with open(out_path+'/uniq_orf.txt', 'w') as uni_seq:
            new_data = list(set(txt_file))
            for uniq in new_data:
                uni_seq.write("%s" %('>' + uniq.replace('\t', '\n')))
    '''
    ##translating ORFs by using EMBOSS transeq
    args_a = ['transeq', '-sequence',out_path+'/uniq_orf.txt', '-outseq',out_path+'/pep_seq.txt', '-frame','F', '-trim','Y']
    try:
        subprocess.run(args_a, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)

    ##removing peptides having * in between the sequences
    handle = open(out_path+'/pep_seq.txt', "r")
    filtered = [record for record in SeqIO.parse(handle, "fasta") if record.seq.count('*') == 0]
    output_handle = open(out_path+'/pep_seq_nostop.txt', "w")
    SeqIO.write(filtered, output_handle, "fasta")
    output_handle.close()
    handle.close()

    ##multiline fasta to single line fasta
    with open(out_path+'/pep_seq_nostop.txt', 'r') as f_input, open(out_path+'/pep_seq_f.txt', 'w') as f_output:
        block = []

        for line in f_input:
            if line.startswith('>'):
                if block:
                    f_output.write(''.join(block) + '\n')
                    block = []
                f_output.write(line)
            else:
                block.append(line.strip())

        if block:
            f_output.write(''.join(block) + '\n')
    
    ##prepare a file of fasta sequences of length ranges between 5-255aa
    with open(out_path+'/pep_seq_f.tab', 'w') as inf:
        for record in SeqIO.parse(out_path+'/pep_seq_f.txt', 'fasta'):
            inf.write('>{}\t{}'.format(record.description, record.seq + '\n'))

    with open(out_path+'/pep_seq_f.tab', 'r') as inf, open(out_path+'/pep_seq_f2.txt', 'w') as ouf:
        new_d = list(set(inf))
        for line in new_d:
            if ((len(line.split('\t')[1])>= 6) and (len(line.split('\t')[1])<= 256)) :
                ouf.write("%s\n" % (line.strip().replace('\t', '\n')))

    
    
    ##remove signal peptide region from sequences if any
    args_sig = ['signalp6', '-ff',out_path+'/pep_seq_f.txt', '-od',out_path + '/signal_res']
    try:
        subprocess.run(args_sig, check=True)
        sys.exit
    except subprocess.CalledProcessError as e:
        sys.exit(e.stderr)

    data2=[]
    lineco=[]
    lines2 = []
    tabseq=[]
    with open(out_path + '/signal_res/prediction_results.txt', 'r') as fin:
        data = fin.read().splitlines(True)
        data2.append(data[2:])
    for num, line in enumerate(data2[0], 1):
        parts=line.split("\t")
        if parts[1] == 'SP':
            lineco.append(num)
    for record in SeqIO.parse(out_path+'/pep_seq_f.txt', 'fasta'):
        tabseq.append('>{}\t{}'.format(record.description, record.seq))
    for i, line in enumerate(tabseq,1):
        if i not in lineco:
            lines2.append(line.strip())

    with open(out_path + '/signal_res/processed_entries.fasta', 'a') as fp:
        for item in lines2:
            fp.write("%s\n" %(item.replace('\t', '\n')))
    

    MMPATH = sys.argv[0].replace('orfi_res.py','')
    if MMPATH == "":
        MMPATH = os.getcwd()
    else:
        pass
    MPATH=MMPATH.replace("//","/")
    print(MPATH)

    ##to prepare files for model
    #cmd1 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/fasta2sfasta', '-i',out_path+'/signal_res/processed_entries.fasta', '-o',out_path+'/pep_f.sfa']
    cmd1 = ['perl', MPATH+'/util/gpsr_1.0/gpsr/bin/fasta2sfasta', '-i',out_path+'/pep_seq_f2.txt', '-o',out_path+'/pep_f.sfa']
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
    #seq_id, pept = parse_fasta (open(out_path+'/signal_res/processed_entries.fasta', 'r').read().split('\n'))
    seq_id, pept = parse_fasta (open(out_path+'/pep_seq_f2.txt', 'r').read().split('\n'))
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


    if args.dmsearch:
        os.chdir(os.path.abspath(interpro))
        args_ipro = ['./interproscan.sh', '-i', out_path+'/signal_res/processed_entries.fasta', '-f', 'TSV', '-d', out_path]
        try:
            subprocess.run(args_ipro, check=True)
            sys.exit
        except subprocess.CalledProcessError as e:
            sys.exit(e.stderr)
        sid1, silen1, did_n1, des_n1 = dom_search (MPATH+'/util/domain_list', out_path+'/processed_entries.fasta.tsv', seq_id, length_sq)
        with open(out_path+'/domain_result', 'w') as ow:
            for w, x, y, z in zip(sid1, silen1, did_n1, des_n1):
                ow.write(w + '\t' + str(x) + '\t' + y + '\t' + z + '\n')
    
    '''