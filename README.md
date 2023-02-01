# smAMPsTK
__smAMPsTK__ is written in python to detect antimicrobial peptides from plants' transcriptome data. smAMPsTK detects peptides of four different activities i.e. antimicrobial, antibacterial, antifungal, and antiviral.
# Dependencies
1. python3

Python modules: orfipy (v0.0.4), orffinder (v1.8), pandas (v1.1.5), configparser (v5.2.0), biopython (v1.79) 

Python module can be easily installed by using following command:
```
pip3 install < module name > --user
```
2. EMBOSS (v6.6.0) -transeq

3. ncbi-blast (v2.9.0+) 

[Note: Python3, EMBOSS, and ncbi-blast are needed to be globally installed]

4. MiPepid tool 

5. GPSR package 

6. SVM classify

[Note: MiPepid and GPSR package are provided as zip file. SVM classify is provided as a script]

7. interproscan (v5.56-89.0)

# Installation
Download
```
git clone https://github.com/skbinfo/smAMPsTK.git

cd smAMPsTK/
```
Installing MiPepid and GPSR:
```
unzip MiPepid-master.zip 

unzip gpsr.zip 
```
[Note: GPSR will be used from util directory]

# Usage
First change the path of input file, MiPepid, InterProScan, and output folder in the provided `config.ini` file

Once config file is set according to the user's path;

For AMPs prediction from smORFs:
```
python3 smAMPs.py config.ini -th 0 
```
For domain search `-domain` argument is needed to be provided, like this:
```
python3 smAMPs.py config.ini -th 0 -domain 
```
## Regarding the datasets

<a href="https://doi.org/10.5281/zenodo.7592940"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.7592940.svg" alt="DOI"></a>

`analysis_data.tar.gz` file on aforementioned zenodo link contains generated AMPs data of five plant organisms used in this paper. 

## 

If you have any questions, bug reports, or suggestions, please e-mail

Dr. Shailesh Kumar

Staff Scientist, Bioinformatics Laboratory #202

National Institute of Plant Genome Research (NIPGR), New Delhi

shailesh@nipgr.ac.in
