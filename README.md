# smAMPsTK
smAMPsTK is written in python to detect antimicrobial peptides from plants' transcriptome data. smAMPsTK detects peptides of four different activities i.e. antimicrobial, antibacterial, antifungal, and antiviral.
# Dependencies
1. python3

Python modules: orfipy (v0.0.4), orffinder (v1.8), pandas (v1.1.5), configparser (v5.2.0), biopython (v1.79) 

Python module can be easily installed by using following command:

`pip3 install < module name > --user`

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

<code>git clone https://github.com/skbinfo/smAMPsTK.git</code>

`cd smAMPsTK/`

Installing MiPepid and GPSR:

<code>unzip MiPepid-master.zip </code>

<code>unzip gpsr.zip </code>




