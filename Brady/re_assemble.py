import os

gnm_id_txt = '/Users/songweizhi/Desktop/Brady/reassembly/mis-assembly.txt'
#gnm_id_txt = '/home-user/wzsong/Brady/mis-assembly.txt'

for each_gnm in open(gnm_id_txt):
    gnm_id         = each_gnm.strip()
    fastqc_cmd_raw = '/home-user/software/fastqc/v0.11.7/fastqc %s_trim1.fastq %s_trim2.fastq' % (gnm_id, gnm_id)
    trim_cmd       = 'java -jar /home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 %s_trim1.fastq %s_trim2.fastq %s_R1_P.fastq %s_R1_UP.fastq %s_R2_P.fastq %s_R2_UP.fastq ILLUMINACLIP:/home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 CROP:140 HEADCROP:10 SLIDINGWINDOW:5:25 MINLEN:50' % (gnm_id, gnm_id, gnm_id, gnm_id, gnm_id, gnm_id)
    fastqc_cmd     = '/home-user/software/fastqc/v0.11.7/fastqc %s_R1_P.fastq %s_R2_P.fastq' % (gnm_id, gnm_id)

    #spades_cmd = 'spades.py --only-assembler -1 reads/HKCCYLR1010_trim1.fastq -2 reads/HKCCYLR1010_trim2.fastq -o HKCCYLR1010_spade_wd -t 12'
    #spades_cmd = 'spades.py --isolate -1 reads/HKCCYLR1010_trim1.fastq -2 reads/HKCCYLR1010_trim2.fastq -o HKCCYLR1010_spade_wd -t 12'
    #spades_cmd = 'spades.py --isolate -k 55,77,99 -1 %s_R1_P.fastq -2 %s_R2_P.fastq -o %s_spade_wd_k55-99 -t 12' % (gnm_id, gnm_id, gnm_id)
    #spades_cmd = 'spades.py -k 55,77,99 -1 %s_R1_P.fastq -2 %s_R2_P.fastq -o %s_spade_wd_k55-99_no_isolate -t 12' % (gnm_id, gnm_id, gnm_id)
    #spades_cmd = 'spades.py --meta -k 55,77,99 -1 %s_R1_P.fastq -2 %s_R2_P.fastq -o %s_spade_wd_k55-99_meta -t 12' % (gnm_id, gnm_id, gnm_id)
    spades_cmd          = 'spades.py --isolate -k 33,55,77 -1 %s_R1_P.fastq -2 %s_R2_P.fastq -o %s_spades_k33-55-77_isolate -t 5' % (gnm_id, gnm_id, gnm_id)
    submit_spades_cmd   = 'submitHPC.sh --cmd "%s" -n 5 -c %s_spades_k33-55-77_isolate' % (spades_cmd, gnm_id)
    cp_gfa_cmd          =  'cp %s_spades_k33-55-77_isolate/assembly_graph_with_scaffolds.gfa gfa_files/%s_k33-55-77_isolate_assembly_graph_with_scaffolds.gfa' % (gnm_id, gnm_id)

    # print(trim_cmd)
    # os.system(trim_cmd)

    # print(fastqc_cmd)
    # os.system(fastqc_cmd)

    #print(submit_spades_cmd)
    print(cp_gfa_cmd)

'''


java -jar /home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 SZCCHNRI2007_L1_1.fq SZCCHNRI2007_L1_2.fq SZCCHNRI2007_R1_P.fq SZCCHNRI2007_R1_UP.fq SZCCHNRI2007_R2_P.fq SZCCHNRI2007_R2_UP.fq ILLUMINACLIP:/home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 CROP:140 HEADCROP:10 SLIDINGWINDOW:5:25 MINLEN:50
java -jar /home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 SZCCHNS30412_L1_1.fq SZCCHNS30412_L1_2.fq SZCCHNS30412_R1_P.fq SZCCHNS30412_R1_UP.fq SZCCHNS30412_R2_P.fq SZCCHNS30412_R2_UP.fq ILLUMINACLIP:/home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 CROP:140 HEADCROP:10 SLIDINGWINDOW:5:25 MINLEN:50
java -jar /home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 6 -phred33 SZCCHNS3014_L1_1.fq SZCCHNS3014_L1_2.fq SZCCHNS3014_R1_P.fq SZCCHNS3014_R1_UP.fq SZCCHNS3014_R2_P.fq SZCCHNS3014_R2_UP.fq ILLUMINACLIP:/home-user/software/trimmomatic/v0.36/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa:2:30:10 LEADING:30 TRAILING:30 CROP:140 HEADCROP:10 SLIDINGWINDOW:5:25 MINLEN:50

fastqc SZCCHNRI2007_R1_P.fq SZCCHNRI2007_R2_P.fq
fastqc SZCCHNS30412_R1_P.fq SZCCHNS30412_R2_P.fq
fastqc SZCCHNS3014_R1_P.fq SZCCHNS3014_R2_P.fq

submitHPC.sh --cmd "spades.py --isolate -k 33,55,77 -1 SZCCHNRI2007_R1_P.fq -2 SZCCHNRI2007_R2_P.fq -o SZCCHNRI2007_spades_k33-55-77_isolate -t 5" -n 5 -c SZCCHNRI2007
submitHPC.sh --cmd "spades.py --isolate -k 33,55,77 -1 SZCCHNS30412_R1_P.fq -2 SZCCHNS30412_R2_P.fq -o SZCCHNS30412_spades_k33-55-77_isolate -t 5" -n 5 -c SZCCHNS30412
submitHPC.sh --cmd "spades.py --isolate -k 33,55,77 -1 SZCCHNS3014_R1_P.fq -2 SZCCHNS3014_R2_P.fq -o SZCCHNS3014_spades_k33-55-77_isolate -t 5" -n 5 -c SZCCHNS3014


cd /home-user/wzsong/Brady/ITS/SZCCHNRI2007_spades_k33-55-77_isolate
BioSAK gfa2fa -gfa assembly_graph_with_scaffolds.gfa -o assembly_graph_with_scaffolds.gfa.fa

cd /home-user/wzsong/Brady/ITS/SZCCHNS3014_spades_k33-55-77_isolate
BioSAK gfa2fa -gfa assembly_graph_with_scaffolds.gfa -o assembly_graph_with_scaffolds.gfa.fa

cd /home-user/wzsong/Brady/ITS/SZCCHNS30412_spades_k33-55-77_isolate
BioSAK gfa2fa -gfa assembly_graph_with_scaffolds.gfa -o assembly_graph_with_scaffolds.gfa.fa



'''