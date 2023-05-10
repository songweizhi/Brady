import os
from Bio import SeqIO
import multiprocessing as mp


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


########################################################################################################################

# file in
gnm_id_txt      = '/Users/songweizhi/Desktop/get_shc_tree/gnm_id.txt'
ref_seq         = '/Users/songweizhi/Desktop/get_shc_tree/shc_ref.fa'
faa_dir         = '/Users/songweizhi/Desktop/get_shc_tree/protein'
faa_ext         = 'protein'
e_index_cutoff  = 30
num_threads     = 10

# file out
blast_op_dir    = '/Users/songweizhi/Desktop/get_shc_tree/blastp_op_dir'
faa_out         = '/Users/songweizhi/Desktop/get_shc_tree/shc.faa'
aln_out         = '/Users/songweizhi/Desktop/get_shc_tree/shc.aln'
aln_out_trimmed = '/Users/songweizhi/Desktop/get_shc_tree/shc_trimmed.aln'

########################################################################################################################

# prepare blastp cmds
blastp_cmd_list = []
for each_gnm in open(gnm_id_txt):
    gnm_id = each_gnm.strip()
    pwd_faa = '%s/%s.%s' % (faa_dir, gnm_id, faa_ext)
    blast_op = '%s/%s_blastp.txt' % (blast_op_dir, gnm_id)
    blastp_cmd = 'blastp -query %s -subject %s -out %s -outfmt 6 -evalue 0.0001' % (ref_seq, pwd_faa, blast_op)
    blastp_cmd_list.append(blastp_cmd)

# create blast_op_dir
if os.path.isdir(blast_op_dir):
    os.system('rm -r %s' % blast_op_dir)
os.mkdir(blast_op_dir)

# run blastp with multiprocessing
pool = mp.Pool(processes=num_threads)
pool.map(os.system, blastp_cmd_list)
pool.close()
pool.join()

# process blast results and extract sequences
faa_out_handle = open(faa_out, 'w')
for each_gnm in open(gnm_id_txt):

    gnm_id            = each_gnm.strip()
    pwd_faa           = '%s/%s.%s'                  % (faa_dir, gnm_id, faa_ext)
    blast_op          = '%s/%s_blastp.txt'          % (blast_op_dir, gnm_id)
    blast_op_best_hit = '%s/%s_blastp_best_hit.txt' % (blast_op_dir, gnm_id)

    # keep best hit
    best_hit(blast_op, blast_op_best_hit)

    for each_hit in open(blast_op_best_hit):
        each_hit_split = each_hit.strip().split('\t')
        ref_id         = each_hit_split[0]
        best_hit_id    = each_hit_split[1]
        evalue         = each_hit_split[10]

        # check evalue
        qualified_hit = False
        if evalue == '0.0':
            qualified_hit = True
        elif evalue == '0.0e+00':
            qualified_hit = True
        else:
            if '-' in evalue:
                e_index = int(evalue.split('-')[-1])
                if e_index >= e_index_cutoff:
                    qualified_hit = True
            else:
                print(each_hit_split)
                print(evalue)

        # write out sequence of qualified hits
        if qualified_hit is True:
            for each_seq in SeqIO.parse(pwd_faa, 'fasta'):
                if each_seq.id == best_hit_id:
                    faa_out_handle.write('>%s\n' % gnm_id)
                    faa_out_handle.write('%s\n' % str(each_seq.seq))
faa_out_handle.close()

# run mafft-einsi
mafft_cmd = 'mafft-einsi --thread %s --quiet %s  > %s' % (num_threads, faa_out, aln_out)
print(mafft_cmd)
os.system(mafft_cmd)

# trim msa
trimal_cmd = 'trimal -in %s -out %s -automated1' % (aln_out, aln_out_trimmed)
print(trimal_cmd)
os.system(trimal_cmd)

# run iqtree with untrimmed msa
iqtree_cmd = 'iqtree -s %s -m LG+G+F -pre %s -bb 1000 -nt %s -redo' % (aln_out, 'shc_tree_untrimmed', num_threads)
print(iqtree_cmd)
os.system(iqtree_cmd)

# run iqtree with trimmed msa
iqtree_cmd = 'iqtree -s %s -m LG+G+F -pre %s -bb 1000 -nt %s -redo' % (aln_out_trimmed, 'shc_tree_trimmed', num_threads)
print(iqtree_cmd)
os.system(iqtree_cmd)



Notes = '''

# Here are some notes, please ignore!

cd /home-user/wzsong/Brady/SHC_tree
submitHPC.sh --cmd "mafft-einsi --thread 9 --quiet shc.faa  > shc.aln" -n 9 -c mafft

cd /home-user/wzsong/Brady/SHC_tree
submitHPC.sh --cmd "trimal -in shc.aln -out shc_trimmed.aln -automated1" -n 1 -c trimal

cd /home-user/wzsong/Brady/SHC_tree
submitHPC.sh --cmd "iqtree -s shc_trimmed.aln -m LG+G+F -pre shc_tree -bb 1000 -nt 10 -redo" -n 10 -c model3_spp

# iqtree -s shc_trimmed.aln -m MFP -mset LG -mrate G -pre shc_tree -bb 1000 -nt 6 -redo
# iqtree -s shc_trimmed.aln -m MFP -mset LG,WAG -mrate G -pre shc_tree -bb 1000 -nt 6 -redo
# iqtree -s shc_trimmed.aln -m MFP+MERGE -mset LG -mrate G -pre shc_tree -bb 1000 -nt 6 -redo
-spp partitions.txt

conda activate ruby
RUBYLIB=$RUBYLIB:/home-user/wzsong/Software/ruby_lib_sswang/
export RUBYLIB
cd /home-user/wzsong/Brady/SHC_tree/shc_tree_LG
bash /home-user/sswang/LHW-tools/SEQ2TREE/SEQ2TREE.sh --seq_indir seq_dir --seq_suffix faa --outdir shc_tree_LG --iqtree --cpu 6 --force --LG

'''
