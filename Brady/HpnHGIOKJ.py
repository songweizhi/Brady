import os
import glob
import multiprocessing as mp


# ########################################################################################################################
#
# # file in
# ref_seq_fa      = '/home-user/wzsong/Brady/HpnHGIOKJ.fa'
# faa_file_dir    = '/home-user/wzsong/Brady/faa_files_394_PB'
# faa_file_ext    = 'protein'
# threads_num     = 12
# evalue_cutoff   = 0.0001
#
# # file out
# blast_op_dir            = '/home-user/wzsong/Brady/HpnHGIOKJ_vs_394_PB_blastp'
# blast_op_dir_best_hit   = '/home-user/wzsong/Brady/HpnHGIOKJ_vs_394_PB_blastp_best_hit'
#
# ########################################################################################################################
#
# faa_file_re = '%s/*.%s' % (faa_file_dir, faa_file_ext)
# faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
#
# blastp_cmd_list = []
# get_best_hit_cmd_list = []
# for each_faa in faa_file_list:
#     gnm_id                  = each_faa.split('.%s' % faa_file_ext)[0]
#     pwd_faa                 = '%s/%s' % (faa_file_dir, each_faa)
#     pwd_blast_op            = '%s/%s_blastp.txt' % (blast_op_dir, gnm_id)
#     pwd_blast_op_best_hit   = '%s/%s_blastp.txt' % (blast_op_dir_best_hit, gnm_id)
#
#     blastp_cmd              = 'blastp -query %s -subject %s -out %s -outfmt 6 -evalue %s' % (ref_seq_fa, pwd_faa, pwd_blast_op, evalue_cutoff)
#     get_best_hit_cmd        = 'BioSAK BestHit -i %s -o %s' % (pwd_blast_op, pwd_blast_op_best_hit)
#
#     blastp_cmd_list.append(blastp_cmd)
#     get_best_hit_cmd_list.append(get_best_hit_cmd)
#
# if os.path.isdir(blast_op_dir):
#     os.system('rm -r %s' % blast_op_dir)
# os.mkdir(blast_op_dir)
#
# # run blastp with multiprocessing
# pool = mp.Pool(processes=threads_num)
# pool.map(os.system, blastp_cmd_list)
# pool.close()
# pool.join()
#
# if os.path.isdir(blast_op_dir_best_hit):
#     os.system('rm -r %s' % blast_op_dir_best_hit)
# os.mkdir(blast_op_dir_best_hit)
#
# # get_best_hit with multiprocessing
# pool = mp.Pool(processes=threads_num)
# pool.map(os.system, get_best_hit_cmd_list)
# pool.close()
# pool.join()

from Bio import SeqIO

##########################################################################################

# file in
ref_seq_fa              = '/Users/songweizhi/Desktop/HpnHGIOKJ.fa'
blast_op_dir_best_hit   = '/Users/songweizhi/Desktop/HpnHGIOKJ_vs_394_PB_blastp_best_hit'
e_index_cutoff          = 30


# file out
op_df                   = '/Users/songweizhi/Desktop/HpnHGIOKJ_vs_394_PB_blastp_best_hit_df.txt'

##########################################################################################

best_hit_file_re = '%s/*.txt' % blast_op_dir_best_hit
best_hit_file_list = [os.path.basename(file_name) for file_name in glob.glob(best_hit_file_re)]


ref_gene_list = []
for each_gene in SeqIO.parse(ref_seq_fa, 'fasta'):
    ref_gene_list.append(each_gene.id)

ref_gene_list_sorted = sorted(ref_gene_list)


op_df_handle = open(op_df, 'w')
op_df_handle.write('Genome\t' + '\t'.join(ref_gene_list_sorted) + '\n')
for best_hit_file in best_hit_file_list:

    gnm_id                  = best_hit_file.split('_blastp.txt')[0]
    pwd_blast_op_best_hit   = '%s/%s'                               % (blast_op_dir_best_hit, best_hit_file)

    current_gnm_identified_gene = set()
    for each_hit in open(pwd_blast_op_best_hit):
        each_hit_split = each_hit.strip().split('\t')
        ref_id = each_hit_split[0]
        evalue = each_hit_split[10]

        if evalue == '0.0':
            current_gnm_identified_gene.add(ref_id)
        elif evalue == '0.0e+00':
            current_gnm_identified_gene.add(ref_id)
        else:
            if '-' in evalue:
                e_index = int(evalue.split('-')[-1])
                if e_index >= e_index_cutoff:
                    current_gnm_identified_gene.add(ref_id)
            else:
                print(each_hit_split)
                print(evalue)

    anno_op_list = [gnm_id]
    for each_ref in ref_gene_list_sorted:
        if each_ref in current_gnm_identified_gene:
            anno_op_list.append('1')
        else:
            anno_op_list.append('0')
    op_df_handle.write('\t'.join(anno_op_list) + '\n')
    if gnm_id in ['Bradyrhizobium_sp_SZCCHNS2021', 'Bradyrhizobium_sp_SZCCHNS2022', 'Bradyrhizobium_sp_SZCCHNRI3043']:
        print('\t'.join(anno_op_list))

op_df_handle.close()
