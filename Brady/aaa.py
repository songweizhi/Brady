import os
import glob
import multiprocessing as mp


########################################################################################################################

# file in
ref_seq_fa      = '/home-user/wzsong/Brady/HpnHGIOKJ.fa'
faa_file_dir    = '/home-user/wzsong/Brady/faa_files_394_PB'
faa_file_ext    = 'protein'
threads_num     = 12
evalue_cutoff   = 0.0001

# file out
blast_op_dir            = '/home-user/wzsong/Brady/HpnHGIOKJ_vs_394_PB_blastp'
blast_op_dir_best_hit   = '/home-user/wzsong/Brady/HpnHGIOKJ_vs_394_PB_blastp_best_hit'

########################################################################################################################

faa_file_re = '%s/*.%s' % (faa_file_dir, faa_file_ext)
faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]

blastp_cmd_list = []
get_best_hit_cmd_list = []
for each_faa in faa_file_list:
    gnm_id                  = each_faa.split('.%s' % faa_file_ext)[0]
    pwd_faa                 = '%s/%s' % (faa_file_dir, each_faa)
    pwd_blast_op            = '%s/%s_blastp.txt' % (blast_op_dir, gnm_id)
    pwd_blast_op_best_hit   = '%s/%s_blastp.txt' % (blast_op_dir_best_hit, gnm_id)

    blastp_cmd              = 'blastp -query %s -subject %s -out %s -outfmt 6 -evalue %s' % (ref_seq_fa, pwd_faa, pwd_blast_op, evalue_cutoff)
    get_best_hit_cmd        = 'BioSAK BestHit -i %s -o %s' % (pwd_blast_op, pwd_blast_op_best_hit)

    blastp_cmd_list.append(blastp_cmd)
    get_best_hit_cmd_list.append(get_best_hit_cmd)

if os.path.isdir(blast_op_dir):
    os.system('rm -r %s' % blast_op_dir)
os.mkdir(blast_op_dir)

# run blastp with multiprocessing
pool = mp.Pool(processes=threads_num)
pool.map(os.system, blastp_cmd_list)
pool.close()
pool.join()

if os.path.isdir(blast_op_dir_best_hit):
    os.system('rm -r %s' % blast_op_dir_best_hit)
os.mkdir(blast_op_dir_best_hit)

# get_best_hit with multiprocessing
pool = mp.Pool(processes=threads_num)
pool.map(os.system, get_best_hit_cmd_list)
pool.close()
pool.join()


