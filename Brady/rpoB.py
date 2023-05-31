import os
import glob


def sep_path_basename_ext(file_in):
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'
    f_base, f_ext = os.path.splitext(file_name)
    return f_path, f_base, f_ext


query_file  = '/home-user/wzsong/Brady/rpoB/G2026_rpoB2106_2.fas'
gnm_txt     = '/Users/songweizhi/Desktop/genome.txt'
cmd_txt     = '/Users/songweizhi/Desktop/blastn_cmds_2.txt'


query_path, query_base, query_ext = sep_path_basename_ext(query_file)

cmd_txt_handle = open(cmd_txt, 'w')
for each_gnm in open(gnm_txt):
    pwd_gnm = each_gnm.strip()
    gnm_path, gnm_base, gnm_ext = sep_path_basename_ext(pwd_gnm)
    blast_cmd = 'blastn -query %s -subject %s -out %s_vs_%s.txt -outfmt 6' % (query_file, pwd_gnm, query_base, gnm_base)
    cmd_txt_handle.write(blast_cmd + '\n')
cmd_txt_handle.close()

########################################################################################################################

blast_op_dir_1 = '/Users/songweizhi/Desktop/blast_op_G2026_rpoB2106_1'
blast_op_dir_2 = '/Users/songweizhi/Desktop/blast_op_G2026_rpoB2106_2'
file_ext = 'txt'
blast_op_dir = blast_op_dir_2

file_re = '%s/*.%s' % (blast_op_dir, file_ext)
file_list = glob.glob(file_re)
op_txt = blast_op_dir + '_bitscore.txt'

op_txt_handle = open(op_txt, 'w')
for each_file in file_list:
    f_path, f_base, f_ext = sep_path_basename_ext(each_file)
    gnm_id = f_base.split('_vs_')[1]
    best_hit_line = open(each_file).readline()
    best_hit_line_split = best_hit_line.strip().split()
    if len(best_hit_line_split) == 0:
        bit_score = 0
    else:
        bit_score = int(best_hit_line_split[-1])
    op_txt_handle.write('%s\t%s\n' % (gnm_id, bit_score))
op_txt_handle.close()
