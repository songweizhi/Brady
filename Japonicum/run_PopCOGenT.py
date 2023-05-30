import os

# uclust_16s_op_txt       = '/Users/songweizhi/Desktop/genome_16S_dereplicated_0.995.uc.reorganised.txt'
# gnm_to_cluster_txt      = '/Users/songweizhi/Desktop/gnm_to_cluster_0.995.txt'
# gnm_to_cluster_txt_itol = '/Users/songweizhi/Desktop/gnm_to_cluster_0.995_iTOL.txt'
#
# cluster_to_genome_dict = dict()
# genome_to_cluster_dict = dict()
# for each_c in open(uclust_16s_op_txt):
#     each_c_split = each_c.strip().split()
#     c_id = each_c_split[0]
#     sequence_list = each_c_split[1].split(',')
#     gnm_set = set()
#     for each_seq in sequence_list:
#         gnm_set.add(each_seq.split('_16S_')[0])
#     print('%s\t%s' % (c_id, len(gnm_set)))
#     cluster_to_genome_dict[c_id] = gnm_set
#
#     for each_gnm in gnm_set:
#         if each_gnm not in genome_to_cluster_dict:
#             genome_to_cluster_dict[each_gnm] = set()
#         genome_to_cluster_dict[each_gnm].add(c_id)
#
#
# gnm_to_cluster_txt_handle = open(gnm_to_cluster_txt, 'w')
# for each_gnm in genome_to_cluster_dict:
#     gnm_c_set = genome_to_cluster_dict[each_gnm]
#     if len(gnm_c_set) == 1:
#         gnm_to_cluster_txt_handle.write('%s\t%s\n' % (each_gnm, list(gnm_c_set)[0]))
#     else:
#         print('%s\t%s' % (each_gnm, gnm_c_set))
# gnm_to_cluster_txt_handle.close()


'''
cd /home-user/wzsong/Japonicum/PopCOGenT
BioSAK iTOL -ColorStrip -lg PopCOGenT_MCs.txt -lt PopCOGenT_MC -o PopCOGenT_MCs_iTOL_Strip.txt
BioSAK iTOL -ColorRange -lg PopCOGenT_MCs.txt -lt PopCOGenT_MC -o PopCOGenT_MCs_iTOL_Range.txt
'''
#
# PopC_txt    = '/Users/songweizhi/Desktop/PopC.txt'
# PopC_txt    = 'PopC.txt'
# gnm_dir     = '/home-user/wzsong/Japonicum/data/genome'
# gnm_ext     = 'genome'
# copy_to_dir = '/home-user/wzsong/Japonicum/PopCOGenT/genome_sep'
#
#
# if os.path.isdir(copy_to_dir) is True:
#     os.system('rm -r %s' % copy_to_dir)
# os.system('mkdir %s' % copy_to_dir)
#
# for each_gnm in open(PopC_txt):
#     each_gnm_split = each_gnm.strip().split('\t')
#     gnm_id = each_gnm_split[0]
#     cluster_id = each_gnm_split[1]
#
#     copy_to_subdir = '%s/%s' % (copy_to_dir, cluster_id)
#     if os.path.isdir(copy_to_subdir) is False:
#         os.system('mkdir %s' % copy_to_subdir)
#
#     pwd_gnm_from = '%s/%s.%s' % (gnm_dir, gnm_id, gnm_ext)
#     pwd_gnm_to = '%s/%s.%s' % (copy_to_subdir, gnm_id, gnm_ext)
#
#     os.system('cp %s %s' % (pwd_gnm_from, pwd_gnm_to))
#
#

