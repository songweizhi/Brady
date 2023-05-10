# import os
#
# file_1 = 0
# file_0 = 0
# file_not_found = []
# for each_gnm in open('/home-user/wzsong/Brady/PB_gnm_list.txt'):
#     gnm_id = each_gnm.strip()
#     pwd_faa = '/home-user/sswang/project/Brady/data/zong/protein/%s.protein' % gnm_id
#     pwd_faa = '/home-user/sswang/project/Brady/data/LingHK/old/protein/%s.protein' % gnm_id
#     pwd_faa = '/home-user/wzsong/Brady/faa/%s.protein' % gnm_id
#     if os.path.isfile(pwd_faa) is False:
#         file_not_found.append(gnm_id)
#         file_0 += 1
#     else:
#         file_1 += 1
#
# print('\n'.join(file_not_found))
# print(len(file_not_found))
#
# print('not found %s' % file_0)
# print('found %s' % file_1)





# import os
# import glob
# import multiprocessing as mp
#
#
# ########################################################################################################################
#
# # file in
# #faa_file_dir    = '/home-user/sswang/project/Brady/data/zong/protein'
# faa_file_dir    = '/home-user/wzsong/Brady/faa'
# faa_file_ext    = 'protein'
# threads_num     = 12
# evalue_cutoff   = 0.001
# kegg_db_diamond = '/home-user/wzsong/resource/db/kegg/v2017/blastdb/kegg_prot_bacteria_201704.dmnd'
#
# # file out
# blast_op_dir    = 'blast_op_dir_diamond'
#
# ########################################################################################################################
#
# faa_file_re = '%s/*.%s' % (faa_file_dir, faa_file_ext)
# faa_file_list = [os.path.basename(file_name) for file_name in glob.glob(faa_file_re)]
#
# if os.path.isdir(blast_op_dir):
#     os.system('rm -r %s' % blast_op_dir)
# os.mkdir(blast_op_dir)
#
#
# diamond_cmd_list = []
# for each_faa in faa_file_list:
#     gnm_id       = each_faa.split('.%s' % faa_file_ext)[0]
#     pwd_faa      = '%s/%s' % (faa_file_dir, each_faa)
#     pwd_blast_op = '%s/%s_.txt' % (blast_op_dir, gnm_id)
#     diamond_cmd = 'diamond blastp -q %s --db %s --out %s --outfmt 6 --evalue %s --block-size 1 --threads 1 --quiet' % (pwd_faa, kegg_db_diamond, pwd_blast_op, evalue_cutoff)
#     #print(diamond_cmd)
#     #os.system(diamond_cmd)
#     diamond_cmd_list.append(diamond_cmd)
#
# # run blastp with multiprocessing
# pool = mp.Pool(processes=threads_num)
# pool.map(os.system, diamond_cmd_list)
# pool.close()
# pool.join()


import os

blast_op_dir            = '/Users/songweizhi/Desktop/lipid/blast_op_dir_diamond'
blast_op_dir_best_hit   = '/Users/songweizhi/Desktop/lipid/blast_op_dir_diamond_394_PB_best_hit'
gene2ko_txt             = '/Users/songweizhi/DB/KEGG_v2017_luo/kegg_bacteria.gene2ko'
interested_ko_txt       = '/Users/songweizhi/Desktop/lipid/interested_ko.txt'
op_txt                  = '/Users/songweizhi/Desktop/lipid/interested_ko_df.txt'

# read in interested ko
interested_ko_to_des_dict = dict()
for each_ko in open(interested_ko_txt):
    each_ko_split = each_ko.strip().split('\t')
    interested_ko_to_des_dict[each_ko_split[0]] = each_ko_split[1]

# read in kegg db
gene2ko_dict = dict()
for each_gene in open(gene2ko_txt):
    each_gene_split = each_gene.strip().split('\t')
    gene_id = each_gene_split[0]
    ko_id = each_gene_split[1].split(':')[1]
    gene2ko_dict[gene_id] = ko_id
# print(gene2ko_dict)
# 'zpr:ZPR_1225': 'K22757', 'zpr:ZPR_3938': 'K00265'


gnm_annotation_dod = dict()
for each_gnm in open('/Users/songweizhi/Desktop/lipid/PB_gnm_list.txt'):

    gnm_id                 = each_gnm.strip()
    pwd_blastn_op          = '%s/%s_.txt'                 % (blast_op_dir, gnm_id)
    pwd_blastn_op_best_hit = '%s/%s.txt'                  % (blast_op_dir_best_hit, gnm_id)
    get_best_hit_cmd       = 'BioSAK BestHit -i %s -o %s' % (pwd_blastn_op, pwd_blastn_op_best_hit)
    #print(gnm_id)

    # get best hit
    # os.system(get_best_hit_cmd)

    if gnm_id in ['Bradyrhizobium_sp_SZCCHNS2021', 'Bradyrhizobium_sp_SZCCHNS2022', 'Bradyrhizobium_sp_SZCCHNRI3043']:
        for each_line in open(pwd_blastn_op_best_hit):
            each_line_split = each_line.strip().split('\t')
            ref_id = each_line_split[1]
            ref_ko = gene2ko_dict.get(ref_id, '')
            if ref_ko == 'K06045':
                print('%s\t%s' % (gnm_id, each_line.strip()))

    current_gnm_annotation_dict = dict()
    for each_line in open(pwd_blastn_op_best_hit):
        each_line_split = each_line.strip().split('\t')
        ref_id = each_line_split[1]
        ref_ko = gene2ko_dict.get(ref_id, '')
        if ref_ko in interested_ko_to_des_dict:
            if ref_ko not in current_gnm_annotation_dict:
                current_gnm_annotation_dict[ref_ko] = 1
            else:
                current_gnm_annotation_dict[ref_ko] += 1

    gnm_annotation_dod[gnm_id] = current_gnm_annotation_dict


interested_ko_list_sorted = sorted([i for i in interested_ko_to_des_dict])

interested_ko_list_sorted_with_des = []
for each in interested_ko_list_sorted:
    interested_ko_list_sorted_with_des.append('%s_%s' % (each, interested_ko_to_des_dict[each]))

op_txt_handle = open(op_txt, 'w')
op_txt_handle.write('Genome\t%s\n' % '\t'.join(interested_ko_list_sorted_with_des))
for each_gnm in gnm_annotation_dod:
    anno_dict = gnm_annotation_dod[each_gnm]
    num_list = [each_gnm]
    for each_ko in interested_ko_list_sorted:
        #hits_list = anno_dict.get(each_ko, [])
        #ko_num = len(hits_list)
        ko_num = anno_dict.get(each_ko, '0')
        num_list.append(str(ko_num))
    op_txt_handle.write('%s\n' % ('\t'.join(num_list)))
op_txt_handle.close()

