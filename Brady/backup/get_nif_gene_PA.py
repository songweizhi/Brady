import os
import glob
from Bio import SeqIO

def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '..'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


################################################################################

# file in
blast_op_dir    = '/Users/songweizhi/Desktop/LingHK360_Zong1070_rpsblast_op_nif'
blast_op_ext    = 'blast8'
special_gnm_txt = '/Users/songweizhi/Desktop/special_gnm.txt'

# in file
genes_to_check_txt  = '/Users/songweizhi/Desktop/genes_to_check_nif.txt'
# or as string
genes_to_check_txt  = 'nifA,nifB,nifD,nifE,nifH,nifK,nifN,nifX'
gene_to_count_num   = 'nifH'

# file out
op_df               = '/Users/songweizhi/Desktop/nif_PA_df.txt'
nifH_copy_num_txt   = '/Users/songweizhi/Desktop/nifH_copy_num.txt'

################################################################################

# file in
blast_op_dir    = '/Users/songweizhi/Desktop/nif_seq_updated'
blast_op_ext    = 'blast8'
special_gnm_txt = ''

genes_to_check_txt  = 'nifA,nifB,nifD,nifE,nifH,nifK,nifN,nifX'
gene_to_count_num   = 'nifH'

# file out
op_df               = '/Users/songweizhi/Desktop/nif_seq_updated_nif_PA_df.txt'
nifH_copy_num_txt   = '/Users/songweizhi/Desktop/nif_seq_updated_nifH_copy_num.txt'

################################################################################

special_gnm_list = []
if os.path.isfile(special_gnm_txt) is True:
    for each_gnm in open(special_gnm_txt):
        special_gnm_list.append(each_gnm.strip())

################################################################################

# get gnm list
blast_op_re = '%s/*.%s' % (blast_op_dir, blast_op_ext)
blast_op_list = glob.glob(blast_op_re)

gene_list_to_check = []
if os.path.isfile(genes_to_check_txt) is True:
    for each_g in open(genes_to_check_txt):
        gene_list_to_check.append(each_g.strip())
else:
    gene_list_to_check = genes_to_check_txt.split(',')

op_df_handle = open(op_df, 'w')
op_df_handle.write('Genome\t' + '\t'.join(gene_list_to_check) + '\n')
nifH_copy_num_txt_handle = open(nifH_copy_num_txt, 'w')
for each_op in blast_op_list:
    file_path, file_basename, file_ext = sep_path_basename_ext(each_op)

    identified_nif_list = []
    if file_basename not in special_gnm_list:
        for each_hit in open(each_op):
            each_hit_split = each_hit.strip().split('\t')
            ref_id = each_hit_split[1]
            e_value = each_hit_split[-2]
            qualified_evalue = False
            if '-' in e_value:
                e_value_index = int(e_value.split('-')[-1])

                if ref_id == 'nifX':
                    if e_value_index >= 50:
                        qualified_evalue = True
                    else:
                        if e_value_index >= 100:
                            qualified_evalue = True
            elif float(e_value) == 0:
                qualified_evalue = True

            if qualified_evalue is True:
                identified_nif_list.append(ref_id)
    else:
        pass
        # for the 13 special strains
        # print(file_basename)
        # for each_hit in open(each_op):
        #     each_hit_split = each_hit.strip().split('\t')
        #     ref_id = each_hit_split[1]
        #     e_value = each_hit_split[-2]
        #     qualified_evalue = False
        #     if '-' in e_value:
        #         e_value_index = int(e_value.split('-')[-1])
        #
        #         if ref_id == 'nifX':
        #             if e_value_index >= 50:
        #                 qualified_evalue = True
        #             else:
        #                 if e_value_index >= 100:
        #                     qualified_evalue = True
        #     elif float(e_value) == 0:
        #         qualified_evalue = True
        #     if qualified_evalue is True:
        #         print(each_hit_split)
        # print()

    # write out nifH_copy_num
    nifH_copy_num = identified_nif_list.count(gene_to_count_num)
    nifH_copy_num_txt_handle.write('%s\t%s\n' % (file_basename, nifH_copy_num))

    # write out nif pa
    pa_list = []
    for nif_gene in gene_list_to_check:
        if nif_gene in identified_nif_list:
            pa_list.append('1')
        else:
            pa_list.append('0')
    op_df_handle.write('%s\t%s\n' % (file_basename, '\t'.join(pa_list)))
op_df_handle.close()
nifH_copy_num_txt_handle.close()
