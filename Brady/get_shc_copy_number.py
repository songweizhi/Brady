import os
from Bio import SeqIO
import multiprocessing as mp


########################################################################################################################

# file in
gnm_id_txt          = '/Users/songweizhi/Desktop/Brady/SHC/gnm_id.txt'
blast_op_dir        = '/Users/songweizhi/Desktop/Brady/SHC/blastp_op_dir'
e_index_cutoff      = 30
shc_copy_num_txt    = '/Users/songweizhi/Desktop/Brady/SHC/shc_copy_num.txt'

########################################################################################################################

shc_copy_num_txt_handle = open(shc_copy_num_txt, 'w')
copy_num_dict = dict()
for each_gnm in open(gnm_id_txt):

    gnm_id            = each_gnm.strip()
    blast_op          = '%s/%s_blastp.txt'          % (blast_op_dir, gnm_id)
    copy_num_dict[gnm_id] = 0

    for each_hit in open(blast_op):
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
            copy_num_dict[gnm_id] += 1

    if copy_num_dict[gnm_id] == 2:
        print('%s\t%s' % (gnm_id, copy_num_dict[gnm_id]))

    shc_copy_num_txt_handle.write('%s\t%s\n' % (gnm_id, copy_num_dict[gnm_id]))

shc_copy_num_txt_handle.close()


'''

cd /Users/songweizhi/Desktop/Brady/SHC
BioSAK iTOL -SimpleBar -lv shc_copy_num.txt -scale 0-1-2 -lt SHC -out shc_copy_num_SimpleBar.txt

'''