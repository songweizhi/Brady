import os
import glob
import argparse
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def rm_problematic_gnms(args):

    og_dir              = args['i']
    og_ext              = args['x']
    problematic_gnm_txt = args['g']
    og_dir_updated      = args['o']
    force_mkdir         = args['f']

    if os.path.isdir(og_dir_updated) is True:
        if force_mkdir is False:
            print('%s detected, program exited!' % og_dir_updated)
            exit()
        else:
            os.system('rm -r %s' % og_dir_updated)
    os.system('mkdir %s' % og_dir_updated)

    problematic_gnm_set = set()
    for each_gnm in open(problematic_gnm_txt):
        problematic_gnm_set.add(each_gnm.strip())

    og_file_re = '%s/*.%s' % (og_dir, og_ext)
    og_file_list = glob.glob(og_file_re)

    for og_file in og_file_list:
        file_path, file_basename, file_ext = sep_path_basename_ext(og_file)
        og_file_updated = '%s/%s%s' % (og_dir_updated, file_basename, file_ext)
        og_file_updated_handle = open(og_file_updated, 'w')
        for each_seq in SeqIO.parse(og_file, 'fasta'):
            seq_id = each_seq.id
            if seq_id not in problematic_gnm_set:
                og_file_updated_handle.write('>%s\n' % each_seq.id)
                og_file_updated_handle.write('%s\n' % str(each_seq.seq))
        og_file_updated_handle.close()


########################################################################################################################

# # file in
# og_dir              = '/Users/songweizhi/Desktop/LingHK360_Zong1070_123_OG'
# og_ext              = 'fas'
# problematic_gnm_txt = '/Users/songweizhi/Desktop/strains_to_ignore.txt'
# force_mkdir         = True
#
# # file out
# og_dir_updated      = '/Users/songweizhi/Desktop/LingHK360_Zong1070_123_OG_problematic_gnm_removed'

########################################################################################################################

# # file in
# og_dir              = '/home-user/wzsong/Brady/subsampling_nif_tree/pep-clean_subsampled'
# og_ext              = 'fas'
# gnms_to_remove_txt  = '/home-user/wzsong/Brady/subsampling_nif_tree/distant_outgroup.txt'
# force_mkdir         = True
#
# # file out
# og_dir_updated      = '/home-user/wzsong/Brady/subsampling_nif_tree/pep-clean_subsampled_distant_outgroup_removed2'

########################################################################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',   required=True,                          help='sequence file dir')
    parser.add_argument('-x',   required=False, default='faa',          help='sequence file extension, default:faa')
    parser.add_argument('-g',   required=True,                          help='genome id txt')
    parser.add_argument('-o',   required=True,                          help='output dir')
    parser.add_argument('-f',   required=False, action="store_true",    help='force overwrite existing output')
    args = vars(parser.parse_args())
    rm_problematic_gnms(args)

'''
cd /home-user/wzsong/Brady/nif_tree_2023_05_09_tree1
python3 /home-user/wzsong/Scripts/rm_problematic_gnms.py -i nif_seq_updated -x fas -g gnms_to_ignore.txt -o nif_seq_updated2 -f 
BioSAK get_gnm_size -i nif_seq_updated -x fas -total
BioSAK get_gnm_size -i nif_seq_updated2 -x fas -total

cd /home-user/wzsong/Brady/nif_tree_2023_05_09_tree2
python3 /home-user/wzsong/Scripts/rm_problematic_gnms.py -i nif_seq_updated -x fas -g gnms_to_ignore.txt -o nif_seq_updated2 -f 

'''