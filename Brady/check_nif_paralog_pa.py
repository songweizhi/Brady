import os
import glob
from Bio import SeqIO


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext

########################################################################################################################

# file in
file_dir = '/Users/songweizhi/Desktop/nif_seq_updated'
file_ext = 'fas'

# file out
file_out_txt = '/Users/songweizhi/Desktop/nif_seq_updated.txt'

########################################################################################################################

file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)

gnm_to_identified_nif_dict= dict()
nif_set = set()
for each_nif in file_list:
    f_path, f_base, f_ext = sep_path_basename_ext(each_nif)
    nif_set.add(f_base)
    identified_gnm_set = set()
    for each_gnm in SeqIO.parse(each_nif, 'fasta'):
        gnm_id = each_gnm.id
        if '|' in gnm_id:
            identified_gnm_set.add(gnm_id)
            if gnm_id not in gnm_to_identified_nif_dict:
                gnm_to_identified_nif_dict[gnm_id] = set()
            gnm_to_identified_nif_dict[gnm_id].add(f_base)


nif_list_sorted = sorted([i for i in nif_set])


file_out_txt_handle = open(file_out_txt, 'w')

file_out_txt_handle.write('%s\t%s\n' % ('Genome', '\t'.join(nif_list_sorted)))
for each_gnm in gnm_to_identified_nif_dict:
    gnm_detected_nif_set = gnm_to_identified_nif_dict[each_gnm]
    nif_pa_list = []
    for each_nif in nif_list_sorted:
        if each_nif in gnm_detected_nif_set:
            nif_pa_list.append('1')
        else:
            nif_pa_list.append('0')
    file_out_txt_handle.write('%s\t%s\n' % (each_gnm, '\t'.join(nif_pa_list)))
file_out_txt_handle.close()

