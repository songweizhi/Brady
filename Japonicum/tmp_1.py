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


file_dir = '/home-user/wzsong/Japonicum/data/genome'
#file_dir = '/Users/songweizhi/Desktop/genome'
file_ext = 'genome'

file_re = '%s/*.%s' % (file_dir, file_ext)
file_list = glob.glob(file_re)
print(file_list)

good_gnm_num = 0
for each_gnm in file_list:

    f_path, f_base, f_ext = sep_path_basename_ext(each_gnm)

    seq_len_dict = dict()
    for each_seq in SeqIO.parse(each_gnm, 'fasta'):
        seq_len_dict[each_seq.id] = len(each_seq.seq)

    if len(seq_len_dict) <= 5:
        good_gnm_num += 1
        print('%s\t%s' % (f_base, seq_len_dict))

print(len(file_list))
print(good_gnm_num)



print(30000*0.8868)
print(30000*0.9017)
