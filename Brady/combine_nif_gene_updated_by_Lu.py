import os
from Bio import SeqIO
import argparse

def combined_nif_gene(args):
    
    nif_seq_from_ss_dir           = args['Sishuo_dir']
    new_nif_seq_dir               = args['New_nif_seq']
    special_strain_nif_seq_dir    = args['Special']
    updated_nif_seq_dir           = args['updated_nif_seq']
    special_strain_id_list        = args['special_list']
    nif_gene_list                 = ['nifA', 'nifB', 'nifD', 'nifE', 'nifH', 'nifK', 'nifN', 'nifX']


    for each_nif in nif_gene_list:

        pwd_nif_ss = '%s/%s.fas' % (nif_seq_from_ss_dir, each_nif)
        pwd_nif_new = '%s/%s.fas' % (new_nif_seq_dir, each_nif)
        pwd_nif_special = '%s/%s.fas' % (special_strain_nif_seq_dir, each_nif)
        pwd_nif_updated = '%s/%s.fas' % (updated_nif_seq_dir, each_nif)

        cp_cmd = 'cp %s %s' % (pwd_nif_ss, pwd_nif_updated)
        os.system(cp_cmd)

        with open(pwd_nif_updated, 'a') as pwd_nif_updated_handle:

            for each_seq in SeqIO.parse(pwd_nif_special, 'fasta'):
                seq_id = each_seq.id
                seq_seq = str(each_seq.seq)
                pwd_nif_updated_handle.write('>%s\n' % seq_id)
                pwd_nif_updated_handle.write('%s\n' % seq_seq)

            for each_seq in SeqIO.parse(pwd_nif_new, 'fasta'):
                seq_id = each_seq.id
                seq_seq = str(each_seq.seq)
                if seq_id not in special_strain_id_list:
                    pwd_nif_updated_handle.write('>%s\n' % seq_id)
                    pwd_nif_updated_handle.write('%s\n' % seq_seq)


if __name__ == '__main__':

    # initialize the options parser
    parser = argparse.ArgumentParser()

    # arguments for combined nif genes
    parser.add_argument('-Sishuo_dir',         required=True,     help='1070 strains with paralog sequence file')
    parser.add_argument('-New_nif_seq',        required=True,     help='newly isolated strains with nif genes')
    parser.add_argument('-Special',            required=True,     help='newly isolated strains with paralog nif genes')
    parser.add_argument('-updated_nif_seq',    required=True,     help='output combined file')
    parser.add_argument('-special_list',       required=False,    help='The special strains with paralog genes ID')
    args = vars(parser.parse_args())
    combined_nif_gene(args)


'''
python3 

'''