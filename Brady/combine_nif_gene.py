import os
from Bio import SeqIO


nif_seq_from_ss_dir           = '/Users/songweizhi/Desktop/555/1070_nif_seq_with_paralog_from_Sishuo'
new_nif_seq_dir               = '/Users/songweizhi/Desktop/555/extracted_nif_gene_seq'
special_strain_nif_seq_dir    = '/Users/songweizhi/Desktop/555/all'
updated_nif_seq_dir           = '/Users/songweizhi/Desktop/555/nif_seq_updated'
nif_gene_list                 = ['nifA', 'nifB', 'nifD', 'nifE', 'nifH', 'nifK', 'nifN', 'nifX']
special_strain_id_list        = ['Bradyrhizobium_iriomotense_HKCCYLS1036', 'Bradyrhizobium_oligotrophicum_HKCCYLR1017', 'Bradyrhizobium_oligotrophicum_HKCCYLR10141']


for each_nif in nif_gene_list:

    pwd_nif_ss = '%s/%s.fas' % (nif_seq_from_ss_dir, each_nif)
    pwd_nif_new = '%s/%s.fas' % (new_nif_seq_dir, each_nif)
    pwd_nif_special = '%s/%s.fas' % (special_strain_nif_seq_dir, each_nif)
    pwd_nif_updated = '%s/%s.fas' % (updated_nif_seq_dir, each_nif)

    cp_cmd = 'cp %s %s' % (pwd_nif_ss, pwd_nif_updated)
    os.system(cp_cmd)

    gnms_in_nif_ss = set()
    for each_seq in SeqIO.parse(pwd_nif_ss, 'fasta'):
        seq_id = each_seq.id
        if '|' not in seq_id:
            gnms_in_nif_ss.add(seq_id)
        else:
            gnm_id = seq_id.split('|')[0]
            gnms_in_nif_ss.add(gnm_id)

    with open(pwd_nif_updated, 'a') as pwd_nif_updated_handle:

        for each_seq in SeqIO.parse(pwd_nif_special, 'fasta'):
            seq_id = each_seq.id
            seq_seq = str(each_seq.seq)
            gnm_id = seq_id.split('|')[0]
            if gnm_id not in gnms_in_nif_ss:
                pwd_nif_updated_handle.write('>%s\n' % seq_id)
                pwd_nif_updated_handle.write('%s\n' % seq_seq)

        for each_seq in SeqIO.parse(pwd_nif_new, 'fasta'):
            seq_id = each_seq.id
            seq_seq = str(each_seq.seq)

            if seq_id not in special_strain_id_list:
                if seq_id not in gnms_in_nif_ss:
                    pwd_nif_updated_handle.write('>%s\n' % seq_id)
                    pwd_nif_updated_handle.write('%s\n' % seq_seq)
