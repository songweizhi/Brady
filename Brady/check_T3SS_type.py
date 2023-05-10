from Bio import SeqIO


def get_t3ss_core_gene_dict(t3ss_core_gene_info_txt):
    t3ss_type_to_core_gene_dict = dict()
    for each_line in open(t3ss_core_gene_info_txt):
        each_line_split = each_line.strip().split('\t')
        gene_id = each_line_split[0]
        t3_type = '_'.join(gene_id.split('_')[:-1])
        if len(each_line_split) == 3:
            if t3_type not in t3ss_type_to_core_gene_dict:
                t3ss_type_to_core_gene_dict[t3_type] = []
            t3ss_type_to_core_gene_dict[t3_type].append(gene_id)
    return t3ss_type_to_core_gene_dict


def get_t3ss_type(identified_gene_dir, core_gene_list, type_to_test, file_out):

    genome_to_identified_core_dict = dict()
    for each_core in core_gene_list:
        identified_seq = '%s/%s.fas' % (identified_gene_dir, each_core)
        for each_seq in SeqIO.parse(identified_seq, 'fasta'):
            seq_id = each_seq.id

            # add identified core to dict
            if seq_id not in genome_to_identified_core_dict:
                genome_to_identified_core_dict[seq_id] = []
            genome_to_identified_core_dict[seq_id].append(each_core)

    # write out
    file_out_handle = open(file_out, 'w')
    for each_gnm in genome_to_identified_core_dict:
        current_gnm_identified_core_list = genome_to_identified_core_dict[each_gnm]
        if len(current_gnm_identified_core_list) >= 4:
            file_out_handle.write('%s\t%s\n' % (each_gnm, type_to_test))
    file_out_handle.close()


# file in
identified_gene_dir = '/Users/songweizhi/Desktop/all'
core_gene_list      = ['S58_rhcC1', 'S58_rhcJ', 'S58_rhcN', 'S58_rhcU', 'S58_rhcV']
type_to_test        = 'S58'

# file out
file_out            = '/Users/songweizhi/Desktop/T3SS_S58.txt'

get_t3ss_type(identified_gene_dir, core_gene_list, type_to_test, file_out)
