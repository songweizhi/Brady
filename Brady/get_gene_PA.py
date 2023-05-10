import os
import argparse


def get_gene_pa(args):

    op_prefix           = args['p']
    blast_op_dir        = args['i']
    blast_op_ext        = args['x']
    genes_to_check_str  = args['g']
    evalue_index_cutoff = args['e']
    evalue_cutoff       = args['e']


    op_df_num = '%s_Num.txt' % op_prefix
    op_df_PA  = '%s_PA.txt'  % op_prefix

    # get sub dir list under blast_op_dir
    sub_dir_list = next(os.walk(blast_op_dir))[1]

    gene_list_to_check = genes_to_check_str.split(',')


    evalue_cutoff_dict = dict()




    op_num_handle = open(op_df_num, 'w')
    op_num_handle.write('Genome\t' + '\t'.join(gene_list_to_check) + '\n')
    op_PA_handle = open(op_df_PA, 'w')
    op_PA_handle.write('Genome\t' + '\t'.join(gene_list_to_check) + '\n')
    for each_gnm in sub_dir_list:

        identified_gene_num_dict = dict()
        for each_gene in gene_list_to_check:
            blast_op = '%s/%s/%s.%s' % (blast_op_dir, each_gnm, each_gene, blast_op_ext)

            gene_found_num = 0
            for each_hit in open(blast_op):
                each_hit_split = each_hit.strip().split('\t')
                e_value = each_hit_split[-2]
                if '-' in e_value:
                    e_value_index = int(e_value.split('-')[-1])
                    if e_value_index >= evalue_index_cutoff:
                        gene_found_num += 1
                elif float(e_value) == 0:
                    gene_found_num += 1

            identified_gene_num_dict[each_gene] = gene_found_num

        num_list = [identified_gene_num_dict[i] for i in gene_list_to_check]

        pa_list = []
        for each in num_list:
            if each > 0:
                pa_list.append('1')
            else:
                pa_list.append('0')

        num_list_as_str = [str(i) for i in num_list]

        num_str_to_write = '%s\t%s\n' % (each_gnm, '\t'.join(num_list_as_str))
        pa_str_to_write = '%s\t%s\n' % (each_gnm, '\t'.join(pa_list))
        op_num_handle.write(num_str_to_write)
        op_PA_handle.write(pa_str_to_write)

    op_num_handle.close()
    op_PA_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-p',     required=True,                    help='output prefix')
    parser.add_argument('-i',     required=True,                    help='blast results')
    parser.add_argument('-x',     required=True,                    help='file extension of blast results')
    parser.add_argument('-g',     required=True,                    help='genes to check')
    parser.add_argument('-e',     required=False, default=1e-100,   help='evalue cutoff, default: 1e-100')
    args = vars(parser.parse_args())
    get_gene_pa(args)



'''

# format of evalue_cutoffs.txt (tab separated)
nifA	1e-100
nifB	1e-50
nifD	0
nifE	0.0001

python3 /home-user/wzsong/Japonicum/get_gene_PA.py -p Japonicum_nod -i nod_blastp_outdir -x blast8 -g nodA,nodB,nodC,nodI,nodJ -e 1e-70
python3 /home-user/wzsong/Japonicum/get_gene_PA.py -p Japonicum_nod -i nod_blastp_outdir -x blast8 -g nodA,nodB,nodC,nodI,nodJ -e evalue_cutoffs.txt

'''