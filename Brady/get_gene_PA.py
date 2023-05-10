import os
import argparse


def numericable_str(my_string):

    is_numericable_str = False
    try:
        my_numeric_value = int(my_string)
        is_numericable_str = True
    except ValueError:
        try:
            my_numeric_value = float(my_string)
            is_numericable_str = True
        except ValueError:
            is_numericable_str = False

    return is_numericable_str


def get_gene_pa(args):

    op_prefix           = args['p']
    blast_op_dir        = args['i']
    blast_op_ext        = args['x']
    genes_to_check_str  = args['g']
    evalue_cutoff       = args['e']

    op_df_num = '%s_Num.txt' % op_prefix
    op_df_PA  = '%s_PA.txt'  % op_prefix

    # get sub dir list under blast_op_dir
    sub_dir_list = next(os.walk(blast_op_dir))[1]
    gene_list_to_check = genes_to_check_str.split(',')

    # get evalue_cutoff_dict
    evalue_cutoff_dict_tmp = dict()
    if os.path.isfile(evalue_cutoff) is False:
        if numericable_str(evalue_cutoff) is False:
            print('invalid e-value cutoff or wrong file path/name, please double check')
            exit()
        else:
            for each_g in gene_list_to_check:
                evalue_cutoff_dict_tmp[each_g] = float(evalue_cutoff)
    else:
        for each_cutoff in open(evalue_cutoff):
            each_cutoff_split = each_cutoff.strip().split('\t')
            if numericable_str(each_cutoff_split[1]) is False:
                print('invalid e-value cutoff, please double check:\n%s' % each_cutoff.strip())
                exit()
            else:
                evalue_cutoff_dict_tmp[each_cutoff_split[0]] = float(each_cutoff_split[1])

    evalue_cutoff_dict = dict()
    for each_gene in gene_list_to_check:
        evalue_cutoff_dict[each_gene] = evalue_cutoff_dict_tmp.get(each_gene, 1e-100)

    print('The following evalue cutoffs will be applied:')
    for each_gene in evalue_cutoff_dict:
        if each_gene not in evalue_cutoff_dict_tmp:
            print('%s\t%s (did not provide, default value will be applied)' % (each_gene, evalue_cutoff_dict[each_gene]))
        else:
            print('%s\t%s' % (each_gene, evalue_cutoff_dict[each_gene]))

    op_num_handle = open(op_df_num, 'w')
    op_num_handle.write('Genome\t' + '\t'.join(gene_list_to_check) + '\n')
    op_PA_handle = open(op_df_PA, 'w')
    op_PA_handle.write('Genome\t' + '\t'.join(gene_list_to_check) + '\n')
    for each_gnm in sub_dir_list:
        identified_gene_num_dict = dict()
        for each_gene in gene_list_to_check:
            current_evalue_cutoff =  evalue_cutoff_dict[each_gene]
            blast_op = '%s/%s/%s.%s' % (blast_op_dir, each_gnm, each_gene, blast_op_ext)

            gene_found_num = 0
            for each_hit in open(blast_op):
                each_hit_split = each_hit.strip().split('\t')
                e_value = float(each_hit_split[-2])
                if e_value <= current_evalue_cutoff:
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

    print('Results exported to %s and %s.' % (op_df_num, op_df_PA))
    print('Done!')


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
cd /Users/songweizhi/Desktop/get_PA_wd
python3 get_gene_PA.py -p nod -i nod_blastp_outdir -x blast8 -g nodA,nodB,nodC,nodI,nodJ -e 1e-70
python3 get_gene_PA.py -p nod -i nod_blastp_outdir -x blast8 -g nodA,nodB,nodC,nodI,nodJ -e evalue_cutoffs.txt

# format of evalue_cutoffs.txt (tab separated)
nodA	1e-100
nodB	1e-50
nodC	0
nodI	0.0001
nodJ	1e-100

'''
