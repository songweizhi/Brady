import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def perform_Benjamini_H_correction(p_values_list):

    p_value_array = np.array(p_values_list)
    _, adjusted_p_values, _, _ = multipletests(p_value_array, method='fdr_bh')

    adjust_p_value_dict = dict()
    for p_value, adj_p_value in zip(p_value_array, adjusted_p_values):
        if p_value not in adjust_p_value_dict:
            adjust_p_value_dict[p_value] = []
        adjust_p_value_dict[p_value] = adj_p_value

    return adjust_p_value_dict


def perform_fisher_exact_test(pa_table_txt, gnm_cate_txt, interested_gnm_txt):

    # read in interested_gnm
    interested_gnm_set = set()
    for each_gnm in open(interested_gnm_txt):

        gnm_id = each_gnm.strip()
        if '\t' in gnm_id:
            gnm_id = gnm_id.split('\t')[0]
        elif ' ' in gnm_id:
            gnm_id = gnm_id.split(' ')[0]
        elif ',' in gnm_id:
            gnm_id = gnm_id.split(',')[0]
        interested_gnm_set.add(gnm_id)

    # read in gnm_cate_txt
    gnm_cate_set = set()
    gnm_to_cate_dict = dict()
    for each_gnm in open(gnm_cate_txt):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_cate_set.add(each_gnm_split[1])
        gnm_to_cate_dict[each_gnm_split[0]] = each_gnm_split[1]

    gnm_cate_list_sorted = sorted([i for i in gnm_cate_set])

    df_in = pd.read_csv(pa_table_txt, sep='\t', header=0, index_col=0)

    p_value_dict = dict()
    for (columnName, columnData) in df_in.iteritems():

        count_dict = dict()
        for (genome_id, pa_value) in zip(df_in.index, columnData.values):
            if genome_id in interested_gnm_set:
                gnm_cate = gnm_to_cate_dict[genome_id]

                if gnm_cate not in count_dict:
                    count_dict[gnm_cate] = dict()
                    count_dict[gnm_cate]['p'] = 0
                    count_dict[gnm_cate]['a'] = 0

                if pa_value == 1:
                    count_dict[gnm_cate]['p'] += 1

                if pa_value == 0:
                    count_dict[gnm_cate]['a'] += 1

        # get contingency table
        obs_table = [[count_dict[gnm_cate_list_sorted[0]]['p'],
                      count_dict[gnm_cate_list_sorted[0]]['a']],
                     [count_dict[gnm_cate_list_sorted[1]]['p'],
                      count_dict[gnm_cate_list_sorted[1]]['a']]]
        fisher_res = fisher_exact(obs_table)
        p_value = fisher_res[1]

        # add to dict
        p_value_dict[columnName] = p_value

    return p_value_dict


def perform_differential_presence_analysis(pa_table_txt, gnm_cate_txt, interested_gnm_txt, p_value_cutoff, output_txt):

    # perform_fisher_exact_test
    p_value_dict = perform_fisher_exact_test(pa_table_txt, gnm_cate_txt, interested_gnm_txt)

    # perform_Benjamini_H_correction
    adjust_p_value_dict = perform_Benjamini_H_correction(list(p_value_dict.values()))

    # write out results
    op_txt_handle = open(output_txt, 'w')
    for each_item in p_value_dict:
        p_value = p_value_dict[each_item]
        p_value_adjusted = adjust_p_value_dict[p_value]
        if p_value_adjusted <= p_value_cutoff:
            op_txt_handle.write('%s\t%s\t%s\n' % (each_item, p_value, p_value_adjusted))
    op_txt_handle.close()


# file in
pathways_pa_energy_txt  = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc/Pathways_PA_energy.txt'
gnm_cate_txt            = '/Users/songweizhi/Desktop/Japonicum/genome_cate_by_nod.txt'
interested_gnm_txt      = '/Users/songweizhi/Desktop/Japonicum/genome_cate_by_nod.txt'
p_value_cutoff          = 0.05

# file out
op_txt                  = '/Users/songweizhi/Desktop/differential_gene_presence_analysis.txt'


perform_differential_presence_analysis(pathways_pa_energy_txt, gnm_cate_txt, interested_gnm_txt, p_value_cutoff, op_txt)

