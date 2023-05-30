import os
import glob
import numpy as np
import pandas as pd
from ete3 import Tree
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def subset_tree(tree_file_in, genomes_to_keep, tree_file_out):
    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(genomes_to_keep, preserve_branch_length=True)
    subset_tree.write(outfile=tree_file_out)


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def subset_df(file_in, rows_to_keep, cols_to_keep, sep_symbol, column_name_pos, row_name_pos, file_out):
    print(file_in)
    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)

    if len(rows_to_keep) == 0:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[:, :]
        else:
            subset_df = df.loc[:, cols_to_keep]
    else:
        if len(cols_to_keep) == 0:
            subset_df = df.loc[rows_to_keep, :]
        else:
            subset_df = df.loc[rows_to_keep, cols_to_keep]

    subset_df.to_csv(file_out, sep=sep_symbol)


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


########################################################################################################################

pwd_pathway_pa_txt      = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc/Pathway_PA.txt'
gapseq_db_meta_pwy_tbl  = '/Users/songweizhi/DB/gapseq/meta_pwy.tbl'
gnm_cate_txt            = '/Users/songweizhi/Desktop/Japonicum/differential_gene_presence_analysis_with_name/genome_cate_by_nod.txt'
pwy_cate_color_txt      = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc/pwy_cate_color.txt'
p_value_cutoff          = 0.05
add_pathway_name        = True
analysis_clades_txt     = '/Users/songweizhi/Desktop/Japonicum/Japo_analysis_clades.txt'
tree_file               = '/Users/songweizhi/Desktop/Japonicum/Japonicum_121_OG_tree_LG.treefile'
gnm_to_ignore_list      = ['GCA_024171065.1']

# op dir
op_dir                  = '/Users/songweizhi/Desktop/Japonicum/differential_gene_presence_analysis_by_clade'
pathway_color_txt       = '%s/pathway_color.txt' % op_dir

########################################################################################################################

# read in analysis_clades_txt
analysis_clade_to_gnm_dict = dict()
for each_japo in open(analysis_clades_txt):
    each_japo_split = each_japo.strip().split('\t')
    gnm_id = each_japo_split[0]
    ac_id = each_japo_split[1]
    if ac_id not in analysis_clade_to_gnm_dict:
        analysis_clade_to_gnm_dict[ac_id] = set()
    analysis_clade_to_gnm_dict[ac_id].add(gnm_id)

# read in gapseq's meta_pwy.tbl
meta_pwy_id2name_dict = dict()
for each_pwy in open(gapseq_db_meta_pwy_tbl):
    if not each_pwy.startswith('id\t'):
        each_pwy_split = each_pwy.strip().split('\t')
        pwy_id = each_pwy_split[0]
        pwy_name = each_pwy_split[1]
        meta_pwy_id2name_dict[pwy_id.replace('|', '')] = pwy_name.replace(' ', '_').replace(',', '_')

# read in color
pwy_cate_color_dict = dict()
if os.path.isfile(pwy_cate_color_txt) is True:
    for each_pwy_cate in open(pwy_cate_color_txt):
        each_pwy_cate_split = each_pwy_cate.strip().split('\t')
        pwy_cate_color_dict[each_pwy_cate_split[0]] = each_pwy_cate_split[1]

for analysis_clade in analysis_clade_to_gnm_dict:

    print('Processing %s' % analysis_clade)

    ac_gnm_set          = analysis_clade_to_gnm_dict[analysis_clade]
    tree_file_subset    = '%s/Genome_tree_%s.treefile' % (op_dir, analysis_clade)
    interested_gnm_txt  = '%s/gnm_id_%s.txt'           % (op_dir, analysis_clade)

    # write out genome id to txt
    with open(interested_gnm_txt, 'w')  as gnm_id_txt_handle:
        gnm_id_txt_handle.write('\n'.join(ac_gnm_set))

    # subset genome tree
    subset_tree(tree_file, ac_gnm_set, tree_file_subset)

    pwd_analysis_result_txt           = '%s/%s_pathway_PA_analysis_result.txt'      % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff           = '%s/%s_pathway_PA_diff.txt'                 % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff_with_name = '%s/%s_pathway_PA_diff_with_name.txt'       % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff_0as1      = '%s/%s_pathway_PA_diff_with_name_0as1.txt'  % (op_dir, analysis_clade)
    pwd_pathway_color_txt             = '%s/%s_pathway_color.txt'                   % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff_0as1_itol = '%s/%s_pathway_PA_diff_with_name_iTOL.txt'  % (op_dir, analysis_clade)

    # perform_differential_presence_analysis
    perform_differential_presence_analysis(pwd_pathway_pa_txt, gnm_cate_txt, interested_gnm_txt, p_value_cutoff, pwd_analysis_result_txt)

    # get pathways with differential presence
    diff_pwy_id_set = set()
    for each_line in open(pwd_analysis_result_txt):
        each_line_split = each_line.strip().split('\t')
        diff_pwy_id_set.add(each_line_split[0])

    diff_pwy_id_list_sorted = sorted([i for i in diff_pwy_id_set])

    # get pathways with differential presence
    interested_gnm_set = set()
    for each_line in open(interested_gnm_txt):
        if each_line.strip() not in gnm_to_ignore_list:
            interested_gnm_set.add(each_line.strip())

    interested_gnm_list_sorted = sorted([i for i in interested_gnm_set])

    if len(diff_pwy_id_set) > 0:

        # Subset PA df by selected columns
        subset_df(pwd_pathway_pa_txt, interested_gnm_list_sorted, diff_pwy_id_list_sorted, '\t', 0, 0, pwd_pathway_pa_txt_diff)

        ####################################################################################################################

        # add name to pathway id
        if add_pathway_name is True:
            pwd_pathway_pa_txt_diff_with_name_handle = open(pwd_pathway_pa_txt_diff_with_name, 'w')
            for each_line in open(pwd_pathway_pa_txt_diff):
                each_line_split = each_line.strip().split('\t')
                if each_line.startswith('\t'):
                    header_list_with_name = ['']
                    for each_pwy in each_line_split:
                        pwy_id_no_cate = each_pwy.split('__')[1]
                        pwy_name = meta_pwy_id2name_dict.get(pwy_id_no_cate, '')
                        if pwy_name == '':
                            pwy_with_name = each_pwy
                        else:
                            pwy_with_name = '%s__%s' % (each_pwy, pwy_name)
                        header_list_with_name.append(pwy_with_name)
                    pwd_pathway_pa_txt_diff_with_name_handle.write('%s\n' % '\t'.join(header_list_with_name))
                else:
                    pwd_pathway_pa_txt_diff_with_name_handle.write(each_line)
            pwd_pathway_pa_txt_diff_with_name_handle.close()

        # write out pathway_color_txt
        all_pwy_list = open(pwd_pathway_pa_txt_diff_with_name).readline().strip().split('\t')
        pwd_pathway_color_txt_handle = open(pwd_pathway_color_txt, 'w')
        for each_pwy in all_pwy_list:
            if '__' in each_pwy:
                pwy_cate = each_pwy.split('__')[0]
                pwy_cate_color = pwy_cate_color_dict[pwy_cate]
                pwd_pathway_color_txt_handle.write('%s\t%s\n' % (each_pwy, pwy_cate_color))
        pwd_pathway_color_txt_handle.close()


        # turn 0 to -1
        pwd_pathway_pa_txt_diff_0as1_handle = open(pwd_pathway_pa_txt_diff_0as1, 'w')
        for each_line in open(pwd_pathway_pa_txt_diff_with_name):
            pwd_pathway_pa_txt_diff_0as1_handle.write(each_line.replace('\t0', '\t-1'))
        pwd_pathway_pa_txt_diff_0as1_handle.close()

        # iTOL
        itol_cmd = 'BioSAK iTOL -Binary -lm %s -lt Pathway -cc %s -o %s' % (pwd_pathway_pa_txt_diff_0as1, pwd_pathway_color_txt, pwd_pathway_pa_txt_diff_0as1_itol)
        os.system(itol_cmd)

        # remove tmp files
        os.system('rm %s' % pwd_pathway_color_txt)
        os.system('rm %s' % pwd_pathway_pa_txt_diff)
        os.system('rm %s' % pwd_pathway_pa_txt_diff_0as1)
