import os
import glob
import numpy as np
import pandas as pd
from ete3 import Tree
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def rm_single_value_cols(df_in):

    # Check if there is only one unique value in each column
    unique_values = df_in.nunique()
    same_value_cols = unique_values[unique_values == 1].index

    # Drop columns with the same value
    df_out = df_in.drop(same_value_cols, axis=1)

    return df_out


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


def perform_fisher_exact_test(perform_stats, pa_table_txt, gnm_to_cate_dict, interested_gnm_txt):

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

    gnm_cate_list_sorted = sorted(list(gnm_to_cate_dict.keys()))

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


def perform_binaryPGLMM_test(perform_stats, tree_file, pa_table_txt, gnm_to_cate_dict, interested_gnm_txt, PhyloBiAssoc_R):

    f_path, f_base, f_ext       = sep_path_basename_ext(interested_gnm_txt)
    binaryPGLMM_input_txt_tmp   = '%s/binaryPGLMM_%s_input_tmp.txt'  % (f_path, f_base)
    binaryPGLMM_input_txt       = '%s/binaryPGLMM_%s_input.txt'      % (f_path, f_base)
    binaryPGLMM_output_txt      = '%s/binaryPGLMM_%s_output.txt'     % (f_path, f_base)

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

    if perform_stats is True:
        binaryPGLMM_input_txt_handle = open(binaryPGLMM_input_txt_tmp, 'w')
        for each_gnm in open(pa_table_txt):
            each_gnm_split = each_gnm.strip().split('\t')
            if each_gnm.startswith('\t'):
                binaryPGLMM_input_txt_handle.write('ID\tcate\t%s\n' % '\t'.join(each_gnm_split))
            gnm_id = each_gnm_split[0]
            if gnm_id in interested_gnm_set:
                gnm_cate = gnm_to_cate_dict[gnm_id]
                binaryPGLMM_input_txt_handle.write('%s\t%s\t%s\n' % (gnm_id, gnm_cate, '\t'.join(each_gnm_split[1:])))
        binaryPGLMM_input_txt_handle.close()

        df = pd.read_csv(binaryPGLMM_input_txt_tmp, sep='\t', header=0, index_col=0)
        df_without_single_value_cols = rm_single_value_cols(df)
        df_without_single_value_cols.to_csv(binaryPGLMM_input_txt, sep='\t')
        os.system('rm %s' % binaryPGLMM_input_txt_tmp)
        binaryPGLMM_cmd = 'Rscript %s -t %s -d %s > %s' % (PhyloBiAssoc_R, tree_file, binaryPGLMM_input_txt, binaryPGLMM_output_txt)
        os.system(binaryPGLMM_cmd)

    p_value_dict = dict()
    for each_test in open(binaryPGLMM_output_txt):
        each_test_split = each_test.strip().split('\t')
        pwy_id = each_test_split[0].replace('.', '-')
        p_value = float(each_test_split[-1])
        p_value_dict[pwy_id] = p_value

    return p_value_dict


def perform_differential_presence_analysis(perform_stats, test_to_perform, tree_file_subset, pa_table_txt, gnm_to_cate_dict, interested_gnm_txt, p_value_cutoff, output_txt, PhyloBiAssoc_R, perform_multiple_test_correction):

    if test_to_perform == 'FisherExact':
        p_value_dict = perform_fisher_exact_test(perform_stats, pa_table_txt, gnm_to_cate_dict, interested_gnm_txt)
    elif test_to_perform == 'BinaryPGLMM':
        p_value_dict = perform_binaryPGLMM_test(perform_stats, tree_file_subset, pa_table_txt, gnm_to_cate_dict, interested_gnm_txt, PhyloBiAssoc_R)
    else:
        print('Please specify either FisherExact or BinaryPGLMM, program exited!')
        exit()

    # perform_Benjamini_H_correction
    adjust_p_value_dict = perform_Benjamini_H_correction(list(p_value_dict.values()))

    # write out results
    op_txt_handle = open(output_txt, 'w')
    for each_item in p_value_dict:
        p_value = p_value_dict[each_item]
        p_value_adjusted = adjust_p_value_dict[p_value]

        if perform_multiple_test_correction is False:
            p_value_to_use = p_value
        else:
            p_value_to_use = p_value_adjusted

        if p_value_to_use <= p_value_cutoff:
            op_txt_handle.write('%s\t%s\n' % (each_item, p_value_to_use))
    op_txt_handle.close()


########################################################################################################################

pwd_transporter_pa_txt              = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc_transporter/transport_PA.txt'
analysis_clades_txt                 = '/Users/songweizhi/Desktop/Japonicum/Japo_analysis_clades.txt'
tree_file                           = '/Users/songweizhi/Desktop/Japonicum/Japonicum_121_OG_tree_LG.treefile'
gnm_cate_txt                        = '/Users/songweizhi/Desktop/Japonicum/genome_cate_by_nod.txt'
p_value_cutoff                      = 0.05
test_to_perform                     = 'BinaryPGLMM' # FisherExact, BinaryPGLMM
PhyloBiAssoc_R                      = '/Users/songweizhi/PycharmProjects/BioSAK/BioSAK/PhyloBiAssoc.R'
perform_PhyloBiAssoc_stats          = True
perform_multiple_test_correction    = True
gnm_to_ignore_list                  = []

# op dir
op_dir                              = '/Users/songweizhi/Desktop/Japonicum/differential_presence_analysis_by_clade_transporter'
pathway_color_txt                   = '%s/pathway_color.txt' % op_dir

########################################################################################################################

# read in gnm_cate_txt
gnm_cate_set = set()
gnm_to_cate_dict = dict()
for each_gnm in open(gnm_cate_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_cate = each_gnm_split[1]
    if gnm_id not in gnm_to_ignore_list:
        gnm_cate_set.add(gnm_cate)
        gnm_to_cate_dict[gnm_id] = gnm_cate

# read in analysis_clades_txt
analysis_clade_to_gnm_dict = dict()
for each_japo in open(analysis_clades_txt):
    each_japo_split = each_japo.strip().split('\t')
    gnm_id = each_japo_split[0]
    ac_id = each_japo_split[1]
    if ac_id not in analysis_clade_to_gnm_dict:
        analysis_clade_to_gnm_dict[ac_id] = set()
    if gnm_id not in gnm_to_ignore_list:
        analysis_clade_to_gnm_dict[ac_id].add(gnm_id)

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

    pwd_analysis_result_txt           = '%s/%s_transporter_PA_analysis_result.txt'  % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff           = '%s/%s_transporter_PA_diff.txt'             % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff_0as1      = '%s/%s_transporter_PA_diff_0as1.txt'        % (op_dir, analysis_clade)
    pwd_pathway_pa_txt_diff_0as1_itol = '%s/%s_transporter_PA_diff_iTOL.txt'        % (op_dir, analysis_clade)
    pwd_pathway_color_txt             = '%s/%s_transporter_color.txt'               % (op_dir, analysis_clade)

    # perform_differential_presence_analysis
    perform_differential_presence_analysis(perform_PhyloBiAssoc_stats, test_to_perform, tree_file_subset, pwd_transporter_pa_txt, gnm_to_cate_dict, interested_gnm_txt, p_value_cutoff, pwd_analysis_result_txt, PhyloBiAssoc_R, perform_multiple_test_correction)

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
        subset_df(pwd_transporter_pa_txt, interested_gnm_list_sorted, diff_pwy_id_list_sorted, '\t', 0, 0, pwd_pathway_pa_txt_diff)

        ####################################################################################################################

        # write out pathway_color_txt
        all_pwy_list = open(pwd_pathway_pa_txt_diff).readline().strip().split('\t')
        pwd_pathway_color_txt_handle = open(pwd_pathway_color_txt, 'w')
        for each_pwy in all_pwy_list:
            if '__' in each_pwy:
                pwy_cate = each_pwy.split('__')[0]
                pwd_pathway_color_txt_handle.write('%s\t%s\n' % (each_pwy, 'lightgreen'))
        pwd_pathway_color_txt_handle.close()

        # turn 0 to -1
        pwd_pathway_pa_txt_diff_0as1_handle = open(pwd_pathway_pa_txt_diff_0as1, 'w')
        for each_line in open(pwd_pathway_pa_txt_diff):
            pwd_pathway_pa_txt_diff_0as1_handle.write(each_line.replace('\t0', '\t-1'))
        pwd_pathway_pa_txt_diff_0as1_handle.close()

        # iTOL
        itol_cmd = 'BioSAK iTOL -Binary -lm %s -lt Pathway -cc %s -o %s' % (pwd_pathway_pa_txt_diff_0as1, pwd_pathway_color_txt, pwd_pathway_pa_txt_diff_0as1_itol)
        print(itol_cmd)
        os.system(itol_cmd)
