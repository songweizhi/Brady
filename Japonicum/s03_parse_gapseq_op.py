import os
import glob
import pandas as pd


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def rm_single_value_cols(df_in):

    # Check if there is only one unique value in each column
    unique_values = df_in.nunique()
    same_value_cols = unique_values[unique_values == 1].index

    # Drop columns with the same value
    df_out = df_in.drop(same_value_cols, axis=1)

    return df_out


def filter_col(df_in, min_value, min_value_min_num):

    df_in_copy = df_in.copy(deep=True)

    for col in df_in_copy.columns:
        count = (df_in_copy[col] >= min_value).sum()
        if count < min_value_min_num:
            df_in_copy.drop(col, axis=1, inplace=True)

    return df_in_copy


########################################################################################################################

'''
all			Pathways|seed|kegg
amino		Amino-Acid-Biosynthesis	
nucl		Nucleotide-Biosynthesis
cofactor	Cofactor-Biosynthesis
carbo		CARBO-BIOSYNTHESIS
carbo-deg	Carbohydrates-Degradation
polyamine	Polyamine-Biosynthesis
fatty		Fatty-acid-biosynthesis
energy		Energy-Metabolism
terpenoid	Terpenoid-Biosynthesis
degradation	Degradation
'''


# file in
wd                                  = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc'
key_list                            = ['amino', 'carbo', 'carbo-deg', 'cofactor', 'degradation', 'energy', 'fatty', 'nucl', 'polyamine', 'terpenoid']
file_ext                            = 'tbl'
absent_as_minus_one                 = True
min_value                           = 1
min_value_min_num                   = 10


for pathway_key in key_list:

    print('Processing %s' % pathway_key)

    #pathway_key                         = 'degradation'  # amino, nucl, cofactor, carbo, polyamine
    pathways_tbl_dir                    = '%s/Pathways_%s'                                              % (wd, pathway_key)

    # file out
    op_df                               = '%s/Pathways_PA_%s.txt'                                       % (wd, pathway_key)
    op_df_no_boring_cols                = '%s/Pathways_PA_%s_no_single_value_cols.txt'                  % (wd, pathway_key)
    op_df_no_boring_cols_itol           = '%s/Pathways_PA_%s_no_single_value_cols_iTOL.txt'             % (wd, pathway_key)
    op_df_no_boring_cols_filtered       = '%s/Pathways_PA_%s_no_single_value_cols_min%s_num%s.txt'      % (wd, pathway_key, min_value, min_value_min_num)
    op_df_no_boring_cols_filtered_itol  = '%s/Pathways_PA_%s_no_single_value_cols_min%s_num%s_iTOL.txt' % (wd, pathway_key, min_value, min_value_min_num)

    ########################################################################################################################

    pathways_tbl_file_re   = '%s/*.%s' % (pathways_tbl_dir, file_ext)
    pathways_tbl_file_list = glob.glob(pathways_tbl_file_re)

    # read in prediction
    all_pwy_set = set()
    gnm_to_detected_pwy_dict = dict()
    for pathways_tbl_file in pathways_tbl_file_list:
        f_path, f_base, f_ext = sep_path_basename_ext(pathways_tbl_file)
        gnm_id = f_base.split('-%s-' % pathway_key)[0]
        detected_pwy_set = set()
        for each_pwy in open(pathways_tbl_file):
            if (not each_pwy.startswith('#')) and (not each_pwy.startswith('ID')):
                each_pwy_split = each_pwy.strip().split('\t')
                all_pwy_set.add(each_pwy_split[0])
                if each_pwy_split[2] == 'true':
                    detected_pwy_set.add(each_pwy_split[0])
        gnm_to_detected_pwy_dict[gnm_id] = detected_pwy_set

    gnm_list_sorted = sorted(list(gnm_to_detected_pwy_dict.keys()))
    pwy_list_sorted = sorted([i for i in all_pwy_set])

    # write out prediction
    op_df_handle = open(op_df, 'w')
    op_df_handle.write('\t%s\n' % ('\t'.join(pwy_list_sorted)))
    for each_gnm in gnm_list_sorted:
        detected_pwy_set = gnm_to_detected_pwy_dict[each_gnm]
        detected_pwy_pa_list = []
        for each_pwy in pwy_list_sorted:
            if each_pwy in detected_pwy_set:
                detected_pwy_pa_list.append('1')
            else:
                if absent_as_minus_one is True:
                    detected_pwy_pa_list.append('-1')
                else:
                    detected_pwy_pa_list.append('0')
        op_df_handle.write('%s\t%s\n' % (each_gnm, '\t'.join(detected_pwy_pa_list)))
    op_df_handle.close()

    # remove single value columns
    df = pd.read_csv(op_df, sep='\t', header=0, index_col=0)
    df_without_single_value_cols = rm_single_value_cols(df)
    df_without_single_value_cols.to_csv(op_df_no_boring_cols, sep='\t')
    df_without_single_value_cols_filtered = filter_col(df_without_single_value_cols, min_value, min_value_min_num)
    df_without_single_value_cols_filtered.to_csv(op_df_no_boring_cols_filtered, sep='\t')

    # iTOL
    itol_cmd          = 'BioSAK iTOL -Binary -lm %s -lt %s -out %s' % (op_df_no_boring_cols, pathway_key, op_df_no_boring_cols_itol)
    itol_cmd_filtered = 'BioSAK iTOL -Binary -lm %s -lt %s -out %s' % (op_df_no_boring_cols_filtered, pathway_key, op_df_no_boring_cols_filtered_itol)
    os.system(itol_cmd)
    os.system(itol_cmd_filtered)

