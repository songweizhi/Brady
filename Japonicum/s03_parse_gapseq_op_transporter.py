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

wd                              = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc_transporter'
gap_seq_transport_stdout_dir    = '%s/transport_stdout'                                             % wd
gap_seq_transport_stdout_ext    = 'txt'
min_value                       = 1
min_value_min_num               = 10

# file out
op_df                           = '%s/transport_PA.txt'                                             % wd
op_df_no_boring_cols            = '%s/transport_PA_no_single_value_cols.txt'                        % wd
op_df_no_boring_cols_filtered   = '%s/transport_PA_no_single_value_cols_min%s_num%s.txt'            % (wd, min_value, min_value_min_num)

########################################################################################################################

stdout_file_re = '%s/*.%s' % (gap_seq_transport_stdout_dir, gap_seq_transport_stdout_ext)
stdout_file_list = glob.glob(stdout_file_re)

overall_found_transporter_set = set()
gnm_to_found_transporter_dict = dict()
for each_stdout in stdout_file_list:
    _, f_base, _ = sep_path_basename_ext(each_stdout)
    gnm_id = f_base.split('_find-transport_stdout')[0]
    to_count = False
    found_transporter_set = set()
    for each_line in open(each_stdout):
        if each_line.strip() == 'Found transporter for substances:':
            to_count = True
        if each_line.strip() == 'Found transporter without reaction in DB:':
            to_count = False
        if to_count is True:
            if each_line.strip() != '':
                found_transporter_set.add(each_line.strip())
                overall_found_transporter_set.add(each_line.strip())
    gnm_to_found_transporter_dict[gnm_id] = found_transporter_set

overall_found_transporter_list_sorted = sorted([i for i in overall_found_transporter_set])
gnm_list_sorted = sorted(list(gnm_to_found_transporter_dict.keys()))

op_df_txt_handle = open(op_df, 'w')
op_df_txt_handle.write('\t' + '\t'.join(overall_found_transporter_list_sorted) + '\n')
for each_gnm in gnm_list_sorted:
    pa_list = [each_gnm]
    current_found_transporter_set = gnm_to_found_transporter_dict[each_gnm]
    for each_transporter in overall_found_transporter_list_sorted:
        if each_transporter in current_found_transporter_set:
            pa_list.append('1')
        else:
            pa_list.append('0')
    pa_str_to_write = '\t'.join(pa_list)
    op_df_txt_handle.write(pa_str_to_write + '\n')
op_df_txt_handle.close()

# remove single value columns
df = pd.read_csv(op_df, sep='\t', header=0, index_col=0)
df_without_single_value_cols = rm_single_value_cols(df)
df_without_single_value_cols.to_csv(op_df_no_boring_cols, sep='\t')
df_without_single_value_cols_filtered = filter_col(df_without_single_value_cols, min_value, min_value_min_num)
df_without_single_value_cols_filtered.to_csv(op_df_no_boring_cols_filtered, sep='\t')

