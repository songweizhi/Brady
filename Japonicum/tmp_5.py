import os
import glob
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests


def subset_df_by_cols(file_in, rows_to_keep, cols_to_keep, sep_symbol, column_name_pos, row_name_pos, file_out):

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


df_file     = '/Users/songweizhi/Desktop/demo.txt'
df_file_out = '/Users/songweizhi/Desktop/demo_subset.txt'
to_keep_col_list = ['PWY0-1578', 'P283-PWY']
to_keep_row_list = ['a1', 'a3', 'a5']

subset_df_by_cols(df_file, df_file_out, to_keep_col_list, to_keep_row_list, '\t', 0, 0)
