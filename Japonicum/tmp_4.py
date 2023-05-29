import numpy as np
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


p_values_list = [0.01, 0.05, 0.1, 0.05, 0.2, 0.3, 0.01, 0.05, 0.1, 0.05, 0.2, 0.3]
adjust_p_value_dict = perform_Benjamini_H_correction(p_values_list)

print(adjust_p_value_dict)


gapseq_db_meta_pwy_tbl              = '/Users/songweizhi/DB/gapseq/meta_pwy.tbl'








