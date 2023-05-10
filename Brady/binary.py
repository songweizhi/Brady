
fd_txt          = '/Users/songweizhi/Desktop/HpnHGIOKJ_vs_394_PB_blastp_best_hit_df_fmted.txt'
itol_binary_txt = '/Users/songweizhi/Desktop/HpnHGIOKJ_vs_394_PB_blastp_best_hit_df_iTOL.txt'


itol_binary_txt_handle = open(itol_binary_txt, 'w')

col_name_list = []
line_index = 0
for each_line in open(fd_txt):
    each_line_split = each_line.strip().split('\t')
    if line_index == 0:
        col_name_list = each_line_split[1:]
        itol_binary_txt_handle.write('DATASET_BINARY\n\nSEPARATOR TAB\nDATASET_LABEL\tlabel1\nCOLOR\t#85C1E9\n')
        itol_binary_txt_handle.write('SHOW_LABELS\t1\nLABEL_ROTATION\t45\nLABEL_SHIFT\t5\n')
        itol_binary_txt_handle.write('FIELD_LABELS\t%s\n' % '\t'.join(col_name_list))
        itol_binary_txt_handle.write('FIELD_SHAPES\t%s\n' % '\t'.join(['1'] * len(col_name_list)))
        itol_binary_txt_handle.write('\nDATA\n')
    else:
        itol_binary_txt_handle.write(each_line)
    line_index += 1
itol_binary_txt_handle.close()
