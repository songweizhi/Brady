
name_correlation_txt = '/Users/songweizhi/Desktop/Brady/hkcct2species.tbl'

rename_dict_short_to_long = dict()
rename_dict_long_to_short = dict()
for each_isolate in open(name_correlation_txt):
    each_isolate_split = each_isolate.strip().split('\t')
    name_short = each_isolate_split[0]
    name_long = each_isolate_split[1]
    rename_dict_short_to_long[name_short] = name_long
    rename_dict_long_to_short[name_long] = name_short


only_id_txt  = '/Users/songweizhi/Desktop/Brady/rename/nifH_1copies_HKLL_itol_only_id.txt'
renamed_file = '/Users/songweizhi/Desktop/Brady/rename/nifH_1copies_HKLL_itol_only_id_renamed.txt'

only_id_txt  = '/Users/songweizhi/Desktop/Brady/rename/nifH_2copies_itol_only_id.txt'
renamed_file = '/Users/songweizhi/Desktop/Brady/rename/nifH_2copies_itol_only_id_renamed.txt'

renamed_file_handle = open(renamed_file, 'w')
for each in open(only_id_txt):
    each_split = each.strip().split(',')
    isolate_id = each_split[0]
    new_name = rename_dict_short_to_long.get(isolate_id, isolate_id)
    renamed_file_handle.write('%s,%s\n' % (new_name, ','.join(each_split[1:])))
renamed_file_handle.close()
