import os
import glob


KEGG_DB_seq2ko              = '/Users/songweizhi/DB/KEGG_v2017_luo/kegg_db_seq2ko.txt'
lipid_metabolism_pwys_txt   = '/Users/songweizhi/Desktop/Brady/Lipid_metabolism_pwys.txt'
best_hit_dir_394_PB_KEGG    = '/Users/songweizhi/Desktop/Brady/394_PB_KEGG_best_hit_dir'
bacterial_ko_txt            = '/Users/songweizhi/DB/KEGG_v2017_luo/bacterial_KO.txt'
e_index_cutoff              = 30

op_df                       = '/Users/songweizhi/Desktop/Brady/394_PB_KEGG_lipid_pwy.txt'


bacterial_ko_set = set()
for each_ko in open(bacterial_ko_txt):
    bacterial_ko_set.add(each_ko.strip())

# get db_seq_to_KO_dict
db_seq_to_KO_dict = {}
for each_hit in open(KEGG_DB_seq2ko):
    each_hit_split = each_hit.strip().split('\t')
    db_seq = each_hit_split[0]
    hit_id_KO = each_hit_split[1]
    if hit_id_KO != '':
        db_seq_to_KO_dict[db_seq] = hit_id_KO

Cs_description_dict = {}
Ds_description_dict = {}
D2ABCD_dict = {}
current_B = ''
current_C = ''
for each_line in open(lipid_metabolism_pwys_txt):
    if 'fungi type' in each_line:
        print(each_line.strip())

    if each_line[0] in ['B', 'C', 'D']:
        each_line_split = each_line.strip().split(' ')

        if each_line[0] == 'B':
            if len(each_line_split) > 1:
                current_B_id = each_line_split[2]
                current_B_description = ' '.join(each_line_split[3:])
                current_B = current_B_id

        elif each_line[0] == 'C':
            current_C_id = each_line_split[4]
            current_C_description = ' '.join(each_line_split[5:])
            current_C = current_C_id
            Cs_description_dict[current_C_id] = current_C_description

        elif each_line[0] == 'D':
            current_D_id = each_line_split[6]
            current_D_description = ' '.join(each_line_split[7:])
            Ds_description_dict[current_D_id] = current_D_description
            ABCD_value = 'B_%s|C_%s' % (current_B, current_C)
            if current_D_id not in D2ABCD_dict:
                D2ABCD_dict[current_D_id] = [ABCD_value]
            elif (current_D_id in D2ABCD_dict) and (ABCD_value not in D2ABCD_dict[current_D_id]):
                D2ABCD_dict[current_D_id].append(ABCD_value)

ko_to_path_dict = dict()
path_to_ko_dict = dict()
for each_d in D2ABCD_dict:
    if each_d in bacterial_ko_set:
        d_group_list = D2ABCD_dict[each_d]
        for each_group in d_group_list:
            group_split = each_group.split('|')
            c_id = group_split[1]
            if c_id not in path_to_ko_dict:
                path_to_ko_dict[c_id] = {each_d}
            else:
                path_to_ko_dict[c_id].add(each_d)

for each_path in path_to_ko_dict:
    path_gene_set = path_to_ko_dict[each_path]
    path_desc = Cs_description_dict.get(each_path[2:], '')
    #print('%s\t%s (%s)\t%s' % (each_path, path_desc, len(path_gene_set), path_gene_set))

path_list_sorted = sorted([i[2:] for i in path_to_ko_dict])
path_list_sorted_with_desc = []
for i in path_list_sorted:
    id_with_desc = '%s_%s(%s)' % (i, Cs_description_dict[i], len(path_to_ko_dict['C_'+ i]))
    path_list_sorted_with_desc.append(id_with_desc)

best_hit_file_re = '%s/*.tab' % best_hit_dir_394_PB_KEGG
best_hit_file_list = [os.path.basename(file_name) for file_name in glob.glob(best_hit_file_re)]

op_df_handle = open(op_df, 'w')
op_df_handle.write('Genome\t%s\n' % '\t'.join(path_list_sorted_with_desc))
for each_file in best_hit_file_list:

    pwd_best_hit = '%s/%s' % (best_hit_dir_394_PB_KEGG, each_file)
    gnm_id = each_file.split('_blast_best_hits.tab')[0]

    current_gnm_anno_dict = dict()
    for each_c_id in Cs_description_dict:
        current_gnm_anno_dict[each_c_id] = set()

    for each_line in open(pwd_best_hit):
        each_file_split = each_line.strip().split('\t')
        ref_id = each_file_split[1]
        evalue = each_file_split[10]

        # only look at good hit
        good_hit = False
        if evalue == '0.0':
            good_hit = True
        elif evalue == '0.0e+00':
            good_hit = True
        else:
            if '-' in evalue:
                e_index = int(evalue.split('-')[-1])
                if e_index >= e_index_cutoff:
                    good_hit = True
            else:
                pass
                #print(each_file_split)

        if good_hit is True:
            ref_ko = db_seq_to_KO_dict.get(ref_id, 'NA')
            if ref_ko in D2ABCD_dict:
                ko_pwy_list = D2ABCD_dict[ref_ko]

                for each_pwy in ko_pwy_list:
                    pwy_id = each_pwy.split('|')[1][2:]
                    current_gnm_anno_dict[pwy_id].add(ref_ko)

    anno_result_list = [gnm_id]
    for each_pwy in path_list_sorted:
        full_ko_list = path_to_ko_dict['C_' + each_pwy]
        identified_ko_list = current_gnm_anno_dict[each_pwy]
        identified_ko_pct = float("{0:.2f}".format(len(identified_ko_list)*100/len(full_ko_list)))
        anno_result_list.append(str(identified_ko_pct))

    op_df_handle.write('\t'.join(anno_result_list) + '\n')



print(path_list_sorted)
print(len(path_list_sorted))