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


def transpose_csv(file_in, file_out, sep_symbol, column_name_pos, row_name_pos):

    csv = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    df_csv = pd.DataFrame(data=csv)
    transposed_csv = df_csv.T
    transposed_csv.to_csv(file_out, sep=sep_symbol)


def subset_df_by_cols(file_in, file_out, to_keep_col_list, sep_symbol, column_name_pos, row_name_pos):

    df = pd.read_csv(file_in, sep=sep_symbol, header=column_name_pos, index_col=row_name_pos)
    subset_df = df.loc[:, to_keep_col_list]
    subset_df.to_csv(file_out, sep=sep_symbol)


########################################################################################################################

key_list                                    = ['amino', 'carbo', 'carbo-deg', 'cofactor', 'degradation', 'energy', 'fatty', 'nucl', 'polyamine', 'terpenoid']
pathway_pa_dir                              = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc'
dgpa_wd                                     = '/Users/songweizhi/Desktop/Japonicum/differential_gene_presence_analysis_with_name'
gnm_cate_txt                                = '/Users/songweizhi/Desktop/Japonicum/differential_gene_presence_analysis_with_name/genome_cate_by_nod.txt'
pwy_cate_color_txt                          = '/Users/songweizhi/Desktop/Japonicum/gapseq_metacyc/pwy_cate_color.txt'

'''

cd /Users/songweizhi/Desktop/Japonicum/differential_gene_presence_analysis
scoary -g Pathways_PA_energy_transposed_for_Scoary.txt -t genome_cate_by_nod.txt -o Scoary_op_dir

'''

########################################################################################################################

fixed_col_header = 'Gene,Non-unique Gene name,Annotation,No. isolates,No. sequences,Avg sequences per isolate,Genome fragment,Order within fragment,Accessory Fragment,Accessory Order with Fragment,QC,Min group size nuc,Max group size nuc,Avg group size nuc'

gnm_cate_f_path, gnm_cate_f_base, gnm_cate_f_ext = sep_path_basename_ext(gnm_cate_txt)

# read in color
pwy_cate_color_dict = dict()
if os.path.isfile(pwy_cate_color_txt) is True:
    for each_pwy_cate in open(pwy_cate_color_txt):
        each_pwy_cate_split = each_pwy_cate.strip().split('\t')
        pwy_cate_color_dict[each_pwy_cate_split[0]] = each_pwy_cate_split[1]


for each_key in key_list:

    pwd_pathway_pa_txt                  = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10.txt'                % (pathway_pa_dir, each_key)
    pwd_pathway_pa_txt_t                = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_T.txt'              % (dgpa_wd, each_key)
    pwd_pathway_pa_txt_t_Scoary         = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_T_for_Scoary.txt'   % (dgpa_wd, each_key)
    pwd_scoary_op_dir                   = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_scoary_op_%s'       % (dgpa_wd, each_key, gnm_cate_f_base)
    pwd_pathway_pa_txt_diff             = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_diff.txt'           % (dgpa_wd, each_key)
    pwd_pathway_pa_txt_diff_0as1        = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_diff_0as1.txt'      % (dgpa_wd, each_key)
    pwd_pathway_pa_txt_diff_0as1_itol   = '%s/Pathways_PA_%s_no_single_value_cols_min1_num10_diff_0as1_iTOL.txt' % (dgpa_wd, each_key)

    # get color
    current_color                   = pwy_cate_color_dict.get(each_key, 'lightblue')

    # transpose data matrix
    transpose_csv(pwd_pathway_pa_txt, pwd_pathway_pa_txt_t, '\t', 0, 0)

    # change to Scoary format
    pwd_pathway_pa_txt_t_Scoary_handle = open(pwd_pathway_pa_txt_t_Scoary, 'w')
    for each_line in open(pwd_pathway_pa_txt_t):
        each_line_split = each_line.strip().split('\t')
        if each_line.startswith('\t'):
            pwd_pathway_pa_txt_t_Scoary_handle.write('%s,%s\n' % (fixed_col_header, ','.join(each_line_split)))
        else:
            pwd_pathway_pa_txt_t_Scoary_handle.write('%s,,pathway,100,100,1,1,1,,,,1000,1000,1000,%s\n' % (each_line_split[0], ','.join(each_line_split[1:])))
    pwd_pathway_pa_txt_t_Scoary_handle.close()

    # run Scoary
    scoary_cmd = 'scoary -t %s -g %s -o %s' % (gnm_cate_txt, pwd_pathway_pa_txt_t_Scoary, pwd_scoary_op_dir)
    print(scoary_cmd)
    #os.system(scoary_cmd)

    # get Scoary op txt
    pwd_scoary_op_txt_re = '%s/*.results.csv' % pwd_scoary_op_dir
    pwd_scoary_op_txt   = glob.glob(pwd_scoary_op_txt_re)[0]

    print('pwd_scoary_op_txt')
    print(pwd_scoary_op_txt)

    # get genes with differential presence
    diff_pwy_id_set = set()
    for each_line in open(pwd_scoary_op_txt):
        if not each_line.startswith('"Gene",'):
            each_line_split = each_line.strip().split(',')
            pwy_id          = each_line_split[0].replace('"', '')
            Benjamini_H_p   = float(each_line_split[12].replace('"', ''))
            if Benjamini_H_p <= 0.05:
                diff_pwy_id_set.add(pwy_id)

    print('diff_pwy_id_set')
    print(diff_pwy_id_set)

    # Subset PA df by selected columns
    subset_df_by_cols(pwd_pathway_pa_txt, pwd_pathway_pa_txt_diff, diff_pwy_id_set, '\t', 0, 0)

    # turn 0 to -1
    pwd_pathway_pa_txt_diff_0as1_handle = open(pwd_pathway_pa_txt_diff_0as1, 'w')
    for each_line in open(pwd_pathway_pa_txt_diff):
        pwd_pathway_pa_txt_diff_0as1_handle.write(each_line.replace('\t0', '\t-1'))
    pwd_pathway_pa_txt_diff_0as1_handle.close()

    # iTOL
    itol_cmd = 'BioSAK iTOL -Binary -lm %s -lt %s_diff -gc "%s" -o %s' % (pwd_pathway_pa_txt_diff_0as1, each_key, current_color, pwd_pathway_pa_txt_diff_0as1_itol)
    print(itol_cmd)
    os.system(itol_cmd)

