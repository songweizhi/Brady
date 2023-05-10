import os


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


########################################################################################################################

'''
columns in metadata

Genome_ID	        
Genome_name	        
Genome_name_alias   
Genome_size
Completeness	
Contamination	
Taxonomy            
Bradyrhizobium      # yes/no
Clade               # Basal_1/Basal_2/Basal_3/non_basal
Location	        # full path to gbk file
GBK_by	            # id/name/alias
nif_copy_number     # 1/2
Description         # other descriptions
assemly_status      # succeed/failed
outgroup            # yes/no
'''

# file in
japonicum_jj_txt            = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_JJ.txt'
japonicum_ll_zzc_lj_txt     = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_LL_ZZC_LJ.txt'
japonicum_outgroup_txt      = '/Users/songweizhi/Desktop/Japonicum/get_gnm/japonicum_outgroup.txt'

gbk_dir_jj729_download      = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_jj729_download'
gbk_dir_jj731_sequencing    = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_jj731_sequencing'
gbk_dir_ling                = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_ling'
gbk_dir_jieliu              = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_jieliu'
gbk_dir_zong                = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_zong'
gbk_dir_zzc                 = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_zzc'

id2bspe_txt                 = '/Users/songweizhi/Desktop/Japonicum/get_gnm/id2bspe.txt'
specified_gnm_to_id_txt     = '/Users/songweizhi/Desktop/Japonicum/get_gnm/specified_id_to_gnm.txt'
genome_quality_txt          = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_genome_quality_reformatted.txt'
cpl_cutoff                  = 95
ctm_cutoff                  = 5
copy_file                   = True

# file out
japo_metadata_txt_all                    = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_all.txt'
japo_metadata_txt_qualified              = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_qualified.txt'
japo_metadata_txt_qualified_no_redundant = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_qualified_no_redundant.txt'
qualified_gbk_dir_no_redundant           = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_qualified_nonredundant'

'''
missing: Bradyrhizobium_liaoningense_HKCCYLR1011    # 要不是组装问题，要不就是污染和完整度的问题
missing: Bradyrhizobium_oligotrophicum_HKCCYLS3021  # 要不是组装问题，要不就是污染和完整度的问题
missing: Bradyrhizobium_AXAI_s_HKCCYLS3084          # 要不是组装问题，要不就是污染和完整度的问题
'''

########################################################################################################################

gnm_cpl_dict = dict()
gnm_ctm_dict = dict()
for each_q in open(genome_quality_txt):
    if not each_q.startswith('Genome	Completeness	Contamination	Heterogeneity'):
        each_q_split = each_q.strip().split('\t')
        gnm_id = each_q_split[0]
        gnm_cpl = float(each_q_split[1])
        gnm_ctm = float(each_q_split[2])
        gnm_cpl_dict[gnm_id] = gnm_cpl
        gnm_ctm_dict[gnm_id] = gnm_ctm

if copy_file is True:
    if os.path.isdir(qualified_gbk_dir_no_redundant) is True:
        os.system('rm -r %s' % qualified_gbk_dir_no_redundant)
    os.system('mkdir %s' % qualified_gbk_dir_no_redundant)

specified_gnm_to_id_dict = dict()
for each_gnm in open(specified_gnm_to_id_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_name = each_gnm_split[1]
    specified_gnm_to_id_dict[gnm_name] = gnm_id

id_to_species_dict = dict()
species_to_id_dict = dict()
for each_gnm in open(id2bspe_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gnm_name = each_gnm_split[1]
    if ('_nif' in gnm_id) or ('_sym' in gnm_id) or ('_chromosome' in gnm_id) or ('_plasmid' in gnm_id):
        pass
    else:
        id_to_species_dict[gnm_id] = gnm_name
        species_to_id_dict[gnm_name] = gnm_id

japonicum_set_jj = set()
for each_gnm in open(japonicum_jj_txt):
    japonicum_set_jj.add(each_gnm.strip())

japonicum_set_ll_zzc_lj = set()
for each_gnm in open(japonicum_ll_zzc_lj_txt):
    japonicum_set_ll_zzc_lj.add(each_gnm.strip())

japonicum_set_outgroup = set()
for each_gnm in open(japonicum_outgroup_txt):
    japonicum_set_outgroup.add(each_gnm.strip())

japonicum_set_all                       = japonicum_set_jj.union(japonicum_set_ll_zzc_lj)
japonicum_set_all_with_outgroup         = japonicum_set_all.union(japonicum_set_outgroup)
japonicum_set_all_sorted                = sorted([i for i in japonicum_set_all])
japonicum_set_all_with_outgroup_sorted  = sorted([i for i in japonicum_set_all_with_outgroup])


gbk_by_id_set = set()
gnm_to_id_dict = dict()
found_gnm_set = set()
missing_gnm_set = set()
found_in_dict= dict()
for each_japo in japonicum_set_all_with_outgroup_sorted:

    if each_japo not in found_in_dict:
        found_in_dict[each_japo] = set()

    source_gbk_jj729_download   = '%s/%s.gbk' % (gbk_dir_jj729_download, each_japo)
    source_gbk_jj731_sequencing = '%s/%s.gbk' % (gbk_dir_jj731_sequencing, each_japo)
    source_gbk_ling             = '%s/%s.gbk' % (gbk_dir_ling, each_japo)
    source_gbk_zong             = '%s/%s.gbk' % (gbk_dir_zong, each_japo)
    source_gbk_jieliu           = '%s/%s.gbk' % (gbk_dir_jieliu, each_japo)
    source_gbk_zzc              = '%s/%s.gbk' % (gbk_dir_zzc, each_japo)

    if os.path.isfile(source_gbk_zzc) is True:
        found_in_dict[each_japo].add('zzc')
    if os.path.isfile(source_gbk_jj729_download) is True:
        found_in_dict[each_japo].add('jj729_download')
    if os.path.isfile(source_gbk_jj731_sequencing) is True:
        found_in_dict[each_japo].add('jj731_sequencing')
    if os.path.isfile(source_gbk_ling) is True:
        found_in_dict[each_japo].add('ling')
    if os.path.isfile(source_gbk_zong) is True:
        found_in_dict[each_japo].add('zong')
    if os.path.isfile(source_gbk_jieliu) is True:
        found_in_dict[each_japo].add('jieliu')

    if (os.path.isfile(source_gbk_zzc) is True) or (os.path.isfile(source_gbk_jj729_download) is True) or (os.path.isfile(source_gbk_jj731_sequencing) is True) or (os.path.isfile(source_gbk_ling) is True) or (os.path.isfile(source_gbk_zong) is True) or (os.path.isfile(source_gbk_jieliu) is True):
        found_gnm_set.add(each_japo)
    else:
        if each_japo in species_to_id_dict:
            japo_by_id      = species_to_id_dict[each_japo]
            gbk_by_id_jj729 = '%s/%s.gbk' % (gbk_dir_jj729_download, japo_by_id)
            if (os.path.isfile(gbk_by_id_jj729) is True):
                found_gnm_set.add(each_japo)
                found_in_dict[each_japo].add('jj729_download')
                gbk_by_id_set.add(each_japo)
            else:
                missing_gnm_set.add(each_japo)
        else:
            missing_gnm_set.add(each_japo)

qualified_gnm_num = 0
japo_metadata_txt_all_handle = open(japo_metadata_txt_all, 'w')
japo_metadata_txt_qualified_handle = open(japo_metadata_txt_qualified, 'w')
japo_metadata_txt_all_handle.write('ID\tGenome\tCompleteness\tContamination\tLocation\tOutgroup\tGBK_by\n')
japo_metadata_txt_qualified_handle.write('ID\tGenome\tCompleteness\tContamination\tLocation\tOutgroup\tGBK_by\n')
for japo_gnm in found_in_dict:

    # get is_outgroup info
    is_outgroup = 'n'
    if japo_gnm in japonicum_set_outgroup:
        is_outgroup = 'y'

    # get japo_id info
    japo_id = ''
    if japo_gnm in specified_gnm_to_id_dict:
        japo_id = specified_gnm_to_id_dict[japo_gnm]
    elif japo_gnm in species_to_id_dict:
        japo_id = species_to_id_dict[japo_gnm]
    elif 'HKCC' in japo_gnm:
        japo_id = 'HKCC' + japo_gnm.split('HKCC')[-1]
    elif 'SZCC' in japo_gnm:
        japo_id = 'SZCC' + japo_gnm.split('SZCC')[-1]
    else:
        print(japo_gnm)

    # get from_dir info
    from_dir = 'missing'
    dir_list_sorted = sorted([i for i in found_in_dict[japo_gnm]])
    if len(dir_list_sorted) == 1:
        from_dir = dir_list_sorted[0]
    elif len(dir_list_sorted) > 1:
        if dir_list_sorted == ['ling', 'zong']:
            from_dir = 'zong'
        elif dir_list_sorted == ['jieliu', 'jj731_sequencing']:
            from_dir = 'jj731_sequencing'
        else:
            print('missing: %s\t%s' % (japo_gnm, dir_list_sorted))
    else:
        print('missing:\t%s' % japo_gnm)

    # gbk by id
    gbk_by = 'name'
    if japo_gnm in gbk_by_id_set:
        gbk_by = 'id'

    # get quality
    gbk_cpl = gnm_cpl_dict.get(japo_id, 'na')
    gbk_ctm = gnm_ctm_dict.get(japo_id, 'na')

    japo_metadata_txt_all_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (japo_id, japo_gnm, gbk_cpl, gbk_ctm, from_dir, is_outgroup, gbk_by))

    if (gbk_cpl != 'na') and (gbk_ctm != 'na'):
        if (from_dir != 'missing') and (gbk_cpl >= 95) and (gbk_ctm <= 5):
            japo_metadata_txt_qualified_handle.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (japo_id, japo_gnm, gbk_cpl, gbk_ctm, from_dir, is_outgroup, gbk_by))
            qualified_gnm_num += 1

japo_metadata_txt_all_handle.close()
japo_metadata_txt_qualified_handle.close()

# get gnm_name_to_id_dict
gnm_name_to_id_dict= dict()
for gnm_meta in open(japo_metadata_txt_qualified):
    if not gnm_meta.startswith('ID	Genome'):
        gnm_meta_split = gnm_meta.strip().split('\t')
        gnm_id       = gnm_meta_split[0]
        gnm_name     = gnm_meta_split[1]
        gnm_name_to_id_dict[gnm_name] = gnm_id

# get duplicated genome list
duplicated_gnm_set = set()
dot_to_no_dot_dict = dict()
for each_gnm in gnm_name_to_id_dict:
    if '.' in each_gnm:
        gnm_name_not_dot = each_gnm
        gnm_name_not_dot = gnm_name_not_dot.replace('._', '_').replace('.', '_')
        if gnm_name_not_dot in gnm_name_to_id_dict:
            duplicated_gnm_set.add(each_gnm)
            duplicated_gnm_set.add(gnm_name_not_dot)
            dot_to_no_dot_dict[each_gnm] = gnm_name_not_dot


# get no redundant qualified gnms
no_redundant_gnm_num = 0
japo_metadata_txt_qualified_no_redundant_handle = open(japo_metadata_txt_qualified_no_redundant, 'w')
japo_metadata_txt_qualified_no_redundant_handle.write('ID\tGenome\tCompleteness\tContamination\tLocation\tOutgroup\tGBK_by\tOther_name\n')
for gnm_meta in open(japo_metadata_txt_qualified):
    if not gnm_meta.startswith('ID	Genome'):
        gnm_meta_split = gnm_meta.strip().split('\t')
        gnm_id       = gnm_meta_split[0]
        gnm_name     = gnm_meta_split[1]
        gbk_dir      = gnm_meta_split[4]
        gbk_named_by = gnm_meta_split[6]
        if gnm_name not in duplicated_gnm_set:
            japo_metadata_txt_qualified_no_redundant_handle.write('%s\tno\n' % gnm_meta.strip())
            no_redundant_gnm_num += 1
        else:
            if '.' in gnm_name:
                name_without_dot = dot_to_no_dot_dict[gnm_name]
                japo_metadata_txt_qualified_no_redundant_handle.write('%s\t%s\n' % (gnm_meta.strip().replace(gnm_name, name_without_dot), gnm_name))
                no_redundant_gnm_num += 1
japo_metadata_txt_qualified_no_redundant_handle.close()


print('all_japonicum\t%s'               % len(japonicum_set_all))
print('all_japonicum_with_outgroup\t%s' % len(japonicum_set_all_with_outgroup))
print('found genome\t%s'                % len(found_gnm_set))
print('missing genome\t%s'              % len(missing_gnm_set))
print('qualified genome\t%s'            % qualified_gnm_num)
print('nonredundant genome\t%s'         % no_redundant_gnm_num)

# get gbk of qualified gnms
print('copying %s qualified gbk files to %s' % (no_redundant_gnm_num, qualified_gbk_dir_no_redundant))
for gnm_meta in open(japo_metadata_txt_qualified_no_redundant):
    if not gnm_meta.startswith('ID	Genome'):
        gnm_meta_split = gnm_meta.strip().split('\t')
        gnm_id       = gnm_meta_split[0]
        gnm_name     = gnm_meta_split[1]
        gbk_dir      = gnm_meta_split[4]
        gbk_named_by = gnm_meta_split[6]
        pwd_gbk_dir  = '/Users/songweizhi/Desktop/Japonicum/get_gnm/gbk_%s' % gbk_dir

        pwd_gbk_file = ''
        if gbk_named_by == 'id':
            pwd_gbk_file = '%s/%s.gbk' % (pwd_gbk_dir, gnm_id)
        if gbk_named_by == 'name':
            pwd_gbk_file = '%s/%s.gbk' % (pwd_gbk_dir, gnm_name)

        if os.path.isfile(pwd_gbk_file) is False:
            print(gnm_meta.strip())
        else:
            pwd_gbk_file_copied = '%s/%s.gbk' % (qualified_gbk_dir_no_redundant, gnm_id)
            if copy_file is True:
                os.system('cp %s %s' % (pwd_gbk_file, pwd_gbk_file_copied))


