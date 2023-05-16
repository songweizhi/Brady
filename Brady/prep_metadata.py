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

# file in
combined_gbk_txt            = '/Users/songweizhi/Desktop/prep_metadata/combined_gbk.txt'
id2spe_txt                  = '/Users/songweizhi/Desktop/prep_metadata/id2spe.txt'
outgroup_gnm_txt            = '/Users/songweizhi/Desktop/prep_metadata/outgroup_gnms.txt'
duplicated_genomes_txt      = '/Users/songweizhi/Desktop/prep_metadata/duplicated_genomes.txt'
gnm_quality_txt             = '/Users/songweizhi/Desktop/prep_metadata/combined_quality_report.tsv'

# file out
metadata_txt                = '/Users/songweizhi/Desktop/prep_metadata/metadata.txt'

########################################################################################################################

# read in outgroup_gnm_txt
outgroup_gnm_set = set()
for each_outgroup in open(outgroup_gnm_txt):
    outgroup_gnm_set.add(each_outgroup.strip())

# read in outgroup_gnm_txt
gnm_stats_dict = dict()
for each_qual in open(gnm_quality_txt):
    each_qual_split = each_qual.strip().split('\t')
    gnm_name = each_qual_split[0]
    gnm_cpl  = each_qual_split[1]
    gnm_ctm  = each_qual_split[2]
    gnm_cd   = each_qual_split[5]
    gnm_n50  = each_qual_split[6]
    gnm_size = each_qual_split[8]
    gnm_gc   = each_qual_split[9]
    gnm_stats_dict[gnm_name] = [gnm_cpl, gnm_ctm, gnm_size, gnm_gc, gnm_cd, gnm_n50]

########## get the location of gbk file ##########

# read in combined_gbk_txt
gbk_name_to_all_location_dict = dict()
for each_gbk in open(combined_gbk_txt):
    pwd_gbk = each_gbk.strip()
    f_path, gbk_name, f_ext = sep_path_basename_ext(pwd_gbk)
    if gbk_name not in gbk_name_to_all_location_dict:
        gbk_name_to_all_location_dict[gbk_name] = set()
    gbk_name_to_all_location_dict[gbk_name].add(pwd_gbk)

# get genomes with multiple gbk file
gbk_name_to_location_dict = dict()
gbk_name_to_alternative_location_dict = dict()
for each_gbk in gbk_name_to_all_location_dict:

    file_location = gbk_name_to_all_location_dict[each_gbk]

    the_good_one = ''
    the_alternative_one = ''
    if len(file_location) == 1:
        the_good_one = [i for i in file_location][0]
    else:
        for each_location in file_location:
            if 'zong_Sishuo' in each_location:
                the_good_one = each_location
            else:
                the_alternative_one = each_location

    # add to dict
    gbk_name_to_location_dict[each_gbk] = the_good_one
    if the_alternative_one != '':
        gbk_name_to_alternative_location_dict[each_gbk] = the_alternative_one

########## get duplicated genomes ##########

# read in id2bspe_txt
id_to_species_dict = dict()
species_to_id_dict = dict()
for each_gnm in open(id2spe_txt):
    each_gnm_split = each_gnm.strip().split('\t')
    gnm_id = each_gnm_split[0]
    gbk_name = each_gnm_split[1]
    gbk_name_not_dot = gbk_name
    gbk_name_not_dot = gbk_name_not_dot.replace('._', '_').replace('.', '_')

    if ('_nif' in gnm_id) or ('_sym' in gnm_id) or ('_chromosome' in gnm_id) or ('_plasmid' in gnm_id):
        pass
    else:
        id_to_species_dict[gnm_id] = gbk_name
        species_to_id_dict[gbk_name] = gnm_id
        if '.' in gbk_name:
            species_to_id_dict[gbk_name_not_dot] = gnm_id


# read in provided duplicated genome first
duplicated_gnm_set = set()
gnm_to_duplicated_gnm_dict = dict()
for each_pair in open(duplicated_genomes_txt):
    each_pair_split = each_pair.strip().split('\t')
    duplicated_gnm_set.add(each_pair_split[1])
    gnm_to_duplicated_gnm_dict[each_pair_split[0]] = each_pair_split[1]


for each_gnm in gbk_name_to_location_dict:

    if ('_nif' in each_gnm) or ('_sym' in each_gnm) or ('_chromosome' in each_gnm) or ('_plasmid' in each_gnm):
        pass
    elif ('GCA_' in each_gnm):
        gca_gnm_name = id_to_species_dict.get(each_gnm, 'na')

        gca_gnm_name_not_dot = gca_gnm_name
        gca_gnm_name_not_dot = gca_gnm_name_not_dot.replace('._', '_').replace('.', '_')

        if gca_gnm_name in gbk_name_to_location_dict:
            duplicated_gnm_set.add(each_gnm)
            gnm_to_duplicated_gnm_dict[gca_gnm_name] = each_gnm

        elif gca_gnm_name_not_dot in gbk_name_to_location_dict:
            duplicated_gnm_set.add(each_gnm)
            gnm_to_duplicated_gnm_dict[gca_gnm_name_not_dot] = each_gnm

        else:
            pass  # these are the good ones


gbk_name_set_without_duplication = set()
for each_gbk in gbk_name_to_location_dict:
    if ('_nif' in each_gbk) or ('_sym' in each_gbk) or ('_chromosome' in each_gbk) or ('_plasmid' in each_gbk):
        pass
    else:
        if each_gbk not in duplicated_gnm_set:
            gbk_name_set_without_duplication.add(each_gbk)


metadata_txt_handle = open(metadata_txt, 'w')
metadata_txt_handle.write('ID\tName\tAlias\tCompleteness\tContamination\tSize\tGC\tCoding_density\tN50\tDuplicated\tOutgroup\tPath_to_gbk\tDuplicated_gbk\n')
for uniq_gbk in gbk_name_set_without_duplication:

    ##### column 1: get gnm id #####
    gnm_id = 'NA'
    if (uniq_gbk in outgroup_gnm_set):
        gnm_id = uniq_gbk
    elif 'GCA_' in uniq_gbk:
        gnm_id = uniq_gbk
    elif 'SZCC' in uniq_gbk:
        gnm_id = 'SZCC' + uniq_gbk.split('SZCC')[-1]
    elif 'HKCC' in uniq_gbk:
        gnm_id = 'HKCC' + uniq_gbk.split('HKCC')[-1]
    elif 'MSLE' in uniq_gbk:
        gnm_id = 'MSLE' + uniq_gbk.split('MSLE')[-1]
    elif uniq_gbk in species_to_id_dict:
        gnm_id = species_to_id_dict[uniq_gbk]
    else:
        # use the id of the duplicated one
        if uniq_gbk in gnm_to_duplicated_gnm_dict:
            duplicated_gnm = gnm_to_duplicated_gnm_dict[uniq_gbk]
            gnm_id = species_to_id_dict.get(duplicated_gnm, 'NA')
        else:
            gnm_id = uniq_gbk

    ##### column 2: get gnm name #####
    genome_name = uniq_gbk
    if 'GCA_' in genome_name:
        genome_name = id_to_species_dict.get(genome_name, genome_name)

    ##### column 3: get alias #####

    alias = 'na'
    is_duplicated = 'no'

    if uniq_gbk in gnm_to_duplicated_gnm_dict:
        alias = gnm_to_duplicated_gnm_dict[uniq_gbk]
        is_duplicated = 'yes'

    if 'GCA_' in alias:
        alias = id_to_species_dict.get(alias, alias)

        if 'GCA_' not in gnm_id:
            gnm_id = species_to_id_dict.get(alias, alias)

    ##### column 4: get out group #####
    is_outgroup = 'no'
    if uniq_gbk in outgroup_gnm_set:
        is_outgroup = 'yes'

    ##### column get path to gbk #####

    path_to_gbk = gbk_name_to_location_dict.get(uniq_gbk, 'missing')

    path_to_duplicated_gbk = 'na'
    if is_duplicated == 'yes':

        path_to_duplicated_gbk = gbk_name_to_all_location_dict.get(gnm_id, 'na')
        if len(path_to_duplicated_gbk) == 1:
            path_to_duplicated_gbk = [i for i in path_to_duplicated_gbk][0]

    ##### get genome stats #####

    gnm_stats_list = ['na', 'na', 'na', 'na', 'na', 'na']
    if gnm_id in gnm_stats_dict:
        gnm_stats_list = gnm_stats_dict[gnm_id]
    elif genome_name in gnm_stats_dict:
        gnm_stats_list = gnm_stats_dict[genome_name]
    elif alias in gnm_stats_dict:
        gnm_stats_list = gnm_stats_dict[alias]

    ##### write out #####
    if genome_name == alias:
        alias = 'na'
    str_to_write = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s' % (gnm_id, genome_name, alias, '\t'.join(gnm_stats_list), is_duplicated, is_outgroup, path_to_gbk, path_to_duplicated_gbk)
    metadata_txt_handle.write(str_to_write + '\n')

metadata_txt_handle.close()
