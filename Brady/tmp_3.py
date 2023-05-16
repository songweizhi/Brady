
gnm_35_txt = '/Users/songweizhi/Desktop/35.txt'
id2spe_txt                  = '/Users/songweizhi/Desktop/prep_metadata/id2spe.txt'

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


for each_gnm in open(gnm_35_txt):
    gnm_name = each_gnm.strip()
    strain = gnm_name.split('_')[-1]


    for each in species_to_id_dict:
        if strain in each:
            print('%s\t%s' % (gnm_name, each))

    print()