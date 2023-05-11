
japo_metadata_txt = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_qualified_no_redundant.txt'
id_to_name_itol_txt = '/Users/songweizhi/Desktop/Japonicum/get_gnm/id_to_name_iTOL.txt'

id_to_name_itol_txt_handle = open(id_to_name_itol_txt, 'w')
for each_japo in open(japo_metadata_txt):

    if not each_japo.startswith('ID\t'):
        each_japo_split = each_japo.strip().split('\t')
        gnm_id = each_japo_split[0]
        gnm_name = each_japo_split[1]

        to_write = ''
        if ('GCA_' in gnm_id):
            to_write = '%s\t%s__%s' % (gnm_id, gnm_id, gnm_name)
        elif ('HKCC' in gnm_id) or ('SZCC' in gnm_id):
            to_write = '%s\t%s' % (gnm_id, gnm_name)
        else:
            to_write = '%s\t%s__%s' % (gnm_id, gnm_id, gnm_name)
        id_to_name_itol_txt_handle.write(to_write + '\n')
id_to_name_itol_txt_handle.close()
