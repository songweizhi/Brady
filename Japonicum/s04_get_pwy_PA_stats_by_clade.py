
'''
cd /Users/songweizhi/Desktop/Japonicum
BioSAK iTOL -ColorStrip -lg Japo_analysis_clades.txt -lt AnalysisClade -o Japo_analysis_clades_iTOL_Strip.txt
'''


analysis_clades_txt = '/Users/songweizhi/Desktop/Japonicum/Japo_analysis_clades.txt'
tree_file           = '/Users/songweizhi/Desktop/Japonicum/Japonicum_121_OG_tree_LG.treefile'


analysis_clade_to_gnm_dict = dict()
for each_japo in open(analysis_clades_txt):
    each_japo_split = each_japo.strip().split('\t')
    gnm_id = each_japo_split[0]
    ac_id = each_japo_split[1]
    if ac_id not in analysis_clade_to_gnm_dict:
        analysis_clade_to_gnm_dict[ac_id] = set()
    analysis_clade_to_gnm_dict[ac_id].add(gnm_id)
print(analysis_clade_to_gnm_dict)


for each_ac in analysis_clade_to_gnm_dict:
    ac_gnm_set = analysis_clade_to_gnm_dict[each_ac]
    print('%s\t%s\t%s' % (each_ac, len(ac_gnm_set), ac_gnm_set))





