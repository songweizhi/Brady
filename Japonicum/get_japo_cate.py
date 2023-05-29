
nif_pa_txt              = '/Users/songweizhi/Desktop/Japonicum/Japonicum_nif_PA_df.txt'
nod_pa_txt              = '/Users/songweizhi/Desktop/Japonicum/Japonicum_nod_PA.txt'
genome_cate_txt         = '/Users/songweizhi/Desktop/Japonicum/genome_cate.txt'

genome_cate_txt_by_nif  = '/Users/songweizhi/Desktop/Japonicum/genome_cate_by_nif.txt'
genome_cate_txt_by_nod  = '/Users/songweizhi/Desktop/Japonicum/genome_cate_by_nod.txt'


gnm_set = set()
nif_pa_dict = dict()
for each_gnm in open(nif_pa_txt):
    if not each_gnm.startswith('Genome'):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        pa_list = [int(i) for i in each_gnm_split[1:]]
        present_pct = sum(pa_list)/len(pa_list)
        if present_pct >= 0.75:
            nif_pa_dict[gnm_id] = 1
        else:
            nif_pa_dict[gnm_id] = 0

        gnm_set.add(gnm_id)


nod_pa_dict = dict()
for each_gnm in open(nod_pa_txt):
    if not each_gnm.startswith('Genome'):
        each_gnm_split = each_gnm.strip().split('\t')
        gnm_id = each_gnm_split[0]
        pa_list = [int(i) for i in each_gnm_split[1:]]
        present_pct = sum(pa_list)/len(pa_list)
        if present_pct == 1:
            nod_pa_dict[gnm_id] = 1
        elif present_pct == 0:
            nod_pa_dict[gnm_id] = 0
        else:
            if sum(pa_list[:3]) == 1:
                nod_pa_dict[gnm_id] = 0
            else:
                if present_pct >= 0.75:
                    nod_pa_dict[gnm_id] = 1
                else:
                    nod_pa_dict[gnm_id] = 0

        gnm_set.add(gnm_id)



genome_cate_txt_by_nif_handle = open(genome_cate_txt_by_nif, 'w')
genome_cate_txt_by_nod_handle = open(genome_cate_txt_by_nod, 'w')
genome_cate_txt_by_nif_handle.write(',nif\n')
genome_cate_txt_by_nod_handle.write(',nod\n')
gnm_cate_dict = dict()
for each_gnm in gnm_set:

    if each_gnm != 'GCA_024171065.1':

        nif_pa = nif_pa_dict[each_gnm]
        nod_pa = nod_pa_dict[each_gnm]

        if nif_pa == 1:
            genome_cate_txt_by_nif_handle.write('%s,1\n' % each_gnm)
        else:
            genome_cate_txt_by_nif_handle.write('%s,0\n' % each_gnm)

        if nod_pa == 1:
            genome_cate_txt_by_nod_handle.write('%s,1\n' % each_gnm)
        else:
            genome_cate_txt_by_nod_handle.write('%s,0\n' % each_gnm)




        if (nif_pa == 1) and (nod_pa == 1):
            gnm_cate_dict[each_gnm] = 'BothFound'
        elif (nif_pa == 0) and (nod_pa == 0):
            gnm_cate_dict[each_gnm] = 'BothAbsent'
        elif (nif_pa == 1) and (nod_pa == 0):
            gnm_cate_dict[each_gnm] = 'OnlyNif'
        else:
            print('%s\t%s\t%s' % (each_gnm, nif_pa, nod_pa))

genome_cate_txt_by_nif_handle.close()
genome_cate_txt_by_nod_handle.close()


genome_cate_txt_handle = open(genome_cate_txt, 'w')
for each_gnm in gnm_cate_dict:
    gnm_cate = gnm_cate_dict[each_gnm]
    genome_cate_txt_handle.write('%s\t%s\n' % (each_gnm, gnm_cate))
genome_cate_txt_handle.close()



'''
GCA_024171065.1 Bradyrhizobium_japonicum_USDA_5

# nif
GCA_024171065.1	1	0	0	0	0	0	0	0

# nod
GCA_024171065.1	1	1	1	1	1

'''


