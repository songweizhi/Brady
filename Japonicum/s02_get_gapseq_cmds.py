import os
import pandas as pd


def divide_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


'''
all			Pathways|seed|kegg

amino		Amino-Acid-Biosynthesis	
nucl		Nucleotide-Biosynthesis
cofactor	Cofactor-Biosynthesis
carbo		CARBO-BIOSYNTHESIS
carbo-deg	Carbohydrates-Degradation
polyamine	Polyamine-Biosynthesis
fatty		Fatty-acid-biosynthesis
energy		Energy-Metabolism
terpenoid	Terpenoid-Biosynthesis
degradation	Degradation

'''

pathway_key_list = ['amino', 'nucl', 'cofactor', 'carbo', 'carbo-deg', 'polyamine', 'fatty', 'energy', 'terpenoid', 'degradation']
#pathway_key_list = ['carbo-deg']

# file in
japo_metadata_txt   = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_qualified_no_redundant.txt'
prot_dir            = '/home-user/wzsong/Japonicum/data/protein'
gapseq_exe          = '/home-user/wzsong/Software/gapseq/gapseq'
gapseq_tmp_dir      = '/home-user/wzsong/Japonicum/gapseq_metacyc/gapseq_tmp_dir'
pathway_db          = 'metacyc'
js_cpu              = 1
cmd_per_js          = 50
submit_js           = False
js_dir              = '/Users/songweizhi/Desktop/Japonicum_gapseq_js'


if os.path.isdir(js_dir) is True:
    os.system('rm -r %s' % js_dir)
os.mkdir(js_dir)


print('conda activate gapseq-dev')
for pathway_key in pathway_key_list:

    print('mkdir %s' % pathway_key)

    submit_cmds_txt = '%s/gapseq_%s_%s_cmds_%s.txt' % (js_dir, pathway_db, pathway_key, cmd_per_js)

    # get gapseq_cmd_list
    gapseq_cmd_list = []
    for gnm_meta in open(japo_metadata_txt):
        if not gnm_meta.startswith('ID\t'):
            gnm_meta_split = gnm_meta.strip().split('\t')
            if gnm_meta_split[2] != 'missing':
                gnm_id     = gnm_meta_split[0]
                pwd_prot   = '%s/%s.protein' % (prot_dir, gnm_id)
                gapseq_cmd = '%s find -p %s -t Bacteria -l %s -K %s -M prot -T %s %s > %s_%s_%s_stdout.txt' % (gapseq_exe, pathway_key, pathway_db, js_cpu, gapseq_tmp_dir, pwd_prot, gnm_id, pathway_db, pathway_key)
                gapseq_cmd_list.append(gapseq_cmd)

    # write out gapseq cmds
    submit_cmds_txt_handle = open(submit_cmds_txt, 'w')
    js_index =1
    for each_sub_list in list(divide_list(gapseq_cmd_list, cmd_per_js)):
        each_sub_list_as_str = ';'.join(each_sub_list)
        submit_cmd = 'submitHPC.sh --cmd "%s" -n %s -c %s_%s_%s' % (each_sub_list_as_str, js_cpu, pathway_db, pathway_key, js_index)
        submit_cmds_txt_handle.write(submit_cmd + '\n')

        if submit_js is True:
            print(submit_cmd)
            #os.system(submit_cmd)

        js_index += 1
    submit_cmds_txt_handle.close()
