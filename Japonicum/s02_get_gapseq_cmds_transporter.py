import os
import pandas as pd


def divide_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


# file in
japo_metadata_txt   = '/Users/songweizhi/Desktop/Japonicum/get_gnm/Japonicum_metadata_qualified_no_redundant.txt'
prot_dir            = '/home-user/wzsong/Japonicum/data/protein'
gapseq_exe          = '/home-user/wzsong/Software/gapseq/gapseq'
gapseq_tmp_dir      = '/home-user/wzsong/Japonicum/gapseq_transporter'
js_cpu              = 1
cmd_per_js          = 25
submit_js           = False
js_dir              = '/Users/songweizhi/Desktop/Japonicum_gapseq_transport_js'


if os.path.isdir(js_dir) is True:
    os.system('rm -r %s' % js_dir)
os.mkdir(js_dir)


print('conda activate gapseq-dev')
print('cd /home-user/wzsong/Japonicum/gapseq_transporter')

submit_cmds_txt = '%s/gapseq_find-transport_cmds_%s.txt' % (js_dir, cmd_per_js)

# get gapseq_cmd_list
gapseq_cmd_list = []
for gnm_meta in open(japo_metadata_txt):
    if not gnm_meta.startswith('ID\t'):
        gnm_meta_split = gnm_meta.strip().split('\t')
        if gnm_meta_split[2] != 'missing':
            gnm_id     = gnm_meta_split[0]
            pwd_prot   = '%s/%s.protein' % (prot_dir, gnm_id)
            gapseq_cmd = '%s find-transport -b 100 -c 70 -K %s -M prot -f %s %s > %s_find-transport_stdout.txt' % (gapseq_exe, js_cpu, gapseq_tmp_dir, pwd_prot, gnm_id)
            gapseq_cmd_list.append(gapseq_cmd)

# write out gapseq cmds
submit_cmds_txt_handle = open(submit_cmds_txt, 'w')
js_index =1
for each_sub_list in list(divide_list(gapseq_cmd_list, cmd_per_js)):
    each_sub_list_as_str = ';'.join(each_sub_list)
    submit_cmd = 'submitHPC.sh --cmd "%s" -n %s -c transport_%s' % (each_sub_list_as_str, js_cpu, js_index)
    submit_cmds_txt_handle.write(submit_cmd + '\n')

    if submit_js is True:
        print(submit_cmd)
        #os.system(submit_cmd)

    js_index += 1
submit_cmds_txt_handle.close()
