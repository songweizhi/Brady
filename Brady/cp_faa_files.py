import os
import glob
import multiprocessing as mp


faa_id_txt = '/home-user/wzsong/Brady/PB_gnm_list.txt'

for each_gnm in open(faa_id_txt):

    gnm_id = each_gnm.strip()
    pwd_faa = 'faa_files/%s.protein' % gnm_id
    cp_cmd = 'cp %s faa_files_394_PB/' % pwd_faa
    # os.system(cp_cmd)

