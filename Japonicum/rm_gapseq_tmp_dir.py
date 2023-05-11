import os
import time


def get_time_since_last_modification(target_folder_or_file):

    last_modified_time = os.path.getmtime(target_folder_or_file)
    current_time = time.time()
    tslm_sec  = current_time - last_modified_time
    tslm_min  = tslm_sec / 60
    tslm_hour = tslm_sec / (60 * 60)
    tslm_day  = tslm_sec / (60 * 60 * 24)
    tslm_sec  = float("{0:.2f}".format(tslm_sec))
    tslm_min  = float("{0:.2f}".format(tslm_min))
    tslm_hour = float("{0:.2f}".format(tslm_hour))
    tslm_day  = float("{0:.2f}".format(tslm_day))

    return tslm_sec, tslm_min, tslm_hour, tslm_day


gapseq_tmp_dir = '/home-user/wzsong/Japonicum/gapseq_metacyc/gapseq_tmp_dir'

sub_dir_list = next(os.walk(gapseq_tmp_dir))[1]

for each_sub_dir in sub_dir_list:
    pwd_sub_dir = '%s/%s' % (gapseq_tmp_dir, each_sub_dir)
    subdir_tslm_sec, subdir_tslm_min, subdir_tslm_hour, subdir_tslm_day = get_time_since_last_modification(pwd_sub_dir)
    if subdir_tslm_min >= 30:
        print('Removing %s' % each_sub_dir)
        os.system('rm -r %s' % pwd_sub_dir)
