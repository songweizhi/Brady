import os
from time import sleep

n = 0
while n < 1:
    sleep(600)
    cmd = 'python3 /home-user/wzsong/Scripts/rm_tmp_dir.py -d /home-user/wzsong/Japonicum/gapseq_metacyc/gapseq_tmp_dir -l 20m'
    os.system(cmd)
