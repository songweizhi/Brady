import random


file_in = '/Users/songweizhi/Desktop/AAA.txt'


mc_list = []
color_list = []
for each_line in open(file_in):
    each_line_split = each_line.strip().split('\t')
    mc_list.append(each_line_split[0])
    color_list.append(each_line_split[1])


mc_list_r = random.sample(mc_list, len(mc_list))
color_list_r = random.sample(color_list, len(color_list))


for mc, color in zip(mc_list_r, color_list_r):
    print('%s\t%s' % (mc, color))

