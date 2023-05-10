import os
import glob
import argparse


def sep_path_basename_ext(file_in):

    # separate path and file name
    file_path, file_name = os.path.split(file_in)
    if file_path == '':
        file_path = '.'

    # separate file basename and extension
    file_basename, file_extension = os.path.splitext(file_name)

    return file_path, file_basename, file_extension


def check_paralog(args):

    blast_op_dir        = args['i']
    blast_op_ext        = args['x']
    marker_gene         = args['m']
    evalue_index_cutoff = args['e']

    blast_op_file_re = '%s/*.%s' % (blast_op_dir, blast_op_ext)
    blast_op_file_list = glob.glob(blast_op_file_re)

    gnm_set_with_paralog = set()
    qualified_hit_dict = dict()
    for each_blast_op_file in blast_op_file_list:

        file_path, file_basename, file_extension = sep_path_basename_ext(each_blast_op_file)

        identified_nifB_num = 0
        qualified_hit_list = []
        for each_hit in open(each_blast_op_file):
            each_hit_split = each_hit.strip().split('\t')
            ref_id = each_hit_split[1]
            e_value = each_hit_split[-2]

            qualified_evalue = False
            if '-' in e_value:
                e_value_index = int(e_value.split('-')[-1])
                if e_value_index >= evalue_index_cutoff:
                    qualified_evalue = True
            elif float(e_value) == 0:
                qualified_evalue = True

            if qualified_evalue is True:
                if 'nif' in ref_id:
                    qualified_hit_list.append(each_hit.strip())

                if ref_id == marker_gene:
                    identified_nifB_num += 1

        qualified_hit_dict[file_basename] = qualified_hit_list

        if identified_nifB_num >= 2:
            gnm_set_with_paralog.add(file_basename)

    for_check_txt = blast_op_dir + '.txt'
    for_check_txt_handle = open(for_check_txt, 'w')
    for gnm in gnm_set_with_paralog:
        current_gnm_qulified_hit_list = qualified_hit_dict[gnm]
        for_check_txt_handle.write('\n'.join(current_gnm_qulified_hit_list))
        for_check_txt_handle.write('\n\n')
    for_check_txt_handle.close()

    # report
    if len(gnm_set_with_paralog) == 0:
        print('No genome with paralog of nif cluster was identified')
    else:
        print('nif cluster paralog was found in the following genomes:')
        print('\n'.join(gnm_set_with_paralog))
        print('nif cluster paralog was found in above genomes')


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',     required=True,                         help='blast results ')
    parser.add_argument('-x',     required=True,                         help='file extension')
    parser.add_argument('-m',     required=True,                         help='marker gene, e.g., nifB')
    parser.add_argument('-e',     required=False, type=int, default=100, help='evalue index cutoff, default:100')
    args = vars(parser.parse_args())
    check_paralog(args)


'''
python3 identify_paralog.py -h
python3 identify_paralog.py -i LingHK360_Zong1070_rpsblast_op_nif -x blast8 -m nifB -e 100
'''

