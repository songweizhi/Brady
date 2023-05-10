import os
import glob
import random
from Bio import SeqIO
from ete3 import Tree


def sep_path_basename_ext(file_in):

    # separate path and file name
    f_path, file_name = os.path.split(file_in)
    if f_path == '':
        f_path = '.'

    # separate file basename and extension
    f_base, f_ext = os.path.splitext(file_name)

    return f_path, f_base, f_ext


def subset_tree_by_treecluster(tree_file_in, TreeCluster_op_txt, tips_to_keep_txt, tree_file_out):

    specified_to_keep_node_set = set()
    if os.path.isfile(tips_to_keep_txt) is True:
        for each in open(tips_to_keep_txt):
            specified_to_keep_node_set.add(each.strip())

    gnm_to_keep_set = set()
    cluster_to_gnm_dict = dict()
    for each_line in open(TreeCluster_op_txt):
        if not each_line.startswith('SequenceName	ClusterNumber'):
            each_line_split = each_line.strip().split('\t')
            gnm_id = each_line_split[0]
            cluster_id = each_line_split[1]
            if cluster_id == '-1':
                gnm_to_keep_set.add(gnm_id)
            else:
                if cluster_id not in cluster_to_gnm_dict:
                    cluster_to_gnm_dict[cluster_id] = set()
                cluster_to_gnm_dict[cluster_id].add(gnm_id)

    for each_c in cluster_to_gnm_dict:
        c_gnm_set = cluster_to_gnm_dict[each_c]

        current_c_gnm_to_keep_set = set()
        for each_g in c_gnm_set:
            if each_g in specified_to_keep_node_set:
                current_c_gnm_to_keep_set.add(each_g)

        if len(current_c_gnm_to_keep_set) == 0:
            current_c_gnm_to_keep_set.add(random.choice([i for i in c_gnm_set]))

        for each in current_c_gnm_to_keep_set:
            gnm_to_keep_set.add(each)

    # subsample tree
    input_tree = Tree(tree_file_in)
    subset_tree = input_tree.copy()
    subset_tree.prune(gnm_to_keep_set, preserve_branch_length=True)
    subset_tree.write(outfile=tree_file_out)


def rm_problematic_gnms(og_dir, og_ext, to_keep_gnm_set, og_dir_updated, force_mkdir):

    if os.path.isdir(og_dir_updated) is True:
        if force_mkdir is False:
            print('%s detected, program exited!' % og_dir_updated)
            exit()
        else:
            os.system('rm -r %s' % og_dir_updated)
    os.system('mkdir %s' % og_dir_updated)

    og_file_re = '%s/*.%s' % (og_dir, og_ext)
    og_file_list = glob.glob(og_file_re)

    for og_file in og_file_list:
        file_path, file_basename, file_ext = sep_path_basename_ext(og_file)
        og_file_updated = '%s/%s%s' % (og_dir_updated, file_basename, file_ext)
        og_file_updated_handle = open(og_file_updated, 'w')
        for each_seq in SeqIO.parse(og_file, 'fasta'):

            gnm_id = each_seq.id

            if gnm_id in to_keep_gnm_set:
                og_file_updated_handle.write('>%s\n' % each_seq.id)
                og_file_updated_handle.write('%s\n' % str(each_seq.seq))
        og_file_updated_handle.close()


########################################################################################################################

# file in
full_tree_file              = '/Users/songweizhi/Desktop/subsample_gnm_wd/iqtree.treefile'
tree_file_pb                = '/Users/songweizhi/Desktop/subsample_gnm_wd/iqtree_PB.treefile'
tree_file_sym1              = '/Users/songweizhi/Desktop/subsample_gnm_wd/iqtree_sym1.treefile'
nodes_to_keep_txt           = '/Users/songweizhi/Desktop/subsample_gnm_wd/nodes_to_keep.txt'
TreeCluster_op_txt_pb       = '/Users/songweizhi/Desktop/subsample_gnm_wd/iqtree_PB_clusters_0.00001_by_len.txt'
TreeCluster_op_txt_sym1     = '/Users/songweizhi/Desktop/subsample_gnm_wd/iqtree_sym1_clusters_0.0001_by_len.txt'
og_dir                      = '/Users/songweizhi/Desktop/subsample_gnm_wd/pep-clean'
og_ext                      = 'fas'
force_mkdir                 = True

# file out
op_dir                      = '/Users/songweizhi/Desktop/subsample_gnm_wd/output'
tree_file_subsampled_pb     = '%s/iqtree_PB_subsampled_0.00001_by_len.treefile'     % op_dir
tree_file_subsampled_sym1   = '%s/iqtree_sym1_subsampled_0.0001_by_len.treefile'    % op_dir
og_dir_updated              = '%s/pep-clean_subsampled'                             % op_dir

########################################################################################################################

if os.path.isdir(op_dir) is True:
    if force_mkdir is True:
        os.system('rm -r %s' % op_dir)
    else:
        print('%s found, program exited!' % op_dir)
        exit()
os.system('mkdir %s' % op_dir)

# subsample tree
subset_tree_by_treecluster(tree_file_pb, TreeCluster_op_txt_pb, nodes_to_keep_txt, tree_file_subsampled_pb)
subset_tree_by_treecluster(tree_file_sym1, TreeCluster_op_txt_sym1, nodes_to_keep_txt, tree_file_subsampled_sym1)

# get genomes on unsubsampled trees
leaves_overall  = [i.name for i in Tree(full_tree_file).get_leaves()]
leaves_pb       = [i.name for i in Tree(tree_file_pb).get_leaves()]
leaves_sym1     = [i.name for i in Tree(tree_file_sym1).get_leaves()]

# get genomes on subsampled trees
leaves_pb_subsampled   = [i.name for i in Tree(tree_file_subsampled_pb).get_leaves()]
leaves_sym1_subsampled = [i.name for i in Tree(tree_file_subsampled_sym1).get_leaves()]

# get all genomes to keep in the main input tree
gnm_set_to_keep = set()
for each_gnm in leaves_overall:
    if (each_gnm not in leaves_pb) and (each_gnm not in leaves_sym1):
        gnm_set_to_keep.add(each_gnm)
    elif each_gnm in leaves_pb:
        if each_gnm in leaves_pb_subsampled:
            gnm_set_to_keep.add(each_gnm)
    elif each_gnm in leaves_sym1:
        if each_gnm in leaves_sym1_subsampled:
            gnm_set_to_keep.add(each_gnm)

# rm_problematic_gnms from sequence
rm_problematic_gnms(og_dir, og_ext, gnm_set_to_keep, og_dir_updated, force_mkdir)

