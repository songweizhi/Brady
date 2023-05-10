
isolate_metadata_txt = '/Users/songweizhi/PycharmProjects/Brady/metadata/isolate_metadata.txt'


# read in isolate_metadata_txt
duplicated_isolate_set = set()
isolate_to_source_dict = dict()
source_to_isolate_dict = dict()
for each_isolate in open(isolate_metadata_txt):

    each_isolate_split = each_isolate.strip().split('\t')
    isolate_name = each_isolate_split[0]
    isolate_source = each_isolate_split[1]

    if isolate_name in isolate_to_source_dict:
        duplicated_isolate_set.add(isolate_name)

    isolate_to_source_dict[isolate_name] = isolate_source
    if isolate_source not in source_to_isolate_dict:
        source_to_isolate_dict[isolate_source] = {isolate_name}
    else:
        source_to_isolate_dict[isolate_source].add(isolate_name)

# check duplication
if len(duplicated_isolate_set) > 0:
    print('Duplicated isolates:')
    print('n'.join(duplicated_isolate_set))
    print('The above isolates were found multiple times in metadata')
    print('Program exited!')

# print out summary
total_num = 0
for each_source in source_to_isolate_dict:
    current_isolate_set = source_to_isolate_dict[each_source]
    total_num += len(current_isolate_set)
    print('%s\t%s' % (each_source, len(current_isolate_set)))
print('Total\t%s' % total_num)


ref_seq         = '/home-user/wzsong/Brady/metadata_and_ref_Seq/shc.fa'
faa_dir         = '/home-user/wzsong/Brady/faa_files'
blast_op_dir    = 'op_dir'


for each_gnm in isolate_to_source_dict:
    pwd_faa = '%s/%s.protein' % (faa_dir, each_gnm)
    print(pwd_faa)
    blast_op = '%s/' % blast_op_dir

    blastp_cmd = 'blastp -query %s -subject %s -out %s -outfmt 6 -evalue %s' % (ref_seq, pwd_faa, blast_op, evalue_cutoff)
