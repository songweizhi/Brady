import argparse


def best_hit(file_in, file_out):

    file_out_handle = open(file_out, 'w')
    best_hit_line = ''
    best_hit_query_id = ''
    best_hit_score = 0
    for blast_hit in open(file_in):
        blast_hit_split = blast_hit.strip().split('\t')
        query_id = blast_hit_split[0]
        bit_score = float(blast_hit_split[11])

        if best_hit_query_id == '':
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

        elif (query_id == best_hit_query_id) and (bit_score > best_hit_score):
            best_hit_score = bit_score
            best_hit_line = blast_hit

        elif query_id != best_hit_query_id:
            file_out_handle.write(best_hit_line)
            best_hit_query_id = query_id
            best_hit_line = blast_hit
            best_hit_score = bit_score

    file_out_handle.write(best_hit_line)
    file_out_handle.close()


def hkcct2speciesName(args):

    blast_results = args['i']
    op_table = args['o']

    # file out
    blast_results_best_hit = blast_results + '.besthit'

    best_hit(blast_results, blast_results_best_hit)

    hkcct2speciesName_tbl_handle = open(op_table, 'w')
    processed_gnm = set()
    for each_line in open(blast_results_best_hit):
        each_line_split = each_line.strip().split('\t')
        gnm_id = each_line_split[0].split('.genome_16S_')[0]
        tax_str = each_line_split[1]
        if gnm_id not in processed_gnm:
            tax_str_split = tax_str.split(';')
            species_name = tax_str_split[-1]
            if 'Bradyrhizobium_' not in species_name:
                species_name = 'Bradyrhizobium_' + species_name

            processed_gnm.add(gnm_id)
            id_with_taxon = species_name + '_' + gnm_id
            hkcct2speciesName_tbl_handle.write('%s\t%s\n' % (gnm_id, id_with_taxon))

    hkcct2speciesName_tbl_handle.close()


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',     required=True,  help='blast results ')
    parser.add_argument('-o',     required=True,  help='hkcct2speciesName.tbl')
    args = vars(parser.parse_args())
    hkcct2speciesName(args)
