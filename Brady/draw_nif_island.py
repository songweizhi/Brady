import os
from collections import defaultdict
import pandas as pd
from Bio import SeqIO
from dna_features_viewer import GraphicFeature, GraphicRecord
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


# get nifH locus, several genomes have more than one nifH gene, their genoem was named like species_name|1
def get_nifH_id(df, spe_order):
    nifH_id = defaultdict(lambda: list)
    for s in spe_order:
        ns = s.replace('.', '')
        if ns == 'Bradyrhizobium_sp_R22-H':
            ns = 'Bradyrhizobium_sp_R2.2-H'
        if ns == 'Bradyrhizobium_sp_NAS962':
            ns = 'Bradyrhizobium_sp_NAS96.2'
        if ns not in df.columns.tolist():
            print(s, 'not in merged_hmm_info')
        if pd.isna(df.loc['K02588', ns]):
            print(s, 'doesn\'t have nifH gene')
        else:
            locus = df.loc['K02588', ns].split(',')
            dup = 0
            for locu in locus:
                if len(locus) > 1:
                    dup = dup + 1
                    spedup = s + '|' + str(dup)
                else:
                    spedup = s
                nifH_id[spedup] = locu.split('|')[-1]
    return nifH_id


# get gene name
def get_locas2ko(ko_info_tab, df):
    kotbl = [line.strip() for line in open(ko_info_tab, 'r')]
    koname = {}
    for line in kotbl:
        ko = line.split('\t')[0].strip('ko:')
        fullname = line.split('\t')[1]
        name = fullname.split(';')[0]
        koname[ko] = name
    locas2ko = defaultdict()
    for i in df.index.tolist():
        row = df.loc[i]
        for info in row:
            if not pd.isna(info):
                if ',' in info:
                    para = info.split(',')
                    for a in para:
                        locus = a.split('|')[1]
                        locas2ko[locus] = i
                else:
                    locus = info.split('|')[1]
                    locas2ko[locus] = i
    return koname, locas2ko


def get_correct_name(key, indir):
    filelist = os.listdir(indir)
    if key + '.gbk' in filelist:
        fh = indir + key + '.gbk'
        return fh
    else:
        key = key.replace('sp', 'sp.')
        if key + '.gbk' in filelist:
            fh = indir + key + '.gbk'
            return fh
        else:
            key = '_'.join(key.split('_')[:-1])
            fh = indir + key + '.gbk'
            return fh


def trans_features_to_location(CDS, locas2ko, koname):
    '''transform features to location'''
    location = defaultdict(lambda: defaultdict())
    for feature in CDS:
        locus_tag = feature.qualifiers['locus_tag'][0]
        location[locus_tag]['start'] = int(feature.location.start)
        location[locus_tag]['end'] = int(feature.location.end)
        location[locus_tag]['strand'] = int(feature.location.strand)
        if locus_tag in locas2ko.keys():
            location[locus_tag]['label'] = koname[locas2ko[locus_tag]]
        else:
            location[locus_tag]['label'] = 'unknown'
    return location


def select_contig(nifH_id, indir, feature_type, locas2ko, koname):
    locations = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe, nifhid in nifH_id.items():
        contig = []
        filename = spe.split('|')[0]
        fh = get_correct_name(filename, indir)
        record = list(SeqIO.parse(fh, format='genbank'))
        for n in range(len(record)):
            cds_num = 0
            for (index, feature) in enumerate(record[n].features):
                if feature.type == feature_type:
                    cds_num = cds_num + 1
                    if feature.qualifiers['locus_tag'][0] == nifhid:  # find the contig
                        contig = record[n]
                        locate = cds_num
        CDS = []
        for (index_, feature_) in enumerate(contig.features):
            if feature_.type == feature_type:
                CDS.append(feature_)
        if len(CDS) <= 180:
            location = trans_features_to_location(CDS, locas2ko, koname)
        else:
            if locate < 80:
                CDS = CDS[:locate + 80]
            elif locate > (len(CDS) - 80):
                CDS = CDS[locate - 80:]
            else:
                CDS = CDS[locate - 80:locate + 80]
            location = trans_features_to_location(CDS, locas2ko, koname)
        locations[spe] = location
        print(spe, len(locations[spe]))
    return locations


# reverse chromosome and redefine the starting point
def reverse_contig(locations, nifH_id):
    locations_sort = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations.keys():
        nifhid = nifH_id[spe]
        origin = list(locations[spe].values())[0]['start']
        terminal = list(locations[spe].values())[-1]['end']
        if locations[spe][nifhid]['strand'] == 1:
            print(spe, 'forward', origin)
            for locus_tag in locations[spe].keys():
                locations_sort[spe][locus_tag]['start'] = locations[spe][locus_tag]['start'] - origin
                locations_sort[spe][locus_tag]['end'] = locations[spe][locus_tag]['end'] - origin
                locations_sort[spe][locus_tag]['strand'] = locations[spe][locus_tag]['strand']
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
        else:
            print(spe, 'backward', terminal)
            for locus_tag in sorted(locations[spe].keys(), reverse=True):
                locations_sort[spe][locus_tag]['start'] = terminal - locations[spe][locus_tag]['end']
                locations_sort[spe][locus_tag]['end'] = terminal - locations[spe][locus_tag]['start']
                locations_sort[spe][locus_tag]['strand'] = -(locations[spe][locus_tag]['strand'])
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
    return locations_sort


# cut genome region based on gene numbers
def cut_contig_by_genenumber(locations_sort, lf, rf, nifH_id):
    cut_n = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations_sort.keys():
        nifhid = nifH_id[spe]
        count = 0
        locate = 0
        for locus_tag in locations_sort[spe].keys():
            count = count + 1
            if locus_tag == nifhid:
                locate = count
        for (index, locus_tag) in enumerate(locations_sort[spe]):
            if locate < lf:
                left = 0
                right = locate + rf
            elif locate > (len(locations_sort[spe].keys()) - rf):
                left = locate - lf
                right = len(locations_sort[spe].keys())
            else:
                left = locate - lf
                right = locate + rf
            for k, v in list(locations_sort[spe].items())[left:right]:
                cut_n[spe][k] = v
    return cut_n


# draw picture
def plot(start_sites, cut, length, plot_out, fig_size, color):
    xdh = ['BRADO5451', 'BRAD285_1784', 'SE92_10885', 'BDOA9_0161380']
    fig, ax = plt.subplots(len(start_sites), 1, figsize=(fig_size, len(start_sites) * 1.5))
    i = 0
    for spe in start_sites:
        fea = []
        if start_sites[spe] == 0:
            origin = list(cut[spe].values())[0]['end']
        else:
            origin = start_sites[spe]
        for locus_tag in cut[spe].keys():
            s = cut[spe][locus_tag]['end']
            t = cut[spe][locus_tag]['end']
            if t - origin > length or s < origin:
                pass
            else:
                a = cut[spe][locus_tag]['start'] - start_sites[spe]
                z = cut[spe][locus_tag]['end'] - start_sites[spe]
                d = cut[spe][locus_tag]['strand']
                if locus_tag in xdh:
                    l = 'xdh'
                else:
                    l = cut[spe][locus_tag]['label'].split(',')[0]
                if l == 'unknown':
                    gf = GraphicFeature(start=a, end=z, strand=d, color="#FFFFCC")
                elif l not in color.keys():
                    gf = GraphicFeature(start=a, end=z, strand=d, color="#FFFFCC", label=l)
                else:
                    c = color[l]
                    gf = GraphicFeature(start=a, end=z, strand=d, color=c, label=l)
                fea.append(gf)
        record = GraphicRecord(sequence_length=length, features=fea)
        for f in record.features:
            f.data['fixed_level'] = 0
        ax[i].set_xlabel(spe.replace('_', ' '))
        record.plot(ax=ax[i])
        i = i + 1
    plt.savefig(plot_out)


def plot_gbk(gene_color_txt, annotation, ko_info_tab, gnm_id_txt, gbk_dir, gbk_ext, feature_type, lf, rf, flnif_dict, plot_out):

    gnm_list = []
    for each_gnm in open(gnm_id_txt):
        gnm_list.append(each_gnm.strip())

    gene_color_dict = dict()
    if os.path.isfile(gene_color_txt):
        for each_gene in open(gene_color_txt):
            each_gene_split = each_gene.strip().split('\t')
            gene_color_dict[each_gene_split[0]] = each_gene_split[1]
    else:
        print('color file not found, keep plotting anyway!')

    protein_info = pd.read_csv(annotation, sep='\t', header=0, index_col=0, low_memory=False)
    koname, locas2ko = get_locas2ko(ko_info_tab, protein_info)
    nifH_id = get_nifH_id(protein_info, gnm_list)
    locations = select_contig(nifH_id, gbk_dir, feature_type, locas2ko, koname)
    locations_sort = reverse_contig(locations, nifH_id)
    cut_n = cut_contig_by_genenumber(locations_sort, lf, rf, nifH_id)

    plot(flnif_dict, cut_n, 50000, plot_out, 56, gene_color_dict)

    print('Done!')


############################################# main #############################################

# file in
wd                  = '/Users/songweizhi/Desktop/Brady/draw_nif_island_wd'
gbk_dir             = '%s/gbk_files/'           % wd
gbk_ext             = 'gbk'
gnm_id_txt          = '%s/gnm_id.txt'           % wd
gene_color_txt      = '%s/gene_color.txt'       % wd
feature_type        = 'CDS'
gene_num_to_plot_l  = 65
gene_num_to_plot_r  = 50

# file out
plot_out            = '%s/flnif.pdf'            % wd

# the other files
annotation          = '%s/merged_hmm_info.tab'  % wd
ko_info_tab         = '%s/ko_info.tab'          % wd
flnif_dict          = {'Bradyrhizobium_sp._S23321|1': 31195, 'Bradyrhizobium_sp._S23321|2': 62979, 'Bradyrhizobium_sp._AT1': 46651}

plot_gbk(gene_color_txt, annotation, ko_info_tab, gnm_id_txt, gbk_dir, gbk_ext, feature_type, gene_num_to_plot_l, gene_num_to_plot_r, flnif_dict, plot_out)
