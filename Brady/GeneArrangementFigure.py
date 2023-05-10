#! /usr/bin/env python

color = {'xdh':'#FFFF00',
         'fabG' : '#FF3399', #Catalyzes the NADPH-dependent reduction of beta-ketoacyl-ACP substrates to beta-hydroxyacyl-ACP products, the first reductive step in the elongation cycle of fatty acid biosynthesis.
         'nifA' : '#33CC99', #synthesize nitrogenase
         'nifB' : '#33CC99',
         'nifD' : '#33CC99',
         'nifK' : '#33CC99',
         'nifH' : '#33CC99',
         'nifQ' : '#33CC99',
         'nifE' : '#33CC99',
         'nifN' : '#33CC99',
         'nifT' : '#33CC99',
         'nifU-like protein' : '#33CC99',
         'nifV' : '#33CC99',
         'nifW' : '#33CC99',
         'nifX' : '#33CC99',
         'nifZ' : '#33CC99',
         'fixA' : '#99CC00', # transfer electron to nitrogenase
         'fixB' : '#99CC00',
         'fixC' : '#99CC00',
         'fixX' : '#99CC00',
         'fixL' : '#99CC00',
         'fixJ' : '#99CC00',
         'sufB' : '#FFD700', # FeS cluster assembly protein
         'sufC' : '#FFD700',
         'sufS' : '#FFD700',
         'sufD' : '#FFD700',
         'sufE' : '#FFD700',
         'sufX' : '#FFD700',
         'iscA' : '#FFD700',
         'iscS' : '#FFD700',
         'hcaC' : '#FFFF00', # like ferredoxin
         'fdx' : '#FFFF00', # ferredoxin  transfer electron
         'per' : '#FFC0CB', # synthesize sugar
         'cysE' : '#FF7F50', # synthesizes L-cysteine
         'paaD' : '#CD853F', # phenylacetate degradation, part of Aromatic compound metabolism
         'bolA': '#66CCFF', # DNA-binding transcriptional regulator BolA
         'grxD' : '#FFD700', # Monothiol glutaredoxin involved in the biogenesis of iron-sulfur clusters.
         'irr' : '#66CCFF', #Acts as a global negative controlling element, employing Fe2+ as a cofactor to bind the operator of the repressed genes.
         'hyaA' : '#99FFFF',
         'hyaB' : '#99FFFF',
         'hyaC' : '#99FFFF',
         'hyaD' : '#99FFFF',
         'hyaF' : '#99FFFF',
         'hypA' : '#99FFFF',
         'hypB' : '#99FFFF',
         'hypBA1': '#99FFFF',
         'hypC' : '#99FFFF',
         'hypD' : '#99FFFF',
         'hypE' : '#99FFFF',
         'hypF' : '#99FFFF',
         'hypV' : '#99FFFF',
         'hupV' : '#99FFFF',
         'glbN' : '#CC0000',
         'glbO' : '#CC0000',
         'hspQ' : '#1E90FF',
         'modA' : '#CC99FF',
         'modB' : '#CC99FF',
         'modC' : '#CC99FF',
         'modD' : '#CC99FF',
         'modE' : '#CC99FF',
         'nodA' : '#F08080',
         'nodB' : '#F08080',
         'nodC' : '#F08080',
         'nodD' : '#F08080',
         'nodI' : '#F08080',
         'nodJ' : '#F08080',
          }

#####################################################################################
# The script is originally written by Jinjin Tao from CUHK.
# It is further modified by Sishuo Wang (sishuowang@hotmail.ca) to add a few more functions allowing more easily specifying parameters.
# All rights are reserved.


#####################################################################################
import os
import sys
import getopt
import os.path
import re

import collections
from collections import defaultdict
from tqdm import tqdm
import pandas as pd

from Bio import SeqIO

from matplotlib.backends.backend_pdf import PdfPages
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

from dna_features_viewer import GraphicFeature, GraphicRecord,CircularGraphicRecord


#####################################################################################
# user-defined
#from GeneColorDict import color

sys.path.append('lib')
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), './lib')))
import user_defined
from user_defined import get_gene_id, create_flnif, get_list


#####################################################################################
#get gene name
def get_locas2ko(df):
    #kotbl = [line.strip() for line in open("/home-user/thliao/data/protein_db/kegg/ko_info.tab", 'r')]
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
                        spe = a.split('|')[0]
                        locus = a.split('|')[-1]
                        locas2ko[locus]= i
                else:
                    spe = info.split('|')[0]
                    locus = info.split('|')[-1]
                    locas2ko[locus]= i
    return koname,locas2ko

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
        
def trans_features_to_location(CDS):
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

def select_contig(nifH_id):
    locations = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    tqdm.write("reading gbk files")
    for spe,nifhid in nifH_id.items():
            contig = []
            print(spe)
            #spe = 'Bradyrhizobium_sp._AT1'
            #nifHid = 'SE92_10630'
            filename = spe.split('|')[0]
            fh = get_correct_name(filename, indir) # brady_sp_123.gbk
            if not os.path.exists(fh):
                continue
            record = list(SeqIO.parse(fh,format='genbank'))
            for n in range(len(record)):
                cds_num = 0         
                for (index, feature) in enumerate(record[n].features):
                    if feature.type == feature_type:
                        cds_num = cds_num +1
                        print(feature.qualifiers['locus_tag'][0], nifhid)
                        if feature.qualifiers['locus_tag'][0] == nifhid: #find the contig
                            #contig.append record[n]
                            contig = record[n]
                            print(record)
                            locate = cds_num 
            CDS = []                                   
            print(spe, nifhid)
            for (index_, feature_) in enumerate(contig.features):
                if feature_.type == feature_type:
                    CDS.append(feature_)
            if len(CDS) <= 180:
                location = trans_features_to_location(CDS)
            else:
                if locate < 80:
                    CDS = CDS[:locate+80]
                elif locate > (len(CDS)-80):
                    CDS = CDS[locate-80:]
                else:
                    CDS = CDS[locate-80:locate+80]
                location = trans_features_to_location(CDS)   
            locations[spe] = location
            #print(spe, len(locations[spe]))
    return locations

#reverse chromosome and redefine the starting point
def reverse_contig(locations):
    locations_sort = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations.keys():
        #spe = 'Bradyrhizobium_sp._S23321'
        nifhid = nifH_id[spe]
        origin = list(locations[spe].values())[0]['start']
        terminal = list(locations[spe].values())[-1]['end']
        if locations[spe][nifhid]['strand'] == 1:
            #print(spe,'forward', origin)
            for locus_tag in locations[spe].keys():
                locations_sort[spe][locus_tag]['start'] = locations[spe][locus_tag]['start']- origin
                locations_sort[spe][locus_tag]['end'] = locations[spe][locus_tag]['end'] - origin
                locations_sort[spe][locus_tag]['strand'] = locations[spe][locus_tag]['strand']
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
        else:   
            #print(spe,'backward',terminal)
            for locus_tag in sorted(locations[spe].keys(),reverse = True):
                locations_sort[spe][locus_tag]['start'] = terminal - locations[spe][locus_tag]['end']
                locations_sort[spe][locus_tag]['end'] = terminal - locations[spe][locus_tag]['start'] 
                locations_sort[spe][locus_tag]['strand'] = -(locations[spe][locus_tag]['strand'])
                locations_sort[spe][locus_tag]['label'] = locations[spe][locus_tag]['label']
    return locations_sort

#cut genome region based on gene numbers 
def cut_contig_by_genenumber(locations_sort,lf,rf):
    cut_n = defaultdict(lambda: defaultdict(lambda: defaultdict()))
    for spe in locations_sort.keys(): 
        nifhid = nifH_id[spe]
        count = 0
        for locus_tag in locations_sort[spe].keys():
            count = count +1
            if locus_tag == nifhid:
                locate = count                                         
        for (index, locus_tag) in enumerate(locations_sort[spe]):  
            if locate < lf:
                left = 0
                right = locate+rf
            elif locate > (len(locations_sort[spe].keys())-rf):
                left = locate-lf
                right = len(locations_sort[spe].keys())
            else:
                left = locate-lf
                right = locate+rf      
            for k,v in list(locations_sort[spe].items())[left:right]:
                cut_n[spe][k] = v
    return cut_n

#draw picture 
def plot(start_sites,cut,length,group_name,fig_size):
    xdh = ['BRADO5451','BRAD285_1784','SE92_10885','BDOA9_0161380']
    fig, ax = plt.subplots(len(start_sites),1,figsize=(fig_size,len(start_sites)*2))
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
            if t-origin > length or s < origin:
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
                    c = "#FFFFCC"
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c)          
                elif l not in color.keys():            
                    c = "#FFFFCC"
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c,
                                        label = l 
                                        )
                
                else:
                    c = color[l]           
                    gf = GraphicFeature(start=a,
                                        end=z,
                                        strand=d,
                                        color=c,
                                        label=l)
                fea.append(gf)
        record = GraphicRecord(sequence_length=length, features=fea)
        for f in record.features:
            f.data['fixed_level'] = 0
        ax[i].set_xlabel(spe.replace('_',' '))
        record.plot(ax = ax[i])
        i = i+1
    plt.savefig('flnif.pdf')
    #plt.savefig(f'{group_name}.pdf')
    return('done')


############################################# main #############################################
# default
indir = '/home-user/jjtao/Rhizobiales/data/Rhizobiales/gbk/' 
annotation = '/home-user/jjtao/Rhizobiales/kegg_hmmsearch/Brady345/after_rename_3Div/1e10/merged_hmm_info.tab'
kotbl = [line.strip() for line in open("/home-user/thliao/data/protein_db/kegg/ko_info.tab", 'r')]
species_list_file = None
length = 50000
is_one = False

nifH_KO = 'K02588'
nifD_KO = 'K02586'

options, args = getopt.getopt(sys.argv[1:], "", ['gbk_indir=', 'kotbl=', 'ko_tbl=', 'ko_info=', 'annotation=', 'species_list=', 'length=', 'one'])
for opt, arg in options:
    if opt in ('--gbk_indir'):
        indir = arg
    elif opt in ('--kotbl', '--ko_tbl', '--ko_info'):
        kotbl = arg
    elif opt in ('--annotation'):
        annotation = arg
    elif opt in ('--species_list'):
        species_list_file = arg
    elif opt in ('--length'):
        length = int(arg)
    elif opt in ('--one'):
        is_one = True
    elif opt in ('--gene1'):
        nifH_KO = arg
    elif opt in ('--gene2'):
        nifD_KO = arg


#############################################
## species list
#my = ['Bradyrhizobium_jicamae_SZCCHNW2005']
my = get_list(species_list_file)

#spe_order = PB + flnif + flwithsym + SymBasal + sym 
spe_order = my

## nifH id
protein_info = pd.read_csv(annotation,sep='\t',header=0,index_col=0,low_memory=False)
#print(protein_info['Bradyrhizobium_jicamae_SZCCHNW2005'])

nifH_id = get_gene_id(protein_info, spe_order, nifH_KO) # gene name is actually KO id
nifD_id = get_gene_id(protein_info, spe_order, nifD_KO) # gene name is the one used to locate the genomic context where nifH is nearby
#nifH_id = get_gene_id(protein_info, spe_order, 'K02586') # SSW
#nifD_id = get_gene_id(protein_info, spe_order, 'K02586') # SSW

koname,locas2ko = get_locas2ko(protein_info)
#print(locas2ko); sys.exit()

##read gbk and get contig 
#indir = '/home-user/jjtao/Rhizobiales/data/Rhizobiales/gbk/' 
feature_type = 'CDS'
locations = select_contig(nifH_id)

locations_sort = reverse_contig(locations)
lf = 200 #gene numbers in left flank of nifH
rf = 200 #gene numbers in right flank of nifH
cut_n = cut_contig_by_genenumber(locations_sort,lf,rf)


################################### draw picture ####################################
'''
flnif = {
	'Bradyrhizobium_sp._S23321|1':31195,
	'Bradyrhizobium_sp._S23321|2':62979,        
	#'Bradyrhizobium_japonicum_22':-3225,
	'Bradyrhizobium_japonicum_22':28505,
	'Bradyrhizobium_sp_DOA9_chromosome':35549,
	'Bradyrhizobium_guangxiense_CCBAU_53363_CP022219':29168,      
	'Bradyrhizobium_iriomotense_SZCCT0007':-3225,
	'Bradyrhizobium_sacchari_p9-20_LWIG01000010':-3225,
	'Bradyrhizobium_sp._AT1':46651,
	'Bradyrhizobium_sp._BM-T':35390,
	'Bradyrhizobium_sp_DOA9_chromosome|2':39546,
    }
'''

flnif = create_flnif(cut_n, nifD_id, 12000)
#flnif = create_flnif(cut_n, nifH_id, 12000) #SSW

if is_one:
    print(flnif);
    for k in list(flnif.keys()):
        v = flnif[k]
        if re.search('\|[2-9]+$', k):
            flnif.pop(k)

plot(flnif, cut_n, length, 'flnif', 56)


