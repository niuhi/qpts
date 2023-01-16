import bisect # used for finding highest value lower than threshold
import collections # used for OrderedDict(), see https://pymotw.com/2/collections/ordereddict.html
import csv
import glob
import os, re, sys
import shutil
import time
from timeit import default_timer as timer
from pprint import pprint
import pandas as pd
import numpy as np

start = timer()
now = time.strftime('%Y%m%d%H%M%S', time.localtime())

# starting from 'v8.0'
# tree_topology_sorter_v8.0.py n s
  # n for nested, s for sister
  # m for MTA, s for STA

""" before v8.0
thread = int(sys.argv[1]) # 0..all_threads-1
print('''thread={}'''.format(thread))
"""
thread = 0
all_threads = 1

copy_tree = False
# copy_tree = True
"""
if only_nested:
    bs_treshold = 0 # standard value: 50
    relationship = 'nested'
else:
    bs_treshold = 75
    relationship = 'sister'
"""# You can change for testing purposes
only_nested = None
if sys.argv[1] == 'm' or sys.argv[2] == 'm' or sys.argv[3] == 'm':#if 'MeEGT' in os.getcwd():
    eEGT = 'MeEGT'
    #test_size = 8644 # without single MTA
    test_size = 11777 # all MTA; correct as of v8.2
    in_tree_dir =  'Q:/STEPA42db/MTA/'
    # result_dir =  'Q:/____STEPA41db_results/MeEGT/'
else:
    eEGT = 'SeEGT' # Q:\Dropbox\university\q_papers\eEGT_Q2020\STEPA41db\sorting_v8.0\input
    test_size = 76912
    in_tree_dir =  'Q:/STEPA42db/256409212138/STA/'
    #result_dir = 'Q:/____STEPA41db_results/SeEGT/'
if sys.argv[1] == 'n' or sys.argv[2] == 'n' or sys.argv[3] == 'n':
    only_nested = True
    bs_treshold = 50 # standard value: 50
    relationship = 'nested'
    relationship = 'NR' + str(bs_treshold).zfill(2)
else:
    only_nested = False
    bs_treshold = 75 # standard value: 75
    relationship = 'sister'
    relationship = 'SR' + str(bs_treshold).zfill(2)
for i in range(1, len(sys.argv)):
    try:
        bins = int(sys.argv[i])
    except ValueError:
        pass
else:
    pass
script_name = sys.argv[0].split('\\')[-1]
version = script_name[1+script_name.find('_v'):-3]
#lbl = 'dip_nomix_nested_node'
lbl = '' # 'thread_' + str(thread)
settings = eEGT[:1] + str(bins).zfill(2) + relationship
label = '_' + version + '_' + settings
print(label)
#input()
tree_ext = 'bip'

# You are not allowed to change anything below this line!
contamination_dir = 'Q:/Dropbox/university/q_papers/eEGT_Q2020/STEPA42db/contamination/target/' # used in count_algae(), but to slow
trees_done = int(test_size/100+1)
# print(script_name)
# print('''version: {}'''.format(version))
# print('''label: {}'''.format(label))
print('''{}: {} analysis of {}; time stamp: {}'''.format(script_name, relationship, eEGT, now))
#input()
"""
**************************************************************
*******************  TREE OPTIONS  ***************************
**************************************************************
"""



target_tree_file = 'Q:/Dropbox/STEPA41db/scripts/trees/DExx_all_alt.tre'
full_Euglenoidea_tree = '(((((((DEEZ,DEEG),DEEH),DEEE),DERV),(DErc,DEnp)),DEpt),DEpv);'
    # Target tree must be fully resolved, i.e. contains only bipartitions, not multipartitions!!!!
    # i.e. the number of '(' and ')' must be the same as the number of ','
    # does not true for R analysis, only for opening  the tree in treeviewer
# To be correctly processed in R, the full tree of Euglenids have form:
"""((((((((Euglena_longa:1,Euglena_gracilis:1):1,Euglena_hiemalis:1):1,Eutreptiella_gymnastica:1):1,Rapaza_viridis:1):1,Peranema_trichophorum:1):1,(Rhabdomonas_costata:1,Neometanema_parovale:1):1):1,Ploeotia_vitrela:1):1);"""
# in near future the tree should be in this file and not external!!!!!!!!!
# the string must be in form:
# '(((((((DEEZ,DEEG),DEEH),DEEE),DERV),DEpt),DEnp),DEpv);\n'
# with open(target_tree_file) as f:
    # target_tree = f.read().strip('\n')
target_tree = full_Euglenoidea_tree # _no_RCo

"""
**************************************************************
****************  END OF TREE OPTIONS  ***********************
**************************************************************
"""

code_length = 4 # 4-letter codes
Q000_absence = 1

outdir = 'Q:/STEPA42db/256409212138/output/' + settings + '/'
#outdir = '../output/' + settings + '/'
path = outdir + now + '_' + str(test_size) + '/' # + label
result_dir = path # historical reasons
if not os.path.isdir(outdir):
    os.mkdir(outdir)
if not os.path.isdir(path):
    os.mkdir(path)
#pth = now + '_' + str(test_size)
#print(pth)
#print(label)
#input()
#path = result_dir + pth + '/'
#print(path)
log_path = path + 'log/'
path_selected = path + '_selected/'
log_path_selected = path_selected + 'log/' # not_used
nested_trees = path + 'nested_target/'
len_intree = len(in_tree_dir)
trees = glob.glob(in_tree_dir + '*.bip')
# print(trees)
"""
# Now I am sorting only the trees which were preselected; decontamination takes very long time
with open('preselected_trees.q') as f:
    trees = list()
    for line in f:
        trees.append(in_tree_dir[:-1] + '\\' + line.strip())
print(trees[0])
#"""
# print(trees)
# trees = glob.glob(in_tree_dir + 'q1080830.bip')
# trees = glob.glob(in_tree_dir + 'q1000443.bip')

# results_file      = result_dir + now + '_' + str(test_size) + label + '.csv' # eEGT_EEF_'
# results_node_file = result_dir + now + '_' + str(test_size) + label + '_node.csv'
# summary_file      = result_dir + now + '_' + str(test_size) + label + '_summary.csv'
# decontaminated    = results_file[:-4] + '_decon.csv'

# results_nested_file      = result_dir + now + '_' + str(test_size) + label + '_nested.csv'
# results_node_nested_file = result_dir + now + '_' + str(test_size) + label + '_node_nested.csv'
# decontaminated_nested    = results_nested_file[:-4] + '_decon.csv'


# Script produces several output files; starting 'v6.03' - they are all listed here (for nested only version)
results_file                              = result_dir + now + '_' + str(test_size) + label + '.csv' # 1. id_label.csv - all trees meeting BS condition and sisterhood relationship eugl+algae, incl. mix; without 'node #'
results_nested_file                       = results_file[:-4] + '_nested.csv'                        # 2. id_label_nested.csv - only nested trees with nested relationship eugl inside algae, incl. mix
results_nested_unique_file                = results_nested_file[:-4] + '_unique_rows.csv'            # 3. without duplicite lines
results_nested_unique_decontaminated_file = results_nested_unique_file[:-4] + '_decon.csv'           # 4. decontaminated (using external files), without mix; but multiple occurence!!!! with mix
summary_nested_unique_decontaminated_file = result_dir + settings + '_' + now + '_' + str(test_size) + '.csv' # 'summary_nested_' + str(thread) + '.csv'                                     # 5. summary

results_unique_file                = results_file[:-4] + '_unique_rows.csv'            # 3. without duplicite lines
results_unique_decontaminated_file = results_unique_file[:-4] + '_decon.csv'           # 4. decontaminated (using external files), without mix; but multiple occurence!!!! with mix
summary_unique_decontaminated_file = result_dir + settings + '_' + now + '_' + str(test_size) + '.csv' # 'summary_sister_' + str(thread) + '.csv'                                     # 5. summary

results_node_file                         = results_file[:-4] + '_node.csv'
results_node_nested_file                  = results_nested_file[:-4] + '_node.csv'
r_input_file = summary_nested_unique_decontaminated_file
# r_input_file = result_dir + now + '_' + str(test_size) + label + '_unique_rows_summary.csv'
r_script_name = 'meegt_' + version + '.r'
#eliminated_mix_file                       = 

#decontaminated_nested    = results_nested_file[:-4] + '_decon.csv'
#results_node_nested_unique_file         = results_nested_unique_file[:-4] + '_node_nested.csv'
#results_node_nested_file                = results_nested_unique_file[:-4] + '_node_nested.csv'


def tree2lst(tree):
    '''
    input: tree in newick format as file or string
    output: list of all taxa in the same order as given in the tree
    '''
    # t = ''
    try:
        with open(tree) as f:
            t = f.read().strip('\n')
            # print('''reading file...''')
    except FileNotFoundError:
        t = tree
        # print('''cannot find file: t = {}'''.format(t))
    else:
        # print('''what???''')
        pass
    finally:
        t = t.replace('(', '')
        t = t.replace(')', '')
        t = t.replace(';', '')
        t = t.replace(':', '')
        return t.split(',')

def first_nonzero_element_index(lst):
    # print(lst)
    return [i for i in range(len(lst)) if lst[i] != 0][0]

full_target = tree2lst(full_Euglenoidea_tree)
target = tree2lst(target_tree)
# ignored_taxa = list(set(full_target) - set(target))
ignored_taxa = [] # ['DErc']
    # should contain only target for now (v6.02)
# print('''full_target:{}\ntarget: {}\nignored_taxa:{}'''.format(full_target, target, ignored_taxa))

def create_folders(path):
    '''
    changed order in 'v7.10'
    '''
    # global target
    global dic_algae
    # green:
    viridiplantae = ['ACCR', 'ACCH', 'ACMP', 'ACNP', 'ACOL', 'ACOT', 'ACPC', 'ACPO', 'ACPP', 'ACTA', 'ACVC', 'ASAT', 'ASGM', 'ASPP', 'ASPT', 'ASSM', 'ASSB'] # 17/17
    chlorarachniophyta = ['RCBN', 'RCGY', 'RCCR', 'RCLA', 'RCLG'] # 5/5
    
    
    # red:
    ochrophyta  = ['KOBP', 'KOCC', 'KOFJ', 'KOCS', 'KONI', 'KOOC', 'KOPC', 'KOPT', 'KOTP'] # 9/9; KO..
    heterokonta = ochrophyta + ['klap', 'kpcr', 'kppi'] # + Aplanochytrium, Cafeteria, Phytophthora
    cryptophyta = ['PCCP', 'PCGT', 'PCHP', 'PCHA', 'PCRA'] # 5/5; PC.
    cryptomonads = cryptophyta + ['pggo', 'pkrt'] # + Goniomonas + Roombia
    haptophyta  = ['HPCL', 'HPEH', 'HPCR', 'HPCP', 'HPPH', 'HPPC', 'HPPP'] # 7/7; HP..; no heterotrophic relative included

    rhodophyta  = ['ARCM', 'ARGS', 'ARCC', 'ARPP', 'ARPY', 'ARRC', 'ARRM'] # 7/7; AR..
    chromerida  = ['VCCV', 'VCVB'] # new
    dinophyta   = ['VDSM'] # new
    
    glaucophyta = ['AGCP', 'AGGW'] # AG..
    paulinella = ['RCPC'] # RCPC
    rhizaria = chlorarachniophyta + ['rfso'] + paulinella

    # Euglenoidea = ['DEEZ', 'DEEG', 'DEEH', 'DEEE', 'DERV', 'DEpt', 'DErc', 'DEnp', 'DEpv'] # DE..
    Heterolobosea = ['dhnf', 'dhng', 'dhpc'] # dh..
    Kinetoplastea = ['dkbs', 'dkem', 'dklp', 'dknd', 'dktb', 'dktc', 'dktg'] # 7/7; dk..
    discoba = Heterolobosea + Kinetoplastea
    Apicomplexa   = ['VABb', 'VACh', 'VACm', 'VACp', 'VALa', 'VANc', 'VAPb', 'VAPc', 'VAPy', 'VATp', 'VATg', 'Vggn'] # 12/12; Vggn = Gragerina niphandrodes
    ciliophora    = ['vceh', 'vcot', 'vcpt', 'vctt'] # 4/4
    Alveolata     = Apicomplexa + chromerida + dinophyta + ciliophora + ['vppm'] # + Perkinsus
    Metazoa       = ['omce', 'omeg', 'omhm', 'ommm', 'omnv', 'omsj', 'omtc'] # 7/7
    Amoebozoa     = ['naed', 'naeh', 'ndpe', 'ndva', 'nmdd', 'nmdp', 'nvfn'] # 7/7
    Fungi         = ['famg', 'fand', 'fasc', 'fats', 'fatm', 'faur', 'fava', 'fbcc', 'fblb', 'fbpp', 'fmec', 'fmeb', 'fmnc'] # 13/13
    Choanozoa     = ['ocmb', 'ocsa'] # 2/2
    Metamonada    = ['mfgl', 'mfss', 'mptv', 'mptf', 'mome', 'mppp'] # 6/6

    hbacteria = ['bzsu', 'baao', 'baas', 'bacp', 'baks', 'baml', 'bamp', 'bpci', 'bphp', 'bppb', 'bprl', 'bprb', 'bprd', 'bptp', 'bptx', 'bbbf', 'bbca', 'bbcl', 'bbdf', 'bbft', 'bbmp', 'bbpp', 'bfal', 'bfam', 'bfao', 'bflb', 'bfoi', 'bufn', 'bgcj', 'bgea', 'bgmm', 'bgsb', 'bgta', 'bgvc', 'bmct', 'bmpa', 'bmwc', 'bxrc', 'bsbp', 'bsli', 'bssa', 'bssc', 'bsss', 'bsst', 'btmc', 'btml', 'btmm', 'bvcf', 'bvvs']
    cyanophyta = ['BYAV', 'BYCW', 'BYFI', 'BYGL', 'BYLY', 'BYSE', 'BYSY']
    bacteria = hbacteria + cyanophyta
    archaea = ['ecmc', 'ecms', 'ecmy', 'eeaf', 'eeap', 'eeav', 'eehs', 'eemp', 'eepa', 'eepf', 'eeph', 'etnk', 'etnl'] # 13/13
    # target = Euglenoidea # + Kinetoplastea + Heterolobosea
    dic_algae = collections.OrderedDict()
    
    if bins == 3: # sort into archaea and bacteria or hbact and cyanophyta; 2 runs
        dic_algae['archaea'] = archaea
        dic_algae['bacteria'] = bacteria
    else:
        if bins == 5: # primary colours only; single run
            dic_algae['green'] = viridiplantae + chlorarachniophyta
            dic_algae['red'] = ochrophyta + cryptophyta + haptophyta + rhodophyta + chromerida + dinophyta
            dic_algae['glaucophyta'] = glaucophyta
            dic_algae['paulinella'] = paulinella
        elif bins == 12:
            dic_algae['viridiplantae'] = viridiplantae
            dic_algae['chlorarachniophyta'] = chlorarachniophyta
            dic_algae['ochrophyta'] = ochrophyta
            dic_algae['cryptophyta'] = cryptophyta
            dic_algae['haptophyta'] = haptophyta
            dic_algae['rhodophyta'] = rhodophyta
            dic_algae['chromerida'] = chromerida
            dic_algae['dinophyta'] = dinophyta
            dic_algae['apicomplexa'] = Apicomplexa
            dic_algae['glaucophyta'] = glaucophyta
            dic_algae['paulinella'] = paulinella
            dic_algae['cyanophyta'] = cyanophyta
        if bins == 11:
            dic_algae['heterokonta'] = heterokonta
            dic_algae['cryptomonads'] = cryptomonads
            dic_algae['rhizaria'] = rhizaria
            dic_algae['discoba'] = discoba
            dic_algae['alveolata'] = Alveolata
            dic_algae['metazoa'] = Metazoa
            dic_algae['amoebozoa'] = Amoebozoa
            dic_algae['fungi'] = Fungi
            dic_algae['choanozoa'] = Choanozoa
            dic_algae['metamonada'] = Metamonada
    dic_algae['mix'] = []
    #print(dic_algae['red'])
    #input()
    # dic_algae['cyanophyta'] = ['BYAV', 'BYCW', 'BYFI', 'BYGL', 'BYLY', 'BYSE', 'BYSY'] # BY..
        # unfortunately, there is one error: KOCR is Cafeteria - heterotroph
    # dic_algae['cryptophyta'] = ['PCCP', 'PCGT', 'PCHP', 'PCHA', 'PCRA', 'pggo', 'pkrt'] # including Roombia and Goniomonas
    # dic_algae['apicomplexa'] = ['VABb', 'VACh', 'VACm', 'VACp', 'VALa', 'VANc', 'VAPb', 'VAPc', 'VAPy', 'VATp', 'VATg', 'Vggn'] # Va.
    # dic_algae['alveolata'] = ['VABb', 'VACh', 'VACm', 'VACp', 'VALa', 'VANc', 'VAPb', 'VAPc', 'VAPy', 'VATp', 'VATg', 'vceh', 'vcot', 'vcpt', 'vctt', 'VDSM', 'Vggn', 'VCCV', 'VCVB', 'vppm'] # VA../v..
    #dic_algae['chlorophyta'] = ['ACCR', 'ACCH', 'ACMP', 'ACNP', 'ACOL', 'ACOT', 'ACPC', 'ACPO', 'ACPP', 'ACTA', 'ACVC'] # AC
    #dic_algae['streptophyta'] = ['ASAT', 'ASGM', 'ASPP', 'ASPT', 'ASSM', 'ASSB'] # AS..
    # dic_algae['kinetoplastea'] = Kinetoplastea
    # dic_algae['heterolobosea'] = Heterolobosea
    # dic_algae['metazoa'] = ['omce', 'omeg', 'omhm', 'ommm', 'omnv', 'omsj', 'omtc']
    # dic_algae['amoebozoa'] = ['naed', 'naeh', 'ndpe', 'ndva', 'nmdd', 'nmdp', 'nvfn'] 
    # dic_algae['fungi'] = ['famg', 'fand', 'fasc', 'fats', 'fatm', 'faur', 'fava', 'fbcc', 'fblb', 'fbpp', 'fmec', 'fmeb', 'fmnc']
    # dic_algae['choanozoa'] = ['ocmb', 'ocsa']
    # dic_algae['excavata'] = ['ehn', 'eha', 'ehw','ehf', 'ept','epg'] # , 'ekt', 'eke', 'ekb', 'ekn', 'ekp', 'ekg', 'ekc'] #  'epg', 'ept', 
    # dic_algae['metamonada'] = ['mfgl', 'mfss', 'mptv', 'mptf', 'mome', 'mppp'] # Metamonada
    # dic_algae['prokaryota'] = ['ecmc', 'ecms', 'ecmy', 'eeaf', 'eeap', 'eeav', 'eehs', 'eemp', 'eepa', 'eepf', 'eeph', 'etnk', 'etnl', 'bzsu', 'baao', 'baas', 'bacp', 'baks', 'baml', 'bamp', 'bpci', 'bphp', 'bppb', 'bprl', 'bprb', 'bprd', 'bptp', 'bptx', 'bbbf', 'bbca', 'bbcl', 'bbdf', 'bbft', 'bbmp', 'bbpp', 'bfal', 'bfam', 'bfao', 'bflb', 'bfoi', 'bufn', 'bgcj', 'bgea', 'bgmm', 'bgsb', 'bgta', 'bgvc', 'bmct', 'bmpa', 'bmwc', 'bxrc', 'bsbp', 'bsli', 'bssa', 'bssc', 'bsss', 'bsst', 'btmc', 'btml', 'btmm', 'bvcf', 'bvvs']
        # contamination of RCo (=DErc) with 'bgta'
        # Bacteria	Gammaproteobacteria	Tolumonas	auensis	bgt	bgta
    # dic_algae['prokaryota'] = ['bgta'] # checking for contamination
    # print(dic_algae)
    if not os.path.isdir(path):
        os.mkdir(path)
    for key in dic_algae:
        # print(key)
        folder = path + key
        if not os.path.isdir(folder):
            os.mkdir(folder)
    folder = path + 'mix'
    if not os.path.isdir(folder):
        os.mkdir(folder)
    if not os.path.isdir(log_path):
        os.mkdir(log_path)
    if not os.path.isdir(path_selected):
        os.mkdir(path_selected)
    if not os.path.isdir(nested_trees):
        os.mkdir(nested_trees)
    for key in dic_algae:
        # print(key)
        folder = path_selected + key
        if not os.path.isdir(folder):
            os.mkdir(folder)
    folder = path_selected + 'mix'
    if not os.path.isdir(folder):
        os.mkdir(folder)
    if not os.path.isdir(log_path_selected):
        os.mkdir(log_path_selected)

def m_array_index(array, searchItem):
    '''source:
    http://stackoverflow.com/questions/6518291/using-index-on-multidimensional-lists
    '''
    for i,x in enumerate(array):
        for j,y in enumerate(x):
            if y == searchItem:
                return i,j
    return -1,-1 # not found

def csv2matrix(filename):
    '''
    source: https://gist.github.com/dettmering/3767366
    256101181532 not used in aar
    '''
    ifile = open(filename, "rU")
    reader = csv.reader(ifile, delimiter=";")
    rownum = 0
    a = []
    for row in reader:
        a.append(row)
        rownum += 1
    ifile.close()
    return a

def bs_condition(b):
    '''
    What is the bootstrap treshold?
    '''
    if b >= bs_treshold:
        return True
    else:
        return False

def lst2str(lst):
    '''
    Convert list of int into one simple string;
    i. e. [1,2,3] -> '1.2.3'
    '''
    # https://stackoverflow.com/questions/5618878/how-to-convert-list-to-string
    return '.'.join(str(e) for e in lst)

class TreeError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Newick(str):
    # The same conditions as in single-euglenid analysis, i.e. cyano, apico excluded, strepto+chlorophyta = viridiplantae; ehf excluded
    add_treshold = 0  # 0 -> not used; how many foreign taxa can be in monophyletic clade - absolute value
    mult_treshold = 1 # 1 -> not used; how many foreign taxa can be in monophyletic clade - fraction, i.e. 0<x<1
    epsilon = 0
    prev_organism = ''
    prev_file = ''
    prev_algal_group = ''
    max_bs = 100
    b_one_alga_only = 102
    pure_tree_bs = 111
    len_max_bs = len(str(max_bs))

    def __init__(self, tree_string):
        self.tree_string = tree_string

    def find_all_targets(self,replacement=True):
        '''
        Return list of all (usually replaced) target taxa.
        Exception! It finds also 'bptx004183' as x004
        '''
        if replacement:
            # Xnnn contains always 1 target organism
            # Qnnn contains always more than 1 target organism - have to be counted!
            Qnnn = r'Q\d{3}'
            Xnnn = r'X\d{3}'
            QXnnn = r'[QX]\d{3}'
            # q = re.findall(QXnnn, self, re.IGNORECASE) # why ignorecase? The replacement is always with CAPITAL X or Q
            q = re.findall(QXnnn, self)
            # print(q)
            # if len(q) == 3:
                # print(self)
        else:
            pass
        return q

    def count_target(self, replacement=False):
        '''
        Returns number of targets in the (sub)tree.
        replacement = False - counting EEExxxxxx
        replacement = True - counting Qxxx
        '''
        nTarget = 0
        if replacement:
            nTarget = len(self.find_all_targets(replacement=True))
            # if nTarget == 2:
                # print('''count_target: {nTarget}'''.format(**locals()))
        else:
            for t in target:
                nTarget += self.count(t)
        return nTarget

    def count_ignored(self):
        '''
        count all taxa included in 'ignored_taxa' list
        count taxa, not sequences!
        '''
        nIgnored = 0
        for taxon in ignored_taxa:
            nIgnored += self.count(taxon)
        return nIgnored

    def analyse_target_composition(self):
        '''
        Return list representing presented target taxa in the group.
        '''
        target_frequency = [0] * len(target)
        for index, item in enumerate(target):
            #print('''subtree = {self}'''.format(**locals()))
            target_frequency[index] = self.count(item)
            #print(target_frequency)
        return target_frequency

    def count_taxa(self):
        '''
        Returns number of taxa in the (sub)tree. But it is a little cheating - counting only ','
        Now it seems that counting of actual taxa is neccessary?
        '''
        if False:
            return 1+self.count(',') # very dangerous
        else:
            '''
            taxon id:
                AAAnnnnnn  -> [A-Z]{3}\d{6,7}
                AAAnnnnnnn -> 
                Qnnn       -> Q\d{3}
            '''
            # AAAnnnnnn = r'[A-Z]{3}\d{6,7}' # Qnnn = r'Q\d{3}' 3-letter codes
            AAAAnnnnnn = r'[A-Z]{4}\d{6}' # 4-letter codes in STEPA40db
            
            k = re.findall(AAAAnnnnnn, self, re.IGNORECASE)
            q = self.count_target(replacement=True)
            #print(k)
            #print(q)
            return len(k) + q

    def count_algae(self):
        '''
        Count algae in the tree (each algal group as one item in list).
        key ~ phylum
        alga ~ taxon
        '''
        
        phyla_frequency = [0] * len(dic_algae)
        # contamination = dict()
        #phyla_frequency_cnt = [0] * len(dic_algae)
        t_0 = timer()
        for index, key in enumerate(dic_algae):
            # contamination[key] = set()
            # name = 'algae_' + phylum
            # print(name)
            for alga in dic_algae[key]:
                #print(alga)
                """ do you want do deal with contamination? It takes quite a long time :-(
                cnt_file = contamination_dir + 'algae/' + alga + '.cnt'
                #print(cnt_file)
                if os.path.isfile(cnt_file):
                    # print('''opening {}'''.format(cnt_file))
                    with open(cnt_file) as f:
                        for line in f:
                            row = line.strip()
                            if row in self:
                                # print('''cnt row: {}'''.format(row))
                                contamination[key].add(row)
                #phyla_frequency_cnt[index] += self.count(alga)
                # print('''{}: {} - {}'''.format(key, self.count(alga), len(contamination[key])))
                #"""
                n_cont = self.count(alga + '9')
                if n_cont != 0:
                    print('''The {} of {} contaminations not counted.'''.format(n_cont, alga))
                    n_alga = self.count(alga) - self.count(alga + '9')
                else:
                    n_alga = self.count(alga)
                phyla_frequency[index] += self.count(alga)
            # phyla_frequency[index] -= len(contamination[key])
        # print('''contamination: {}'''.format(contamination))
        # print(phyla_frequency_cnt)
        # print(phyla_frequency)
        #print('''decon takes {} s.'''.format(timer() - t_0))
        return phyla_frequency

    def count_brackets(self, type='left'):
        if type == 'left':
            return self.count('(')
        elif type == 'right':
            return self.count(')')

    def is_newick(self):
        '''
        '''
        # print(self)
        if len(self) > 0:
            if self[-1] != ';':
                return False
            else:
                left_bracket = self.count_brackets() #self.count('(') # number of open bracket :-)
                right_bracket = self.count_brackets(type='right') #self.count(')') # number of close bracket
                if left_bracket != right_bracket:
                    return False
                else:
                    return True
        else:
            return False

    def topology(self, nodes=False):
        '''
        input: string in Newick format, i. e. ((a,b),(c,(d,e)));
        output: two lists with indexes of left and right brackets
            the index of first char in the tree is 1!! not 0!!
        '''
        nBracket = self.count('(')
        # print(nBracket)
        nTaxa = 1 + self.count(',') # Number of Taxa in self; each taxon, except the last one, ends with ","
                                    # Does the number of "," correspond number of branches? - Yes.
        #"""
        try:
            nOpen = [0 for x in range(nBracket)]  # indexes of open brackets
            mClose = [0 for x in range(nBracket)] # indexes of close brackets
            nClose = [0 for x in range(nBracket)] # indexes of close brackets in right order
            mNode = [0] * (nTaxa - 1) # get indexes of ',' weird trees (A,(B,C), D)
            nNode = [0] * (nTaxa - 1) # transform to right order
            # print(mNode)
            n = 0
            m = 0
            i = 0
            o = 0
            for char in self: # index of '(' and ')'
                i += 1
                if char == '(':
                    nOpen[n] = i
                    n += 1
                elif char == ')':
                    mClose[m] = i
                    m += 1
                elif char==',':
                    # print(o)
                    mNode[o] = i
                    o += 1
            for n in range(nBracket-1, -1, -1): # indexes of ')' in right order
                for m in range(0,nBracket):
                    if nOpen[n] < mClose[m]:
                        nClose[n] = mClose[m]
                        mClose[m] = -1
                        break
            if nodes:
                for n in range(nBracket-1, -1, -1): # backward for cycle of brackets
                    # print('''n = {}: {}, {}'''.format(n, nOpen[n], nClose[n]))
                    inside = []
                    for node in mNode: # for cycle of ','
                        if nOpen[n] < node < nClose[n]:
                            inside.append(node)
                    # print('''   inside: {}'''.format(inside))
                    if len(inside) == 1:
                        nNode[n] = inside[0]
                    else:
                        err = 0
                        for node in inside:
                            if node not in nNode:
                                err += 1
                                nNode[n] = node
                        if err > 1:
                            pass
                            # raise('''More than 1 most outer node present''')
                return nOpen, nClose, nNode
            else:
                return nOpen, nClose #, None
        except:
            raise TreeError('''Topology not determinable''')
       #"""

    def get_bootstraps(self, right_brackets):
        '''
        input: tree string + list of position of right bracket
        output: list of bootstraps
        '''
        nBracket = self.count('(')
        bootstraps = [0 for x in range(nBracket)]
        # print('''get_bootstraps():\n''')
        # print('''    bootstraps: {}'''.format(bootstraps))
        # print('''    nBracket: {}'''.format(nBracket))
        # print('''    right_brackets: {}'''.format(right_brackets))
        j = -1
        for bracket in right_brackets:
            # print(bracket)
            j += 1
            i = self.find(':', bracket)
            # print(i)
            if i < bracket and self[-2] == ')': # j == 0: 
                # Have to work also for subtrees!!!
                # i < bracket only in full tree!!! in the one case (i.e. the bootstrap for whole tree)
                # self[-2] != ')' i.e. additional b = 101 added in analyse_root_epsilon()
                bootstraps[j] = Newick.pure_tree_bs
            else:
                try:
                    bootstraps[j] = int(self[bracket:i])
                    # print('''    bootstraps: {}'''.format(bootstraps))
                except ValueError:
                    # print(self[bracket:i])
                    print('''The tree is corupted.''')
        return bootstraps

    # def append_bootstrap(self,right): # ??? in get_all_sisterhoods()
        # '''
        # Append bootstrap of given branch.
        # '''
        # self + str(b[right.index(in_sister_branch_right)])

    def get_id(self, taxon, start = 0):
        '''
        Get id of protein (of taxon) in the tree, i.e. the position of the last char (of the id) of the first occurrence in the tree.
        Based on position of next colon ':'.
        Applicable only for trees with branch length.
        '''
        taxon_index = self.find(taxon, start)
        colon_index = self.find(':', taxon_index)
        if colon_index == -1:
            raise TreeError('''No branch tree!''')
        else:
            # print(self[taxon_index:colon_index])
            return self[taxon_index:colon_index], colon_index

    def get_target_ids(self):
        '''
        '''
        target_ids_lst = [''] * self.count_target(replacement=False)
        i = 0
        for t in target:
            start = 0
            max_target = self.count(t)
            # print('''taxon: {} - number: {}'''.format(t, max_target))
            for j in range(max_target):
                # print('''j: {}'''.format(j))
                target_ids_lst[i], start = self.get_id(t, start)
                i += 1
        return target_ids_lst

    def branch2tree(self, r, inner_branch):
        '''
        Adds appendices to the branch, so if replaced, the rest is regular tree.
        '''
        #print('''{:x^50}'''.format(' branch2tree() '))
        #print('''branch2tree()''')
        #print('''{:*^50}'''.format(' INPUT '))
        #print('''self =\n{self}\nr= {r}\ninner_branch = {inner_branch}'''.format(**locals()))
        colon = self.find(':', r)
        comma = self.find(',', r)
        bracket = self.find(')',r)
        i = min(bracket,comma)
        length_appendices = self[r:i]
        #print('''{:+^50}'''.format(' intermezzo '))
        #print('''i = min(comma = {comma}, bracket = {bracket}) = {i}; colon = {colon}; length_appendices = {length_appendices}'''.format(**locals()))
        #print('''{:=^50}'''.format(' OUTPUT '))
        return inner_branch + length_appendices

    def tiniest_target_group_determination(self): # Blind way? :-( but cute
        '''
        '''
        target_groups = []
        left, right = self.topology()
        for t in self.get_target_ids():
            # tindex for index in newick tree string
            # lindex for index in list of brackets, e.g. left.index(l_tindex_l)
            t_tindex = self.find(t)
            r_tindex_r = self.find(')', t_tindex) + 1
            l_tindex_r = left[right.index(r_tindex_r)]
            l_tindex_l = self.rfind('(', 0, t_tindex) + 1
            r_tindex_l = right[left.index(l_tindex_l)]
            l_tindex = min(l_tindex_l,l_tindex_r)
            r_tindex = max(r_tindex_l,r_tindex_r)
            # print('''target: {} - left: {}   right: {}'''.format(t, l_tindex, r_tindex))
            tiniest_subtree = Newick(self[l_tindex:r_tindex])
            if t not in tiniest_subtree:
                raise TreeError('''The target is not in subtree!''')
            nTarget = tiniest_subtree.count_target(replacement = False)
            nTaxa = tiniest_subtree.count_taxa()
            # print('''Number of target / taxa in the tree: {} / {}'''.format(tiniest_subtree.count_target(), tiniest_subtree.count_taxa()))
        return None

    def analyse_root_epsilon(self, left = False, right = False, b = False):
        '''
        The target organisms can be on the 'base' of the tree in newick format.
        The function is looking for such sitation and creates one target group as big as possible out of them.
        Tree (E,(A,E),E) is not solved properly -> two target groups instead of one
        '''
        #print('''analyse_root_epsilon()''')
        length = 0.1 # I do not care the branch length!
        nTaxa = self.count_taxa()
        nTarget = self.count_target(replacement=False)
        if nTaxa - nTarget == 1:
            ids = self.get_target_ids()
            comp = self.analyse_target_composition()
            AAAAnnnnnn = r'[A-Z]{4}\d{6}'
            ids_all = re.findall(AAAAnnnnnn, self, re.IGNORECASE)
            #print(ids_all)
            for id in ids_all:
                print('''id[:{}]: {}'''.format(code_length, id[:code_length]))
                if id[:code_length] not in target: # There is always one which is not in target!
                    # print('''id {id} not in target!'''.format(**locals))
                    code = 'Q000'
                    msg = code + ': ' + str(comp) + ' --- ' + str(ids) + '\n'
                    tout.write(msg)
                    # print(msg)
                    return Newick('(Q000,' + id + ')' + str(Newick.b_one_alga_only)  + ';'), comp, code, True
        else:
            if not (left and right and b):
                left, right = self.topology()
                b = self.get_bootstraps(right)
            # print(right)
            for i in reversed(range(len(right))):
                # print('''i = {i}'''.format(**locals()))
                inner_branch = Newick(self[left[i]-1:right[i]])
                inner = self.branch2tree(right[i], inner_branch)
                #print(inner)
                outer_branch = Newick(self.replace(inner, 'X'))
                if outer_branch == self:
                    raise TreeError('''The branch was not cutted!!!''')
                #print(outer_branch)
                oTarget = outer_branch.count_target(replacement=False)
                oTaxa = outer_branch.count_taxa()
                # print('''oTarget = {oTarget}, oTaxa = {oTaxa}'''.format(**locals()))
                if oTarget == oTaxa and oTaxa > 0:
                    if oTarget == 1:
                        letter = 'X'
                        print('''X000''')
                        # There should not be any and there is not any.
                    else:
                        letter = 'Q'
                    code = letter + '000'
                    # print('''outer branch: {outer_branch}\ninner branch = {inner}'''.format(**locals()))
                    b_outer = b[i]
                    ids = outer_branch.get_target_ids()
                    comp = outer_branch.analyse_target_composition()
                    msg = code + ': ' + str(comp) + ' --- ' + str(ids) + '\n'
                    tout.write(msg)
                    # print(comp)
                    #print('''b = {b_outer}'''.format(**locals()))
                    # print('..................................\n......,,,,,,,,,,,..............')
                    return Newick('(' + code + ',' + inner_branch + str(b_outer) + ':' + str(length) + ');'), comp, code, True
                    #break
            else:
                return self, None, None, False

    def target_group_determination(self):
        '''
        Finds all branches containing only target organisms, i.e. does not catch 'branches' with only one target (incl. special condition - in near future).
        Returns list of the branches and new tree, where the branches are replacet by CODE for easier searching.
        
        Returns four values:
            target_groups - list of ids of target group organisms, only for debbuging; as in 256007231139 excluding Q000 (since not needed)
            new_tree - tree with all target groups replaced
            q_list - list of replaced target groups, i.e. Q000, Q001, ..., Qxxx
            target_composition_list - list of i-element lists, where i is number of target organisms
        '''
        # print('''target_group_determination()\n''')
        global Q000_absence
        temp_group = ''
        target_groups = []
        number_target_groups = 0
        q_index = 0
        q_list = []
        target_composition_list = []
        left, right = self.topology()
        b = self.get_bootstraps(right)
        #print('''self = {self}\n    left: {left}\n    right: {right}\n    b = {b}'''.format(**locals()))
        self, composition, code, change = self.analyse_root_epsilon(left = left, right = right, b = b)
        if change: # it means that Q000 exists. i. e. Q000_absence must be set to 0
            Q000_absence = 0
            #print('''CHANGE: Q000_absence is set to {}'''.format(Q000_absence))
            left, right = self.topology()
            b = self.get_bootstraps(right)
            print('''q_list + code: {} + {}'''.format(q_list, code))
            q_list.append(code)
            # print('''q_list: {}'''.format(q_list))
            target_composition_list.append(composition)
            # print(target_composition_list)
            # print('''target_groups added: {}'''.format(target_groups))
            # target_groups.append(inner_branch.get_target_ids()) # original line, commented in aaq
            # target_groups.append(self.get_target_ids()) # I do not know why it is here
            # print('''target_groups added: {}'''.format(target_groups))
            number_target_groups +=1
        else:
            Q000_absence = 1
        new_tree = self
        # print('''    inner_branch\n{}'''.format(self))
        # print('''    left: {}\n    right: {}'''.format(left, right))
        # print('''    b = {}'''.format(b))
        #self.analyse_root_epsilon()
        for i in range(len(right)):
            #print(i)
            inner_branch = Newick(self[left[i]-1:right[i]])
            nTarget = inner_branch.count_target(replacement=False)
            nIgnored = inner_branch.count_ignored()
            nTaxa = inner_branch.count_taxa()
            # print('''nTarget = {nTarget}'''.format(**locals()))
            if nTarget+nIgnored == nTaxa: # Does not catch one-target branches! Add '2nd run'! Done, see abandoned.
                if inner_branch in temp_group:
                    pass
                    # print('''temp: {}\ntemp_group: {}'''.format(inner_branch, temp_group))
                else:
                    number_target_groups +=1
                    q_index += 1
                    temp_group = inner_branch
                    b_inner_branch = inner_branch + str(b[i])
                    ids = inner_branch.get_target_ids()
                    target_groups.append(ids)
                    #print(b_inner_branch)
                    q = 'Q' + str(q_index).zfill(self.len_max_bs)
                    new_tree = Newick(new_tree.replace(b_inner_branch, q))
                    q_list.append(q)
                    comp = inner_branch.analyse_target_composition()
                    target_composition_list.append(comp)
                    #print('''q = {q}: {comp} --- {ids}'''.format(**locals()))
                    msg = q + ': ' + str(comp) + ' --- ' + str(ids) + '\n'
                    tout.write(msg)
            else:
                pass
        abandoned = new_tree.get_target_ids()
        # print('''abandoned: {abandoned}'''.format(**locals()))
        a_index = -1
        for k in abandoned:
            number_target_groups += 1
            q_index += 1
            a_index += 1
            #q = 'Q' + str(q_index).zfill(self.len_max_bs)
            q = 'X' + str(q_index).zfill(self.len_max_bs)
            target_groups.append(k)
            new_tree = Newick(new_tree.replace(k, q))
            q_list.append(q)
            comp = Newick(k).analyse_target_composition()
            target_composition_list.append(comp)
            msg = q + ': ' + str(comp) + ' --- ' + str(abandoned[a_index]) + '\n'
            tout.write(msg)
        # print('''Number of target groups is {}.'''.format(number_target_groups))
        # print('''list of q: {q_list}'''.format(**locals()))
        # print('''Tree with replaced targets:\n{new_tree}'''.format(**locals()))
        return target_groups, new_tree, q_list, target_composition_list

    def analyse_sisterhood_composition(self):
        '''
        Analyse composition of given branch (sisterhood).
        If the branch contain taxa from very one phylum, returns the phylum, False.
        If the branch contain taxa from very one phylum + 1 additional algal taxa, returns phylum, True. (Does not make sense to me; the additional taxa should be any taxa)
        '''
        # print('''analyse_sisterhood_composition()''')
        nTaxa = self.count_taxa()
        nAlga = self.count_algae() # starting v5.34 does not contain MMETSP contamination
        nAlgae = sum(nAlga)
        nIgnored = self.count_ignored()
        nQ = self.count_target(replacement=True) # count groups of Q, not Qs themselves!!!!
        #print('''nTaxa = {nTaxa}, nQ = {nQ}, nAlgae = {nAlgae}: {nAlga}'''.format(**locals()))
        delta = (nTaxa - nQ - self.add_treshold - nAlgae - nIgnored) * self.mult_treshold
        #delta = nTaxa - nAlgae - nQ
        # print('''delta = {delta}'''.format(**locals()))
        if  delta <= self.epsilon and nAlgae > 0:
        #if  delta == 0 and nAlgae > 0: # if delta = 1 -> bs = 101 ??? for cycle????
            # print('''delta = {delta}: nTaxa = {nTaxa}, nQ = {nQ}, nAlgae = {nAlgae}: {nAlga}'''.format(**locals()))
            '''
            In near future change to 'relaxed' condition.
            Should return the branch 'color' if pure - working well.
            '''
            if nQ == 1:
                # print('''nQ = 1''')
                multi = False
                multi_target = False
                concatenate_target_group_id = self.find_all_targets(replacement=True)[0]
                # print('''concatenate_target_group_id: {concatenate_target_group_id}'''.format(**locals()))
            else:
                # print('''The concatenation of more monophyletic groups of Q is in progress!''')
                multi_target = self.find_all_targets()
                # print('''multitarget = {multi_target}!'''.format(**locals()))
                concatenation = []
                concatenate_target_group_id = lst2str(multi_target) # e.g. 'Q004.X005'
                # print('''concatenate_target_group_id: {concatenate_target_group_id}'''.format(**locals()))
                for i in multi_target:
                    j = int(i[1:])-1
                    # print('''j = {j}'''.format(**locals()))
                    # print('''l[j] = {}'''.format(l[j]))
                    if isinstance(l[j], list): # check if l[j] is a list; l is list of indeces of left brackets; see "for i, (l, r) in enumerate(zip(left, right)):"
                        concatenation += l[j]
                    else:
                        concatenation.append(l[j])
                    multi = concatenation
            for index, key in enumerate(dic_algae):
                #print('''Monoalga?\nindex = {index}, key = {key}'''.format(**locals()))
                if nAlga[index] == nAlgae:
                    #print('''The sistergroup of Q is one particular group of algae.''')
                    return key, multi, concatenate_target_group_id, False
                    break
            else: # else after for is executed only if the for cycle is not interrupted by break! Exactly what I need :-)
                #print('''MonoAlga!!! +-1''')
                for index, key in enumerate(dic_algae): # This for cycle does not do anything.
                    #print('''index = {index}, key = {key}'''.format(**locals()))
                    # if nAlga[index] > 1 and nAlga[index] + 1 == nAlgae:
                    if False: # nAlga[index] > 1 and nAlga[index] + nQ + 1 == nTaxa:
                        #print('''Find something!!! key = {key}'''.format(**locals()))
                        return key, multi, concatenate_target_group_id, True
                        break
                else:
                    #print('''PolyAlga!!!''')
                    return 'mix', multi, concatenate_target_group_id, False # 'mix'
                    # return False, multi, concatenate_target_group_id, False # 'mix'
        else: # print('''black-and-white''')
            return False, False, False, False

    def get_closest_sisterhoods(self, organism): # not used
        '''
        Find both sister branch of an organism.
        '''
        #print('''get_sisterhood()\n''')
        #print('''    subtree:\n        {}'''.format(self))
        left, right = self.topology()
        b = self.get_bootstraps(right)
        #print('''    left: {}\n    right: {}\n    b = {}'''.format(left, right, b))
        if self.find('(' + organism) > -1:
            # print('''Found FIRST''')
            in_sister_branch_left  = self.find('(', self.find(organism))
            in_sister_branch_right = right[left.index(in_sister_branch_left+1)]
            in_sister_branch = Newick(self[in_sister_branch_left:in_sister_branch_right])
            b_in_sister_branch =  in_sister_branch + str(b[right.index(in_sister_branch_right)])
            out_sister_branch = self.replace(self[self.find(organism)-1:in_sister_branch_right], 'QIN')
                # The result is not valid tree, i.e. 
            b_out_sister_branch = self.replace(b_in_sister_branch, 'QIN')
            # print('''in_sister_branch_left: {}'''.format(in_sister_branch_left))
            # print('''in_sister_branch_right: {}'''.format(in_sister_branch_right))
            # print('''in_sister_branch: {}'''.format(in_sister_branch))
            # print('''b_in_sister_branch: {}'''.format(b_in_sister_branch))
            # print('''out_sister_branch: {}'''.format(out_sister_branch))
            # print('''b_out_sister_branch: {}'''.format(b_out_sister_branch))
            # print('''left.index(in_sister_branch_left): {}'''.format(left.index(in_sister_branch_left+1)))
        elif self.find(',' + organism) > -1:
            pass
            #print('''found LAST''')
            ## in <-> out_sister_branch
            # pass
        # print(in_sister_branch)
        # print(out_sister_branch)
        return None

    def analyse_sisterhoods(self, organism):
        '''
        'Mixed-concatenation' of get_all_sisterhoods() and analyse_sisterhood_composition().
        Better? Easier?
        Since I do not care about best bootstrap support, but only if the bootstrap support is high enough , I keep the branch with bs high enough only once.
        Neccessary for last step of analysis
        
        # is q present OK
        # get b ok
        # count alga OK
        # count q OK
        # count taxa OK
        # is algal? OK
        # copy tree
        '''
        #print(self)
        global Q000_absence
        # print('''analyse_sisterhoods()''')
        left, right = self.topology()
        b = self.get_bootstraps(right)
        sorted_right = sorted(right)
        # print('''    left: {}\n    right: {}\n    b = {}'''.format(left, right, b))
        for i, (l, r) in enumerate(zip(left, right)):
            # run for indexes and values of 2 lists
            #print('''{0:-^25}'''.format('subtree'))
            #print('''l = {}, r = {}'''.format(l, r))
            inner_branch = Newick(self[l-1:r])# + str(b[i])) # trouble line
            if bs_condition(b[i]):
                if b[i] == Newick.pure_tree_bs:
                    outer_branch = ''
                else:
                    '''
                    in order to properly replace inner branch by empty string
                    needed for ',' counting - count taxa changed to count taxa not ',' in aag.py
                    '''
                    inner = self.branch2tree(r, inner_branch)
                    #print('''inner: {inner}'''.format(**locals()))
                    if False: #l-1 == left[i-1]: # OK
                        '''
                        In order to superproperly replace inner branch by empty string; looking for resting double bootstrap (successful) and keeping only the highest (unsuccessful)and summation of branch length (unsuccessful)
                        NOT USED
                        '''
                        pass
                        """
                        #if self[l-1] == '(': !!! does not work for unknown reason [256007191223] !!!
                        b_next = b[i-1]
                        #b_prev = b[i+1] # right[i-1] najdi nejvyssi nizsi cislo; experimentalne zjisteno
                        #print('''right[i-1] = {}'''.format(right[i-1]))
                        m = bisect.bisect_left(sorted_right, right[i-1])
                        n = sorted_right[m-1]
                        o = right.index(n)
                        b_prev_better = b[o]
                        #print('''sorted right = {}'''.format(sorted_right))
                        #print('''#####################################m = {}'''.format(m))
                        #print('''m = {}, right value = {}, right index = {}, b = {}'''.format(m, n, o, b_prev_better))
                        #bisect.bisect(a, x, lo=0, hi=len(a))
                        #print('''nextb = {}, prevb = {}; better prevb = {}'''.format(b_next, b_prev, b_prev_better))
                        text_l = '''l-1 = ''' + str(l-1)
                        text_r = '''r[l-1] = ''' + str(right[left.index(l-1)])
                        #print(text_l)
                        #texta = '''self[l-1] = ''' + str(self[l-1])
                        #print(texta)
                        #print(text_r)
                        #print('''left: {}\nright: {}'''.format(l-1, right[left.index(l-1)]))
                     # else:
                        # print('''{0:*<515}'''.format(1))
                     #print('''length_appendices: {}'''.format(length_appendices))
                     #print('''inner:\n {}'''.format(inner))
                     """
                    outer_branch = Newick(self.replace(inner, ''))
                    #if outer_branch == self:
                    #    raise TreeError('''The branch was not cutted!!!''')
                # print('''Subtree:\n    inner branch: {inner_branch}\n    outer branch: {outer_branch}'''.format(**locals()))
                # print('''Target group: {}'''.format(q))
                if q in inner_branch:
                    algal_group, multi, concatenate_target_group_id, b101 = inner_branch.analyse_sisterhood_composition()
                    # print('''q present in inner branch {}'''.format(inner_branch))
                    # print('''algal group: {algal_group}'''.format(**locals()))
                elif q in outer_branch:
                    algal_group, multi, concatenate_target_group_id, b101 = outer_branch.analyse_sisterhood_composition()
                    # print('''q present in outer branch {}'''.format(outer_branch))
                    # print('''algal group: {algal_group}'''.format(**locals()))
                #print('''Q000_absence (in analyze_sisterhood()) is {}.'''.format(Q000_absence))
                if algal_group:
                    # print('''algal group found ({})'''.format(algal_group))
                    organism = concatenate_target_group_id
                    if multi:
                        numstr_target_composition = Newick(multi).analyse_target_composition()
                        # print('''target_composition = {numstr_target_composition}'''.format(**locals()))
                        # print('''\n{:~^50}'''.format('~'))
                        # print('''{:~^50}'''.format('#'))
                        # print('''{:~^50}'''.format('~'))
                        # print(multi)
                    else:
                        #print('''target_composition: {}\n  organism: {}'''.format(target_composition, organism))
                        #print('''Q000 absence is {}'''.format(Q000_absence))
                        # print('''index is {}'''.format(int(organism[1:])-1*Q000_absence))
                        numstr_target_composition = target_composition[int(organism[1:])-1*Q000_absence]
                        """
                        the original line:
                        numstr_target_composition = target_composition[int(organism[1:])-1]
                        why the hell was there '-1'?????
                        There are two types of trees: with Q000 and without Q000 (starting from 1)
                        """
                        # print('''numstr_target_composition: {}'''.format(numstr_target_composition))
                        # cause problems with active ignored_taxa
                            # e.g. ValueError: invalid literal for int() with base 10: '001.Q002.Q003.Q004'
                    src = file
                    #print('src:  {src}'.format(**locals())) # prints 'path to file'
                    # print('file: {file}'.format(**locals())) # returns key error
                    #print('file: {}'.format(file)) # prints 'path to file'
                    bootstrap = b[i]
                    if b101:
                        # print('101')
                        bootstrap = 101
                    id_tree = file[file.rfind('\\')+1:file.rfind('.')]
                    print(id_tree)
                    # dst = path + algal_group + '/' + str(bootstrap).zfill(3) + '-' + file[-10:-4] + '-' + organism + '-' + lst2str(numstr_target_composition) + file[-4:]
                    # dst = path + algal_group + '/' + str(bootstrap).zfill(3) + '-' + id_tree + '-' + organism + '-' + lst2str(numstr_target_composition) + file[-4:]
                    dst = path + algal_group + '/' + id_tree + '-' + str(bootstrap).zfill(3) + '-' + organism + '-' + lst2str(numstr_target_composition) + file[-4:]
                    # print(dst)
                    if copy_tree:
                        shutil.copy(src, dst)
                    msg = dst + '\n'
                    tout.write(msg)
                    # write into csv #
                    print('''{organism} == {self.prev_organism}'''.format(**locals()))
                    try:
                        print('''{} == {}'''.format(file, self.prev_file))
                    except KeyError:
                        print('''no file specified''')
                    print('''{algal_group} == {self.prev_algal_group}'''.format(**locals()))
                    # important conditions for nested sorting
                    algal_group_is_same = algal_group == self.prev_algal_group
                    algal_group_is_mix = algal_group == 'mix'
                    prev_algal_group_is_mix = self.prev_algal_group == 'mix'
                    if organism == self.prev_organism and file == self.prev_file and (algal_group_is_same or algal_group_is_mix or prev_algal_group_is_mix): # originally: organism in self.prev_organism ...
                        # checking if the taxon is nested, i.e. it is the second check of sisterhood
                            # 1. check if the target containing branch is the same (important in MeEGT)
                            # 2. check if the tree is still the same
                            # 3. check if last two algal groups are the same OR
                            # 4. if one of them is 'mix' (the final group is then 'mix')
                        pass
                        # first condition is the main change in comparison to aaq.py
                        #print('''   Condition 1007 met: prev_organim - {} (organism = {}), prev_algal_group - {}; prev_file: {}.'''.format(self.prev_organism, organism, self.prev_algal_group, self.prev_file))
                        if False:#'.' in self.prev_organism:
                            pass
                            # print('''Target {} is not monophyletic, excluding'''.format(self.prev_organism))
                        else:
                            print('''writing into nested result file''')
                            #'''
                            if prev_algal_group_is_mix or algal_group_is_mix: # label: 2mix
                                algal_group2 = 'mix'
                            else:
                                algal_group2 = algal_group
                            # Interstingly, the above if clausule does not produce correct result, but without it, the sorting seems to be correct.
                            #'''
                            nested_csv.writerow([algal_group2] + [id_tree] + [tree_size] + [organism] + [str(bootstrap).zfill(3)] + numstr_target_composition)
                            #nested_csv.writerow([algal_group] + [id_tree] + [tree_size] + [organism] + [str(bootstrap).zfill(3)] + numstr_target_composition)
                            dst = nested_trees + id_tree + '.dip'
                            if copy_tree:
                                shutil.copy(src, dst)
                    # else: # writing into the result file
                    print('''writing all into result file''')
                    # id_tree = file[file.rfind('\\')+1:file.rfind('.')]
                    # print(id_tree)
                    results_csv.writerow([algal_group] + [id_tree] + [tree_size] + [organism] + [str(bootstrap).zfill(3)] + numstr_target_composition)
                    self.prev_organism = organism
                    self.prev_file = file
                    self.prev_algal_group = algal_group
                    # print(dst)
                else:
                    pass
                    # print('''No algal group found''')
        return None

    def reverse_target_counting_in_subtree(self):
        '''
        Looks into the subtree and counts all target and compare it with the initial target counting.
        I would like to do:
            1. input: all subtrees with sisterhoods of euglenids and algae.
            2. do: look into all of them and count targets
            3. do: decide if the target is the same
            4. do: if not concatenate targets
            '''
        pass
        return None

def is_enter(left_branch, right_branch): # of target tree in file 'target_tree'
    '''
    input: arrays of presence
    trying to find the first occurence of 'gene' in target tree
    'basal' is a single taxon - the content is TRUE or FALSE
    terminal is any number of taxa (tuple)
    basal AND [OR(terminal)] = b AND (t1 OR t2 OR t3 OR ... OR tn)
    upgrade: basal is also any number of taxa
    
    '''
    left = 0
    right = 0
    for i in left_branch:
        left = left or i
    for i in right_branch:
        right = right or i
    return left and right

def find_enter(presence_of_gene):
    '''
    input: array of leaves of target tree: 0 - absence of protein, 1 - presence of protein; target tree
    output: id of node of enter of the gene
    new: recursively looks in all branches and returns the index of most ancient single entrance of a gene
    
    every ',' is a node
    
    '''
    # print('''{:#^50}'''.format('#'))
    # print('''{:#^50}'''.format('#'))
    # print('''{:#^50}'''.format('#'))
    # print('''{:#^50}'''.format('#'))
    # print('''{:#^50}'''.format('#'))
    # print('''{:#^50}'''.format('#'))
    # print('''find_enter()''')
    left, right, commas = Newick.topology(target_tree, True)
    tree_list = tree2lst(target_tree)
    # nodes = [0] * len(tree_list) # terminal nodes, i.e. taxa; left could have one element less!!
    # print('''number of taxa: {}; tree list: {}'''.format(len(tree_list), tree_list))
    # print('''nodes: {}'''.format(nodes))
    # print('''len(target_tree) = {}'''.format(len(target_tree)))
    entrance = [0] * len(commas) # every comma is an internal node, could be an enter of gene
    # print('''left: {}, right: {}, commas: {}; len(commas) = {}, len(left) = {}'''.format(left, right, commas, len(commas), len(left)))
    # print(tree)
    # print('''presence_of_gene {}'''.format(presence_of_gene))
    if len(presence_of_gene) != len(tree_list):
        raise ValueError('''The number of taxa in presence list is different from the number of taxa in target tree.''')
    if len(presence_of_gene) - presence_of_gene.count(0) > 1: # how many non-zero
        # print('''The protein is presented in more than one taxon, i. e. the common ancestor of all the taxa obtained the protein.''')
        # print(range(len(entrance)))
        # for  i in range(len(entrance)): # nodes
            # print(i)
        for node in range(len(entrance)): # nodes
            # print('''node: {}'''.format(node))
            # print('''node: {}; left: '{}'; right: '{}' '''.format(commas[node], target_tree[commas[node]-2], target_tree[commas[node]]))
            # find composition of left branch
            if target_tree[commas[node]-2] == ')':
                i = 0
                while True:
                    try:
                        first_left_taxon = target_tree[left[node]-2+i:left[node]-2+4+i]
                    except IndexError:
                        print('''Be carefull, something is wrong.'''.format())
                        # print('''{:~^50}'''.format('#'))
                        # print('''{}\ni = {}, left = {}, node = {}'''.format(target_tree, i, left, node))
                    i += 1
                    if first_left_taxon in tree_list:
                        break
                i = 0
                while True:
                    last_left_taxon = target_tree[commas[node]-4-i:commas[node]-i]
                    i += 1
                    if last_left_taxon in tree_list:
                        break
                lpresence = presence_of_gene[tree_list.index(first_left_taxon):tree_list.index(last_left_taxon)+1]
                # print('''    first_left_taxon: {}; last_left_taxon: {}'''.format(first_left_taxon, last_left_taxon))
            else: # i.e. CODE,....
                lcode = target_tree[commas[node]-code_length-1:commas[node]-1]
                # print('''lcode(-1): {}'''.format(lcode))
                lindex = tree_list.index(lcode)
                lpresence = [presence_of_gene[lindex]]
                # print('''    lcode(lindex): '{}({})-{}'' '''.format(lcode, lindex, lpresence))
            # find composition of right branch
            if target_tree[commas[node]] == '(':
                
                i = 0
                while True:
                    first_right_taxon = target_tree[commas[node]+i:commas[node]+4+i]
                    i += 1
                    if first_right_taxon in tree_list:
                        break
                i = 0
                while True:
                    #last_right_taxon = target_tree[right[node-2]-4-i:right[node-2]-i]
                    last_right_taxon = target_tree[right[node]-4-i:right[node]-i]
                    i += 1
                    if last_right_taxon in tree_list:
                        break
                # print('''    first_right_taxon: {}; last_right_taxon: {}'''.format(first_right_taxon, last_right_taxon))
                rpresence = presence_of_gene[tree_list.index(first_right_taxon):tree_list.index(last_right_taxon)+1]                
            else:
                rcode = target_tree[commas[node]:commas[node]+code_length]
                rindex = tree_list.index(rcode)
                rpresence = [presence_of_gene[rindex]]
                # print('''    rcode(rindex): '{}({})-{}'' '''.format(rcode, rindex, rpresence))
            # print('''lcode(lindex): '{}({})-{}'; rcode(rindex): '{}({})-{}' '''.format(lcode, lindex, lpresence, rcode, rindex, rpresence))
            enter = is_enter(lpresence, rpresence)
            # print('''n = {}: l = {}, r = {}, c = {}; ENTER = {}'''.format(node, left[node], right[node], commas[node], enter))
            # print(enter)
            entrance[node] = enter
        # print(entrance)
        try:
            # first_enter_of_gene = entrance.index(1) # wrong: looking for first non-zero
            first_enter_of_gene = first_nonzero_element_index(entrance)
            first_enter_of_gene = first_nonzero_element_index(entrance) + len(left)
        except ValueError:
            # print('''Cannot establish entrance (id: {})'''.format(tree))
            raise ValueError('''in find_enter()''')
        else: # +2 is neede because there is 1 empty element and we must start with 1 (not to overwrite last element in taxa -> +1+1 = +2
            return first_enter_of_gene + 2
    else: # print('''Single target in tree only''')
        # return presence_of_gene.index(1) + len(entrance) + 1 # it could be larger number than 1!!!!!
        # return first_nonzero_element_index(presence_of_gene) + len(entrance) # +1 is wrong; wrong order, i.e. internal nodes before taxa
        return first_nonzero_element_index(presence_of_gene) # +1 is wrong

def remove_duplicite_line(csv_file, ignore_bs=True):
    '''
    https://stackoverflow.com/questions/15741564/removing-duplicate-rows-from-a-csv-file-using-a-python-script
    Ignore bs for comparison added in 'v5.54' - it was troublesome for nested summary
    returns file without duplicite lines
    '''
    if '_unique_rows' not in csv_file:
        out_file = csv_file[:-4] + '_unique_rows.csv'
        with open(csv_file,'r') as csv, open(out_file,'w') as out:
            seen = set() # set for fast O(1) amortized lookup
            for line in csv:
                line2compare = line
                if ignore_bs:
                    list2compare = line.split('\t')
                    bs2 = int(list2compare[4]) # first attempt to save the highest values of BS (v7.40 - not trying very much)
                    del list2compare[4] # delete BS
                    line2compare = '\t'.join(list2compare)
                if line2compare in seen: continue # skip duplicate
                seen.add(line2compare)
                out.write(line)
        return out_file
    return csv_file

def remove_duplicite_line_mix(csv_file): # tag: "RDLM"
    '''
    https://stackoverflow.com/questions/15741564/removing-duplicate-rows-from-a-csv-file-using-a-python-script
    Ignore bs for comparison added in 'v5.54' - it was troublesome for nested summary
    returns file without duplicite lines and without mix, where specific algal group present
    ### edit Q256509051551: unfortunately return only single color for every tree, even if there should be multiple colors !!!!!!! - solved in 'v8.2'
    # header: algal group	tree id	tree size	subgroup id	b	DEEZ	DEEG	DEEH	DEEE	DERV	DErc	DEnp	DEpt	DEpv
    '''
    print('''uniquarum: {}'''.format(csv_file))
    if '_unique_rows' not in csv_file:
        out_file = csv_file[:-4] + '_unique_rows.csv'
        with open(csv_file,'r') as csv:
            seen_mix = dict() # set for fast O(1) amortized lookup
            for line in csv:
                #print(line)
                #input()
                #"""
                list2compare = line.split('\t', 4)
                #print(list2compare)
                #input()
                algal_group = list2compare[0]
                tree_id = list2compare[1]
                subgroup_id = list2compare[3]
                tree_size = list2compare[2]
                taxon_composition = list2compare[4]
                # key = tree_id + subgroup_id 
                key = tree_id + subgroup_id + algal_group # differentiate the colors
                value = line
                if key not in seen_mix:
                    # print('''key "{}" not found'''.format(key))
                    seen_mix[key] = value
                elif algal_group == 'mix':
                    # print('''found mix''')
                    pass
                else:
                    seen_mix[key] = value
                #"""
        with open(out_file,'w', newline='\n') as out:
            for key in seen_mix:
                out.write(seen_mix[key])
        return out_file
    else:
        return csv_file

def row2csv(row):
    '''
    input: pandas row
    output: line to be written into csv
    '''
    row['b'] = str(row['b']).zfill(3)
    new_row = ''
    for cell in row:
        new_row += str(cell) + '\t'
    return new_row[:-1] + '\n'

def target_decontamination(unique_rows_file, log_dir, output):
    '''
    vyřešit dvojitý output
    '''
    # log_dir = 'Q:/____STEPA40db_results/_important/Florida/20190327185420_13977_5.32_Florida_b=75/log/'
    #unique_rows_file = 'Q:/____STEPA40db_results/_important/Florida/Hollywood_v0.2.csv'
    #unique_rows_file = 'Q:/____STEPA40db_results/_important/Florida/test18.csv'
    #unique_rows_onlySL = 'Q:/____STEPA40db_results/_important/Florida/Hollywood_v0.3_nomix.csv'
    contamination = 'Q:/____STEPA40db_input/contamination/target/'
    cnts = glob.glob(contamination + '*.cnt')
    # header = 'algal group	tree id	tree size	subgroup id	b	DEEZ	DEEG	DEEH	DEEE	DERV	DErc	DEnp	DEpt	DEpv	node\n'
    header = 'algal group	tree id	tree size	subgroup id	b	DEEZ	DEEG	DEEH	DEEE	DERV	DErc	DEnp	DEpt	DEpv\n'
    tmp = unique_rows_file
    i = 0
    for cnt in cnts:
        taxon = cnt[1+cnt.rfind('\\'):cnt.rfind('.cnt')]
        print(taxon)
        with open(cnt) as f:
            # sl = f.readlines() # contains '\n's
            sl = f.read().splitlines()
            #print(sl)
        print('''panda's starting...''')
        df = pd.read_csv(tmp, sep='\t')
        print('''panda read''')
        # tmp = 'tmp_' + str(i) + '.csv'
        tmp = 'tmp_' + taxon + '.csv'
        with open(tmp, 'w', newline='\n') as out:
            out.write(header)
            i = 0
            new_row = ''
            for index, row in df.iterrows():
                #print('''index: {};\nrow: {};'''.format(index, row))
                #print('''row[algal group] = {}'''.format(row['algal group']))
                # if row['algal group'] != 'mix':
                if True: # starting 'v6.03' I do care about mix!
                    # Starting 'v8.00' I do care about mix
                    if row[taxon] > 0: # I have to deal with 2090 trees/logs
                        sub_ids = row['subgroup id'].split('.')
                        # print(sub_ids)
                        # print('''tree_id'''.format(row['tree id']))
                        logf = log_dir + row['tree id'] + '.log'
                        # print(logf)
                        with open(logf) as log:
                            difference = dict()
                            for line in log:
                                if line.startswith('X0') or line.startswith('Q0'):
                                    total = line.count(taxon)
                                    cnt_count = 0
                                    for item in sl:
                                        cnt_count += line.count(item)
                                    if cnt_count == 0:
                                        difference[line.split(':')[0]] = 0
                                    else:
                                        difference[line.split(':')[0]] = cnt_count
                                else: # lines which do not start with 'Q0'/'X0' are not important
                                    pass
                        total_difference = 0
                        #print(difference)
                        for sub_id in sub_ids:
                            total_difference += difference[sub_id]
                        row[taxon] -= total_difference
                        out.write(row2csv(row))
                    else: # no problematic taxon present in the group (DERV/DEEE/DErc)
                        out.write(row2csv(row))
                else: # I do not care about mix in MeEGT
                    pass
                    # out.write(row2csv(row))
    shutil.copy(tmp, output)

def summarize_results(csv_tab, csv_results, copy):
    '''
    The function takes csv file containing all possibilities of all trees (and group) and takes only the best ones (highest BS).
    Creates data for each node of target tree.
    '''
    # print('''summarize_results()''')
    nNodes = 1 + 2 * len(target) - 1 # number of nodes in the tree is dependent only on number of leaves, i.e. terminal nodes; i.e. is independent on the topology!
    # 1+ is for the empty node (otherwaise root doubled in tree-plot)
    # field = [0] * nNodes # cannot be used - it creates only multiple shortcuts to the same field
    # node_composition = {k:field for k,v in dic_algae.items()} # items are not indempendent
    node_composition = {k:[0] * (nNodes) for k,v in dic_algae.items()} # creates an initialized dictionary with the very same keys as dic_algae
    # print(node_composition)
    # pprint(node_composition)
    header = next(csv_tab)
    results_dict = dict()
    with open(results_node_file, 'w', newline='\n') as out_node:
        out_node.write('\t'.join(header) + '\t' + 'node' + '\n')
    with open(results_node_nested_file, 'w', newline='\n') as out_node_nested:
        out_node_nested.write('\t'.join(header) + '\t' + 'node' + '\n')
    print(header)
    for algal_group, tree_id, tree_size, subgroup_id, b, *leaves in csv_tab: # *leaves should contain all (=any number of) terminal taxa
        # leaves = [int(i) for i in leaves] # '0.0' sometimes present
        for i,j in enumerate(leaves):
            if j == '0.0':
                leaves[i] = 0
            else:
                leaves[i] = int(j)
        # print('''leaves: {leaves}'''.format(**locals()))
        # print('''tree_id: {tree_id}; algal_group: {algal_group}; subgroup_id: {subgroup_id}'''.format(**locals()))
        if False: #algal_group == 'mix': #or int(b) == 102 or int(b) == 101 or int(b) == 111:
            pass
        elif sum(leaves) == 0:
            print('''Tree {} does not have any targets in some branch (cause by ignored taxa).'''.format(tree_id))
            # If some taxa from target tree are ignored, they still appear in the sorting; i.e. they have their own Qxxx number, but only zero color
        elif int(b) == 102:
            pass
            print('''{}: {} ({})'''.format(tree_id, b, algal_group))
        elif '.' in subgroup_id:
            pass
            print('''{}: {}'''.format(tree_id, subgroup_id))
        else:
            '''
            Since the target decontamination is done quite late in pipe line, composition of target can be different, i.e. lst2str(leaves)
            '''
            if copy_tree: # if copy:
                sorted_id_path = algal_group + '/' + b + '-' + tree_id + '-' + subgroup_id + '-' + lst2str(leaves) + '.' + tree_ext
                #id_path = algal_group + '/' + tree_id + '-' + b + '-' + subgroup_id + '-' + + '.' + tree_ext
                id_path = algal_group + '/' + tree_id + '-' + b + '-' + subgroup_id + '*.' + tree_ext
                id_path = algal_group + '/' + tree_id + '*.' + tree_ext
                #print('''sorted id path: {}'''.format(sorted_id_path))
                #print('''id path: {}'''.format(id_path))
                src = glob.glob(path + id_path)[0]
                #src = path + sorted_id_path
                dst = path_selected + sorted_id_path
                shutil.copy(src,dst)
            # print('''find_enter(leaves): {}'''.format(find_enter(leaves)))
            found_enter = find_enter(leaves)
            ##### writing into file #####
            # print('''{} - {} - {} - {} - {} -- {}'''.format(algal_group, tree_id, subgroup_id, b, leaves, found_enter))
            # print('''{}'''.format(leaves))
            # print('''{}'''.format('\t'.join(leaves)))
            # ','.join(str(v) for v in value_list)
            node_line = algal_group + '\t' + tree_id + '\t' + subgroup_id + '\t' + b + '\t' + '\t'.join(str(v) for v in leaves) + '\t' + str(found_enter) + '\n'
#            node_line = algal_group + '\t' + tree_id + '\t' + subgroup_id + '\t' + b + '\t' + found_enter
            key = tree_id + '_' + subgroup_id
            value = algal_group + '\t' + tree_id + '\t' + subgroup_id + '\t' + b + '\t' + '\t'.join(str(v) for v in leaves) + '\t' + str(found_enter) + '\n'
            if key not in results_dict:
                results_dict[key] = value
            elif algal_group == 'mix':
                pass
            else:
                results_dict[key] = value
            # with open(results_node_file, 'a+', newline='\n') as out_node:
                # out_node.write(node_line)
            # with open(results_node_nested_file, 'a+', newline='\n') as out_node_nested:
                # out_node_nested.write(node_line)
            ##### writing into file end #####
            node_composition[algal_group][found_enter] += 1 # SUMMARY UPDATE
            dst2 = path_selected + str(found_enter) + '/' + algal_group + '/'
            os.makedirs(dst2, exist_ok=True)
            # if not os.path.isdir(dst2):
                # os.mkdir(dst2)
            # dst3 = dst2 + tree_id + '.bip'
            if copy_tree:
                shutil.copy(src,dst2)
            # print('''alg. group: {}'''.format(node_composition[algal_group]))
    with open(results_node_nested_file, 'a+', newline='\n') as eliminated:
        for key in results_dict:
            eliminated.write(results_dict[key])
    # print(node_composition)
    # categories, i. e. nodes of the tree are in order specified by R
    categories = [i for i in range(nNodes)] # 0..n
    # 0..nTaxa-1, root, nTaxa+1..2nTaxa-1
    categories = target + ['empty', 'root'] + [i+len(target)+2 for i in range(nNodes-len(target)-2)]
    #categories = target + ['empty', 'root'] + ['n' + str(i+len(target)+2) for i in range(nNodes-len(target)-2)]
    # print(categories)
    #['L', 'LxG', 'G', 'LGxE', 'E', 'LGExR', 'R', 'LGERxc', 'C', 'LGERcxt', 'T'] # how to universaly change? node id?
    csv_results.writerow(['algal group'] + categories)
    # print(categories)
    for key in dic_algae:
        if False: # key == 'prokaryota':
            node_composition[key][6] -= 520 # bacterial contamination in RCo - index 6
        csv_results.writerow([key] + node_composition[key])
        # print('''{}: {}'''.format(key, node_composition[key]))
    return None

"""
*************************************************
----------------------  R  ----------------------
*************************************************
"""

def list2c(name, lst, norm):
    '''
    python list to R column
    '''
    s = ''
    i = 0
    for n in lst:
        i += 1
        if i == len(lst)+1:
            s += '0' + ','
        s += str(round(int(n) / float(norm), 5))  + ', '
    s = s[:-2] + ')\n'
    # print(str(name) + ' = c(' + s)
    return str(name) + ' = c(' + s

def csv2r(file):
    maximum = 0
    rdata = ''
    names = []
    with open(file, 'r') as f:
        data = f.readlines()
        # print(data)
    for row in data:
        rowd = row.strip().split('\t')
        print('''rowd: {}'''.format(rowd))
        name = rowd[0]
        if name != 'algal group': # i.e. skip the header line
            names.append(name)
            freq = rowd[1:]
            print('''name: {}'''.format(name))
            print('''freq: {}'''.format(freq))
            nfreq = [int(i) for i in freq]
            # print(name)
            # print(freq)
            try:
                m = int(max(nfreq))
                # print(m)
            except ValueError:
                pass
            else:
                if m > maximum:
                    maximum = m
    for row in data:
        rowd = row.strip().split('\t')
        if rowd[0] != 'algal group':
            rdata +=  list2c(rowd[0], rowd[1:], maximum)
    # print(maximum)
    bars = '''bars <- nodebar(dat, cols=1:14, color = c(a=uni, '''
    ph = ''
    for code, phylum in zip(range(ord('c'), ord('g') + 1), names):
        ph = ph + str(code) + '=' + str(phylum[:3]) + ', '
        # b=gla, c=rho, d=hap, e=och, f=cry, g=din, h=chr, i=vir, j=pau, k=chl, l=pro, m=cya, n=kin
    bars += ph + '''), position = "dodge")\n'''
    return rdata, bars # bars not used (v5.3)

def lst2r(lst):
    '''
    takes list and returns R-like column as string
    '''
    sep = ', '
    col = 'c('
    for item in lst:
        col += str(item) + sep
    return col[:-len(sep)] + ')'

def create_meegt_R(file, taxon_number, out_name):
    '''r_input_file, len(target), r_script_name
    input variables:
        file         ... input file in 2D matrix form (color x node)
        taxon_number ... number of taxa in the target tree, i.e. euglenids
        out_name     ... name of vreated R script
    output:
        none; creates R script, which is capable of visualize tree and bar plot for each node
    '''
    ruller = [1] * (2*taxon_number) # ruller is needed to show all bar plots in the same scale
    ruller[taxon_number] = 0
    unary = 'unary = ' + lst2r(ruller) + '\n'
    # print(ruller, unary)
        # the length of unary (=ruller) depends on number of taxa in target tree
        # len(unary) = 2*len(tree_taxa) + 1
        # the middle element is always set to ZERO, all others are ONES
    with open(out_name, 'w', newline='\n') as r:
        preambule = \
'''#!/usr/bin/env Rscript
# This script was automaticaly created by ''' + "'" + script_name + "'" + ''' on SABETHES with python 3.7
library(treeio)
library(ggplot2)
library(ggimage)
library(ggtree)

rm(list=ls()) # clear workspace
# setwd("Q://Dropbox/STEPA41db/_presentation/data")

uni <- "black"
#uni <- "white"
gla <- "#00B0F0" # blue
rho <- "#FF0000" # red
hap <- "#FAC090"
och <- "#984807"
cry <- "#E46C0A"
din <- "#FFCCAA"
chr <- "#D10F32"
vir <- "#00B050" # green
pau <- "#FFFF00" # yellow
chl <- "#92D050"
#cya <- "#008B8B" # darkcyan
#pro <- "#8B7355" # burlywood4
#kin <- "#A9A9A9" # darkgrey
mixc <- "black"

tr <- read.newick(file="''' + target_tree_file + '''") # should be readed from the string (beggining of the py script)
tree <- ggtree(tr, layout="rectangular") + geom_tiplab(color='blue') + ggplot2::xlim(0, 10)
ggtree(tr) + geom_text2(aes(subset=!isTip, label=node), hjust=-.3) + geom_tiplab()

#unary = c(1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1) # for -RCo dele the first 1
#unary = c(1, 1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1)

'''

        data, bars = csv2r(file)

        end = \
'''
#dat <- data.frame(a = unary, b = glaucophyta, c = rhodophyta, d = haptophyta, e = ochrophyta, f = cryptophyta, g = dinophyta, h = chromerida, i = viridiplantae, j = paulinella, k = chlorarachniophyta, l = prokaryota, m = cyanophyta, n = kinetoplastea)
dat <- data.frame(a = unary, b = glaucophyta, c = rhodophyta, d = haptophyta, e = ochrophyta, f = cryptophyta, g = dinophyta, h = chromerida, i = viridiplantae, j = paulinella, k = chlorarachniophyta, l = mix)
dat$node <- 1:''' + str(2*len(target)) + ''' # number of all nodes (rows), i.e. 18 

# cols - number of colorful categories
# bars <- nodebar(dat, cols=1:14, color = c(a=uni, b=gla, c=rho, d=hap, e=och, f=cry, g=din, h=chr, i=vir, j=pau, k=chl, l=pro, m=cya, n=kin), position = "dodge")
bars <- nodebar(dat, cols=1:12, color = c(a=uni, b=gla, c=rho, d=hap, e=och, f=cry, g=din, h=chr, i=vir, j=pau, k=chl, l=mixc), position = "dodge")
inset(tree, bars, x="branch", width=1, height=1, vjust=-0.35)

ggsave("meegt_automat_''' + version + '''.png")
'''
        r.writelines(preambule)
        r.write(unary)
        r.write(data)
        r.writelines(end)
    return None

"""
*************************************************
------------------- END OF R  -------------------
*************************************************
"""

if __name__ == "__main__":
    threader = 0
    only_R = False
    if not only_R:
        create_folders(path)
        # summarize_results(1,2)
        analyzing = True
        #print('''################################''')
        if analyzing:
            #print('''######################''')
            with open(results_file, 'w', newline='') as results, open(results_nested_file, 'w', newline='') as nested:
                header = ['algal group', 'tree id', 'tree size', 'subgroup id', 'b'] + target #+ ['node']
                results_csv = csv.writer(results, delimiter = '\t') # better using Dictwriter???
                results_csv.writerow(header)
                nested_csv = csv.writer(nested, delimiter = '\t')
                nested_csv.writerow(header)
                w = 0
                #print('''################################''')
                for file in trees:
                    #print('''################################''')
                    print(file)
                    threader += 1
                    #print('''t%t=t: {}%{}={}'''.format(threader, all_threads, threader%all_threads))
                    zbytek = threader%all_threads
                    #print('''zbytek({})==thread({})'''.format(zbytek, thread))
                    if zbytek == thread:
                        #print(file)
                        # file = '000117.bip'
                        log_file = log_path + file[len_intree:-3] + 'log'
                        # print(log_file)
                        tree_size = os.stat(file).st_size
                        # print('''tree-size: {}'''.format(tree_size))
                        w += 1
                        if w > test_size:
                            break
                        else:
                            # print('''\n{:~^50}'''.format('~'))
                            # print('''{:~^50}'''.format(file))
                            # print('''{:~^50}'''.format('~'))
                            # print('''{:~^50}'''.format(file[file.rfind('\\')+1:-4]))
                            with open(file) as tree_file, open(log_file, 'w') as tout:
                                tree = Newick(tree_file.read().strip()) # strip better in class declaration!?
                                if tree.is_newick():
                                    # left, right = tree.topology()
                                    nIgnored = Newick.count_ignored(tree)
                                    nTaxa = Newick.count_taxa(tree)
                                    nTarget = Newick.count_target(tree, replacement=False)
                                    if nTaxa - nIgnored > 3 and nTarget - nIgnored > 0:
                                        l, t, q_list, target_composition = tree.target_group_determination()
                                        # print('''q_list = {q_list}'''.format(**locals()))
                                        # print('''target_composition = {target_composition}'''.format(**locals()))
                                        for q in q_list:
                                            # print('''q: {}'''.format(q))
                                            t.analyse_sisterhoods(q)
                                """
                                            # print('''{:~^50}'''.format('~'))
                                            # print('''{:~^50}'''.format(q))
                                            # print('''{:~^50}'''.format('~'))
                                            # # t.get_sisterhoods(q)
                                            # # Newick.analyse_sisterhoods(t, q)
                                            #t.get_closest_sisterhoods(q)
                                            #break
                                        #print('''tree:\n {tree}'''.format(**locals()))
                                        # print('''count_taxa() = {}'''.format(tree.count_taxa()))
                                        left, right = tree.topology()
                                        b = tree.get_bootstraps(right)
                                        nTarget = tree.count_target()
                                        s = 0
                                        for g in target_composition:
                                            s += sum(g)
                                        # print('''left brackets: {left}'''.format(**locals()))
                                        # print('''right brackets: {right}'''.format(**locals()))
                                        # print('''bootstraps: {}'''.format(b))
                                        # print('''Number of algae in the tree: {}'''.format(tree.count_algae()))
                                        # print('''Number of target in the tree: {}'''.format(tree.count_target()))
                                        # print('''Number of taxa in the tree: {}'''.format(tree.count_taxa()))
                                        # print(tree.get_id('EEG')) # position of first occurrence in the tree
                                        # print(tree.get_target_ids())
                                        #print(tree)
                                        #print('''l = {l}\nt = {t}\nq_list = {q_list}\ntarget_composition = {target_composition}'''.format(**locals()))
                                        # print('''l = {l}\nq_list = {q_list}\ntarget_composition = {target_composition}'''.format(**locals()))
                                        # print('''l = {l}\nq_list = {q_list}\ntarget_composition = {target_composition}'''.format(**locals()))
                                        # print(target_composition)
                                        # print(l)
                                        # print('''count_target() = {nTarget}, groups_target = {s}'''.format(**locals()))
                                        #print(tree)
                                        #print(t)
                                        # print('''q_list = {q_list}'''.format(**locals()))
                                        # print(tree.analyse_root_epsilon())
                                    # else:
                                        # print('''{} is not newick tree!'''.format(file))
                                """
                            if w % trees_done == 0:
                                t = int(timer())
                                p = int(100 * w/test_size)
                                print('''I've already sorted {w} trees ({p} %) in {t} sec.'''.format(**locals()))
        print('''{:~^50}'''.format('#'))
        print('''{:~^50}'''.format('Summary starts'))
        print('''{:~^50}'''.format('#'))
        if not analyzing:
            results_file = result_dir + '20190330195803_13977_v5.34_dip_decon_b=75_unique_rows.csv'
            r_input_file = results_file[:-4] + '_unique_rows_summary.csv'
        # decontaminated = 'Q:/____STEPA40db_results/_important/Florida/Hollywood_v5.33.csv'
        if only_nested:
            #print('''Let's go -- only nested''')
            results_nested_unique_file = remove_duplicite_line_mix(results_nested_file)
            # results_nested_unique_file = remove_duplicite_line(results_nested_file)
            target_decontamination(results_nested_unique_file, log_path, results_nested_unique_decontaminated_file)
            decontaminated = results_nested_unique_decontaminated_file
            summary_file = summary_nested_unique_decontaminated_file

            # summary_nested_file = results_nested_file[:-4] + '_summary.csv'
            #r_input_file = summary_nested_file
            #results_node_file = results_node_nested_file
            # target_decontamination(results_nested_unique_file, log_path, results_nested_unique_decontaminated_file)
#            target_decontamination(results_nested_file, log_path, results_nested_unique_decontaminated_file) # why not UNIQUE??????? used till v7.00
        else:
            results_unique_file = remove_duplicite_line_mix(results_file)
            summary_file = summary_unique_decontaminated_file
            target_decontamination(results_unique_file, log_path, results_unique_decontaminated_file)
            decontaminated = results_unique_decontaminated_file
            
            # summary_file = results_file[:-4] + '_summary' + str(thread) + '.csv' # create file

        with open(decontaminated, 'r', newline='') as results, open(summary_file, 'w', newline='\n') as summary:
            result_csv = csv.reader(results, delimiter = '\t')
            summary_csv = csv.writer(summary, delimiter = '\t')
            # summary_csv.writerow(['glaucophyta'] + ['rhodophyta'] + ['haptophyta'] + ['ochrophyta'] + ['cryptophyta']+ ['viridiplantae']+ ['paulinella'] + ['chlorarachniophyta'] + ['kinetoplastea'])
            summarize_results(result_csv, summary_csv, copy=analyzing) # input, output
        stop = timer()
        try:
            print('''Sorting of {} trees took {} s.'''.format(w, stop-start))
        except NameError:
            print('''Summarizing is done.''')
    else:
        pass
        #file = 'Q:/____STEPA40db_results/20190321210138_139_cryptophyta-DErc_b=75_unique_rows_summary.csv'
        #print(csv2r(file))
    #"""    
        print('''{:~^50}'''.format('~'))
        print('''{:~^50}'''.format('R'))
        print('''{:~^50}'''.format('~'))
        create_meegt_R(r_input_file, len(target), r_script_name)
        os.system(r_script_name)
        print('''Analysis done in {} s.'''.format(timer() - start))
    #"""
