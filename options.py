"""General options"""
import os
import logging
data_folder = os.path.expanduser("~/data/")

cache_folder = data_folder+"cache/"

colors = {"refseq":"#B31B1B",
          "tara":"#32C6A6",
          "Myoviridae": "#4daf4a",
          "Siphoviridae" :"#984ea3",
          "Podoviridae":"#ff7f00",
          "red":"#e41a1c",
          "blue":"#377eb8",
          "green":"#4daf4a",
          "violet":"#984ea3",
          "orange":"#ff7f00",
          "yellow":"#ffff33",
          "brown":"#a65628",
          "pink":"#f781bf",
          "grey":"#999999",
          "Non affiliated":"#999999",
          "Non-affiliated":"#999999",
          "Inoviridae":"#f781bf",}


folders = {}
for f in ("cache","proteins","contigs","log","results","modules"):
    folders[f] = data_folder+f+"/"
    if not os.path.exists(folders[f]):
        os.mkdir(folders[f])
        logging.info("Created {}/ directory into {}".format(f,data_folder))
        
contig_name_regex = '([NY][CP]_[0-9]*|^[0-9]{2,3}[A-Z]{3}[0-9_]*)'

overwrite = False

refseq_faa = data_folder+"refseq/phage.protein.faa"

# Blast options 
blast_cmd = "blastall -p blastp  -m 8 -a 24 -e 0.00001"


keywords =  ['ATPase', 'wedge', 'junction', 'assembly',  'activator',  'baseplate',  'capsid',  'chaperone',  'diphosphate',  'endonuclease',  'exonuclease',  'fiber',  'head',  'helicase',  'helix',  'homing',  'hydrolase',  'inhibitor',  'injection',  'integrase',  'kinase',  'ligase',  'lysis',  'lysozyme',  'membrane',  'methylase',  'methyltransferase',  'neck',  'nuclease',  'polymerase',  'portal',  'plate',  'scaffold',  'primase',  'prohead',  'protease',  'recombination',  'recombinase',  'transposase',  'reductase',  'repair',  'regulator',  'replication',  'repressor',  'ribonucleoside',  'ribonucleotide',  'structural',  'synthase',  'tail',  'tube',  'tRNA',  'terminase',  'transcription',  'transcriptional',  'transferase',  'virion']

pelagiphages = ["NC_020481","NC_020482","NC_020483","NC_020484"]
names = ["HTVC010P", "HTVC011P", "HTVC019P", "HTVC008M"]

"""
refseq_fi = {}
files = {"refseq":refseq_fi, 
         "tara":tara_fi}
"""
