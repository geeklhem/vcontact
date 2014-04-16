"""General options"""
import os
import logging
data_folder = os.path.expanduser("~/data/")

cache_folder = data_folder+"cache/"

folders = {}
for f in ("cache","proteins","contigs","log"):
    folders[f] = data_folder+f+"/"
    if not os.path.exists(folders[f]):
        os.mkdir(folders[f])
        logging.info("Created {}/ directory into {}".format(data_folder,f))
        
contig_name_regex = '([NY][CP]_[0-9]*|^[0-9]{2,3}[A-Z]{3}[0-9_]*)'

overwrite = False

refseq_faa = data_folder+"refseq/phage.protein.faa"

# Blast options 
blast_cmd = "blastall -p blastp  -m 8 -a 24 -e 0.00001"
