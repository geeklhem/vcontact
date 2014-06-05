import genome_clusters 
import pc_matrix 
import cPickle as pickle
import numpy as np
import pandas
import subprocess
import options
import exports


def extract_pelagiphages(contigs,pelagiphages,names,ntw):
    contigs = contigs.reset_index()
    refseq_pos = contigs.query('origin == "refseq_jan14"').pos.values
    refpel_pos = [list(refseq_pos).index(x) for x in contigs.query("name in pelagiphages").pos.values]
    c_pel = [int(x) for x in contigs.query("name in pelagiphages").mcl_cluster.drop_duplicates().values]
    candidates = contigs.query('mcl_cluster in c_pel').pos.values
    names_pel = dict([(n,[]) for n in names])
    
    a = ntw.tocsr()
    found = dict(zip(pelagiphages,[[],[],[],[]]))
    contigs = contigs.set_index("pos")
    for c in candidates:
        filtered = a[c,refseq_pos]
        if filtered.getnnz():
            p = filtered.indices[filtered.data.argmax()]
            for n,order in zip(names,refpel_pos):
                if p == order:
                    contigs.ix[c,"genus"] = n + '-like'
                    names_pel[n].append(contigs.ix[c,"name"])
                    
    return contigs,names_pel 
