"""main.py""" 

import genomes
import protein_clusters
import pc_matrix
import options
import cPickle as pickle 

g = genomes.Genomes(True,True)
p = protein_clusters.ProteinClusters(options.data_folder+"tara/tara_contigs10_and_refseq_pblast_mcl20.clusters")

pcm = pc_matrix.PCMatrix(p.data.proteins,p.data.clusters,
                         g.data.proteins, g.data.contigs)
                                    
pcm.ntw = pcm.network(pcm.matrix)

with open(options.data_folder+"pcm.pkle", 'wb') as f:
        pickle.dump(pcm,f)

pcm.to_mcl(pcm.ntw,"{0}.ntwk".format(fi))
