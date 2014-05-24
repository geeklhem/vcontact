import pandas
from Bio import SeqIO
import logging
import os
import options
import numpy as np
import re

class Genomes(object):
    """
    Data about reference genomes, object built from the refseq gbff file.
    """
    
    def __init__(self, refseq=True,tara=False,refseq_contigs=False):
        """
        
        Arguments:
        - `fi`:
        """
        if refseq and not tara:
            self.data = self.load_refseq()
        if tara and not refseq:
            self.data = self.load_tara()
        if refseq and tara and not refseq_contigs:
            self.data = pandas.HDFStore(options.cache_folder+"all_contigs.h5")
            if not "contigs" in self.data: 
                refseq = self.load_refseq()
                tara = self.load_tara()
                self.data.append('contigs',pandas.concat((tara.contigs,refseq.contigs)),format='table', data_columns=True)
                self.data.append('proteins',pandas.concat((tara.proteins,refseq.proteins)),format='table', data_columns=True)
                self.data.append('metadata',tara.metadata,format='table')
        if refseq and refseq_contigs and not tara:
            self.data = pandas.HDFStore(options.cache_folder+"all_refseq_contigs.h5")
            if not "contigs" in self.data: 
                refseq = self.load_refseq()
                contigs = self.load_contigs_refseq()
                self.data.append('contigs',pandas.concat((contigs.contigs,refseq.contigs)),format='table', data_columns=True)
                self.data.append('proteins',pandas.concat((contigs.proteins,refseq.proteins)),format='table', data_columns=True)
        else:
            logging.info("Nothing loaded")
        
    def load_tara(self,h5_name = "tara.h5",fi_faa=None,fi_fna=None,fi_meta=None):
        """
        Load TARA metadata
        """
        store =  pandas.HDFStore(options.cache_folder+h5_name)


        if not "metadata" in store:
            # Metadata 
            fi_meta = os.path.expanduser(options.data_folder+"tara/tara_station_metadata.csv") if fi_meta == None else os.path.expanduser(fi_meta)
            meta = pandas.read_csv(fi_meta)
            meta["Station"] = ["".join(k.split("_")) for k in meta["Station"]]
            meta = meta.set_index("Station")
            store.append('metadata',meta,format='table')
            
        # Contigs
        if not "contigs" in store:
            fi_fna = os.path.expanduser(options.data_folder+"tara/contigs10.txt") if fi_fna == None else os.path.expanduser(fi_fna) 
            contigs = {"name":[],
                      "origin":[]}

            with open(fi_fna,"rb") as fna:
                for l in fna:
                    contigs["name"].append(l.rstrip())
                    contigs["origin"].append(l.split("_")[0])

            contigs = pandas.DataFrame(contigs).set_index("name")
            contigs["family"] = None
            contigs["genus"] = None

            store.append('contigs',contigs,format='table',data_columns=True)

        if not "proteins" in store:
            fi_faa = os.path.expanduser(options.data_folder+"tara/tara_c10_prots.txt") if fi_faa == None else os.path.expanduser(fi_faa)
            proteins = {"protein_id":[],
                        "contig":[]}

            with open(fi_faa,"rb") as faa:
                for l in faa:

                    proteins["protein_id"].append(l.split("#")[0].rstrip()[1:])
                    proteins["contig"].append("_".join(l.split("_")[:2])[1:])
                
            proteins = pandas.DataFrame(proteins).set_index("protein_id")
            proteins["function"] = None
            store.append('proteins',proteins,format='table',data_columns=True)
        return store
            
    def load_refseq(self,fi=None,fi_taxa = None,contig_id_eq_prot_id=True,h5_name="refseq.h5"):
        """
        INPUT :
        - fi (str), a genebak file.
        OUTPUT : 
        - names (list) The genomes names (Index2genome)
        - taxonomy (dict of list) map genome name on the list of taxonomic levels
        - prot2genome (dict) map protein name to genome name
        """

        # Loading taxonomic classes 
        fi = os.path.expanduser(options.data_folder+"refseq/phage.genomic.gbff") if not fi else fi 
        fi_taxa = os.path.expanduser(options.data_folder+"ncbi_taxdump/refseqphage_taxa.pkl") if not fi_taxa else fi_taxa  


        store =  pandas.HDFStore(options.cache_folder+h5_name)

        if "contigs" not in store or "proteins" not in store:
            taxa = pandas.read_pickle(fi_taxa)
            families = frozenset(taxa[taxa["Rank"]=="family"]["Name"])
            genera = frozenset(taxa[taxa["Rank"]=="genus"]["Name"])


            genome = {"name":[],
                      "family":[],
                      "genus":[],
                      "size":[],}

            proteins = {"protein_id":[],
                        "contig":[],
                        "function":[],}

            # Reading the gbff file thanks to BioPython
            rec = SeqIO.parse(fi, "gb")
            for r in rec:
                name = r.id.split(".")[0]
                genome["name"].append(name)

                # Taxonomy
                genome["size"].append(len(r.seq))
                genome["family"].append(None)
                genome["genus"].append(None)
                for t in r.annotations["taxonomy"]:
                    if t in families:
                        genome["family"][-1] = t
                    elif t in genera:
                        genome["genus"][-1]  = t
                #Proteins 
                for feature in r.features:
                    if "protein_id" in feature.qualifiers:
                        proteins["protein_id"].append(feature.qualifiers["protein_id"][0].split(".")[0])
                        proteins["contig"].append(name)
                        proteins["function"].append(feature.qualifiers["product"][0])

                 

            genome = pandas.DataFrame(genome).set_index("name")
            genome["origin"] = "refseq_jan14"
            proteins = pandas.DataFrame(proteins).set_index("protein_id")
            store.append('contigs',genome,format='table',data_columns=True)
            store.append('proteins',proteins,format='table',data_columns=True)

            nf = len(genome) - genome.family.count() 
            ng = len(genome) - genome.genus.count()
            if nf!=0 or ng!=0 :
                logging.warning("{0} genomes without family and {1} genomes without genus.".format(nf,ng))

        return store 

        

    
    def load_contigs_refseq(self,fi=None,fi_taxa=None,h5_name="refseq_contigs.h5"):
        """
        INPUT :
        - fi (str), a genebak file.
        OUTPUT : 
        - names (list) The genomes names (Index2genome)
        - taxonomy (dict of list) map genome name on the list of taxonomic levels
        - prot2genome (dict) map protein name to genome name
        """
        logging.warning("Loading Refseq contigs")
        # Loading taxonomic classes 
        fi = os.path.expanduser(options.data_folder+"contigs_refseq/proteins_names") if not fi else fi 
        fi_taxa = os.path.expanduser(options.data_folder+"ncbi_taxdump/refseqphage_taxa.pkl") if not fi_taxa else fi_taxa  


        store =  pandas.HDFStore(options.cache_folder+h5_name)

        if "contigs" not in store or "proteins" not in store:
            taxa = pandas.read_pickle(fi_taxa)
            families = frozenset(taxa[taxa["Rank"]=="family"]["Name"])
            genera = frozenset(taxa[taxa["Rank"]=="genus"]["Name"])


            genome = {"name":[],
                      "ref_genome":[],
                      "origin":[]}

            proteins = {"protein_id":[],
                        "contig":[],}

            # Reading the gbff file thanks to BioPython
            with open(fi,'rb') as f:
                for r in f:
                    prot_name = pid(r)
                    splited = prot_name.split("_")
                    contig_name = "_".join(splited[:-1])
                    genome_name = "_".join(splited[:2])

                    if len(splited) != 4:
                        logging.error("line: {}\n pid: {}\n contig {} \n genome {}".format(r,prot_name,contig_name,genome_name))
                        
                    
                    genome["origin"].append("contigs_{}_refseq_jan14".format(r.split("-")[0]))
                    genome["name"].append(contig_name)
                    genome["ref_genome"].append(genome_name)
                    proteins["protein_id"].append(prot_name)
                    proteins["contig"].append(contig_name)

            genome = pandas.DataFrame(genome).set_index("name").drop_duplicates()
            refseq = pandas.HDFStore(options.cache_folder+"refseq.h5")
            genome = pandas.merge(genome,refseq.contigs.drop("origin",1),left_on="ref_genome",right_index="True")
            
            proteins = pandas.DataFrame(proteins).set_index("protein_id")
            
            store.append('contigs',genome,format='table',data_columns=True)
            store.append('proteins',proteins,format='table',data_columns=True)
            
            nf = len(genome) - genome.family.count() 
            ng = len(genome) - genome.genus.count()
            if nf!=0 or ng!=0 :
                logging.warning("{0} contigs without family and {1} contigs without genus.".format(nf,ng))

        return store 


def pid(prot):
    pid = re.search(options.contig_name_regex, prot).group(1)
    simul = re.search('_(full|[0-9]*)-gene(_[0-9]*)',prot)
    if simul:
        pid = "{}_{}".format(pid,simul.group(1)+simul.group(2)) 
    return pid  


if __name__ == "__name__":
    for opt in [(True,False,True)]:
        print("Caching genomes with options {}".format(opt))
        g = Genomes(*opt)
        del g
