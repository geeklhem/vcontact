import pandas
import re
import os

if __name__ == "__main__":
    store = pandas.HDFStore(os.path.expanduser("~/data/cache/tara_reads_per_contigs.h5"))
    if not "normalised" in store:
        reads = pandas.read_csv(os.path.expanduser("~/data/tara/READS_PER_CONTIGS.TAB"),sep="\t")

        reads.columns = [re.search("/([0-9A-Z]*).LIST",x).groups()[0]
                         if bool(re.search("/([0-9A-Z]*).LIST",x)) else x
                         for i,x in enumerate(reads.columns)]
        reads = reads.drop(["Unnamed: 44","Unnamed: 45"], axis=1,)
        reads = reads.set_index("name")
        print reads.head()

        colums_sum = reads.values.sum(0) 

        for i,cs in enumerate(colums_sum):
            reads.ix[:,i] = reads.ix[:,i] / float(cs) 


        from Bio import SeqIO
        for seq_record in SeqIO.parse("ALL_CONTIGS_NO_QC.FNA", "fasta"):
            reads.ix[seq_record.id,:] = reads.ix[seq_,i] / float(len(seq_record))



        store.append('normalised',reads,format='table')

    if not "c10" in store :

        with open("~/data/tara/contigs10.txt","rb") as f:
            c10 = [x.strip for x in f]

        filtered = store.normalised.loc[c10,:]
        store.append("c10",
                     filtered,
                     format="table")
