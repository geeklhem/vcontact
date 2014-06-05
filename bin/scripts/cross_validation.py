""" Cross validation """
import pc_matrix
import genome_clusters
import options
import os
import logging
import pandas 

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("cv_routine")


output = options.data_folder+"crossvalidations/"

def generate_data(inflations=(2,3,4,5),sig_list=(1,5,10),levels=("family","genus")):
    gcs = {}
    parameters = []
    for sig in sig_list:
        parameters.append({"inflation_pc":2,"sig":sig,"inflation_cc":2})
    for inflation in inflations:
        parameters.append({"inflation_pc":2,"sig":1,"inflation_cc":inflation})
        parameters.append({"inflation_pc":inflation,"sig":1,"inflation_cc":2})

    for p in parameters:
        pcm = pc_matrix.read_pickle("../data/pcm0{}_140526.pkle".format(p["inflation_pc"]))
        name = "gc_mcl{}_sig{}_mcl{}".format(p["inflation_pc"],p["sig"],p["inflation_cc"])

        logger.info(name)
        gc = genome_clusters.GenomeCluster(pcm,
                                           p["inflation_cc"],
                                           p["sig"],
                                           name=name)
        gc.cv_results = {}
        for level in levels:
            gc.cv_results[level] = gc.cross_validation_affiliation(level)
        gcs[name] = gc
    return gcs 

def analyse_data(gcs):

    mean = []
    std = []
    key = []
    for name,gc in gcs.items():
        for level in gc.cv_results.keys():
            for set_ in ["cv_set","train_set","test_set"]:
                key.append(name+"_{}_{}".format(level,set_))
                sname = name.split("_")
                mdf = gc.cv_results[level][set_].mean()
                mdf["inflation_pc"] = int(sname[1][3:])
                mdf["inflation_cc"] = int(sname[3][3:])
                mdf["sig"] = int(sname[2][3:])
                mdf["level"] = level
                mdf["dataset"] = set_
                mean.append(mdf)
                sdf = gc.cv_results[level][set_].std()
                std.append(sdf)

    mean = pandas.DataFrame(mean,key)
    std = pandas.DataFrame(std,key)
    print len(std)
    dataframe = pandas.merge(mean,std,
                             left_index=True, right_index=True,
                             suffixes=["_mean","_std"])
    #print dataframe.columns
    #dataframe.set_index(["level","inflation_pc","sig","inflation_cc"])
    return dataframe


def analyse_data_learning_curve(lc):
    dataframes = []
    for set_,df in lc.items():
        mean = df.groupby("train_size").mean()
        mean["dataset"] = set_
        std = df.groupby("train_size").std()
        dataframes.append(pandas.merge(mean,std,
                                 left_index=True, right_index=True,
                                 suffixes=["_mean","_std"]))
    dataframe = pandas.concat(dataframes)
    return dataframe
