""" Loading tools """
import pandas
import logging
import scipy.sparse as sparse
import os
import cPickle as pickle


logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

def csv(contigs_fi, pcs_fi, pcprofiles_fi, folder, name, force=False):
    """
    Import csv file and build the info table and pc-profiles matrix.
    Save everything into h5 and pickle files. 

    Args:
        contigs_fi (str): path to the csv file containing the contigs
        pcs_fi (str): path to the csv file containing the pcs
        pcprofiles_fi (str): path to the csv file containing the pc_profiles
        folder (str): output folder path 
        name (str): experiment name 
        force (str): overwrite existing files
    """
    store = pandas.HDFStore(folder+name+".h5")

    if "contigs" not in store or force:
        if contigs_fi != None:
            contigs = pandas.read_csv(contigs_fi, sep=None)
            logger.debug("Read {} entries from {}".format(len(contigs),
                                                           contigs_fi))
            contigs.index.name = "pos"
            contigs.reset_index(inplace=True)
            store.append('contigs', contigs, format='table')
        else:
            raise ValueError("Need contig info file")

    pcs = None 
    if ("pcs" not in store or force) and pcs_fi != None:
        pcs = pandas.read_csv(pcs_fi, sep=None)
        logger.debug("Read {} entries from {}".format(len(pcs),
                                                       pcs_fi))
    elif "pcs" in store:
        pcs = store.pcs
        
    if not os.path.exists(folder+"profiles.pkle") or force:
        if pcprofiles_fi is not None:
            profiles = pandas.read_csv(pcprofiles_fi,sep=None)

            if pcs is None:
                pcs = pandas.DataFrame(profiles.pc_id.drop_duplicates())
                pcs.columns = ["id"]

            # Filtering the PC profiles that appears only once
            before_filter = len(profiles)
            cont_by_pc = profiles.groupby("pc_id").count().contig_id.reset_index()

            # get the number of contigs for each pcs and add it to the dataframe
            cont_by_pc.columns = ["pc_id","nb_contigs"]
            pcs = pandas.merge(pcs,cont_by_pc,left_on="id",right_on="pc_id",how="left")
            pcs.fillna({"no_contigs":0},inplace=True)

            # Drop the pcs that <= 1 contig from the profiles.
            pcs = pcs.query("nb_contigs>1")
            at_least_a_cont = cont_by_pc.query("nb_contigs>1")
            profiles = profiles.query("pc_id in at_least_a_cont.pc_id")

            logger.info(("Read {} entries (dropped {} "
                          "singletons) from {}").format(len(profiles),
                                                        before_filter-len(profiles),
                                                        pcprofiles_fi))
        
            pcs = pcs.reset_index(drop=True)
            pcs.index.name = "pos"
            pcs = pcs.reset_index()
            
            matrix,singletons = _matrix(profiles, store.contigs, pcs)

            #save
            profiles = {"matrix":matrix,"singletons":singletons}
            store.append('pcs', pcs, format='table')
            with open(folder+"profiles.pkle","w") as f:
                pickle.dump(profiles,f)
        else:
            raise ValueError("Need profiles file")
        
    else: #If the pickle file exist
        with open(folder+"profiles.pkle","r") as f:
            profiles = pickle.load(f)
        matrix,singletons = profiles["matrix"],profiles["singletons"]

        
    logger.info(("{} contains : \n {:10} contigs, \n {:10} protein-"
                 "clusters, \n {:10} singletons \n {:10} profile-"
                 "entries.").format(folder, len(store.contigs),
                                    len(store.pcs), singletons.sum(),
                                    matrix.getnnz()))

        
    if (len(store.contigs) != matrix.shape[0]
        or len(store.contigs) != singletons.shape[0]
        or len(store.pcs) != matrix.shape[1]):
        logger.error("profile matrix: {}, singletons matrix {}".format(matrix.shape,
                                                                       singletons.shape))
        logger.debug(matrix.todense())
        logger.debug(singletons.todense())
        raise ValueError("Number of contigs or pc non consistent with the profiles data")

    return store.pcs, store.contigs, profiles 

def _matrix(profiles,contigs,pcs):
    """
    Build the pc profiles matrices (shared & singletons) from dataframes.

    Args:
        profiles (dataframe): required fields are contig_id and pc_id.
        contigs (dataframe): contigs info, required field are proteins, pos and id.
        pcs (dataframe): pcs info, required field are pos and id. 

    Returns:
        (tuple of sparse matrix): Shared PCs and singletons matrix.
    """
    
    pc_by_cont = profiles.groupby("contig_id").count().pc_id
    pc_by_cont = pc_by_cont.reset_index()
    pc_by_cont = pandas.merge(contigs.sort("pos").loc[:,["pos","id","proteins"]],
                              pc_by_cont, how="left",
                              left_on="id",right_on="contig_id").fillna(0)
    singletons = (pc_by_cont.proteins - pc_by_cont.pc_id).values
    singletons = sparse.lil_matrix(singletons).transpose()

    # Matrix 
    profiles.index.name = "pos"
    profiles.reset_index(inplace=True)

    profiles = pandas.merge(profiles,pcs.loc[:,["id","pos"]],
                            left_on="pc_id", right_on="id",
                            how="left",suffixes=["","_pc"])

    profiles = pandas.merge(profiles,contigs.loc[:,["id","pos"]],
                            left_on="contig_id", right_on="id",
                            how="left",suffixes=["","_contig"])
    
    profiles = profiles.loc[:,["pos_contig", "pos_pc"]]
    matrix = sparse.coo_matrix(([1]*len(profiles), (zip(*profiles.values))),
                               shape=(len(contigs), len(pcs)),
                               dtype="bool")

    return matrix.tocsr(), singletons.tocsr()
