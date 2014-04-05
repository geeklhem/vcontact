import pc_matrix
class GenomeCluster(object):
    
    def __init__(self,pcm ):
        """
        pcm : PCMatrix object or tuble features (pandas df),contigs (pandas df),network (sparse matriex).
        """
        if isinstance(pcm,pc_matrix.PCMatrix):
            self.features = pcm.features
            self.contigs = pcm.contigs
            self.network = pcm.network
        else:
            self.features,self.contigs,self.network = pcm

    def membership_matrix(self,network=None,contgs=None,features=None,):
        """Build the membership matrix
        INPUT:
        - network (sparse matrix)
        OUTPUT:
        - (sparse matrix) Membership matrix #target X #clusters, B(g,c) is the proportion 
        of edges weight linking the node g to the cluster C"""

        # Parameters default values
        network = self.network if network == None else network
        
        network = (network - np.min(network)) / (np.max(network) - np.min(network))
        B_sum = np.hstack([network.sum(1)]*len(self.features)) #Sum of weights linking to a given target.

        
        # A list giving the columns of the i-th clusters members :
        clusters = [[targets.index(m) for m in members] for members in clusters2target]
        B_clust = np.hstack([network[:,cols].sum(1) for cols in clusters])


        B = B_clust/B_sum
        B = np.nan_to_num(B)
        return B

        
