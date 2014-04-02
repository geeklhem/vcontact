class Genomes(object):
    """Data about reference genomes, object built from the refseq gbff file.
    """
    
    def __init__(self, fi ):
        """
        
        Arguments:
        - `fi`:
        """
        self._fi = fi
        self.names = None
        self.prot2genome = None
        self.taxonomy = None
        self.prot2func = None

    def load_geome(self,fi=None):
        """
        Load the genomes names and associated proteins
        Arguments:
        - `self`:
        - `fi`:
        """

    def load_taxonomy(self,fi=None):
        """
        
        Arguments:
        - `self`:
        - `fi`:
        """
        
