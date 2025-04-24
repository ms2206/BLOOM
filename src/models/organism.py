""" Contains Organism class.

"""

import services.bloom_functions as bf


class Organism:
    """Organism to be studied by the user.

    This class stores information about a organism to be studied as well as the list of barcodes
    for each type of barcode found in NCBI. 

    Attributes:
        name (str): scientific name of the organism
        taxid (str): taxonomical ID of the organism
        taxonomy (dict): taxonomical lineage of the organism. Keys: Generic name of the rank.
            Values: [Taxonomical ID, Scientific name, Name of the rank]
    """
    def __init__(self, name):
        """Initialises the instance based on the Organism name."""
        self.name = name
        self.taxid = self._set_taxid()
        self.taxonomy = self._set_taxonomy()
    
    def _set_taxid(self):
        """Sets the taxonomical ID given the organism's name"""
        return bf.get_taxa_id(self.name)

    def _set_taxonomy(self):
        """Sets the taxonomy given the taxonomical ID."""
        return bf.get_taxonomy(self.taxid)
                
                    
        