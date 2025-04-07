""" Contains Organism class.

"""

import services.bloom_functions as bf

from models.barcode import Barcode


class Organism:
    """Organism to be studied by the user.

    This class stores information about a organism to be studied as well as the list of barcodes
    for each type of barcode found in NCBI. 

    Attributes:
        name (str): scientific name of the organism
        taxid (str): taxonomical ID of the organism
        taxonomy (dict): taxonomical lineage of the organism. Keys: Generic name of the rank.
            Values: [Taxonomical ID, Scientific name, Name of the rank]
        barcodes (dict): list of barcode sequences classified by type found for the organism.
            Keys: type of barcode. Values: instances of the class Barcode
    """
    def __init__(self, name):
        """Initialises the instance based on the Organism name."""
        self.name = name
        self.taxid = self.set_taxid()
        self.taxonomy = self.set_taxonomy()
        self.barcodes = {}
    
    def set_taxid(self):
        """Sets the taxonomical ID given the organism's name"""
        return bf.get_taxa_id(self.name)

    def set_taxonomy(self):
        """Sets the taxonomy given the taxonomical ID."""
        return bf.get_taxonomy(self.taxid)
    
    def eval_barcode(self, header, rank):
        for key in self.barcodes:
            for barcode in self.barcodes[key]:
                if header == barcode.get_header():
                    print(barcode)
                    return bf.blast(barcode, rank)
                
    def get_taxonomy(self):
        return self.taxonomy
    
    def get_taxid(self):
        return self.taxid
                    
        