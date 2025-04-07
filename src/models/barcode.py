"""Contains Barcode Class and dictionary for barcode queries.
 
"""

from models.alignment import find_primers 

class Barcode:
    """Barcode card which the user can choose to BLAST.

    This class stores information about a barcode. 

    Attributes:
        type (str): type of barcode
        query (str): query for the gene used to search in NCBI
        header (str): FASTA header from the NCBI entry where the sequence comes from
        primers (tuple): starting primer and ending primer sequences 
        sequence (str): clean sequence with nucleotides before and after limit primers removed
    """

    def __init__(self, type, header, sequence, primers):
        """ Initialises the instance based on a barcode sequence and its characteristics

        Args:
            type (str): defines the type of barcode 
            header (str): defines the header of the NCBI entry 
            sequence (str): raw nucleotide sequence 
            primers (tuple): tuple with the starting and ending primer (in that order) of the barcode
        """
        self.type = type
        self.headers = [header]
        self.primers = primers
        self.primer_coords = []
        self.sequence = self.clean_sequence(sequence)  

    def clean_sequence(self, sequence):
        """Trims a raw sequence to just the primers and the nucleotides in between."""
        coordinates = find_primers(sequence, self.primers)
        starting_offset = coordinates[0][0]
        self.primer_coords = [(a - starting_offset, b - starting_offset) for a, b in coordinates]
        return  sequence[coordinates[0][0]:coordinates[-1][1]+1]
   
    def __str__(self):
        """Returns the sequence so that when converted to string or printed the sequences is passed"""
        return self.sequence
    
    def __int__(self):
        return len(self.headers)
    
    def add_header(self, new_header):
        self.headers.append(new_header)

    def get_headers(self):
        return self.headers
    
    def get_primers_intervals(self):
        return self.primer_coords


"""Dictionary that stores the NCBI query to use for each barcode"""
BARCODE_QUERIES = {
    "trnL": "trnL",
    "trnL-UAA": "trnL-UAA",
    "trnLP6": "trnL",           # The P6 is found inside the trnL loop
    "matK": "matK"
}

    
    



    
    
