"""
Contains Barcode Class and dictionary for barcode queries. 
"""

from services.alignment import find_primers 

class Barcode:
    """Barcode card which the user can choose to BLAST.

    This class stores information about a barcode being studied: the trimmed nucleotide sequence with 
    all nucleotides outside of the primers deleted, a list with all the fasta headers of barcodes with
    the same nucleotide sequence, the barcode type, the forward and reverse primers and the coordinates
    of the primers so that they can be used to highlight matching bases. 

    Attributes:
        type (str): type of barcode
        headers (str): list of string FASTA headers from the NCBI entries of all barcodes with the 
            same nucleotide sequence. Originally it is only one element and it increases when detecting
            for duplicates.
        primers (tuple): starting primer and ending primer sequences 
        sequence (str): clean sequence with nucleotides before and after limit primers removed
    """

    def __init__(self, type, header, sequence, primers):
        self.type = type
        self.headers = [header]
        self.primers = primers
        self.primer_coords = []
        self.sequence = self._clean_sequence(sequence)  

    def _get_trim_limits(self, intervals):
        """
        Returns the coordinates where the sequence must be trimmed to delete all nucleotides outside
        the primers.
        """
        all_coordinates  = [num for sublist in intervals for tup in sublist for num in tup]
        return [min(all_coordinates), max(all_coordinates)]
    
    def __str__(self):
        """Returns the sequence so that when converted to string or printed the sequences is passed."""
        return self.sequence
    
    def __int__(self):
        """Returns the number of fasta headers stored inside the object when converted to integer."""
        return len(self.headers)
    
    def _clean_sequence(self, sequence):
        """Trims a raw sequence to just the primers and the nucleotides in between."""
        # Get intervals where primers are located
        intervals = find_primers(sequence, self.primers)
        # Get trim limits from primers positions
        trim_limits = self._get_trim_limits(intervals)
        # Get trimming offset to accurately locate primer coordinates in trimmed sequence
        start_off = trim_limits[0]
        # Get primer coordinates in trimmed sequence
        self.primer_coords = [[(x - start_off, y - start_off) for (x, y) in sublist] for sublist in intervals]
        # Return trimmed sequence
        return  sequence[trim_limits[0]:trim_limits[1]+1]
   
    def add_header(self, new_header):
        """Adds a header to the list of headers stored in the object"""
        self.headers.append(new_header)