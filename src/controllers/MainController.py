"""Contains the main controller which manages the GUI and the back-end"""

import services.bloom_functions as bf
from models.organism import Organism
from models.barcode import Barcode, BARCODE_QUERIES
from models.alignment import get_seqs_alignment
from pathlib import Path


class MainController:
    """Class that manages the actions in the GUI with the functions and classes in the back-end.

    Attributes (widgets):
        studied_organism(Organism): organism to study
        barcodes(dict): barcodes being studied grouped by type to study
        main_window(MainWindow): main window of the app
        seqs_to_blast(dict): sequences grouped by barcode type to blast
    """
    def __init__(self):
        """Initialises and instance of the class."""
        self.studied_organism = None  # Track the active organism
        self.barcodes = {}
        self.main_window = None
        self.seqs_to_blast = {}

    def set_main_window(self, main_window):
        """Set the main window"""
        self.main_window = main_window

    def create_new_organism(self, name):
        """Creates a new organism and stores it as the active instance."""
        self.studied_organism = Organism(name)

    def add_barcodes(self, barcode_type, primers):
        """Creates the list of instances of barcodes with their primers for each barcode type."""
        # Get the NCBI query given the barcode
        barcode_query = BARCODE_QUERIES[barcode_type]
        # Get the list of data from NCBI which includes fasta header and sequence
        sequence_list = bf.get_NCBI_Sequences(barcode_query, self.studied_organism.get_taxid())
        # Create a list to store the barcodes
        barcode_list = []
        for sequence in sequence_list:
            # Create a barcode for each list given the barcode type, fasta header, the sequence and the primers
            barcode = Barcode(barcode_type, sequence[0], sequence[1], primers)
            barcode_list.append(barcode)
        # Filter the barcodes to avoid repetitions
        self.filter_barcodes(barcode_list, barcode_type)

    def filter_barcodes(self, barcode_list, barcode_type):
        """Filters the list of barcodes so that all sequences are unique."""
        # Dictionary to ster the barcodes with sequnce as key
        unique_barcodes = {}
        for barcode in barcode_list:
            # Check if the sequence is already in the dict
            if str(barcode) in unique_barcodes:
                for header in barcode.get_headers():
                    # Append the header of the repeated barcode to the unique barcode
                    unique_barcodes[str(barcode)].add_header(header)
            else:
                # Add the new unique barcode
                unique_barcodes[str(barcode)] = barcode
            self.barcodes[barcode_type] = list(unique_barcodes.values())
    
    def search_for_barcodes(self, organism, barcode_type, primers):
        """Goes throught the process of creating a new organism and the barcode sequences requested by the user."""
        self.create_new_organism(organism)
        self.add_barcodes(barcode_type, primers)

    def align_seqs(self, seq1, seq2):
        """Returns the alignment of two sequences"""
        return get_seqs_alignment(seq1, seq2)
    
    def add_new_barcode_tab(self, tab_name):
        """Gets the information to create a new tab with barcodes"""
        # Get the list of barcodes given the barcode name which is the tab name
        barcodes_to_add = self.barcodes[tab_name]
        barcodes_data = []
        for barcode in barcodes_to_add:
            # Append data to list to send to the barcodes tab
            data = (str(barcode), barcode.get_headers()[0], int(barcode), barcode.get_primers_intervals())
            barcodes_data.append(data)
        # Call main window to create the tab
        self.main_window.add_new_barcode_tab(tab_name, barcodes_data)
    
    def write_in_logbook(self, text):
        """Calls the main window to write in the logbook"""
        self.main_window.write_in_logbook(text)
    
    def load_species_list(self):
        """Loads the list of species names from data/species_list.txt."""
        # Get the path to the current script (MainController.py)
        current_file = Path(__file__)
        # Navigate to the target file relative to the current file
        file_path = current_file.parent.parent / "assets" / "data" / "species_list.txt"
        # Read the contents
        try:
            with file_path.open("r", encoding="utf-8") as f:
                return [line.strip() for line in f if line.strip()]
        except FileNotFoundError:
            print(f"File not found: {file_path}")
            return None
    
    def get_icon(self):
        """Returns the path to the BLOOM icon"""
        # Get the path to the current script (MainController.py)
        current_file = Path(__file__)
        # Navigate to the target file relative to the current file
        file_path = current_file.parent.parent / "assets" / "icons" / "bloom_icon.ico"
        return str(file_path)

    def get_logo(self):
        """Returns the path to the BLOOM png logo"""
        # Get the path to the current script (MainController.py)
        current_file = Path(__file__)
        # Navigate to the target file relative to the current file
        file_path = current_file.parent.parent / "assets" / "images" / "bloom_logo.png"
        return str(file_path)
    
    def load_style(self):
        """Loads the style sheet."""
        # Get the path to the current script (MainController.py)
        current_file = Path(__file__)
        # Navigate to the target file relative to the current file
        file_path = current_file.parent.parent / "assets" / "styles" / "style.qss"
        with open(file_path, "r") as f:
            return f.read()
    
    def add_sequence_to_blast(self, sequence, type):
        self.seqs_to_blast[type] = sequence
        print(self.seqs_to_blast)

    

