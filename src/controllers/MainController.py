"""Contains the main controller which manages the GUI and the back-end"""

import csv
import math
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
        self.blast_results = {}
        self.sequence = ""

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
        self.barcodes = self.filter_barcodes(barcode_list, barcode_type)

    def filter_barcodes(self, barcode_list, barcode_type):
        """Filters the list of barcodes so that all sequences are unique."""
        # Dictionary to ster the barcodes with sequnce as key
        unique_barcodes = {}
        filtered_barcodes = {}
        for barcode in barcode_list:
            # Check if the sequence is already in the dict
            if str(barcode) in unique_barcodes:
                for header in barcode.get_headers():
                    # Append the header of the repeated barcode to the unique barcode
                    unique_barcodes[str(barcode)].add_header(header)
            else:
                # Add the new unique barcode
                unique_barcodes[str(barcode)] = barcode
            filtered_barcodes[barcode_type] = list(unique_barcodes.values())
        return filtered_barcodes
    
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
    
    def remove_sequence_to_blast(self, type):
        if type in self.seqs_to_blast.keys():
            del self.seqs_to_blast[type]
    
    def start_blast(self, blast_mode, rank):
        results = ()
        if len(self.seqs_to_blast.keys()) == 1:
            # Get barcode query from tab name and list of queries
            barcode_type = list(self.seqs_to_blast.keys())[0]
            barcode_query = BARCODE_QUERIES[barcode_type]
            # Get sequence
            self.sequence = self.seqs_to_blast[barcode_type]
            # Check if user wants to use megablast
            megablast_use = True if blast_mode == "Megablast" else False
            # Get taxonomy rank name
            rank_name = ""
            for level in self.studied_organism.taxonomy:
                if level['Rank'] == rank:
                    rank_name = level['ScientificName']
                    break
            # Write in logbook
            self.write_in_logbook(f'Starting BLAST search for a {barcode_type} barcode in the {rank} of {self.studied_organism.get_name()}.')
            # Run one_blast
            self.blast_results = bf.blast2(self.sequence, barcode_query, rank_name, megablast_use)
            # Write in logbook
            self.write_in_logbook("BLAST search has finished.")
            # Analyse results
            results = self.analyse_results_by_diffs("Real")
        self.main_window.show_results(results, "diffs")
            
    def analyse_results_by_diffs(self, mode):
        target_key = "Real num. diffs" if mode == "Real" else "Aligned num. diffs"
        # Get differences for all of BLAST hits
        diffs_all_hits = [0, 0, 0, 0, 0, 0, 0] # [0, 1, 2, 3, 4, 5, more than 5] differences
        for datapoint in self.blast_results:
            index = datapoint[target_key] if datapoint[target_key] <= 5 else 6
            diffs_all_hits[index] += 1
        # Get differences taking into account unique species (lowest number of diffs)
        results_species = bf.filter_data(self.blast_results, target_key)   
        diffs_unique_species = [0, 0, 0, 0, 0, 0, 0] # [0, 1, 2, 3, 4, 5, more than 5] differences
        for key in results_species.keys():
            index = results_species[key] if results_species[key] <= 5 else 6    
            diffs_unique_species[index] += 1
        return (diffs_all_hits, diffs_unique_species)

    def analyse_results_by_pct(self, mode):
        target_key = "Real sim. pct" if mode == "Real" else "Aligned sim. pct"
        # Analyse data for all hits
        pcts_all_hits = [0, 0, 0, 0, 0, 0, 0] # [100, 100-99, 99-98, 98-95, 95-90, 90-80, <80] %
        for datapoint in self.blast_results:
            diffs = datapoint[target_key]
            index = 6
            if diffs == 100:
                index = 0
            elif 99 <= diffs < 100:
                index = 1
            elif 98 <= diffs < 99:
                index = 2
            elif 95 <= diffs < 98:
                index = 3
            elif 90 <= diffs < 95:
                index = 4
            elif 80 <= diffs < 90:
                index = 5
            pcts_all_hits[index] += 1   
        # Analyse data for unique species
        results_species = bf.filter_data(self.blast_results, target_key)   
        pcts_species = [0, 0, 0, 0, 0, 0, 0] 
        for key in results_species.keys():
            diffs = results_species[key]
            index = 6
            if diffs == 100:
                index = 0
            elif 99 <= diffs < 100:
                index = 1
            elif 98 <= diffs < 99:
                index = 2
            elif 95 <= diffs < 98:
                index = 3
            elif 90 <= diffs < 95:
                index = 4
            elif 80 <= diffs < 90:
                index = 5
            pcts_species[index] += 1   
        return (pcts_all_hits, pcts_species)

    def update_results(self, data_type, data_group):
        results = None
        mode = "Real" if data_group == "Complete sequences" else "Aligned"
        if data_type == "Number of differences":
            results = self.analyse_results_by_diffs(mode)
            type = "diffs"
        else:
            results = self.analyse_results_by_pct(mode)
            type = "pcts"
        self.main_window.show_results(results, type)
        
    def write_csv(self, file_path):
        column_names = self.blast_results[0].keys()
        try:
            # Write to the file
            with open(file_path, mode='w', newline='', encoding='utf-8') as file:
                writer = csv.DictWriter(file, fieldnames=column_names)
                writer.writeheader()
                writer.writerows(self.blast_results)
            print(f"Data successfully saved to {file_path}")
        except Exception as e:
            print(f"Error writing file: {e}")

    