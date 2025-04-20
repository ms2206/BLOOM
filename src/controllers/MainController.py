"""Contains the main controller which manages the GUI and the back-end"""

import csv
import math
import services.bloom_functions as bf
from models.organism import Organism
from models.barcode import Barcode
from models.alignment import get_seqs_alignment
from config.config_manager import ConfigManager
from pathlib import Path


"""Configuration constants."""
config = ConfigManager()
# BLAST percentage threshold
PCT_THRESHOLD = config.get('percentage_threshold')
# Barcode queries 
QUERIES_DICT = config.get_barcodes_queries()


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
        self.blast_results = []
        self.blast_results_unique = []
        self.results_diffs_stats = []
        self.results_pcts_stats = []
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
        barcode_query = QUERIES_DICT[barcode_type]
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
        if len(self.seqs_to_blast.keys()) == 1:
            # Get barcode query from tab name and list of queries
            barcode_type = list(self.seqs_to_blast.keys())[0]
            barcode_query = QUERIES_DICT[barcode_type]
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
            # Run blast
            self.blast_results = bf.blast(self.sequence, barcode_query, rank_name, megablast_use)
            self.blast_results_unique = bf.filter_data(self.blast_results)
            # Write in logbook
            self.write_in_logbook("BLAST search has finished.")
            # Analyse results
            self.results_pcts_stats = self.analyse_results_by_pcts()
            self.results_diffs_stats = self.analyse_results_by_diffs()
            self.main_window.show_results(self.results_diffs_stats)
        else:
            self.error_pop_up("Select ONE barcode to BLAST")

    def update_results(self, data_type, show_dissimilar):
        if not self.results_diffs_stats or not self.results_pcts_stats:
            self.error_pop_up("Missing BLAST results")
            return
        # Create empty results holder
        results = None
        # Check if the user wants differences or percentages
        if data_type == "Number of differences":
            # Don't pass the last element of each tuple in results if user does not want
            # to see dissimilar hits
            if show_dissimilar:
                results = self.results_diffs_stats
            else:
                results = tuple(element[:-1] for element in self.results_diffs_stats)
        else:
            if show_dissimilar:
                results = self.results_pcts_stats
            else:
                results = tuple(element[:-1] for element in self.results_pcts_stats)
        # Return results
        self.main_window.show_results(results)
        
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

    def analyse_results_by_pcts(self):
        # Create array of stats for identity percentage
        pcts_all_hits = [0]*(100 - PCT_THRESHOLD + 2)
        pcts_unique = [0]*(100 - PCT_THRESHOLD + 2)
        # Create an array of labels
        labels = ["100%"]
        for i in range(1, len(pcts_all_hits)-1):
            labels.append(str(100-i) + "-" + str(100-i-1) + "%")
        labels.append("<" +  str(PCT_THRESHOLD) + "%")
        # Calculate stats for all hits
        for datapoint in self.blast_results:
            data_pct = datapoint["Identity percentage"]
            index = len(pcts_all_hits) - 1
            for pct in range(PCT_THRESHOLD, 101):
                if data_pct >= pct:
                    index -= 1
                else:
                    break
            pcts_all_hits[index] += 1
        # Calculate stats for unique species
        for datapoint in self.blast_results_unique:
            data_pct = datapoint["Identity percentage"]
            index = len(pcts_unique) - 1
            for pct in range(PCT_THRESHOLD, 101):
                if data_pct >= pct:
                    index -= 1
                else:
                    break
            pcts_unique[index] += 1
        # Return results
        return (labels, pcts_all_hits, pcts_unique)
        
    def analyse_results_by_diffs(self):
        # Calculate the maximum number of differences given threshold percentage
        max_diffs = int(math.ceil(len(self.sequence)*(100 - PCT_THRESHOLD)/100))
        # Create array to count the number of hits or species with a specific number of hits
        diffs_all_hits = [0]*(max_diffs + 2)
        diffs_unique = [0]*(max_diffs + 2)
        # Create an array of labels
        labels = list(str(i) for i in range(max_diffs+1))
        labels.append(">" + str(max_diffs))
        # Calculate stats for all BLAST hits data
        for datapoint in self.blast_results:
            index = datapoint["Differences"] if datapoint["Differences"] <= max_diffs else (max_diffs+1)
            diffs_all_hits[index] += 1
        # Calculate stats for unique species data
        for datapoint in self.blast_results_unique:
            index = datapoint["Differences"] if datapoint["Differences"] <= max_diffs else (max_diffs+1)
            diffs_unique[index] += 1  
        # Return results
        return (labels, diffs_all_hits, diffs_unique) 
    
    def error_pop_up(self, message):
        self.main_window.show_error(message)
    
    def clear(self):
        # Clear own data
        self.studied_organism = None  # Track the active organism
        self.barcodes = {}
        self.seqs_to_blast = {}
        self.blast_results = []
        self.blast_results_unique = []
        self.results_diffs_stats = []
        self.results_pcts_stats = []
        self.sequence = ""
        # Clear gui data
        self.main_window.clear()
