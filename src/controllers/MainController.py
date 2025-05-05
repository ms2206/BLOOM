"""
Contains the main controller which manages the GUI and the models/services.
"""

import csv
import math
import services.bloom_functions as bf
from models.organism import Organism
from models.barcode import Barcode
from services.alignment import get_seqs_alignment
from models.taxo_tree import TaxoTree
from config.config_manager import ConfigManager
from pathlib import Path


"""Configuration constants."""
config = ConfigManager()
# BLAST percentage threshold
PCT_THRESHOLD = config.get('percentage_threshold')
# Barcode queries 
QUERIES_DICT = config.get_barcodes_queries()


class MainController:
    """Class that manages the actions in the GUI with the functions and classes in the back.

    This class is the most important of the entire app. It communicates the GUI with the models and
    services that operate the calculations and APIs. Only one instance of the class is created when
    declaring the app so that the data stored in its attributes can be accessed by other widgets. 

    Attributes:
        main_window (MainWindow): main window which communicates with this controller
        studied_organism (Organism): contains the biological information of the organism to study  
        barcodes (dict): barcodes obtained from an NCBI search grouped by their types. Keys: barcode types, 
            values: list with instances of the Barcode class                   
        seqs_to_blast (dict): sequences to blast grouped by type of barcode so that there is only one sequence
            per barcode type. Keys: barcode type, values: tuple with the sequence and the accession number. This
            implementation facilitates multiple barcode blast but it is not implemented in the current version.
        blast_results (list): dat obtained by BLAST
        blast_results_unique (list): data obtained by BLAST filtered to get only one datapoint per species
        results_diffs_stats (list): list of frequencies of hits/species with a specific number of differences 
        results_pcts_stats (list): list of frequencies of hits/species within a specific range of identity percentage
        sequence (str): nucleotide sequence selected to BLAST
    """
    def __init__(self):
        """Initialises and instance of the class."""
        self.main_window = None
        self.studied_organism = None       
        self.barcodes = {}                            
        self.seqs_to_blast = {}
        self.blast_results = []
        self.blast_results_unique = []
        self.results_diffs_stats = []
        self.results_pcts_stats = []
        self.sequence = ""

    # PRIVATE METHODS
    def _filter_barcodes(self, barcode_list, barcode_type):
        """Filters the list of barcodes so that all sequences are unique.
        
        Args: 
            barcode_list (list): list of barcodes to filter.
            barcode_type (str): type of barcode to organise them."""
        # Dictionary to store the barcodes with sequnce as key
        unique_barcodes = {}
        filtered_barcodes = {}
        for barcode in barcode_list:
            # Check if the sequence is already in the dict
            if str(barcode) in unique_barcodes:
                for header in barcode.headers:
                    # Append the header of the repeated barcode to the unique barcode
                    unique_barcodes[str(barcode)].add_header(header)
            else:
                # Add the new unique barcode
                unique_barcodes[str(barcode)] = barcode
            filtered_barcodes[barcode_type] = list(unique_barcodes.values())
        return filtered_barcodes
    
    def _analyse_results_by_pcts(self):
        """
        Analyses results to return a list of numbers with the frequency of hits within a specific
        range of identity percentage.
        """
        # Create array of stats for identity percentage
        pcts_all_hits = [0]*(100 - PCT_THRESHOLD + 2)
        pcts_unique = [0]*(100 - PCT_THRESHOLD + 2)
        # Create an array of labels
        labels = ["100%"]
        for i in range(0, len(pcts_all_hits)-2):
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
        
    def _analyse_results_by_diffs(self):
        """
        Analyses results to return a list of numbers with the frequency of hits with a specific
        number of differences.
        """
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

    # PUBLIC METHODS
    ## Searching barcodes 
    def search_barcodes(self, organism, barcode_type, primers):
        """Goes throught the process of creating a new organism and the barcode sequences requested by the user.
        
        Args:
            organism (str): name of the organism to study
            barcode_type (str): type of barcode to search
            primers (tuple): tuple with the forward and reverse primers of the barcode
        """
        try:
            self.studied_organism = Organism(organism)
        except Exception as e:
            self.error_pop_up(f'Error creating the organism: {e}')
            return 1
        else:
            # Get the NCBI query given the barcode
            barcode_query = QUERIES_DICT[barcode_type]
            # Get the list of data from NCBI which includes fasta header and sequence
            sequence_list = bf.get_NCBI_Sequences(barcode_query, self.studied_organism.taxid)
            # Check if there has been an error getting the barcodes
            if not sequence_list[0]:
                # Display pop up error emssage
                self.error_pop_up(f'Error when searching barcodes: {sequence_list[1]}')
                return 1
            else:
                # Create a list to store the barcodes
                barcode_list = []
                for sequence in sequence_list:
                    # Create a barcode for each list given the barcode type, fasta header, the sequence and the primers
                    barcode = Barcode(barcode_type, sequence[0], sequence[1], primers)
                    barcode_list.append(barcode)
                # Filter the barcodes to avoid repetitions
                self.barcodes = self._filter_barcodes(barcode_list, barcode_type)

    def add_new_barcode_tab(self, tab_name):
        """Gets the information to create a new tab with barcodes given the barcode type as tab name."""
        # Get the list of barcodes given the barcode name which is the tab name
        barcodes_to_add = self.barcodes[tab_name]
        barcodes_data = []
        for barcode in barcodes_to_add:
            # Append data to list to send to the barcodes tab
            data = (str(barcode), barcode.headers[0], int(barcode), barcode.primer_coords)
            barcodes_data.append(data)
        # Call main window to create the tab
        self.main_window.add_new_barcode_tab(tab_name, barcodes_data)
    
    ## Selecting barcodes
    def add_sequence_to_blast(self, sequence, barcode_type, acc_number):
        """Adds a new sequence to BLAST to the list given its barcode type.
        
        Args:
            sequence (str): nucleotide sequence to BLAST
            barcode_type (str): type of barcode to organise the sequence to BLAST so that
                there is only one sequence ber type
            acc_number (str): accession number of the sequence selected to keep track of 
                which one was selected"""
        self.seqs_to_blast[barcode_type] = (sequence, acc_number)
    
    def remove_sequence_to_blast(self, type):
        """Remove the sequence from the unchecked barcode from the sequences to BLAST given the barcode type."""
        if type in self.seqs_to_blast.keys():
            del self.seqs_to_blast[type]
    
    ## Getting results from BLAST
    def start_blast(self, blast_mode, rank):
        """BLASTs the sequence selected by the user with the parameters inputted
        
        Args:
            blast_mode (str): option to use Megablast or discontiguous blast
            rank (str): taxonomy rank to narrow down BLAST search
        """
        # check if there is only one sequence selected
        if len(self.seqs_to_blast.keys()) == 1:
            # Get barcode query from tab name and list of queries
            barcode_type = list(self.seqs_to_blast.keys())[0]
            barcode_query = QUERIES_DICT[barcode_type]
            # Get sequence
            self.sequence = self.seqs_to_blast[barcode_type][0]
            # Check if user wants to use megablast
            megablast_use = True if blast_mode == "Megablast" else False
            # Get taxonomy rank name
            rank_name = ""
            for level in self.studied_organism.taxonomy:
                if level['Rank'].lower() == rank.lower():
                    rank_name = level['ScientificName']
                    break
            print(rank_name)
            # Write in logbook
            self.write_in_logbook(f'Starting {blast_mode} search for a {barcode_type} barcode in the {rank} of {self.studied_organism.name}.')
            # Run blast
            blast_results = bf.blast(self.sequence, barcode_query, rank_name, megablast_use)
            if not blast_results:
                self.error_pop_up(f'BLAST error: the hit-list is empty')
                return 1
            if not blast_results[0]:
                self.error_pop_up(f'BLAST error: {blast_results[1]}')
                return 1
            self.blast_results = blast_results
            self.blast_results_unique = bf.filter_data(self.blast_results)
            # Write in logbook
            self.write_in_logbook("BLAST search has finished.")
            # Analyse results
            self.results_pcts_stats = self._analyse_results_by_pcts()
            self.results_diffs_stats = self._analyse_results_by_diffs()
            # Get title for results display
            acc_number = self.seqs_to_blast[barcode_type][1]
            title = f'{blast_mode} search for the {rank} of {self.studied_organism.name} [{barcode_type} gene: {acc_number}]'
            # Display results
            self.main_window.show_results(self.results_diffs_stats, title)
            return None
        else:
            self.error_pop_up("Select ONE barcode to BLAST")
            return None

    def update_results(self, data_type, show_dissimilar):
        """Updates the result visualisatcion based on the selected parameters.
        
        Args:
            data_type (str): type of data to display (identity percentages or number of
                differences)
            show_dissimilar (bool): if true, dissimilar results are shown, otherwise they 
                are not taken into account
        """
        # Check if there are available BLAST results
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

    ## Using taxonomy tree
    def create_tree(self):
        """Creates and displays a taxonomy tree."""
        # Check if there are results available
        if self.blast_results_unique:
            try:
                # Create tree
                tree = TaxoTree(species_info=self.blast_results_unique, 
                                target=self.studied_organism.name,
                                controller=self)
                # Display tree
                tree.show_tree()
                self.write_in_logbook(f'Taxonomy tree created for {self.studied_organism.name}')
            except Exception as e:
                self.error_pop_up(f'Error creating taxonomy tree: {e}')
        else:
            self.error_pop_up("Missing BLAST results")
        
    def show_alignment_from_tree(self, name):
        """
        Gets the alignment from the name of the selected organism on the tree and the sequence blasted.
        Then calls the main window to create a pop-up window with the alignment.
        """
        # Initialise alignment text to be displayed
        alignment_text = ""
        # Get sequence and name from the studied organism
        target = self.studied_organism.name
        target_seq = self.sequence
        # Get all sequences from the name selected in the tree
        for element in self.blast_results:
            if element["Scientific name"].lower() == name.lower():
                sequence = element["Subject sequence"]
                # Sequences in blast results contain gaps as "-" and need to be eliminated
                sequence = sequence.replace("-", "")
                # Add accession number differentiate alignments 
                acc_num = element["Accession Number"]
                alignment = self.align_seqs(target_seq, sequence)
                # Add alignment text to text to display
                alignment_text += acc_num + "\n" + str(alignment) + "\n\n\n"
        self.main_window.create_popup_from_tree(alignment_text, target, name)

    ## Creating CSV
    def write_csv(self, file_path):
        """Creates and writes a csv file given a path to the file."""
        # Get column names from keys in dictionary with results
        column_names = self.blast_results[0].keys()
        try:
            # Write to the file
            with open(file_path, mode='w', newline='', encoding='utf-8') as file:
                writer = csv.DictWriter(file, fieldnames=column_names)
                writer.writeheader()
                writer.writerows(self.blast_results)
            # Write message in the logbook
            self.write_in_logbook(f"Data successfully saved to {file_path}")
        except Exception as e:
            # Error message
            self.error_pop_up(f"Error writing csv file: {e}") 
    
    ## Messaging
    def error_pop_up(self, message):
        """Calls the main window to display a pop-up error message."""
        self.main_window.show_error(message)
    
    def write_in_logbook(self, text):
        """Calls the main window to write a text in the logbook"""
        self.main_window.write_in_logbook(text)

    ## Others
    def clear(self):
        """Resets own data and calls main window to clear all information."""
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
    
    def align_seqs(self, seq1, seq2):
        """Returns the alignment of two sequences.
        
        Args:
            seq1 (str): target or reference sequence.
            seq2 (str): aligned or subejct sequence
        """
        return get_seqs_alignment(seq1, seq2)    

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
            self.error_pop_up("Error getting list of species. File not found")
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
        try:
            with open(file_path, "r") as f:
                return f.read()
        except FileNotFoundError:
            self.error_pop_up("Error getting style sheet. File not found")
            return None