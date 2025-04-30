"""
Contaings path to config file and ConfigManager class.
"""

import os
import yaml

# Path to config file
CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config.yaml')


class ConfigManager:
    """Class to load and get values stored in config.yaml file."""

    def __init__(self):
        """Creates instance by loading data from file"""
        self.config = self._load_config()

    def _load_config(self):
        """Loads the user configuration."""
        with open(CONFIG_PATH, 'r') as file:
            return yaml.safe_load(file)

    def get(self, key_path, default=None):
        """Access config values using dot notation."""
        keys = key_path.split('.')
        value = self.config
        try:
            for key in keys:
                value = value[key]
            return value
        except KeyError:
            return default
    
    def get_barcodes_names(self):
        """Acces the barcode structure to retrieve the names of all barcodes stored."""
        barcodes = self.get('barcodes', [])
        return list(barcode['name'] for barcode in barcodes)
    
    def get_barcodes_queries(self):
        """Access the barcode structure to return a dict with the queries for each barcode"""
        barcodes = self.get('barcodes', [])
        queries_dict = {}
        for barcode in barcodes:
            name = barcode['name']
            primers = barcode['query']
            queries_dict[name] = primers
        return queries_dict

    def get_barcodes_primers(self):
        """Access the barcode structure to return a dict with the primers for each barcode"""
        barcodes = self.get('barcodes', [])
        primers_dict = {}
        for barcode in barcodes:
            name = barcode['name']
            primers = (barcode['fprimer'], barcode['rprimer'])
            primers_dict[name] = primers
        return primers_dict