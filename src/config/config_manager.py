"""Contaings path to config file and ConfigManager class.

"""

import os
import yaml

# Path to config file
CONFIG_PATH = os.path.join(os.path.dirname(__file__), 'config.yaml')


class ConfigManager:
    """Class to load and get values stored in config.yaml file."""

    def __init__(self):
        """Creates instance by loading data from file"""
        self.config = self.load_config()

    def load_config(self):
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