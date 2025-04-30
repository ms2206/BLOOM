"""
Contains the main window of the app.
"""
from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QSplitter, QMessageBox
from PyQt5.QtCore import Qt
from gui.widgets import CustomMenuBar, InputModule, OutputModule, LogBook, Header, AlignmentPopup
import os, signal

class MainWindow(QMainWindow):
    """It is the main window of the app where all the widgets are contained.

    Attributes (and important widgets):
        controller (MainController): controller used to manage everything happening behind the GUI.
        app (QApplication): application where this main window belongs.
        menu_bar (CustomMenuBar): menu bar on top of the window with some functionalities.
        header (Header): title of the app for aesthetic purposes with a button to hide the input module.
        input_module (InputModule): section in the app where the user inputs data.
        output_module (OutputModule): section in the app where the data obtained is displayed.
        logbook (LogBook): section in the app where every action is recorded to provide the user with a
            record of what has been done for reproducibility.
    """
    def __init__(self, controller, app):
        """ Initialises an instance of the class.

        Args:
            controller(MainController)
            app (QApplication)
        """
        super().__init__()
        # Add controller
        self.controller = controller
        self.controller.main_window = self
        # App
        self.app = app
        # Title
        self.setWindowTitle("BLOOM")        
        # Geometry
        self.setGeometry(100, 100, 1500, 900)
        """Add Widgets."""
        # Set up the menu bar
        self.menu_bar = CustomMenuBar(self, self.controller)
        self.setMenuBar(self.menu_bar)   
        # Add header
        self.header = Header(self, self.controller)
        # Add Input Module
        self.input_module = InputModule(self.controller)
        # Add logbook
        self.log_book = LogBook()
        # Add Output module
        self.output_module = OutputModule(self.controller)
        """Splitter to encapsulate widgets."""
        # Create a splitter for logbook and output module
        self.first_splitter = QSplitter(Qt.Vertical)
        ## Add modules to splitter
        self.first_splitter.addWidget(self.output_module)
        self.first_splitter.addWidget(self.log_book)
        ## Set splitter dimensions
        self.first_splitter.setSizes([800, 200])
        self.first_splitter.setStretchFactor(1, 1)
        # Create a second splitter for the input module and the first splitter
        self.second_splitter = QSplitter(Qt.Horizontal)
        ## Add modules to splitter
        self.second_splitter.addWidget(self.input_module)
        self.second_splitter.addWidget(self.first_splitter)
        ## Set splitter dimenstions
        self.second_splitter.setSizes([300, 700])
        self.second_splitter.setStretchFactor(1, 1)
        """Widgets' Layout."""
        # Central widget and layout
        self.central_widget = QWidget()
        self.layout = QVBoxLayout()     
        # Add Header to layout
        self.layout.addWidget(self.header)   
        # Add splitter to layout
        self.layout.addWidget(self.second_splitter)
        # Set layout
        self.central_widget.setLayout(self.layout)
        # Set central widget as center of main window
        self.setCentralWidget(self.central_widget)

    def add_new_barcode_tab(self, tab_name, barcodes_data):
        """Adds a new tab of barcode cards given the name and the list with data for the cards."""
        self.output_module.add_new_barcode_tab(tab_name, barcodes_data)
    
    def write_in_logbook(self, text):
        """Writes in the logbook."""
        self.log_book.log(text)
        # Update the app to show the message immediately
        self.app.processEvents()
    
    def toggle_input_module(self):
        """Hides or shows the input module."""
        if self.input_module.isVisible():
            self.input_module.hide()
        else:
            self.input_module.show()
        self.header.toggle_button_text()
    
    def show_results(self, results):
        """Shows results obtained by BLAST in the results tab."""
        self.output_module.show_results(results)
    
    def show_error(self, message):
        """Creates a pop up window with the error message."""
        QMessageBox.critical(self, "Error", message)
    
    def clear(self):
        """Fully clears the app's widgets"""
        self.input_module.clear()
        self.output_module.clear()
        self.log_book.clear()
    
    def create_popup_from_tree(self, alignment_text, target_name, aligned_name):
        """
        Creates a pop up window with the alignment between the species selected on the taxonomy
        tree and the target species studied.
        """
        # Create popup
        popup = AlignmentPopup(alignment_text, target_name, aligned_name)
        popup.show()
        # Keep a reference to avoid garbage collection
        self._popup = popup
    
    def closeEvent(self, event):
        """
        Kills the app when closing the window. It is not a "delicate" waty to handle it but it 
        is the only thing that kills the taxonomy tree processes running in the back.
        """
        os.kill(os.getpid(), signal.SIGTERM)