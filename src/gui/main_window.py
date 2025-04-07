"""Contains the main window of the app."""
from PyQt5.QtWidgets import QMainWindow, QVBoxLayout, QWidget, QSplitter
from PyQt5.QtCore import Qt
from gui.widgets import CustomMenuBar, InputModule, OutputModule, LogBook, Header


class MainWindow(QMainWindow):
    """Main window of the app.

    Attributes (widgets):
        
    """
    def __init__(self, controller):
        """ Initialises an instance of the class.

        Args:
            controller(MainController)
        """
        super().__init__()
        # Add controller
        self.controller = controller
        self.controller.set_main_window(self)
        # Title
        self.setWindowTitle("BLOOM")        
        # Geometry
        self.setGeometry(100, 100, 1500, 900)
        """Add Widgets."""
        # Set up the menu bar
        self.menu_bar = CustomMenuBar(self)
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
    
    def toggle_input_module(self):
        """Hides or shows the input module."""
        if self.input_module.isVisible():
            self.input_module.hide()
        else:
            self.input_module.show()
        self.header.toggle_button_text()