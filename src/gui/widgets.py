"""Contains the declaration of all custom widgets as classes."""
from datetime import datetime
from PyQt5.QtCore import Qt, QStringListModel
from PyQt5.QtGui import QTextCursor, QTextCharFormat, QColor
from PyQt5.QtWidgets import *


"""Constants."""
BARCODE_LIST = ["trnL", "trnL-UAA", "trnLP6"]
BLAST_MODES = ["Megablast", "Discontiguous megablast", "Blastn"]
TAXONOMY_RANKS = ["Genus", "Family", "Order", "Kingdom"]
MAX_THRESHOLD = 100 # For similarity
MIN_THRESHOLD = 97


class CustomMenuBar(QMenuBar):
    """Menu bar at the top of the app

    This class contains the menus and submenus at the top of the app 

    Attributes (widgets):
        file_menu: contains the actions to close the app and generate report (future)
        settings_menu: contains the actions to change default parameters (future)
        help_menu: contains the action to provide the user with help (future)
    """
    def __init__(self, parent):
        super().__init__(parent)
        """ Initialises the instance 

        Args:
            parent (QMainWindow): main window of the app where the menu belongs
        """
        # File Menu
        self.file_menu = self.addMenu("File")
        exit_action = QAction("Exit", self)         # Close the app
        exit_action.triggered.connect(parent.close)
        self.file_menu.addAction(exit_action)
        # Settings Menu
        self.settings_menu = self.addMenu("Settings")
        undo_action = QAction("Undo", self)     # Currently does nothing
        redo_action = QAction("Redo", self)     # Currently does nothing
        self.settings_menu.addAction(undo_action)
        self.settings_menu.addAction(redo_action)
        # Help Menu
        self.help_menu = self.addMenu("Help")


class Header(QWidget):
    """Header with the title of the app and button to close input module.

    Attributes (widgets):
        toggle_button(QToolButton): hides and shows the input module
        title(QLabel): full title of the app
    """
    def __init__(self, main_window, controller):
        """ Initialises the instance 

        Args:
            main_window (QMainWindow): main window of the app where the title belongs
            controller (MainController): controller to take the actions of the signals
        """
        super().__init__()
        # Add main window
        self.main_window = main_window
        # Add controller
        self.controller = controller
        """ Add widgets """
        # Add Button
        self.toggle_button = QToolButton()
        self.toggle_button.setText("☰")
        self.toggle_button.setFixedSize(40, 40) 
        self.toggle_button.clicked.connect(self.toggle_input_module)
        # Add Label
        self.title = QLabel("Barcode Linker for Organism and Ontology Mapping")
        self.title.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
        self.title.setObjectName("title")   # Sets a custom style
        """ Widget's layout """
        # Define layout
        self.setFixedHeight(70)
        self.layout = QHBoxLayout()
        self.layout.setSpacing(10)
        # Add toggle button
        self.layout.addWidget(self.toggle_button)
        # Add title
        self.layout.addWidget(self.title)
        # Set layout
        self.setLayout(self.layout)

    def toggle_input_module(self):
        """Commands the main window to hide or show the input module."""
        self.main_window.toggle_input_module()
    
    def toggle_button_text(self):
        """Changes the aspect of the button."""
        if self.toggle_button.text() == "☰":
            self.toggle_button.setText(">")
        else:
            self.toggle_button.setText("☰")


class SearchTab(QFrame):
    """Tab in the input module to define the barcode and organism to search.

    This class inherits from QFrame. It incorporates the entry widgets to select which
    organism which autocompletes with NCBi data and a barcode dropdown to select which
    barcode to search for. There are also boxes to input the primers' sequences.

    Attributes (widgets):
        species_input(QLineEdit): to enter the name of the organism
        species_model(QStringListModel): list of organism names to choose from
        completer(QCompleter): completer for the organism name
        barcode_dropdown(QComboBox): to choose a barcode
        f_primer_input(QTextEdit): to enter the forward primer sequence
        r_primer_input(QTextEdit): to enter the reverse primer sequence
        search_button(QPushButton): button to search the selected parameters
    """
    def __init__(self, controller):
        """ Initialises the instance 

        Args:
            controller (MainController): controller to take the actions of the signals
        """
        super().__init__()
        # Add controller
        self.controller = controller
        # Load list of species
        self.full_species_list = self.controller.load_species_list()
        """ Add Widgets """
        # Input for name of species
        self.species_input = QLineEdit()
        self.species_input.setPlaceholderText("Enter species name...")
        # Completer for the list of species
        self.species_model = QStringListModel()
        self.completer = QCompleter(self.species_model, self)
        self.completer.setCaseSensitivity(Qt.CaseInsensitive)
        self.completer.setFilterMode(Qt.MatchContains)
        self.completer.setCompletionMode(QCompleter.PopupCompletion)
        # Popup for the list of species
        popup = QListView()
        self.completer.setPopup(popup)
        # Add completer and popup to species_input and add signal to suggest names
        self.species_input.setCompleter(self.completer)
        self.species_input.textEdited.connect(self.update_species_completer)
        # Barcode selection dropdown
        self.barcode_dropdown = QComboBox()
        self.barcode_dropdown.addItems(BARCODE_LIST)
        # Input for forward primer
        self.f_primer_input = QTextEdit()
        self.f_primer_input.setPlaceholderText("Enter forward primer sequence...")
        self.f_primer_input.setFixedHeight(50)
        self.f_primer_input.setObjectName("forwardPrimer")
        # Input for reverse primer
        self.r_primer_input = QTextEdit()
        self.r_primer_input.setPlaceholderText("Enter reverse primer sequence...")
        self.r_primer_input.setFixedHeight(50)
        self.r_primer_input.setObjectName("reversePrimer")
        # Search Button
        self.search_button = QPushButton("🔍 Search")
        self.search_button.clicked.connect(self.search_button_action)
        # Labels
        self.species_label = QLabel("Species")
        self.barcode_label = QLabel("Barcode")
        self.f_primer_label = QLabel("Forward primer")
        self.r_primer_label = QLabel("Reverse primer")
        """ Widgets' layout """
        # Define layour
        self.layout = QVBoxLayout()
        # Add species input
        self.layout.addSpacing(30) 
        self.layout.addWidget(self.species_label)
        self.layout.addWidget(self.species_input)
        self.layout.addSpacing(30) 
        # Add barcode input
        self.layout.addWidget(self.barcode_label)
        self.layout.addWidget(self.barcode_dropdown)
        self.layout.addSpacing(30) 
        # Add primers input
        self.layout.addWidget(self.f_primer_label)
        self.layout.addWidget(self.f_primer_input)
        self.layout.addSpacing(10) 
        self.layout.addWidget(self.r_primer_label)
        self.layout.addWidget(self.r_primer_input)
        # Add space between button and other widgets
        self.layout.addStretch()
        # Add search button
        self.layout.addWidget(self.search_button)
        self.layout.addSpacing(30)     
        # Set layout
        self.setLayout(self.layout)
    
    def search_button_action(self):
        """Gathers the data input by the user and commands the controller to search."""
        # Block organism name button so user can't change it
        self.species_input.setEnabled(False)
        # Get organism, barcode type and primers chosen by the user
        organism = self.species_input.text()
        barcode_type = BARCODE_LIST[self.barcode_dropdown.currentIndex()]
        f_primer = self.f_primer_input.toPlainText()
        r_primer = self.r_primer_input.toPlainText()
        primers = (f_primer, r_primer)
        # Logbook
        self.controller.write_in_logbook(f'Fetching {barcode_type} sequences for {organism}.')
        # Call controller to get barcodes
        self.controller.search_for_barcodes(organism, barcode_type, primers)
        # Call controller to create new tab with barcodes
        self.controller.add_new_barcode_tab(barcode_type)
    
    def update_species_completer(self):
        """Updates the list of suggested species given what the user has written."""
        if len(self.species_input.text()) < 3:
            self.species_model.setStringList([])  # Wait for user to write at least 3 characters
            return
        # Filter list and use the nuew filtered list instead
        filtered = [s for s in self.full_species_list if self.species_input.text().lower() in s.lower()]
        self.species_model.setStringList(filtered)


class BlastTab(QFrame):
    """Tab in the input module to define the parameters for BLAST and to generate the tree.

    This class inherits from QFrame. It incorporates the entry widgets to select which
    BLAST mode to use, the taxonomy rank to search and the percentage of identity from 
    which two sequences are considered identical. It also incorporates a button to use
    BLAST and a button to create a phylogenetic tree.

    Attributes (widgets):
        blast_mode_dropdown(QComboBox): to select the blast mode
        rank_dropdown(QComboBox): to select the taxonomy rank to analyse
        identity_slider(QSlider): to select the similarity threshold
        blast_button(QPushButton): to use BLAST with the parameters selected
        tree_button(QPushButton): to generate the tree
    """
    def __init__(self):
        """Creates an instance of the class."""
        super().__init__()
        """ Add Widgets """
        # Input blast mode
        self.blast_mode_dropdown = QComboBox()
        self.blast_mode_dropdown.addItems(BLAST_MODES)
        # Input taxonomy rank
        self.rank_dropdown = QComboBox()
        self.rank_dropdown.addItems(TAXONOMY_RANKS)
        # Input identity threshold
        self.identity_slider = QSlider(Qt.Horizontal)
        self.identity_slider.setMinimum(MIN_THRESHOLD)
        self.identity_slider.setMaximum(MAX_THRESHOLD)
        self.identity_slider.setTickInterval(1)
        self.identity_slider.setTickPosition(QSlider.TicksBelow)
        # Blast button
        self.blast_button = QPushButton("🚀 BLAST")
        self.tree_button = QPushButton("🌳 Generate Tree")
        # Labels
        self.blast_mode_label = QLabel("BLAST mode")
        self.rank_label = QLabel("Taxonomy rank")
        self.identity_label = QLabel("Identity threshold")
        # Create layou to add labels for numbers below the slider
        self.label_row = QHBoxLayout()
        for i in range(4):  # inclusive
            lbl = QLabel(str(MIN_THRESHOLD + i))
            if i == 0:
                lbl.setAlignment(Qt.AlignLeft)
            elif i == 3:
                lbl.setAlignment(Qt.AlignRight)
            else:
                lbl.setAlignment(Qt.AlignCenter)
            self.label_row.addWidget(lbl)
        """ Widget's layout """
        # Define layour
        self.layout = QVBoxLayout()
        # Add blast mode input
        self.layout.addSpacing(30)
        self.layout.addWidget(self.blast_mode_label)
        self.layout.addWidget(self.blast_mode_dropdown)
        self.layout.addSpacing(30) 
        # Add taxonomy rank input
        self.layout.addWidget(self.rank_label)
        self.layout.addWidget(self.rank_dropdown)
        self.layout.addSpacing(30) 
        # Add identity threshold input
        self.layout.addWidget(self.identity_label)
        self.layout.addWidget(self.identity_slider)
        self.layout.addLayout(self.label_row)
        # Add space between buttons and other widgets
        self.layout.addStretch()
        # Add search button
        self.layout.addWidget(self.blast_button)
        self.layout.addSpacing(30)     
        self.layout.addWidget(self.tree_button)
        self.layout.addSpacing(30)     
        # Set layout
        self.setLayout(self.layout)


class LogBook(QFrame):
    """Logbook to write all the events that have occurred in the app.

    This class inherits from QFrame. It incorporates a text box where every time
    and action is carried out, it is recorded with the time it happened.

    Attributes (widgets):
        status_box(QTextEdit): box where the actions are written.
    """
    def __init__(self):
        """Creates an instance of the class."""
        super().__init__()
        """Add widgets."""
        # Status box
        self.status_box = QTextEdit()
        self.status_box.setReadOnly(True)
        self.status_box.setObjectName("logBook")        # To apply custom style
        # Labels
        self.log_label = QLabel("📓 LOGBOOK")
        self.log_label.setObjectName("logBookTitle")    # To apply custom style
        """Widgets' layout."""
        # Define layour
        self.layout = QVBoxLayout()
        # Add status box
        self.layout.addWidget(self.log_label)
        self.layout.addWidget(self.status_box)
        # Set layout
        self.setLayout(self.layout)
    
    def log(self, message):
        """Writes in status_bow the action with the time when the action happened."""
        timestamp = datetime.now().strftime("[%H:%M:%S]")
        self.status_box.append(f"{timestamp} {message}")


class InputModule(QTabWidget):
    """Tab widget with that incorporates a SearchTab and a BlastTab.

    Attributes (widgets):
        search_tab(SearchTab): see SearchTab class
        blast_tab(BlastTab): see BlasTab class
    """
    def __init__(self, controller):
        """ Initialises the instance 

        Args:
            controller (MainController): controller to take the actions of the signals
        """
        super().__init__()
        # Add controller
        self.controller = controller
        # Set position tabs to horizontal on top
        self.setTabPosition(QTabWidget.North)
        """ Add tabs """
        self.search_tab = SearchTab(self.controller)
        self.blast_tab = BlastTab()
        # Add tabs to tab list widget
        self.insertTab(0, self.search_tab, "SEARCH")
        self.insertTab(1, self.blast_tab, "BLAST")


class Sequence(QTextEdit):
    """Text widget to display the barcode sequence and highlight it.

    Attributes (widgets):
        sequence(str): nucleotide string
        primers_intervals(list of tuples): intervals marking the positions where the
            nucleotide sequence aligns with the primers.
    """
    def __init__(self, sequence, primers_intervals):
        """ Initialises an instance of the class.

        Args:
            sequence(str): nucleotide string
            primers_intervals(list of tuples): intervals marking the positions where the
                nucleotide sequence aligns with the primers.
        """
        super().__init__()
        # Add sequence and intervals
        self.sequence = sequence
        self.primers_intervals = primers_intervals
        # Set the mouse tracking for hover effects
        self.setMouseTracking(True)
        # Display the sequence
        self.setPlainText(self.sequence)
        
    def enterEvent(self, event):
        """Enters when the mouse is on top of the sequence. Calls function to highlight."""
        self.apply_highlight(self.primers_intervals, (QColor("#FAC898"), QColor("#B3EBF2")))
        super().enterEvent(event)
        
    def leaveEvent(self, event):
        """Enters when the mouse is NOT on top of the sequence. Calls function to delete highlight."""
        self.apply_highlight(self.primers_intervals, (Qt.transparent, Qt.transparent))
        super().leaveEvent(event)

    def apply_highlight(self, intervals, colours):
        """Applies highlighter to a set of intervals with a given color."""
        for i, primer_list in enumerate(intervals):
            cursor = self.textCursor()
            fmt = QTextCharFormat()
            fmt.setBackground(colours[i])

            for start, end in primer_list:
                if start > end:
                    continue  # skip invalid intervals
                length = end - start + 1
                cursor.setPosition(start)
                cursor.movePosition(QTextCursor.Right, QTextCursor.KeepAnchor, length)
                cursor.mergeCharFormat(fmt)


class AlignmentPopup(QWidget):
    """Pop up window to display the alignment of two sequences.

    Attributes (widgets):
        alignment_text(QTextEdit): text box to display the alignment
    """
    def __init__(self, alignment, seq1_id, seq2_id):
        """ Initialises an instance of the class.

        Args:
            alignment(Alignment): alignment object between two sequences
            seq1_id(str): string of nucleotides of the reference sequence (the checked one)
            seq2_id(str): string of nucleotides of the selected sequence (right-clisked one)
        """
        super().__init__()
        self.setObjectName("alignmentPopup") # To apply custom style
        # Set name of the window
        self.setWindowTitle(f'Alignment between {seq1_id} and {seq2_id}')
        """Add widgets."""
        self.alignment_text = QTextEdit()
        self.alignment_text.setText(str(alignment))
        self.alignment_text.setObjectName("alignmentTextPopup")
        self.alignment_text.setReadOnly(True)
        """Define layout"""
        layout = QVBoxLayout()
        layout.addWidget(self.alignment_text)
        self.setLayout(layout)
        # Set size of the window
        self.resize(800, 600)


class BarcodeCard(QFrame):
    """Card to display the searched barcode sequence.

    It inherits from the QFrame class. It is a selectable card with a radio button
    that displays a sequence of a specific barcode and organism retrieved from 
    a fasta file uploaded to NCBI. It also shows the number of duplicates, i.e. 
    the number of fasta files uploaded to NCBI where the sequence is the same.

    Attributes (widgets):
        sequence(str): nucleotide sequence
        select_radio(QRadioButton): radio button to select a sequence
        sequence_text(Sequence): see Sequence class
        id(str): fasta file id for the barcode. Serves as id for the card to differentiate
            it from other cards in the same tab.
    """
    def __init__(self, parent, sequence, header, duplicates, primers_intervals):
        """ Initialises an instance of the class.

        Args:
            parent(BarcodesTab): see BarcodesTab class
            sequence(str): nucleotide sequence
            header(str): first header of the barcode sequence
            duplicates(int): number of duplicates
            primers_intervals(list of tuples): intervals marking the positions where the
                nucleotide sequence aligns with the primers.
        """
        super().__init__()
        # Add parent (the tab widget)
        self.parent = parent
        # Add sequence 
        self.sequence = sequence
        self.setObjectName("barcodeCard")   # To apply custom style
        """Add widgets."""
        # Fasta header label
        self.header_label = QLabel(header)
        self.header_label.setObjectName("barcodeLabel")
        # Set the identifier of the card as the fasta file id
        self.id = header.split(" ")[0][1:]
        # Number of duplicates label
        self.duplicates_label = QLabel(f'Duplicates: {duplicates}')
        self.duplicates_label.setObjectName("barcodeLabel")
        # Radio button to select card
        self.select_radio = QRadioButton()
        self.select_radio.clicked.connect(self.check)
        # Sequence text class
        self.sequence_text = Sequence(sequence, primers_intervals)
        self.sequence_text.setReadOnly(True)
        # Adjust size of widget
        self.sequence_text.document().adjustSize()
        height = self.sequence_text.document().size().height()
        self.sequence_text.setFixedHeight(int(height/3)+10)
        self.sequence_text.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.sequence_text.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        """Widget's layout."""
        # Define layout
        self.layout = QVBoxLayout()
        self.top_layout = QHBoxLayout()
        # Add fasta header to top_layout
        self.top_layout.addWidget(self.header_label)
        # Add stretch between fasta header and radio button
        self.top_layout.addStretch()
        # Add radio button
        self.top_layout.addWidget(self.select_radio)
        # Add top layout to main layout
        self.layout.addLayout(self.top_layout)
        # Add duplicates label
        self.layout.addWidget(self.duplicates_label)
        # Add sequence text
        self.layout.addWidget(self.sequence_text) 
        # Set layout
        self.setLayout(self.layout)
    
    def check(self): 
        """Checks if the card is selected and acts accordingly."""
        if self.select_radio.isChecked():
            # The barcode was previously unchecked so this is the new checked card
            self.parent.uncheck_barcode_cards(self.id)  # Uncheck other cards     
            self.parent.check_barcode_card(self.id)     # Select this one as checked
        else:
            # The barcode was previously checked
            self.parent.uncheck_barcode_cards(None) # Uncheck all cards
                    
    def uncheck(self):
        """Unchecks the card."""
        self.select_radio.setChecked(False)
    
    def mousePressEvent(self, event):
        "Function to crete pop up window with the alignment."
        # Only works if a card is selected and the user right clicks a card
        if event.button() == Qt.RightButton and self.parent.checked_card:
            # The card can be the same as the selected one
            if self.parent.checked_card.id != self.id:
                # Get sequence from reference card
                ref_seq = self.parent.checked_card.sequence
                # Call controller to obtain alignment
                alignment = self.parent.controller.align_seqs(ref_seq, self.sequence)
                # Create popup
                popup = AlignmentPopup(alignment, self.parent.checked_card.id, self.id)
                popup.show()
                # Keep a reference to avoid garbage collection
                self._popup = popup
    

class ResultsTab(QScrollArea):
    """Tab to store the results obtained after the BLAST.

    Attributes (widgets):
        
    """
    def __init__(self):
        """Creates an instance of the class."""
        super().__init__()
        # Allow to be resized
        self.setWidgetResizable(True)
        """Add widgets."""
        self.result_container = QWidget()
        self.result_container.setObjectName("outputWidget")
        """Widgets' layout"""
        # Define layout
        self.result_layout = QVBoxLayout()
        # Set layout
        self.result_container.setLayout(self.result_layout)
        self.setWidget(self.result_container)


class BarcodesTab(QScrollArea):
    """Tab to display the barcode cards.

    Attributes (widgets):
        barcode_card_list(list): list of BarcodeCard instances contained in the tab.
        checked_card(str): id of the selected card
        barcodes_data(list): list of lists containing the data for each barcode card. The
            list includes: sequence, fasta header, number of duplicates and primer intervals.
        name(str): name of the tab which corresponds to the type of barcode
    """
    def __init__(self, controller, tab_name, barcodes_data):
        """ Initialises an instance of the class.

        Args:
            controller(MainController): controller to call in functions.
            barcodes_data(list): list of lists containing the data for each barcode card. The
                list includes: sequence, fasta header, number of duplicates and primer intervals.
        """
        super().__init__()
        # Checked card
        self.checked_card = None
        # Add controller
        self.controller = controller
        # Add data
        self.barcodes_data = barcodes_data
        # Add name of the tab (barcode type)
        self.tab_name = tab_name
        # Allow to be resized
        self.setWidgetResizable(True)
        # List of barcode cards
        self.barcode_card_list = []
        # Add barcode data to barcode list
        self.populate_barcode_tab()
        # Allocate widgets in the tab
        self.result_container = QWidget()
        self.result_container.setObjectName("outputWidget")
        self.result_layout = QVBoxLayout()
        for barcode_card in self.barcode_card_list:
            barcode_card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
            self.result_layout.addWidget(barcode_card)
        # Add expanding spacer to push widgets up and prevent stretching
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
        self.result_layout.addWidget(spacer)       
        self.result_container.setLayout(self.result_layout)
        self.setWidget(self.result_container)
    
    def populate_barcode_tab(self):
        """Populates the tab with barcode cards"""
        # Extract data for each element in list
        for data_point in self.barcodes_data:
            sequence = data_point[0]
            header = data_point[1]
            duplicates = data_point[2]
            primers_intervals = data_point[3]
            # Create barcode card
            new_barcode_card = BarcodeCard(self, sequence, header, duplicates, primers_intervals)
            # Add barcode card to list
            self.barcode_card_list.append(new_barcode_card)
    
    def check_barcode_card(self, barcode_id):
        """Adds a barcode card that has been checked to the variable checked_card given its id."""
        for barcode_card in self.barcode_card_list:
            if barcode_card.id == barcode_id:
                self.checked_card = barcode_card
                self.controller.add_sequence_to_blast(barcode_card.sequence, self.tab_name)
        self.controller.write_in_logbook(barcode_id)

    def uncheck_barcode_cards(self, barcode_id):
        """Deletes a barcode card that has been unchecked from the variable checked_card given its id."""
        for barcode_card in self.barcode_card_list:
            if barcode_card.id != barcode_id:
                barcode_card.uncheck()
        self.checked_card = None

    
class OutputModule(QTabWidget):
    """Module for the output data: barcode cards and results.

    Attributes (widgets):
        results_tab(ResultsTab): see ResultsTab
    """
    def __init__(self, controller):
        """ Initialises an instance of the class.

        Args:
            controller(MainController): see BarcodesTab class
        """
        super().__init__()
        # Add controller
        self.controller = controller
        # Set position tabs to horizontal on top
        self.setTabPosition(QTabWidget.North)
        """Add tabs."""
        self.results_tab = ResultsTab()
        # Add tabs to tab widget
        self.insertTab(0, self.results_tab, "Results")
        self.setCurrentIndex(0)
    
    def add_new_barcode_tab(self, barcode_name, barcodes_data):
        """Adds a barcode tab given a barcod name and a list of data for the barcode cards"""
        # Create tab
        self.new_barcode_tab = BarcodesTab(self.controller, barcode_name, barcodes_data)
        # Check if a tab already exists for this barcode
        for index in range(self.count()):
            if self.tabText(index) == barcode_name:
                # Replace the widget in the tab
                self.removeTab(index)
                self.insertTab(index, self.new_barcode_tab, barcode_name)
                self.setCurrentIndex(index)  # Focus the updated tab
                return
        # If no match found, add new tab
        self.addTab(self.new_barcode_tab, barcode_name)
        self.setCurrentIndex(self.count() - 1)  
          