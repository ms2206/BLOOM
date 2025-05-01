"""
Contains the declaration of all custom widgets as classes for the GUI.
"""

from datetime import datetime
from PyQt5.QtCore import Qt, QStringListModel, QTimer, QRectF
from PyQt5.QtGui import QTextCursor, QTextCharFormat, QColor, QPainter, QBrush, QFont
from PyQt5.QtWidgets import *
from config.config_manager import ConfigManager


"""Configuration constants."""
config = ConfigManager()
# Barcode info
BARCODE_LIST = config.get_barcodes_names()
PRIMERS_DICT = config.get_barcodes_primers()
# BLAST configuration
BLAST_MODES = config.get('blast_modes')
TAXONOMY_RANKS = config.get('taxonomy_ranks')
# Results visualisation configuration
DATA_TYPES = config.get('data_types')
ANIMATION_DURATION = config.get('animation_duration')
# Github link for help
GITHUB_LINK = config.get('github_link')


class PieChartWidget(QWidget):
    """Widget used to display a donut chart with the stats for the results obtained.

    Displays a donut chart with sectors representing how many hits/species have a specificic number
    of differences or fall in a range of identity percentages. In the middle there is a text label with 
    the total number of hits/species recorded. Each sector's size is dependent on the value linked to it.
    Each sector has a colour that ranges from red (for similar hits) to green (dissimilar hits).
    The chart is animated with a pan movement so that every time the results tab the chart is drawn again for
    aesthetic purposes.

    Arguments:
        values (list): list of numerical values to display. Each sector in the chart has a value
        colours (list): list of colours for each value. Each sector in the chart has a colour
        animation_progress (double): used to track how much has been animated
        start_angle_offset (int): starting position of the animation.
        total_angle (int): total angle that is goint go be covered by the chart
    """
    def __init__(self, values, colours):
        """Creates an instance of the class.

        Args:
            values (list): list of numerical values to display. Each sector in the chart has a value
            colours (list): list of colours for each value. Each sector in the chart has a colour
        """
        super().__init__()
        # Set data and colours
        self.values = values
        self.colours = colours
        # Animation parameters
        self.animation_progress = 0.0
        self.start_angle_offset = -90 # Start drawing from the top
        self.total_angle = 360

    def set_progress(self, progress):
        """Sets the current progress of the animation."""
        self.animation_progress = progress
        self.update()

    def paintEvent(self, event):
        """Paints the chart for the current progress."""
        # Define painter object
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        # Define size of drawing
        side = min(self.width(), self.height()) - 20
        rect = QRectF((self.width() - side) / 2, (self.height() - side) / 2, side, side)
        # Calculate current angle 
        progressAngle = self.animation_progress * self.total_angle
        startAngle = self.start_angle_offset
        total_percent = sum(self.values)
        # Draw the pie chart
        for i, value in enumerate(self.values):
            # Scale each slice to fit 300° instead of 360°
            sliceAngle = (value / total_percent) * self.total_angle
            # Check if drawing is complete
            if startAngle - self.start_angle_offset >= progressAngle:
                break
            # Check which slice is set to draw
            if (startAngle - self.start_angle_offset) + sliceAngle <= progressAngle:
                drawAngle = sliceAngle
            else:
                drawAngle = progressAngle - (startAngle - self.start_angle_offset)
            # Paint with the coulours given
            painter.setBrush(QBrush(self.colours[i]))
            painter.setPen(Qt.NoPen)
            painter.drawPie(rect, int(-startAngle * 16), int(-drawAngle * 16))
            # Define new start angle for next instant drawing
            startAngle += sliceAngle
        # Draw the center hole
        innerRatio = 0.5
        innerRectSize = side * innerRatio
        innerRect = QRectF(
            (self.width() - innerRectSize) / 2,
            (self.height() - innerRectSize) / 2,
            innerRectSize,
            innerRectSize,
        )
        painter.setBrush(QBrush(QColor("#ffffff")))
        painter.drawEllipse(innerRect)
        # Animated number in the center
        target_number = sum(self.values)
        current_number = int(target_number * self.animation_progress)
        painter.setPen(QColor(0, 0, 0))
        font = QFont("Bahnschrift Semibold", int(innerRectSize * 0.2))
        painter.setFont(font)
        painter.drawText(innerRect, Qt.AlignCenter, str(current_number))


class BarChartWidget(QWidget):
    """Widget used to display a bar chart with the stats for the results obtained.

    Displays a histogram with bars representing how many hits/species have a specificic number
    of differences or fall in a range of identity percentages. The X axis has a label to interpret
    what bars represent. Each bar's height is dependent on the value linked to it.
    Each bar has a colour that ranges from red (for similar hits) to green (dissimilar hits).
    The chart is animated by increasing gradually the bars' heights so that every time the results 
    tab the chart is drawn again for aesthetic purposes.

    Arguments:
        labels (list): list of labels for the x-axis.
        values (list): list of numerical values to display. Each sector in the chart has a value.
        colours (list): list of colours for each value. Each sector in the chart has a colour.
        title (str): title for the x-axis
        animation_progress (double): used to track how much has been animated.
    """
    def __init__(self, values, labels, colours, data_type):
        """Creates an instance of the class.

        Args:
            values (list): list of numerical values to display
            labels (list): list of labels for the X-Axis in the bar chart
            colours (list): list of colours for each value. Each sector in the chart has a colour.
            data_type (str): type of data to be represented (percentages of number of differences)
        """
        super().__init__()
        # Define data numbers, colours
        self.labels = labels
        self.values = values
        self.colours = colours
        # Define the title given the type of data to represent
        self.title = "Number of differences" if data_type == "differences" else "Percentage of identity"
        # Animation parameters
        self.animation_progress = 0.0

    def set_progress(self, progress):
        """Sets the current progress of the animation."""
        self.animation_progress = progress
        self.update()

    def paintEvent(self, event):
        """Paints the chart for the current progress."""
        # Define painter
        painter = QPainter(self)
        painter.setRenderHint(QPainter.Antialiasing)
        # Set drawing parameters
        width = self.width()
        height = self.height()
        num_bars = len(self.values)
        spacing = 30
        bar_width = (width - (num_bars + 1) * spacing) / num_bars
        max_value = max(self.values)
        bottom_margin = 50
        top_margin = 30
        bar_area_height = height - top_margin - bottom_margin
        # Set font for labels
        font = QFont("Bahnschrift Semibold", 10)
        painter.setFont(font)
        # Draw bars
        for i, value in enumerate(self.values):
            bar_height = (value / max_value) * bar_area_height * self.animation_progress # Margin for labels
            x = spacing + i * (bar_width + spacing)
            y = top_margin + bar_area_height - bar_height
            # Add colours
            painter.setBrush(QBrush(self.colours[i]))
            painter.setPen(Qt.NoPen)
            painter.drawRect(QRectF(x, y, bar_width, bar_height))
            # Value label above
            painter.setPen(QColor(0, 0, 0))
            painter.drawText(QRectF(x, y - 20, bar_width, 20), Qt.AlignCenter, str(int(value)))
            # Index label
            label = self.labels[i]
            painter.drawText(QRectF(x, top_margin + bar_area_height + 5, bar_width, 20), Qt.AlignCenter, label)
        # Axis label
        painter.setFont(QFont("Bahnschrift Semibold", 12, QFont.Bold))
        painter.drawText(QRectF(0, height - 20, width, 20), Qt.AlignCenter, self.title)


class CombinedChartWindow(QWidget):
    """Widget that contains one line of results to be displayed in the results tab.

    This widget contains a title that indicates if the results are for all hits or just unique
    species, the donut plot and the bar plot. The parameters needed to set up the animation 
    process are defined in this class.

    Arguments: 
        data (list): numerical data to be displayed in the graphs
        labels (list): list of labels for the X-Axis in the barc chart
        title (str): title of the results
        timer (QTimer): timer object used for the animations
        elapsed (double): elapsed time that has passed since the animation started
        animation_interval (int): time interval used to calculate the time that has to
            pass for the animation to progress
    
    Widgets:
        title_label (QLabel): label with the title of the results
        pie_chart (PieChartWidget): donut chart
        bar_chart (BarChartWidget): bar plot
    """
    def __init__(self, data, labels, title, data_type):
        """Creates an instance of the class.

        Args:
            data (list): numerical data to be displayed in the graphs
            labels (list): list of labels for the X-Axis in the bar chart
            title (str): title of the results
            data_type (str): type of data to be represented (percentages of number of differences)
        """
        super().__init__()
        self.data = data
        self.labels = labels
        self.title = title
        """Widgets."""
        # Label with the title
        self.title_label = QLabel(self.title)
        self.title_label.setObjectName("resultsLabel")
        # Colour palette
        palette = self.get_color_palette(len(self.data))
        # Pie chart
        self.pie_chart = PieChartWidget(values=self.data, colours=palette)
        # Bar chart
        self.bar_chart = BarChartWidget(values=self.data, labels=self.labels, colours=palette, data_type=data_type)
        # Splitter
        main_splitter = QSplitter(Qt.Horizontal)
        graphs_splitter = QSplitter(Qt.Horizontal)
        """Splitters layout"""
        # Splitter for graphs
        graphs_splitter.addWidget(self.pie_chart)
        graphs_splitter.addWidget(self.bar_chart)
        graphs_splitter.setSizes([500, 500])
        graphs_splitter.setStretchFactor(1, 1)
        # Main splitter
        main_splitter.addWidget(self.title_label)
        main_splitter.addWidget(graphs_splitter)
        main_splitter.setSizes([100, 900])
        main_splitter.setStretchFactor(1, 1)
        """Widgets' layout"""
        # Layout
        layout = QHBoxLayout()
        layout.addWidget(main_splitter)
        self.setLayout(layout)
        """Animation parameters-"""
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_animation)
        self.elapsed = 0
        self.animation_interval = 16  # ~60 FPS
        self.restart_animation()

    def restart_animation(self):
        """Start the animation again from the beginning."""
        self.elapsed = 0
        self.timer.start(self.animation_interval)

    def update_animation(self):
        """Updates the animation given the time elapsed."""
        self.elapsed += self.animation_interval
        progress = min(self.elapsed / ANIMATION_DURATION, 1.0)
        self.pie_chart.set_progress(progress)
        self.bar_chart.set_progress(progress)
        # Stop when progress is completed
        if progress >= 1.0:
            self.timer.stop()
    
    def get_color_palette(self, n):
        """
        Creates a list of colours given the number of values to represent. The colour palette
        goes from red to green passing trhough yellow.
        """
        # Return a single colour in case there is only one value
        if n == 1:
            return [QColor(*(255, 0, 97))]
        # Calculate step between colours
        step = int(510/(n-1))
        f_palette = []  # First half of the palette
        r_palette = []  # Second half of the palette
        # Get colour progression 
        for i in range(int(n/2)):
            f_colour = (255, 0+i*step, 97) # Starts in red and goes to yellow
            r_colour = (0+i*step, 255, 97) # Starts in green and goes to yellow
            f_palette.append(QColor(*f_colour))
            r_palette.append(QColor(*r_colour))
        # Compensate with a yellow colour if the number of values is odd
        if n%2 == 1:
            colour = (255, 255, 97)
            f_palette.append(QColor(*colour))    
        # Append reversed second half of the palette
        r_palette.reverse()
        palette = f_palette + r_palette
        return palette
        

class CustomMenuBar(QMenuBar):
    """Menu bar at the top of the app

    Attributes:
        controller (MainController): controller to perform actions when options are pressed

    Wdgets:
        file_menu (QMenu): contains the actions to close the app and generate report 
        help_menu (QMenu): contains the action to provide the user with help 
    """
    def __init__(self, parent, controller):
        super().__init__(parent)
        """Creates an instance of the class. 

        Args:
            parent (QMainWindow): main window of the app where the menu belongs
            controller (MainController): controller to perform actions when options are pressed
        """
        self.controller = controller
        # File Menu
        self.file_menu = self.addMenu("File")
        exit_action = QAction("Exit", self)             # Close the app
        exit_action.triggered.connect(parent.close)
        self.file_menu.addAction(exit_action)
        csv_action = QAction("Create .csv file", self)
        csv_action.triggered.connect(self.write_csv)
        self.file_menu.addAction(csv_action)
        clear_action = QAction("Clear workspace", self) # Clear the workspace
        clear_action.triggered.connect(self.clear_workspace)
        self.file_menu.addAction(clear_action)
        # Help Menu
        self.help_menu = self.addMenu("Help")
        help_action = QAction('Help', self) # Show a pop up message with the link to the github
        help_action.triggered.connect(self.show_help_message)
        self.help_menu.addAction(help_action)
    
    def write_csv(self):
        """
        Creates a pop up window to select name and location of csv file and passes the
        info to the controller to create and save the csv file.
        """
        # Open the Save File dialog
        options = QFileDialog.Options()
        file_path, _ = QFileDialog.getSaveFileName(None,
                                                   "Save CSV File",
                                                   "",
                                                   "CSV Files (*.csv);;All Files (*)",
                                                   options=options)
        # Exit if the user cancels
        if not file_path:
            print("Save cancelled.")
            return
        # Ensure file ends with .csv
        if not file_path.endswith('.csv'):
            file_path += '.csv'
        self.controller.write_csv(file_path)
    
    def clear_workspace(self):
        """Commands the controller to clean the entire workspace."""
        self.controller.clear()
    
    def show_help_message(self):
        """Shows a pop-up window with a message with the link to the github."""
        message = f'Check this link for more info:\n{GITHUB_LINK}'
        QMessageBox.information(self, "Help", message)


class Header(QWidget):
    """Header with the title of the app and button to hide the input module.

    Attributes:
        main_window (QMainWindow): main window of the app where the title belongs

    Widgets:
        toggle_button (QToolButton): hides and shows the input module
        title (QLabel): full title of the app
    """
    def __init__(self, main_window):
        """Creates an instance of the class.

        Args:
            main_window (QMainWindow): main window of the app where the title belongs
        """
        super().__init__()
        # Add main window
        self.main_window = main_window
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
        layout = QHBoxLayout()
        layout.setSpacing(10)
        # Add toggle button
        layout.addWidget(self.toggle_button)
        # Add title
        layout.addWidget(self.title)
        # Set layout
        self.setLayout(layout)

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

    This incorporates all widgets to input the paramaeters for the barcode search. This includes 
    the organism to study, the type of barcode and the sequence for the forward and reverse primers.
    The tab has an autocompleter for the organism name so that the user can get help to write it. The 
    barcode selection has a drop down to select from a list of barcode types. Finally the primer sequences
    are autocompleted depending on the barcode selected but they can be later edited by the user.

    Attributes:
        controller (MainController): controller to perform actions when buttons are pressed
        full_species_list (list): list of all species in NCBI for the autocomplete
    
    Widgets:
        species_input(QLineEdit): line widget to enter the name of the organism
        species_model(QStringListModel): list of organism names to choose from
        completer(QCompleter): completer for the organism name
        barcode_dropdown(QComboBox): to choose a barcode
        f_primer_input(QTextEdit): to enter the forward primer sequence
        r_primer_input(QTextEdit): to enter the reverse primer sequence
        search_button(QPushButton): button to search the selected parameters
    """
    def __init__(self, controller):
        """Creates an instance of the class.

        Args:
            controller (MainController): controller to perform actions when buttons are pressed
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
        completer = QCompleter(self.species_model, self)
        completer.setCaseSensitivity(Qt.CaseInsensitive)
        completer.setFilterMode(Qt.MatchContains)
        completer.setCompletionMode(QCompleter.PopupCompletion)
        # Popup for the list of species
        popup = QListView()
        completer.setPopup(popup)
        # Add completer and popup to species_input and add signal to suggest names
        self.species_input.setCompleter(completer)
        self.species_input.textEdited.connect(self.update_species_completer)
        # Barcode selection dropdown
        self.barcode_dropdown = QComboBox()
        self.barcode_dropdown.addItem("Choose barcode...")      # Add an empty option
        self.barcode_dropdown.model().item(0).setEnabled(False)
        self.barcode_dropdown.model().item(0).setForeground(QColor('gray'))
        self.barcode_dropdown.addItems(BARCODE_LIST)
        self.barcode_dropdown.activated.connect(self.autocomplete_primers)
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
        self.search_button.setObjectName("regularButton")
        self.search_button.clicked.connect(self.search_button_action)
        # Labels
        species_label = QLabel("Species")
        species_label.setObjectName("regularLabel")
        barcode_label = QLabel("Barcode")
        barcode_label.setObjectName("regularLabel")
        f_primer_label = QLabel("Forward primer")
        f_primer_label.setObjectName("regularLabel")
        r_primer_label = QLabel("Reverse primer")
        r_primer_label.setObjectName("regularLabel")
        """ Widgets' layout """
        # Define layour
        layout = QVBoxLayout()
        # Add species input
        layout.addSpacing(30) 
        layout.addWidget(species_label)
        layout.addWidget(self.species_input)
        layout.addSpacing(30) 
        # Add barcode input
        layout.addWidget(barcode_label)
        layout.addWidget(self.barcode_dropdown)
        layout.addSpacing(30) 
        # Add primers input
        layout.addWidget(f_primer_label)
        layout.addWidget(self.f_primer_input)
        layout.addSpacing(10) 
        layout.addWidget(r_primer_label)
        layout.addWidget(self.r_primer_input)
        # Add space between button and other widgets
        layout.addStretch()
        # Add search button
        layout.addWidget(self.search_button)
        layout.addSpacing(30)     
        # Set layout
        self.setLayout(layout)
    
    def search_button_action(self):
        """Gathers the data inputted by the user and commands the controller to search."""
        # Get organism and primers chosen by the user
        organism = self.species_input.text()
        index = self.barcode_dropdown.currentIndex()
        f_primer = self.f_primer_input.toPlainText()
        r_primer = self.r_primer_input.toPlainText()
        primers = (f_primer, r_primer)
        if not organism or not f_primer or not r_primer or index == 0:
            message = "One or more fields are incomplete."
            self.controller.error_pop_up(message)
            return
        # Get barcode type taking into account the blank option
        barcode_type = BARCODE_LIST[index-1]
        # Logbook
        self.controller.write_in_logbook(f'Fetching {barcode_type} sequences for {organism}.')
        # Call controller to get barcodes
        searching = self.controller.search_barcodes(organism, barcode_type, primers)
        if searching:
            self.controller.write_in_logbook(f'Error when searching for barcodes')
            return
        else:
            # Block organism name button so user can't change it
            self.species_input.setEnabled(False)
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
    
    def autocomplete_primers(self):
        """Autocompletes the primers boxes with the primers of the selected barcode."""
        # Get index
        index = self.barcode_dropdown.currentIndex()
        # If index is first option return (blank option)
        if index == 0:
            return
        # Get barcode type from list to get primers
        barcode_type = BARCODE_LIST[index-1]
        primers = PRIMERS_DICT[barcode_type]
        self.f_primer_input.setText(primers[0])
        self.r_primer_input.setText(primers[1])
    
    def clear(self):
        """Deletes any info written by the user in the module and enables the organism name to be changed."""
        self.species_input.setEnabled(True)
        self.species_input.setText("")
        self.f_primer_input.setText("")
        self.r_primer_input.setText("")


class BlastTab(QFrame):
    """Tab in the input module to define the parameters for BLAST and to generate the tree.

    Thi tab incorporates the entry widgets to analyse results. It is divided intro three sections:
        - BLAST section: to select parameters to analyse barcodes with the BLAST tool and the button to 
        start it
        - DATA DISPLAY section: to select how to show the results obtained for better interpretation
        - TREE: button to create and show the taxonomy tree from the results obtained

    Attributes:
        controller (MainController): controller to perform actions when buttons are pressed
      
    Widgets:
        blast_mode_dropdown(QComboBox): to select the blast mode
        rank_dropdown(QComboBox): to select the taxonomy rank to analyse
        data_type_dropdown (QComboBox): to select wchich of type of data to show the results
        dissimilars_checkbox (QCheckBox): allows the user to display or not dissimilar species as they 
            usually minimise other stats.
        update_button (QPushButton): to redraw the stats with the selected parameters
        blast_button(QPushButton): to use BLAST with the parameters selected
        tree_button(QPushButton): to generate and display the taxonomy tree
    """
    def __init__(self, controller):
        """Initialises an instance of the class.
            
        Args:
            controller (MainController): controller to perform actions when buttons are pressed
        """
        super().__init__()
        # Add controller
        self.controller = controller
        """ Add Widgets """
        # Input blast mode
        self.blast_mode_dropdown = QComboBox()
        self.blast_mode_dropdown.addItems(BLAST_MODES)
        # Input taxonomy rank
        self.rank_dropdown = QComboBox()
        self.rank_dropdown.addItems(TAXONOMY_RANKS)
        # Select the type of data to display
        self.data_type_dropdown = QComboBox()
        self.data_type_dropdown.addItems(DATA_TYPES)
        # Select which data to represent
        self.dissimilars_checkbox = QCheckBox('Show dissimilar hits', self)
        # Update results button
        self.update_button = QPushButton("Update results")
        self.update_button.setObjectName("regularButton")
        self.update_button.clicked.connect(self.update_button_action)
        # Blast button
        self.blast_button = QPushButton("🚀 BLAST")
        self.blast_button.setObjectName("regularButton")
        self.blast_button.clicked.connect(self.blast_button_action)
        # Tree button
        self.tree_button = QPushButton("🌳 Generate Tree")
        self.tree_button.setObjectName("regularButton")
        self.tree_button.clicked.connect(self.show_tree)
        # Labels
        blast_mode_label = QLabel("BLAST mode")
        blast_mode_label.setObjectName("regularLabel")
        rank_label = QLabel("Taxonomy rank")
        rank_label.setObjectName("regularLabel")
        data_type_label = QLabel("Type of data")
        data_type_label.setObjectName("regularLabel") 
        """ Widget's layout """
        # Define layour
        layout = QVBoxLayout()
        # Add blast mode input
        layout.addSpacing(30)
        layout.addWidget(blast_mode_label)
        layout.addWidget(self.blast_mode_dropdown)
        layout.addSpacing(20) 
        # Add taxonomy rank input
        layout.addWidget(rank_label)
        layout.addWidget(self.rank_dropdown)
        layout.addSpacing(20) 
        # Add blast button
        layout.addWidget(self.blast_button)
        layout.addSpacing(50)  
        # Add a line to divide sections
        line1 = QFrame()
        line1.setFrameShape(QFrame.HLine)
        line1.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line1)
        layout.addSpacing(50)  
        # Add section to update results
        layout.addWidget(data_type_label)
        layout.addWidget(self.data_type_dropdown)
        layout.addSpacing(20)
        layout.addWidget(self.dissimilars_checkbox)
        layout.addSpacing(20)
        layout.addWidget(self.update_button)
        layout.addSpacing(50)
        # Add a line to divide sections
        line2 = QFrame()
        line2.setFrameShape(QFrame.HLine)
        line2.setFrameShadow(QFrame.Sunken)
        layout.addWidget(line2)
        layout.addSpacing(40)  
        # Add tree button
        layout.addWidget(self.tree_button)
        # Add space after all widgets
        layout.addStretch()   
        # Set layout
        self.setLayout(layout)
    
    def blast_button_action(self):
        """Calls the main controller to do a BLAST search."""
        blast_mode = BLAST_MODES[self.blast_mode_dropdown.currentIndex()]
        taxonomy_rank = TAXONOMY_RANKS[self.rank_dropdown.currentIndex()]
        blasting = self.controller.start_blast(blast_mode, taxonomy_rank)
        if blasting:
            self.controller.write_in_logbook("Error when using BLAST")
    
    def update_button_action(self):
        """Updates the results display with the new parameters selected."""
        data_type = DATA_TYPES[self.data_type_dropdown.currentIndex()]
        show_dissimilars = self.dissimilars_checkbox.isChecked()
        self.controller.update_results(data_type, show_dissimilars)
    
    def show_tree(self):
        """Calls controller to create the taxonomy tree and display it"""
        self.controller.create_tree()


class LogBook(QFrame):
    """Logbook to write all the events that have occurred in the app.

    This module contains a text box where every time and action is carried out, it is recorded 
    with the time it happened. This aids the user to check what actions have been previously performed
    to keep track of them and the time an action has taked to be carried out.

    Widgets:
        status_box (QTextEdit): box where the actions are written.
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
        log_label = QLabel("📓 LOGBOOK")
        log_label.setObjectName("logBookTitle")    # To apply custom style
        """Widgets' layout."""
        # Define layout
        layout = QVBoxLayout()
        # Add status box
        layout.addWidget(log_label)
        layout.addWidget(self.status_box)
        # Set layout
        self.setLayout(layout)
    
    def log(self, message):
        """Writes in status_bow the action with the time when the action happened."""
        # Get time when action happened
        timestamp = datetime.now().strftime("[%H:%M:%S]")
        # Write message in the box
        self.status_box.append(f"{timestamp} {message}")
    
    def clear(self):
        """Clears all lines written in the log book"""
        self.status_box.setText("")


class InputModule(QTabWidget):
    """Tab widget with that incorporates a SearchTab and a BlastTab.

    This module includes all sections where the user can input data to perform searches, analyse 
    barcodes and see results. The module is divided into two tabs: one to search barcodes and another
    one to analyse barcodes and get visual results.

    Attributes:
        controller (MainController): controller to pass down to each tab
    
    Widgets:
        search_tab(SearchTab): tab where the user can select the parameters to search for barcodes
        blast_tab(BlastTab): tab where the user can select the parameters to analyse a selected barcode, 
            display results and create a taxonomy tree.
    """
    def __init__(self, controller):
        """Creates an instance of the class.

        Args:
            controller (MainController): controller to pass down to each tab
        """
        super().__init__()
        # Add controller
        self.controller = controller
        # Set position tabs to horizontal on top
        self.setTabPosition(QTabWidget.North)
        """ Add tabs """
        self.search_tab = SearchTab(self.controller)
        self.blast_tab = BlastTab(self.controller)
        # Add tabs to tab list widget
        self.insertTab(0, self.search_tab, "SEARCH")
        self.insertTab(1, self.blast_tab, "BLAST")
    
    def clear(self):
        """Clean the parameters inputted by the user in the search tab"""
        self.search_tab.clear()


class Sequence(QTextEdit):
    """Text widget to display the barcode sequence and highlight it.

    This object is used inside barcode cards to display the sequence of the barcode and
    to highlight the primers whenver the mouse is hovered over the nucleotide sequence.

    Attributes:
        sequence(str): nucleotide string
        primers_intervals (list): intervals marking the positions where the nucleotide 
            sequence aligns with the primers.
    """
    def __init__(self, sequence, primers_intervals):
        """Creates an instance of the class.

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
        # Check all intervals
        for i, primer_list in enumerate(intervals):
            # Prepare cursor to highlight text
            cursor = self.textCursor()
            fmt = QTextCharFormat()
            fmt.setBackground(colours[i])
            # Highlight the selected interval
            for start, end in primer_list:
                if start > end:
                    continue  # skip invalid intervals
                length = end - start + 1
                # Apply highlight
                cursor.setPosition(start)
                cursor.movePosition(QTextCursor.Right, QTextCursor.KeepAnchor, length)
                cursor.mergeCharFormat(fmt)


class AlignmentPopup(QWidget):
    """Pop up window to display the alignment of two sequences.

    This pop-up window displays the alignment between two sequences. Each row contains
    60 bases. The alignment symbols are the following: 
        | : identical bases
        - : gap
        . : mismatch

    Widgets:
        alignment_text (QTextEdit): text box to display the alignment
    """
    def __init__(self, alignment, seq1_id, seq2_id):
        """Creates an instance of the class.

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

    Attributes:
        parent (BarcodesTab): BarcodesTab where the Barcode card belongs
        sequence (str): nucleotide sequence to display
        id (str): accesion number of the sequence to differentiate it from others and keep track
            of the selected one
        sequence_text (Sequence): Sequence object that is used to display highlighted primers

    Widgets:
        select_radio (QRadioButton): radio button to select barcode card
    """
    def __init__(self, parent, sequence, header, duplicates, primers_intervals):
        """Creates an instance of the class.

        Args:
            parent(BarcodesTab): see BarcodesTab class
            sequence(str): nucleotide sequence
            header(str): first header in the list of headers for a specific barcode
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
        header_label = QLabel(header)
        header_label.setObjectName("barcodeLabel")
        # Set the identifier of the card as the fasta file id
        self.id = header.split(" ")[0][1:]
        # Number of duplicates label
        duplicates_label = QLabel(f'Duplicates: {duplicates}')
        duplicates_label.setObjectName("barcodeLabel")
        # Length of sequence label
        length_label = QLabel(f'Length: {len(self.sequence)}')
        length_label.setObjectName("barcodeLabel")
        # Radio button to select card
        self.select_radio = QRadioButton()
        self.select_radio.clicked.connect(self.check)
        # Sequence text class
        self.sequence_text = Sequence(sequence, primers_intervals)
        self.sequence_text.setReadOnly(True)
        # Adjust size of widget
        self.sequence_text.document().adjustSize()
        height = self.sequence_text.document().size().height()
        self.sequence_text.setFixedHeight(int(height/3)+20)
        self.sequence_text.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.sequence_text.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        """Widget's layout."""
        # Define layout
        layout = QVBoxLayout()
        top_layout = QHBoxLayout()
        # Add fasta header to top_layout
        top_layout.addWidget(header_label)
        # Add stretch between fasta header and radio button
        top_layout.addStretch()
        # Add radio button
        top_layout.addWidget(self.select_radio)
        # Add top layout to main layout
        layout.addLayout(top_layout)
        # Add duplicates label
        layout.addWidget(duplicates_label)
        # Add length label
        layout.addWidget(length_label)
        # Add sequence text
        layout.addWidget(self.sequence_text) 
        # Set layout
        self.setLayout(layout)
    
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
        "Function to crete a pop up window with the alignment."
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
    """Tab to display the results obtained after the BLAST.

    This tab contains the initialisation for the chart widgets that will display the results obtained
    through the barcode analysis. Whenever stats are created or modified for the display the function
    add_stats updates the widgets with the new information.

    Widgets:
        result_container (QWidget): Widget that will sever as a container to other widgets
        layout (QVBoxLayout): layout to organise the widgets inside result_container
        all_hits_chart (CombinedChartWindow): widget that displays the stats for all hits
        species_chart (CombinedChartWindow): widget that displays the stats for only unique species

    """
    def __init__(self):
        """Creates an instance of the class."""
        super().__init__()
        # Allow to be resized
        self.setWidgetResizable(True)
        """Add widgets."""
        self.result_container = QWidget()
        self.result_container.setObjectName("outputWidget")
        self.all_hits_chart = None
        self.species_chart = None
        """Widgets' layout"""
        # Define layout
        self.layout = QVBoxLayout()
        # Set layout
        self.result_container.setLayout(self.layout)
        self.setWidget(self.result_container)
    
    def restart_animation(self):
        """Restarts the animation of the graphs."""
        if self.all_hits_chart and self.species_chart:
            self.all_hits_chart.restart_animation()
            self.species_chart.restart_animation()
    
    def clear_layout(self):
        """Deletes everything contained in the layout."""
        while self.layout.count():
            item = self.layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
    
    def add_stats(self, results, data_type, title):
        """Creates the widgets to display the stats obtained along with the title."""
        # Clean layout
        self.clear_layout()
        # Add title if it is first time displaying results
        if title:
            self.title_label = QLabel(title)
            self.title_label.setObjectName("resultsTitleLabel")
            self.title_label.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
            self.title_label.setAlignment(Qt.AlignHCenter)
        # Create widgets with charts
        labels = results[0]
        all_hits = results[1]
        unique_species = results[2]
        self.all_hits_chart = CombinedChartWindow(data=all_hits, 
                                                  labels=labels, 
                                                  title="BLAST"+"\n"+"HITS", 
                                                  data_type=data_type)
        self.species_chart = CombinedChartWindow(data=unique_species,
                                                 labels=labels,
                                                 title="UNIQUE"+"\n"+"SPECIES", 
                                                 data_type=data_type)
        # Add widgets
        self.layout.addWidget(self.title_label)
        self.layout.addWidget(self.all_hits_chart)
        self.layout.addWidget(self.species_chart)
        # Set layout
        self.result_container.setLayout(self.layout)
        self.setWidget(self.result_container)
           

class BarcodesTab(QScrollArea):
    """Tab to display the barcode cards.

    Attributes:
        checked_card (BarcodeCard): card selected by the user
        controller (MainController): controller to pass the selected and unselected barcodes
        barcodes_data (list): list containing the data for all barcode cards to be displayed in the tab
        tab_name (str): title of the tab
        barcode_card_list (list): list of barcodes contained in the tab    
    """
    def __init__(self, controller, tab_name, barcodes_data):
        """Creates an instance of the class.

        Args:
            controller(MainController): controller to pass the selected and unselected barcodes.
            tab_name (str): name of the tab, which is the type of barcode selected
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
        result_container = QWidget()
        result_container.setObjectName("outputWidget")
        result_layout = QVBoxLayout()
        for barcode_card in self.barcode_card_list:
            barcode_card.setSizePolicy(QSizePolicy.Preferred, QSizePolicy.Fixed)
            result_layout.addWidget(barcode_card)
        # Add expanding spacer to push widgets up and prevent stretching
        spacer = QWidget()
        spacer.setSizePolicy(QSizePolicy.Minimum, QSizePolicy.Expanding)
        result_layout.addWidget(spacer)       
        result_container.setLayout(result_layout)
        self.setWidget(result_container)
    
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
        # Looks for the card that has been selected
        for barcode_card in self.barcode_card_list:
            if barcode_card.id == barcode_id:
                # Check the barcode and add that sequence for further analysis
                self.checked_card = barcode_card
                self.controller.add_sequence_to_blast(sequence=barcode_card.sequence, 
                                                      barcode_type=self.tab_name, 
                                                      acc_number=barcode_id)
        self.controller.write_in_logbook("The barcode " + barcode_id + " has been selected.")

    def uncheck_barcode_cards(self, barcode_id):
        """Deletes a barcode card that has been unchecked from the variable checked_card given its id."""
        # Look for the barcodes card that has been deselected
        for barcode_card in self.barcode_card_list:
            if barcode_card.id != barcode_id:
                # Uncheck the card
                barcode_card.uncheck()
        # Remove card from checked_card and remove sequence from the list of sequences to analyse
        self.checked_card = None
        self.controller.remove_sequence_to_blast(self.tab_name)

    
class OutputModule(QTabWidget):
    """Module for the output data: barcode cards and results.

    This module contains the tabs with the barcode cards found and the results tab where the
    statistical results are shown. Initially there is only an empty results tab and barcode tabs
    are added as the user searches specific barcodes. 

    Attributes:
        controller(MainController): controller to pass down when addint tabs to the module
    
    Widgets:
        results_tab(ResultsTab): see ResultsTab
    """
    def __init__(self, controller):
        """Creates an instance of the class.

        Args:
            controller(MainController): controller to pass down when addint tabs to the module
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
        self.currentChanged.connect(self.on_tab_changed)
    
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
    
    def show_results(self, results, title):
        "Passes down the results and title from the analysis to display them in the results tab."
        # Check type of results
        labels = results[0]
        data_type = "percentages" if '%' in labels[0] else 'differences'
        # Write action in logbook
        self.controller.write_in_logbook(f'Displaying results by {data_type}')
        # Change tab to results tab
        self.setCurrentIndex(0) 
        # Add results in results tab
        self.results_tab.add_stats(results, data_type, title)     
    
    def on_tab_changed(self, index):
        """Resstarts the animation when the tab is changed to results tab for aesthetic purposes."""
        if index == 0:
            self.results_tab.restart_animation()
    
    def clear(self):
        """Clears al the information in the module to start a new search"""
        # Set tab to results tab
        self.setCurrentIndex(0)
        # Delete all barcode tabs
        while self.count() > 1:
            self.removeTab(1)
        # Delete results tab info
        self.results_tab.clear_layout()
          