"""Contains the main creation of the main window and the controller. It has the main function to run the app"""
import sys
from PyQt5.QtWidgets import QApplication
from PyQt5.QtGui import QIcon
from gui.main_window import MainWindow
from controllers.MainController import MainController


def main():
    # Initialise controller
    bloom_controller = MainController()
    # Initialise application
    app = QApplication([])
    # Set icon
    app.setWindowIcon(QIcon(bloom_controller.get_icon()))
    # Load style sheet
    app.setStyleSheet(bloom_controller.load_style())
    # Initialise main window
    window = MainWindow(bloom_controller, app)
    # Run app
    window.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
    