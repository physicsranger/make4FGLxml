#!/usr/bin/env python

import sys
from PyQt6.QtWidgets import (
    QApplication,
    QMainWindow,
    QVBoxLayout,
    QWidget)

#custom imports
from make4FGLxml_Components.control import ControlWidget
from make4FGLxml_Components.terminal import TerminalWidget

#class for the main window of the GUI
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle('make4FGLxml')

        self.control_widget=ControlWidget(self)

        self.terminal_widget=TerminalWidget(self)

        self.vbox_layout=QVBoxLayout()

        self.control_widget.pack_widgets()

        self.vbox_layout.addWidget(self.control_widget)
        self.vbox_layout.addWidget(self.terminal_widget)

        self.main_widget=QWidget()
        self.main_widget.setLayout(self.vbox_layout)

        self.setCentralWidget(self.main_widget)

        #now we want to redirect output from the 'terminal'
        #or any print statement to our terminal widget
        self.old_stdout=sys.stdout
        self.old_stderr=sys.stderr

        sys.stdout=self.terminal_widget
        sys.stderr=self.terminal_widget


#I'll likely need to tweak things a bit, but this should work
def main():
    app=QApplication()
    
    window=MainWindow()
    window.show()
    
    app.exec()

if __name__=='main':
    main()
