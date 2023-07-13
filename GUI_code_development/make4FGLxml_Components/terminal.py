import sys

from PyQt6.QtWidgets import (
    QWidget,
    QTextEdit,
    QHBoxLayout)

from PyQt6.QtGui import QTextOption,QTextCursor

class TerminalWidget(QWidget):
    '''
    A class for creating and managing a "terminal" widget in the
    make4FGLxml GUI (this class can, actually, be applied more
    generally and is not limited to the make4FGL GUI).
    This captures the stdout and stderr and prints
    to a window in the GUI.

    ...

    Attributes
    ----------
    <inherits some properties from the PyQt6 QWidget class,
    see that documentation for more details>
    main_window : PyQt6 QMainWindow object
        the main window in which this widget will be contained
    terminal_layout : PyQt6 QHBoxLayout object
        layout object for the terminal
    terminal_out : PyQt6 QTextEdit object
        the widget handing the display of text from stdout and stderr

    Methods
    -------
    flush()
        method for stdout and stderr compatability
    write(text)
        method to insert text received from stdout and stderr
        and move the cursor to the end of the line
    '''
    
    def __init__(self,main_window):
        '''
        Parameters
        ----------
        main_window - PyQt6 QMainWindow object
            the main window within which this widget is contained
        '''

        #call init for the base class
        super().__init__(main_window)

        self.main_window=main_window

        self.terminal_out=QTextEdit(parent=self)
        
        
        #set the word wrap to try and wrap at word boundaries
        #but allow it to wrap anywere in a line if necessary
        self.terminal_out.setWordWrapMode(QTextOption().WrapMode(4))

        #start things out as readOnly after adding some intro text
        #figure out a way to access the version number without hardcoding it here
        self.terminal_out.setMarkdown('## Welcome to the make4FGLxml GUI\n---')
        self.terminal_out.append('')
        self.terminal_out.moveCursor(QTextCursor().MoveOperation(11))
        self.terminal_out.setReadOnly(True)

        #set the layout, ensure that widget extends to edge of window
        self.terminal_layout=QHBoxLayout()
        self.terminal_layout.setContentsMargins(0,0,0,0)

        self.terminal_layout.addWidget(self.terminal_out)

        self.setLayout(self.terminal_layout)

    def write(self,text):
        '''
        method to insert text into the QTextEdit instance and display

        Parameters
        ----------
        text - str
            text received from stderr or stdout to display
        '''
        
        #turn off the read only property
        self.terminal_out.setReadOnly(False)

        #insert text, hopefully at the end of the line, whic is what the
        #scrollToAnchor command hopefully does
        self.terminal_out.insertPlainText(text)

        #we'll try and make sure we can get to the last by
        #splitting on the white space character and grabbing the last few
        self.terminal_out.moveCursor(QTextCursor().MoveOperation(11))
        
        #then turn readonly back on
        self.terminal_out.setReadOnly(True)

    def flush(self):
        '''
        method for compatability when redirecting stderr and stdout
        '''

        #do nothing
        pass
        
        
