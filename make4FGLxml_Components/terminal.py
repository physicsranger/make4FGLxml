import sys, os

from PyQt6.QtWidgets import (
    QWidget,
    QTextEdit,
    QHBoxLayout)
from PyQt6.QtGui import QTextOption,QTextCursor

class TerminalWidget(QWidget):
    def __init__(self,main_window):
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

        self.terminal_layout=QHBoxLayout()
        self.terminal_layout.setContentsMargins(0,0,0,0)

        self.terminal_layout.addWidget(self.terminal_out)

        self.setLayout(self.terminal_layout)

    def write(self,text):
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
        pass
        
        
