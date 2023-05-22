#!/usr/bin/env python

#using script example from pyinstaller website
import PyInstaller.__main__
import sys

#honestly, not sure if this will work here, but if you
#get output messages from pyinstaller, add the recursion limit
#line below to your .spec file and follow the instructions
#to create the exe from the modified .spec file
sys.setrecursionlimit(5000)

PyInstaller.__main__.run([
    'make4FGLxml_GUI.py',
    '--onefile',
    '--windowed',
    '--noconsole'])
