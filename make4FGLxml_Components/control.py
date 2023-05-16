import sys,os

from PyQt6.QtWidgets import (
    QWidget,
    QPushButton,
    QCheckBox,
    QSpinBox,
    QLineEdit,
    QVBoxLayout,
    QHBoxLayout,
    QLabel,
    QFileDialog,
    QApplication)

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QDoubleValidator

class ControlWidget(QWidget):
    def __init__(self,main_window):
        super().__init__(main_window)

        self.main_window=main_window
        
        self.fermi_dir=os.environ.get('FERMI_DIR')

        self.create_widgets()

    def create_widgets(self):
        self.create_first_row()
        self.create_second_row()
        self.create_third_row()
        self.create_fourth_row()
        self.create_fifth_row()
        self.create_sixth_row()
        self.create_seventh_row()
        self.create_eighth_row()
        self.create_ninth_row()
        self.create_tenth_row()
        self.create_eleventh_row()
        self.create_twelvth_row()

    def create_first_row(self):
        self.first_row=QWidget(parent=self)
        #create widgets for inputing catalog file
        self.catalog_button=QPushButton("Catalog",parent=self.first_row)
        self.catalog_button.setCheckable(True)
        self.catalog_button.clicked.connect(self.choose_catalog)
        self.catalog_button.setFixedWidth(100)
        
        self.catalog_entry=QLineEdit(parent=self.first_row)
        self.catalog_entry.setFixedWidth(400)

        #create widgets for the DR version
        self.DR_label=QLabel('DR:',parent=self.first_row)
        self.DR_spinbox=QSpinBox(parent=self.first_row)
        self.DR_spinbox.setValue(3)
        self.DR_spinbox.setFixedWidth(30)

        #set reasonable limits on the spinbox
        self.DR_spinbox.setMinimum(1)
        self.DR_spinbox.setMaximum(3)

    def create_second_row(self):
        self.second_row=QWidget(parent=self)

        #create widgets relating to output file
        self.save_directory_button=QPushButton("Save Directory",parent=self.second_row)
        self.save_directory_button.setCheckable(True)
        self.save_directory_button.clicked.connect(self.choose_save_directory)
        self.save_directory_entry=QLineEdit(os.getcwd(),parent=self.second_row)
        self.save_directory_button.setFixedWidth(100)
        self.save_directory_entry.setFixedWidth(400)

        self.output_file_label=QLabel('Output File:',parent=self.second_row)
        self.output_file_entry=QLineEdit('my_model.xml',parent=self.second_row)

    def create_third_row(self):
        self.third_row=QWidget(parent=self)
        
        #create widgets for the region file
        self.region_check=QCheckBox("Make ds9 region file?",parent=self.third_row)
        self.region_entry=QLineEdit('my_region.reg',parent=self.third_row)
        ##set up to tie to default above but do we need to do this?
        ##do we want to have this always match the output_file text? might annoy
        ##users who edit this and then make a change to the xml name

    def create_fourth_row(self):
        self.fourth_row=QWidget(parent=self)
        
        #create widgets for ROI center
        self.ROI_center_label=QLabel("ROI information [deg.]:",parent=self.fourth_row)
        self.ROI_center_RA_label=QLabel("RA =",parent=self.fourth_row)
        self.ROI_center_DEC_label=QLabel("DEC =",parent=self.fourth_row)
        self.ROI_radius_label=QLabel("radius =",parent=self.fourth_row)
        
        #set up the entry widgets with constrained input
        self.RA_validator=QDoubleValidator(bottom=0,top=360)
        self.DEC_validator=QDoubleValidator(bottom=-90,top=90)
        self.radius_validator=QDoubleValidator(bottom=0)
        
        self.ROI_center_RA_entry=QLineEdit(parent=self.fourth_row)
        self.ROI_center_RA_entry.setValidator(self.RA_validator)
        self.ROI_center_RA_entry.setFixedWidth(50)
        
        self.ROI_center_DEC_entry=QLineEdit(parent=self.fourth_row)
        self.ROI_center_DEC_entry.setValidator(self.DEC_validator)
        self.ROI_center_DEC_entry.setFixedWidth(50)

        self.ROI_radius_entry=QLineEdit(parent=self.fourth_row)
        self.ROI_radius_entry.setValidator(self.radius_validator)
        self.ROI_radius_entry.setFixedWidth(25)
        #now need self.RA and self.DEC tied to the QLineEdit widgets above

        #allow user to instead choose an event file to get the ROI center information
        self.or_label=QLabel(' or ',parent=self.fourth_row)
        self.use_event_file_button=QPushButton("Get From Event File",parent=self.fourth_row)
        self.use_event_file_button.setCheckable(True)
        self.use_event_file_button.clicked.connect(self.get_ROI_center)
        self.use_event_file_button.setFixedWidth(125)
        
        self.event_file_label=QLabel("Event File:",parent=self.fourth_row)
        self.event_file_entry=QLineEdit(parent=self.fourth_row)

    def create_fifth_row(self):
        self.fifth_row=QWidget(parent=self)

        #now create widgets for optional parameters        
        self.free_radius_label=QLabel("Radius Limit For Most Free Sources:",parent=self.fifth_row)
        self.free_radius_validator=QDoubleValidator(bottom=-1)
        self.free_radius_entry=QLineEdit('-1',parent=self.fifth_row)
        self.free_radius_entry.setFixedWidth(25)
        self.free_radius_entry.setValidator(self.free_radius_validator)

        self.max_free_radius_label=QLabel("Radius Limit For All Free Sources:",parent=self.fifth_row)
        self.max_free_radius_validator=QDoubleValidator(bottom=-1)
        self.max_free_radius_entry=QLineEdit('-1',parent=self.fifth_row)
        self.max_free_radius_entry.setValidator(self.max_free_radius_validator)
        self.max_free_radius_entry.setFixedWidth(25)

        self.extra_radius_label=QLabel("Radial Distance Beyond ROI To Include Fixed Sources:",parent=self.fifth_row)
        self.extra_radius_validator=QDoubleValidator(bottom=0,top=30)
        self.extra_radius_entry=QLineEdit('10',parent=self.fifth_row)
        self.extra_radius_entry.setValidator(self.extra_radius_validator)
        self.extra_radius_entry.setFixedWidth(25)

    def create_sixth_row(self):
        self.sixth_row=QWidget(parent=self)
        
        #would be good to have a 'trace' set on this to modify based on XML or FITS
        #version of the catalog (e.g., avg sig vs avg TS)
        self.significance_label=QLabel("Minimum Average Significance In Catalog For Free Sources:",parent=self.sixth_row)
        self.significance_validator=QDoubleValidator(bottom=0)
        self.significance_entry=QLineEdit(parent=self.sixth_row)
        self.significance_entry.setValidator(self.significance_validator)
        self.significance_entry.setFixedWidth(25)

    def create_seventh_row(self):
        self.seventh_row=QWidget(parent=self)
        
        #now for option check boxes
        self.normalizations_only_check=QCheckBox("Only Free Normalization Parameters?",parent=self.seventh_row)
        self.variable_sources_free_check=QCheckBox("Free Variable Sources?",parent=self.seventh_row)
        self.force_point_source_check=QCheckBox("Change Extended Sources To Point Sources?",parent=self.seventh_row)

    def create_eighth_row(self):
        self.eighth_row=QWidget(parent=self)
        
        self.use_extended_catalog_names=QCheckBox("Use Catalog Names For Extended Sources?",parent=self.eighth_row)
        self.free_galactic_index_check=QCheckBox("Modify Galactic Spectrum By Power Law?",parent=self.eighth_row)
        self.use_old_name_convention_check=QCheckBox("Use '_4FGLJXXXX.X+XXXX' Naming Convention?",parent=self.eighth_row)

    def create_ninth_row(self):
        self.ninth_row=QWidget(parent=self)
        
        #now for diffuse files/directories and names
        if self.fermi_dir is None:
            self.galactic_model_file_entry=QLineEdit(os.sep.join(['$(FERMI_DIR)','refdata','fermi',
                    'galdiffuse','gll_iem_v07.fits']),parent=self.ninth_row)
        
        else:
            self.galactic_model_file_entry=QineEdit(os.path.join(self.fermi_dir,'refdata','fermi',
                    'galdiffuse','gll_iem_v07.fits'),parent=self.ninth_row)

        self.galactic_model_file_entry.setFixedWidth(400)

        self.galactic_model_file_label=QLabel('File Path For Galactic Diffuse Model:',parent=self.ninth_row)
        self.galactic_model_file_button=QPushButton("Change File",parent=self.ninth_row)
        self.galactic_model_file_button.setCheckable(True)
        self.galactic_model_file_button.clicked.connect(self.choose_galactic_model_file)
        self.galactic_model_file_button.setFixedWidth(100)

    def create_tenth_row(self):
        self.tenth_row=QWidget(parent=self)

        #now for diffuse files/directories and names
        if self.fermi_dir is None:
            self.isotropic_template_file_entry=QLineEdit(os.sep.join(['$(FERMI_DIR)','refdata','fermi',
                    'galdiffuse','iso_P8R3_SOURCE_V3_v1.txt']),parent=self.tenth_row)
        
        else:
            self.isotropic_template_file_entry=QLineEdit(os.path.join(self.fermi_dir,'refdata','fermi',
                    'galdiffuse','iso_P8R3_SOURCE_V3_v1.txt'),parent=self.tenth_row)

        self.isotropic_template_file_entry.setFixedWidth(400)

        self.isotropic_template_file_label=QLabel('File Path For Isotropic Template:',parent=self.tenth_row)
        self.isotropic_template_file_button=QPushButton("Change File",parent=self.tenth_row)
        self.isotropic_template_file_button.setCheckable(True)
        self.isotropic_template_file_button.clicked.connect(self.choose_isotropic_template_file)
        self.isotropic_template_file_button.setFixedWidth(100)

    def create_eleventh_row(self):
        self.eleventh_row=QWidget(parent=self)
        
        #now for diffuse files/directories and names
        if self.fermi_dir is None:
            self.extended_template_directory_entry=\
                    QLineEdit(os.sep.join(['$(FERMI_DIR)','data','pyBurstAnalysisGUI','templates']),
                              parent=self.eleventh_row)
        
        else:
            self.extended_template_directory_entry=\
                    QineEdit(os.path.join(self.fermi_dir,'data','pyBurstAnalysisGUI','template'),
                             parent=self.eleventh_row)

        self.extended_template_directory_entry.setFixedWidth(400)

        self.extended_template_directory_label=QLabel('Directory With Extended Source Templates:',
                                                      parent=self.eleventh_row)
        self.extended_template_directory_button=QPushButton("Change Directory",
                                                            parent=self.eleventh_row)
        self.extended_template_directory_button.setCheckable(True)
        self.extended_template_directory_button.clicked.connect(self.choose_extended_directory)
        self.extended_template_directory_button.setFixedWidth(125)

    def create_twelvth_row(self):
        self.twelvth_row=QWidget(parent=self)
        
        #now we need a button to actually make the files
        self.make_model_button=QPushButton('Make Model',parent=self.twelvth_row)
        self.make_model_button.setCheckable(True)
        self.make_model_button.clicked.connect(self.create_model)
        
        #and finish it off with a 'Quit' button
        self.quit_button=QPushButton("Quit",parent=self.twelvth_row)
        self.quit_button.setCheckable(True)
        self.quit_button.clicked.connect(QApplication.instance().quit)

        self.make_model_button.setFixedWidth(100)
        self.quit_button.setFixedWidth(100)

    
    def pack_widgets(self):
        #We're going to stack rows vertically, so make a main layout
        #and then pack the widgets in 'rows', with only related widgets in each row
        #so each row will need to be a horizontal layout
        self.main_layout=QVBoxLayout()
        #self.main_layout.setSpacing(0)

        self.pack_first_row()
        self.pack_second_row()
        self.pack_third_row()
        self.pack_fourth_row()
        self.pack_fifth_row()
        self.pack_sixth_row()
        self.pack_seventh_row()
        self.pack_eighth_row()
        self.pack_ninth_row()
        self.pack_tenth_row()
        self.pack_eleventh_row()
        self.pack_twelvth_row()
        
        self.main_layout.addWidget(self.first_row)
        self.main_layout.addWidget(self.second_row)
        self.main_layout.addWidget(self.third_row)
        self.main_layout.addWidget(self.fourth_row)
        self.main_layout.addWidget(self.fifth_row)
        self.main_layout.addWidget(self.sixth_row)
        self.main_layout.addWidget(self.seventh_row)
        self.main_layout.addWidget(self.eighth_row)
        self.main_layout.addWidget(self.ninth_row)
        self.main_layout.addWidget(self.tenth_row)
        self.main_layout.addWidget(self.eleventh_row)
        self.main_layout.addWidget(self.twelvth_row)

        #self.main_layout.setSpacing(0)
        
        self.setLayout(self.main_layout)

    def pack_first_row(self):
        self.first_row_layout=QHBoxLayout()
        self.first_row_layout.setContentsMargins(0,0,0,0)
        
        #pack the first row, catalog file and DR stuff
        self.first_row_layout.addWidget(self.catalog_button)
        self.first_row_layout.addWidget(self.catalog_entry)
        self.first_row_layout.addWidget(self.DR_label)
        self.first_row_layout.addWidget(self.DR_spinbox)
        self.first_row_layout.addStretch(1)

        self.first_row.setLayout(self.first_row_layout)

    def pack_second_row(self):
        self.second_row_layout=QHBoxLayout()
        self.second_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the second row, save directory and output fil stuff
        self.second_row_layout.addWidget(self.save_directory_button)
        self.second_row_layout.addWidget(self.save_directory_entry)
        self.second_row_layout.addWidget(self.output_file_label)
        self.second_row_layout.addWidget(self.output_file_entry)

        self.second_row.setLayout(self.second_row_layout)
        
    def pack_third_row(self):
        self.third_row_layout=QHBoxLayout()
        self.third_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the third row, region file stuff
        self.third_row_layout.addWidget(self.region_check)
        self.third_row_layout.addWidget(self.region_entry)
        self.third_row_layout.addStretch(1)

        self.third_row.setLayout(self.third_row_layout)
        
    def pack_fourth_row(self):
        self.fourth_row_layout=QHBoxLayout()
        self.fourth_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the fourth row, ROI information stuff
        self.fourth_row_layout.addWidget(self.ROI_center_label)

        self.fourth_row_layout.addWidget(self.ROI_center_RA_label)
        self.fourth_row_layout.addWidget(self.ROI_center_RA_entry)

        self.fourth_row_layout.addWidget(self.ROI_center_DEC_label)
        self.fourth_row_layout.addWidget(self.ROI_center_DEC_entry)

        self.fourth_row_layout.addWidget(self.ROI_radius_label)
        self.fourth_row_layout.addWidget(self.ROI_radius_entry)

        self.fourth_row_layout.addStretch(1)
        self.fourth_row_layout.addWidget(self.or_label)
        self.fourth_row_layout.addStretch(1)
        
        self.fourth_row_layout.addWidget(self.use_event_file_button)
        self.fourth_row_layout.addWidget(self.event_file_label)
        self.fourth_row_layout.addWidget(self.event_file_entry)

        self.fourth_row.setLayout(self.fourth_row_layout)

    def pack_fifth_row(self):
        self.fifth_row_layout=QHBoxLayout()
        self.fifth_row_layout.setContentsMargins(0,0,0,0)

        #grid the fifth row, the free source radius limit stuff
        self.fifth_row_layout.addWidget(self.free_radius_label)
        self.fifth_row_layout.addWidget(self.free_radius_entry)
        self.fifth_row_layout.addWidget(self.max_free_radius_label)
        self.fifth_row_layout.addWidget(self.max_free_radius_entry)
        self.fifth_row_layout.addWidget(self.extra_radius_label)
        self.fifth_row_layout.addWidget(self.extra_radius_entry)
        self.fifth_row_layout.addStretch(1)
        
        self.fifth_row.setLayout(self.fifth_row_layout)

    def pack_sixth_row(self):
        self.sixth_row_layout=QHBoxLayout()
        self.sixth_row_layout.setContentsMargins(0,0,0,0)

        #grid the sixth row, minimum significance stuff
        self.sixth_row_layout.addWidget(self.significance_label)
        self.sixth_row_layout.addWidget(self.significance_entry)
        self.sixth_row_layout.addStretch(1)

        self.sixth_row.setLayout(self.sixth_row_layout)

    def pack_seventh_row(self):
        self.seventh_row_layout=QHBoxLayout()
        self.seventh_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the seventh row, normalization and variability checkboxes
        self.seventh_row_layout.addWidget(self.normalizations_only_check)
        self.seventh_row_layout.addWidget(self.variable_sources_free_check)
        self.seventh_row_layout.addWidget(self.force_point_source_check)
        self.seventh_row_layout.addStretch(1)

        self.seventh_row.setLayout(self.seventh_row_layout)

    def pack_eighth_row(self):
        self.eighth_row_layout=QHBoxLayout()
        self.eighth_row_layout.setContentsMargins(0,0,0,0)

        #grid the eighth row, extended source related checkboxes
        self.eighth_row_layout.addWidget(self.use_extended_catalog_names)
        self.eighth_row_layout.addWidget(self.free_galactic_index_check)
        self.eighth_row_layout.addWidget(self.use_old_name_convention_check)

        self.eighth_row.setLayout(self.eighth_row_layout)

    def pack_ninth_row(self):
        self.ninth_row_layout=QHBoxLayout()
        self.ninth_row_layout.setContentsMargins(0,0,0,0)

        #grid the ninth row, Galactic diffuse model stuff
        self.ninth_row_layout.addWidget(self.galactic_model_file_label)
        self.ninth_row_layout.addWidget(self.galactic_model_file_entry)
        self.ninth_row_layout.addWidget(self.galactic_model_file_button)
        self.ninth_row_layout.addStretch(1)

        self.ninth_row.setLayout(self.ninth_row_layout)

    def pack_tenth_row(self):
        self.tenth_row_layout=QHBoxLayout()
        self.tenth_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the tenth row, Isotropic diffuse stuff
        self.tenth_row_layout.addWidget(self.isotropic_template_file_label)
        self.tenth_row_layout.addWidget(self.isotropic_template_file_entry)
        self.tenth_row_layout.addWidget(self.isotropic_template_file_button)
        self.tenth_row_layout.addStretch(1)

        self.tenth_row.setLayout(self.tenth_row_layout)

    def pack_eleventh_row(self):
        self.eleventh_row_layout=QHBoxLayout()
        self.eleventh_row_layout.setContentsMargins(0,0,0,0)

        #grid the eleventh row, Extended source template directory stuff
        #and quit button
        self.eleventh_row_layout.addWidget(self.extended_template_directory_label)
        self.eleventh_row_layout.addWidget(self.extended_template_directory_entry)
        self.eleventh_row_layout.addWidget(self.extended_template_directory_button)
        self.eleventh_row_layout.addStretch(1)

        self.eleventh_row.setLayout(self.eleventh_row_layout)

    def pack_twelvth_row(self):
        self.twelvth_row_layout=QHBoxLayout()
        self.twelvth_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the twelvth row, buttons to generate the file(s) and to quit
        self.twelvth_row_layout.addWidget(self.make_model_button)
        self.twelvth_row_layout.addStretch(1)
        self.twelvth_row_layout.addWidget(self.quit_button)

        self.twelvth_row.setLayout(self.twelvth_row_layout)

    def choose_catalog(self):
        #do stuff to open a file browser so the user can choose the catalog file
        #and then assign the choice to the catalog_entry and catlog value
        return
        
    def choose_save_directory(self):
        #do stuff to open a file browser and choose a directory to save the output
        #xml model (and possible .reg file) in
        return

    def choose_extended_directory(self):
        #do stuff to open a file browser and specify the directory with the extended source files
        return

    def choose_galactic_model_file(self):
        #do stuff to open a file browser and choose a file for the Galactic diffuse model
        return

    def choose_isotropic_template_file(self):
        #do stuff to open a file browser and choose a file for the isotropic diffuse template
        return

    def get_ROI_center(self):
        #do stuff to open a file browser and choose an event file to get the ROI
        #center information from, and get that info and fill the LineEdit widgets
        return

    def create_model(self):
        #do stuff to make the XML file (and possible .reg file)
        return

