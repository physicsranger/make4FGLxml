import sys,os,warnings

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

from PyQt6.QtCore import QTimer

from PyQt6.QtGui import QDoubleValidator

from build_model.utilities import get_ROI_from_event_file

from make4FGLxml import SourceList

class ControlWidget(QWidget):
    '''
    A class to encompass and manage the interactive widgets
    of the make4FGLxml GUI

    ...

    Attributes
    ----------
    <inherits some properties from the PyQt6 QWidget class,
    see that documentation for more details>
    catalog_button : PyQt6 QPushButton
        button to open a file browser and select 4FGL catalog
        file to select sources from
    catalog_entry : PyQt6 QLineEdit
        entry widget for typing in catalog file name, will also
        display full path to file selected via catalog_button
    DEC_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the value given to
        ROI_center_DEC_entry
    DR_label : PyQt6 QLabel
        label for DR_spinbox
    DR_spinbox : PyQt6 QSpinBox
        spinbox to select 4FGL data release version
    eighth_row: PyQt6 Qwidget object
        eighth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    eighth_row_layout : PyQt6 QHBoxLayout
        layout for the eighth row
    eleventh_row : PyQt6 Qwidget object
        eleventh row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    eleventh_row_layout : PyQt6 QHBoxLayout
        layout for the eleventh row
    event_file_entry : PyQt6 QLineEdit
        widget to enter the name of the LAT event file from which
        the region of interest information will be derived, will
        also display the path to the file chosen via the
        use_event_file_button
    event_file_label : PyQt6 QLabel
        label for the event_file_entry widget
    extended_template_directory_button : PyQt6 QPushButton
        button to open a directory browser and specify the directory
        with extended source templates
    extended_template_directory_entry : PyQt6 QLineEdit
        widget to enter the path to the directory where the extended source
        template files are stored
    extended_template_directory_label : PyQt6 QLabe
        label for the extended_template_directory_entry
    extra_radius_entry : PyQt6 QLineEdit
        widget to enter the radius beyond the region of interest out to
        which sources are included to account for the size of the
        low-energy point-spread function
    extra_radius_label : PyQt6 QLabel
        label for the extra_radius_entry widget
    extra_radius_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the value given to
        the extra_radius_entry
    fermi_dir : str
        path to the fermitools installation, accessed through
        the FERMI_DIR environment variable
    fifth_row : PyQt6 Qwidget object
        fifth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    fifth_row_layout : PyQt6 QHBoxLayout
        layout for the fifth row
    first_row : PyQt6 Qwidget object
        first row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    first_row_layout : PyQt6 QHBoxLayout
        layout for the first row
    force_point_source_check : PyQt6 QCheckBox
        checkbox to indicate if extended sources should be cast as
        point sources
    fourth_row : PyQt6 Qwidget object
        fourth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    fourth_row_layout : PyQt6 QHBoxLayout
        layout for the fourth row
    free_galactic_index_check : PyQt6 QCheckBox
        check box to indicate if the Galactic diffuse component spectrum
        should be modified by a power law with free spectral index
    free_radius_entry : PyQt6 QLineEdit
        widget to enter the radius limit value for sources
        to be free
    free_radius_label : PyQt6 QLabel
        label for the free_radius_entry
    free_radius_validator : PyQt6 QDoubleValidator
        validator to enforce constraints on the value provided
        to the free_radius_entry widget
    galactic_model_file_button : PyQt6 QPushButton
        button to open a file browser and select the file for the
        Galactic diffuse emission
    galactic_model_file_entry : PyQt6 QLineEdit
        widget for the path to the file for the Galactic diffuse emission,
        will display file chosen by galactic_model_file_button
    galactic_model_file_label : PyQt 6 QLabel
        label for the galactic_model_file_entry widget
    galactic_name_entry : PyQt 6QLineEdit
        widget for the name of the Galactic diffuse emission
        component in the model
    galactic_name_label : PyQt6 QLabel
        label for the galactic_name_entry
    isotropic_name_entry : PyQt6 QLineEdit
        widget to enter the name of the isotropic diffuse component
        in the model
    isotropic_name_label : PyQt6 QLabel
        label for the isotropic_name_entry
    isotropic_template_file_button : PyQt6 QPushButton
        button to open a file browser and select the template file for
        the isotropic diffuse component
    isotropic_template_file_entry : PyQt6 QLineEdit
        widget for the isotropic diffuse component template file
    isotropic_template_file_label : PyQ6 QLabel
        label for the isotropic_template_file_entry
    main_layout : PyQt6 QVBoxLayout
        overall layout for the ControlWidget
    main_window : PyQt6 QMainWindow object
        the main window in which this widget will be contained
    make_model_button : PyQt6 QPushButton
        button to use the specified inputs to generate a source model file
    max_free_radius_entry : PyQt6 QLineEdit
        widget to enter the absoulte maximum radius for free sources
    max_free_radius_label : PyQt6 QLabel
        label for max_free_radius_entry widget
    max_free_radius_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the value given to
        max_free_radius_entry
    normalizations_only_check : PyQt6 QCheckBox
        checkbox to indicate if only source normalizations can possibly
        be set free or not
    ninth_row : PyQt6 Qwidget object
        ninth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    ninth_row_layout : PyQt6 QHBoxLayout
        layout for the ninth row
    or_label : PyQt6 QLabel
        label which simply says 'or' with white space around to
        separate direct region of interest widgets from widgets
        to get the region of interest information from a LAT
        event file
    output_file_entry : PyQt6 QLineEdit
        widget to enter name of output model file
    output_file_label : PyQt6 QLabel
        label widget for the output_file_entry widget
    quit_button : PyQt6 QPushButton
        button to restore stderr and stdout and quit the GUI
    RA_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the value given to
        the ROI_center_RA_entry
    radius_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the value given to
        the ROI_radius_entry
    region_check : PyQt6 QCheckBox
        check box to make or not make a ds9-style .reg file
    region_entry : PyQt6 QLineEdit
        widget to enter name of output ds9-style .reg file,
        if region_check is checked
    ROI_center_DEC_entry : PyQt6 QLineEdit
        widget to enter declination for region of interest center
        if not using a LAT event file
    ROI_center_DEC_label : PyQt6 QLabel
        label for the ROI_center_DEC_entry
    ROI_center_label : PyQt6 QLabel
        label for all region of interest related widgets
    ROI_center_RA_entry : PyQt6 QLineEdit
        widget to enter right ascension for region of interest center
        if not using a LAT event file
    ROI_center_RA_label : PyQt6 QLabel
        label for the ROI_center_RA_entry widget
    ROI_radius_entry : PyQt6 QLineEdit
        widget to enter the radius of the region of interest, if
        not using a LAT event file
    ROI_radius_label : PyQt6 QLabel
        label for the ROI_radius_entry widget
    save_directory_button : PyQt6 QPushButton
        button to open a directory browser and select the target
        directory to save source model file in
    save_directory_entry : PyQt6 QLineEdit
        widget to enter target directory for files, will also
        display path chosen with save_directory_button
    second_row : PyQt6 Qwidget object
        second row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    second_row_layout : PyQt6 QHBoxLayout
        layout for the second row
    seventh_row : PyQt6 Qwidget object
        seventh row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    seventh_row_layout : PyQt6 QHBoxLayout
        layout for the seventh row
    significance_entry : PyQt6 QLineEdit
        widget to enter the significance threshold for free sources
    significance_label : PyQt6 QLabel
        label for significance_entry widget
    significance_validator : PyQt6 QDoubleValidator
        validator to enforce some constraints on the values given to
        significance_entry
    sixth_row : PyQt6 Qwidget object
        sixth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    sixth_row_layout : PyQt6 QHBoxLayout
        layout for the sixth row
    tenth_row : PyQt6 Qwidget object
        tenth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    tenth_row_layout : PyQt6 QHBoxLayout
        layout for the tenth row
    third_row : PyQt6 Qwidget object
        third row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    third_row_layout : PyQt6 QHBoxLayout
        layout for the third row
    twelvth_row : PyQt6 Qwidget object
        twelvth row of widgets, packed horizontally, combined
        with other 'rows' stacked vertically to achieve
        desired layout
    twelvth_row_layout : PyQt6 QHBoxLayout
        layout for the twelvth row
    use_event_file_button : PyQt6 QPushButton
        button to open a file browser and choose a LAT event
        file from which the region of interest information will
        be obtained
    use_extended_catalog_names_check : PyQt6 QCheckBox
        check box to indicate if 4FGL catalog names should be used for
        extended sources instead of multi-wavelength association names
    use_old_name_convention_check : PyQt6 QCheckBox
        check box to indicate if the old naming convention from earlier
        #FGL scripts should be used (prepend an underscore and remove
        all spaces from names)
    variable_sources_free_check : PyQt6 QCheckBox
        check box to indicate if variable sources should be have free
        normalization parameters even if they do not satisfy other
        criteria

    Methods
    -------
    choose_catalog()
        handles opening file browser to select the 4FGL catalog file
    choose_galactic_model_file()
        handles opening file browser to select the Galactic diffuse
        emission file
    choose_extended_directory()
        handles opening the directory browser to select the directory
        with the extended source template files
    choose_isotropic_template_file()
        handles opening file browser to select the isotropic diffuse
        component template file
    choose_save_directory()
        handles opening the directory browser to select the target
        directory for the source model
    create_eighth_row()
        create the eighth_row widget and all associated widgets
    create_eleventh_row()
        create the eleventh_row widget and all associated widgets
    create_fifth_row()
        create the fifth_row widget and all associated widgets
    create_first_row()
        create the first_row widget and all associated widgets
    create_fourth_row()
        create the fourth_row widget and all associated widgets
    create_model()
        invoked by the make_model_button, use the supplied inputs
        in a SourceList object from make4FGLxml to create a
        source model XML file, and possibly a ds9-style .reg file
    create_ninth_row()
        create the ninth_row widget and all associated widgets
    create_second_row()
        create the second_row widget and all associated widgets
    create_seventh_row()
        create the seventh_row widget and all associated widgets
    create_sixth_row()
        create the sixth_row widget and all associated widgets
    create_tenth_row()
        create the tenth_row widget and all associated widgets
    create_third_row()
        create the third_row widget and all associated widgets
    create_twelvth_row()
        create the twelvth_row widget and all associated widgets
    create_widgets()
        create all necessary widgets and assign them as object
        attributes
    get_ROI_center()
        handles opening file browser to select LAT event file
        from which the method then extracts the region of
        interest information
    pack_eighth_row()
        pack the eighth_row widget and all associated widgets
    pack_eleventh_row()
        pack the eleventh_row widget and all associated widgets
    pack_fifth_row()
        pack the fifth_row widget and all associated widgets
    pack_first_row()
        pack the first_row widget and all associated widgets
    pack_fourth_row()
        pack the fourth_row widget and all associated widgets
    pack_ninth_row()
        pack the ninth_row widget and all associated widgets
    pack_second_row()
        pack the second_row widget and all associated widgets
    pack_seventh_row()
        pack the seventh_row widget and all associated widgets
    pack_sixth_row()
        pack the sixth_row widget and all associated widgets
    pack_tenth_row()
        pack the tenth_row widget and all associated widgets
    pack_third_row()
        pack the third_row widget and all associated widgets
    pack_twelvth_row()
        pack the twelvth_row widget and all associated widgets
    pack_widgets()
        pack all the widgets in the individual rows and within
        the ControlWidget
    quit_application()
        restore the old stdout and stderror and quit the GUI
    set_make_model_state()
        check the widgets set to be watched to determine if
        the make_model_button is active or not
    set_significance_type()
        determine if the significance threshold is for TS or
        average significance based on the extension of the
        4FGL catalog file chosen
    '''
    
    def __init__(self,main_window):
        '''
        create the ControlWidget object with all child widgets
        and pack everything in a layout of vertically stacked
        rows in which widgets are horizontally stacked
        '''

        #init for the base class
        super().__init__(main_window)

        #pointer to the main window this widget is in
        self.main_window=main_window

        #get the FERMI_DIR variable if is set
        self.fermi_dir=os.environ.get('FERMI_DIR')

        #create all the widgets
        self.create_widgets()

    def create_widgets(self):
        '''
        create all the widgets for each row of the GUI
        '''

        #create the 'row' widgets and all associated widgets
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
        '''
        create the first_row widget and all associated widgete=s
        '''
        
        self.first_row=QWidget(parent=self)
        
        #create widgets for inputing catalog file
        self.catalog_button=QPushButton("Catalog",parent=self.first_row)
        self.catalog_button.clicked.connect(self.choose_catalog)
        self.catalog_button.setFixedWidth(100)
        
        self.catalog_entry=QLineEdit(parent=self.first_row)
        self.catalog_entry.setFixedWidth(400)

        #set a trace on the catalog_entry to detect
        #the file type and change the significance label
        #between average sigma and test statistic
        self.catalog_entry.textEdited.connect(self.set_significance_type)

        #create widgets for the DR version
        self.DR_label=QLabel('DR:',parent=self.first_row)
        self.DR_spinbox=QSpinBox(parent=self.first_row)
        self.DR_spinbox.setValue(3)
        self.DR_spinbox.setFixedWidth(30)

        #set reasonable limits on the spinbox
        self.DR_spinbox.setMinimum(1)
        self.DR_spinbox.setMaximum(4)

    def create_second_row(self):
        '''
        create the second_row widget and all associated widgete=s
        '''
        
        self.second_row=QWidget(parent=self)

        #create widgets relating to output file
        self.save_directory_button=QPushButton("Save Directory",parent=self.second_row)
        self.save_directory_button.clicked.connect(self.choose_save_directory)
        self.save_directory_button.setFixedWidth(100)
        
        self.save_directory_entry=QLineEdit(os.getcwd(),parent=self.second_row)
        self.save_directory_entry.setFixedWidth(400)

        self.output_file_label=QLabel('Output File:',parent=self.second_row)
        self.output_file_entry=QLineEdit('my_model.xml',parent=self.second_row)

    def create_third_row(self):
        '''
        create the third_row widget and all associated widgete=s
        '''
        
        self.third_row=QWidget(parent=self)
        
        #create widgets for the region file
        self.region_check=QCheckBox("Make ds9 region file?",parent=self.third_row)
        self.region_entry=QLineEdit('my_region.reg',parent=self.third_row)
        ##set up to tie to default above but do we need to do this?
        ##do we want to have this always match the output_file text? might annoy
        ##users who edit this and then make a change to the xml name

    def create_fourth_row(self):
        '''
        create the fourth_row widget and all associated widgete=s
        '''
        
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
        self.ROI_radius_entry.setFixedWidth(35)
        
        #now need to monitor for when those values are changed
        #either by the user or programmatically
        self.ROI_center_RA_entry.textEdited.connect(self.set_make_model_state)
        self.ROI_center_RA_entry.textChanged.connect(self.set_make_model_state)

        self.ROI_center_DEC_entry.textEdited.connect(self.set_make_model_state)
        self.ROI_center_DEC_entry.textChanged.connect(self.set_make_model_state)

        self.ROI_radius_entry.textEdited.connect(self.set_make_model_state)
        self.ROI_radius_entry.textChanged.connect(self.set_make_model_state)

        #allow user to instead choose an event file to get the ROI center information
        self.or_label=QLabel(' or ',parent=self.fourth_row)
        self.use_event_file_button=QPushButton("Get From Event File",parent=self.fourth_row)
        self.use_event_file_button.clicked.connect(self.get_ROI_center)
        self.use_event_file_button.setFixedWidth(125)
        
        self.event_file_label=QLabel("Event File:",parent=self.fourth_row)
        self.event_file_entry=QLineEdit(parent=self.fourth_row)

    def create_fifth_row(self):
        '''
        create the fifth_row widget and all associated widgete=s
        '''
        
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

        #put checks on these values being edited to see if the
        #make model button can be active
        self.free_radius_entry.textEdited.connect(self.set_make_model_state)
        self.max_free_radius_entry.textEdited.connect(self.set_make_model_state)
        self.extra_radius_entry.textEdited.connect(self.set_make_model_state)

    def create_sixth_row(self):
        '''
        create the sixth_row widget and all associated widgete=s
        '''
        
        self.sixth_row=QWidget(parent=self)
        
        #would be good to have a 'trace' set on this to modify based on XML or FITS
        #version of the catalog (e.g., avg sig vs avg TS)
        self.significance_label=QLabel("Minimum Average Significance In Catalog For Free Sources:",parent=self.sixth_row)
        self.significance_validator=QDoubleValidator(bottom=0)
        self.significance_entry=QLineEdit(parent=self.sixth_row)
        self.significance_entry.setValidator(self.significance_validator)
        self.significance_entry.setFixedWidth(25)

        #check if the make_model_button can be active based on this value
        self.significance_entry.textEdited.connect(self.set_make_model_state)

    def create_seventh_row(self):
        '''
        create the seventh_row widget and all associated widgete=s
        '''
        
        self.seventh_row=QWidget(parent=self)
        
        #now for option check boxes
        self.normalizations_only_check=QCheckBox("Only Free Normalization Parameters?",parent=self.seventh_row)
        self.variable_sources_free_check=QCheckBox("Free Variable Sources?",parent=self.seventh_row)
        self.force_point_source_check=QCheckBox("Change Extended Sources To Point Sources?",parent=self.seventh_row)

        self.variable_sources_free_check.setChecked(True)

    def create_eighth_row(self):
        '''
        create the eighth_row widget and all associated widgete=s
        '''
        
        self.eighth_row=QWidget(parent=self)
        
        self.use_extended_catalog_names_check=QCheckBox("Use Catalog Names For Extended Sources?",parent=self.eighth_row)
        self.free_galactic_index_check=QCheckBox("Modify Galactic Spectrum By Power Law?",parent=self.eighth_row)
        self.use_old_name_convention_check=QCheckBox("Use '_4FGLJXXXX.X+XXXX' Naming Convention?",parent=self.eighth_row)

        self.free_galactic_index_check.setChecked(True)

    def create_ninth_row(self):
        '''
        create the ninth_row widget and all associated widgete=s
        '''
        
        self.ninth_row=QWidget(parent=self)
        
        #now for diffuse files/directories and names
        if self.fermi_dir is None:
            self.galactic_model_file_entry=QLineEdit(os.sep.join(['$(FERMI_DIR)','refdata','fermi',
                    'galdiffuse','gll_iem_v07.fits']),parent=self.ninth_row)
        
        else:
            self.galactic_model_file_entry=QLineEdit(os.path.join(self.fermi_dir,'refdata','fermi',
                    'galdiffuse','gll_iem_v07.fits'),parent=self.ninth_row)

        self.galactic_model_file_entry.setFixedWidth(400)

        self.galactic_model_file_label=QLabel('File Path For Galactic Diffuse Model:',parent=self.ninth_row)
        self.galactic_model_file_button=QPushButton("Change File",parent=self.ninth_row)
        self.galactic_model_file_button.clicked.connect(self.choose_galactic_model_file)
        self.galactic_model_file_button.setFixedWidth(100)

        self.galactic_name_label=QLabel('Galactic Diffuse Name:',parent=self.ninth_row)
        self.galactic_name_entry=QLineEdit('gll_iem_v07',parent=self.ninth_row)

        #link these entries to make_model_button state as well
        #essentially just making sure they're not empty
        self.galactic_model_file_entry.textEdited.connect(self.set_make_model_state)
        self.galactic_model_file_entry.textChanged.connect(self.set_make_model_state)
        
        self.galactic_name_entry.textEdited.connect(self.set_make_model_state)

    def create_tenth_row(self):
        '''
        create the tenth_row widget and all associated widgete=s
        '''
        
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
        self.isotropic_template_file_button.clicked.connect(self.choose_isotropic_template_file)
        self.isotropic_template_file_button.setFixedWidth(100)

        self.isotropic_name_label=QLabel('Isotropic Component Name:',parent=self.tenth_row)
        self.isotropic_name_entry=QLineEdit('iso_P8R3_SOURCE_V3_v1',parent=self.tenth_row)

        #use edits here to check state of make_model_button
        self.isotropic_template_file_entry.textEdited.connect(self.set_make_model_state)
        self.isotropic_template_file_entry.textChanged.connect(self.set_make_model_state)

        self.isotropic_name_entry.textEdited.connect(self.set_make_model_state)

    def create_eleventh_row(self):
        '''
        create the eleventh_row widget and all associated widgete=s
        '''
        
        self.eleventh_row=QWidget(parent=self)
        
        #now for diffuse files/directories and names
        if self.fermi_dir is None:
            self.extended_template_directory_entry=\
                    QLineEdit(os.sep.join(['$(FERMI_DIR)','data','pyBurstAnalysisGUI','templates']),
                              parent=self.eleventh_row)
        
        else:
            self.extended_template_directory_entry=\
                    QLineEdit(os.path.join(self.fermi_dir,'data','pyBurstAnalysisGUI','template'),
                             parent=self.eleventh_row)

        self.extended_template_directory_entry.setFixedWidth(400)

        self.extended_template_directory_label=QLabel('Directory With Extended Source Templates:',
                                                      parent=self.eleventh_row)
        self.extended_template_directory_button=QPushButton("Change Directory",
                                                            parent=self.eleventh_row)
        self.extended_template_directory_button.clicked.connect(self.choose_extended_directory)
        self.extended_template_directory_button.setFixedWidth(125)

    def create_twelvth_row(self):
        '''
        create the twelvth_row widget and all associated widgete=s
        '''
        
        self.twelvth_row=QWidget(parent=self)
        
        #now we need a button to actually make the files
        self.make_model_button=QPushButton('Make Model',parent=self.twelvth_row)
        self.make_model_button.clicked.connect(self.create_model)

        #we will start off with the button disabled
        #need to have valid inputs before we can push the button
        self.make_model_button.setEnabled(False)
        
        #and finish it off with a 'Quit' button
        self.quit_button=QPushButton("Quit",parent=self.twelvth_row)
        self.quit_button.clicked.connect(self.quit_application)

        self.make_model_button.setFixedWidth(100)
        self.quit_button.setFixedWidth(100)

    
    def pack_widgets(self):
        '''
        pack all of the widgets with each 'row' widget and
        within the ControlWidget
        '''
        
        #We're going to stack rows vertically, so make a main layout
        #and then pack the widgets in 'rows', with only related widgets in each row
        #so each row will need to be a horizontal layout
        self.main_layout=QVBoxLayout()
        #self.main_layout.setSpacing(0)

        #peack all widgets within each 'row'
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

        #add all 'row' widgets to the ControlWidget Layout
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
        '''
        pack the first_row widget and all associated widgete=s
        '''
        
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
        '''
        pack the second_row widget and all associated widgete=s
        '''
        
        self.second_row_layout=QHBoxLayout()
        self.second_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the second row, save directory and output fil stuff
        self.second_row_layout.addWidget(self.save_directory_button)
        self.second_row_layout.addWidget(self.save_directory_entry)
        self.second_row_layout.addWidget(self.output_file_label)
        self.second_row_layout.addWidget(self.output_file_entry)

        self.second_row.setLayout(self.second_row_layout)
        
    def pack_third_row(self):
        '''
        pack the third_row widget and all associated widgete=s
        '''
        
        self.third_row_layout=QHBoxLayout()
        self.third_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the third row, region file stuff
        self.third_row_layout.addWidget(self.region_check)
        self.third_row_layout.addWidget(self.region_entry)
        self.third_row_layout.addStretch(1)

        self.third_row.setLayout(self.third_row_layout)
        
    def pack_fourth_row(self):
        '''
        pack the fourth_row widget and all associated widgete=s
        '''
        
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
        '''
        pack the fifth_row widget and all associated widgete=s
        '''
        
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
        '''
        pack the sixth_row widget and all associated widgete=s
        '''
        
        self.sixth_row_layout=QHBoxLayout()
        self.sixth_row_layout.setContentsMargins(0,0,0,0)

        #grid the sixth row, minimum significance stuff
        self.sixth_row_layout.addWidget(self.significance_label)
        self.sixth_row_layout.addWidget(self.significance_entry)
        self.sixth_row_layout.addStretch(1)

        self.sixth_row.setLayout(self.sixth_row_layout)

    def pack_seventh_row(self):
        '''
        pack the seventh_row widget and all associated widgete=s
        '''
        
        self.seventh_row_layout=QHBoxLayout()
        self.seventh_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the seventh row, normalization and variability checkboxes
        self.seventh_row_layout.addWidget(self.normalizations_only_check)
        self.seventh_row_layout.addWidget(self.variable_sources_free_check)
        self.seventh_row_layout.addWidget(self.force_point_source_check)
        self.seventh_row_layout.addStretch(1)

        self.seventh_row.setLayout(self.seventh_row_layout)

    def pack_eighth_row(self):
        '''
        pack the eighth_row widget and all associated widgete=s
        '''
        
        self.eighth_row_layout=QHBoxLayout()
        self.eighth_row_layout.setContentsMargins(0,0,0,0)

        #grid the eighth row, extended source related checkboxes
        self.eighth_row_layout.addWidget(self.use_extended_catalog_names_check)
        self.eighth_row_layout.addWidget(self.free_galactic_index_check)
        self.eighth_row_layout.addWidget(self.use_old_name_convention_check)

        self.eighth_row.setLayout(self.eighth_row_layout)

    def pack_ninth_row(self):
        '''
        pack the ninth_row widget and all associated widgete=s
        '''
        
        self.ninth_row_layout=QHBoxLayout()
        self.ninth_row_layout.setContentsMargins(0,0,0,0)

        #grid the ninth row, Galactic diffuse model stuff
        self.ninth_row_layout.addWidget(self.galactic_model_file_label)
        self.ninth_row_layout.addWidget(self.galactic_model_file_entry)
        self.ninth_row_layout.addWidget(self.galactic_model_file_button)
        self.ninth_row_layout.addStretch(1)
        self.ninth_row_layout.addWidget(self.galactic_name_label)
        self.ninth_row_layout.addWidget(self.galactic_name_entry)

        self.ninth_row.setLayout(self.ninth_row_layout)

    def pack_tenth_row(self):
        '''
        pack the tenth_row widget and all associated widgete=s
        '''
        
        self.tenth_row_layout=QHBoxLayout()
        self.tenth_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the tenth row, Isotropic diffuse stuff
        self.tenth_row_layout.addWidget(self.isotropic_template_file_label)
        self.tenth_row_layout.addWidget(self.isotropic_template_file_entry)
        self.tenth_row_layout.addWidget(self.isotropic_template_file_button)
        self.tenth_row_layout.addStretch(1)
        self.tenth_row_layout.addWidget(self.isotropic_name_label)
        self.tenth_row_layout.addWidget(self.isotropic_name_entry)

        self.tenth_row.setLayout(self.tenth_row_layout)

    def pack_eleventh_row(self):
        '''
        pack the eleventh_row widget and all associated widgete=s
        '''
        
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
        '''
        pack the twelvth_row widget and all associated widgete=s
        '''
        
        self.twelvth_row_layout=QHBoxLayout()
        self.twelvth_row_layout.setContentsMargins(0,0,0,0)
        
        #grid the twelvth row, buttons to generate the file(s) and to quit
        self.twelvth_row_layout.addWidget(self.make_model_button)
        self.twelvth_row_layout.addStretch(1)
        self.twelvth_row_layout.addWidget(self.quit_button)

        self.twelvth_row.setLayout(self.twelvth_row_layout)

    def set_make_model_state(self):
        '''
        determine if the make_model_button should be enabled or not
        '''
        
        #lots of things to check so that the make_model_button can be enabled
        if os.path.exists(self.catalog_entry.text()) and\
           self.ROI_center_RA_entry.hasAcceptableInput() and\
           self.ROI_center_DEC_entry.hasAcceptableInput() and\
           self.ROI_radius_entry.hasAcceptableInput() and\
           self.free_radius_entry.hasAcceptableInput() and\
           self.max_free_radius_entry.hasAcceptableInput() and\
           self.extra_radius_entry.hasAcceptableInput() and\
           self.significance_entry.hasAcceptableInput() and\
           self.galactic_model_file_entry.text()!='' and\
           self.galactic_name_entry.text()!='' and\
           self.isotropic_template_file_entry.text()!='' and\
           self.isotropic_name_entry.text()!='':
            self.make_model_button.setEnabled(True)

        else:
            self.make_model_button.setEnabled(False)

    def choose_catalog(self):
        '''
        open a QFileDialog for selecting the 4FGL catalog file and then
        assign that to the appropriate entry
        '''
        
        #first, determine the starting directory
        start_dir=(os.getcwd() if self.fermi_dir is None else self.fermi_dir)

        #open the file dialog to select the catalog file
        catalog_file_name=QFileDialog.getOpenFileName(self,'Choose Catalog File',start_dir,
                        "Fermi LAT catlog files (gll_psc*);;FITS files (*.fit*);;XML files (*.xml);;All Files(*)")[0]

        #set the returned value to the corresponding entry
        self.catalog_entry.setText(catalog_file_name)
        return

    def set_significance_type(self):
        '''
        change the text of the significance_label to reflect what value
        the significance limit applies to, TS or average significance, based
        on the extension of the 4FGL catalog file
        '''
        
        ##if os.path.dirname(self.catalog_entry.text()).split('.')[-1].lower()=='xml':
        if self.catalog_entry.text().lower().endswith('xml'):
            self.significance_label.setText("Minimum Test Statistic Value In Catalog For Free Sources:")
        else:
            self.significance_label.setText("Minimum Average Significance In Catalog For Free Sources:")

        #also check the make_model_button
        self.set_make_model_state()

    def choose_save_directory(self):
        '''
        open a QFileDialog to allow the user to select the target directory for
        the source model file and set the selection to the appropriate entry
        '''
        
        #first, determine the starting directory
        start_dir=(self.save_directory_entry.text()\
                   if os.path.exists(self.save_directory_entry.text()) else os.getcwd())
        
        #open file dialog to select a directory
        save_directory=QFileDialog.getExistingDirectory(self,'Choose Save Directory',start_dir)

        #set the returned value to the corresponding entry
        self.save_directory_entry.setText(save_directory)
        return

    def choose_extended_directory(self):
        '''
        open a QFileDialog to allow the user to select the directory with the
        extended source template files and then set the selection to
        the appropriate entry
        '''
        
        #first, determine the starting directory for the file dialog
        start_dir=(self.extended_template_directory_entry.text()\
                   if os.path.exists(self.extended_template_directory_entry.text()) else os.getcwd())

        #open up the dialog
        extended_directory=QFileDialog.getExistingDirectory(self,
                                'Choose Extended Source Templates Directory',start_dir)

        #set this as the corresponding entry text
        self.extended_template_directory_entry.setText(extended_directory)
        return

    def choose_galactic_model_file(self):
        '''
        open a QFileDialog to allow the user to select the Galactic diffuse model file
        and then set the selection to the appropriate entry
        '''
        
        #first, determine the starting directory for the file dialog
        start_dir=(self.galactic_model_file_entry.text()\
                   if os.path.exists(self.galactic_model_file_entry.text()) else os.getcwd())

        #open up the dialog
        galactic_file=QFileDialog.getOpenFileName(self,'Choose Galactic Diffuse Model File',
                                        start_dir,'FITS files (*.fit*);;All Files (*)')[0]

        #set the returned value as the corresponding entry text
        if galactic_file!='' and galactic_file is not None:
            self.galactic_model_file_entry.setText(galactic_file)
        return

    def choose_isotropic_template_file(self):
        '''
        open a QFileDialog to allow the user to select the isotropic diffuse component
        template file and set the selection to the appropriate entry
        '''
        
        #first, determine the starting directory
        start_dir=(self.isotropic_template_file_entry.text()\
                   if os.path.exists(self.isotropic_template_file_entry.text()) else os.getcwd())

        #open up the dialot
        isotropic_file=QFileDialog.getOpenFileName(self,'Choose Isotropic Component Template File',
                            start_dir,'Text files (*.txt);;All Files (*)')[0]

        #set the returned value as the corresponding entry text
        if isotropic_file!='' and isotropic_file is not None:
            self.isotropic_template_file_entry.setText(isotropic_file)

    def get_ROI_center(self):
        '''
        open a QFileDialog to allow the user to select a LAT event file from which the
        region of interest information is then extracted
        '''
        
        #first, we need to determine if the current value in the event file entry is valid
        if not os.path.exists(self.event_file_entry.text()):
            #if we get in here, open a file dialog to choose an event file, start at the current dir
            event_file=QFileDialog.getOpenFileName(self,'Choose Fermi LAT Event File',os.getcwd(),
                            "FITS files (*.fit*);;All Files (*)")[0]
            self.event_file_entry.setText(event_file)
        else:
            event_file=self.event_file_entry.text()

        #check that the event entry is not empty
        if event_file!='' and event_file is not None:
            #now we feed this to one of the make4FGLxml utility functions
            RA,DEC,radius=get_ROI_from_event_file(event_file)
    
            #finish by setting the corresponding entry
            self.ROI_center_RA_entry.setText(f'{RA:.3f}')
            self.ROI_center_DEC_entry.setText(f'{DEC:.3f}')
            self.ROI_radius_entry.setText(f'{radius:.1f}')

        else:
            warnings.warn('Note, must select a valid event file')

    def create_model(self):
        '''
        use the inputs specified in the GUI to create a SourceList object
        and creat a source model XML file
        '''
        
        #do stuff to make the XML file (and possible .reg file)
        
        source_list=SourceList(self.catalog_entry.text(),
                               [float(self.ROI_center_RA_entry.text()),
                                  float(self.ROI_center_DEC_entry.text()),
                                  float(self.ROI_radius_entry.text())],
                               self.output_file_entry.text(),
                               self.DR_spinbox.value(),
                               self.save_directory_entry.text())

        #print to verify
        source_list.Print()

        #now make the model and we're done
        source_list.make_model(self.galactic_model_file_entry.text(),
                               self.galactic_name_entry.text(),
                               self.isotropic_template_file_entry.text(),
                               self.isotropic_name_entry.text(),
                               self.normalizations_only_check.isChecked(),
                               self.extended_template_directory_entry.text(),
                               float(self.free_radius_entry.text()),
                               float(self.max_free_radius_entry.text()),
                               float(self.extra_radius_entry.text()),
                               float(self.significance_entry.text()),
                               self.variable_sources_free_check.isChecked(),
                               self.force_point_source_check.isChecked(),
                               self.use_extended_catalog_names_check.isChecked(),
                               self.region_check.isChecked(),
                               self.region_entry.text(),
                               self.free_galactic_index_check.isChecked(),
                               self.use_old_name_convention_check.isChecked())

    def quit_application(self):
        '''
        print a polite message to the terminal before resetting
        stdout and stderr and then quitting the GUI
        '''
        
        self.main_window.terminal_widget.write('Goodbye!')
        
        sys.stdout=self.main_window.old_stdout
        sys.stderr=self.main_window.old_stderr

        QTimer.singleShot(1000,QApplication.instance().quit)

