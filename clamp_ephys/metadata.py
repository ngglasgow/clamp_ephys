from collections import OrderedDict
import pandas as pd
import os


def get_metadata(file, data_path):
    '''Takes a filename and parses it for metadata, and returns metadata in an
    orderedDict as a pandas DataFrame for saving later
    Also takes information from the cell spreadsheet in data_notes'''

    # pull out cell id, cell number, date and condition
    file_split = file.split('_')
    cell_id = file_split[0]+'_'+file_split[1]
    cell_num = cell_id[-1:]
    date = '20'+cell_id[2:4]+'-'+cell_id[4:6]+'-'+cell_id[6:8]

    if 'light' in file:
        condition = 'light'
    else:
        condition = 'spontaneous'

    # grab metadata from data notes spreadsheet
    data_notes = pd.read_csv(data_path, index_col=[0])
    file_data = data_notes[data_notes['File Path'] == file]
    cell_path = file_data['File Path'].tolist()[0]
    cell_type = file_data['Cell Type'].tolist()[0]

    # save metadate into orderedDict pandas DataFrame
    dict = OrderedDict()
    dict['Date'] = date
    dict['Cell ID'] = cell_id
    dict['Cell Number'] = cell_num
    dict['Cell Path'] = cell_path
    dict['Condition'] = condition
    dict['Cell Type'] = cell_type
    metadata = pd.DataFrame(dict, index=range(1))

    return metadata


def create_data_notes(timepoint, summary_file, ibw_file_paths):
    '''
    Create data_notes summary spreadsheets
    
    Parameters
    ----------
    timepoint: str
        Name of the injection timepoint used in the analysis
    summary_file: .csv file
        Manually-generated summary file of all cells in the dataset of a timepoint
    ibw_file_paths:  list
        List of ibw files found for a timepoint

    Returns:
    file_name_list: list
        List of all the file names in the timepoint data set
    data_notes: pandas.DataFrame
        DataFrame of parsed notes for each cell from manually-inputted summary_file
    '''
    # Open the notes spreadsheet and parse for what we want to analyze ## '''
    # open metadata file
    data_notes = pd.read_csv(os.path.join(paths.tables, summary_file))

    # pull out cell_id for directory, file name, and make the full path
    file_name_list = data_notes['Cell name'].tolist()

    data_notes = pd.concat([pd.DataFrame({'File Path': ibw_file_paths}), 
        data_notes], axis=1)

    light_data_notes = data_notes[data_notes['Cell name'].str.contains("light")]
    spontaneous_data_notes = data_notes[data_notes['Cell name'].str.contains("spontaneous")]

    data_notes.to_csv(os.path.join(paths.tables, '{}_data_notes.csv'.format(timepoint)))
    light_data_notes.to_csv(os.path.join(paths.tables, '{}_light_data_notes.csv'.format(timepoint)))
    spontaneous_data_notes.to_csv(os.path.join(paths.tables, '{}_spontaneous_data_notes.csv'.format(timepoint)))

    return file_name_list, data_notes


class file_structure:
    def __init__(self, location, project_path):
        '''
        Creates an object with paths as attributes:
        location:   str value 'local' or 'server' only, refers to where you are
                    doing the actual work, 'local' by default.
        project_path:   str of the root project path from wherever your home dir is
        '''
        machine = platform.uname()[0]

        if location == 'local':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"C:\Users\jhuang\Documents\phd_projects"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        elif location == 'server':
            if machine == 'Darwin':
                home_dir = '/Volumes/Urban'

            elif machine == 'Linux':
                home_dir = os.path.join(os.path.expanduser('~'), 'urban/neurobio/Huang')

            elif machine == 'Windows':
                home_dir = r"N:\Huang"

            else:
                print("OS not recognized. \nPlease see Nate for correction.")

        self.project = os.path.join(home_dir, project_path)
        self.figures = os.path.join(self.project, 'figures')
        self.tables = os.path.join(self.project, 'tables')
        self.p2 = os.path.join(self.project, 'data/p2')
        self.p2_files = os.listdir(self.p2)
        self.p14 = os.path.join(self.project, 'data/p14')
        self.p14_files = os.listdir(self.p14)

    def __repr__(self):
        return 'Project file structure and file lists for {}'.format(self.project)
