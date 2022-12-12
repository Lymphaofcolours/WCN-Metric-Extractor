import WCN
import glob
import os
import numpy as np
import warnings
import statistics
import pandas as pd
from collections import defaultdict
import re
import json


def input_file_path():
    '''
    :return: str: Input filepath. Detect invalid inputs and ask for path again.
    '''
    try:
        filepath: str = input('Please, Insert Your Input Directory') or '/home/lympha/RAS/isoraspdb/'
        print('This process is computationally intensive. Make yourself a coffee. A big one.')
        return str(filepath)
    except ValueError:
        print('Wrong Input. Please, Try Again:')


def output_file_path():
    '''
    :return: str: Input filepath. Detect invalid inputs and ask for path again.
    '''
    try:
        filepath: str = input('Please, Insert Your Output Directory') or '/home/lympha/RAS/csv/'
        return str(filepath)
    except ValueError:
        print('Wrong Input. Please, Try Again:')


class Extractor:
    def name_extractor(self, pdb): -> str
        '''
        Extract name of protein from file name.
        '''
        pdb_file_name: str = pdb.split('/')[-1]
        pdb_protein: str = pdb_file_name.split('.')[0]
        return pdb_protein

    def classification_extractor(self, pdb): -> str
        '''
        Generate string containing 'active' or 'inactive' state from file name.
        '''
        pdb_file_name: str = pdb.split('/')[-1]
        pdb_protein: str = pdb_file_name.split('.')[0]
        if '_a' in pdb_protein:
            return 'Active'
        elif '_i' in pdb_protein:
            return 'Inactive'


class Calculator:
    '''
    Calculate means for current sequence
    '''
    def arithmetic_mean(self, array: np.ndarray) -> float:
        mean = np.mean(array)
        return mean

    def geometric_mean(self, array: np.ndarray) -> float:
        array = np.absolute(array)
        mean = statistics.geometric_mean(array)
        return mean

    def harmonic_mean(self, array: np.ndarray) -> float:
        array = np.absolute(array)
        mean = statistics.harmonic_mean(array)
        return mean


def emptydict(): -> list or ndarray
        '''
        Return value for use in defaultdict for later dataframe creation. Can be easily modified to return empty ndarrays for array calculation,
        but as the heavy workload is carried already in numpy arrays generated through the WCN module, it´s currently working on lists of n float values,
        where n is the number of files in the input directory.
        '''
        # return np.empty([len(glob.glob(os.path.join(input_file_path(), '*.pdb'))), 1], dtype=float)
        return []


class Listdb:
    pdb_names: list = []
    pdb_classification: list = []
    
    # Available metrics in the WCN library. Makes it very easy to add new study variables, or apply them to specific domains of the peptide
    available_metrics: list = ['getBfactors',
                     'calculateCAlpha',
                     'calculateAllAtom',
                     'calculateAveragePerResidue'
                         ]
    
    # Available methods in the Calculator class, dunder methods filtered out. Makes it very easy to add new operations to the Calculator class
    available_ops: list = [method for method in dir(Calculator) if method.startswith('__') is False]
    
    # Default dictionary, generated through the later constructor. Defaults to empty list or empty numpy array for memory preallocation
    ras_data = defaultdict(emptydict)


if __name__ == '__main__':
    # .pdb files are usually full of discontinuities in their protein chains, which raise warnings. As there is
    # currently no known method to block them specifically, warnings are filtered off the console
    warnings.filterwarnings('ignore')
    
    # Loops over the files in the input directory. Later versions of this script could open them all, extract available metrics and operate directly over
    # constructor-generated dataframes
    for pdb in glob.glob(os.path.join(input_file_path(), '*.pdb')):
        #Instantiates an object for each object in the input directory
        loc = WCN.WCNObject(pdb)
        
        
        # Exports name and classification to list, so they can be added as a column once the dataframe is created
        Listdb.pdb_names.append(Extractor.name_extractor(loc, pdb))
        Listdb.pdb_names.append(Extractor.classification_extractor(loc, pdb))
        print(f'Processing {pdb.split("/")[-1].split(".")[0]}. Please Wait...')
        
        
        for metric in Listdb.available_metrics:
            # Obtains metric from the list of available metrics and call to the attribute of the WCN object from it
            tmpmetric = getattr(WCN.WCNObject, metric)
            # Extracts name of the metric based on capital letter ocurrence in the method name string
            tmpmetric_name = ''.join(re.findall('[A-Z][^A-Z]*', metric))
            
            
            # There is high recursivity in this code, so it features a quadratic time complexity or even more. Future versions of this script
            # could extract raw features in dataframes (reducing complexity in a factor of n) and normalize or apply functions through applymaps.
            for op in Listdb.available_ops:
                # Calls to the available operations from the list of methods in the Calculator class
                tmpop = getattr(Calculator, op)
                Listdb.ras_data[f'{tmpmetric_name}_{op}'].append(tmpop(tmpmetric(loc), tmpmetric(loc)))


    print(Listdb.ras_data)
    ras_df: pd.DataFrame = pd.DataFrame(Listdb.ras_data)
    # Pickles preserve data types for later analysis
    ras_df.to_pickle(f'{output_file_path()}/rasoop.pkl')
    print(ras_df)

    # JSONs allow us to understand the hierarchies of data in the generated defauldict, in case we have nested arrays
    json_object = json.dump(ras_df, f'{output_file_path()}result.json')

