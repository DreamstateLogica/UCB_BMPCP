#!/usr/bin/env python3
"""
A script to quickly pull desired data from ORCA .out files within a series of
subdirectories and compile the results into a single CSV file.

The script searches for subdirectories starting with 'mol_' in the execution
directory, finds the first .out file in each, and processes it.
"""
__author__ = "Peter Waddell"
__copyright__ = "Copyright 2024"
__credits__ = ["Peter Waddell"]
__version__ = "0.2.0"
__date__ = "2024/07/11"
__maintainer__ = "Peter Waddell"
__email__ = "pmwaddell9@gmail.com"
__status__ = "Prototype"

import os
import sys
import pandas as pd

# --- FIX FOR ModuleNotFoundError ---
# This block correctly identifies the project's root directory and adds it to
# the system path, allowing absolute imports to work correctly.
try:
    src_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(src_dir))
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
except NameError:
    # Fallback for environments where __file__ is not defined
    sys.path.insert(0, os.path.abspath('../..'))
# --- END FIX ---

from orca_data_extraction.src.structure_data_builder import StructureDataBuilder
from orca_data_extraction.src.orca_out_to_json import make_json_list

def create_csv_from_sds(sd_list, csv_name):
    """
    Writes the data in a list of StructureData instances to a CSV file.

    Parameters
    ----------
    sd_list : list
        List containing the set of StructureData instances that each come
        from the ORCA .out files.
    csv_name : str
        Name of the CSV file where the data will be stored.
    """
    if not sd_list:
        print("No data was successfully extracted. CSV file will not be created.")
        return

    json_list = make_json_list(sd_list)
    df = pd.json_normalize(json_list)

    # Make the 'origin' column the first column for clarity
    if 'origin' in df.columns:
        cols = ['origin'] + [col for col in df.columns if col != 'origin']
        df = df[cols]

    df.to_csv(csv_name + '.csv', index=False)
    print(f'\nProcess complete! Results saved as "{csv_name}.csv"')


def main():
    # --- MODIFIED ARGUMENT PARSING for directory parsing ---
    if len(sys.argv) != 3:
        print("Usage: python orca_out_to_csv.py <input_json_file> <output_csv_name>")
        print("Example: python orca_out_to_csv.py gen_input.json my_results")
        sys.exit(1)

    inputs_name = sys.argv[1]
    csv_name = sys.argv[2]
    # --- END MODIFICATION ---

    # The input JSON file is expected to be in the execution directory
    execution_path = os.getcwd()
    json_input_path = os.path.join(execution_path, inputs_name)

    if not os.path.isfile(json_input_path):
        print(f"Error: Input JSON file '{json_input_path}' not found.")
        quit()

    print(f"Using input file: {json_input_path}\n")
    
    all_sd_objects = []
    structure_data_builder = StructureDataBuilder(json_input_path)
    
    # Find all subdirectories starting with 'mol_'
    subdirs = sorted([d for d in os.listdir(execution_path) if os.path.isdir(os.path.join(execution_path, d)) and d.startswith('mol_')])

    if not subdirs:
        print("No subdirectories with prefix 'mol_' found in the current directory.")
        return

    for dirname in subdirs:
        dirpath = os.path.join(execution_path, dirname)
        print(f"--- Searching in directory: {dirpath} ---")
        
        found_out_file = False
        for filename in os.listdir(dirpath):
            if filename.endswith(".out"):
                out_filepath = os.path.join(dirpath, filename)
                print(f"Found output file: {filename}")
                
                try:
                    # Build the StructureData object and add the origin directory name
                    sd_object = structure_data_builder.build(out_filepath)
                    sd_object.origin = dirname  # Add origin attribute to the object
                    all_sd_objects.append(sd_object)
                    
                except Exception as e:
                    print(f"  -> An error occurred while processing {filename}: {e}\n")
                
                found_out_file = True
                break  # Process only the first .out file found in the directory
        
        if not found_out_file:
            print(f"Warning: No .out file found in {dirname}")
        print("-" * (len(dirpath) + 24))


    create_csv_from_sds(all_sd_objects, csv_name)


if __name__ == '__main__':
    main()
