#!/usr/bin/env python3
"""
A DataSection subclass that extracts thermochemical data (Gibbs Free Energy,
Enthalpy, Entropy, ZPE, and Electronic Energy) from an ORCA output file.
"""
__author__ = "Peter Waddell"
__copyright__ = "Copyright 2024"
__credits__ = ["Peter Waddell"]
__version__ = "0.1.3"
__date__ = "2024/07/10"
__maintainer__ = "Peter Waddell"
__email__ = "pmwaddell9@gmail.com"
__status__ = "Prototype"

import re
from data_section import DataSection


class ThermoData(DataSection):
    """
    Extracts key thermochemical values from an ORCA output file.

    This class searches for the GIBBS FREE ENERGY block and works backwards
    to efficiently find all related thermochemical data.

    Attributes
    ----------
    _section_name : str
        The name of the section to be extracted.

    Methods
    -------
    _find_data(self)
        Extracts the data from the ORCA output file.
    """
    def __init__(self, out_filename, outfile_contents):
        """
        Parameters
        ----------
        out_filename : str
            The name of the ORCA output file.
        outfile_contents : str
            The full text of the ORCA output file.
        """
        super().__init__(out_filename, outfile_contents)
        self._section_name = 'Thermochemistry'

    def _find_data(self):
        """
        Extracts the thermochemical values by searching for each one
        individually for robustness.

        Returns
        -------
        dict
            A dictionary containing all the extracted thermochemical data.
            Returns 'N/A' for values that are not found and prints detailed
            error messages.
        """
        reversed_contents = self._outfile_contents[::-1]
        data = {}
        errors = []

        # Define regex patterns for each thermochemical value (in reverse)
        patterns = {
            'Final Gibbs Free Energy': re.compile(r"hE\s+([\d\.-]+)\s+\.\.\.\s+ygnere\s+eerf\s+sbbG\s+laniF"),
            'Final Entropy Term': re.compile(r"hE\s+([\d\.-]+)\s+\.\.\.\s+mret\s+yportne\s+laniF"),
            'Total Enthalpy': re.compile(r"hE\s+([\d\.-]+)\s+\.\.\.\s+yplahtnE\s+latoT"),
            'Zero Point Energy': re.compile(r"hE\s+([\d\.-]+)\s+\.\.\.\s+ygrene\s+tniop\s+oreZ"),
            'Electronic Energy': re.compile(r"hE\s+([\d\.-]+)\s+\.\.\.\s+ygrene\s+cinortcelE")
        }

        found_any = False
        for key, pattern in patterns.items():
            match = pattern.search(reversed_contents)
            if match:
                # Reverse the captured value back to its correct orientation
                data[key] = match.group(1)[::-1]
                found_any = True
            else:
                data[key] = 'N/A'
                errors.append(f"Could not find '{key}'")
        
        if not found_any:
             print(f"Error in {self._out_filename}: Could not find any thermochemistry data blocks.")
        elif errors:
            print(f"Warning in {self._out_filename}: The following thermochemistry values were not found:")
            for error in errors:
                print(f"  - {error}")
        else:
            print(f"Successfully extracted all thermochemistry data from {self._out_filename}.")

        return data
