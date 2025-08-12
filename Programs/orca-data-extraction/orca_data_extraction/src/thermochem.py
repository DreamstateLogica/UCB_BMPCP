#!/usr/bin/env python3
"""
A DataSection subclass that extracts and converts thermochemical data from an
ORCA output file.
"""
_author__ = "Nabil Mouhajir, Peter Waddell"
__copyright__ = "Copyright 2024"
__credits__ = ["Nabil Mouhajir, Peter Waddell"]
__version__ = "1.0"
__date__ = "2025/07/22"
__maintainer__ = "Nabil Mouhajir"
__email__ = "mouhajirni@hotmail.com"
__status__ = "Prototype"

import re
from data_section import DataSection

# Define physical constants for unit conversions
HARTREE_TO_JOULE = 4.3597447222071e-18
HARTREE_TO_EV = 27.211386245988
TEMPERATURE_K = 298.15


class ThermoData(DataSection):
    """
    Extracts and converts key thermochemical values from an ORCA output file.
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
        Extracts the thermochemical values, converts them to standard units,
        and returns them in a dictionary with updated labels.
        """
        raw_data = {}
        errors = []
        contents = self._outfile_contents

        # Define forward-facing regex patterns for each value in Hartrees.
        patterns = {
            'Final Gibbs Free Energy': re.compile(r"Final\s+Gibbs\s+free\s+energy\s*\.*\s*(-?[\d\.]+)\s*Eh", re.IGNORECASE),
            'Final Entropy Term': re.compile(r"Final\s+entropy\s+term\s*\.*\s*(-?[\d\.]+)\s*Eh", re.IGNORECASE),
            'Total Enthalpy': re.compile(r"Total\s+Enthalpy\s*\.*\s*(-?[\d\.]+)\s*Eh", re.IGNORECASE),
            'Zero Point Energy': re.compile(r"Zero\s+point\s+energy\s*\.*\s*(-?[\d\.]+)\s*Eh", re.IGNORECASE),
            'Electronic Energy': re.compile(r"Electronic\s+energy\s*\.*\s*(-?[\d\.]+)\s*Eh", re.IGNORECASE)
        }

        print(f"--- Analyzing {self._out_filename} for thermochemistry ---")
        for key, pattern in patterns.items():
            matches = pattern.findall(contents)
            if matches:
                # The last match in the file is the final value.
                raw_data[key] = matches[-1]
            else:
                raw_data[key] = 'N/A'
                errors.append(f"Could not find '{key}'")

        # --- Perform Unit Conversions ---
        converted_data = {}

        # Gibbs Free Energy (Eh -> J)
        try:
            val_eh = float(raw_data['Final Gibbs Free Energy'])
            converted_data['Final Gibbs Free Energy (J)'] = val_eh * HARTREE_TO_JOULE
        except (ValueError, TypeError):
            converted_data['Final Gibbs Free Energy (J)'] = 'N/A'

        # Final Entropy Term (T*S in Eh -> S in J/K)
        try:
            val_eh = float(raw_data['Final Entropy Term'])
            # ORCA outputs T*S in Eh. Convert to Joules, then divide by T.
            converted_data['Final Entropy (J/K)'] = (val_eh * HARTREE_TO_JOULE) / TEMPERATURE_K
        except (ValueError, TypeError):
            converted_data['Final Entropy (J/K)'] = 'N/A'

        # Total Enthalpy (Eh -> J)
        try:
            val_eh = float(raw_data['Total Enthalpy'])
            converted_data['Total Enthalpy (J)'] = val_eh * HARTREE_TO_JOULE
        except (ValueError, TypeError):
            converted_data['Total Enthalpy (J)'] = 'N/A'

        # Zero Point Energy (Eh -> eV)
        try:
            val_eh = float(raw_data['Zero Point Energy'])
            converted_data['Zero Point Energy (eV)'] = val_eh * HARTREE_TO_EV
        except (ValueError, TypeError):
            converted_data['Zero Point Energy (eV)'] = 'N/A'
            
        # Electronic Energy (Eh -> eV)
        try:
            val_eh = float(raw_data['Electronic Energy'])
            converted_data['Electronic Energy (eV)'] = val_eh * HARTREE_TO_EV
        except (ValueError, TypeError):
            converted_data['Electronic Energy (eV)'] = 'N/A'
        
        # Report errors if any values were not found
        if errors:
            print(f"Warning in {self._out_filename}: The following thermochemistry values were not found:")
            for error in errors:
                print(f"  - {error}")
        elif 'N/A' not in raw_data.values():
            print(f"Successfully extracted and converted all thermochemistry data from {self._out_filename}.")
        
        return converted_data