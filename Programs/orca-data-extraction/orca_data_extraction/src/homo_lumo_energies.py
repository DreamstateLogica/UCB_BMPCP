#!/usr/bin/env python3
"""
A DataSection subclass that extracts HOMO/LUMO energies and calculates
related conceptual DFT descriptors from an ORCA output file.
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

# Conversion from Hartrees to electronvolts
HARTREE_TO_EV = 27.211386245988

class HOMOLUMOEnergies(DataSection):
    """
    Finds HOMO/LUMO energies and calculates electronegativity, hardness,
    softness, and electrophilicity.
    """
    def __init__(self, out_filename, outfile_contents):
        """
        Parameters
        ----------
        out_filename : str
            Name of the ORCA .out file that will be searched.
        outfile_contents : str
            String containing the full text of the ORCA .out file.
        """
        super().__init__(out_filename, outfile_contents)
        self._section_name = 'Conceptual DFT Descriptors'

    def _find_data(self):
        """
        Extracts HOMO/LUMO, calculates descriptors, and returns them as a dict.

        Returns
        -------
        dict
            Dictionary containing all calculated values.
        """
        # Regex to find the last "ORBITAL ENERGIES" block in the file.
        orbital_blocks = re.findall(
            r"ORBITAL ENERGIES\n-+\n\n\s+NO\s+OCC\s+E\(Eh\)\s+E\(eV\)\s*\n(.*?)\n\n",
            self._outfile_contents,
            re.DOTALL
        )
        
        if not orbital_blocks:
            print(f"Warning in {self._out_filename}: Could not find any ORBITAL ENERGIES block.")
            return self._return_na_dict()

        last_block = orbital_blocks[-1]
        lines = last_block.strip().split('\n')

        homo_eh = None
        lumo_eh = None

        # Find the HOMO-LUMO transition
        for i, line in enumerate(lines):
            parts = line.split()
            try:
                occ = float(parts[1])
                if occ > 0:
                    # This is an occupied orbital, potentially the HOMO
                    homo_eh = float(parts[2])
                else:  # occ == 0.0
                    # This is the first unoccupied orbital, so it's the LUMO
                    lumo_eh = float(parts[2])
                    break  # Found both, no need to continue
            except (ValueError, IndexError):
                continue  # Skip malformed lines

        if homo_eh is None or lumo_eh is None:
            print(f"Warning in {self._out_filename}: Could not determine HOMO and/or LUMO energies.")
            return self._return_na_dict()
            
        # --- Calculate Descriptors ---
        try:
            # Convert to eV for conceptual DFT calculations
            homo_ev = homo_eh * HARTREE_TO_EV
            lumo_ev = lumo_eh * HARTREE_TO_EV

            electronegativity = -(lumo_ev + homo_ev) / 2.0
            hardness = (lumo_ev - homo_ev) / 2.0
            
            if hardness == 0:
                softness = float('inf')
                electrophilicity = float('inf')
            else:
                softness = 1.0 / hardness
                electrophilicity = (electronegativity**2) / (2.0 * hardness)
        except Exception as e:
            print(f"Warning in {self._out_filename}: Could not calculate conceptual DFT descriptors. Error: {e}")
            return self._return_na_dict()

        print(f"Successfully extracted conceptual DFT data from {self._out_filename}.")
        return {
            'HOMO (eV)': homo_ev,
            'LUMO (eV)': lumo_ev,
            'Electronegativity (eV)': electronegativity,
            'Hardness (eV)': hardness,
            'Softness (1/eV)': softness,
            'Electrophilicity (eV)': electrophilicity
        }

    def _return_na_dict(self):
        """Returns a dictionary with 'N/A' for all descriptor values."""
        return {
            'HOMO (eV)': 'N/A',
            'LUMO (eV)': 'N/A',
            'Electronegativity (eV)': 'N/A',
            'Hardness (eV)': 'N/A',
            'Softness (1/eV)': 'N/A',
            'Electrophilicity (eV)': 'N/A'
        }
