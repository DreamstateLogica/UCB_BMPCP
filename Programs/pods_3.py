#!/usr/bin/env python3


## Author ##

# 👤 Maintainer and Author **Nabil Mouhajir**

# * Website: https://github.com/DreamstateLogica
# * Github: [@DreamstateLogica] (https://github.com/DreamstateLogica)
# * LinkedIn: https://www.linkedin.com/in/nabil-mouhajir/

# 👤 Maintainer and Co-Author **Rohan Wallace**

# * Website: https://github.com/Fierypigz
# * Github: [@Fierypigz] (https://github.com/Fierypigz)
# * LinkedIn: https://www.linkedin.com/in/rohan-wallace-099283234/

# 👤 Maintainer and Co-Author **Sean Treacy**

# * LinkedIn: https://www.linkedin.com/in/sean-treacy

# 👤 Original Author, Licenser **Peter Waddell**

# * Website: https://github.com/pmwaddell
# * Github: [@pmwaddell](https://github.com/pmwaddell)
# * LinkedIn: [@https://www.linkedin.com/in/pmwaddell-ph-d/](https://www.linkedin.com/in/pmwaddell-ph-d/)

## END ##

import os
import subprocess
import re
import argparse
import time
import shutil # For shutil.which

# Need os module for path operations etc.
import os

# --- RDKit Import ---
try:
    # Ensure RDKit modules are imported for helper functions
    from rdkit import Chem
    from rdkit.Chem import AllChem
    RDKIT_AVAILABLE = True
except ImportError:
    print("ERROR: RDKit module not found. Please install RDKit (e.g., via Conda).")
    RDKIT_AVAILABLE = False
    # Define dummy classes simply to allow script parsing if RDKit is missing.
    class Chem: pass
    class AllChem: pass
# --- End RDKit Import ---

# --- Configuration ---
# Updated DEFAULT_ORCA_PATH for Windows
DEFAULT_ORCA_PATH = "C:\\ORCA\\orca.exe" # Typical Windows ORCA path, adjust as needed
# Base keywords define the calculation type (Opt or Opt+Freq)
OPT_FREQ_BASE = "Opt Freq" # Keyword part for Gibbs energy calc
OPT_BASE = "Opt"           # Keyword part for Electronic energy calc
# Standard options for optimization runs (SlowConv added via flag)
OTHER_OPTS = "XYZFile"

HARTREE_TO_EV = 27.211386
DEFAULT_TIMEOUT_SECONDS = 7200 * 2 # Increased default timeout slightly for Opt(+Freq)

# --- Constants for Redox Potential Calculation (User Literature Method) ---
HARTREE_TO_KCALMOL = 627.5095   # Hartree to kcal/mol
FARADAY_KCAL_MOL_V = 23.0605    # Faraday constant in kcal mol^-1 V^-1
NE = 1.0                       # Number of electrons transferred (n_e)
E_ABS_SHE_V = 4.281            # Absolute SHE potential [V] from user's formula context
# ** WARNING: Dictionary below contains APPROXIMATE values - VERIFY FROM LITERATURE **
# ** Values are E_ref(SCE) vs SHE in that solvent [V] ("solvent correction" in user formula) **
# ** Ensure the SIGN convention matches the user's formula: E = Term - E_abs - E_ref **
SCE_REF_POTENTIALS_VS_SHE = {
    # Format: 'solvent_name_lowercase': E_ref(SCE) vs SHE in that solvent [V]
    'acetonitrile': -0.141, # Value provided previously by user
    'methanol':    -0.21,  # APPROXIMATE - CHECK LITERATURE! Often Ag/AgCl used instead.
    'water':       +0.244, # APPROXIMATE (vs SHE(aq)) - CHECK LITERATURE!
    'dmso':        -0.28,  # APPROXIMATE - CHECK LITERATURE!
    'thf':         -0.40,  # APPROXIMATE - CHECK LITERATURE!
    # Add more solvents and *verified* literature values here
}
# --- End Configuration ---


# --- Helper Functions --- (Defined Globally)

def smiles_to_xyz(smiles_string, molecule_name="molecule"):
    """Converts SMILES to XYZ coordinates using RDKit."""
    if not RDKIT_AVAILABLE: print(f"RDKit not available..."); return None, None
    try:
        from rdkit import Chem; from rdkit.Chem import AllChem # Ensure loaded
        mol = Chem.MolFromSmiles(smiles_string);
        if not mol: print(f"Error: Could not parse SMILES: {smiles_string}"); return None, None
        mol = Chem.AddHs(mol); cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, pruneRmsThresh=0.5, randomSeed=42)
        if not cids: cids = AllChem.EmbedMultipleConfs(mol, numConfs=10, useRandomCoords=True, pruneRmsThresh=0.5, randomSeed=42)
        if not cids:
             print(f"Warning: Could not generate 3D conformers for {smiles_string}. Trying basic UFF.")
             try: AllChem.EmbedMolecule(mol, randomSeed=42); AllChem.UFFOptimizeMolecule(mol)
             except Exception as e_uff: print(f"Error: UFF fallback failed for {smiles_string}: {e_uff}"); return None, None
             if mol.GetNumConformers() == 0: print(f"Error: UFF fallback also failed (no conformers) for {smiles_string}"); return None, None
             results = []
        else:
             try: results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=200)
             except ValueError: print(f"Warning: MMFF optimization failed for {smiles_string} (ValueError), trying UFF."); results = []
             except Exception as e_mmff: print(f"Warning: MMFF optimization failed for {smiles_string} ({type(e_mmff).__name__}), trying UFF."); results = []
             if not results:
                  try: results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=200)
                  except Exception as e_uff_conf: print(f"Error: UFF conformer optimization failed for {smiles_string}: {e_uff_conf}"); results = []
        if mol.GetNumConformers() == 0: print(f"Error: No conformers available for {smiles_string}"); return None, None
        if not results: print(f"Warning: Using first embedded conformer for {smiles_string} due to optimization failures."); conf_id = mol.GetConformers()[0].GetId()
        else:
             converged_results = [(i, res[1]) for i, res in enumerate(results) if len(res)>1 and res[0] == 0 and isinstance(res[1], (int, float))]
             if not converged_results:
                  print(f"Warning: No conformers converged for {smiles_string}. Using lowest energy found.")
                  all_energies = [(i, res[1]) for i, res in enumerate(results) if len(res)>1 and isinstance(res[1], (int, float))]
                  if not all_energies: print(f"Error: No valid conformer energies found for {smiles_string}. Using first."); min_energy_idx = 0
                  else: min_energy_idx = min(all_energies, key=lambda item: item[1])[0]
             else: min_energy_idx = min(converged_results, key=lambda item: item[1])[0]
             try: conf_id = mol.GetConformers()[min_energy_idx].GetId()
             except IndexError: print(f"Error: Conformer index {min_energy_idx} out of range. Using first."); conf_id = mol.GetConformers()[0].GetId()
        conf = mol.GetConformer(id=conf_id)
        xyz_coords = f"{mol.GetNumAtoms()}\n{molecule_name}\n"
        for i, atom in enumerate(mol.GetAtoms()): pos = conf.GetAtomPosition(i); xyz_coords += f"{atom.GetSymbol()} {pos.x:.8f} {pos.y:.8f} {pos.z:.8f}\n"
        return xyz_coords, mol # Return mol object
    except Exception as e: print(f"Error during RDKit processing for SMILES {smiles_string}: {e}"); return None, None


def generate_orca_input(full_input_path, xyz_coords_str, charge, multiplicity, keywords, num_orca_cores=None):
    """Creates an ORCA input file at the specified full path."""
    if not keywords.strip().startswith('!'): keywords = "! " + keywords.strip()
    try:
        os.makedirs(os.path.dirname(full_input_path), exist_ok=True)
        with open(full_input_path, 'w') as f:
            f.write(f"{keywords}\n");
            if num_orca_cores and num_orca_cores > 1: f.write(f"%pal nprocs {num_orca_cores} end\n")
            f.write(f"* xyz {charge} {multiplicity}\n");
            # Write coordinates - expects xyz_coords_str to be multiline string starting with num_atoms
            f.write("\n".join(xyz_coords_str.strip().split('\n')[2:]) + "\n") # Skip num_atoms and title line
            f.write("*\n")
    except Exception as e: print(f"ERROR: Failed to write ORCA input file {full_input_path}: {e}"); raise


def run_orca(orca_executable_path, full_input_path, full_base_name, working_dir, timeout=DEFAULT_TIMEOUT_SECONDS):
    """Runs ORCA calculation in a specified working directory, passes environment."""
    input_basename = os.path.basename(full_input_path); output_basename = os.path.basename(full_base_name) + ".out"; full_output_path = full_base_name + ".out"; command = [orca_executable_path, input_basename]
    env = os.environ.copy() # Pass environment copy
    print(f"Running ORCA in '{working_dir}': {' '.join(command)} > {output_basename}")
    try:
        # On Windows, shell=True might be needed for some executables if they are in PATH but not found directly,
        # or if ORCA relies on batch scripts/environment variables.
        # However, for direct executable paths, it's generally not needed and can be less secure.
        # Keeping shell=False as it's typically more robust if the full path is correct.
        result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True, timeout=timeout, cwd=working_dir, env=env) # Added env=env
        with open(full_output_path, 'w') as outfile: outfile.write(result.stdout)
        if result.stderr: print(f"ORCA STDERR for {output_basename}:\n{result.stderr}")
        print(f"ORCA calculation finished successfully for {output_basename}.")
        return full_output_path
    except FileNotFoundError: print(f"ERROR: ORCA executable not found at {orca_executable_path}"); return None
    except subprocess.TimeoutExpired: print(f"ERROR: ORCA calculation timed out (> {timeout}s) for {output_basename}"); try_write_output_on_error(full_output_path, f"TIMEOUT ERROR after {timeout}s"); return None
    except subprocess.CalledProcessError as e: print(f"ERROR: ORCA process failed for {output_basename} with exit code {e.returncode}."); try_write_output_on_error(full_output_path, e.stdout, e.stderr); return None
    except Exception as e: print(f"An unexpected error occurred during ORCA run for {output_basename}: {e}"); try_write_output_on_error(full_output_path, f"Unexpected Python error: {e}"); return None


def try_write_output_on_error(full_output_path, stdout_content=None, stderr_content=None):
    """Helper to write ORCA output/error even if run failed"""
    try:
        output_dir = os.path.dirname(full_output_path); os.makedirs(output_dir, exist_ok=True)
        with open(full_output_path, 'w') as f_err:
            f_err.write(f"--- ORCA RUN FAILED OR TIMED OUT ({time.strftime('%Y-%m-%d %H:%M:%S')}) ---\n")
            if stdout_content: f_err.write("\n--- STDOUT ---\n"); f_err.write(stdout_content)
            if stderr_content: f_err.write("\n--- STDERR ---\n"); f_err.write(stderr_content)
            f_err.write("\n--- END OF CAPTURED OUTPUT ---")
    except Exception as e_write: print(f"Could not write error details to {os.path.basename(full_output_path)}: {e_write}")


def parse_orca_energy(full_output_path):
    """
    Parses the final energy (Gibbs Free Energy if available [e.g., from Opt+Freq],
    else Final Single Point E_el [e.g., from Opt or SP])
    from ORCA output file given its full path.
    """
    energy = None; gibbs_energy = None; final_sp_energy = None; normal_termination = False; convergence_info = "Convergence status unknown"
    if not os.path.isfile(full_output_path): print(f"Error: Output file for parsing not found: {os.path.basename(full_output_path)}"); return None
    try:
        content = "";
        with open(full_output_path, 'r') as f: content = f.read()
        if "****ORCA TERMINATED NORMALLY****" in content: normal_termination = True
        # Prioritize Gibbs Free Energy (relevant for Opt+Freq outputs)
        gibbs_match = re.search(r"Final Gibbs free energy\s+\.+\s+([-+]?\d+\.\d+)", content)
        if gibbs_match: gibbs_energy = float(gibbs_match.group(1))
        # Fallback to Final Single Point Energy (relevant for SP or Opt outputs)
        sp_match = re.search(r"FINAL SINGLE POINT ENERGY\s+([-+]?\d+\.\d+)", content)
        if sp_match: final_sp_energy = float(sp_match.group(1))

        # Check convergence status
        if "OPTIMIZATION HAS CONVERGED" in content: convergence_info = "Optimization Converged"
        elif "OPTIMIZATION DID NOT CONVERGE" in content: convergence_info = "Optimization NOT Converged"
        elif "SCF NOT CONVERGED" in content: convergence_info = "SCF NOT CONVERGED"

        energy = gibbs_energy if gibbs_energy is not None else final_sp_energy
        return energy
    except Exception as e: print(f"Error parsing {os.path.basename(full_output_path)}: {e}"); return None


def check_orca_completion(full_output_path):
    """Checks if an ORCA output file exists and contains the normal termination string."""
    if not os.path.isfile(full_output_path): return False
    try:
        file_size = os.path.getsize(full_output_path); read_size = min(file_size, 4096)
        if file_size == 0: return False
        with open(full_output_path, 'r') as f:
            if read_size > 0: f.seek(max(0, file_size - read_size), os.SEEK_SET); tail = f.read()
            else: f.seek(0); tail = f.read()
            if "****ORCA TERMINATED NORMALLY****" in tail: return True
    except Exception as e: return False
    return False


def check_orca_convergence(full_output_path):
    """Checks if an ORCA optimization output file contains the convergence string."""
    if not os.path.isfile(full_output_path): return False
    try:
        with open(full_output_path, 'r') as f_check: content = f_check.read()
        if "OPTIMIZATION HAS CONVERGED" in content or \
           "The optimization did not converge but the final gradient is very small" in content:
             if "OPTIMIZATION HAS CONVERGED" in content:
                  return True
             else:
                  return False # Treat as not converged
    except Exception as e: print(f"Warning: could not read output file {os.path.basename(full_output_path)} for convergence check: {e}"); return False
    return False


def check_imaginary_freq(full_output_path):
    """
    Checks ORCA output file (from Opt Freq) for imaginary frequencies.
    Returns True if imaginary frequencies are found, False otherwise or if check fails.
    """
    if not os.path.isfile(full_output_path):
        print(f"Warning: Cannot check for imaginary frequencies, file not found: {os.path.basename(full_output_path)}")
        return False

    try:
        with open(full_output_path, 'r') as f:
            in_freq_section = False; imag_freq_count = 0
            for line in f:
                if line.strip() == "VIBRATIONAL FREQUENCIES" or \
                   line.strip() == "*** Vibrational Frequencies":
                    in_freq_section = True
                    continue
                if in_freq_section:
                    match = re.match(r"\s*\d+\s*:\s*(-?\d+\.\d+)\s+cm", line)
                    if match:
                        frequency = float(match.group(1));
                        if frequency < -1.0:
                           imag_freq_count += 1
                    if "NORMAL MODES" in line or \
                       "Thermochemistry" in line or \
                       "THERMOCHEMISTRY AT" in line or \
                       "ZERO POINT ENERGY" in line or \
                       "Writing the Hessian" in line or \
                       "End of Frequency Analysis" in line or \
                       "IR SPECTRUM" in line or \
                       line.strip() == "-----------------------" or \
                       line.strip() == "=========================":
                         break

            if imag_freq_count > 0:
                 print(f"    Warning: Found {imag_freq_count} imaginary frequency/frequencies in {os.path.basename(full_output_path)}.")
                 return True

    except Exception as e:
        print(f"Warning: Could not read/check {os.path.basename(full_output_path)} for imaginary frequencies: {e}")
        return False
    return False


def read_xyz_file(xyz_filepath):
    """Reads coordinates from an XYZ file."""
    if not os.path.isfile(xyz_filepath): raise ValueError(f"XYZ file not found: {os.path.basename(xyz_filepath)}")
    try:
        with open(xyz_filepath, 'r') as f: lines = f.readlines()
        if len(lines) < 2: raise ValueError(f"XYZ file is too short (must have count and title lines): {os.path.basename(xyz_filepath)}")
        try: num_atoms = int(lines[0].strip())
        except ValueError: raise ValueError(f"Cannot read atom count from first line of {os.path.basename(xyz_filepath)}")
        if len(lines) != num_atoms + 2: raise ValueError(f"Line count ({len(lines)}) in {os.path.basename(xyz_filepath)} does not match atom count ({num_atoms}) + 2.")
        coord_pattern = re.compile(r"^\s*[A-Za-z]{1,3}\s+-?\d+(\.\d*)?([eE][-+]?\d+)?\s+-?\d+(\.\d*)?([eE][-+]?\d+)?\s+-?\d+(\.\d*)?([eE][-+]?\d+)?")
        coord_lines = lines[2:];
        if not all(coord_pattern.match(line) for line in coord_lines): bad_lines_info = [f"L{k+3}: {line.strip()}" for k, line in enumerate(coord_lines) if not coord_pattern.match(line)]; raise ValueError(f"Invalid coordinate line format found in {os.path.basename(xyz_filepath)}.\nProblematic line(s):\n" + "\n".join(bad_lines_info))
        return "".join(lines)
    except ValueError as ve: raise ve
    except Exception as e: raise ValueError(f"Error reading or parsing XYZ file {os.path.basename(xyz_filepath)}: {e}")


# --- Main Prediction Function (Handles Mode and Energy Type) ---
def predict_redox(smiles_list, functional, basis_set, solvent_name, orca_path, main_output_dir,
                  mode='both', energy_type='electronic', dispersion='D3BJ',
                  use_tightscf=False, use_slowconv=False,
                  num_orca_cores=None):
    """
    Predicts redox potentials sequentially based on specified mode and energy type.
    Includes options for TightSCF and SlowConv keywords. Thank 
    MODIFIED: Returns results and a list of summary_lines for file output.
    """
    if not RDKIT_AVAILABLE:
        print("Cannot proceed without RDKit.")
        return [], ["ERROR: Cannot proceed without RDKit."] # MODIFICATION: Return empty summary
    if not smiles_list:
        print("Error: No SMILES strings provided to process.")
        return [], ["ERROR: No SMILES strings provided to process."] # MODIFICATION: Return empty summary

    # MODIFICATION START: Initialize summary_lines list
    summary_lines = []
    # MODIFICATION END

    results = []; start_time = time.time()
    # Use os.path.normpath to ensure consistent path separators for display
    output_dir_name = os.path.normpath(os.path.basename(main_output_dir))

    # MODIFICATION START: Add initial settings to summary_lines
    # Terminal prints for these are handled in __main__ before calling this function
    summary_lines.append(f"--- Redox Potential Calculation Summary ---")
    summary_lines.append(f"Timestamp: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    summary_lines.append(f"Calculation Mode: {mode}")
    summary_lines.append(f"Energy Type Used for Delta: {energy_type.capitalize()} {'(Runs Freq!)' if energy_type == 'gibbs' else '(Opt Only)'}")
    summary_lines.append(f"Optional ORCA Keywords: TightSCF={use_tightscf}, SlowConv={use_slowconv}")
    # MODIFICATION END

    print(f"--- Starting Sequential Processing ---") # Keep this as a console print
    print(f"Main output directory: {main_output_dir}") # Keep as console print

    solvent_name_cleaned = solvent_name.strip(); solvation_keyword = f"CPCM({solvent_name_cleaned})"
    # MODIFICATION START: Add to summary
    summary_lines.append(f"Using Solvation Model: {solvation_keyword}")
    # MODIFICATION END
    print(f"Using Solvation Model: {solvation_keyword}") # Keep console print

    solvent_key = solvent_name_cleaned.lower(); e_ref_rel_v = None; ref_electrode_name = "SCE"; ref_solvent_used = None; reference_potential_found = False
    if solvent_key in SCE_REF_POTENTIALS_VS_SHE:
        e_ref_rel_v = SCE_REF_POTENTIALS_VS_SHE[solvent_key]; ref_solvent_used = solvent_name_cleaned; reference_potential_found = True
        ref_info_line = f"Using reference potential for {ref_electrode_name} in {solvent_name_cleaned}: E_ref={e_ref_rel_v:.3f} V vs SHE"
        print(ref_info_line) # Keep console print
        summary_lines.append(ref_info_line) # MODIFICATION: Add to summary
    else:
        warn_line1 = f"Warning: Reference potential data for {ref_electrode_name} in '{solvent_name_cleaned}' not found in script's dictionary."
        print(warn_line1) # Keep console print
        summary_lines.append(warn_line1) # MODIFICATION: Add to summary
        fallback_solvent = 'acetonitrile'; e_ref_rel_v = SCE_REF_POTENTIALS_VS_SHE.get(fallback_solvent)
        if e_ref_rel_v is not None:
            ref_solvent_used = fallback_solvent; reference_potential_found = True;
            warn_line2 = f"         Using reference values for SCE in Acetonitrile as fallback (E_ref={e_ref_rel_v:.3f} V)."
            warn_line3 = f"         >>> Calculated potentials vs SCE will be inaccurate for {solvent_name_cleaned}! <<<"
            print(warn_line2); print(warn_line3) # Keep console prints
            summary_lines.append(warn_line2); summary_lines.append(warn_line3) # MODIFICATION: Add to summary
        else:
            err_line = f"         Acetonitrile reference data also missing! Cannot calculate potentials vs {ref_electrode_name}."
            print(err_line); summary_lines.append(err_line) # MODIFICATION: Add to summary
    reference_correction_v = (E_ABS_SHE_V - e_ref_rel_v) if reference_potential_found else None

    for i, smiles in enumerate(smiles_list):
        index = i + 1; smiles = smiles.strip()
        if not smiles: print(f"Skipping empty line at index {i}."); continue

        mol_subdir_name = f"mol_{index}"; mol_dir_path = os.path.join(main_output_dir, mol_subdir_name)
        # Use os.path.normpath for the relative path in summary
        relative_mol_dir = os.path.normpath(os.path.join(output_dir_name, mol_subdir_name))
        print(f"\n--- Processing Molecule {index}/{len(smiles_list)}: {smiles} ---")
        print(f"    Calculations in subdirectory: ./{relative_mol_dir}/")
        result_data = {"smiles": smiles, "index": index, "Energy_Neutral": None, "Energy_Cation": None, "Energy_Anion": None, "E_ox_vs_ref_V": None, "E_red_vs_ref_V": None, "Prediction": "Processing Started", "status": "Failed", "directory": relative_mol_dir, "energy_type": energy_type, "mode": mode}

        try:
            os.makedirs(mol_dir_path, exist_ok=True)
            mol_name_in_dir = "molecule"

            print(f"    Step 0: Generating initial 3D structure...")
            if not RDKIT_AVAILABLE: raise ImportError("RDKit cannot be imported")
            from rdkit import Chem; from rdkit.Chem import AllChem
            initial_xyz_coords_str, rdkit_mol = smiles_to_xyz(smiles, mol_name_in_dir)
            if not initial_xyz_coords_str or not rdkit_mol: raise ValueError("Structure Generation Failed")
            num_atoms = rdkit_mol.GetNumAtoms();
            if num_atoms <= 0: raise ValueError("RDKit molecule has zero atoms.")

            states_to_calculate = ["neutral"]
            if mode in ['both', 'oxidation']: states_to_calculate.append("cation")
            if mode in ['both', 'reduction']: states_to_calculate.append("anion")
            print(f"    Mode '{mode}': Will calculate states: {', '.join(states_to_calculate)}")

            neutral_charge = Chem.GetFormalCharge(rdkit_mol); num_radical_electrons = sum(atom.GetNumRadicalElectrons() for atom in rdkit_mol.GetAtoms())
            total_electrons_neutral = sum(atom.GetAtomicNum() for atom in rdkit_mol.GetAtoms()) - neutral_charge; neutral_multiplicity = 1 if (total_electrons_neutral + num_radical_electrons) % 2 == 0 else 2
            cation_charge = neutral_charge + 1; cation_multiplicity = 1 if (total_electrons_neutral - 1) % 2 == 0 else 2
            anion_charge = neutral_charge - 1; anion_multiplicity = 1 if (total_electrons_neutral + 1) % 2 == 0 else 2
            all_states_params = {
                 "neutral": {"charge": neutral_charge, "mult": neutral_multiplicity, "start_geom": initial_xyz_coords_str},
                 "cation":  {"charge": cation_charge,  "mult": cation_multiplicity,  "start_geom": None},
                 "anion":   {"charge": anion_charge,   "mult": anion_multiplicity,   "start_geom": None}
            }

            energies = {}; opt_xyz_files = {}; calculation_ok = True

            for state_name in states_to_calculate:
                 state_params = all_states_params[state_name]
                 job_type_str = "Opt+Freq" if energy_type == 'gibbs' else "Opt"
                 print(f"    Step {job_type_str}: Processing {state_name} state (Q={state_params['charge']}, M={state_params['mult']})...")

                 job_suffix = "opt_freq" if energy_type == 'gibbs' else "opt"
                 opt_freq_keyword = OPT_FREQ_BASE if energy_type == 'gibbs' else OPT_BASE
                 calc_input_file = os.path.join(mol_dir_path, f"{mol_name_in_dir}_{state_name}_{job_suffix}.inp")
                 calc_base_name = os.path.join(mol_dir_path, f"{mol_name_in_dir}_{state_name}_{job_suffix}")
                 calc_output_file = calc_base_name + ".out"
                 calc_xyz_file = calc_base_name + ".xyz"

                 keyword_list = ['!']
                 keyword_list.append(functional)
                 if dispersion and dispersion.lower() != 'none':
                     keyword_list.append(dispersion)
                 if use_tightscf:
                     keyword_list.append("TightSCF")
                 keyword_list.append(basis_set)
                 keyword_list.append(opt_freq_keyword)
                 if use_slowconv:
                     keyword_list.append("SlowConv")
                 keyword_list.append(OTHER_OPTS)
                 keyword_list.append(solvation_keyword)
                 keywords = " ".join(filter(None, keyword_list))

                 converged = False
                 has_imag_freq = (energy_type == 'gibbs')

                 start_geom_str = state_params["start_geom"]
                 if start_geom_str is None:
                      neutral_xyz_path = opt_xyz_files.get("neutral")
                      if not neutral_xyz_path: raise ValueError(f"Cannot start {state_name} {job_type_str}, missing neutral optimized xyz path")
                      try: start_geom_str = read_xyz_file(neutral_xyz_path)
                      except ValueError as e_read: raise ValueError(f"Failed to read neutral opt xyz '{neutral_xyz_path}' for {state_name}: {e_read}")

                 if check_orca_completion(calc_output_file):
                      print(f"    Skipping {state_name} {job_type_str}: Found completed file.")
                      converged = check_orca_convergence(calc_output_file)
                      if not converged: result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Skipped but Not Converged"; calculation_ok = False; break
                      print(f"    Existing {state_name} output file shows convergence.")
                      if energy_type == 'gibbs':
                           has_imag_freq = check_imaginary_freq(calc_output_file)
                           if has_imag_freq: result_data["Prediction"] = f"{state_name.capitalize()} Opt+Freq Found Imaginary Frequencies"; calculation_ok = False; break
                           else: print(f"    Existing {state_name} Opt+Freq file has no imaginary frequencies.")
                      energy = parse_orca_energy(calc_output_file)
                      if energy is None: result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Skipped but Failed (Energy Parse Error)"; calculation_ok = False; break
                      energies[state_name] = energy
                      if os.path.exists(calc_xyz_file): opt_xyz_files[state_name] = calc_xyz_file
                      print(f"    Parsed {'G' if energy_type=='gibbs' else 'E_el'} ({state_name}) from existing file: {energy:.6f} Hartree")

                 else:
                      generate_orca_input(calc_input_file, start_geom_str, state_params['charge'], state_params['mult'], keywords, num_orca_cores)
                      output_file_run = run_orca(orca_path, calc_input_file, calc_base_name, working_dir=mol_dir_path, timeout=DEFAULT_TIMEOUT_SECONDS)
                      if output_file_run is None or not os.path.exists(output_file_run): result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Failed (Run Error/No Output)"; calculation_ok = False; break
                      calc_output_file = output_file_run
                      if not check_orca_completion(calc_output_file): result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Failed (Abnormal Termination)"; calculation_ok = False; break
                      converged = check_orca_convergence(calc_output_file)
                      if not converged: result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Failed (Did Not Converge)"; calculation_ok = False; break
                      print(f"    {state_name.capitalize()} optimization converged.")
                      if energy_type == 'gibbs':
                           has_imag_freq = check_imaginary_freq(calc_output_file)
                           if has_imag_freq: result_data["Prediction"] = f"{state_name.capitalize()} Opt+Freq Found Imaginary Frequencies"; calculation_ok = False; break
                           else: print(f"    {state_name.capitalize()} Opt+Freq file has no imaginary frequencies.")
                      energy = parse_orca_energy(calc_output_file)
                      if energy is None: result_data["Prediction"] = f"{state_name.capitalize()} {job_suffix.upper()} Failed (Energy Parse Error)"; calculation_ok = False; break
                      energies[state_name] = energy
                      if os.path.exists(calc_xyz_file): opt_xyz_files[state_name] = calc_xyz_file
                      print(f"    Parsed {'G' if energy_type=='gibbs' else 'E_el'} ({state_name}) from new file: {energy:.6f} Hartree")

                 if state_name == "neutral":
                     if not opt_xyz_files.get("neutral"):
                          raise ValueError("Neutral optimization finished but corresponding XYZ file path was not stored or found!")
                 if not calculation_ok: break

            if calculation_ok:
                result_data["status"] = "Success (Energies Calculated)"
                result_data["Energy_Neutral"] = energies.get("neutral")
                if mode in ['both', 'oxidation']: result_data["Energy_Cation"] = energies.get("cation")
                if mode in ['both', 'reduction']: result_data["Energy_Anion"] = energies.get("anion")

                if reference_potential_found:
                    e_ox_calculated = False
                    if mode in ['both', 'oxidation'] and "neutral" in energies and "cation" in energies:
                        delta_E_half_ox_hartree = energies["neutral"] - energies["cation"]
                        delta_E_half_ox_kcal = delta_E_half_ox_hartree * HARTREE_TO_KCALMOL
                        term_ox_V = - delta_E_half_ox_kcal / (NE * FARADAY_KCAL_MOL_V)
                        e_ox_vs_ref = term_ox_V - E_ABS_SHE_V - e_ref_rel_v
                        result_data["E_ox_vs_ref_V"] = e_ox_vs_ref
                        e_ox_calculated = True

                    e_red_calculated = False
                    if mode in ['both', 'reduction'] and "neutral" in energies and "anion" in energies:
                        delta_E_half_red_hartree = energies["neutral"] - energies["anion"]
                        delta_E_half_red_kcal = delta_E_half_red_hartree * HARTREE_TO_KCALMOL
                        term_red_V = - delta_E_half_red_kcal / (NE * FARADAY_KCAL_MOL_V)
                        e_red_vs_ref = term_red_V - E_ABS_SHE_V - e_ref_rel_v
                        result_data["E_red_vs_ref_V"] = e_red_vs_ref
                        e_red_calculated = True

                    if e_ox_calculated or e_red_calculated:
                         result_data["Prediction"] = f"Potentials Calculated ({energy_type.capitalize()})"
                         result_data["status"] = "Success"
                         ref_label = f"{ref_electrode_name} in {ref_solvent_used}"
                         if e_ox_calculated: print(f"    Calculated E_ox = {e_ox_vs_ref:.3f} V vs {ref_label} (using {energy_type} energy)")
                         if e_red_calculated: print(f"    Calculated E_red = {e_red_vs_ref:.3f} V vs {ref_label} (using {energy_type} energy)")
                    elif calculation_ok: result_data["Prediction"] = f"Energies Calculated (Mode={mode}, No Potentials)"

                else:
                    result_data["Prediction"] = f"Energies Calculated (No Ref Potential Data for {solvent_name_cleaned})"
                    print(f"    Skipping potential calculation vs {ref_electrode_name} (no reference data for {solvent_name_cleaned}).")
                    result_data["status"] = "Success (No Potential)"
            else:
                 if result_data["Prediction"] == "Processing Started": result_data["Prediction"] = "Calculation Failed (Energy Calculation Error)"
                 print(f"    Failed: {result_data['Prediction']}")
                 result_data["status"] = "Failed"
            results.append(result_data)
        except Exception as e:
             print(f"    ERROR processing Molecule {index} ({smiles}): {e}")
             import traceback; traceback.print_exc()
             if result_data["Prediction"] == "Processing Started": result_data["Prediction"] = f"Calculation Failed ({type(e).__name__})"
             result_data["status"] = "Failed"; results.append(result_data); print(f"    Failed: {result_data['Prediction']}")
             continue

    end_time = time.time()
    total_time_line = f"\nTotal processing time: {end_time - start_time:.2f} seconds"
    print(total_time_line) # Keep console print
    summary_lines.append(total_time_line) # MODIFICATION: Add to summary

    # --- Analysis and Final Output (MODIFIED to append to summary_lines) ---
    summary_lines.append("\n--- Overall Results ---")
    valid_potential_results = [r for r in results if r.get("status").startswith("Success") and (r.get("E_ox_vs_ref_V") is not None or r.get("E_red_vs_ref_V") is not None)]
    if valid_potential_results:
        try:
            ox_results = [r for r in valid_potential_results if r.get("E_ox_vs_ref_V") is not None]
            if ox_results:
                min_eox_mol = min(ox_results, key=lambda x: x.get("E_ox_vs_ref_V", float('inf')))
                ref_label_summary = f"{ref_electrode_name} in {ref_solvent_used}"
                summary_lines.append(f"Molecule most easily OXIDIZED (Lowest E_ox): {min_eox_mol['smiles']} ({min_eox_mol['E_ox_vs_ref_V']:.3f} V vs {ref_label_summary})")
            else:
                summary_lines.append("No oxidation potentials were successfully calculated and referenced.")

            red_results = [r for r in valid_potential_results if r.get("E_red_vs_ref_V") is not None]
            if red_results:
                max_ered_mol = max(red_results, key=lambda x: x.get("E_red_vs_ref_V", float('-inf')))
                ref_label_summary = f"{ref_electrode_name} in {ref_solvent_used}"
                summary_lines.append(f"Molecule most easily REDUCED (Highest E_red): {max_ered_mol['smiles']} ({max_ered_mol['E_red_vs_ref_V']:.3f} V vs {ref_label_summary})")
            else:
                summary_lines.append("No reduction potentials were successfully calculated and referenced.")
        except (ValueError, TypeError):
            summary_lines.append("Could not determine min E_ox / max E_red.")
    elif any(r.get("status").startswith("Success") for r in results):
        summary_lines.append("Calculations finished for some molecules but potentials vs reference could not be determined or calculated.")
    else:
        summary_lines.append("No successful calculations were completed.")

    summary_lines.append("\n--- Detailed Results ---")
    summary_lines.append(f"* Potentials calculated using Delta {energy_type.capitalize()} and Literature Formula: E = - (E_N - E_Ion) / F - E_abs_SHE - E_ref_vs_SHE")
    if energy_type == 'electronic':
        summary_lines.append(f"* WARNING: Using Delta E_el (Opt only) is an approximation for Delta G_free and lacks frequency validation for minima.")
    if energy_type == 'gibbs':
        summary_lines.append(f"* Frequency calculation performed on required states for minimum validation.")
    if reference_potential_found and e_ref_rel_v is not None:
         ref_label_detail = f"{ref_electrode_name} in {ref_solvent_used}"
         summary_lines.append(f"* Reference Electrode: vs {ref_label_detail} (E_abs_SHE={E_ABS_SHE_V}V, E_ref={e_ref_rel_v:.3f}V)")
         if ref_solvent_used.lower() != solvent_name_cleaned.lower():
             summary_lines.append(f"* WARNING: Using fallback {ref_solvent_used} reference for {solvent_name_cleaned} solvent!")
    else:
        summary_lines.append(f"* Reference Electrode: Potentials vs reference could not be calculated for {solvent_name_cleaned}.")
    summary_lines.append(f"* WARNING: Ensure values in SCE_REF_POTENTIALS_VS_SHE dictionary are accurate and appropriate!")

    results.sort(key=lambda x: x.get('index', float('inf')))
    for res in results:
        summary_lines.append(f"Index: {res.get('index', 'N/A')} | SMILES: {res['smiles']}")
        # Use os.path.normpath when displaying directory for Windows style
        summary_lines.append(f"  Directory: .{os.sep}{os.path.normpath(res.get('directory', 'N/A'))}{os.sep}")
        energy_label = "G" if res.get('energy_type') == 'gibbs' else "E_el"
        if res.get('status').startswith("Success"):
            summary_lines.append(f"  {energy_label}_Neutral: {res.get('Energy_Neutral'):.6f} H" if res.get('Energy_Neutral') is not None else f"  {energy_label}_Neutral: N/A")
            if res.get('mode') in ['both', 'oxidation']:
                summary_lines.append(f"  {energy_label}_Cation:  {res.get('Energy_Cation'):.6f} H" if res.get('Energy_Cation') is not None else f"  {energy_label}_Cation:  N/A (Not Calc. or Failed)")
            if res.get('mode') in ['both', 'reduction']:
                summary_lines.append(f"  {energy_label}_Anion:   {res.get('Energy_Anion'):.6f} H" if res.get('Energy_Anion') is not None else f"  {energy_label}_Anion:   N/A (Not Calc. or Failed)")
            if res.get("E_ox_vs_ref_V") is not None:
                summary_lines.append(f"  E_ox vs Ref: {res['E_ox_vs_ref_V']:.3f} V")
            elif res.get('mode') in ['both', 'oxidation'] and res.get('status') != "Success (No Potential)" and res.get('status').startswith("Success"):
                summary_lines.append(f"  E_ox vs Ref: Not Calculated (Error)")
            if res.get("E_red_vs_ref_V") is not None:
                summary_lines.append(f"  E_red vs Ref: {res['E_red_vs_ref_V']:.3f} V")
            elif res.get('mode') in ['both', 'reduction'] and res.get('status') != "Success (No Potential)" and res.get('status').startswith("Success"):
                summary_lines.append(f"  E_red vs Ref: Not Calculated (Error)")
            if res.get("status") == "Success (No Potential)":
                summary_lines.append(f"  Potentials vs Ref: Not calculated (missing ref data for solvent)")
        else:
            error_detail = res.get("Prediction", "Unknown Failure")
            summary_lines.append(f"  Result: Calculation Failed ({error_detail})")
        summary_lines.append("-" * 20)
    # --- End of MODIFIED Detailed Results for summary_lines ---

    # MODIFICATION: Return summary_lines along with results
    return results, summary_lines


# --- Main Execution Block ---
if __name__ == "__main__":
    if not RDKIT_AVAILABLE: exit(1)

    parser = argparse.ArgumentParser(
        description="Predict redox potentials using ORCA (SEQUENTIAL). Allows selecting calculation mode and energy type.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("-s", "--smiles", nargs='+', help="SMILES strings (creates main dir 'smiles_run_TIMESTAMP/' in CWD).")
    input_group.add_argument("-i", "--input-file", dest='input_filename', help="Path to SMILES file (creates main dir 'FILENAME_BASE/' in CWD).")

    parser.add_argument("-f", "--functional", default="B3LYP", help="DFT functional (e.g., M062X, PBE0, wB97X-D).")
    parser.add_argument("-b", "--basis", default="def2-SVP", help="Basis set (e.g., def2-SVP, def2-TZVPD, ma-def2-SVP).")
    parser.add_argument("--solvent", default="Acetonitrile", help="Solvent name for CPCM model (e.g., Water, DMSO). Case-insensitive. Used for ORCA calc & selects Ref. Potential if available.")
    parser.add_argument(
        "--dispersion", default="D3BJ", choices=['None', 'D3BJ', 'D4'],
        help="Dispersion correction method ('None' for no explicit dispersion or if functional includes it like wB97X-D)."
    )
    parser.add_argument(
        "--mode", choices=['both', 'oxidation', 'reduction'], default='both',
        help="Which potential(s) to calculate: 'oxidation' (N+C only), 'reduction' (N+A only), or 'both'."
    )
    parser.add_argument(
        "--energy-type", choices=['gibbs', 'electronic'], default='electronic',
        help="Energy type: 'gibbs' (Opt+Freq, Delta G, expensive) or 'electronic' (Opt only, Delta E_el, faster approx.)."
    )
    parser.add_argument(
        "--tightscf", action='store_true',
        help="Use ORCA's TightSCF keyword for stricter SCF convergence."
    )
    parser.add_argument(
        "--slowconv", action='store_true',
        help="Use ORCA's SlowConv keyword for potentially more robust geometry optimization."
    )
    parser.add_argument("--orca_path", default=None, help=f"Full path to ORCA executable (default: auto-detect or '{DEFAULT_ORCA_PATH}').") # Updated default help
    parser.add_argument("--orca-cores", type=int, default=None, help="Cores for ORCA job (default: ORCA default).")

    args = parser.parse_args()

    ORCA_EXECUTABLE_PATH = None;
    if args.orca_path: ORCA_EXECUTABLE_PATH = args.orca_path; print(f"Using ORCA path provided by user: {ORCA_EXECUTABLE_PATH}")
    else: print("Attempting to auto-detect ORCA path..."); ORCA_EXECUTABLE_PATH = shutil.which("orca") # shutil.which works on Windows
    if ORCA_EXECUTABLE_PATH: print(f"Auto-detected ORCA path: {ORCA_EXECUTABLE_PATH}")
    else: ORCA_EXECUTABLE_PATH = DEFAULT_ORCA_PATH; print(f"Could not auto-detect ORCA. Falling back to default: {ORCA_EXECUTABLE_PATH}")
    # Check if the path ends with .exe and make sure it's executable
    if not ORCA_EXECUTABLE_PATH.lower().endswith(".exe"):
        # For Windows, if ORCA path is provided without .exe, try appending it.
        # This is a common issue on Windows where users might provide just "C:\ORCA\orca"
        # but the executable is actually "orca.exe".
        temp_orca_path = ORCA_EXECUTABLE_PATH + ".exe"
        if os.path.isfile(temp_orca_path) and os.access(temp_orca_path, os.X_OK):
            ORCA_EXECUTABLE_PATH = temp_orca_path
            print(f"Adjusted ORCA path to include .exe: {ORCA_EXECUTABLE_PATH}")

    if not os.path.isfile(ORCA_EXECUTABLE_PATH) or not os.access(ORCA_EXECUTABLE_PATH, os.X_OK): print("="*60 + "\nERROR: Invalid ORCA executable path." + f"\n Path determined: '{ORCA_EXECUTABLE_PATH}'"); exit(1)

    smiles_list = []; main_output_dir_name = None; original_dir = os.getcwd();
    if args.input_filename:
        abs_input_path = os.path.abspath(args.input_filename);
        if not os.path.isfile(abs_input_path): print(f"Error: Input file not found: {abs_input_path}"); exit(1)
        main_output_dir_name = os.path.splitext(os.path.basename(abs_input_path))[0]; print(f"Reading SMILES strings from file: {abs_input_path}")
        try:
            with open(abs_input_path, 'r') as f: smiles_list = [line.strip() for line in f if line.strip() and not line.strip().startswith(('#', ';'))]
            if not smiles_list: print(f"Error: No valid SMILES strings found in file: {args.input_filename}"); exit(1)
            print(f"Found {len(smiles_list)} SMILES strings in file.")
        except Exception as e: print(f"Error reading input file {args.input_filename}: {e}"); exit(1)
    elif args.smiles:
        smiles_list = args.smiles; timestamp = time.strftime("%Y%m%d_%H%M%S"); main_output_dir_name = f"smiles_run_{timestamp}"; print(f"Processing {len(smiles_list)} SMILES strings provided via command line.")
    if not smiles_list: print("Error: No SMILES strings to process."); exit(1)
    if not main_output_dir_name: print("Error: Could not determine main output directory name."); exit(1)
    main_output_path = os.path.join(original_dir, main_output_dir_name)

    try:
        print(f"Creating main output directory: {main_output_path}"); os.makedirs(main_output_path, exist_ok=True)

        # MODIFICATION START: Define summary file path and write initial settings
        summary_file_path = os.path.join(main_output_path, "calculation_summary.txt")
        with open(summary_file_path, 'w') as summary_file:
            summary_file.write(f"--- Calculation Settings Summary ---\n")
            summary_file.write(f"Timestamp: {time.strftime('%Y%m%d_%H%M%S')}\n")
            summary_file.write(f"Input source: {'File - ' + args.input_filename if args.input_filename else 'Direct SMILES input'}\n")
            summary_file.write(f"Number of SMILES processed: {len(smiles_list)}\n")
            summary_file.write(f"ORCA Executable: {ORCA_EXECUTABLE_PATH}\n")
            if args.orca_cores: summary_file.write(f"ORCA Cores per Job: {args.orca_cores}\n")
            else: summary_file.write("ORCA Cores per Job: ORCA default (usually 1)\n")
            summary_file.write(f"Functional: {args.functional}\n")
            summary_file.write(f"Basis Set: {args.basis}\n")
            summary_file.write(f"Dispersion: {args.dispersion}\n")
            summary_file.write(f"Solvent: {args.solvent}\n")
            summary_file.write(f"Redox Calculation Mode: {args.mode}\n")
            summary_file.write(f"Energy Type Used for Delta: {args.energy_type} {'(Uses Freq!)' if args.energy_type == 'gibbs' else '(Opt Only)'}\n")
            summary_file.write(f"Optional ORCA Keywords: TightSCF={args.tightscf}, SlowConv={args.slowconv}\n")
            summary_file.write(f"Main Output Directory: {os.path.normpath(main_output_path)}\n") # Normalize path for summary
            summary_file.write("-" * 40 + "\n\n")
        # MODIFICATION END

        # Print settings to console (kept from original script)
        try: cpu_cores = os.cpu_count()
        except NotImplementedError: cpu_cores = None; print("Warning: Could not determine number of CPU cores.")
        if args.orca_cores and cpu_cores and args.orca_cores > cpu_cores: print(f"WARNING: --orca-cores ({args.orca_cores}) > available cores ({cpu_cores}).")
        elif args.orca_cores and args.orca_cores > 1: print(f"Note: Each ORCA job will attempt to use {args.orca_cores} cores.")
        print(f"\n--- Starting Sequential Calculations ---") # This part is more for console
        # The following lines are already effectively covered by the summary file preamble for predict_redox
        # and also the console prints within predict_redox.
        # No need to duplicate them here in console if predict_redox also prints them.

        # MODIFICATION START: Call predict_redox and get summary_lines
        # Note: predict_redox itself will print its own progress to console.
        # The summary_lines it returns will be appended to the summary_file.
        results, summary_output_lines = predict_redox(
                                        smiles_list,
                                        args.functional,
                                        args.basis,
                                        args.solvent,
                                        ORCA_EXECUTABLE_PATH,
                                        main_output_path,
                                        mode=args.mode,
                                        energy_type=args.energy_type,
                                        dispersion=args.dispersion,
                                        use_tightscf=args.tightscf,
                                        use_slowconv=args.slowconv,
                                        num_orca_cores=args.orca_cores)

        # Append results from predict_redox to the summary file
        with open(summary_file_path, 'a') as summary_file: # Append mode
            for line in summary_output_lines:
                summary_file.write(line + "\n")

        print(f"\nSummary of calculations written to: {summary_file_path}")
        # MODIFICATION END

    except Exception as e: print(f"\nFATAL ERROR in main execution block: {e}"); import traceback; traceback.print_exc()
    finally:
        print(f"\nFinished calculations. Main output directory: {main_output_path}")
        print(f"Current directory is still: {os.getcwd()}")