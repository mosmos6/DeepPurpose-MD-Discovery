# utils.py
# Author: Iori Mochizuki
# Created: 2025-05-13
# Shared helper functions for DeepPurpose-MD-Discovery pipeline

import numpy as np
import re

def extract_coordinates(pdb_file):
    """
    Extract atomic coordinates and full lines from a PDB or PDBQT file.

    Args:
        pdb_file (str): Path to the file.

    Returns:
        coords (np.ndarray): Nx3 array of atomic positions.
        atoms (list): List of tuples (atom_name, x, y, z, original_line).
    """
    atoms = []
    coords = []
    with open(pdb_file, "r") as infile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    atom_name = line[12:16].strip()
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    atoms.append((atom_name, x, y, z, line))
                    coords.append([x, y, z])
                except ValueError:
                    print(f"⚠️ Skipping malformed line: {line.strip()}")
                    continue
    return np.array(coords, dtype=float), atoms


def kabsch(P, Q):
    """
    Kabsch algorithm to find optimal rotation matrix to align P onto Q.

    Args:
        P (np.ndarray): Nx3 coordinate matrix (source).
        Q (np.ndarray): Nx3 coordinate matrix (target).

    Returns:
        R (np.ndarray): 3x3 rotation matrix.
        centroid_P, centroid_Q: Mean of input coordinates.
    """
    centroid_P = np.mean(P, axis=0)
    centroid_Q = np.mean(Q, axis=0)
    P_centered = P - centroid_P
    Q_centered = Q - centroid_Q

    H = P_centered.T @ Q_centered
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T

    # Avoid reflection
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T

    return R, centroid_P, centroid_Q


def fix_pdb_element_column(input_pdb, output_pdb):
    """
    Fix the element column (77–78) in a PDB file if it's missing or incorrect.

    Args:
        input_pdb (str): Path to PDB input.
        output_pdb (str): Destination path.
    """
    corrected_lines = []

    element_map = {
        "C": " C", "N": " N", "O": " O", "H": " H", "S": " S", "P": " P"
    }

    with open(input_pdb, "r") as infile:
        for line in infile:
            if line.startswith(("ATOM", "HETATM")):
                atom_type = line[12:14].strip()[0]
                atom_name = element_map.get(atom_type.upper(), " X")
                corrected_line = line[:76] + f"{atom_name:>2}" + line[78:]
                corrected_lines.append(corrected_line)
            else:
                corrected_lines.append(line)

    with open(output_pdb, "w") as outfile:
        outfile.writelines(corrected_lines)

    print(f"✅ Fixed PDB saved as {output_pdb}")


def calculate_centroid_from_pdb(pdb_path, atom_type="ATOM"):
    """
    Calculate the centroid (x, y, z) of all ATOM records in a PDB file.

    Args:
        pdb_path (str): Path to the receptor PDB.
        atom_type (str): "ATOM" or "HETATM".

    Returns:
        dict: {'x': float, 'y': float, 'z': float}
    """
    x_sum = y_sum = z_sum = 0.0
    atom_count = 0

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith(atom_type):
                x_sum += float(line[30:38].strip())
                y_sum += float(line[38:46].strip())
                z_sum += float(line[46:54].strip())
                atom_count += 1

    if atom_count == 0:
        raise ValueError(f"No {atom_type} entries found in {pdb_path}.")

    return {
        "x": round(x_sum / atom_count, 3),
        "y": round(y_sum / atom_count, 3),
        "z": round(z_sum / atom_count, 3),
    }
