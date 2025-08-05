# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 18:33:24 2025

@author: georg
"""

import spglib
import numpy as np
from pymatgen.core import Structure
import sympy as sp

tolerance = 1e-3

def extract_cif_data(file_path):
    """
    Extract lattice vectors, atomic fractional coordinates, and atomic numbers from a CIF file.
    
    Args:
        file_path (str): Path to the CIF file.
        
    Returns:
        lattice_vectors (numpy.ndarray): 3×3 array of lattice vectors.
        atomic_fractional_coords (numpy.ndarray): N×3 array of fractional atomic positions.
        atomic_numbers (list[int]): List of atomic numbers for each site.
    """
    structure = Structure.from_file(file_path)
    lattice_vectors = structure.lattice.matrix
    atomic_fractional_coords = structure.frac_coords
    atomic_numbers = structure.atomic_numbers
    return lattice_vectors, atomic_fractional_coords, atomic_numbers

def convert_lattice_vectors_to_NumPy_array(lattice_vectors):
    """
    Convert lattice vectors to a NumPy array.
    
    Args:
        lattice_vectors (array-like): Sequence of three 3‐element lattice vectors.
        
    Returns:
        lattice_vectors_np (numpy.ndarray): 3×3 array of lattice vectors.
    """
    lattice_vectors_np = np.array([np.array(lv) for lv in lattice_vectors])
    return lattice_vectors_np

def convert_atomic_fractional_coords_to_NumPy_array(atomic_fractional_coords):
    """
    Convert fractional atomic coordinates to a NumPy array.
    
    Args:
        atomic_fractional_coords (array-like): Sequence of N fractional coordinate triplets.
    
    Returns:
        atomic_fractional_coords_np (numpy.ndarray): N×3 array of fractional coordinates.
    """
    atomic_fractional_coords_np = np.array([np.array(pos) for pos in atomic_fractional_coords])
    return atomic_fractional_coords_np

def get_spacegroup(symmetry_dataset):
    """
    Get the international space‐group symbol from a symmetry dataset.

    Args:
        symmetry_dataset (SpglibDataset): Output of spglib.get_symmetry_dataset.

    Returns:
        str: International space‐group symbol (e.g., 'Im-3m').
    """
    return symmetry_dataset.international

def find_potential_point_group_operations(symmetry_dataset):
    """
    Retrieve all rotation operations (potential point‐group ops) from a symmetry dataset.

    Args:
        symmetry_dataset (SpglibDataset): Output of spglib.get_symmetry_dataset.

    Returns:
        rotations (numpy.ndarray): Array of shape (M, 3, 3) containing rotation matrices.
    """
    rotations = symmetry_dataset.rotations
    return rotations

def convert_fractional_coordinates_to_cartesian(atomic_fractional_coords, lattice_vectors):
    """
    Convert fractional coordinates to Cartesian coordinates.

    Args:
        atomic_fractional_coords (array-like): N×3 fractional coordinates.
        lattice_vectors (array-like): 3×3 lattice‐vector matrix.

    Returns:
        atomic_cartesian_coords (numpy.ndarray): N×3 array of Cartesian coordinates.
    """
    lattice_matrix = np.array(lattice_vectors)
    fractional_coords_matrix = np.array(atomic_fractional_coords)
    atomic_cartesian_coords = np.dot(fractional_coords_matrix, lattice_matrix)
    return atomic_cartesian_coords

def generate_lattice_positions(atomic_cartesian_coords, atomic_numbers, lattice_vectors, lattice_size):
    """Generate positions and atomic numbers for all periodic images in a super‐lattice.

    Args:
        atomic_cartesian_coords (numpy.ndarray): N×3 Cartesian coords of atoms in the unit cell.
        atomic_numbers (array-like): Length‐N sequence of atomic numbers.
        lattice_vectors (array-like): 3×3 matrix of lattice vectors.
        lattice_size (int): Range of supercell images in each direction (±lattice_size).

    Returns:
        all_positions (numpy.ndarray): (N*(2*lattice_size+1)^3)×3 array of all Cartesian positions.
        all_numbers (numpy.ndarray): Corresponding array of atomic numbers.
    """
    range_of_lattice = np.arange(-lattice_size, lattice_size + 1)
    grid_points = np.array(np.meshgrid(range_of_lattice, range_of_lattice, range_of_lattice)).transpose().reshape(-1, 3)
    lattice_vectors_np = np.array(lattice_vectors)
    lattice_points = np.dot(grid_points, lattice_vectors_np)
    all_positions = []
    all_numbers = []
    for i, coord in enumerate(atomic_cartesian_coords):
        current_positions = lattice_points + coord
        all_positions.extend(current_positions.tolist())
        all_numbers.extend([atomic_numbers[i]] * len(current_positions))
    return np.array(all_positions), np.array(all_numbers)

def calculate_vectors_from_point(ref_point, lattice, lattice_atomic_numbers):
    """Compute vectors and distances from a reference point to every lattice point.

    Args:
        ref_point (array-like): 3‐element Cartesian reference point.
        lattice (numpy.ndarray): M×3 array of lattice point positions.
        lattice_atomic_numbers (array-like): Length‐M sequence of atomic numbers at those points.

    Returns:
        vectors (numpy.ndarray): M×3 array of displacement vectors (lattice – ref_point).
        magnitudes (numpy.ndarray): Length‐M array of vector norms.
        atomic_numbers (numpy.ndarray): Copy of input atomic numbers array.
    """
    lattice_array = np.array(lattice)
    vectors = lattice_array - ref_point
    magnitudes = np.linalg.norm(vectors, axis=1)
    atomic_numbers = np.array(lattice_atomic_numbers)
    return vectors, magnitudes, atomic_numbers


def dict_to_matrix(relations):
    """
    Build a 3×3 sympy matrix from a dict mapping symbols Jᵢⱼ to expressions.

    Args:
        relations (dict): Mapping from sympy.Symbol J{i}{j} to sympy expressions.

    Returns:
        matrix (list[list[sympy.Expr]]): 3×3 nested list representing the tensor matrix.
    """
    matrix = [[None for _ in range(3)] for _ in range(3)]
    for i in range(1, 4):
        for j in range(1, 4):
            key_sym = sp.symbols(f"J{i}{j}")
            matrix[i-1][j-1] = relations.get(key_sym, key_sym)
    return matrix

def endpoints_equivalent_with_common_translation(r1p, r1, r2p, r2):
    """
    Check if (r1′,r2′) = (r1,r2) up to the same integer translation.

    Args:
        r1p (array-like): Transformed first endpoint.
        r1 (array-like): Original first endpoint.
        r2p (array-like): Transformed second endpoint.
        r2 (array-like): Original second endpoint.

    Returns:
        bool: True if both endpoints match up to one common lattice lattice translation.
    """
    diff1 = r1p - r1
    diff2 = r2p - r2
    T1 = np.round(diff1)
    T2 = np.round(diff2)
    if not (np.allclose(r1p, r1 + T1, atol=tolerance) and np.allclose(r2p, r2 + T2, atol=tolerance)):
        return False
    return np.allclose(T1, T2, atol=tolerance)

def endpoints_equivalent_with_common_translation_inverted(r1p, r1, r2p, r2):
    """
    Check if (r1′,r2′) = (r2,r1) (inverted) up to the same integer lattice translation.

    Args:
        r1p (array-like): Transformed first endpoint.
        r1 (array-like): Original first endpoint.
        r2p (array-like): Transformed second endpoint.
        r2 (array-like): Original second endpoint.

    Returns:
        bool: True if the pair is inverted and matches by a common lattice translation.
    """
    diff1 = r1p - r2
    diff2 = r2p - r1
    T1 = np.round(diff1)
    T2 = np.round(diff2)
    if not (np.allclose(r1p, r2 + T1, atol=tolerance) and np.allclose(r2p, r1 + T2, atol=tolerance)):
        return False
    return np.allclose(T1, T2, atol=tolerance)

def find_specialising_rotations(r1, r2, rotations, translations):
    """
    Identify rotations that either preserve or invert a bond under symmetry.

    Args:
        r1 (array-like): Fractional coords of the first endpoint.
        r2 (array-like): Fractional coords of the second endpoint.
        rotations (array-like): List/array of 3×3 rotation matrices.
        translations (array-like): Corresponding list/array of 3‐element translation vectors.

    Returns:
        preserving_rotations (list[numpy.ndarray]): Ops mapping (r1,r2)→(r1,r2).
        inverting_rotations (list[numpy.ndarray]): Ops mapping (r1,r2)→(r2,r1).
    """
    preserving_rotations = []
    inverting_rotations = []
    for R, t in zip(rotations, translations):
        r1p = np.dot(R, r1) + t
        r2p = np.dot(R, r2) + t
        if endpoints_equivalent_with_common_translation(r1p, r1, r2p, r2):
            preserving_rotations.append(R)
        elif endpoints_equivalent_with_common_translation_inverted(r1p, r1, r2p, r2):
            inverting_rotations.append(R)
    return preserving_rotations, inverting_rotations

def specialize_tensor(J, preserving_rotations, inverting_rotations):
    """
    Impose symmetry constraints on a general 3×3 tensor J.

    For preserving ops: J = R·J·Rᵀ; for inverting ops: J = (R·J·Rᵀ)ᵀ.

    Args:
        J (sympy.Matrix): 3×3 symbolic tensor.
        preserving_rotations (list[numpy.ndarray]): Rotations that preserve the bond.
        inverting_rotations (list[numpy.ndarray]): Rotations that invert the bond.

    Returns:
        solutions (list[dict]): Solution dict(s) mapping tensor symbols to relations, if any.
    """
    eqs = []
    for R in preserving_rotations:
        R_sym = sp.Matrix(R)
        constraint = R_sym * J * R_sym.T
        for i in range(3):
            for j in range(3):
                eqs.append(sp.Eq(constraint[i, j], J[i, j]))
    for R in inverting_rotations:
        R_sym = sp.Matrix(R)
        constraint = (R_sym * J * R_sym.T).T
        for i in range(3):
            for j in range(3):
                eqs.append(sp.Eq(constraint[i, j], J[i, j]))
    variables_J = list(J)
    solutions = sp.solve(eqs, variables_J, dict=True, simplify=True)
    return solutions


def dict_to_matrix_bond(J, solution):
    """
    Substitute a symmetry solution into J and simplify to get the specialized tensor.

    Args:
        J (sympy.Matrix): 3×3 symbolic tensor.
        solution (dict): Mapping from symbols to expressions (one result of specialize_tensor).

    Returns:
        sympy.Matrix: Simplified 3×3 tensor after applying symmetry relations.
    """
    return sp.simplify(J.subs(solution))


def group_bonds_by_magnitude_and_symmetry(bond_dicts, site_sym_ops, atom_cart):
    """
    Groups the bonds for an atom (each provided as a dictionary with a 'bond_vector_frac' and 'bond_index')
    first by unique bond magnitude (within tolerance), computed in Cartesian space, and then further groups bonds that
    are symmetry-equivalent under the site's symmetry operations.
   
    Args:
        bond_dicts (list[dict]): Each has 'bond_vector_frac', 'bond_endpoint_cartesian', etc.
        site_sym_ops (list[numpy.ndarray]): Rotation matrices of the site symmetry group.
        atom_cart (array-like): 3‐element Cartesian coordinates of the central atom.

    Returns:
        list[tuple]:
            Each tuple is (bond_length (float), groups (list[list[int]])),
            where each inner list lists bond indices that are symmetry‐equivalent.
    """
    # Compute bond magnitude in Cartesian coordinates for each bond.
    for bond in bond_dicts:
        bond['bond_magnitude'] = np.linalg.norm(bond['bond_endpoint_cartesian'] - atom_cart)
   
    # Group bonds by bond magnitude.
    mag_groups = {}
    for bond in bond_dicts:
        mag = bond['bond_magnitude']
        found_key = None
        for key in mag_groups:
            if np.isclose(mag, key, atol=tolerance):
                found_key = key
                break
        if found_key is None:
            mag_groups[mag] = [bond]
        else:
            mag_groups[found_key].append(bond)
   
    # Sort groups by magnitude.
    mag_groups_list = sorted(mag_groups.items(), key=lambda x: x[0])
   
    grouped_output = []
    for mag, bonds in mag_groups_list:
        vecs = [bond['bond_vector_frac'] for bond in bonds]
        indices = [bond['bond_index'] for bond in bonds]
        classified = [False] * len(vecs)
        symmetry_groups = []
        for i in range(len(vecs)):
            if not classified[i]:
                rep = vecs[i]
                group = [indices[i]]
                classified[i] = True
                for j in range(i+1, len(vecs)):
                    if not classified[j]:
                        candidate = vecs[j]
                        equivalent = False
                        for op in site_sym_ops:
                            transformed = np.dot(op, rep)
                            if np.allclose(transformed, candidate, atol=tolerance):
                                equivalent = True
                                break
                        if equivalent:
                            group.append(indices[j])
                            classified[j] = True
                symmetry_groups.append(group)
        grouped_output.append((mag, symmetry_groups))
    return grouped_output


# Functions for finding point group of site


def determine_candidate_symmetry_operations(symmetry_dataset, atomic_fractional_coord):
    """
    Find all space‐group ops that leave a given site invariant (up to lattice translation).

    Args:
        symmetry_dataset (dict): Output of spglib.get_symmetry_dataset.
        atomic_fractional_coord (array-like): 3‐element fractional coord of the atom.

    Returns:
        candidate_symmetry_operations (list[tuple]):
            Each tuple is (rotation (3×3 numpy.ndarray), translation (3‐element array)).
    
    Reference:
        Grosse-Kunstleve, R. W., & Adams, P. D. (2002). Algorithms for deriving crystallographic space-group information.
        II. Treatment of special positions. Acta Crystallographica Section A, 58(1), 60-65.
    """
    rotations = symmetry_dataset.rotations
    translations = symmetry_dataset.translations
    candidate_symmetry_operations = []
    for rotation, translation in zip(rotations, translations):
        transformed_atomic_fractional_coord = np.dot(rotation, atomic_fractional_coord) + translation
        displacement_vector = transformed_atomic_fractional_coord - atomic_fractional_coord
        normalised_displacement_vector = np.mod(displacement_vector + 0.5, 1) - 0.5
        distance_squared = np.sum(np.square(normalised_displacement_vector))
        if distance_squared < np.square(tolerance):
            candidate_symmetry_operations.append((rotation, translation))
    return candidate_symmetry_operations

def sort_candidate_operations_by_squared_distance(candidate_symmetry_operations, atomic_fractional_coord):
    """
    Sort site‐symmetry ops by how closely they map the site onto itself.

    Args:
        candidate_symmetry_operations (list[tuple]): (rotation, translation) pairs.
        atomic_fractional_coord (array-like): 3‐element fractional coord of the atom.

    Returns:
        candidate_symmetry_operations (list[tuple]): Same format as input, but sorted by minimal displacement².
        
    Reference:
        Grosse-Kunstleve, R. W., & Adams, P. D. (2002). Algorithms for deriving crystallographic space-group information.
        II. Treatment of special positions. Acta Crystallographica Section A, 58(1), 60-65.
    """
    candidate_symmetry_operations.sort(
        key=lambda op: np.sum(np.square(np.mod(np.dot(op[0], atomic_fractional_coord) + op[1] - atomic_fractional_coord + 0.5, 1) - 0.5))
    )
    return candidate_symmetry_operations

def generate_symmetry_group_by_multiplication(candidate_symmetry_operations):
    """
    Build the point‐group via closure under multiplication of candidate ops to ensure a full point group.

    Args:
        candidate_symmetry_operations (list[tuple]): (rotation, translation) pairs for the site.

    Returns:
        point_group_name (str): Hermann–Mauguin symbol from spglib.
        final_symmetry_operations (numpy.ndarray): Array of unique rotation matrices.
    
    Reference:
        Grosse-Kunstleve, R. W., & Adams, P. D. (2002). Algorithms for deriving crystallographic space-group information.
        II. Treatment of special positions. Acta Crystallographica Section A, 58(1), 60-65.
    """
    final_symmetry_operations = []
    for rotation_1, translation_1 in candidate_symmetry_operations:
        for rotation_2, translation_2 in candidate_symmetry_operations:
            resultant_rotation = np.dot(rotation_1, rotation_2)
            if not any(np.allclose(resultant_rotation, rot, atol=tolerance) for rot in final_symmetry_operations):
                final_symmetry_operations.append(resultant_rotation)
    final_symmetry_operations = np.array(final_symmetry_operations)
    final_symmetry_operations_point_group_name = spglib.get_pointgroup(final_symmetry_operations)[0]
    return final_symmetry_operations_point_group_name, final_symmetry_operations


# Functions for processing and printing specialized tensors


def process_atoms(lattice_vectors, atomic_fractional_coords, atomic_numbers, lattice_size):
    """
    Compute per‐atom symmetry data, bond lists, and specialized exchange tensors.

    Args:
        lattice_vectors (array-like): 3×3 lattice matrix.
        atomic_fractional_coords (array-like): N×3 fractional coords.
        atomic_numbers (array-like): Length‐N atomic numbers.
        lattice_size (int): Supercell range (± lattice_size) for bond search.

    Returns:
        dict: Comprehensive data including
            - 'spacegroup' (str)
            - 'lattice_vectors' (numpy.ndarray)
            - 'atomic_fractional_coords' (numpy.ndarray)
            - 'atoms' (dict): atom‐index → subdict with site symmetry, bond tensors, etc.
    """
    lattice_vectors_np = convert_lattice_vectors_to_NumPy_array(lattice_vectors)
    atomic_fractional_coords_np = convert_atomic_fractional_coords_to_NumPy_array(atomic_fractional_coords)
    atomic_cartesian_coords = convert_fractional_coordinates_to_cartesian(atomic_fractional_coords, lattice_vectors)
   
    cell = (lattice_vectors_np, atomic_fractional_coords_np, np.array(atomic_numbers))
    symmetry_dataset = spglib.get_symmetry_dataset(cell, symprec=tolerance)
    spacegroup = symmetry_dataset.international
   
    supercell_positions_cartesian, supercell_atomic_numbers = generate_lattice_positions(
        atomic_cartesian_coords, atomic_numbers, lattice_vectors, lattice_size)
   
    results = {}
    inv_lattice = np.linalg.inv(lattice_vectors_np)
   
    for atom_index, (atom_frac, atom_cart, atom_num) in enumerate(zip(atomic_fractional_coords, atomic_cartesian_coords, atomic_numbers)):
        atom_data = {}
        atom_data["atomic_number"] = atom_num
        atom_data["atom_fractional_coord"] = atom_frac
        atom_data["atom_cartesian_coord"] = atom_cart
       
        # Compute site symmetry group for the atom.
        candidate_sym_ops = determine_candidate_symmetry_operations(symmetry_dataset, atom_frac)
        candidate_sym_ops = sort_candidate_operations_by_squared_distance(candidate_sym_ops, atom_frac)
        site_sym_name, site_sym_ops = generate_symmetry_group_by_multiplication(candidate_sym_ops)
        atom_data["site_symmetry_group_name"] = site_sym_name
        atom_data["site_symmetry_group_operations"] = site_sym_ops
       
        # Calculate bond vectors from the atom.
        bonds_cartesian, bonds_magnitudes, bond_atomic_numbers = calculate_vectors_from_point(atom_cart, supercell_positions_cartesian, supercell_atomic_numbers)
        atom_data["bond_vectors_cartesian"] = bonds_cartesian
        atom_data["bond_magnitudes"] = bonds_magnitudes
        atom_data["bond_endpoint_atomic_numbers"] = bond_atomic_numbers
       
        # Specializing exchange tensor using bond endpoints
        bond_specialized_tensors = []
        for i, bond_vector in enumerate(bonds_cartesian):
            endpoint_cartesian = atom_cart + bond_vector
            # Retain the full endpoint (no modulo applied).
            endpoint_fractional = np.dot(endpoint_cartesian, inv_lattice)
            r1 = np.array(atom_frac)
            r2 = np.array(endpoint_fractional)
            preserving_rot, inverting_rot = find_specialising_rotations(r1, r2, symmetry_dataset.rotations, symmetry_dataset.translations)
           
            # Define a general symbolic 3x3 exchange tensor J.
            Jxx, Jxy, Jxz, Jyx, Jyy, Jyz, Jzx, Jzy, Jzz = sp.symbols('Jxx Jxy Jxz Jyx Jyy Jyz Jzx Jzy Jzz')
            J = sp.Matrix([[Jxx, Jxy, Jxz],
                           [Jyx, Jyy, Jyz],
                           [Jzx, Jzy, Jzz]])
            solutions = specialize_tensor(J, preserving_rot, inverting_rot)
            if solutions:
                specialized_tensor = dict_to_matrix_bond(J, solutions[0])
            else:
                specialized_tensor = J
            # Compute the bond vector in fractional coordinates without normalization.
            bond_vector_frac = endpoint_fractional - atom_frac
            bond_specialized_tensors.append({
                "bond_index": i,
                "bond_endpoint_cartesian": endpoint_cartesian,
                "bond_endpoint_fractional": endpoint_fractional,
                "bond_vector_frac": bond_vector_frac,
                "bond_endpoint_atomic_number": bond_atomic_numbers[i],
                "preserving_rotations": preserving_rot,
                "inverting_rotations": inverting_rot,
                "specialized_tensor_solution": solutions[0] if solutions else {},
                "specialized_tensor_matrix": specialized_tensor
            })
        atom_data["bond_specialized_tensors"] = bond_specialized_tensors
       
        # Group bonds by magnitude (Cartesian) and symmetry equivalence using the site symmetry operations.
        grouped_bonds = group_bonds_by_magnitude_and_symmetry(bond_specialized_tensors, site_sym_ops, atom_cart)
        atom_data["grouped_bond_specialized_tensors"] = grouped_bonds
       
        results[atom_index] = atom_data
   
    overall_data = {
        "spacegroup": spacegroup,
        "lattice_vectors": lattice_vectors_np,
        "atomic_fractional_coords": atomic_fractional_coords_np,
        "atoms": results
    }
    return overall_data



def print_grouped_bond_tensors(atom_index, overall_data):
    """
    Pretty‐print specialized bond tensors for a given atom, grouped by length and symmetry.

    Args:
        atom_index (int): Index of the atom within overall_data['atoms'].
        overall_data (dict): Output from process_atoms().

    Returns:
        None
    """
    atom_data = overall_data["atoms"][atom_index]
    print(f"Atom {atom_index}: Atomic number: {atom_data['atomic_number']}, Cartesian coordinates: {atom_data['atom_cartesian_coord']}")
    print("Site symmetry group:", atom_data["site_symmetry_group_name"])
    groups = atom_data["grouped_bond_specialized_tensors"]
    for mag, sym_groups in groups:
        print(f"\nBond magnitude (Cartesian): {mag:.4f}")
        for grp in sym_groups:
            print("  Symmetry-equivalent bond indices:")
            for idx in grp:
                bond_info = atom_data["bond_specialized_tensors"][idx]
                print(f"    Bond index: {idx}")
                print(f"    Bond endpoint atomic number: {bond_info['bond_endpoint_atomic_number']}")
                print(f"    Bond endpoint fractional coordinates: {bond_info['bond_endpoint_fractional']}")
                print(f"    Bond endpoint Cartesian coordinates: {bond_info['bond_endpoint_cartesian']}")
                print("    Specialized tensor:")
                sp.pretty_print(bond_info["specialized_tensor_matrix"])
                print("    ---")


# Main section

# lattice_vectors, atomic_fractional_coords, atomic_numbers = extract_cif_data(r"C:\Users\georg\OneDrive - University of Leeds\MPhys Project\Code\Materials\MnF2.cif")
lattice_vectors, atomic_fractional_coords, atomic_numbers = extract_cif_data(r"C:\Users\georg\OneDrive - University of Leeds\MPhys Project\Code\Materials\Y3Fe5O12_symmetrised.cif")

# Example Diamond lattice (Lax)
# lattice_vectors = np.array([[1,0,0],
#                             [0,1,0],
#                             [0,0,1]])

# atomic_numbers = (1,1,1,1,1,1,1,1)

# atomic_fractional_coords = np.array([[0,0,0],
#                                      [0.5,0.5,0],
#                                      [0.5,0,0.5],
#                                      [0,0.5,0.5],
#                                      [0.25,0.25,0.25],
#                                      [0.75,0.75,0.25],
#                                      [0.75,0.25,0.75],
#                                      [0.25,0.75,0.75]])

# lattice_vectors = np.array([[1,0,0],
#                             [0,1,0],
#                             [0,0,1]])

# atomic_numbers = (1,1)

# atomic_fractional_coords = np.array([[0,0,0],
#                                      [0.5,0.5,0.5]])

overall_data = process_atoms(lattice_vectors, atomic_fractional_coords, atomic_numbers, lattice_size=0)

i = 62
print(f"Specialized bond tensors for bonds connected to atom {i}:")
print_grouped_bond_tensors(i, overall_data)

# Note: The bond which 'connects the atom to itself' and which has zero magnitude is also given (i.e. the self interaction tensor, I don't think this has any physical meaning)