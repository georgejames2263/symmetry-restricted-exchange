# symmetry-restricted-exchange

Code developed for my MPhys project to automate symmetry analysis of magnetic exchange interactions in crystals, replacing manual tensor specialization using Neumann’s principle.

---

## Overview

This Python tool determines the symmetry-allowed form of exchange tensors between magnetic atoms in a crystal lattice. By analyzing the local bond symmetry (rather than the full space group), it automates the process of constraining tensor components based on Neumann’s principle — a task typically done manually in theoretical magnetism research.

---

## Features

- Parses CIF files and extracts atomic structure
- Computes site symmetry groups for individual atoms
- Finds all bond vectors and groups symmetry-equivalent bonds
- Applies preserving/inverting symmetry operations to constrain a symbolic 3×3 exchange tensor
- Outputs specialized tensors grouped by bond length and symmetry

---

## Dependencies

You’ll need Python 3.9+ and the following packages:

- `numpy`
- `sympy`
- `spglib`
- `pymatgen`

Install with:

```bash
pip install numpy sympy spglib pymatgen

---

## 🔍 Example

### Input

Suppose you want to analyze the crystal structure in `Materials/Y3Fe5O12_symmetrised.cif` and examine the exchange tensors associated with atom 62 (the 63rd atom in the structure). Modify the bottom of the script like this:

```python
file_path = "Materials/Y3Fe5O12_symmetrised.cif"
atom_index = 62

lattice_vectors, atomic_fractional_coords, atomic_numbers = extract_cif_data(file_path)
overall_data = process_atoms(lattice_vectors, atomic_fractional_coords, atomic_numbers, lattice_size=0)

print(f"Specialized bond tensors for bonds connected to atom {atom_index}:")
print_grouped_bond_tensors(atom_index, overall_data)

Atom 62: Atomic number: 26, Cartesian coordinates: [3.216 2.955 1.750]
Site symmetry group: m-3m

Bond magnitude (Cartesian): 2.9821
  Symmetry-equivalent bond indices:
    Bond index: 15
    Bond endpoint atomic number: 26
    Bond endpoint fractional coordinates: [0.25 0.5  0.25]
    Bond endpoint Cartesian coordinates: [5.865 2.955 1.750]
    Specialized tensor:
⎡ Jxx   0    0  ⎤
⎢  0   Jxx   0  ⎥
⎣  0    0   Jxx ⎦
    ---

Bond magnitude (Cartesian): 3.1004
  Symmetry-equivalent bond indices:
    Bond index: 23
    Bond endpoint atomic number: 8
    Bond endpoint fractional coordinates: [0.13 0.61 0.36]
    Bond endpoint Cartesian coordinates: [3.825 3.318 2.465]
    Specialized tensor:
⎡ Jxx   0   Jxz ⎤
⎢  0   Jyy   0  ⎥
⎣ Jzx   0   Jzz ⎦
    ---

