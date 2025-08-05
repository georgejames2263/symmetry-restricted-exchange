# symmetry-restricted-exchange

Code developed for my MPhys project to automate symmetry analysis of magnetic exchange interactions in crystals, replacing manual tensor specialization using Neumann’s principle.

---

## 🧠 Overview

This Python tool determines the symmetry-allowed form of exchange tensors between magnetic atoms in a crystal lattice. By analyzing the local bond symmetry (rather than the full space group), it automates the process of constraining tensor components based on Neumann’s principle — a task typically done manually in theoretical magnetism research.

---

## ⚙️ Features

- Parses CIF files and extracts atomic structure
- Computes site symmetry groups for individual atoms
- Finds all bond vectors and groups symmetry-equivalent bonds
- Applies preserving/inverting symmetry operations to constrain a symbolic 3×3 exchange tensor
- Outputs specialized tensors grouped by bond length and symmetry

---

## 📦 Dependencies

You’ll need Python 3.9+ and the following packages:

- `numpy`
- `sympy`
- `spglib`
- `pymatgen`

Install with:

```bash
pip install numpy sympy spglib pymatgen
```
---

## 🔍 Example

### Input

Suppose you want to analyze the crystal structure in `Materials/Mn5Ge3.cif` and examine the exchange tensors associated with atom 0 (the 1st atom in the .cif file). Modify the bottom of the script like this:

```python
file_path = "Materials/Mn5Ge3.cif"
atom_index = 0

lattice_vectors, atomic_fractional_coords, atomic_numbers = extract_cif_data(file_path)
overall_data = process_atoms(lattice_vectors, atomic_fractional_coords, atomic_numbers, lattice_size=1)

print(f"Specialized bond tensors for bonds connected to atom {atom_index}:")
print_grouped_bond_tensors(atom_index, overall_data)
```
### Output (partial)
```yaml
Specialized bond tensors for bonds connected to atom 0:
Atom 0: Atomic number: 25, Cartesian coordinates: [0.86905913 1.50525457 3.74845551]
Site symmetry group: mm2

Bond magnitude (Cartesian): 0.0000
  Symmetry-equivalent bond indices:
    Bond index: 13
    Bond endpoint atomic number: 25
    Bond endpoint fractional coordinates: [0.24141583 0.24141583 0.75      ]
    Bond endpoint Cartesian coordinates: [0.86905913 1.50525457 3.74845551]
    Specialized tensor:
⎡Jyy  Jyx   0 ⎤
⎢             ⎥
⎢Jyx  Jyy   0 ⎥
⎢             ⎥
⎣ 0    0   Jzz⎦
    ---

Bond magnitude (Cartesian): 2.4799
  Symmetry-equivalent bond indices:
    Bond index: 364
    Bond endpoint atomic number: 32
    Bond endpoint fractional coordinates: [2.77555756e-17 3.94443570e-01 7.50000000e-01]
    Bond endpoint Cartesian coordinates: [-1.419935    2.45939956  3.74845551]
    Specialized tensor:
⎡Jxx  Jxy   0 ⎤
⎢             ⎥
⎢Jyx  Jyy   0 ⎥
⎢             ⎥
⎣ 0    0   Jzz⎦
    ---
    Bond index: 418
    Bond endpoint atomic number: 32
    Bond endpoint fractional coordinates: [0.39444357 0.         0.75      ]
    Bond endpoint Cartesian coordinates: [2.83986999 0.         3.74845551]
    Specialized tensor:
⎡Jxx  Jxy   0 ⎤
⎢             ⎥
⎢Jyx  Jyy   0 ⎥
⎢             ⎥
⎣ 0    0   Jzz⎦
    ---

Bond magnitude (Cartesian): 2.6217
  Symmetry-equivalent bond indices:
    Bond index: 310
    Bond endpoint atomic number: 32
    Bond endpoint fractional coordinates: [0.60555643 0.60555643 0.75      ]
    Bond endpoint Cartesian coordinates: [2.17990818 3.77571173 3.74845551]
    Specialized tensor:
⎡Jyy  Jyx   0 ⎤
⎢             ⎥
⎢Jyx  Jyy   0 ⎥
⎢             ⎥
⎣ 0    0   Jzz⎦
    ---

Bond magnitude (Cartesian): 2.7311
  Symmetry-equivalent bond indices:
    Bond index: 391
    Bond endpoint atomic number: 32
    Bond endpoint fractional coordinates: [0.39444357 0.39444357 0.25      ]
    Bond endpoint Cartesian coordinates: [1.419935   2.45939956 1.24948517]
    Specialized tensor:
⎡Jyy  Jyx  Jyz⎤
⎢             ⎥
⎢Jyx  Jyy  Jyz⎥
⎢             ⎥
⎣Jzy  Jzy  Jzz⎦
    ---
    Bond index: 400
    Bond endpoint atomic number: 32
    Bond endpoint fractional coordinates: [0.39444357 0.39444357 1.25      ]
    Bond endpoint Cartesian coordinates: [1.419935   2.45939956 6.24742585]
    Specialized tensor:
⎡Jyy  Jyx  Jyz⎤
⎢             ⎥
⎢Jyx  Jyy  Jyz⎥
⎢             ⎥
⎣Jzy  Jzy  Jzz⎦
    ---

Bond magnitude (Cartesian): 3.0105
  Symmetry-equivalent bond indices:
    Bond index: 39
    Bond endpoint atomic number: 25
    Bond endpoint fractional coordinates: [ 0.         -0.24141583  0.75      ]
    Bond endpoint Cartesian coordinates: [ 0.86905913 -1.50525457  3.74845551]
    Specialized tensor:
⎡Jxx  -Jyx + Jyy   0 ⎤
⎢                    ⎥
⎢Jyx     Jyy       0 ⎥
⎢                    ⎥
⎣ 0       0       Jzz⎦
    ---
    Bond index: 64
    Bond endpoint atomic number: 25
    Bond endpoint fractional coordinates: [-0.24141583  0.          0.75      ]
    Bond endpoint Cartesian coordinates: [-1.73811826  0.          3.74845551]
    Specialized tensor:
⎡Jxy + Jyx  Jxy   0 ⎤
⎢                   ⎥
⎢   Jyx     Jyy   0 ⎥
⎢                   ⎥
⎣    0       0   Jzz⎦
    ---
```

### 💬 Output Explanation

Each block of output corresponds to a group of bonds connected to the selected atom. Bonds are grouped first by their **length (magnitude)** and then by **symmetry equivalence** under the site's local point group. For each symmetry-equivalent bond, the output includes the atomic number and coordinates of the bonded atom, as well as the **specialized exchange tensor** that describes allowed magnetic interactions along that bond.

The tensor is expressed symbolically using SymPy, and its structure reflects the constraints imposed by the bond’s symmetry. For example, diagonal-only tensors indicate isotropic exchange, while off-diagonal terms may correspond to anisotropic or antisymmetric (e.g. DMI-like) interactions.

⚠️ **Note:** The bond group with magnitude `0.0000` represents a self-bond (i.e. the atom bonded to itself) and does **not** have physical meaning in this context — it appears due to how the code includes all nearby lattice points, including the origin.



