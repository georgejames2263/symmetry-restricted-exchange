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
