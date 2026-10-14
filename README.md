# Molecular_Structure_Toolkit

A Streamlit-based cheminformatics application for molecular visualization, descriptor calculation, structure-file conversion, chemical-name lookup, and 3D conformer analysis.

The application combines **RDKit**, **Open Babel**, **PubChem**, **py3Dmol**, and **Streamlit** in a single graphical interface.

## Weblink

https://molecularstructuretoolkit-fbxpsgpepzgz7go6udznll.streamlit.app/

---

## Features

The application is organized into four main panels:

### 1. 2D Structure

Enter one or more SMILES strings in the shared SMILES workspace.

For each valid molecule, the application provides:

- input SMILES
- canonical isomeric SMILES
- publication-quality 2D molecular structure
- adjustable atom-label font size
- adjustable bond thickness
- adjustable image dimensions
- high-resolution image export

Calculated molecular descriptors include:

- molecular formula
- molecular weight
- exact mass
- LogP
- TPSA
- hydrogen-bond donors
- hydrogen-bond acceptors
- rotatable bonds
- ring count
- heavy-atom count
- formal charge

Available image downloads:

- SVG
- PNG
- JPG
- high-resolution PNG

The descriptors are calculated locally with **RDKit**.

### 2D drawing settings

Tune your parameter of choice. 
For publication-grade image: Atom label font size = 45, Bond thickness = 5, Image width = 600, Image height = 600, High-res PNG = 3x

---

### 2. 3D Conformers

The 3D workflow is separated into two stages.

#### Stage 1 — Initial ETKDGv3 structure

As soon as a valid SMILES molecule is available, the application generates a single 3D structure using:

- RDKit ETKDGv3
- explicit hydrogens
- stereochemistry preservation
- small-ring torsion handling
- macrocycle torsion handling

This initial structure is **not force-field minimized**.

It can be viewed immediately and downloaded without running a conformational search.

Available ETKDGv3 downloads:

- XYZ
- MOL
- SDF
- PDB
- MOL2

#### Stage 2 — Optional conformer search

The user may optionally perform a conformational search.

The force field is chosen explicitly:

- MMFF94s
- UFF

The conformer-search workflow is:

```text
SMILES
  ↓
ETKDGv3 multiple-conformer generation
  ↓
embedding RMSD pruning
  ↓
MMFF94s or UFF optimization
  ↓
energy ranking
  ↓
RMSD diversity filtering
  ↓
lowest-energy conformer
```

User-adjustable settings include:

- number of conformers
- ETKDG embedding RMSD pruning threshold
- final RMSD diversity threshold
- maximum number of retained conformers
- MMFF94s or UFF force field

The application displays:

- number of conformers generated
- number successfully optimized
- number retained after RMSD filtering
- selected force field
- lowest-energy conformer
- conformer ID
- energy
- full conformer energy table
- relative energy, ΔE
- RMSD-diverse conformer table

The lowest-energy conformer can be downloaded as:

- XYZ
- MOL
- SDF
- PDB
- MOL2
- interactive 3D HTML

> MMFF94s and UFF energies should only be compared within the same force field. Absolute energies from different force fields are not directly comparable.

---

### 3. File → Structure

The application currently accepts:

- MOL
- MOL2
- SDF

PDB input is intentionally disabled.

Small-molecule PDB files often do not contain reliable bond-order information, which can lead to incorrect reconstruction of:

- C=O bonds
- C=C bonds
- aromatic systems
- formal charges
- protonation states

For reliable structure reconstruction, MOL, MOL2, or SDF is preferred because these formats normally retain explicit bond-order information.

For uploaded structure files, the application provides:

- canonical isomeric SMILES
- calculated RDKit descriptors
- reconstructed 2D structure
- SVG download
- PNG download
- SMILES `.smi` download
- option to add the generated SMILES to the shared workspace

#### Explicit hydrogen control

The file-analysis panel includes:

```text
Hide explicit hydrogens in 2D structure and SMILES analysis
```

This option is enabled by default.

It removes unnecessary explicit hydrogens from:

- CH3 groups
- CH2 groups
- carbon rings
- ordinary carbon atoms

while preserving chemically meaningful hydrogens where required by valence or aromaticity.

---

### 4. Name → Structure

The application can resolve chemical names through the **PubChem PUG REST** service.

Supported input can include:

- IUPAC names
- systematic chemical names
- common names
- drug names
- other names recognized by PubChem

Multiple names can be entered, one per line.

Workflow:

```text
Chemical name
  ↓
PubChem lookup
  ↓
SMILES
  ↓
RDKit validation
  ↓
canonical isomeric SMILES
  ↓
RDKit descriptors
  ↓
2D structure
```

For each successfully resolved molecule, the application displays:

- PubChem CID
- PubChem title
- PubChem IUPAC name
- SMILES returned by PubChem
- RDKit canonical isomeric SMILES
- calculated molecular descriptors
- reconstructed 2D structure

Available downloads:

- `.smi`
- SVG
- PNG
- JPG
- high-resolution PNG

The resolved SMILES can also be added directly to the shared SMILES workspace.

> The Name → Structure panel requires an internet connection because it queries PubChem.

---

## Shared SMILES Workspace

The application uses a shared SMILES workspace at the top of the interface.

Molecules can enter the workspace from:

- manually entered SMILES
- MOL/MOL2/SDF conversion
- PubChem name lookup

Once a molecule is added to the workspace, the same structure can be used in:

- 2D analysis
- descriptor calculation
- ETKDGv3 3D generation
- conformer search

---

## Main Dependencies

The project uses:

- Python
- Streamlit
- RDKit
- py3Dmol
- NumPy
- pandas
- Pillow
- requests
- Open Babel

---

## Installation

A Conda environment is recommended.

Example:

```bash
conda create -n molecular_structure python=3.12 -y
conda activate molecular_structure
```

Install RDKit:

```bash
conda install -c conda-forge rdkit -y
```

Install the Python dependencies:

```bash
pip install streamlit py3Dmol numpy pandas Pillow requests
```

Install Open Babel:

```bash
conda install -c conda-forge openbabel -y
```

Check Open Babel:

```bash
obabel -V
```

---

## Running the Application

From the project directory:

```bash
conda activate molecular_structure
streamlit run app.py
```

Then open:

```text
http://localhost:8501
```

Streamlit will normally open the application automatically in the default browser.

---

## Suggested Project Structure

```text
molecular_structure_toolkit/
│
├── app.py
├── chemistry_tools.py
├── conformer_tools.py
├── file_tools.py
├── name_tools.py
├── requirements.txt
├── packages.txt
├── README.md
└── .gitignore
```

---

## Module Overview

### `app.py`

Main Streamlit graphical interface.

Responsibilities include:

- shared SMILES workspace
- tab navigation
- 2D visualization
- 3D visualization
- descriptor display
- structure-file upload
- name lookup
- download controls
- session-state management

### `chemistry_tools.py`

Core RDKit/Open Babel chemistry utilities.

Includes functionality for:

- SMILES parsing
- canonical SMILES generation
- molecular descriptors
- SVG drawing
- PNG drawing
- JPG conversion
- XYZ export
- MOL export
- SDF export
- PDB export
- MOL2 export

### `conformer_tools.py`

3D structure and conformational-search utilities.

Includes:

- single ETKDGv3 structure generation
- multiple ETKDGv3 conformer generation
- MMFF94s optimization
- UFF optimization
- conformer energy evaluation
- relative-energy calculation
- RMSD calculation
- RMSD-based diversity filtering
- lowest-energy conformer extraction

### `file_tools.py`

Molecular structure-file import utilities.

Currently supports:

- MOL
- MOL2
- SDF

Also handles:

- RDKit parsing
- Open Babel MOL2 fallback
- structure sanitization
- explicit-hydrogen removal
- canonical SMILES export

### `name_tools.py`

Chemical-name lookup utilities using PubChem PUG REST.

Handles:

- single-name resolution
- multiple-name resolution
- PubChem CID retrieval
- IUPAC name retrieval
- title retrieval
- SMILES retrieval
- per-name error handling

---

## Scientific Notes

### ETKDGv3

ETKDGv3 is used to generate physically reasonable initial 3D geometries from molecular topology.

The first displayed ETKDGv3 structure is intentionally kept separate from force-field optimization.

This makes it possible to:

- inspect the initial RDKit-generated geometry
- download the unoptimized structure
- decide independently whether a conformational search is necessary

### MMFF94s

MMFF94s is generally suitable for many organic and drug-like molecules.

If parameters are unavailable for a molecule, the MMFF94s search may fail.

### UFF

UFF is available as an explicit alternative force field.

It supports a broader range of elements than MMFF94s, although its parameterization and energetic behavior differ.

### Conformer energies

Conformer energies are most useful as **relative energies within one force field**.

Do not directly compare:

```text
MMFF94s energy
```

with:

```text
UFF energy
```

as though they were on the same absolute energy scale.

### PDB input

PDB is not currently accepted as an input structure format because bond orders can be missing or ambiguous.

PDB export from a structure already defined in RDKit is still available in the 3D section.

---

## External Services

### PubChem

The Name → Structure panel uses:

**PubChem PUG REST**

https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

PubChem is used only for chemical-name resolution.

Molecular descriptor calculation and 2D rendering are performed locally with RDKit.

Applications using PubChem programmatically should respect PubChem's request-rate guidance.

---

## Third-Party Software

This project uses several open-source scientific software packages.

### RDKit

Used for:

- molecular graph handling
- SMILES parsing
- molecular descriptors
- 2D depiction
- ETKDGv3
- MMFF94s
- UFF
- conformer analysis

Project:

https://www.rdkit.org/

### Open Babel

Used primarily for MOL2 export and fallback MOL2 conversion.

Project:

https://openbabel.org/

### py3Dmol

Used for interactive 3D molecular visualization.

Project:

https://3dmol.csb.pitt.edu/

### Streamlit

Used for the graphical web interface.

Project:

https://streamlit.io/

---

## Current Scope

The current version focuses on:

- molecular structure handling
- molecular visualization
- classical cheminformatics
- descriptor calculation
- conformational analysis
- chemical-name resolution

No machine-learning molecular-property prediction model is currently used for the descriptor panel.

Properties such as:

- molecular weight
- LogP
- TPSA
- H-bond donors
- H-bond acceptors
- rotatable bonds

are calculated directly from the molecular graph using RDKit.

---

## Possible Future Extensions

Potential future additions include:

- solubility prediction
- pKa prediction
- permeability prediction
- toxicity prediction
- CYP inhibition prediction
- ADMET models
- molecular fingerprints
- similarity searching
- substructure searching
- ligand preparation
- conformer ensemble SDF export
- molecular database integration
- structure-based molecular analysis

These features should remain separate from the current calculated-descriptor section so that calculated descriptors and machine-learning predictions are clearly distinguished.

---

## Disclaimer

This application is intended primarily for research, education, and molecular-data exploration.

Automatically generated structures, conformers, protonation states, stereochemistry, and chemical-name mappings should be verified before use in critical scientific calculations.

In particular:

- name resolution depends on the corresponding PubChem record
- MOL2 interpretation can depend on atom typing and charge conventions
- conformer energies depend on the selected force field
- an ETKDGv3 geometry is an initial generated geometry, not necessarily the global-energy minimum or biologically active conformation

---

## License and Attribution

Before public redistribution, include the appropriate licenses and notices for all third-party software used by the project.

If this application is used in scientific work, cite the relevant software and methods according to the requirements or recommendations of their respective authors.
