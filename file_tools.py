from rdkit import Chem

import io
import os
import subprocess
import tempfile


SUPPORTED_EXTENSIONS = ["mol2", "mol", "sdf"]


def _sanitize_mol(mol):
    if mol is None:
        return None, "RDKit could not parse this structure."

    try:
        Chem.SanitizeMol(mol)
        return mol, None
    except Exception as exc:
        return mol, (
            "The structure was loaded, but RDKit sanitization was incomplete. "
            "Please verify bond orders, aromaticity, protonation, formal "
            f"charges and valence. Details: {exc}"
        )


def _load_mol(text):
    mol = Chem.MolFromMolBlock(
        text,
        sanitize=True,
        removeHs=False,
        strictParsing=False,
    )

    if mol is not None:
        return [(mol, None)]

    mol = Chem.MolFromMolBlock(
        text,
        sanitize=False,
        removeHs=False,
        strictParsing=False,
    )

    mol, warning = _sanitize_mol(mol)
    return [(mol, warning)]


def _load_sdf(data):
    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(data),
        sanitize=True,
        removeHs=False,
        strictParsing=False,
    )

    items = []

    for i, mol in enumerate(supplier, start=1):
        if mol is None:
            items.append((None, f"SDF record {i} could not be parsed."))
        else:
            items.append((mol, None))

    return items


def _mol2_via_openbabel(text):
    with tempfile.TemporaryDirectory() as temp_dir:
        mol2_path = os.path.join(temp_dir, "input.mol2")
        sdf_path = os.path.join(temp_dir, "output.sdf")

        with open(mol2_path, "w", encoding="utf-8") as handle:
            handle.write(text)

        result = subprocess.run(
            [
                "obabel",
                "-imol2",
                mol2_path,
                "-osdf",
                "-O",
                sdf_path,
            ],
            capture_output=True,
            text=True,
        )

        if result.returncode != 0 or not os.path.exists(sdf_path):
            raise RuntimeError(
                "Open Babel MOL2 conversion failed: "
                + (result.stderr or result.stdout or "unknown error")
            )

        with open(sdf_path, "rb") as handle:
            sdf_bytes = handle.read()

    supplier = Chem.ForwardSDMolSupplier(
        io.BytesIO(sdf_bytes),
        sanitize=True,
        removeHs=False,
        strictParsing=False,
    )

    return [mol for mol in supplier if mol is not None]


def _load_mol2(text):
    mol = Chem.MolFromMol2Block(
        text,
        sanitize=True,
        removeHs=False,
    )

    if mol is not None:
        return [(mol, None)]

    try:
        mols = _mol2_via_openbabel(text)

        if mols:
            warning = (
                "The MOL2 file was parsed through Open Babel because direct "
                "RDKit MOL2 parsing failed. Verify SYBYL atom types, "
                "aromaticity, formal charges and protonation."
            )
            return [(mols[0], warning)]

    except Exception as exc:
        return [
            (
                None,
                "MOL2 parsing failed in both RDKit and Open Babel: "
                f"{exc}",
            )
        ]

    return [(None, "MOL2 file could not be parsed.")]


def load_molecules_from_bytes(filename, data):
    """
    Parse MOL2, MOL, or SDF bytes.

    PDB is intentionally not supported in this version because PDB files
    commonly lack reliable bond-order information for small molecules.
    """
    if "." not in filename:
        raise ValueError("The uploaded file has no extension.")

    extension = filename.rsplit(".", 1)[-1].lower()

    if extension not in SUPPORTED_EXTENSIONS:
        raise ValueError(
            f"Unsupported file type '.{extension}'. "
            f"Supported types: {', '.join(SUPPORTED_EXTENSIONS)}"
        )

    if extension == "sdf":
        loaded = _load_sdf(data)
    else:
        text = data.decode("utf-8", errors="replace")

        if extension == "mol":
            loaded = _load_mol(text)
        elif extension == "mol2":
            loaded = _load_mol2(text)
        else:
            loaded = []

    records = []

    for record_index, (mol, warning) in enumerate(loaded, start=1):
        records.append(
            {
                "source_name": filename,
                "record_index": record_index,
                "format": extension.upper(),
                "mol": mol,
                "warning": warning,
            }
        )

    return records


def prepare_mol_for_analysis(mol, remove_explicit_h=True):
    """
    Create a cleaned copy for SMILES, descriptors and 2D depiction.

    When remove_explicit_h=True, explicit hydrogens that can be represented
    implicitly are removed. This hides H atoms on CH3/CH2/rings and produces
    a conventional publication-style 2D drawing.

    Chemically required hydrogens such as [nH] and protonated nitrogens remain
    chemically encoded where appropriate.
    """
    if mol is None:
        return None

    clean_mol = Chem.Mol(mol)

    if remove_explicit_h:
        try:
            clean_mol = Chem.RemoveHs(
                clean_mol,
                sanitize=True,
            )
        except Exception:
            try:
                clean_mol = Chem.RemoveHs(
                    clean_mol,
                    sanitize=False,
                )
                Chem.SanitizeMol(clean_mol)
            except Exception:
                clean_mol = Chem.Mol(mol)

    try:
        Chem.AssignStereochemistry(
            clean_mol,
            cleanIt=True,
            force=True,
        )
    except Exception:
        pass

    return clean_mol


def mol_to_smiles_text(mol, name="molecule"):
    smiles = Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
    )

    return f"{smiles}\t{name}\n"
