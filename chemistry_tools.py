from rdkit import Chem
from rdkit.Chem import (
    AllChem,
    Descriptors,
    Crippen,
    Lipinski,
    rdMolDescriptors,
    rdDepictor,
)
from rdkit.Chem.Draw import rdMolDraw2D

from PIL import Image

import io
import os
import subprocess
import tempfile


# ============================================================
# BASIC MOLECULE HANDLING
# ============================================================

def parse_smiles(smiles):
    try:
        return Chem.MolFromSmiles(smiles)
    except Exception:
        return None


def canonical_smiles(mol):
    return Chem.MolToSmiles(
        mol,
        canonical=True,
        isomericSmiles=True,
    )


# ============================================================
# PROPERTIES
# ============================================================

def get_properties(mol):
    return {
        "Formula": rdMolDescriptors.CalcMolFormula(mol),
        "Molecular Weight": round(Descriptors.MolWt(mol), 3),
        "Exact Mass": round(Descriptors.ExactMolWt(mol), 4),
        "LogP": round(Crippen.MolLogP(mol), 3),
        "TPSA": round(rdMolDescriptors.CalcTPSA(mol), 3),
        "H-bond Donors": Lipinski.NumHDonors(mol),
        "H-bond Acceptors": Lipinski.NumHAcceptors(mol),
        "Rotatable Bonds": Lipinski.NumRotatableBonds(mol),
        "Ring Count": Lipinski.RingCount(mol),
        "Heavy Atoms": mol.GetNumHeavyAtoms(),
        "Formal Charge": Chem.GetFormalCharge(mol),
    }


# ============================================================
# 2D SVG
# ============================================================

def mol_to_svg(
    mol,
    width=800,
    height=550,
    font_size=26,
    bond_width=3.0,
):
    mol2d = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol2d)

    drawer = rdMolDraw2D.MolDraw2DSVG(
        width,
        height,
    )

    opts = drawer.drawOptions()

    opts.fixedFontSize = font_size
    opts.minFontSize = max(12, font_size - 6)
    opts.maxFontSize = font_size + 6

    opts.bondLineWidth = bond_width
    opts.padding = 0.05

    drawer.DrawMolecule(mol2d)
    drawer.FinishDrawing()

    return drawer.GetDrawingText()


# ============================================================
# 2D PNG
# ============================================================

def mol_to_png(
    mol,
    width=1200,
    height=800,
    font_size=32,
    bond_width=3.0,
):
    """
    High-resolution PNG suitable for PowerPoint.
    """

    mol2d = Chem.Mol(mol)
    rdDepictor.Compute2DCoords(mol2d)

    drawer = rdMolDraw2D.MolDraw2DCairo(
        width,
        height,
    )

    opts = drawer.drawOptions()

    opts.fixedFontSize = font_size
    opts.minFontSize = max(14, font_size - 8)
    opts.maxFontSize = font_size + 8

    opts.bondLineWidth = bond_width
    opts.padding = 0.05

    drawer.DrawMolecule(mol2d)
    drawer.FinishDrawing()

    return drawer.GetDrawingText()


# ============================================================
# PNG -> JPG
# ============================================================

def png_to_jpg(png_data, quality=95):
    image = Image.open(
        io.BytesIO(png_data)
    )

    # JPEG does not support transparency
    if image.mode in ("RGBA", "LA"):
        background = Image.new(
            "RGB",
            image.size,
            "white"
        )

        background.paste(
            image,
            mask=image.getchannel("A")
        )

        image = background

    else:
        image = image.convert("RGB")

    buffer = io.BytesIO()

    image.save(
        buffer,
        format="JPEG",
        quality=quality,
    )

    return buffer.getvalue()


# ============================================================
# 3D GENERATION
# ============================================================

def generate_3d(mol):

    mol3d = Chem.AddHs(
        Chem.Mol(mol)
    )

    params = AllChem.ETKDGv3()
    params.randomSeed = 42

    result = AllChem.EmbedMolecule(
        mol3d,
        params,
    )

    if result != 0:
        params.useRandomCoords = True

        result = AllChem.EmbedMolecule(
            mol3d,
            params,
        )

    if result != 0:
        return None

    try:

        if AllChem.MMFFHasAllMoleculeParams(mol3d):

            AllChem.MMFFOptimizeMolecule(
                mol3d,
                maxIters=1000,
            )

        else:

            AllChem.UFFOptimizeMolecule(
                mol3d,
                maxIters=1000,
            )

    except Exception:

        try:
            AllChem.UFFOptimizeMolecule(
                mol3d,
                maxIters=1000,
            )
        except Exception:
            pass

    return mol3d


# ============================================================
# XYZ
# ============================================================

def mol_to_xyz(mol):

    conf = mol.GetConformer()

    lines = [
        str(mol.GetNumAtoms()),
        "Generated by Molecular Structure Toolkit",
    ]

    for atom in mol.GetAtoms():

        pos = conf.GetAtomPosition(
            atom.GetIdx()
        )

        lines.append(
            f"{atom.GetSymbol():<3s}"
            f"{pos.x:15.7f}"
            f"{pos.y:15.7f}"
            f"{pos.z:15.7f}"
        )

    return "\n".join(lines)


# ============================================================
# MOL / SDF / PDB
# ============================================================

def mol_to_mol(mol):

    return Chem.MolToMolBlock(
        mol,
        includeStereo=True,
    )


def mol_to_sdf(mol, name="Molecule"):

    mol_copy = Chem.Mol(mol)

    mol_copy.SetProp(
        "_Name",
        name
    )

    return (
        Chem.MolToMolBlock(
            mol_copy,
            includeStereo=True,
        )
        + "\n$$$$\n"
    )


def mol_to_pdb(mol):

    return Chem.MolToPDBBlock(
        mol
    )


# ============================================================
# MOL2 VIA OPEN BABEL
# ============================================================

def mol_to_mol2(
    mol,
    molecule_name="Molecule"
):

    with tempfile.TemporaryDirectory() as temp_dir:

        mol_file = os.path.join(
            temp_dir,
            "input.mol"
        )

        mol2_file = os.path.join(
            temp_dir,
            "output.mol2"
        )

        Chem.MolToMolFile(
            mol,
            mol_file
        )

        command = [
            "obabel",
            "-imol",
            mol_file,
            "-omol2",
            "-O",
            mol2_file,
        ]

        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                result.stderr
            )

        with open(
            mol2_file,
            "r",
            encoding="utf-8",
            errors="ignore",
        ) as f:
            text = f.read()

    text = text.replace(
        "input",
        molecule_name,
        1
    )

    return text