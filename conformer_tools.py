from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign

import math


# ============================================================
# SINGLE ETKDGv3 STRUCTURE — NO FORCE-FIELD OPTIMIZATION
# ============================================================

def generate_etkdg_structure(
    mol,
    random_seed=42,
    use_small_ring_torsions=True,
    use_macrocycle_torsions=True,
):
    """
    Generate ONE ETKDGv3 3D structure without force-field minimization.

    This is intended as the immediate 3D structure shown before the optional
    conformer search.

    Returns
    -------
    RDKit Mol or None
    """
    mol3d = Chem.AddHs(
        Chem.Mol(mol)
    )

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.enforceChirality = True
    params.useSmallRingTorsions = use_small_ring_torsions
    params.useMacrocycleTorsions = use_macrocycle_torsions

    result = AllChem.EmbedMolecule(
        mol3d,
        params,
    )

    if result != 0:
        # A fallback embedding attempt, still ETKDGv3 and still unoptimized.
        params.useRandomCoords = True
        result = AllChem.EmbedMolecule(
            mol3d,
            params,
        )

    if result != 0:
        return None

    return mol3d


# ============================================================
# FORCE-FIELD ENERGY
# ============================================================

def _get_mmff_energy(mol, conf_id):
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        return None

    properties = AllChem.MMFFGetMoleculeProperties(
        mol,
        mmffVariant="MMFF94s",
    )

    force_field = AllChem.MMFFGetMoleculeForceField(
        mol,
        properties,
        confId=conf_id,
    )

    if force_field is None:
        return None

    return float(force_field.CalcEnergy())


def _get_uff_energy(mol, conf_id):
    try:
        force_field = AllChem.UFFGetMoleculeForceField(
            mol,
            confId=conf_id,
        )

        if force_field is None:
            return None

        return float(force_field.CalcEnergy())

    except Exception:
        return None


# ============================================================
# OPTIMIZE ONE CONFORMER WITH EXPLICIT USER CHOICE
# ============================================================

def _optimize_conformer(
    mol,
    conf_id,
    force_field="MMFF94s",
    max_iters=1000,
):
    """
    Optimize one conformer with the selected force field.

    force_field must be:
        - 'MMFF94s'
        - 'UFF'

    No automatic switching is performed. This keeps MMFF94s and UFF
    experiments scientifically distinct.
    """
    force_field = force_field.upper()

    if force_field == "MMFF94S":
        if not AllChem.MMFFHasAllMoleculeParams(mol):
            return (
                "MMFF94s",
                -1,
                None,
                "MMFF94s parameters are unavailable for this molecule.",
            )

        try:
            properties = AllChem.MMFFGetMoleculeProperties(
                mol,
                mmffVariant="MMFF94s",
            )

            ff = AllChem.MMFFGetMoleculeForceField(
                mol,
                properties,
                confId=conf_id,
            )

            if ff is None:
                return (
                    "MMFF94s",
                    -1,
                    None,
                    "Could not construct the MMFF94s force field.",
                )

            status = ff.Minimize(
                maxIts=max_iters
            )

            energy = float(
                ff.CalcEnergy()
            )

            return (
                "MMFF94s",
                status,
                energy,
                None,
            )

        except Exception as exc:
            return (
                "MMFF94s",
                -1,
                None,
                str(exc),
            )

    if force_field == "UFF":
        try:
            ff = AllChem.UFFGetMoleculeForceField(
                mol,
                confId=conf_id,
            )

            if ff is None:
                return (
                    "UFF",
                    -1,
                    None,
                    "Could not construct the UFF force field.",
                )

            status = ff.Minimize(
                maxIts=max_iters
            )

            energy = float(
                ff.CalcEnergy()
            )

            return (
                "UFF",
                status,
                energy,
                None,
            )

        except Exception as exc:
            return (
                "UFF",
                -1,
                None,
                str(exc),
            )

    raise ValueError(
        "force_field must be either 'MMFF94s' or 'UFF'."
    )


# ============================================================
# GENERATE MULTIPLE CONFORMERS
# ============================================================

def generate_conformers(
    mol,
    num_conformers=50,
    random_seed=42,
    prune_rms_thresh=0.5,
    max_iters=1000,
    use_small_ring_torsions=True,
    use_macrocycle_torsions=True,
    force_field="MMFF94s",
):
    """
    Generate a conformational ensemble and optimize every conformer with the
    user-selected force field.

    Workflow
    --------
    molecule
      -> add explicit H
      -> ETKDGv3 multiple-conformer embedding
      -> embedding RMSD pruning
      -> MMFF94s OR UFF minimization (explicit user choice)
      -> energy ranking

    Important
    ---------
    Absolute energies from different force fields should not be compared
    directly. Compare relative energies within one selected force field.
    """
    mol3d = Chem.AddHs(
        Chem.Mol(mol)
    )

    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    params.pruneRmsThresh = prune_rms_thresh
    params.useSmallRingTorsions = use_small_ring_torsions
    params.useMacrocycleTorsions = use_macrocycle_torsions
    params.enforceChirality = True
    params.numThreads = 0

    conf_ids = list(
        AllChem.EmbedMultipleConfs(
            mol3d,
            numConfs=num_conformers,
            params=params,
        )
    )

    results = []

    for conf_id in conf_ids:

        method, status, energy, error = _optimize_conformer(
            mol3d,
            conf_id,
            force_field=force_field,
            max_iters=max_iters,
        )

        results.append(
            {
                "conf_id": int(conf_id),
                "method": method,
                "optimization_status": int(status),
                "energy": energy,
                "error": error,
            }
        )

    valid_energies = [
        r["energy"]
        for r in results
        if r["energy"] is not None
        and math.isfinite(r["energy"])
    ]

    if valid_energies:

        e_min = min(valid_energies)

        for r in results:
            if (
                r["energy"] is not None
                and math.isfinite(r["energy"])
            ):
                r["delta_energy"] = (
                    r["energy"] - e_min
                )
            else:
                r["delta_energy"] = None

    else:
        for r in results:
            r["delta_energy"] = None

    results.sort(
        key=lambda x: (
            float("inf")
            if x["energy"] is None
            else x["energy"]
        )
    )

    for rank, r in enumerate(
        results,
        start=1,
    ):
        r["rank"] = rank

    return mol3d, results


# ============================================================
# PAIRWISE RMSD
# ============================================================

def conformer_rmsd(
    mol,
    conf_id_1,
    conf_id_2,
):
    return float(
        rdMolAlign.GetBestRMS(
            mol,
            mol,
            prbId=conf_id_1,
            refId=conf_id_2,
        )
    )


# ============================================================
# RMSD DIVERSITY FILTER
# ============================================================

def filter_conformers_by_rmsd(
    mol,
    ranked_results,
    rmsd_threshold=1.0,
    max_keep=None,
):
    """
    Greedy RMSD diversity filtering in energy-ranked order.
    """
    if not ranked_results:
        return []

    # Only conformers with valid energies belong in an energy-ranked search.
    candidates = [
        item
        for item in ranked_results
        if item["energy"] is not None
    ]

    if not candidates:
        return []

    kept = []

    for candidate in candidates:

        candidate_id = candidate[
            "conf_id"
        ]

        if not kept:
            kept.append(candidate)

            if (
                max_keep is not None
                and len(kept) >= max_keep
            ):
                break

            continue

        is_diverse = True

        for selected in kept:

            rmsd = conformer_rmsd(
                mol,
                candidate_id,
                selected["conf_id"],
            )

            if rmsd < rmsd_threshold:
                is_diverse = False
                break

        if is_diverse:
            kept.append(candidate)

        if (
            max_keep is not None
            and len(kept) >= max_keep
        ):
            break

    return kept


# ============================================================
# EXTRACT ONE CONFORMER
# ============================================================

def extract_conformer(
    mol,
    conf_id,
):
    new_mol = Chem.Mol(
        mol
    )

    conf = mol.GetConformer(
        conf_id
    )

    new_conf = Chem.Conformer(
        conf
    )

    new_mol.RemoveAllConformers()

    new_mol.AddConformer(
        new_conf,
        assignId=True,
    )

    return new_mol


# ============================================================
# LOWEST ENERGY CONFORMER
# ============================================================

def get_lowest_energy_conformer(
    mol,
    results,
):
    valid = [
        r
        for r in results
        if r["energy"] is not None
    ]

    if not valid:
        return None, None

    best = min(
        valid,
        key=lambda x: x["energy"],
    )

    best_mol = extract_conformer(
        mol,
        best["conf_id"],
    )

    return best_mol, best
