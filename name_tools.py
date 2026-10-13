import requests
from urllib.parse import quote


PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"


def name_to_pubchem(name, timeout=20):
    """
    Resolve an IUPAC/common chemical name using PubChem PUG REST.

    Returns a dictionary containing:
    - query_name
    - cid
    - title
    - iupac_name
    - smiles
    """

    name = name.strip()

    if not name:
        raise ValueError("Chemical name is empty.")

    encoded_name = quote(name, safe="")

    url = (
        f"{PUBCHEM_BASE}/compound/name/{encoded_name}/property/"
        "SMILES,CanonicalSMILES,IsomericSMILES,IUPACName,Title/JSON"
    )

    response = requests.get(
        url,
        timeout=timeout,
        headers={
            "User-Toolkit": "Molecular-Structure-Toolkit/1.0"
        },
    )

    if response.status_code == 404:
        raise ValueError(
            f"No PubChem compound was found for '{name}'."
        )

    response.raise_for_status()

    data = response.json()

    properties = (
        data
        .get("PropertyTable", {})
        .get("Properties", [])
    )

    if not properties:
        raise ValueError(
            f"PubChem returned no structure for '{name}'."
        )

    item = properties[0]

    smiles = (
        item.get("SMILES")
        or item.get("IsomericSMILES")
        or item.get("CanonicalSMILES")
    )

    if not smiles:
        raise ValueError(
            f"PubChem returned a record for '{name}' but no SMILES."
        )

    return {
        "query_name": name,
        "cid": item.get("CID"),
        "title": item.get("Title"),
        "iupac_name": item.get("IUPACName"),
        "smiles": smiles,
    }


def resolve_name_list(names, timeout=20):
    """
    Resolve multiple names while keeping failures isolated to each name.
    """

    results = []

    for name in names:

        name = name.strip()

        if not name:
            continue

        try:
            result = name_to_pubchem(
                name,
                timeout=timeout,
            )

            results.append(
                {
                    "ok": True,
                    "result": result,
                }
            )

        except Exception as exc:

            results.append(
                {
                    "ok": False,
                    "query_name": name,
                    "error": str(exc),
                }
            )

    return results
