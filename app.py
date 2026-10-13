import streamlit as st
import streamlit.components.v1 as components
import py3Dmol

from chemistry_tools import (
    parse_smiles,
    canonical_smiles,
    get_properties,
    mol_to_svg,
    mol_to_png,
    png_to_jpg,
    mol_to_xyz,
    mol_to_mol,
    mol_to_sdf,
    mol_to_pdb,
    mol_to_mol2,
)

from conformer_tools import (
    generate_etkdg_structure,
    generate_conformers,
    filter_conformers_by_rmsd,
    get_lowest_energy_conformer,
)

from file_tools import (
    SUPPORTED_EXTENSIONS,
    load_molecules_from_bytes,
    mol_to_smiles_text,
    prepare_mol_for_analysis,
)

from name_tools import (
    resolve_name_list,
)



st.set_page_config(
    page_title="Molecular Structure Toolkit",
    page_icon="🧪",
    layout="wide",
)

st.title("🧪 Molecular Structure Toolkit")
st.caption(
    "SMILES analysis, publication-quality 2D depictions, "
    "3D conformer generation, MOL/MOL2/SDF conversion, and chemical-name lookup."
)

# Make the four main navigation panels very large and prominent.
st.markdown(
    """
    <style>
    div[data-testid="stTabs"] div[role="tablist"],
    div[data-baseweb="tab-list"] {
        gap: 26px !important;
        margin-top: 20px !important;
        margin-bottom: 30px !important;
        flex-wrap: wrap !important;
    }

    div[data-testid="stTabs"] button[role="tab"],
    button[data-baseweb="tab"] {
        min-height: 96px !important;
        padding: 22px 44px !important;
        border-radius: 14px 14px 0 0 !important;
    }

    div[data-testid="stTabs"] button[role="tab"] *,
    button[data-baseweb="tab"] * {
        font-size: 32px !important;
        font-weight: 800 !important;
        line-height: 1.15 !important;
        letter-spacing: 0.2px !important;
    }

    div[data-testid="stTabs"] button[role="tab"][aria-selected="true"],
    button[data-baseweb="tab"][aria-selected="true"] {
        border-bottom-width: 5px !important;
    }

    @media (max-width: 1500px) {
        div[data-testid="stTabs"] button[role="tab"],
        button[data-baseweb="tab"] {
            min-height: 86px !important;
            padding: 18px 32px !important;
        }

        div[data-testid="stTabs"] button[role="tab"] *,
        button[data-baseweb="tab"] * {
            font-size: 28px !important;
        }
    }
    </style>
    """,
    unsafe_allow_html=True,
)

if "smiles_list" not in st.session_state:
    st.session_state.smiles_list = []

if "conformer_cache" not in st.session_state:
    st.session_state.conformer_cache = {}

if "file_records" not in st.session_state:
    st.session_state.file_records = []

if "name_results" not in st.session_state:
    st.session_state.name_results = []

if "etkdg_cache" not in st.session_state:
    st.session_state.etkdg_cache = {}



default_smiles = """CC1=CN=C(C(=C1OC)C)C[S@@](=O)C2=NC3=C(N2)C=C(C=C3)OC
CC1=CN=C(C(=C1OC)C)C[S@](=O)C2=NC3=C(N2)C=C(C=C3)OC
"""

with st.expander("SMILES workspace", expanded=True):
    smiles_text = st.text_area(
        "SMILES list",
        value=default_smiles,
        height=170,
        placeholder="Enter one SMILES per line...",
        key="smiles_input",
    )

    c_load, c_clear = st.columns(2)

    with c_load:
        if st.button("Load SMILES", type="primary", key="load_smiles_button"):
            st.session_state.smiles_list = [
                x.strip() for x in smiles_text.splitlines() if x.strip()
            ]
            st.session_state.conformer_cache = {}
            st.session_state.etkdg_cache = {}

    with c_clear:
        if st.button("Clear loaded molecules", key="clear_smiles_button"):
            st.session_state.smiles_list = []
            st.session_state.conformer_cache = {}

if st.session_state.smiles_list:
    st.success(
        f"{len(st.session_state.smiles_list)} molecule(s) loaded into the workspace."
    )


tab_2d, tab_3d, tab_file, tab_name = st.tabs(
    [
        "2D Structure",
        "3D Conformers",
        "File → Structure",
        "Name → Structure",
    ]
)


# ============================================================
# TAB 1 — 2D STRUCTURE
# ============================================================

with tab_2d:
    st.header("2D Structure and Calculated Descriptors")

    if not st.session_state.smiles_list:
        st.info("Load one or more SMILES in the workspace above.")
    else:
        with st.expander("2D drawing settings", expanded=False):
            set1, set2, set3, set4, set5 = st.columns(5)

            with set1:
                font_size = st.slider(
                    "Atom label font size", 14, 45, 28, 1, key="font_size_2d"
                )
            with set2:
                bond_width = st.slider(
                    "Bond thickness", 1.0, 5.0, 3.0, 0.25, key="bond_width_2d"
                )
            with set3:
                image_width = st.slider(
                    "Image width", 400, 1200, 800, 50, key="image_width_2d"
                )
            with set4:
                image_height = st.slider(
                    "Image height", 300, 900, 550, 50, key="image_height_2d"
                )
            with set5:
                png_scale = st.selectbox(
                    "High-res PNG",
                    options=[1, 2, 3, 4],
                    index=2,
                    format_func=lambda x: f"{x}×",
                    key="png_scale_2d",
                )

        for index, smiles in enumerate(st.session_state.smiles_list, start=1):
            st.divider()
            st.subheader(f"Molecule {index}")

            mol = parse_smiles(smiles)
            if mol is None:
                st.error(f"Invalid SMILES: {smiles}")
                continue

            canon = canonical_smiles(mol)
            properties = get_properties(mol)

            st.caption("Input SMILES")
            st.code(smiles)
            st.caption("Canonical isomeric SMILES")
            st.code(canon)

            left, right = st.columns([1, 2], gap="large")

            with left:
                st.markdown("### Calculated Molecular Descriptors")
                for key, value in properties.items():
                    st.write(f"**{key}:** {value}")

            with right:
                st.markdown("### 2D Structure")

                svg = mol_to_svg(
                    mol,
                    width=image_width,
                    height=image_height,
                    font_size=font_size,
                    bond_width=bond_width,
                )

                components.html(
                    svg,
                    width=image_width,
                    height=image_height + 20,
                    scrolling=False,
                )

                png_data = mol_to_png(
                    mol,
                    width=image_width,
                    height=image_height,
                    font_size=font_size,
                    bond_width=bond_width,
                )

                jpg_data = png_to_jpg(png_data)

                png_highres = mol_to_png(
                    mol,
                    width=image_width * png_scale,
                    height=image_height * png_scale,
                    font_size=font_size * png_scale,
                    bond_width=bond_width * png_scale,
                )

                d1, d2, d3, d4 = st.columns(4)

                with d1:
                    st.download_button(
                        "Download SVG",
                        svg,
                        file_name=f"molecule_{index}_2D.svg",
                        mime="image/svg+xml",
                        key=f"2d_svg_{index}",
                        on_click="ignore",
                    )
                with d2:
                    st.download_button(
                        "Download PNG",
                        png_data,
                        file_name=f"molecule_{index}_2D.png",
                        mime="image/png",
                        key=f"2d_png_{index}",
                        on_click="ignore",
                    )
                with d3:
                    st.download_button(
                        "Download JPG",
                        jpg_data,
                        file_name=f"molecule_{index}_2D.jpg",
                        mime="image/jpeg",
                        key=f"2d_jpg_{index}",
                        on_click="ignore",
                    )
                with d4:
                    st.download_button(
                        "Download PNG (High Res)",
                        png_highres,
                        file_name=f"molecule_{index}_2D_highres.png",
                        mime="image/png",
                        key=f"2d_png_high_{index}",
                        on_click="ignore",
                    )


# ============================================================
# TAB 2 — 3D CONFORMERS
# ============================================================

with tab_3d:

    st.header(
        "3D Structure and Conformer Search"
    )

    st.write(
        "A single **ETKDGv3 structure is generated first without force-field "
        "optimization**. Conformer searching is a separate optional step. "
        "For the conformer search, choose either **MMFF94s** or **UFF**."
    )

    if not st.session_state.smiles_list:

        st.info(
            "Load one or more SMILES in the workspace above."
        )

    else:

        with st.expander(
            "Conformer-search settings",
            expanded=False,
        ):

            s0, s1, s2, s3, s4 = st.columns(5)

            with s0:

                selected_force_field = st.radio(
                    "Force field",
                    ["MMFF94s", "UFF"],
                    horizontal=True,
                    key="conformer_force_field",
                    help=(
                        "The selected force field is used for the complete "
                        "conformer-search run. Energies from MMFF94s and UFF "
                        "should not be compared directly."
                    ),
                )

            with s1:

                num_conformers = st.slider(
                    "Number of conformers",
                    5,
                    200,
                    50,
                    5,
                    key="num_conformers_3d",
                )

            with s2:

                embed_rmsd = st.slider(
                    "Embedding RMSD pruning (Å)",
                    0.0,
                    2.0,
                    0.5,
                    0.1,
                    key="embed_rmsd_3d",
                )

            with s3:

                final_rmsd = st.slider(
                    "Final RMSD diversity threshold (Å)",
                    0.1,
                    3.0,
                    1.0,
                    0.1,
                    key="final_rmsd_3d",
                )

            with s4:

                max_selected = st.slider(
                    "Maximum conformers retained",
                    1,
                    50,
                    10,
                    1,
                    key="max_selected_3d",
                )

        for index, smiles in enumerate(
            st.session_state.smiles_list,
            start=1,
        ):

            st.divider()
            st.subheader(
                f"Molecule {index}"
            )

            mol = parse_smiles(
                smiles
            )

            if mol is None:

                st.error(
                    f"Invalid SMILES: {smiles}"
                )

                continue

            st.caption(
                "Input SMILES"
            )
            st.code(
                smiles
            )

            # ------------------------------------------------
            # STAGE 1 — SINGLE ETKDGv3 STRUCTURE
            # ------------------------------------------------

            st.markdown(
                "### 1. Initial ETKDGv3 Structure"
            )

            st.caption(
                "Generated directly from the molecular graph. "
                "No MMFF94s or UFF minimization has been applied."
            )

            etkdg_key = (
                smiles,
                42,
            )

            cached_etkdg = st.session_state.etkdg_cache.get(
                index
            )

            if (
                cached_etkdg is None
                or cached_etkdg["signature"] != etkdg_key
            ):

                with st.spinner(
                    "Generating initial ETKDGv3 structure..."
                ):

                    etkdg_mol = generate_etkdg_structure(
                        mol,
                        random_seed=42,
                    )

                st.session_state.etkdg_cache[index] = {
                    "signature": etkdg_key,
                    "mol": etkdg_mol,
                }

            etkdg_mol = (
                st.session_state.etkdg_cache[index]["mol"]
            )

            if etkdg_mol is None:

                st.error(
                    "ETKDGv3 could not generate a 3D structure for this molecule."
                )

            else:

                etkdg_block = mol_to_mol(
                    etkdg_mol
                )

                etkdg_viewer = py3Dmol.view(
                    width=900,
                    height=520,
                )

                etkdg_viewer.addModel(
                    etkdg_block,
                    "mol",
                )

                etkdg_viewer.setStyle(
                    {
                        "stick": {
                            "radius": 0.15
                        }
                    }
                )

                etkdg_viewer.addStyle(
                    {},
                    {
                        "sphere": {
                            "scale": 0.28
                        }
                    },
                )

                etkdg_viewer.setBackgroundColor(
                    "white"
                )

                etkdg_viewer.zoomTo()

                etkdg_html = etkdg_viewer._make_html()

                components.html(
                    etkdg_html,
                    height=540,
                )

                st.markdown(
                    "#### Download raw ETKDGv3 structure"
                )

                e1, e2, e3, e4, e5 = st.columns(5)

                with e1:
                    st.download_button(
                        "ETKDGv3 XYZ",
                        mol_to_xyz(etkdg_mol),
                        file_name=f"molecule_{index}_ETKDGv3.xyz",
                        mime="text/plain",
                        key=f"etkdg_xyz_{index}",
                        on_click="ignore",
                    )

                with e2:
                    st.download_button(
                        "ETKDGv3 MOL",
                        mol_to_mol(etkdg_mol),
                        file_name=f"molecule_{index}_ETKDGv3.mol",
                        mime="chemical/x-mdl-molfile",
                        key=f"etkdg_mol_{index}",
                        on_click="ignore",
                    )

                with e3:
                    st.download_button(
                        "ETKDGv3 SDF",
                        mol_to_sdf(
                            etkdg_mol,
                            f"Molecule_{index}_ETKDGv3",
                        ),
                        file_name=f"molecule_{index}_ETKDGv3.sdf",
                        mime="chemical/x-mdl-sdfile",
                        key=f"etkdg_sdf_{index}",
                        on_click="ignore",
                    )

                with e4:
                    st.download_button(
                        "ETKDGv3 PDB",
                        mol_to_pdb(etkdg_mol),
                        file_name=f"molecule_{index}_ETKDGv3.pdb",
                        mime="chemical/x-pdb",
                        key=f"etkdg_pdb_{index}",
                        on_click="ignore",
                    )

                with e5:
                    try:
                        etkdg_mol2 = mol_to_mol2(
                            etkdg_mol,
                            f"Molecule_{index}_ETKDGv3",
                        )

                        st.download_button(
                            "ETKDGv3 MOL2",
                            etkdg_mol2,
                            file_name=f"molecule_{index}_ETKDGv3.mol2",
                            mime="text/plain",
                            key=f"etkdg_mol2_{index}",
                            on_click="ignore",
                        )

                    except Exception as exc:
                        st.caption(
                            f"MOL2 export unavailable: {exc}"
                        )

            # ------------------------------------------------
            # STAGE 2 — OPTIONAL CONFORMER SEARCH
            # ------------------------------------------------

            st.markdown(
                "### 2. Optional Conformer Search"
            )

            st.write(
                f"Current force field: **{selected_force_field}**. "
                "The search does not start until you press the button below."
            )

            current_signature = (
                smiles,
                selected_force_field,
                num_conformers,
                embed_rmsd,
                final_rmsd,
                max_selected,
            )

            if st.button(
                f"Run {selected_force_field} conformer search — Molecule {index}",
                key=f"generate_conf_{index}",
                type="primary",
            ):

                with st.spinner(
                    f"Generating and optimizing up to {num_conformers} "
                    f"conformers with {selected_force_field}..."
                ):

                    ensemble_mol, conformer_results = generate_conformers(
                        mol,
                        num_conformers=num_conformers,
                        prune_rms_thresh=embed_rmsd,
                        force_field=selected_force_field,
                    )

                    selected_results = filter_conformers_by_rmsd(
                        ensemble_mol,
                        conformer_results,
                        rmsd_threshold=final_rmsd,
                        max_keep=max_selected,
                    )

                    lowest_mol, lowest_info = get_lowest_energy_conformer(
                        ensemble_mol,
                        conformer_results,
                    )

                    st.session_state.conformer_cache[index] = {
                        "signature": current_signature,
                        "conformer_results": conformer_results,
                        "selected_results": selected_results,
                        "lowest_mol": lowest_mol,
                        "lowest_info": lowest_info,
                    }

            cached = st.session_state.conformer_cache.get(
                index
            )

            if cached is None:

                st.info(
                    "The ETKDGv3 structure above is ready. "
                    "Run the conformer search only if you need an "
                    "energy-ranked ensemble."
                )

                continue

            if cached["signature"] != current_signature:

                st.warning(
                    "The conformer-search settings have changed since the "
                    "cached run. Press the search button to generate a new ensemble."
                )

            conformer_results = cached[
                "conformer_results"
            ]

            selected_results = cached[
                "selected_results"
            ]

            lowest_mol = cached[
                "lowest_mol"
            ]

            lowest_info = cached[
                "lowest_info"
            ]

            num_generated = len(
                conformer_results
            )

            num_optimized = sum(
                1
                for result in conformer_results
                if result["energy"] is not None
            )

            num_selected = len(
                selected_results
            )

            m1, m2, m3, m4 = st.columns(4)

            with m1:
                st.metric(
                    "Generated",
                    num_generated,
                )

            with m2:
                st.metric(
                    f"{selected_force_field} optimized",
                    num_optimized,
                )

            with m3:
                st.metric(
                    "RMSD-diverse retained",
                    num_selected,
                )

            with m4:
                st.metric(
                    "Force field",
                    selected_force_field,
                )

            if (
                lowest_mol is None
                or lowest_info is None
            ):

                st.error(
                    f"No valid {selected_force_field} energy could be obtained. "
                    "This can happen when the selected force field does not "
                    "support the molecule."
                )

                errors = sorted(
                    {
                        item.get("error")
                        for item in conformer_results
                        if item.get("error")
                    }
                )

                if errors:
                    with st.expander(
                        "Force-field details",
                        expanded=False,
                    ):
                        for error in errors:
                            st.write(error)

                continue

            st.markdown(
                "### Lowest-Energy Conformer"
            )

            i1, i2, i3 = st.columns(3)

            with i1:
                st.metric(
                    "Conformer ID",
                    lowest_info["conf_id"],
                )

            with i2:
                st.metric(
                    "Force Field",
                    lowest_info["method"],
                )

            with i3:
                st.metric(
                    "Energy",
                    f"{lowest_info['energy']:.3f} kcal/mol",
                )

            mol3d = lowest_mol

            mol_block = mol_to_mol(
                mol3d
            )

            viewer = py3Dmol.view(
                width=900,
                height=550,
            )

            viewer.addModel(
                mol_block,
                "mol",
            )

            viewer.setStyle(
                {
                    "stick": {
                        "radius": 0.15
                    }
                }
            )

            viewer.addStyle(
                {},
                {
                    "sphere": {
                        "scale": 0.28
                    }
                },
            )

            viewer.setBackgroundColor(
                "white"
            )

            viewer.zoomTo()

            viewer_html = viewer._make_html()

            components.html(
                viewer_html,
                height=570,
            )

            with st.expander(
                "Advanced conformer energies and RMSD-selected ensemble",
                expanded=False,
            ):

                table_data = []

                for result in conformer_results:

                    table_data.append(
                        {
                            "Rank": result["rank"],
                            "Conformer ID": result["conf_id"],
                            "Force Field": result["method"],
                            "Energy (kcal/mol)": (
                                round(result["energy"], 4)
                                if result["energy"] is not None
                                else None
                            ),
                            "ΔE (kcal/mol)": (
                                round(result["delta_energy"], 4)
                                if result["delta_energy"] is not None
                                else None
                            ),
                            "Optimization Status":
                                result["optimization_status"],
                        }
                    )

                st.markdown(
                    "#### Full energy ranking"
                )

                st.dataframe(
                    table_data,
                    use_container_width=True,
                    hide_index=True,
                )

                selected_table = []

                for result in selected_results:

                    selected_table.append(
                        {
                            "Rank": result["rank"],
                            "Conformer ID": result["conf_id"],
                            "Energy (kcal/mol)": round(
                                result["energy"],
                                4,
                            ),
                            "ΔE (kcal/mol)": round(
                                result["delta_energy"],
                                4,
                            ),
                        }
                    )

                st.markdown(
                    "#### RMSD-diverse conformers"
                )

                st.dataframe(
                    selected_table,
                    use_container_width=True,
                    hide_index=True,
                )

            st.markdown(
                "### Download lowest-energy conformer"
            )

            xyz_data = mol_to_xyz(
                mol3d
            )

            mol_data = mol_to_mol(
                mol3d
            )

            sdf_data = mol_to_sdf(
                mol3d,
                f"Molecule_{index}_{selected_force_field}_lowest",
            )

            pdb_data = mol_to_pdb(
                mol3d
            )

            try:

                mol2_data = mol_to_mol2(
                    mol3d,
                    f"Molecule_{index}_{selected_force_field}_lowest",
                )

                mol2_ok = True

            except Exception as error:

                mol2_data = None
                mol2_ok = False

                st.warning(
                    "MOL2 generation failed. "
                    "Check that Open Babel is installed correctly."
                )

                st.code(
                    str(error)
                )

            d1, d2, d3 = st.columns(3)

            with d1:

                st.download_button(
                    "Download XYZ",
                    xyz_data,
                    file_name=f"molecule_{index}_{selected_force_field}_lowest.xyz",
                    mime="text/plain",
                    key=f"3d_xyz_{index}",
                    on_click="ignore",
                )

                st.download_button(
                    "Download MOL",
                    mol_data,
                    file_name=f"molecule_{index}_{selected_force_field}_lowest.mol",
                    mime="chemical/x-mdl-molfile",
                    key=f"3d_mol_{index}",
                    on_click="ignore",
                )

            with d2:

                st.download_button(
                    "Download SDF",
                    sdf_data,
                    file_name=f"molecule_{index}_{selected_force_field}_lowest.sdf",
                    mime="chemical/x-mdl-sdfile",
                    key=f"3d_sdf_{index}",
                    on_click="ignore",
                )

                st.download_button(
                    "Download PDB",
                    pdb_data,
                    file_name=f"molecule_{index}_{selected_force_field}_lowest.pdb",
                    mime="chemical/x-pdb",
                    key=f"3d_pdb_{index}",
                    on_click="ignore",
                )

            with d3:

                if mol2_ok:

                    st.download_button(
                        "Download MOL2",
                        mol2_data,
                        file_name=f"molecule_{index}_{selected_force_field}_lowest.mol2",
                        mime="text/plain",
                        key=f"3d_mol2_{index}",
                        on_click="ignore",
                    )

                st.download_button(
                    "Download Interactive 3D HTML",
                    viewer_html,
                    file_name=f"molecule_{index}_{selected_force_field}_3D.html",
                    mime="text/html",
                    key=f"3d_html_{index}",
                    on_click="ignore",
                )


# ============================================================
# TAB 3 — FILE -> STRUCTURE
# ============================================================

with tab_file:
    st.header("MOL / MOL2 / SDF → SMILES / 2D Analysis")

    st.write(
        "Upload MOL, MOL2, or SDF files. These formats normally preserve "
        "bond-order information much better than PDB and are therefore "
        "preferred for SMILES reconstruction and descriptor analysis."
    )

    st.info(
        "PDB upload is intentionally disabled in this version because "
        "small-molecule PDB files may not contain reliable bond-order information."
    )

    uploaded_files = st.file_uploader(
        "Upload MOL, MOL2, or SDF files",
        type=SUPPORTED_EXTENSIONS,
        accept_multiple_files=True,
        key="structure_file_uploader",
    )

    hide_explicit_h = st.checkbox(
        "Hide explicit hydrogens in 2D structure and SMILES analysis",
        value=True,
        key="file_hide_explicit_h",
        help=(
            "Recommended. Removes explicit hydrogens on CH3/CH2/carbon rings "
            "and produces a conventional publication-style structure. "
            "Chemically required hydrogens remain represented where needed."
        ),
    )

    if st.button(
        "Analyze uploaded files",
        type="primary",
        key="analyze_uploaded_files",
    ):
        records = []

        for uploaded_file in uploaded_files or []:
            try:
                records.extend(
                    load_molecules_from_bytes(
                        uploaded_file.name,
                        uploaded_file.getvalue(),
                    )
                )
            except Exception as exc:
                records.append(
                    {
                        "source_name": uploaded_file.name,
                        "record_index": 1,
                        "format": (
                            uploaded_file.name.rsplit(".", 1)[-1].upper()
                            if "." in uploaded_file.name
                            else "UNKNOWN"
                        ),
                        "mol": None,
                        "warning": str(exc),
                    }
                )

        st.session_state.file_records = records

    if st.session_state.file_records:
        with st.expander("File-analysis 2D drawing settings", expanded=False):
            fs1, fs2, fs3, fs4 = st.columns(4)

            with fs1:
                file_font_size = st.slider(
                    "Atom font size", 14, 45, 28, key="file_font_size"
                )
            with fs2:
                file_bond_width = st.slider(
                    "Bond thickness", 1.0, 5.0, 3.0, 0.25, key="file_bond_width"
                )
            with fs3:
                file_image_width = st.slider(
                    "Image width", 400, 1200, 800, 50, key="file_image_width"
                )
            with fs4:
                file_image_height = st.slider(
                    "Image height", 300, 900, 550, 50, key="file_image_height"
                )

        for file_index, record in enumerate(
            st.session_state.file_records,
            start=1,
        ):
            st.divider()

            title = record["source_name"]
            if record["record_index"] > 1:
                title += f" — record {record['record_index']}"

            st.subheader(title)

            if record["warning"]:
                st.warning(record["warning"])

            raw_mol = record["mol"]

            if raw_mol is None:
                st.error("This structure could not be converted to an RDKit molecule.")
                continue

            mol = prepare_mol_for_analysis(
                raw_mol,
                remove_explicit_h=hide_explicit_h,
            )

            if mol is None:
                st.error("Failed to prepare molecule for analysis.")
                continue

            try:
                smiles = canonical_smiles(mol)
            except Exception as exc:
                st.error(f"SMILES generation failed: {exc}")
                continue

            try:
                properties = get_properties(mol)
            except Exception as exc:
                st.error(f"Descriptor calculation failed: {exc}")
                continue

            st.caption(
                f"Source format: {record['format']} | "
                "Generated canonical isomeric SMILES"
            )
            st.code(smiles)

            left, right = st.columns([1, 2], gap="large")

            with left:
                st.markdown("### Calculated Molecular Descriptors")

                for key, value in properties.items():
                    st.write(f"**{key}:** {value}")

                smiles_file = mol_to_smiles_text(
                    mol,
                    name=record["source_name"],
                )

                st.download_button(
                    "Download SMILES (.smi)",
                    smiles_file,
                    file_name=f"{record['source_name'].rsplit('.', 1)[0]}.smi",
                    mime="text/plain",
                    key=f"file_smi_{file_index}",
                    on_click="ignore",
                )

                if st.button(
                    "Add SMILES to workspace",
                    key=f"add_file_smiles_{file_index}",
                ):
                    if smiles not in st.session_state.smiles_list:
                        st.session_state.smiles_list.append(smiles)
                        st.session_state.conformer_cache = {}
                        st.success("SMILES added to the shared workspace.")
                    else:
                        st.info("This SMILES is already in the workspace.")

            with right:
                st.markdown("### Reconstructed 2D Structure")

                svg = mol_to_svg(
                    mol,
                    width=file_image_width,
                    height=file_image_height,
                    font_size=file_font_size,
                    bond_width=file_bond_width,
                )

                components.html(
                    svg,
                    width=file_image_width,
                    height=file_image_height + 20,
                    scrolling=False,
                )

                png_data = mol_to_png(
                    mol,
                    width=file_image_width,
                    height=file_image_height,
                    font_size=file_font_size,
                    bond_width=file_bond_width,
                )

                p1, p2 = st.columns(2)

                with p1:
                    st.download_button(
                        "Download reconstructed SVG",
                        svg,
                        file_name=f"file_{file_index}_2D.svg",
                        mime="image/svg+xml",
                        key=f"file_svg_{file_index}",
                        on_click="ignore",
                    )

                with p2:
                    st.download_button(
                        "Download reconstructed PNG",
                        png_data,
                        file_name=f"file_{file_index}_2D.png",
                        mime="image/png",
                        key=f"file_png_{file_index}",
                        on_click="ignore",
                    )


# ============================================================
# TAB 4 — NAME -> STRUCTURE
# ============================================================

with tab_name:

    st.header(
        "Chemical Name → SMILES / 2D Structure"
    )

    st.write(
        "Enter one or more IUPAC names, systematic names, common chemical "
        "names, or drug names recognized by PubChem. Use one name per line."
    )

    st.info(
        "This panel uses the PubChem web service and therefore requires an "
        "internet connection. Molecular descriptors are then calculated "
        "locally with RDKit from the retrieved SMILES."
    )

    names_text = st.text_area(
        "Chemical names",
        height=150,
        placeholder=(
            "acetylsalicylic acid\n"
            "ibuprofen\n"
            "caffeine"
        ),
        key="chemical_names_input",
    )

    if st.button(
        "Resolve Names",
        type="primary",
        key="resolve_names_button",
    ):

        names = [
            line.strip()
            for line in names_text.splitlines()
            if line.strip()
        ]

        if not names:

            st.warning(
                "Enter at least one chemical name."
            )

        else:

            with st.spinner(
                "Querying PubChem..."
            ):

                st.session_state.name_results = resolve_name_list(
                    names
                )

    if st.session_state.name_results:

        with st.expander(
            "Name-analysis 2D drawing settings",
            expanded=False,
        ):

            ns1, ns2, ns3, ns4, ns5 = st.columns(5)

            with ns1:

                name_font_size = st.slider(
                    "Atom font size",
                    14,
                    45,
                    28,
                    key="name_font_size",
                )

            with ns2:

                name_bond_width = st.slider(
                    "Bond thickness",
                    1.0,
                    5.0,
                    3.0,
                    0.25,
                    key="name_bond_width",
                )

            with ns3:

                name_image_width = st.slider(
                    "Image width",
                    400,
                    1200,
                    800,
                    50,
                    key="name_image_width",
                )

            with ns4:

                name_image_height = st.slider(
                    "Image height",
                    300,
                    900,
                    550,
                    50,
                    key="name_image_height",
                )

            with ns5:

                name_png_scale = st.selectbox(
                    "High-res PNG",
                    options=[1, 2, 3, 4],
                    index=2,
                    format_func=lambda x: f"{x}×",
                    key="name_png_scale",
                )

        for name_index, entry in enumerate(
            st.session_state.name_results,
            start=1,
        ):

            st.divider()

            if not entry["ok"]:

                st.subheader(
                    entry["query_name"]
                )

                st.error(
                    entry["error"]
                )

                continue

            result = entry[
                "result"
            ]

            query_name = result[
                "query_name"
            ]

            st.subheader(
                query_name
            )

            smiles = result[
                "smiles"
            ]

            mol = parse_smiles(
                smiles
            )

            if mol is None:

                st.error(
                    "PubChem returned a SMILES string, but RDKit could not parse it."
                )

                continue

            canon = canonical_smiles(
                mol
            )

            properties = get_properties(
                mol
            )

            meta1, meta2 = st.columns(2)

            with meta1:

                st.write(
                    f"**PubChem CID:** {result['cid']}"
                )

                if result["title"]:

                    st.write(
                        f"**PubChem title:** {result['title']}"
                    )

            with meta2:

                if result["iupac_name"]:

                    st.write(
                        f"**PubChem IUPAC name:** {result['iupac_name']}"
                    )

            st.caption(
                "SMILES returned by PubChem"
            )

            st.code(
                smiles
            )

            st.caption(
                "RDKit canonical isomeric SMILES"
            )

            st.code(
                canon
            )

            left, right = st.columns(
                [1, 2],
                gap="large",
            )

            with left:

                st.markdown(
                    "### Calculated Molecular Descriptors"
                )

                for key, value in properties.items():

                    st.write(
                        f"**{key}:** {value}"
                    )

                smiles_file = mol_to_smiles_text(
                    mol,
                    name=query_name,
                )

                st.download_button(
                    "Download SMILES (.smi)",
                    smiles_file,
                    file_name=f"name_{name_index}.smi",
                    mime="text/plain",
                    key=f"name_smi_{name_index}",
                    on_click="ignore",
                )

                if st.button(
                    "Add SMILES to workspace",
                    key=f"name_add_workspace_{name_index}",
                ):

                    if canon not in st.session_state.smiles_list:

                        st.session_state.smiles_list.append(
                            canon
                        )

                        st.session_state.conformer_cache = {}

                        st.success(
                            "SMILES added to the shared workspace."
                        )

                    else:

                        st.info(
                            "This molecule is already in the workspace."
                        )

            with right:

                st.markdown(
                    "### 2D Structure"
                )

                svg = mol_to_svg(
                    mol,
                    width=name_image_width,
                    height=name_image_height,
                    font_size=name_font_size,
                    bond_width=name_bond_width,
                )

                components.html(
                    svg,
                    width=name_image_width,
                    height=name_image_height + 20,
                    scrolling=False,
                )

                png_data = mol_to_png(
                    mol,
                    width=name_image_width,
                    height=name_image_height,
                    font_size=name_font_size,
                    bond_width=name_bond_width,
                )

                jpg_data = png_to_jpg(
                    png_data
                )

                png_highres = mol_to_png(
                    mol,
                    width=name_image_width * name_png_scale,
                    height=name_image_height * name_png_scale,
                    font_size=name_font_size * name_png_scale,
                    bond_width=name_bond_width * name_png_scale,
                )

                n1, n2, n3, n4 = st.columns(4)

                with n1:

                    st.download_button(
                        "Download SVG",
                        svg,
                        file_name=f"name_{name_index}_2D.svg",
                        mime="image/svg+xml",
                        key=f"name_svg_{name_index}",
                        on_click="ignore",
                    )

                with n2:

                    st.download_button(
                        "Download PNG",
                        png_data,
                        file_name=f"name_{name_index}_2D.png",
                        mime="image/png",
                        key=f"name_png_{name_index}",
                        on_click="ignore",
                    )

                with n3:

                    st.download_button(
                        "Download JPG",
                        jpg_data,
                        file_name=f"name_{name_index}_2D.jpg",
                        mime="image/jpeg",
                        key=f"name_jpg_{name_index}",
                        on_click="ignore",
                    )

                with n4:

                    st.download_button(
                        "Download PNG (High Res)",
                        png_highres,
                        file_name=f"name_{name_index}_2D_highres.png",
                        mime="image/png",
                        key=f"name_png_high_{name_index}",
                        on_click="ignore",
                    )

# ============================================================
# FOOTER
# ============================================================

st.divider()

st.caption(
    "Descriptors are calculated with RDKit. "
    "Initial 3D structures use ETKDGv3 without force-field minimization. "
    "Optional conformer searches use the user-selected MMFF94s or UFF force field. "
    "MOL2 export uses Open Babel. "
    "Chemical-name resolution uses PubChem and requires internet access."
)
