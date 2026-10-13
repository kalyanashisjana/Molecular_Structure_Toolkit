# Third-Party Licenses and Notices

This project includes or depends on third-party open-source software and public
web services. The MIT License in the root `LICENSE` file applies only to the
original source code of this project. Third-party components remain governed by
their own licenses and terms.

## RDKit

Project: https://www.rdkit.org/  
Source: https://github.com/rdkit/rdkit  
License: BSD 3-Clause License

RDKit is used for:

- SMILES parsing and canonicalization
- molecular descriptor calculation
- 2D molecular depiction
- ETKDGv3 3D coordinate generation
- MMFF94s and UFF optimization
- conformer comparison and RMSD calculations
- molecular file handling

RDKit is distributed under the BSD 3-Clause License.

Copyright notice from the RDKit license:

> Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum,
> Julie Penzotti and others. All rights reserved.

When redistributing RDKit itself, retain its copyright notice, license
conditions, and disclaimer as required by the BSD 3-Clause License.

## Open Babel

Project: https://openbabel.org/  
Source: https://github.com/openbabel/openbabel  
License: GNU General Public License version 2 (GPL-2.0-only)

Open Babel is used through its command-line executable (`obabel`) for selected
molecular-format conversion tasks, including MOL2 export and fallback
conversion.

Open Babel is not relicensed under this project's MIT License.

The Open Babel project states that Open Babel is distributed under GPL v2.
Its documentation also distinguishes direct linking to the Open Babel library
from use or distribution of the unmodified command-line program.

If you redistribute Open Babel binaries together with this project, you are
responsible for satisfying any GPL obligations that apply to that distribution.

Suggested scientific citation:

Noel M. O'Boyle, Michael Banck, Craig A. James, Chris Morley,
Tim Vandermeersch, and Geoffrey R. Hutchison.
"Open Babel: An open chemical toolbox."
Journal of Cheminformatics 3, 33 (2011).

## Streamlit

Project: https://streamlit.io/  
Source: https://github.com/streamlit/streamlit  
License: Apache License 2.0

Streamlit provides the graphical web interface.

The Streamlit open-source library is licensed under Apache-2.0. The Streamlit
Community Cloud service is additionally governed by Streamlit / Snowflake
service terms.

If redistributing Streamlit itself, retain the notices and license information
required by the Apache License 2.0.

## 3Dmol.js / py3Dmol

3Dmol.js project: https://3dmol.org/  
Source: https://github.com/3dmol/3Dmol.js  
License: BSD 3-Clause License

This project uses `py3Dmol` for interactive molecular visualization.

The 3Dmol.js project states that it is licensed under the BSD 3-Clause License
and asks that use of the software be properly acknowledged and cited.

Copyright notice from the 3Dmol.js license:

> Copyright (c) 2014, University of Pittsburgh and contributors.
> All rights reserved.

When redistributing 3Dmol.js itself, retain its copyright notice, license
conditions, and disclaimer.

## PubChem PUG REST

Service: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest

PubChem PUG REST is used as an external web service for chemical-name lookup.

PubChem is not bundled with this project and is not covered by this project's
MIT License. Use of the service is subject to PubChem / NCBI policies.

PubChem asks applications not to exceed 5 requests per second and may apply
dynamic throttling. The application should therefore avoid unnecessary
repeated or bulk requests.

Recommended citation:

S. Kim, P. A. Thiessen, T. Cheng, B. Yu, and E. E. Bolton.
"An update on PUG-REST: RESTful interface for programmatic access to PubChem."
Nucleic Acids Research 46(W1), W563-W570 (2018).
DOI: 10.1093/nar/gky294.

## Other Python Dependencies

This project also depends on packages such as:

- NumPy
- pandas
- Pillow
- requests

These packages retain their own licenses and copyright notices and are not
relicensed by this project.

For redistribution, review the exact package versions installed in your
environment and preserve any notices required by their corresponding licenses.

## Important Note

This file is an attribution and license-notice summary for the project. It is
not legal advice. When distributing third-party binaries or bundling libraries
directly, review the exact license text shipped with the specific version you
distribute.
