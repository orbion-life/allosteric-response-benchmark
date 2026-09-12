# Data sources and rights

The replay uses bundled, unchanged public scientific records. Original source and copied-module SHA-256 hashes are in `source-receipts.json`; the final package manifest includes later receipts and figures. The label reference was copied before scoring but opened for contact extraction only after the prediction freeze.

| File | Primary source | Attribution / rights |
|---|---|---|
| `raw/1OPL.cif` | https://files.rcsb.org/download/1OPL.cif | PDB archive data, CC0 1.0; Nagar et al. (2003) |
| `raw/5MO4.cif` | https://files.rcsb.org/download/5MO4.cif | PDB archive data, CC0 1.0; Wylie et al. (2017) |
| `raw/1OPL-sifts.json` | https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/1opl | PDBe SIFTS; Dana et al. (2019); Velankar et al. (2013); EMBL-EBI terms |
| `raw/5MO4-sifts.json` | https://www.ebi.ac.uk/pdbe/api/mappings/uniprot/5mo4 | Same SIFTS attribution and terms |
| `raw/P00519.json` | https://rest.uniprot.org/uniprotkb/P00519.json | UniProt Consortium; CC BY 4.0 |
| `raw/ABL-input.json`, domain PDB and residue map | Derived from the bundled 1OPL/SIFTS/UniProt sources | Project Pulsar preparation; retain upstream attribution |

The PDB archive is distributed under CC0 1.0, as described in the [RCSB policies](https://www.rcsb.org/pages/policies). [UniProt copyrightable database content uses CC BY 4.0](https://www.uniprot.org/help/license/); the official licence-help JSON is cached in `raw/UniProt-license.json`. [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/) impose no additional restrictions on contributed data and request resource attribution. The MIT code licence does not replace these data terms. These policies were checked on 12 September 2026. Neither the data providers nor the structure authors endorse this pilot.

Nagar, B. et al. (2003). Structural basis for the autoinhibition of c-Abl tyrosine kinase. *Cell* 112, 859–871. [doi:10.1016/S0092-8674(03)00194-6](https://doi.org/10.1016/S0092-8674(03)00194-6). Structure [1OPL](https://doi.org/10.2210/pdb1OPL/pdb).

Wylie, A. A. et al. (2017). The allosteric inhibitor ABL001 enables dual targeting of BCR-ABL1. *Nature* 543, 733–737. [doi:10.1038/nature21702](https://doi.org/10.1038/nature21702). Structure [5MO4](https://doi.org/10.2210/pdb5MO4/pdb).

Dana et al. (2019), *Nucleic Acids Research* 47, D482; and Velankar et al. (2013), *Nucleic Acids Research* 41, D483. Preferred citations and links are provided by the [SIFTS resource](https://www.ebi.ac.uk/pdbe/docs/sifts/). The bundled mappings are used as data, not copied implementation code.

RSA in the frozen input preparation uses a 1.4 Å Shrake–Rupley probe with 960 points and the theoretical maxima in Tien et al. (2013), [doi:10.1371/journal.pone.0080635](https://doi.org/10.1371/journal.pone.0080635). This pilot verifies the frozen geometry and masks; it does not rerun the upstream RSA calculation. The residue map retains the underlying SASA and RSA values.

The summary structure plot is generated directly from 1OPL Cα coordinates using Matplotlib. It does not reproduce a publisher image, use Mol*, or superpose the AY7 ligand into unmatched input coordinates. Original explanatory text, scripts and generated figures are provided under the package MIT licence, subject to the upstream data attribution above.
