# Data and software attribution

The PDB data files are distributed under [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/) as described by [RCSB PDB](https://www.rcsb.org/pages/policies). Sequence and annotation records retain the UniProt Consortium's [CC BY 4.0 terms](https://www.uniprot.org/help/license). Derived mappings and domain files indicate their input accessions and changes; no source endorsement is implied. PDBe SIFTS provides the structure-to-sequence mappings under [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use/), which add no restrictions to original data owners' terms and request attribution. These provider terms do not automatically relicense third-party components.

Use of this package should cite the relevant structure deposition and its original publication, the [Protein Data Bank (Berman et al., 2000)](https://doi.org/10.1093/nar/28.1.235), [PDBe SIFTS](https://www.ebi.ac.uk/pdbe/docs/sifts/overview.html), and [UniProt](https://www.uniprot.org/). Sources were retrieved on 12 September 2026. The complete files are cached for reproducibility; retained domains and identifiers are derived by the provided preparation source.

## Original structure sources

- [1NKP](https://www.rcsb.org/structure/1NKP): X-ray structures of Myc-Max and Mad-Max recognizing DNA: Molecular bases of regulation by proto-oncogenic transcription factors (2003). [Original source](https://doi.org/10.1016/S0092-8674(02)01284-9).
- [1OPL](https://www.rcsb.org/structure/1OPL): Structural basis for the autoinhibition of c-Abl tyrosine kinase (2003). [Original source](https://doi.org/10.1016/S0092-8674(03)00194-6).
- [4OBE](https://www.rcsb.org/structure/4OBE): In situ selectivity profiling and crystal structure of SML-8-73-1, an active site inhibitor of oncogenic K-Ras G12C. (2014). [Original source](https://doi.org/10.1073/pnas.1404639111).
- [5MO4](https://www.rcsb.org/structure/5MO4): The allosteric inhibitor ABL001 enables dual targeting of BCR-ABL1. (2017). [Original source](https://doi.org/10.1038/nature21702).
- [5TBY](https://www.rcsb.org/structure/5TBY): Effects of myosin variants on interacting-heads motif explain distinct hypertrophic and dilated cardiomyopathy phenotypes. (2017). [Original source](https://doi.org/10.7554/eLife.24634).
- [6C1H](https://www.rcsb.org/structure/6C1H): High-resolution cryo-EM structures of actin-bound myosin states reveal the mechanism of myosin force sensing. (2018). [Original source](https://doi.org/10.1073/pnas.1718316115).
- [6OIM](https://www.rcsb.org/structure/6OIM): The clinical KRAS(G12C) inhibitor AMG 510 drives anti-tumour immunity. (2019). [Original source](https://doi.org/10.1038/s41586-019-1694-1).

## Methods and software

- [Wang et al. (2020)](https://doi.org/10.1038/s41467-020-17618-2) introduced Ohm. Upstream source and this portability patch retain GNU GPL version 3; the license is in `ohm/LICENSE-GPL-3.0`. The original notices and archive accompany the sibling pilot wrapper.
- [Greener and Sternberg (2015)](https://doi.org/10.1186/s12859-015-0771-1) provides the AlloPred stiffness-perturbation precedent and its caution about unannotated negative pockets. AlloPred has not been executed in this package.
- [Tien et al. (2013)](https://doi.org/10.1371/journal.pone.0080635) provides the theoretical amino-acid maximum areas used to normalize exposure. The numerical table is transcribed in `build_validation_inputs.py`; no article PDF is redistributed.
- [Biopython](https://biopython.org/docs/latest/api/Bio.PDB.SASA.html) supplies the parser and Shrake–Rupley implementation; [NumPy](https://numpy.org/) supplies arrays and the documented PCG64 generator. Dependencies retain their own distributions and licenses.
- [Holm (1979)](https://www.ime.usp.br/~abe/lista/pdf4R8xPVzCnX.pdf) supplies the family correction. [Phipson and Smyth (2010)](https://gksmyth.github.io/pubs/PermPValuesPreprint.pdf) discusses valid discrete permutation p-values.

Project-authored source and documentation are distributed under the containing benchmark repository's stated license. That license does not replace the provider terms above or the GPL terms for Ohm and its patch. The public protocol contains no third-party executable binaries.
