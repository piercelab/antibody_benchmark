# Antibody Docking and Affinity Benchmark
This repository contains the antibody-antigen test cases, with unbound and bound structures, for Docking Benchmark 5.5. This dataset is a major update of Docking Benchmark 5.0, which was released in 2015 (Vreven et al. "Updates to the Integrated Protein-Protein Interaction Benchmarks: Docking Benchmark Version 5 and Affinity Benchmark Version 2" J Mol Biol 427(19):3031-41). Users can download this set to test and benchmark their predictive algorithms.

This update contains antibody-antigen structures from 67 test cases, more than doubling the amount in the previous benchmark. The citation for this benchmark is:
Guest JD, Vreven T, Zhou J, Moal I, Jeliazkov JR, Gray JJ, Weng Z, Pierce BG. "An Expanded Benchmark for Antibody-Antigen Docking and Affinity Prediction Reveals Insights into Antibody Recognition Determinants", Under Review.

# Nomenclature
Each test case is represented by four pdb files, the nomenclature for which is as follows:

'*complex-pdb-code*_r_u.pdb' - Unbound antibody structure

'*complex-pdb-code*_l_u.pdb' - Unbound antigen structure

'*complex-pdb-code*_r_b.pdb' - Bound antibody structure

'*complex-pdb-code*_l_b.pdb' - Bound antigen structure

Bound structures originate from the same PDB as the case name, yet unbound structures were taken from separate PDBs that correspond to the bound complex. For instance, 1AHW_r_b.pdb and 1AHW_l_b.pdb are from the antibody-antigen complex in 1AHW, but 1AHW_r_u.pdb is a structure from 1FGN and 1AHW_l_u.pdb is a structure from 1TFH. Bound-unbound pairs come pre-aligned for easy visualization of conformational changes. If testing a docking algorithm that is sensitive to initial positioning of unbound structures, users may randomize unbound structure positions to avoid possible bias in docking results.

# Cases
Information on the cases, including docking difficulties, conformational changes, and binding affinities, is provided in these tab-delimited tables: antibody_antigen_cases.txt and antibody_antigen_affinities.txt.

Additional information for a columns in antibody_antigen_cases.txt is below:

Complex PDB/Antibody PDB/Antigen PDB: PDB code is followed by IDs for antibody and antigen chains. For complexes, antibody chains are listed first and separated from antigen chains by a colon.
Antibody: Trade names for therapeutic antibodies in test cases are shown in parentheses.
I-RMSD: Interface RMSD, which helped to assign docking difficulty level, was calculated by superposition of unbound antibody and antigen structures onto the bound complex structure using root-mean-square fit of interface residues.                                                            
ΔASA: Measured change in accessible surface area upon complex formation.


Shield: [![CC BY 4.0][cc-by-shield]][cc-by]

This work is licensed under a [Creative Commons Attribution 4.0 International
License][cc-by].

[![CC BY 4.0][cc-by-image]][cc-by]

[cc-by]: http://creativecommons.org/licenses/by/4.0/
[cc-by-image]: https://i.creativecommons.org/l/by/4.0/88x31.png
[cc-by-shield]: https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg
