# Introduction
This directory contains scripts and input files used to identify new antibody-antigen test cases from the PDB. The files are:
Scripts:
get_antibody_benchmark.pl: a Perl script to search for two-chain (heavy and light chain) antibody cases
get_antibody_benchmark_camelids.pl: a Perl script to search for single-chain antibody cases
get_cases.pl: a Perl script that processes the output from one of the above scripts to get sets of triplets (bound and unbound PDB codes) representing new cases
Data files:
example/bm5.5.txt: a list of cases in BM5.5, for use by the get_cases.pl script to not output redundant cases with existing (BM5.5) cases
example/all_antibody_structures_8_2020.txt: a text file of all antibody structures, downloaded from SAbDab (http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/) (Dunbar et al. Nucleic Acids Research 2014)
example/pdb_resolutions_8_2020.txt: a text file of PDB resolution information, downloaded from the PDB (https://www.rcsb.org/pages/general/summaries)
example/pdb_seqres_8_2020.fa.tar.gz: a tarball of all PDB sequences, downloaded from the PDB (https://www.rcsb.org/pages/general/summaries)

# Running
Prerequisite: make sure that the NCBI "blastp" executable is installed and available on your system. Check that the paths in the Perl scripts are in accordance wiht the data file paths on your system.
1. Unzip the "pdb_seqres_8_2020.fa.tar.gz" tarball in the example directory.
2. Run the get_antibody_benchmark.pl script, or the get_antibody_benchmark_camelids.pl script. For example "./get_antibody_benchmark.pl > bound_unbound_pdbs.txt".
3. Run the get_cases.pl script on that output. For example: "./get_cases.pl bound_unbound_pdbs.txt > cases.txt".

Example output files from step #2 are provided in the example direcotry for reference (bound_unbound_pdbs.txt, bound_unbound_pdbs.camelids.txt).
