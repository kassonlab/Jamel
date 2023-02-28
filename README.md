# coronavirus-chimera-prediction
Code for generating and predicting coronavirus chimeras.

This repository is described in Simpson and Kasson, "Structural prediction of chimeric immunogens to elicit targeted antibodies against betacoronaviruses", doi:10.1101/2023.01.31.526494.

Here, we provide a set of Python functions and scripts to do the following:
-  Create arbitrary chimeras between protein sequences
-  Predict structure and pLDDT of these chimeras using AlphaFold
-  Score the resulting chimeras for stability using the AlphaFold outputs
-  Utilize the above functions to create S1/S2 chimeras between SARS-CoV-2 and a large set of betacoronaviruses
-  Computationally validate chimera predictions using molecular dynamics simulations.

## Functionality of specific python files is summarized below:

ProductionScript.py creates chimeras between a reference sequence and a user-supplied list of sequences.

MultimerAlphaFold.sh is a shell script to run AlphaFold

ChimeraAnalysis.py parses AlphaFold outputs to create creates pLDDT files and ultimately stability scores.

md_sims/setup.py will set up molecular dynamics simulations using Gromacs.  This requires [Gromacs](https://gitlab.com/gromacs/gromacs) as well as the [CHARMM36](http://mackerell.umaryland.edu/charmm_ff.shtml#gromacs) forcefield package.
As currently formulated, setup.py will take all PDB files in the current directory beginning with "3mer" and prepare production Gromacs run input files (parsing, energy minimization, and equilibration).

## The following files contain utility routines

AccessiontoAlignment.py contains functions necesary for:
- creating fasta files for the sequences attached to the accession numbers you've collected,
- creating a concatenated fasta file for multiple sequence alignment (MSA),
- using that MSA to find homologous sequences useful for splicing.

ChimeraGenerator.py contains functions that splice and then recombine sequence segments of your choice 

Analysis.py contains routines for analysis of AlphaFold outputs

## More detailed information:
**ProductionScript.py** creates chimeras between a reference sequence (ex:SARS-CoV-2) and a user-supplied list of accession numbers and names preferred names for the corresponding sequences.

Example: python3 ProductionScript.py List_of_Coronaviruses.tsv chimera_arguments.json 6vsb_S1.fasta
will create a set of chimeric fasta files along with the fastas for their wild-type parent proteins.  6vsb_S1.fasta is a fasta file containing the sequence in the reference protein that should be replaced by corresponding chimeric sequences.

Naming is specified in the chimera_arguments.json settings file.  Wildcards in the json file will be replaced by the short name of each sequence
given in the right-hand column of List_of_Coronaviruses.tsv.  This is also how the aligned sequences are identified in the multiple sequence alignment. If you input your own MSA, make sure these are identical. AlphaFold will name output folders after the basename of the input fasta file.

In both argument json files some variables are given as lists with a leading empty string ("").
This is designed to facilitate disabling selected commands by changing the "" to "#".

Example:
"emboss_command":["","/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle"] -> "emboss_command":["#","/scratch/jws6pq/EMBOSS-6.6.0/emboss/needle"]
will prevent emboss from running


**ChimeraAnalysis.py** creates pLDDT files from AlphaFold outputs.  These files contain the per residue confidence score arrays generated by the highest scoring prediction model from AlphaFold for each native and chimeric protein sequences.
The script then generates a csv data table where each row is a predicted chimera.
Data columns include:
- Preferred name of protein splice partner
- Emboss sequence similarity of swapped sequences
- Overall AlphaFold confidence score of the wild-type splice partner protein
- Overall confidence score of the resulting chimera with the reference protein
- Averaged relative stability of the chimera as a percentage.
Each of these columns can be turned off if desired.

ChimeraAnalysis.py has similar inputs on the commandline as the ProductionScript.py.

Example: python3 ChimeraAnalysis.py List_of_Coronaviruses.tsv analysis_arguments.json 6vsb_S1.fasta

Similar syntax applies, and the 6vsb_S1.fasta again contains the sequence in the reference protein that should be replaced by chimeric sequences.